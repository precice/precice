#include "acceleration/config/AccelerationConfiguration.hpp"
#include <algorithm>
#include <memory>
#include <ostream>
#include <stdexcept>
#include <utility>
#include <vector>
#include "acceleration/Acceleration.hpp"
#include "acceleration/AitkenAcceleration.hpp"
#include "acceleration/BroydenAcceleration.hpp"
#include "acceleration/ConstantRelaxationAcceleration.hpp"
#include "acceleration/IQNILSAcceleration.hpp"
#include "acceleration/MVQNAcceleration.hpp"
#include "acceleration/impl/ConstantPreconditioner.hpp"
#include "acceleration/impl/QRFactorization.hpp"
#include "acceleration/impl/ResidualPreconditioner.hpp"
#include "acceleration/impl/ResidualSumPreconditioner.hpp"
#include "acceleration/impl/ValuePreconditioner.hpp"
#include "logging/LogMacros.hpp"
#include "mesh/Data.hpp"
#include "mesh/Mesh.hpp"
#include "mesh/config/MeshConfiguration.hpp"
#include "utils/assertion.hpp"
#include "xml/ConfigParser.hpp"
#include "xml/XMLAttribute.hpp"
#include "xml/XMLTag.hpp"

namespace precice {
namespace acceleration {

using namespace precice::acceleration::impl;

AccelerationConfiguration::AccelerationConfiguration(
    const mesh::PtrMeshConfiguration &meshConfig)
    : TAG("acceleration"),
      TAG_RELAX("relaxation"),
      TAG_INIT_RELAX("initial-relaxation"),
      TAG_MAX_USED_ITERATIONS("max-used-iterations"),
      TAG_TIME_WINDOWS_REUSED("time-windows-reused"),
      TAG_DATA("data"),
      TAG_FILTER("filter"),
      TAG_ESTIMATEJACOBIAN("estimate-jacobian"),
      TAG_PRECONDITIONER("preconditioner"),
      TAG_IMVJRESTART("imvj-restart-mode"),
      ATTR_NAME("name"),
      ATTR_MESH("mesh"),
      ATTR_SCALING("scaling"),
      ATTR_VALUE("value"),
      ATTR_ENFORCE("enforce"),
      ATTR_SINGULARITYLIMIT("limit"),
      ATTR_TYPE("type"),
      ATTR_BUILDJACOBIAN("always-build-jacobian"),
      ATTR_IMVJCHUNKSIZE("chunk-size"),
      ATTR_RSLS_REUSED_TIME_WINDOWS("reused-time-windows-at-restart"),
      ATTR_RSSVD_TRUNCATIONEPS("truncation-threshold"),
      ATTR_PRECOND_NONCONST_TIME_WINDOWS("freeze-after"),
      VALUE_CONSTANT("constant"),
      VALUE_AITKEN("aitken"),
      VALUE_IQNILS("IQN-ILS"),
      VALUE_MVQN("IQN-IMVJ"),
      VALUE_BROYDEN("broyden"),
      VALUE_QR1FILTER("QR1"),
      VALUE_QR1_ABSFILTER("QR1-absolute"),
      VALUE_QR2FILTER("QR2"),
      VALUE_CONSTANT_PRECONDITIONER("constant"),
      VALUE_VALUE_PRECONDITIONER("value"),
      VALUE_RESIDUAL_PRECONDITIONER("residual"),
      VALUE_RESIDUAL_SUM_PRECONDITIONER("residual-sum"),
      VALUE_LS_RESTART("RS-LS"),
      VALUE_ZERO_RESTART("RS-0"),
      VALUE_SVD_RESTART("RS-SVD"),
      VALUE_SLIDE_RESTART("RS-SLIDE"),
      VALUE_NO_RESTART("no-restart"),
      _meshConfig(meshConfig),
      _acceleration(),
      _neededMeshes(),
      _preconditioner(),
      _config()
{
  PRECICE_ASSERT(meshConfig.get() != nullptr);
}

void AccelerationConfiguration::connectTags(xml::XMLTag &parent)
{
  using namespace xml;

  // static int recursionCounter = 0;
  // recursionCounter++;

  XMLTag::Occurrence  occ = XMLTag::OCCUR_NOT_OR_ONCE;
  std::vector<XMLTag> tags;
  {
    XMLTag tag(*this, VALUE_CONSTANT, occ, TAG);
    tag.setDocumentation("Accelerates coupling data with constant underrelaxation.");
    addTypeSpecificSubtags(tag);
    tags.push_back(tag);
  }
  {
    XMLTag tag(*this, VALUE_AITKEN, occ, TAG);
    tag.setDocumentation("Accelerates coupling data with dynamic Aitken under-relaxation.");
    addTypeSpecificSubtags(tag);
    tags.push_back(tag);
  }
  {
    XMLTag tag(*this, VALUE_IQNILS, occ, TAG);
    tag.setDocumentation("Accelerates coupling data with the interface quasi-Newton inverse least-squares method.");
    addTypeSpecificSubtags(tag);
    tags.push_back(tag);
  }
  {
    XMLTag tag(*this, VALUE_MVQN, occ, TAG);
    tag.setDocumentation("Accelerates coupling data with the interface quasi-Newton inverse multi-vector Jacobian method.");

    auto alwaybuildJacobian = makeXMLAttribute(ATTR_BUILDJACOBIAN, false)
                                  .setDocumentation("If set to true, the IMVJ will set up the Jacobian matrix"
                                                    " in each coupling iteration, which is inefficient. If set to false (or not set)"
                                                    " the Jacobian is only build in the last iteration and the updates are computed using (relatively) cheap MATVEC products.");
    tag.addAttribute(alwaybuildJacobian);

    addTypeSpecificSubtags(tag);
    tags.push_back(tag);
  }
  {
    XMLTag tag(*this, VALUE_BROYDEN, occ, TAG);
    tag.setDocumentation("Accelerates coupling data with the (single-vector) Broyden method.");
    addTypeSpecificSubtags(tag);
    tags.push_back(tag);
  }

  for (XMLTag &tag : tags) {
    parent.addSubtag(tag);
  }
}

PtrAcceleration AccelerationConfiguration::getAcceleration()
{
  return _acceleration;
}

void AccelerationConfiguration::xmlTagCallback(
    const xml::ConfigurationContext &context,
    xml::XMLTag &                    callingTag)
{
  PRECICE_TRACE(callingTag.getFullName());

  if (callingTag.getNamespace() == TAG) {
    _config.type = callingTag.getName();

    if (_config.type == VALUE_MVQN)
      _config.alwaysBuildJacobian = callingTag.getBooleanAttributeValue(ATTR_BUILDJACOBIAN);
  }

  if (callingTag.getName() == TAG_RELAX) {
    _config.relaxationFactor = callingTag.getDoubleAttributeValue(ATTR_VALUE);
  } else if (callingTag.getName() == TAG_DATA) {
    std::string dataName = callingTag.getStringAttributeValue(ATTR_NAME);
    std::string meshName = callingTag.getStringAttributeValue(ATTR_MESH);
    auto        success  = _uniqueDataAndMeshNames.emplace(dataName, meshName);
    if (not success.second) {
      PRECICE_ERROR("You have provided a subtag "
                    << "<data name=\"" << dataName << "\" mesh=\"" << meshName << "\"/> more than once in your <acceleration:.../>. "
                    << "Please remove the duplicated entry.");
    }
    _meshName      = callingTag.getStringAttributeValue(ATTR_MESH);
    double scaling = 1.0;
    if (_config.type == VALUE_IQNILS || _config.type == VALUE_MVQN || _config.type == VALUE_BROYDEN) {
      scaling = callingTag.getDoubleAttributeValue(ATTR_SCALING);
    }

    for (mesh::PtrMesh mesh : _meshConfig->meshes()) {
      if (mesh->getName() == _meshName) {
        for (mesh::PtrData data : mesh->data()) {
          if (dataName == data->getName()) {
            _config.dataIDs.push_back(data->getID());
            _config.scalings.insert(std::make_pair(data->getID(), scaling));
          }
        }
      }
    }

    if (_config.dataIDs.empty()) {
      std::ostringstream stream;
      stream << "Data with name \"" << dataName << "\" associated to mesh \""
             << _meshName << "\" not found on configuration of acceleration. Add \""
             << dataName << "\" to the \"<mesh name=" << _meshName << ">\" tag, or change the "
             << "data name in the acceleration scheme.";
      throw std::runtime_error{stream.str()};
    }
    _neededMeshes.push_back(_meshName);
  } else if (callingTag.getName() == TAG_INIT_RELAX) {
    _config.relaxationFactor = callingTag.getDoubleAttributeValue(ATTR_VALUE);
    if (callingTag.hasAttribute(ATTR_ENFORCE)) {
      _config.forceInitialRelaxation = callingTag.getBooleanAttributeValue(ATTR_ENFORCE);
    } else {
      _config.forceInitialRelaxation = false;
    }
  } else if (callingTag.getName() == TAG_MAX_USED_ITERATIONS) {
    _config.maxIterationsUsed = callingTag.getIntAttributeValue(ATTR_VALUE);
  } else if (callingTag.getName() == TAG_TIME_WINDOWS_REUSED) {
    _config.timeWindowsReused = callingTag.getIntAttributeValue(ATTR_VALUE);
  } else if (callingTag.getName() == TAG_FILTER) {
    auto f = callingTag.getStringAttributeValue(ATTR_TYPE);
    if (f == VALUE_QR1FILTER) {
      _config.filter = Acceleration::QR1FILTER;
    } else if (f == VALUE_QR1_ABSFILTER) {
      _config.filter = Acceleration::QR1FILTER_ABS;
    } else if (f == VALUE_QR2FILTER) {
      _config.filter = Acceleration::QR2FILTER;
    } else {
      PRECICE_ASSERT(false);
    }
    _config.singularityLimit = callingTag.getDoubleAttributeValue(ATTR_SINGULARITYLIMIT);
  } else if (callingTag.getName() == TAG_PRECONDITIONER) {
    _config.preconditionerType       = callingTag.getStringAttributeValue(ATTR_TYPE);
    _config.precond_nbNonConstTSteps = callingTag.getIntAttributeValue(ATTR_PRECOND_NONCONST_TIME_WINDOWS);
  } else if (callingTag.getName() == TAG_IMVJRESTART) {

    if (_config.alwaysBuildJacobian)
      PRECICE_ERROR("IMVJ cannot be in restart mode while parameter always-build-jacobian is set to true. "
                    << "Please remove 'always-build-jacobian' from the configuration file or do not run in restart mode.");

#ifndef PRECICE_NO_MPI
    _config.imvjChunkSize = callingTag.getIntAttributeValue(ATTR_IMVJCHUNKSIZE);
    auto f                = callingTag.getStringAttributeValue(ATTR_TYPE);
    if (f == VALUE_NO_RESTART) {
      _config.imvjRestartType = MVQNAcceleration::NO_RESTART;
    } else if (f == VALUE_ZERO_RESTART) {
      _config.imvjRestartType = MVQNAcceleration::RS_ZERO;
    } else if (f == VALUE_LS_RESTART) {
      _config.imvjRSLS_reusedTimeWindows = callingTag.getIntAttributeValue(ATTR_RSLS_REUSED_TIME_WINDOWS);
      _config.imvjRestartType            = MVQNAcceleration::RS_LS;
    } else if (f == VALUE_SVD_RESTART) {
      _config.imvjRSSVD_truncationEps = callingTag.getDoubleAttributeValue(ATTR_RSSVD_TRUNCATIONEPS);
      _config.imvjRestartType         = MVQNAcceleration::RS_SVD;
    } else if (f == VALUE_SLIDE_RESTART) {
      _config.imvjRestartType = MVQNAcceleration::RS_SLIDE;
    } else {
      _config.imvjChunkSize = 0;
      PRECICE_ASSERT(false);
    }
#else
    PRECICE_ERROR("Acceleration IQN-IMVJ only works if preCICE is compiled with MPI");
#endif
  }
}

void AccelerationConfiguration::xmlEndTagCallback(
    const xml::ConfigurationContext &context,
    xml::XMLTag &                    callingTag)
{
  PRECICE_TRACE(callingTag.getName());
  if (callingTag.getNamespace() == TAG) {

    //create preconditioner
    if (callingTag.getName() == VALUE_IQNILS || callingTag.getName() == VALUE_MVQN) {

      // if imvj restart-mode is of type RS-SVD, max number of non-const preconditioned time steps is limited by the chunksize
      if (callingTag.getName() == VALUE_MVQN && _config.imvjRestartType > 0)
        if (_config.precond_nbNonConstTSteps > _config.imvjChunkSize)
          _config.precond_nbNonConstTSteps = _config.imvjChunkSize;

      if (_config.preconditionerType == VALUE_CONSTANT_PRECONDITIONER) {
        std::vector<double> factors;
        for (int id : _config.dataIDs) {
          factors.push_back(_config.scalings[id]);
        }
        _preconditioner = PtrPreconditioner(new ConstantPreconditioner(factors));
      } else if (_config.preconditionerType == VALUE_VALUE_PRECONDITIONER) {
        _preconditioner = PtrPreconditioner(new ValuePreconditioner(_config.precond_nbNonConstTSteps));
      } else if (_config.preconditionerType == VALUE_RESIDUAL_PRECONDITIONER) {
        _preconditioner = PtrPreconditioner(new ResidualPreconditioner(_config.precond_nbNonConstTSteps));
      } else if (_config.preconditionerType == VALUE_RESIDUAL_SUM_PRECONDITIONER) {
        _preconditioner = PtrPreconditioner(new ResidualSumPreconditioner(_config.precond_nbNonConstTSteps));
      } else {
        // no preconditioner defined
        std::vector<double> factors;
        for (int id = 0; id < (int) _config.dataIDs.size(); ++id) {
          factors.push_back(1.0);
        }
        _preconditioner = PtrPreconditioner(new ConstantPreconditioner(factors));
      }
    }

    if (callingTag.getName() == VALUE_CONSTANT) {
      _acceleration = PtrAcceleration(
          new ConstantRelaxationAcceleration(
              _config.relaxationFactor, _config.dataIDs));
    } else if (callingTag.getName() == VALUE_AITKEN) {
      _acceleration = PtrAcceleration(
          new AitkenAcceleration(
              _config.relaxationFactor, _config.dataIDs));
    } else if (callingTag.getName() == VALUE_IQNILS) {
      _acceleration = PtrAcceleration(
          new IQNILSAcceleration(
              _config.relaxationFactor,
              _config.forceInitialRelaxation,
              _config.maxIterationsUsed,
              _config.timeWindowsReused,
              _config.filter, _config.singularityLimit,
              _config.dataIDs,
              _preconditioner));
    } else if (callingTag.getName() == VALUE_MVQN) {
#ifndef PRECICE_NO_MPI
      _acceleration = PtrAcceleration(
          new MVQNAcceleration(
              _config.relaxationFactor,
              _config.forceInitialRelaxation,
              _config.maxIterationsUsed,
              _config.timeWindowsReused,
              _config.filter, _config.singularityLimit,
              _config.dataIDs,
              _preconditioner,
              _config.alwaysBuildJacobian,
              _config.imvjRestartType,
              _config.imvjChunkSize,
              _config.imvjRSLS_reusedTimeWindows,
              _config.imvjRSSVD_truncationEps));
#else
      PRECICE_ERROR("Acceleration IQN-IMVJ only works if preCICE is compiled with MPI");
#endif
    } else if (callingTag.getName() == VALUE_BROYDEN) {
      _acceleration = PtrAcceleration(
          new BroydenAcceleration(
              _config.relaxationFactor,
              _config.forceInitialRelaxation,
              _config.maxIterationsUsed,
              _config.timeWindowsReused,
              _config.filter, _config.singularityLimit,
              _config.dataIDs,
              _preconditioner));
    } else {
      PRECICE_ASSERT(false);
    }
  }
}

void AccelerationConfiguration::clear()
{
  _config       = ConfigurationData();
  _acceleration = PtrAcceleration();
  _neededMeshes.clear();
}

void AccelerationConfiguration::addCommonIQNSubtags(xml::XMLTag &tag)
{
  using namespace precice::xml;

  XMLTag tagData(*this, TAG_DATA, XMLTag::OCCUR_ONCE_OR_MORE);
  tagData.setDocumentation("The data used to compute the acceleration.");
  XMLAttribute<std::string> attrName(ATTR_NAME);
  attrName.setDocumentation("The name of the data.");
  XMLAttribute<std::string> attrMesh(ATTR_MESH);
  attrMesh.setDocumentation("The name of the mesh which holds the data.");
  auto attrScaling = makeXMLAttribute(ATTR_SCALING, 1.0)
                         .setDocumentation(
                             "To improve the performance of a parallel or a multi coupling schemes, "
                             "data values can be manually scaled. We recommend, however, to use an automatic scaling via a preconditioner.");
  tagData.addAttribute(attrScaling);
  tagData.addAttribute(attrName);
  tagData.addAttribute(attrMesh);
  tag.addSubtag(tagData);

  XMLTag tagFilter(*this, TAG_FILTER, XMLTag::OCCUR_NOT_OR_ONCE);
  tagFilter.setDocumentation("Type of filtering technique that is used to "
                             "maintain good conditioning in the least-squares system. Possible filters:\n"
                             " - `QR1-filter`: updateQR-dec with (relative) test \\\\(R(i,i) < \\epsilon *\\lVert R\\rVert_F\\\\)\n"
                             " - `QR1_absolute-filter`: updateQR-dec with (absolute) test \\\\(R(i, i) < \\epsilon\\\\)\n"
                             " - `QR2-filter`: en-block QR-dec with test \\\\(\\lVert v_\\text{orth} \\rVert_2 < \\epsilon * \\lVert v \\rVert_2\\\\)\n\n"
                             "Please note that a QR1 is based on Given's rotations whereas QR2 uses "
                             "modified Gram-Schmidt. This can give different results even when no columns "
                             "are filtered out.");
  XMLAttribute<double> attrSingularityLimit(ATTR_SINGULARITYLIMIT, 1e-16);
  attrSingularityLimit.setDocumentation("Limit eps of the filter.");
  tagFilter.addAttribute(attrSingularityLimit);
  auto attrFilterName = XMLAttribute<std::string>(ATTR_TYPE)
                            .setOptions({VALUE_QR1FILTER,
                                         VALUE_QR1_ABSFILTER,
                                         VALUE_QR2FILTER})
                            .setDocumentation("Type of the filter.");
  tagFilter.addAttribute(attrFilterName);
  tag.addSubtag(tagFilter);
}

void AccelerationConfiguration::addTypeSpecificSubtags(
    xml::XMLTag &tag)
{
  using namespace xml;
  if (tag.getName() == VALUE_CONSTANT) {
    XMLTag               tagRelax(*this, TAG_RELAX, XMLTag::OCCUR_ONCE);
    XMLAttribute<double> attrValue(ATTR_VALUE);
    attrValue.setDocumentation("Constant relaxation factor.");
    tagRelax.addAttribute(attrValue);
    tag.addSubtag(tagRelax);
  } else if (tag.getName() == VALUE_AITKEN) {
    XMLTag tagInitRelax(*this, TAG_INIT_RELAX, XMLTag::OCCUR_ONCE);
    tagInitRelax.setDocumentation("Initial relaxation factor.");
    XMLAttribute<double> attrValue(ATTR_VALUE);
    attrValue.setDocumentation("Initial relaxation factor.");
    tagInitRelax.addAttribute(attrValue);
    tag.addSubtag(tagInitRelax);

    XMLTag tagData(*this, TAG_DATA, XMLTag::OCCUR_ONCE_OR_MORE);
    tagData.setDocumentation("The data used to compute the acceleration.");
    XMLAttribute<std::string> attrName(ATTR_NAME);
    attrName.setDocumentation("The name of the data.");
    XMLAttribute<std::string> attrMesh(ATTR_MESH);
    attrMesh.setDocumentation("The name of the mesh which holds the data.");
    tagData.addAttribute(attrName);
    tagData.addAttribute(attrMesh);
    tag.addSubtag(tagData);
  } else if (tag.getName() == VALUE_IQNILS) {
    XMLTag tagInitRelax(*this, TAG_INIT_RELAX, XMLTag::OCCUR_ONCE);
    tagInitRelax.setDocumentation("Initial relaxation factor.");
    XMLAttribute<double> attrDoubleValue(ATTR_VALUE);
    attrDoubleValue.setDocumentation("Initial relaxation factor.");
    tagInitRelax.addAttribute(attrDoubleValue);
    XMLAttribute<bool> attrEnforce(ATTR_ENFORCE, false);
    attrEnforce.setDocumentation("Enforce initial relaxation in every time window.");
    tagInitRelax.addAttribute(attrEnforce);
    tag.addSubtag(tagInitRelax);

    XMLTag tagMaxUsedIter(*this, TAG_MAX_USED_ITERATIONS, XMLTag::OCCUR_ONCE);
    tagMaxUsedIter.setDocumentation("Maximum number of columns used in low-rank approximation of Jacobian.");
    XMLAttribute<int> attrIntValue(ATTR_VALUE);
    attrIntValue.setDocumentation("The number of columns.");
    tagMaxUsedIter.addAttribute(attrIntValue);
    tag.addSubtag(tagMaxUsedIter);

    XMLTag tagTimeWindowsReused(*this, TAG_TIME_WINDOWS_REUSED, XMLTag::OCCUR_ONCE);
    tagTimeWindowsReused.setDocumentation("Number of past time windows from which columns are used to approximate Jacobian.");
    XMLAttribute<int> attrNumTimeWindowsReused(ATTR_VALUE);
    attrNumTimeWindowsReused.setDocumentation("The number of time windows.");
    tagTimeWindowsReused.addAttribute(attrNumTimeWindowsReused);
    tag.addSubtag(tagTimeWindowsReused);

    addCommonIQNSubtags(tag);

    XMLTag tagPreconditioner(*this, TAG_PRECONDITIONER, XMLTag::OCCUR_NOT_OR_ONCE);
    tagPreconditioner.setDocumentation("To improve the performance of a parallel or a multi coupling schemes a preconditioner"
                                       " can be applied. A constant preconditioner scales every acceleration data by a constant value, which you can define as"
                                       " an attribute of data. "
                                       " A value preconditioner scales every acceleration data by the norm of the data in the previous time window."
                                       " A residual preconditioner scales every acceleration data by the current residual."
                                       " A residual-sum preconditioner scales every acceleration data by the sum of the residuals from the current time window.");
    auto attrPreconditionerType = XMLAttribute<std::string>(ATTR_TYPE)
                                      .setOptions({VALUE_CONSTANT_PRECONDITIONER,
                                                   VALUE_VALUE_PRECONDITIONER,
                                                   VALUE_RESIDUAL_PRECONDITIONER,
                                                   VALUE_RESIDUAL_SUM_PRECONDITIONER})
                                      .setDocumentation("The type of the preconditioner.");
    tagPreconditioner.addAttribute(attrPreconditionerType);
    auto nonconstTSteps = makeXMLAttribute(ATTR_PRECOND_NONCONST_TIME_WINDOWS, -1)
                              .setDocumentation(
                                  "After the given number of time steps, the preconditioner weights "
                                  "are frozen and the preconditioner acts like a constant preconditioner.");
    tagPreconditioner.addAttribute(nonconstTSteps);
    tag.addSubtag(tagPreconditioner);

  } else if (tag.getName() == VALUE_MVQN) {
    XMLTag tagInitRelax(*this, TAG_INIT_RELAX, XMLTag::OCCUR_ONCE);
    tagInitRelax.setDocumentation("Initial relaxation factor.");
    tagInitRelax.addAttribute(
        XMLAttribute<double>(ATTR_VALUE).setDocumentation("Initial relaxation factor."));
    tagInitRelax.addAttribute(
        XMLAttribute<bool>(ATTR_ENFORCE, false).setDocumentation("Enforce initial relaxation in every time window."));
    tag.addSubtag(tagInitRelax);

    XMLTag tagIMVJRESTART(*this, TAG_IMVJRESTART, XMLTag::OCCUR_NOT_OR_ONCE);
    auto   attrRestartName = XMLAttribute<std::string>(ATTR_TYPE)
                               .setOptions({VALUE_NO_RESTART,
                                            VALUE_ZERO_RESTART,
                                            VALUE_LS_RESTART,
                                            VALUE_SVD_RESTART,
                                            VALUE_SLIDE_RESTART})
                               .setDefaultValue(VALUE_SVD_RESTART)
                               .setDocumentation("Type of the restart mode.");
    tagIMVJRESTART.addAttribute(attrRestartName);
    tagIMVJRESTART.setDocumentation("Type of IMVJ restart mode that is used:\n"
                                    "- `no-restart`: IMVJ runs in normal mode with explicit representation of Jacobian\n"
                                    "- `RS-ZERO`:    IMVJ runs in restart mode. After M time steps all Jacobain information is dropped, restart with no information\n"
                                    "- `RS-LS`:      IMVJ runs in restart mode. After M time steps a IQN-LS like approximation for the initial guess of the Jacobian is computed.\n"
                                    "- `RS-SVD`:     IMVJ runs in restart mode. After M time steps a truncated SVD of the Jacobian is updated.\n"
                                    "- `RS-SLIDE`:   IMVJ runs in sliding window restart mode.\n");
    auto attrChunkSize = makeXMLAttribute(ATTR_IMVJCHUNKSIZE, 8)
                             .setDocumentation("Specifies the number of time steps M after which the IMVJ restarts, if run in restart-mode. Defaul value is M=8.");
    auto attrReusedTimeWindowsAtRestart = makeXMLAttribute(ATTR_RSLS_REUSED_TIME_WINDOWS, 8)
                                              .setDocumentation("If IMVJ restart-mode=RS-LS, the number of reused time steps at restart can be specified.");
    auto attrRSSVD_truncationEps = makeXMLAttribute(ATTR_RSSVD_TRUNCATIONEPS, 1e-4)
                                       .setDocumentation("If IMVJ restart-mode=RS-SVD, the truncation threshold for the updated SVD can be set.");
    tagIMVJRESTART.addAttribute(attrChunkSize);
    tagIMVJRESTART.addAttribute(attrReusedTimeWindowsAtRestart);
    tagIMVJRESTART.addAttribute(attrRSSVD_truncationEps);
    tag.addSubtag(tagIMVJRESTART);

    XMLTag tagMaxUsedIter(*this, TAG_MAX_USED_ITERATIONS, XMLTag::OCCUR_ONCE);
    tagMaxUsedIter.setDocumentation("Maximum number of columns used in low-rank approximation of Jacobian.");
    XMLAttribute<int> attrIntValue(ATTR_VALUE);
    attrIntValue.setDocumentation("The number of columns.");
    tagMaxUsedIter.addAttribute(attrIntValue);
    tag.addSubtag(tagMaxUsedIter);

    XMLTag tagTimeWindowsReused(*this, TAG_TIME_WINDOWS_REUSED, XMLTag::OCCUR_ONCE);
    tagTimeWindowsReused.setDocumentation("Number of past time windows from which columns are used to approximate Jacobian.");
    tagTimeWindowsReused.addAttribute(attrIntValue);
    tag.addSubtag(tagTimeWindowsReused);

    addCommonIQNSubtags(tag);

    XMLTag tagPreconditioner(*this, TAG_PRECONDITIONER, XMLTag::OCCUR_NOT_OR_ONCE);
    tagPreconditioner.setDocumentation(
        "To improve the performance of a parallel or a multi coupling schemes a preconditioner"
        " can be applied. A constant preconditioner scales every acceleration data by a constant value, which you can define as"
        " an attribute of data.\n"
        "- A value preconditioner scales every acceleration data by the norm of the data in the previous time window.\n"
        "- A residual preconditioner scales every acceleration data by the current residual.\n"
        "- A residual-sum preconditioner scales every acceleration data by the sum of the residuals from the current time window.\n");
    auto attrPreconditionerType = XMLAttribute<std::string>(ATTR_TYPE)
                                      .setOptions({VALUE_CONSTANT_PRECONDITIONER,
                                                   VALUE_VALUE_PRECONDITIONER,
                                                   VALUE_RESIDUAL_PRECONDITIONER,
                                                   VALUE_RESIDUAL_SUM_PRECONDITIONER})
                                      .setDocumentation("Type of the preconditioner.");
    tagPreconditioner.addAttribute(attrPreconditionerType);
    auto nonconstTSteps = makeXMLAttribute(ATTR_PRECOND_NONCONST_TIME_WINDOWS, -1)
                              .setDocumentation("After the given number of time steps, the preconditioner weights are frozen and the preconditioner acts like a constant preconditioner.");
    tagPreconditioner.addAttribute(nonconstTSteps);
    tag.addSubtag(tagPreconditioner);

  } else if (tag.getName() == VALUE_BROYDEN) {
    XMLTag               tagInitRelax(*this, TAG_INIT_RELAX, XMLTag::OCCUR_ONCE);
    XMLAttribute<double> attrDoubleValue(ATTR_VALUE);
    tagInitRelax.addAttribute(attrDoubleValue);
    XMLAttribute<bool> attrEnforce(ATTR_ENFORCE, false);
    tagInitRelax.addAttribute(attrEnforce);
    tag.addSubtag(tagInitRelax);

    XMLTag            tagMaxUsedIter(*this, TAG_MAX_USED_ITERATIONS, XMLTag::OCCUR_ONCE);
    XMLAttribute<int> attrIntValue(ATTR_VALUE);
    tagMaxUsedIter.addAttribute(attrIntValue);
    tag.addSubtag(tagMaxUsedIter);

    XMLTag tagTimeWindowsReused(*this, TAG_TIME_WINDOWS_REUSED, XMLTag::OCCUR_ONCE);
    tagTimeWindowsReused.addAttribute(attrIntValue);
    tag.addSubtag(tagTimeWindowsReused);

    XMLTag                    tagData(*this, TAG_DATA, XMLTag::OCCUR_ONCE_OR_MORE);
    XMLAttribute<std::string> attrName(ATTR_NAME);
    XMLAttribute<std::string> attrMesh(ATTR_MESH);
    auto                      attrScaling = makeXMLAttribute(ATTR_SCALING, 1.0)
                           .setDocumentation(
                               "To improve the performance of a parallel or a multi coupling schemes, "
                               "data values can be manually scaled. We recommend, however, to use an automatic scaling via a preconditioner.");
    tagData.addAttribute(attrScaling);
    tagData.addAttribute(attrName);
    tagData.addAttribute(attrMesh);
    tag.addSubtag(tagData);
  } else {
    PRECICE_ERROR("Acceleration of type \""
                  << tag.getName() << "\" is unknown. Please choose a valid acceleration scheme or check the spelling in the configuration file.");
  }
}
} // namespace acceleration
} // namespace precice
