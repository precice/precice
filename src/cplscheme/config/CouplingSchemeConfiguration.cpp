#include "CouplingSchemeConfiguration.hpp"
#include <algorithm>
#include <memory>
#include <ostream>
#include <stddef.h>
#include <stdexcept>
#include <utility>
#include <vector>
#include "acceleration/Acceleration.hpp"
#include "acceleration/config/AccelerationConfiguration.hpp"
#include "cplscheme/BaseCouplingScheme.hpp"
#include "cplscheme/BiCouplingScheme.hpp"
#include "cplscheme/CompositionalCouplingScheme.hpp"
#include "cplscheme/MultiCouplingScheme.hpp"
#include "cplscheme/ParallelCouplingScheme.hpp"
#include "cplscheme/SerialCouplingScheme.hpp"
#include "cplscheme/SharedPointer.hpp"
#include "cplscheme/impl/AbsoluteConvergenceMeasure.hpp"
#include "cplscheme/impl/MinIterationConvergenceMeasure.hpp"
#include "cplscheme/impl/RelativeConvergenceMeasure.hpp"
#include "cplscheme/impl/ResidualRelativeConvergenceMeasure.hpp"
#include "logging/LogMacros.hpp"
#include "m2n/SharedPointer.hpp"
#include "m2n/config/M2NConfiguration.hpp"
#include "mesh/Data.hpp"
#include "mesh/Mesh.hpp"
#include "mesh/config/MeshConfiguration.hpp"
#include "precice/impl/SharedPointer.hpp"
#include "utils/Helpers.hpp"
#include "utils/assertion.hpp"
#include "xml/ConfigParser.hpp"
#include "xml/XMLAttribute.hpp"
#include "xml/XMLTag.hpp"

namespace precice {
namespace cplscheme {

using precice::impl::PtrParticipant;

CouplingSchemeConfiguration::CouplingSchemeConfiguration(
    xml::XMLTag &                               parent,
    const mesh::PtrMeshConfiguration &          meshConfig,
    const m2n::M2NConfiguration::SharedPointer &m2nConfig)
    : TAG("coupling-scheme"),
      TAG_PARTICIPANTS("participants"),
      TAG_PARTICIPANT("participant"),
      TAG_EXCHANGE("exchange"),
      TAG_MAX_TIME("max-time"),
      TAG_MAX_TIME_WINDOWS("max-time-windows"),
      TAG_TIME_WINDOW_SIZE("time-window-size"),
      TAG_ABS_CONV_MEASURE("absolute-convergence-measure"),
      TAG_REL_CONV_MEASURE("relative-convergence-measure"),
      TAG_RES_REL_CONV_MEASURE("residual-relative-convergence-measure"),
      TAG_MIN_ITER_CONV_MEASURE("min-iteration-convergence-measure"),
      TAG_MAX_ITERATIONS("max-iterations"),
      TAG_EXTRAPOLATION("extrapolation-order"),
      ATTR_DATA("data"),
      ATTR_MESH("mesh"),
      ATTR_PARTICIPANT("participant"),
      ATTR_INITIALIZE("initialize"),
      ATTR_TYPE("type"),
      ATTR_FIRST("first"),
      ATTR_SECOND("second"),
      ATTR_VALUE("value"),
      ATTR_VALID_DIGITS("valid-digits"),
      ATTR_METHOD("method"),
      ATTR_LIMIT("limit"),
      ATTR_MIN_ITERATIONS("min-iterations"),
      ATTR_NAME("name"),
      ATTR_FROM("from"),
      ATTR_TO("to"),
      ATTR_SUFFICES("suffices"),
      ATTR_STRICT("strict"),
      ATTR_CONTROL("control"),
      VALUE_SERIAL_EXPLICIT("serial-explicit"),
      VALUE_PARALLEL_EXPLICIT("parallel-explicit"),
      VALUE_SERIAL_IMPLICIT("serial-implicit"),
      VALUE_PARALLEL_IMPLICIT("parallel-implicit"),
      VALUE_MULTI("multi"),
      VALUE_FIXED("fixed"),
      VALUE_FIRST_PARTICIPANT("first-participant"),
      _config(),
      _meshConfig(meshConfig),
      _m2nConfig(m2nConfig),
      _couplingSchemes(),
      _couplingSchemeCompositions()
{
  using namespace xml;

  XMLTag::Occurrence  occ = XMLTag::OCCUR_ARBITRARY;
  std::vector<XMLTag> tags;
  {
    XMLTag tag(*this, VALUE_SERIAL_EXPLICIT, occ, TAG);
    tag.setDocumentation("Explicit coupling scheme according to conventional serial staggered procedure (CSS).");
    addTypespecifcSubtags(VALUE_SERIAL_EXPLICIT, tag);
    tags.push_back(tag);
  }
  {
    XMLTag tag(*this, VALUE_PARALLEL_EXPLICIT, occ, TAG);
    tag.setDocumentation("Explicit coupling scheme according to conventional parallel staggered procedure (CPS).");
    addTypespecifcSubtags(VALUE_PARALLEL_EXPLICIT, tag);
    tags.push_back(tag);
  }
  {
    XMLTag tag(*this, VALUE_SERIAL_IMPLICIT, occ, TAG);
    tag.setDocumentation("Implicit coupling scheme according to block Gauss-Seidel iterations (S-System). "
                         "Improved implicit iterations are achieved by using a acceleration (recommended!).");
    addTypespecifcSubtags(VALUE_SERIAL_IMPLICIT, tag);
    tags.push_back(tag);
  }
  {
    XMLTag tag(*this, VALUE_PARALLEL_IMPLICIT, occ, TAG);
    tag.setDocumentation("Parallel Implicit coupling scheme according to block Jacobi iterations (V-System). "
                         "Improved implicit iterations are achieved by using a acceleration (recommended!).");
    addTypespecifcSubtags(VALUE_PARALLEL_IMPLICIT, tag);
    tags.push_back(tag);
  }
  {
    XMLTag tag(*this, VALUE_MULTI, occ, TAG);
    tag.setDocumentation("Multi coupling scheme according to block Jacobi iterations. "
                         "Improved implicit iterations are achieved by using a acceleration (recommended!).");
    addTypespecifcSubtags(VALUE_MULTI, tag);
    tags.push_back(tag);
  }

  for (XMLTag &tag : tags) {
    parent.addSubtag(tag);
  }
}

bool CouplingSchemeConfiguration::hasCouplingScheme(
    const std::string &participantName) const
{
  return utils::contained(participantName, _couplingSchemes);
}

const PtrCouplingScheme &CouplingSchemeConfiguration::getCouplingScheme(
    const std::string &participantName) const
{
  PRECICE_CHECK(utils::contained(participantName, _couplingSchemes),
                "No coupling scheme defined for "
                    << "participant \"" << participantName << "\". "
                    << "Please make sure to provide at least one <coupling-scheme:TYPE> in your "
                    << "precice-config.xml that couples this participant using the <participants .../> tag.");
  return _couplingSchemes.find(participantName)->second;
}

void CouplingSchemeConfiguration::xmlTagCallback(
    const xml::ConfigurationContext &context,
    xml::XMLTag &                    tag)
{
  PRECICE_TRACE(tag.getFullName());
  if (tag.getNamespace() == TAG) {
    _config.type = tag.getName();
    _accelerationConfig->clear();
  } else if (tag.getName() == TAG_PARTICIPANTS) {
    std::string first = tag.getStringAttributeValue(ATTR_FIRST);
    _config.participants.push_back(first);
    std::string second = tag.getStringAttributeValue(ATTR_SECOND);
    PRECICE_CHECK(std::find(_config.participants.begin(), _config.participants.end(), second) == _config.participants.end(),
                  "Provided first participant equals second participant in coupling scheme. Please correct the <participants "
                      << "first=\"" << first << "\" "
                      << "second=\"" << second << "\" "
                      << "/> tag in the <coupling-scheme:...> of your precice-config.xml");
    _config.participants.push_back(second);
  } else if (tag.getName() == TAG_PARTICIPANT) {
    PRECICE_ASSERT(_config.type == VALUE_MULTI);
    bool        control         = tag.getBooleanAttributeValue(ATTR_CONTROL);
    std::string participantName = tag.getStringAttributeValue(ATTR_NAME);
    PRECICE_CHECK(std::find(_config.participants.begin(), _config.participants.end(), participantName) == _config.participants.end() && participantName.compare(_config.controller) != 0,
                  "Participant \""
                      << participantName
                      << "\" is provided multiple times to multi coupling scheme. Please make sure that you do not provide the participant multiple times via the <participant name=\""
                      << participantName
                      << "\" /> tag in the <coupling-scheme:...> of your precice-config.xml");
    if (control) {
      PRECICE_CHECK(not _config.setController,
                    "Only one controller per MultiCouplingScheme can be defined. Please check the <participant "
                        << "name=\"" << participantName << "\" "
                        << "control=\"" << control << "\" "
                        << "/> tag in the <coupling-scheme:...> of your precice-config.xml");
      _config.controller    = participantName;
      _config.setController = true;
    } else {
      _config.participants.push_back(participantName);
    }

  } else if (tag.getName() == TAG_MAX_TIME) {
    _config.maxTime = tag.getDoubleAttributeValue(ATTR_VALUE);
    PRECICE_CHECK(_config.maxTime > 0, "Maximum time has to be larger than zero. Please check the <max-time "
                                           << "value=\"" << _config.maxTime << "\" "
                                           << "/> tag in the <coupling-scheme:...> of your precice-config.xml");
  } else if (tag.getName() == TAG_MAX_TIME_WINDOWS) {
    _config.maxTimeWindows = tag.getIntAttributeValue(ATTR_VALUE);
    PRECICE_CHECK(_config.maxTimeWindows > 0, "Maximum number of time windows has to be larger than zero. Please check the <max-time-windows "
                                                  << "value=\"" << _config.maxTimeWindows << "\" "
                                                  << "/> tag in the <coupling-scheme:...> of your precice-config.xml");
  } else if (tag.getName() == TAG_TIME_WINDOW_SIZE) {
    _config.timeWindowSize = tag.getDoubleAttributeValue(ATTR_VALUE);
    _config.validDigits    = tag.getIntAttributeValue(ATTR_VALID_DIGITS);
    PRECICE_CHECK((_config.validDigits >= 1) && (_config.validDigits < 17), "Valid digits of time window size has to be between 1 and 16.");
    _config.dtMethod = getTimesteppingMethod(tag.getStringAttributeValue(ATTR_METHOD));
    if (_config.dtMethod == constants::TimesteppingMethod::FIXED_TIME_WINDOW_SIZE) {
      PRECICE_CHECK(_config.timeWindowSize > 0, "Time window size has to be larger than zero. Please check the <time-window-size "
                                                    << "value=\"" << _config.timeWindowSize << "\" "
                                                    << "valid-digits=\"" << _config.validDigits << "\" "
                                                    << "method=\"" << tag.getStringAttributeValue(ATTR_METHOD) << "\" "
                                                    << "/> tag in the <coupling-scheme:...> of your precice-config.xml");
    } else {
      PRECICE_ASSERT(_config.dtMethod == constants::TimesteppingMethod::FIRST_PARTICIPANT_SETS_TIME_WINDOW_SIZE);
      PRECICE_CHECK(_config.timeWindowSize == -1, "Time window size value has to be equal to -1 (default), if method=\"first-participant\" is used. Please check the <time-window-size "
                                                      << "value=\"" << _config.timeWindowSize << "\" "
                                                      << "valid-digits=\"" << _config.validDigits << "\" "
                                                      << "method=\"" << tag.getStringAttributeValue(ATTR_METHOD) << "\" "
                                                      << "/> tag in the <coupling-scheme:...> of your precice-config.xml");
    }
    PRECICE_CHECK((_config.validDigits >= 1) && (_config.validDigits < 17), "Valid digits of time window size has to be between 1 and 16. Please check the <time-window-size "
                                                                                << "value=\"" << _config.timeWindowSize << "\" "
                                                                                << "valid-digits=\"" << _config.validDigits << "\" "
                                                                                << "method=\"" << tag.getStringAttributeValue(ATTR_METHOD) << "\" "
                                                                                << "/> tag in the <coupling-scheme:...> of your precice-config.xml");
  } else if (tag.getName() == TAG_ABS_CONV_MEASURE) {
    std::string dataName = tag.getStringAttributeValue(ATTR_DATA);
    std::string meshName = tag.getStringAttributeValue(ATTR_MESH);
    double      limit    = tag.getDoubleAttributeValue(ATTR_LIMIT);
    bool        suffices = tag.getBooleanAttributeValue(ATTR_SUFFICES);
    bool        strict   = tag.getBooleanAttributeValue(ATTR_STRICT);
    PRECICE_ASSERT(_config.type == VALUE_SERIAL_IMPLICIT || _config.type == VALUE_PARALLEL_IMPLICIT || _config.type == VALUE_MULTI);
    addAbsoluteConvergenceMeasure(dataName, meshName, limit, suffices, strict);
  } else if (tag.getName() == TAG_REL_CONV_MEASURE) {
    std::string dataName = tag.getStringAttributeValue(ATTR_DATA);
    std::string meshName = tag.getStringAttributeValue(ATTR_MESH);
    double      limit    = tag.getDoubleAttributeValue(ATTR_LIMIT);
    bool        suffices = tag.getBooleanAttributeValue(ATTR_SUFFICES);
    bool        strict   = tag.getBooleanAttributeValue(ATTR_STRICT);
    PRECICE_ASSERT(_config.type == VALUE_SERIAL_IMPLICIT || _config.type == VALUE_PARALLEL_IMPLICIT || _config.type == VALUE_MULTI);
    addRelativeConvergenceMeasure(dataName, meshName, limit, suffices, strict);
  } else if (tag.getName() == TAG_RES_REL_CONV_MEASURE) {
    std::string dataName = tag.getStringAttributeValue(ATTR_DATA);
    std::string meshName = tag.getStringAttributeValue(ATTR_MESH);
    double      limit    = tag.getDoubleAttributeValue(ATTR_LIMIT);
    bool        suffices = tag.getBooleanAttributeValue(ATTR_SUFFICES);
    bool        strict   = tag.getBooleanAttributeValue(ATTR_STRICT);
    PRECICE_ASSERT(_config.type == VALUE_SERIAL_IMPLICIT || _config.type == VALUE_PARALLEL_IMPLICIT || _config.type == VALUE_MULTI);
    addResidualRelativeConvergenceMeasure(dataName, meshName, limit, suffices, strict);
  } else if (tag.getName() == TAG_MIN_ITER_CONV_MEASURE) {
    std::string dataName      = tag.getStringAttributeValue(ATTR_DATA);
    std::string meshName      = tag.getStringAttributeValue(ATTR_MESH);
    int         minIterations = tag.getIntAttributeValue(ATTR_MIN_ITERATIONS);
    bool        suffices      = tag.getBooleanAttributeValue(ATTR_SUFFICES);
    bool        strict        = tag.getBooleanAttributeValue(ATTR_STRICT);
    PRECICE_ASSERT(_config.type == VALUE_SERIAL_IMPLICIT || _config.type == VALUE_PARALLEL_IMPLICIT || _config.type == VALUE_MULTI);
    addMinIterationConvergenceMeasure(dataName, meshName, minIterations, suffices, strict);
  } else if (tag.getName() == TAG_EXCHANGE) {
    std::string   nameData            = tag.getStringAttributeValue(ATTR_DATA);
    std::string   nameMesh            = tag.getStringAttributeValue(ATTR_MESH);
    std::string   nameParticipantFrom = tag.getStringAttributeValue(ATTR_FROM);
    std::string   nameParticipantTo   = tag.getStringAttributeValue(ATTR_TO);
    bool          initialize          = tag.getBooleanAttributeValue(ATTR_INITIALIZE);
    mesh::PtrData exchangeData;
    mesh::PtrMesh exchangeMesh;
    for (mesh::PtrMesh mesh : _meshConfig->meshes()) {
      if (mesh->getName() == nameMesh) {
        for (mesh::PtrData data : mesh->data()) {
          if (data->getName() == nameData) {
            exchangeData = data;
            exchangeMesh = mesh;
            break;
          }
        }
      }
    }
    PRECICE_CHECK(exchangeData.get(), "Mesh \"" << nameMesh << "\" with data \"" << nameData
                                                << "\" not defined. Please check the <exchange "
                                                << "data=\"" << nameData << "\" "
                                                << "mesh=\"" << nameMesh << "\" "
                                                << "from=\"" << nameParticipantFrom << "\" "
                                                << "to=\"" << nameParticipantTo << "\" "
                                                << "/> tag in the <coupling-scheme:... /> of your precice-config.xml.");
    _meshConfig->addNeededMesh(nameParticipantFrom, nameMesh);
    _meshConfig->addNeededMesh(nameParticipantTo, nameMesh);
    _config.exchanges.emplace_back(Config::Exchange{exchangeData, exchangeMesh, nameParticipantFrom, nameParticipantTo, initialize});
  } else if (tag.getName() == TAG_MAX_ITERATIONS) {
    PRECICE_ASSERT(_config.type == VALUE_SERIAL_IMPLICIT || _config.type == VALUE_PARALLEL_IMPLICIT || _config.type == VALUE_MULTI);
    _config.maxIterations = tag.getIntAttributeValue(ATTR_VALUE);
    PRECICE_CHECK(_config.maxIterations > 0,
                  "Maximal iteration limit has to be larger than zero. Please check the <max-iterations "
                      << "value = \"" << _config.maxIterations << "\" /> subtag in the <coupling-scheme:... /> of your precice-config.xml.");
  } else if (tag.getName() == TAG_EXTRAPOLATION) {
    PRECICE_ASSERT(_config.type == VALUE_SERIAL_IMPLICIT || _config.type == VALUE_PARALLEL_IMPLICIT || _config.type == VALUE_MULTI);
    _config.extrapolationOrder = tag.getIntAttributeValue(ATTR_VALUE);
    PRECICE_CHECK((_config.extrapolationOrder == 0) || (_config.extrapolationOrder == 1) || (_config.extrapolationOrder == 2),
                  "Extrapolation order has to be  0, 1, or 2. Please check the <extrapolation-order "
                      << "value=\"" << _config.extrapolationOrder << "\" "
                      << "/> subtag in the <coupling-scheme:... /> of your precice-config.xml.");
  }
}

void CouplingSchemeConfiguration::xmlEndTagCallback(
    const xml::ConfigurationContext &context,
    xml::XMLTag &                    tag)
{
  PRECICE_TRACE(tag.getFullName());
  if (tag.getNamespace() == TAG) {
    if (_config.type == VALUE_SERIAL_EXPLICIT) {
      std::string       accessor(_config.participants[0]);
      PtrCouplingScheme scheme = createSerialExplicitCouplingScheme(accessor);
      addCouplingScheme(scheme, accessor);
      //_couplingSchemes[accessor] = scheme;
      accessor = _config.participants[1];
      scheme   = createSerialExplicitCouplingScheme(accessor);
      addCouplingScheme(scheme, accessor);
      //_couplingSchemes[accessor] = scheme;
      _config = Config();
    } else if (_config.type == VALUE_PARALLEL_EXPLICIT) {
      std::string       accessor(_config.participants[0]);
      PtrCouplingScheme scheme = createParallelExplicitCouplingScheme(accessor);
      addCouplingScheme(scheme, accessor);
      //_couplingSchemes[accessor] = scheme;
      accessor = _config.participants[1];
      scheme   = createParallelExplicitCouplingScheme(accessor);
      addCouplingScheme(scheme, accessor);
      //_couplingSchemes[accessor] = scheme;
      _config = Config();
    } else if (_config.type == VALUE_SERIAL_IMPLICIT) {
      std::string       accessor(_config.participants[0]);
      PtrCouplingScheme scheme = createSerialImplicitCouplingScheme(accessor);
      addCouplingScheme(scheme, accessor);
      //_couplingSchemes[accessor] = scheme;
      accessor = _config.participants[1];
      scheme   = createSerialImplicitCouplingScheme(accessor);
      addCouplingScheme(scheme, accessor);
      //_couplingSchemes[accessor] = scheme;
      _config = Config();
    } else if (_config.type == VALUE_PARALLEL_IMPLICIT) {
      std::string       accessor(_config.participants[0]);
      PtrCouplingScheme scheme = createParallelImplicitCouplingScheme(accessor);
      addCouplingScheme(scheme, accessor);
      accessor = _config.participants[1];
      scheme   = createParallelImplicitCouplingScheme(accessor);
      addCouplingScheme(scheme, accessor);
      _config = Config();
    } else if (_config.type == VALUE_MULTI) {
      PRECICE_CHECK(_config.setController,
                    "One controller per MultiCoupling needs to be defined. "
                    "Please check the <participant name=... /> tags in the <coupling-scheme:... /> of your precice-config.xml. "
                    "Make sure that at least one participant tag provides the attribute <participant name=... control=\"True\"/>.");
      for (std::string &accessor : _config.participants) {
        PtrCouplingScheme scheme = createMultiCouplingScheme(accessor);
        addCouplingScheme(scheme, accessor);
      }
      PtrCouplingScheme scheme = createMultiCouplingScheme(_config.controller);
      addCouplingScheme(scheme, _config.controller);
      _config = Config();
    } else {
      PRECICE_ASSERT(false, _config.type);
    }
  }
}

void CouplingSchemeConfiguration::addCouplingScheme(
    PtrCouplingScheme  cplScheme,
    const std::string &participantName)
{
  PRECICE_TRACE(participantName);
  if (utils::contained(participantName, _couplingSchemes)) {
    PRECICE_DEBUG("Coupling scheme exists already for participant");
    if (utils::contained(participantName, _couplingSchemeCompositions)) {
      PRECICE_DEBUG("Coupling scheme composition exists already for participant");
      // Fetch the composition and add the new scheme.
      PRECICE_ASSERT(_couplingSchemeCompositions[participantName] != nullptr);
      _couplingSchemeCompositions[participantName]->addCouplingScheme(cplScheme);
    } else {
      PRECICE_DEBUG("No composition exists for the participant");
      // No composition exists, thus, the existing scheme is no composition.
      // Create a new composition, add the already existing and new scheme, and
      // overwrite the existing scheme with the composition.
      CompositionalCouplingScheme *composition = new CompositionalCouplingScheme();
      PRECICE_CHECK(nullptr == dynamic_cast<MultiCouplingScheme *>(_couplingSchemes[participantName].get()), "A Multi Coupling Scheme cannot yet be combined with any other coupling scheme. Try to include all participants within one multi coupling scheme instead.");
      composition->addCouplingScheme(_couplingSchemes[participantName]);
      composition->addCouplingScheme(cplScheme);
      _couplingSchemes[participantName] = PtrCouplingScheme(composition);
    }
  } else {
    PRECICE_DEBUG("No coupling scheme exists for the participant");
    // Store the new coupling scheme.
    _couplingSchemes[participantName] = cplScheme;
  }
}

void CouplingSchemeConfiguration::addTypespecifcSubtags(
    const std::string &type,
    //const std::string& name,
    xml::XMLTag &tag)
{
  PRECICE_TRACE(type);
  addTransientLimitTags(type, tag);
  _config.type = type;
  //_config.name = name;

  if (type == VALUE_SERIAL_EXPLICIT) {
    addTagParticipants(tag);
    addTagExchange(tag);
  } else if (type == VALUE_PARALLEL_EXPLICIT) {
    addTagParticipants(tag);
    addTagExchange(tag);
  } else if (type == VALUE_PARALLEL_IMPLICIT) {
    addTagParticipants(tag);
    addTagExchange(tag);
    addTagAcceleration(tag);
    addTagAbsoluteConvergenceMeasure(tag);
    addTagRelativeConvergenceMeasure(tag);
    addTagResidualRelativeConvergenceMeasure(tag);
    addTagMinIterationConvergenceMeasure(tag);
    addTagMaxIterations(tag);
    addTagExtrapolation(tag);
  } else if (type == VALUE_MULTI) {
    addTagParticipant(tag);
    addTagExchange(tag);
    addTagAcceleration(tag);
    addTagAbsoluteConvergenceMeasure(tag);
    addTagRelativeConvergenceMeasure(tag);
    addTagResidualRelativeConvergenceMeasure(tag);
    addTagMinIterationConvergenceMeasure(tag);
    addTagMaxIterations(tag);
    addTagExtrapolation(tag);
  } else if (type == VALUE_SERIAL_IMPLICIT) {
    addTagParticipants(tag);
    addTagExchange(tag);
    addTagAcceleration(tag);
    addTagAbsoluteConvergenceMeasure(tag);
    addTagRelativeConvergenceMeasure(tag);
    addTagResidualRelativeConvergenceMeasure(tag);
    addTagMinIterationConvergenceMeasure(tag);
    addTagMaxIterations(tag);
    addTagExtrapolation(tag);
  } else {
    // If wrong coupling scheme type is provided, this is already caught by the config parser. If the assertion below is triggered, it's a bug in preCICE, not wrong usage.
    PRECICE_ASSERT(false, "Unknown coupling scheme.");
  }
}

void CouplingSchemeConfiguration::addTransientLimitTags(
    const std::string &type,
    xml::XMLTag &      tag)
{
  using namespace xml;
  XMLTag tagMaxTime(*this, TAG_MAX_TIME, XMLTag::OCCUR_NOT_OR_ONCE);
  tagMaxTime.setDocumentation("Defined the end of the simulation as total time.");

  XMLAttribute<double> attrValueMaxTime(ATTR_VALUE);
  attrValueMaxTime.setDocumentation("The value of the maximum simulation time.");
  tagMaxTime.addAttribute(attrValueMaxTime);
  tag.addSubtag(tagMaxTime);

  XMLTag tagMaxTimeWindows(*this, TAG_MAX_TIME_WINDOWS, XMLTag::OCCUR_NOT_OR_ONCE);
  tagMaxTimeWindows.setDocumentation("Defined the end of the simulation as a total count of time windows.");
  XMLAttribute<int> attrValueMaxTimeWindows(ATTR_VALUE);
  attrValueMaxTimeWindows.setDocumentation("The maximum count of time windows.");
  tagMaxTimeWindows.addAttribute(attrValueMaxTimeWindows);
  tag.addSubtag(tagMaxTimeWindows);

  XMLTag tagTimeWindowSize(*this, TAG_TIME_WINDOW_SIZE, XMLTag::OCCUR_ONCE);
  tagTimeWindowSize.setDocumentation("Defines the size of the time window.");
  auto attrValueTimeWindowSize = makeXMLAttribute(ATTR_VALUE, CouplingScheme::UNDEFINED_TIME_WINDOW_SIZE)
                                     .setDocumentation("The maximum time window size.");
  tagTimeWindowSize.addAttribute(attrValueTimeWindowSize);
  XMLAttribute<int> attrValidDigits(ATTR_VALID_DIGITS, 10);
  attrValidDigits.setDocumentation(R"(Precision to use when checking for end of time windows used this many digits. \\(\phi = 10^{-validDigits}\\))");
  tagTimeWindowSize.addAttribute(attrValidDigits);
  std::vector<std::string> allowedMethods;
  if (type == VALUE_SERIAL_EXPLICIT || type == VALUE_SERIAL_IMPLICIT) {
    // method="first-participant" is only allowed for serial coupling schemes
    allowedMethods = {VALUE_FIXED, VALUE_FIRST_PARTICIPANT};
  } else {
    allowedMethods = {VALUE_FIXED};
  }
  auto attrMethod = makeXMLAttribute(ATTR_METHOD, VALUE_FIXED)
                        .setOptions(allowedMethods)
                        .setDocumentation("The method used to determine the time window size. Use `fixed` to fix the time window size for the participants.");
  tagTimeWindowSize.addAttribute(attrMethod);
  tag.addSubtag(tagTimeWindowSize);
}

void CouplingSchemeConfiguration::addTagParticipants(
    xml::XMLTag &tag)
{
  using namespace xml;
  XMLTag tagParticipants(*this, TAG_PARTICIPANTS, XMLTag::OCCUR_ONCE);
  tagParticipants.setDocumentation("Defines the participants of the coupling scheme.");
  XMLAttribute<std::string> attrFirst(ATTR_FIRST);
  attrFirst.setDocumentation("First participant to run the solver.");
  tagParticipants.addAttribute(attrFirst);
  XMLAttribute<std::string> attrSecond(ATTR_SECOND);
  attrSecond.setDocumentation("Second participant to run the solver.");
  tagParticipants.addAttribute(attrSecond);
  tag.addSubtag(tagParticipants);
}

void CouplingSchemeConfiguration::addTagParticipant(
    xml::XMLTag &tag)
{
  using namespace xml;
  XMLTag                    tagParticipant(*this, TAG_PARTICIPANT, XMLTag::OCCUR_ONCE_OR_MORE);
  XMLAttribute<std::string> attrName(ATTR_NAME);
  attrName.setDocumentation("Name of the participant.");
  tagParticipant.addAttribute(attrName);
  XMLAttribute<bool> attrControl(ATTR_CONTROL, false);
  attrControl.setDocumentation("Does this participant control the coupling?");
  tagParticipant.addAttribute(attrControl);
  tag.addSubtag(tagParticipant);
}

void CouplingSchemeConfiguration::addTagExchange(
    xml::XMLTag &tag)
{
  using namespace xml;
  XMLTag tagExchange(*this, TAG_EXCHANGE, XMLTag::OCCUR_ONCE_OR_MORE);
  tagExchange.setDocumentation("Defines the flow of data between meshes of participants.");

  auto attrData = XMLAttribute<std::string>(ATTR_DATA).setDocumentation("The data to exchange.");
  tagExchange.addAttribute(attrData);
  auto attrMesh = XMLAttribute<std::string>(ATTR_MESH).setDocumentation("The mesh which uses the data.");
  tagExchange.addAttribute(attrMesh);
  auto participantFrom = XMLAttribute<std::string>(ATTR_FROM).setDocumentation("The participant sending the data.");
  tagExchange.addAttribute(participantFrom);
  auto participantTo = XMLAttribute<std::string>(ATTR_TO).setDocumentation("The participant receiving the data.");
  tagExchange.addAttribute(participantTo);
  auto attrInitialize = XMLAttribute<bool>(ATTR_INITIALIZE, false).setDocumentation("Should this data be initialized during initializeData?");
  tagExchange.addAttribute(attrInitialize);
  tag.addSubtag(tagExchange);
}

void CouplingSchemeConfiguration::addTagAbsoluteConvergenceMeasure(
    xml::XMLTag &tag)
{
  using namespace xml;
  XMLTag tagConvergenceMeasure(*this, TAG_ABS_CONV_MEASURE, XMLTag::OCCUR_ARBITRARY);
  tagConvergenceMeasure.setDocumentation(
      "Absolute convergence criterion based on the two-norm difference of data values between iterations.\n"
      "\\$$\\left\\lVert H(x^k) - x^k \\right\\rVert_2 < \\text{limit}\\$$");
  addBaseAttributesTagConvergenceMeasure(tagConvergenceMeasure);
  XMLAttribute<double> attrLimit(ATTR_LIMIT);
  attrLimit.setDocumentation("Limit under which the measure is considered to have converged. Must be in \\((0, 1]\\).");
  tagConvergenceMeasure.addAttribute(attrLimit);
  tag.addSubtag(tagConvergenceMeasure);
}

void CouplingSchemeConfiguration::addTagResidualRelativeConvergenceMeasure(
    xml::XMLTag &tag)
{
  using namespace xml;
  XMLTag tagConvergenceMeasure(*this, TAG_RES_REL_CONV_MEASURE,
                               XMLTag::OCCUR_ARBITRARY);
  tagConvergenceMeasure.setDocumentation(
      "Residual relative convergence criterion based on the relative two-norm differences of data values between iterations.\n"
      "\\$$\\frac{\\left\\lVert H(x^k) - x^k \\right\\rVert_2}{\\left\\lVert H(x^{k-1}) - x^{k-1} \\right\\rVert_2} < \\text{limit}\\$$");
  addBaseAttributesTagConvergenceMeasure(tagConvergenceMeasure);
  XMLAttribute<double> attrLimit(ATTR_LIMIT);
  attrLimit.setDocumentation("Limit under which the measure is considered to have converged. Must be in \\((0, 1]\\).");
  tagConvergenceMeasure.addAttribute(attrLimit);
  tag.addSubtag(tagConvergenceMeasure);
}

void CouplingSchemeConfiguration::addTagRelativeConvergenceMeasure(
    xml::XMLTag &tag)
{
  using namespace xml;
  XMLTag tagConvergenceMeasure(*this, TAG_REL_CONV_MEASURE, XMLTag::OCCUR_ARBITRARY);
  tagConvergenceMeasure.setDocumentation(
      "Relative convergence criterion based on the relative two-norm difference of data values between iterations.\n"
      "\\$$\\frac{\\left\\lVert H(x^k) - x^k \\right\\rVert_2}{\\left\\lVert H(x^k) \\right\\rVert_2} < \\text{limit} \\$$");
  addBaseAttributesTagConvergenceMeasure(tagConvergenceMeasure);
  XMLAttribute<double> attrLimit(ATTR_LIMIT);
  attrLimit.setDocumentation(R"(Limit under which the measure is considered to have converged. Must be in \\((0, 1]\\).)");
  tagConvergenceMeasure.addAttribute(attrLimit);
  tag.addSubtag(tagConvergenceMeasure);
}

void CouplingSchemeConfiguration::addTagMinIterationConvergenceMeasure(
    xml::XMLTag &tag)
{
  xml::XMLTag tagMinIterationConvMeasure(*this,
                                         TAG_MIN_ITER_CONV_MEASURE, xml::XMLTag::OCCUR_ARBITRARY);
  tagMinIterationConvMeasure.setDocumentation("Convergence criterion used to ensure a miminimal amount of iterations. Specifying a mesh and data is required for technical reasons and does not influence the measure.");
  addBaseAttributesTagConvergenceMeasure(tagMinIterationConvMeasure);
  xml::XMLAttribute<int> attrMinIterations(ATTR_MIN_ITERATIONS);
  attrMinIterations.setDocumentation("The minimal amount of iterations.");
  tagMinIterationConvMeasure.addAttribute(attrMinIterations);
  tag.addSubtag(tagMinIterationConvMeasure);
}

void CouplingSchemeConfiguration::addBaseAttributesTagConvergenceMeasure(
    xml::XMLTag &tag)
{
  using namespace xml;
  auto attrData = XMLAttribute<std::string>(ATTR_DATA)
                      .setDocumentation("Data to be measured.");
  tag.addAttribute(attrData);
  auto attrMesh = XMLAttribute<std::string>(ATTR_MESH)
                      .setDocumentation("Mesh holding the data.");
  tag.addAttribute(attrMesh);
  auto attrSuffices = makeXMLAttribute(ATTR_SUFFICES, false)
                          .setDocumentation("If true, convergence of this measure is sufficient for overall convergence.");
  tag.addAttribute(attrSuffices);
  auto attrStrict = makeXMLAttribute(ATTR_STRICT, false)
                        .setDocumentation("If true, non-convergence of this measure ends the simulation. \"strict\" overrules \"suffices\".");
  tag.addAttribute(attrStrict);
}

void CouplingSchemeConfiguration::addTagMaxIterations(
    xml::XMLTag &tag)
{
  using namespace xml;
  XMLTag tagMaxIterations(*this, TAG_MAX_ITERATIONS, XMLTag::OCCUR_ONCE);
  tagMaxIterations.setDocumentation("Allows to specify a maximum amount of iterations per time window.");
  XMLAttribute<int> attrValue(ATTR_VALUE);
  attrValue.setDocumentation("The maximum value of iterations.");
  tagMaxIterations.addAttribute(attrValue);
  tag.addSubtag(tagMaxIterations);
}

void CouplingSchemeConfiguration::addTagExtrapolation(
    xml::XMLTag &tag)
{
  using namespace xml;
  XMLTag            tagExtrapolation(*this, TAG_EXTRAPOLATION, XMLTag::OCCUR_NOT_OR_ONCE);
  XMLAttribute<int> attrValue(ATTR_VALUE);
  attrValue.setDocumentation("The extrapolation order to use.");
  tagExtrapolation.addAttribute(attrValue);
  tagExtrapolation.setDocumentation("Sets order of predictor of interface values for first participant.");
  tag.addSubtag(tagExtrapolation);
}

void CouplingSchemeConfiguration::addTagAcceleration(
    xml::XMLTag &tag)
{
  PRECICE_TRACE(tag.getFullName());
  if (_accelerationConfig.get() == nullptr) {
    _accelerationConfig = acceleration::PtrAccelerationConfiguration(
        new acceleration::AccelerationConfiguration(_meshConfig));
  }
  _accelerationConfig->connectTags(tag);
}

void CouplingSchemeConfiguration::addAbsoluteConvergenceMeasure(
    const std::string &dataName,
    const std::string &meshName,
    double             limit,
    bool               suffices,
    bool               strict)
{
  PRECICE_TRACE();
  PRECICE_CHECK(math::greater(limit, 0.0),
                "Absolute convergence limit has to be greater than zero. Please check the "
                "<absolute-convergence-measure "
                    << "limit=\"" << limit << "\" "
                    << "data=\"" << dataName << "\" "
                    << "mesh=\"" << meshName << "\" "
                    << "/> subtag in your <coupling-scheme ... /> in the preCICE configuration file.");
  impl::PtrConvergenceMeasure measure(new impl::AbsoluteConvergenceMeasure(limit));
  ConvergenceMeasureDefintion convMeasureDef;
  convMeasureDef.data        = getData(dataName, meshName);
  convMeasureDef.suffices    = suffices;
  convMeasureDef.strict      = strict;
  convMeasureDef.meshName    = meshName;
  convMeasureDef.measure     = std::move(measure);
  convMeasureDef.doesLogging = true;
  _config.convergenceMeasureDefinitions.push_back(convMeasureDef);
}

void CouplingSchemeConfiguration::addRelativeConvergenceMeasure(
    const std::string &dataName,
    const std::string &meshName,
    double             limit,
    bool               suffices,
    bool               strict)
{
  PRECICE_TRACE();
  PRECICE_CHECK(math::greater(limit, 0.0) && math::greaterEquals(1.0, limit),
                "Relative convergence limit has to be in ]0;1]. Please check the "
                "<relative-convergence-measure "
                    << "limit=\"" << limit << "\" "
                    << "data=\"" << dataName << "\" "
                    << "mesh=\"" << meshName << "\" "
                    << "/> subtag in your <coupling-scheme ... /> in the preCICE configuration file.");
  impl::PtrConvergenceMeasure measure(new impl::RelativeConvergenceMeasure(limit));
  ConvergenceMeasureDefintion convMeasureDef;
  convMeasureDef.data        = getData(dataName, meshName);
  convMeasureDef.suffices    = suffices;
  convMeasureDef.strict      = strict;
  convMeasureDef.meshName    = meshName;
  convMeasureDef.measure     = std::move(measure);
  convMeasureDef.doesLogging = true;
  _config.convergenceMeasureDefinitions.push_back(convMeasureDef);
}

void CouplingSchemeConfiguration::addResidualRelativeConvergenceMeasure(
    const std::string &dataName,
    const std::string &meshName,
    double             limit,
    bool               suffices,
    bool               strict)
{
  PRECICE_TRACE();
  PRECICE_CHECK(math::greater(limit, 0.0) && math::greaterEquals(1.0, limit),
                "Relative convergence limit has to be in ]0;1]. Please check the "
                "<residual-relative-convergence-measure "
                    << "limit=\"" << limit << "\" "
                    << "data=\"" << dataName << "\" "
                    << "mesh=\"" << meshName << "\" "
                    << "/> subtag in your <coupling-scheme ... /> in the preCICE configuration file.");
  impl::PtrConvergenceMeasure measure(new impl::ResidualRelativeConvergenceMeasure(limit));
  ConvergenceMeasureDefintion convMeasureDef;
  convMeasureDef.data        = getData(dataName, meshName);
  convMeasureDef.suffices    = suffices;
  convMeasureDef.strict      = strict;
  convMeasureDef.meshName    = meshName;
  convMeasureDef.measure     = std::move(measure);
  convMeasureDef.doesLogging = true;
  _config.convergenceMeasureDefinitions.push_back(convMeasureDef);
}

void CouplingSchemeConfiguration::addMinIterationConvergenceMeasure(
    const std::string &dataName,
    const std::string &meshName,
    int                minIterations,
    bool               suffices,
    bool               strict)
{
  PRECICE_TRACE();
  impl::PtrConvergenceMeasure measure(new impl::MinIterationConvergenceMeasure(minIterations));
  ConvergenceMeasureDefintion convMeasureDef;
  convMeasureDef.data        = getData(dataName, meshName);
  convMeasureDef.suffices    = suffices;
  convMeasureDef.strict      = strict;
  convMeasureDef.meshName    = meshName;
  convMeasureDef.measure     = std::move(measure);
  convMeasureDef.doesLogging = false;
  _config.convergenceMeasureDefinitions.push_back(convMeasureDef);
}

mesh::PtrData CouplingSchemeConfiguration::getData(
    const std::string &dataName,
    const std::string &meshName) const
{
  for (mesh::PtrMesh mesh : _meshConfig->meshes()) {
    if (meshName == mesh->getName()) {
      for (mesh::PtrData data : mesh->data()) {
        if (dataName == data->getName()) {
          return data;
        }
      }
    }
  }
  PRECICE_ERROR("Data \"" << dataName << "\" used by mesh \""
                          << meshName << "\" is not configured.");
}

mesh::PtrData CouplingSchemeConfiguration::findDataByID(
    int ID) const
{
  for (mesh::PtrMesh mesh : _meshConfig->meshes()) {
    for (mesh::PtrData data : mesh->data()) {
      if (data->getID() == ID) {
        return data;
      }
    }
  }
  return nullptr;
}

PtrCouplingScheme CouplingSchemeConfiguration::createSerialExplicitCouplingScheme(
    const std::string &accessor) const
{
  PRECICE_TRACE(accessor);
  m2n::PtrM2N m2n = _m2nConfig->getM2N(
      _config.participants[0], _config.participants[1]);
  SerialCouplingScheme *scheme = new SerialCouplingScheme(
      _config.maxTime, _config.maxTimeWindows, _config.timeWindowSize,
      _config.validDigits, _config.participants[0], _config.participants[1],
      accessor, m2n, _config.dtMethod, BaseCouplingScheme::Explicit);

  addDataToBeExchanged(*scheme, accessor);

  return PtrCouplingScheme(scheme);
}

PtrCouplingScheme CouplingSchemeConfiguration::createParallelExplicitCouplingScheme(
    const std::string &accessor) const
{
  PRECICE_TRACE(accessor);
  m2n::PtrM2N m2n = _m2nConfig->getM2N(
      _config.participants[0], _config.participants[1]);
  ParallelCouplingScheme *scheme = new ParallelCouplingScheme(
      _config.maxTime, _config.maxTimeWindows, _config.timeWindowSize,
      _config.validDigits, _config.participants[0], _config.participants[1],
      accessor, m2n, _config.dtMethod, BaseCouplingScheme::Explicit);

  addDataToBeExchanged(*scheme, accessor);

  return PtrCouplingScheme(scheme);
}

PtrCouplingScheme CouplingSchemeConfiguration::createSerialImplicitCouplingScheme(
    const std::string &accessor) const
{
  PRECICE_TRACE(accessor);

  const auto first  = _config.participants[0];
  const auto second = _config.participants[1];

  m2n::PtrM2N m2n = _m2nConfig->getM2N(
      first, second);
  SerialCouplingScheme *scheme = new SerialCouplingScheme(
      _config.maxTime, _config.maxTimeWindows, _config.timeWindowSize,
      _config.validDigits, first, second,
      accessor, m2n, _config.dtMethod, BaseCouplingScheme::Implicit, _config.maxIterations);
  scheme->setExtrapolationOrder(_config.extrapolationOrder);

  addDataToBeExchanged(*scheme, accessor);
  PRECICE_CHECK(scheme->hasAnySendData(), "No send data configured. Use explicit scheme for one-way coupling. "
                                              << "Please check your <coupling-scheme ... /> and make sure that you provide at least one <exchange .../> subtag, where "
                                              << "from=\"" << accessor << "\".");

  // Add convergence measures
  PRECICE_CHECK(not _config.convergenceMeasureDefinitions.empty(),
                "At least one convergence measure has to be defined for an implicit coupling scheme. "
                "Please check your <coupling-scheme ... /> and make sure that you provide at least one <...-convergence-measure/> subtag in the precice-config.xml.");
  for (auto &elem : _config.convergenceMeasureDefinitions) {
    _meshConfig->addNeededMesh(second, elem.meshName);
    checkIfDataIsExchanged(elem.data->getID());
    scheme->addConvergenceMeasure(elem.data, elem.suffices, elem.strict, elem.measure, elem.doesLogging);
  }

  // Set relaxation parameters
  if (_accelerationConfig->getAcceleration().get() != nullptr) {
    for (std::string &neededMesh : _accelerationConfig->getNeededMeshes()) {
      _meshConfig->addNeededMesh(second, neededMesh);
    }
    for (const int dataID : _accelerationConfig->getAcceleration()->getDataIDs()) {
      checkSerialImplicitAccelerationData(dataID, first, second);
    }
    scheme->setAcceleration(_accelerationConfig->getAcceleration());
  }
  if (scheme->doesFirstStep() && _accelerationConfig->getAcceleration() && not _accelerationConfig->getAcceleration()->getDataIDs().empty()) {
    int dataID = *(_accelerationConfig->getAcceleration()->getDataIDs().begin());
    PRECICE_CHECK(not scheme->hasSendData(dataID),
                  "In case of serial coupling, acceleration can be defined for "
                      << "data of second participant only!");
  }

  return PtrCouplingScheme(scheme);
}

PtrCouplingScheme CouplingSchemeConfiguration::createParallelImplicitCouplingScheme(
    const std::string &accessor) const
{
  PRECICE_TRACE(accessor);
  m2n::PtrM2N m2n = _m2nConfig->getM2N(
      _config.participants[0], _config.participants[1]);
  ParallelCouplingScheme *scheme = new ParallelCouplingScheme(
      _config.maxTime, _config.maxTimeWindows, _config.timeWindowSize,
      _config.validDigits, _config.participants[0], _config.participants[1],
      accessor, m2n, _config.dtMethod, BaseCouplingScheme::Implicit, _config.maxIterations);
  scheme->setExtrapolationOrder(_config.extrapolationOrder);

  addDataToBeExchanged(*scheme, accessor);
  PRECICE_CHECK(scheme->hasAnySendData(), "No send data configured. Use explicit scheme for one-way coupling. "
                                              << "Please check your <coupling-scheme ... /> and make sure that you provide at least one <exchange .../> subtag, where "
                                              << "from=\"" << accessor << "\".");

  // Add convergence measures
  PRECICE_CHECK(not _config.convergenceMeasureDefinitions.empty(),
                "At least one convergence measure has to be defined for an implicit coupling scheme. "
                "Please check your <coupling-scheme ... /> and make sure that you provide at least one <...-convergence-measure/> subtag in the precice-config.xml.");
  for (auto &elem : _config.convergenceMeasureDefinitions) {
    _meshConfig->addNeededMesh(_config.participants[1], elem.meshName);
    checkIfDataIsExchanged(elem.data->getID());
    scheme->addConvergenceMeasure(elem.data, elem.suffices, elem.strict, elem.measure, elem.doesLogging);
  }

  // Set relaxation parameters
  if (_accelerationConfig->getAcceleration().get() != nullptr) {
    for (std::string &neededMesh : _accelerationConfig->getNeededMeshes()) {
      _meshConfig->addNeededMesh(_config.participants[1], neededMesh);
    }
    for (const int dataID : _accelerationConfig->getAcceleration()->getDataIDs()) {
      checkIfDataIsExchanged(dataID);
    }
    scheme->setAcceleration(_accelerationConfig->getAcceleration());
  }
  return PtrCouplingScheme(scheme);
}

PtrCouplingScheme CouplingSchemeConfiguration::createMultiCouplingScheme(
    const std::string &accessor) const
{
  PRECICE_TRACE(accessor);

  BaseCouplingScheme *scheme;

  if (accessor == _config.controller) {
    std::vector<m2n::PtrM2N> m2ns;
    for (const std::string &participant : _config.participants) {
      m2ns.push_back(_m2nConfig->getM2N(
          _config.controller, participant));
    }

    scheme = new MultiCouplingScheme(
        _config.maxTime, _config.maxTimeWindows, _config.timeWindowSize,
        _config.validDigits, accessor, m2ns, _config.dtMethod,
        _config.maxIterations);
    scheme->setExtrapolationOrder(_config.extrapolationOrder);

    MultiCouplingScheme *castedScheme = dynamic_cast<MultiCouplingScheme *>(scheme);
    PRECICE_ASSERT(castedScheme, "The dynamic cast of CouplingScheme failed.");
    addMultiDataToBeExchanged(*castedScheme, accessor);
  } else {
    m2n::PtrM2N m2n = _m2nConfig->getM2N(
        accessor, _config.controller);

    scheme = new ParallelCouplingScheme(
        _config.maxTime, _config.maxTimeWindows, _config.timeWindowSize,
        _config.validDigits, accessor, _config.controller,
        accessor, m2n, _config.dtMethod, BaseCouplingScheme::Implicit, _config.maxIterations);
    scheme->setExtrapolationOrder(_config.extrapolationOrder);

    BiCouplingScheme *castedScheme = dynamic_cast<BiCouplingScheme *>(scheme);
    addDataToBeExchanged(*castedScheme, accessor);
  }
  PRECICE_CHECK(scheme->hasAnySendData(), "No send data configured. Use explicit scheme for one-way coupling. "
                                              << "Please check your <coupling-scheme ... /> and make sure that you provide at least one <exchange .../> subtag, where "
                                              << "from=\"" << accessor << "\".");

  // Add convergence measures
  PRECICE_CHECK(not _config.convergenceMeasureDefinitions.empty(),
                "At least one convergence measure has to be defined for an implicit coupling scheme. "
                "Please check your <coupling-scheme ... /> and make sure that you provide at least one <...-convergence-measure/> subtag in the precice-config.xml.");
  for (auto &elem : _config.convergenceMeasureDefinitions) {
    _meshConfig->addNeededMesh(_config.controller, elem.meshName);
    checkIfDataIsExchanged(elem.data->getID());
    scheme->addConvergenceMeasure(elem.data, elem.suffices, elem.strict, elem.measure, elem.doesLogging);
  }

  // Set relaxation parameters
  if (_accelerationConfig->getAcceleration().get() != nullptr) {
    for (std::string &neededMesh : _accelerationConfig->getNeededMeshes()) {
      _meshConfig->addNeededMesh(_config.controller, neededMesh);
    }
    for (const int dataID : _accelerationConfig->getAcceleration()->getDataIDs()) {
      checkIfDataIsExchanged(dataID);
    }

    scheme->setAcceleration(_accelerationConfig->getAcceleration());
  }
  PRECICE_INFO(scheme->doesFirstStep());
  if (not scheme->doesFirstStep() && _accelerationConfig->getAcceleration()) {
    if (_accelerationConfig->getAcceleration()->getDataIDs().size() < 3) {
      PRECICE_WARN("Due to numerical reasons, for multi coupling, the number of coupling data vectors should be at least 3, not: "
                   << _accelerationConfig->getAcceleration()->getDataIDs().size() << ". "
                   << "Please check the <data .../> subtags in your <acceleration:.../> and make sure that you have at least 3.");
    }
  }

  return PtrCouplingScheme(scheme);
}

constants::TimesteppingMethod
CouplingSchemeConfiguration::getTimesteppingMethod(
    const std::string &method) const
{
  PRECICE_TRACE(method);
  if (method == VALUE_FIXED) {
    return constants::FIXED_TIME_WINDOW_SIZE;
  } else if (method == VALUE_FIRST_PARTICIPANT) {
    return constants::FIRST_PARTICIPANT_SETS_TIME_WINDOW_SIZE;
  } else {
    PRECICE_ASSERT(false, "Unknown timestepping method \"" << method << "\"!");
  }
}

void CouplingSchemeConfiguration::addDataToBeExchanged(
    BiCouplingScheme & scheme,
    const std::string &accessor) const
{
  PRECICE_TRACE();
  for (const Config::Exchange &exchange : _config.exchanges) {
    const std::string &from     = exchange.from;
    const std::string &to       = exchange.to;
    const std::string &dataName = exchange.data->getName();
    const std::string &meshName = exchange.mesh->getName();

    PRECICE_CHECK(to != from, "You cannot define an exchange from and to the same participant. "
                                  << "Please check the <exchange "
                                  << "data=\"" << dataName << "\" "
                                  << "mesh=\"" << meshName << "\" "
                                  << "from=\"" << from << "\" "
                                  << "to=\"" << to << "\" "
                                  << "/> tag in the <coupling-scheme:... /> of your precice-config.xml.");

    if (not(utils::contained(from, _config.participants) || from == _config.controller)) {
      PRECICE_CHECK(false, "Participant \"" + from + "\" is not configured for coupling scheme. "
                                                     "Please check the <exchange "
                               << "data=\"" << dataName << "\" "
                               << "mesh=\"" << meshName << "\" "
                               << "from=\"" << from << "\" "
                               << "to=\"" << to << "\" "
                               << "/> tag in the <coupling-scheme:... /> of your precice-config.xml.");
    }

    if (not(utils::contained(to, _config.participants) || to == _config.controller)) {
      PRECICE_CHECK(false, "Participant \"" + to + "\" is not configured for coupling scheme. "
                               << "Please check the <exchange "
                               << "data=\"" << dataName << "\" "
                               << "mesh=\"" << meshName << "\" "
                               << "from=\"" << from << "\" "
                               << "to=\"" << to << "\" "
                               << "/> tag in the <coupling-scheme:... /> of your precice-config.xml.");
    }

    const bool requiresInitialization = exchange.requiresInitialization;
    if (from == accessor) {
      scheme.addDataToSend(exchange.data, exchange.mesh, requiresInitialization);
      if (requiresInitialization && (_config.type == VALUE_SERIAL_EXPLICIT || _config.type == VALUE_SERIAL_IMPLICIT)) {
        PRECICE_CHECK(not scheme.doesFirstStep(), "In serial coupling only second participant can initialize data and send it. "
                                                      << "Please check the <exchange "
                                                      << "data=\"" << dataName << "\" "
                                                      << "mesh=\"" << meshName << "\" "
                                                      << "from=\"" << from << "\" "
                                                      << "to=\"" << to << "\" "
                                                      << "initialize=\"" << requiresInitialization << "\" "
                                                      << "/> tag in the <coupling-scheme:... /> of your precice-config.xml.");
      }
    } else if (to == accessor) {
      scheme.addDataToReceive(exchange.data, exchange.mesh, requiresInitialization);
      if (requiresInitialization && (_config.type == VALUE_SERIAL_EXPLICIT || _config.type == VALUE_SERIAL_IMPLICIT)) {
        PRECICE_CHECK(scheme.doesFirstStep(), "In serial coupling only first participant can receive initial data. "
                                                  << "Please check the <exchange "
                                                  << "data=\"" << dataName << "\" "
                                                  << "mesh=\"" << meshName << "\" "
                                                  << "from=\"" << from << "\" "
                                                  << "to=\"" << to << "\" "
                                                  << "initialize=\"" << requiresInitialization << "\" "
                                                  << "/> tag in the <coupling-scheme:... /> of your precice-config.xml.");
      }
    } else {
      PRECICE_ASSERT(_config.type == VALUE_MULTI);
    }
  }
}

void CouplingSchemeConfiguration::addMultiDataToBeExchanged(
    MultiCouplingScheme &scheme,
    const std::string &  accessor) const
{
  PRECICE_TRACE();
  for (const Config::Exchange &exchange : _config.exchanges) {
    const std::string &from = exchange.from;
    const std::string &to   = exchange.to;

    if (not(utils::contained(from, _config.participants) || from == _config.controller)) {
      throw std::runtime_error{"Participant \"" + from + "\" is not configured for coupling scheme"};
    }

    if (not(utils::contained(to, _config.participants) || to == _config.controller)) {
      throw std::runtime_error{"Participant \"" + to + "\" is not configured for coupling scheme"};
    }

    const bool initialize = exchange.requiresInitialization;
    if (from == accessor) {
      size_t index = 0;
      for (const std::string &participant : _config.participants) {
        PRECICE_DEBUG("from: " << from << ", to: " << to << ", participant: " << participant);
        if (to == participant) {
          break;
        }
        index++;
      }
      PRECICE_ASSERT(index < _config.participants.size(), index, _config.participants.size());
      scheme.addDataToSend(exchange.data, exchange.mesh, initialize, index);
    } else {
      size_t index = 0;
      for (const std::string &participant : _config.participants) {
        PRECICE_DEBUG("from: " << from << ", to: " << to << ", participant: " << participant);
        if (from == participant) {
          break;
        }
        index++;
      }
      PRECICE_ASSERT(index < _config.participants.size(), index, _config.participants.size());
      scheme.addDataToReceive(exchange.data, exchange.mesh, initialize, index);
    }
  }
}

void CouplingSchemeConfiguration::checkIfDataIsExchanged(
    int dataID) const
{
  const auto match = std::find_if(_config.exchanges.begin(),
                                  _config.exchanges.end(),
                                  [dataID](const Config::Exchange &exchange) { return exchange.data->getID() == dataID; });
  if (match != _config.exchanges.end()) {
    return;
  }

  // Data is not being exchanged
  std::string dataName = "";
  auto        dataptr  = findDataByID(dataID);
  if (dataptr) {
    dataName = dataptr->getName();
  }

  PRECICE_ERROR("You need to exchange every data that you use for convergence measures "
                << "and/or the iteration acceleration. Data \"" << dataName << "\" is "
                << "currently not exchanged over the respective mesh on which it is used for convergence measures and/or iteration "
                << "acceleration. Please check the <exchange ... /> and "
                << "<...-convergence-measure ... /> tags in the "
                << "<coupling-scheme:... /> of your precice-config.xml.");
}

void CouplingSchemeConfiguration::checkSerialImplicitAccelerationData(
    int                dataID,
    const std::string &first,
    const std::string &second) const
{
  checkIfDataIsExchanged(dataID);
  const auto match = std::find_if(_config.exchanges.begin(),
                                  _config.exchanges.end(),
                                  [dataID](const Config::Exchange &exchange) { return exchange.data->getID() == dataID; });
  PRECICE_ASSERT(match != _config.exchanges.end());
  const auto &exchange = *match;

  // In serial implicit cplschemes, data is only accelerated on the second participant.
  if (second == exchange.from) {
    return;
  }

  std::string dataName = "";
  auto        dataptr  = findDataByID(dataID);
  if (dataptr) {
    dataName = dataptr->getName();
  }

  // clang-format off
  PRECICE_ERROR(
      "You configured acceleration data \"" << dataName << "\" in the serial implicit coupling scheme between participants \"" << first << "\" and \"" << second << "\". "
      "For serial implicit coupling schemes, only data exchanged from the second to the first participant can be used for acceleration. Here, from \"" << second << "\" to \"" << first << "\". "
      "However, you configured data \"" << dataName << "\" for acceleration, which is exchanged from \"" << first << "\" to \"" << second << "\". "
      "Please remove this acceleration data tag or switch to a parallel implicit coupling scheme.");
  // clang-format on
}

} // namespace cplscheme
} // namespace precice
