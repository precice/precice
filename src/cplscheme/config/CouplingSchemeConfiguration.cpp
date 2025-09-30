#include <algorithm>
#include <cstddef>
#include <memory>
#include <ostream>
#include <stdexcept>
#include <utility>
#include <vector>

#include "CouplingSchemeConfiguration.hpp"
#include "acceleration/Acceleration.hpp"
#include "acceleration/AitkenAcceleration.hpp"
#include "acceleration/config/AccelerationConfiguration.hpp"
#include "cplscheme/BaseCouplingScheme.hpp"
#include "cplscheme/BiCouplingScheme.hpp"
#include "cplscheme/CompositionalCouplingScheme.hpp"
#include "cplscheme/CouplingScheme.hpp"
#include "cplscheme/MultiCouplingScheme.hpp"
#include "cplscheme/ParallelCouplingScheme.hpp"
#include "cplscheme/SerialCouplingScheme.hpp"
#include "cplscheme/SharedPointer.hpp"
#include "cplscheme/impl/AbsoluteConvergenceMeasure.hpp"
#include "cplscheme/impl/AbsoluteOrRelativeConvergenceMeasure.hpp"
#include "cplscheme/impl/RelativeConvergenceMeasure.hpp"
#include "cplscheme/impl/ResidualRelativeConvergenceMeasure.hpp"
#include "logging/LogMacros.hpp"
#include "m2n/SharedPointer.hpp"
#include "m2n/config/M2NConfiguration.hpp"
#include "mesh/Data.hpp"
#include "mesh/Mesh.hpp"
#include "mesh/config/MeshConfiguration.hpp"
#include "precice/config/ParticipantConfiguration.hpp"
#include "precice/impl/SharedPointer.hpp"
#include "precice/impl/Types.hpp"
#include "utils/Helpers.hpp"
#include "utils/assertion.hpp"
#include "xml/ConfigParser.hpp"
#include "xml/XMLAttribute.hpp"
#include "xml/XMLTag.hpp"

namespace precice::cplscheme {

const int CouplingSchemeConfiguration::DEFAULT_MIN_ITERATIONS(1);                                        // min 1 iteration
const int CouplingSchemeConfiguration::DEFAULT_MAX_ITERATIONS(CouplingScheme::UNDEFINED_MAX_ITERATIONS); // max infinite iterations

CouplingSchemeConfiguration::CouplingSchemeConfiguration(
    xml::XMLTag                         &parent,
    mesh::PtrMeshConfiguration           meshConfig,
    m2n::M2NConfiguration::SharedPointer m2nConfig,
    config::PtrParticipantConfiguration  participantConfig)
    : TAG("coupling-scheme"),
      TAG_PARTICIPANTS("participants"),
      TAG_PARTICIPANT("participant"),
      TAG_EXCHANGE("exchange"),
      TAG_MAX_TIME("max-time"),
      TAG_MAX_TIME_WINDOWS("max-time-windows"),
      TAG_TIME_WINDOW_SIZE("time-window-size"),
      TAG_ABS_CONV_MEASURE("absolute-convergence-measure"),
      TAG_ABS_OR_REL_CONV_MEASURE("absolute-or-relative-convergence-measure"),
      TAG_REL_CONV_MEASURE("relative-convergence-measure"),
      TAG_RES_REL_CONV_MEASURE("residual-relative-convergence-measure"),
      TAG_MIN_ITERATIONS("min-iterations"),
      TAG_MAX_ITERATIONS("max-iterations"),
      ATTR_DATA("data"),
      ATTR_MESH("mesh"),
      ATTR_PARTICIPANT("participant"),
      ATTR_INITIALIZE("initialize"),
      ATTR_EXCHANGE_SUBSTEPS("substeps"),
      ATTR_TYPE("type"),
      ATTR_FIRST("first"),
      ATTR_SECOND("second"),
      ATTR_VALUE("value"),
      ATTR_METHOD("method"),
      ATTR_LIMIT("limit"),
      ATTR_ABS_LIMIT("abs-limit"),
      ATTR_REL_LIMIT("rel-limit"),
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
      _meshConfig(std::move(meshConfig)),
      _m2nConfig(std::move(m2nConfig)),
      _participantConfig(participantConfig),
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
                "No coupling scheme defined for participant \"{}\". "
                "Please make sure to provide at least one <coupling-scheme:TYPE> in your "
                "precice-config.xml that couples this participant using the <participants .../> tag.",
                participantName);
  return _couplingSchemes.find(participantName)->second;
}

void CouplingSchemeConfiguration::xmlTagCallback(
    const xml::ConfigurationContext &context,
    xml::XMLTag                     &tag)
{
  PRECICE_TRACE(tag.getFullName());
  if (tag.getNamespace() == TAG) {
    _config.type = tag.getName();
    _accelerationConfig->clear();
  } else if (tag.getName() == TAG_PARTICIPANTS) {
    std::string first  = tag.getStringAttributeValue(ATTR_FIRST);
    std::string second = tag.getStringAttributeValue(ATTR_SECOND);

    PRECICE_CHECK(_participantConfig->hasParticipant(first),
                  "First participant in coupling-scheme <participants first=\"{}\" second=\"{}\" /> is unknown. {}",
                  first, second, _participantConfig->hintFor(first));
    PRECICE_CHECK(_participantConfig->hasParticipant(second),
                  "Second participant in coupling-scheme <participants first=\"{}\" second=\"{}\" /> is unknown. {}",
                  first, second, _participantConfig->hintFor(second));
    PRECICE_CHECK(first != second,
                  "First and second participant in coupling scheme are the same. "
                  "Please choose different in the <participants first=\"{}\" second=\"{}\" /> tag in the <coupling-scheme:...> of your precice-config.xml",
                  first, second);
    _config.participants.push_back(first);
    _config.participants.push_back(second);
  } else if (tag.getName() == TAG_PARTICIPANT) {
    PRECICE_ASSERT(_config.type == VALUE_MULTI);
    bool        control         = tag.getBooleanAttributeValue(ATTR_CONTROL);
    std::string participantName = tag.getStringAttributeValue(ATTR_NAME);
    PRECICE_CHECK(_participantConfig->hasParticipant(participantName),
                  "Provided participant in multi coupling-scheme <participant name=\"{}\" ... /> is unknown. {}",
                  participantName, _participantConfig->hintFor(participantName));
    PRECICE_CHECK(std::find(_config.participants.begin(), _config.participants.end(), participantName) == _config.participants.end() && participantName.compare(_config.controller) != 0,
                  "Participant \"{0}\" is provided multiple times to multi coupling scheme. "
                  "Please make sure that you do not provide the participant multiple times via the <participant name=\"{0}\" /> "
                  "tag in the <coupling-scheme:...> of your precice-config.xml",
                  participantName);
    if (control) {
      PRECICE_CHECK(not _config.setController,
                    "Only one controller per MultiCouplingScheme can be defined. "
                    "Please check the <participant name=\"{}\" control=\"{}\" /> tag in the <coupling-scheme:...> of your precice-config.xml",
                    participantName, control);
      _config.controller    = participantName;
      _config.setController = true;
    }
    _config.participants.push_back(participantName);
  } else if (tag.getName() == TAG_MAX_TIME) {
    _config.maxTime = tag.getDoubleAttributeValue(ATTR_VALUE);
    PRECICE_CHECK(_config.maxTime > 0,
                  "Maximum time has to be larger than zero. "
                  "Please check the <max-time value=\"{}\" /> tag in the <coupling-scheme:...> of your precice-config.xml",
                  _config.maxTime);
  } else if (tag.getName() == TAG_MAX_TIME_WINDOWS) {
    _config.maxTimeWindows = tag.getIntAttributeValue(ATTR_VALUE);
    PRECICE_CHECK(_config.maxTimeWindows > 0,
                  "Maximum number of time windows has to be larger than zero. "
                  "Please check the <max-time-windows value=\"{}\" /> tag in the <coupling-scheme:...> of your precice-config.xml",
                  _config.maxTimeWindows);
  } else if (tag.getName() == TAG_TIME_WINDOW_SIZE) {
    _config.timeWindowSize = tag.getDoubleAttributeValue(ATTR_VALUE);
    // Attribute does not exist for parallel coupling schemes as it is always fixed.
    _config.dtMethod = getTimesteppingMethod(tag.getStringAttributeValue(ATTR_METHOD, VALUE_FIXED));

    if (_config.dtMethod == constants::TimesteppingMethod::FIXED_TIME_WINDOW_SIZE) {
      PRECICE_CHECK(_config.timeWindowSize >= math::NUMERICAL_ZERO_DIFFERENCE,
                    "The minimal time window size supported by preCICE is {}. "
                    "Please check the <time-window-size value=\"{}\" /> tag "
                    "in the <coupling-scheme:...> of your precice-config.xml and pick an appropriate time window size.",
                    math::NUMERICAL_ZERO_DIFFERENCE, _config.timeWindowSize);
    } else {
      PRECICE_ASSERT(_config.dtMethod == constants::TimesteppingMethod::FIRST_PARTICIPANT_SETS_TIME_WINDOW_SIZE);
      PRECICE_WARN_IF(_config.timeWindowSize != CouplingScheme::UNDEFINED_TIME_WINDOW_SIZE,
                      "You combined a custom time-window-size of {} with method=\"first-participant\". "
                      "The given time-window-size will be ignored as it is prescribed by the participant.",
                      _config.timeWindowSize);
      _config.timeWindowSize = CouplingScheme::UNDEFINED_TIME_WINDOW_SIZE;
    }
  } else if (tag.getName() == TAG_ABS_CONV_MEASURE) {
    const std::string &dataName = tag.getStringAttributeValue(ATTR_DATA);
    const std::string &meshName = tag.getStringAttributeValue(ATTR_MESH);
    double             limit    = tag.getDoubleAttributeValue(ATTR_LIMIT);
    bool               suffices = tag.getBooleanAttributeValue(ATTR_SUFFICES);
    bool               strict   = tag.getBooleanAttributeValue(ATTR_STRICT);
    PRECICE_ASSERT(_config.type == VALUE_SERIAL_IMPLICIT || _config.type == VALUE_PARALLEL_IMPLICIT || _config.type == VALUE_MULTI);
    addAbsoluteConvergenceMeasure(dataName, meshName, limit, suffices, strict);
  } else if (tag.getName() == TAG_ABS_OR_REL_CONV_MEASURE) {
    const std::string &dataName = tag.getStringAttributeValue(ATTR_DATA);
    const std::string &meshName = tag.getStringAttributeValue(ATTR_MESH);
    double             absLimit = tag.getDoubleAttributeValue(ATTR_ABS_LIMIT);
    double             relLimit = tag.getDoubleAttributeValue(ATTR_REL_LIMIT);
    bool               suffices = tag.getBooleanAttributeValue(ATTR_SUFFICES);
    bool               strict   = tag.getBooleanAttributeValue(ATTR_STRICT);
    PRECICE_ASSERT(_config.type == VALUE_SERIAL_IMPLICIT || _config.type == VALUE_PARALLEL_IMPLICIT || _config.type == VALUE_MULTI);
    addAbsoluteOrRelativeConvergenceMeasure(dataName, meshName, absLimit, relLimit, suffices, strict);
  } else if (tag.getName() == TAG_REL_CONV_MEASURE) {
    const std::string &dataName = tag.getStringAttributeValue(ATTR_DATA);
    const std::string &meshName = tag.getStringAttributeValue(ATTR_MESH);
    double             limit    = tag.getDoubleAttributeValue(ATTR_LIMIT);
    bool               suffices = tag.getBooleanAttributeValue(ATTR_SUFFICES);
    bool               strict   = tag.getBooleanAttributeValue(ATTR_STRICT);
    PRECICE_ASSERT(_config.type == VALUE_SERIAL_IMPLICIT || _config.type == VALUE_PARALLEL_IMPLICIT || _config.type == VALUE_MULTI);
    addRelativeConvergenceMeasure(dataName, meshName, limit, suffices, strict);
  } else if (tag.getName() == TAG_RES_REL_CONV_MEASURE) {
    const std::string &dataName = tag.getStringAttributeValue(ATTR_DATA);
    const std::string &meshName = tag.getStringAttributeValue(ATTR_MESH);
    double             limit    = tag.getDoubleAttributeValue(ATTR_LIMIT);
    bool               suffices = tag.getBooleanAttributeValue(ATTR_SUFFICES);
    bool               strict   = tag.getBooleanAttributeValue(ATTR_STRICT);
    PRECICE_ASSERT(_config.type == VALUE_SERIAL_IMPLICIT || _config.type == VALUE_PARALLEL_IMPLICIT || _config.type == VALUE_MULTI);
    addResidualRelativeConvergenceMeasure(dataName, meshName, limit, suffices, strict);
  } else if (tag.getName() == TAG_EXCHANGE) {
    std::string nameData            = tag.getStringAttributeValue(ATTR_DATA);
    std::string nameMesh            = tag.getStringAttributeValue(ATTR_MESH);
    std::string nameParticipantFrom = tag.getStringAttributeValue(ATTR_FROM);
    std::string nameParticipantTo   = tag.getStringAttributeValue(ATTR_TO);
    bool        initialize          = tag.getBooleanAttributeValue(ATTR_INITIALIZE);
    bool        exchangeSubsteps    = tag.getBooleanAttributeValue(ATTR_EXCHANGE_SUBSTEPS);

    PRECICE_CHECK(_meshConfig->hasMeshName(nameMesh) && _meshConfig->getMesh(nameMesh)->hasDataName(nameData),
                  "Mesh \"{}\" with data \"{}\" not defined. "
                  "Please check the <exchange data=\"{}\" mesh=\"{}\" from=\"{}\" to=\"{}\" /> "
                  "tag in the <coupling-scheme:... /> of your precice-config.xml.",
                  nameMesh, nameData, nameData, nameMesh, nameParticipantFrom, nameParticipantTo);

    mesh::PtrMesh exchangeMesh = _meshConfig->getMesh(nameMesh);
    PRECICE_ASSERT(exchangeMesh);
    mesh::PtrData exchangeData = exchangeMesh->data(nameData);
    PRECICE_ASSERT(exchangeData);

    Config::Exchange newExchange{exchangeData, exchangeMesh, nameParticipantFrom, nameParticipantTo, initialize, exchangeSubsteps};
    PRECICE_CHECK(!_config.hasExchange(newExchange),
                  R"(Data "{}" of mesh "{}" cannot be exchanged multiple times between participants "{}" and "{}". Please remove one of the exchange tags.)",
                  nameData, nameMesh, nameParticipantFrom, nameParticipantTo);

    _meshConfig->addNeededMesh(nameParticipantFrom, nameMesh);
    _meshConfig->addNeededMesh(nameParticipantTo, nameMesh);
    _config.exchanges.emplace_back(std::move(newExchange));
  } else if (tag.getName() == TAG_MIN_ITERATIONS) {
    PRECICE_ASSERT(_config.type == VALUE_SERIAL_IMPLICIT || _config.type == VALUE_PARALLEL_IMPLICIT || _config.type == VALUE_MULTI);
    _config.minIterations = tag.getIntAttributeValue(ATTR_VALUE);
    PRECICE_CHECK(_config.minIterations > 0,
                  "Minimum iteration limit has to be larger than zero. Please check the <min-iterations value = \"{}\" /> subtag in the <coupling-scheme:... /> of your precice-config.xml.",
                  _config.minIterations);
  } else if (tag.getName() == TAG_MAX_ITERATIONS) {
    PRECICE_ASSERT(_config.type == VALUE_SERIAL_IMPLICIT || _config.type == VALUE_PARALLEL_IMPLICIT || _config.type == VALUE_MULTI);
    _config.maxIterations = tag.getIntAttributeValue(ATTR_VALUE);
    PRECICE_CHECK(_config.maxIterations > 0,
                  "Maximal iteration limit has to be larger than zero. "
                  "Please check the <max-iterations value=\"{0}\" /> subtag in the <coupling-scheme:... /> of your precice-config.xml. "
                  "To disable the iteration limit, remove the <max-iterations value=\"{0}\" /> subtag.",
                  _config.maxIterations);
  }

  // Additional consistency checks
  if (_config.minIterations > 0 && _config.maxIterations > 0) {
    PRECICE_CHECK(_config.minIterations <= _config.maxIterations,
                  "Maximum iteration limit {1} has to be larger or equal than the minimum iteration limit {0}. "
                  "Please check the <min-iterations value = \"{0}\" /> and <max-iterations value = \"{1}\" /> subtags in the <coupling-scheme:... /> of your precice-config.xml.",
                  _config.minIterations,
                  _config.maxIterations);
  }
}

void CouplingSchemeConfiguration::xmlEndTagCallback(
    const xml::ConfigurationContext &context,
    xml::XMLTag                     &tag)
{
  PRECICE_TRACE(tag.getFullName());
  if (tag.getNamespace() == TAG) {
    if (_config.type == VALUE_SERIAL_EXPLICIT) {
      PRECICE_CHECK(!_allowRemeshing, "Remeshing is currently incompatible with serial coupling schemes. Try using a parallel or a multi coupling scheme instead.");
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
      PRECICE_CHECK(!_allowRemeshing, "Remeshing is currently incompatible with serial coupling schemes. Try using a parallel or a multi coupling scheme instead.");
      updateConfigForImplicitCoupling();
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
      PRECICE_INFO_IF(_allowRemeshing, "Remeshing for implicit coupling schemes is in development. Currently, the acceleration data is deleted on remeshing.");
      updateConfigForImplicitCoupling();
      std::string       accessor(_config.participants[0]);
      PtrCouplingScheme scheme = createParallelImplicitCouplingScheme(accessor);
      addCouplingScheme(scheme, accessor);
      accessor = _config.participants[1];
      scheme   = createParallelImplicitCouplingScheme(accessor);
      addCouplingScheme(scheme, accessor);
      _config = Config();
    } else if (_config.type == VALUE_MULTI) {
      /// TODO test multi coupling scheme
      PRECICE_CHECK(!_allowRemeshing, "Remeshing is currently incompatible with multi coupling schemes. Try using a parallel coupling scheme instead.");
      PRECICE_INFO_IF(_allowRemeshing, "Remeshing for implicit coupling schemes is in development. Currently, the acceleration data is deleted on remeshing.");
      updateConfigForImplicitCoupling();
      PRECICE_CHECK(_config.setController,
                    "One controller per MultiCoupling needs to be defined. "
                    "Please check the <participant name=... /> tags in the <coupling-scheme:... /> of your precice-config.xml. "
                    "Make sure that at least one participant tag provides the attribute <participant name=... control=\"True\"/>.");
      for (const std::string &accessor : _config.participants) {
        PtrCouplingScheme scheme = createMultiCouplingScheme(accessor);
        addCouplingScheme(scheme, accessor);
      }
      _config = Config();
    } else {
      PRECICE_ASSERT(false, _config.type);
    }
  }
}

void CouplingSchemeConfiguration::addCouplingScheme(
    const PtrCouplingScheme &cplScheme,
    const std::string       &participantName)
{
  PRECICE_TRACE(participantName);
  if (!utils::contained(participantName, _couplingSchemes)) {
    PRECICE_DEBUG("No coupling scheme exists for the participant");
    // Store the new coupling scheme.
    _couplingSchemes[participantName] = cplScheme;
    return;
  }
  PRECICE_ASSERT(_couplingSchemes.count(participantName) > 0);

  PRECICE_CHECK(!_allowRemeshing, "Remeshing is currently incompatible with compositional coupling schemes. If you need remeshing, try using a multi coupling scheme to compose your participants.");

  // Create a composition to add the new cplScheme to
  if (!utils::contained(participantName, _couplingSchemeCompositions)) {
    PRECICE_DEBUG("Creating a compositional coupling scheme for the participant");
    auto composition = std::make_shared<CompositionalCouplingScheme>();
    composition->addCouplingScheme(_couplingSchemes[participantName]);
    _couplingSchemeCompositions[participantName] = composition.get();
    _couplingSchemes[participantName]            = std::move(composition);
  }

  PRECICE_ASSERT(_couplingSchemeCompositions.count(participantName) > 0);

  // Add the new scheme to the composition
  auto composition = _couplingSchemeCompositions.at(participantName);
  PRECICE_CHECK(!cplScheme->isImplicitCouplingScheme() || !composition->isImplicitCouplingScheme(),
                "You attempted to define a second implicit coupling-scheme for the participant \"{}\", which is not allowed. "
                "Please use a multi coupling-scheme for true implicit coupling of multiple participants.",
                participantName);
  _couplingSchemeCompositions[participantName]->addCouplingScheme(cplScheme);
}

void CouplingSchemeConfiguration::addTypespecifcSubtags(
    const std::string &type,
    xml::XMLTag       &tag)
{
  PRECICE_TRACE(type);
  addTransientLimitTags(type, tag);
  _config.type = type;
  //_config.name = name;

  if (type == VALUE_SERIAL_EXPLICIT) {
    addTagParticipants(tag);
    addTagExchange(tag, false);
  } else if (type == VALUE_PARALLEL_EXPLICIT) {
    addTagParticipants(tag);
    addTagExchange(tag, false);
  } else if (type == VALUE_PARALLEL_IMPLICIT) {
    addTagParticipants(tag);
    addTagExchange(tag, true);
    addTagAcceleration(tag);
    addTagAbsoluteConvergenceMeasure(tag);
    addTagAbsoluteOrRelativeConvergenceMeasure(tag);
    addTagRelativeConvergenceMeasure(tag);
    addTagResidualRelativeConvergenceMeasure(tag);
    addTagMinIterations(tag);
    addTagMaxIterations(tag);
  } else if (type == VALUE_MULTI) {
    addTagParticipant(tag);
    addTagExchange(tag, true);
    addTagAcceleration(tag);
    addTagAbsoluteConvergenceMeasure(tag);
    addTagAbsoluteOrRelativeConvergenceMeasure(tag);
    addTagRelativeConvergenceMeasure(tag);
    addTagResidualRelativeConvergenceMeasure(tag);
    addTagMinIterations(tag);
    addTagMaxIterations(tag);
  } else if (type == VALUE_SERIAL_IMPLICIT) {
    addTagParticipants(tag);
    addTagExchange(tag, true);
    addTagAcceleration(tag);
    addTagAbsoluteConvergenceMeasure(tag);
    addTagAbsoluteOrRelativeConvergenceMeasure(tag);
    addTagRelativeConvergenceMeasure(tag);
    addTagResidualRelativeConvergenceMeasure(tag);
    addTagMinIterations(tag);
    addTagMaxIterations(tag);
  } else {
    // If wrong coupling scheme type is provided, this is already caught by the config parser. If the assertion below is triggered, it's a bug in preCICE, not wrong usage.
    PRECICE_ASSERT(false, "Unknown coupling scheme.");
  }
}

void CouplingSchemeConfiguration::addTransientLimitTags(
    const std::string &type,
    xml::XMLTag       &tag)
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
  if (type == VALUE_SERIAL_EXPLICIT || type == VALUE_SERIAL_IMPLICIT) {
    // method="first-participant" is only allowed for serial coupling schemes
    auto attrMethod = makeXMLAttribute(ATTR_METHOD, VALUE_FIXED)
                          .setOptions({VALUE_FIXED, VALUE_FIRST_PARTICIPANT})
                          .setDocumentation("The method used to determine the time window size. Use `fixed` to fix the time window size for the participants.");
    tagTimeWindowSize.addAttribute(attrMethod);
  } else {
    tagTimeWindowSize.addAttributeHint(ATTR_METHOD, "This feature is only available for serial coupling schemes.");
  }
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
    xml::XMLTag &tag, bool substepsDefault)
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
  auto attrInitialize = XMLAttribute<bool>(ATTR_INITIALIZE, false).setDocumentation("Should this data be initialized during initialize?");
  tagExchange.addAttribute(attrInitialize);
  auto attrExchangeSubsteps = XMLAttribute<bool>(ATTR_EXCHANGE_SUBSTEPS, substepsDefault).setDocumentation("Should this data exchange substeps?");
  tagExchange.addAttribute(attrExchangeSubsteps);
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

void CouplingSchemeConfiguration::addTagAbsoluteOrRelativeConvergenceMeasure(
    xml::XMLTag &tag)
{
  using namespace xml;
  XMLTag tagConvergenceMeasure(*this, TAG_ABS_OR_REL_CONV_MEASURE, XMLTag::OCCUR_ARBITRARY);
  tagConvergenceMeasure.setDocumentation(
      "Absolute or relative convergence, which is the disjunction of an absolute criterion based on the two-norm difference of data values between iterations and a relative criterion based on the relative two-norm difference of data values between iterations,i.e. convergence is reached as soon as one of the both criteria is fulfilled. "
      "\\$$\\left\\lVert H(x^k) - x^k \\right\\rVert_2 < \\text{abs-limit}\\quad\\text{or}\\quad\\frac{\\left\\lVert H(x^k) - x^k \\right\\rVert_2}{\\left\\lVert H(x^k) \\right\\rVert_2} < \\text{rel-limit} \\$$  ");
  addBaseAttributesTagConvergenceMeasure(tagConvergenceMeasure);
  XMLAttribute<double> attrAbsLimit(ATTR_ABS_LIMIT);
  attrAbsLimit.setDocumentation(R"(Absolute limit under which the measure is considered to have converged.)");
  tagConvergenceMeasure.addAttribute(attrAbsLimit);
  XMLAttribute<double> attrRelLimit(ATTR_REL_LIMIT);
  attrAbsLimit.setDocumentation(R"(Relative limit under which the measure is considered to have converged. Must be in \\((0, 1]\\).)");
  tagConvergenceMeasure.addAttribute(attrRelLimit);
  tag.addSubtag(tagConvergenceMeasure);
}

void CouplingSchemeConfiguration::addTagResidualRelativeConvergenceMeasure(
    xml::XMLTag &tag)
{
  using namespace xml;
  XMLTag tagConvergenceMeasure(*this, TAG_RES_REL_CONV_MEASURE,
                               XMLTag::OCCUR_ARBITRARY);
  tagConvergenceMeasure.setDocumentation(
      "Relative convergence criterion comparing the currently measured residual to the residual of the first iteration in the time window.\n"
      "\\$$\\frac{\\left\\lVert H(x^k) - x^k \\right\\rVert_2}{\\left\\lVert H(x^0) - x^0 \\right\\rVert_2} < \\text{limit}\\$$");
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

void CouplingSchemeConfiguration::addTagMinIterations(
    xml::XMLTag &tag)
{
  using namespace xml;
  XMLTag tagMinIterations(*this, TAG_MIN_ITERATIONS, XMLTag::OCCUR_NOT_OR_ONCE);
  tagMinIterations.setDocumentation("Allows to specify a minimum amount of iterations that must be performed per time window.");
  XMLAttribute<int> attrValue(ATTR_VALUE);
  attrValue.setDocumentation("The minimum amount of iterations.");
  tagMinIterations.addAttribute(attrValue);
  tag.addSubtag(tagMinIterations);
}

void CouplingSchemeConfiguration::addTagMaxIterations(
    xml::XMLTag &tag)
{
  using namespace xml;
  XMLTag tagMaxIterations(*this, TAG_MAX_ITERATIONS, XMLTag::OCCUR_NOT_OR_ONCE);
  tagMaxIterations.setDocumentation("Allows to specify a maximum amount of iterations per time window.");
  XMLAttribute<int> attrValue(ATTR_VALUE);
  attrValue.setDocumentation("The maximum value of iterations.");
  tagMaxIterations.addAttribute(attrValue);
  tag.addSubtag(tagMaxIterations);
}

void CouplingSchemeConfiguration::addTagAcceleration(
    xml::XMLTag &tag)
{
  PRECICE_TRACE(tag.getFullName());
  if (!_accelerationConfig) {
    _accelerationConfig = std::make_shared<acceleration::AccelerationConfiguration>(
        _meshConfig);
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
                "Absolute convergence limit has to be greater than zero. "
                "Please check the <absolute-convergence-measure limit=\"{}\" data=\"{}\" mesh=\"{}\" /> subtag "
                "in your <coupling-scheme ... /> in the preCICE configuration file.",
                limit, dataName, meshName);
  ConvergenceMeasureDefintion convMeasureDef;
  convMeasureDef.data     = getData(dataName, meshName);
  convMeasureDef.suffices = suffices;
  convMeasureDef.strict   = strict;
  convMeasureDef.meshName = meshName;
  convMeasureDef.measure  = std::make_shared<impl::AbsoluteConvergenceMeasure>(limit);
  _config.convergenceMeasureDefinitions.push_back(convMeasureDef);
}

void CouplingSchemeConfiguration::addAbsoluteOrRelativeConvergenceMeasure(
    const std::string &dataName,
    const std::string &meshName,
    double             absLimit,
    double             relLimit,
    bool               suffices,
    bool               strict)
{
  PRECICE_TRACE();
  PRECICE_CHECK(math::greater(absLimit, 0.0),
                "Absolute convergence limit has to be greater than zero. "
                "Please check the <absolute-or-relative-convergence-measure abs-limit=\"{}\" rel-limit=\"{}\" data=\"{}\" mesh=\"{}\" /> subtag "
                "in your <coupling-scheme ... /> in the preCICE configuration file.",
                absLimit, relLimit, dataName, meshName);
  PRECICE_CHECK(math::greater(relLimit, 0.0) && math::greaterEquals(1.0, relLimit),
                "Relative convergence limit has to be in ]0;1]. "
                "Please check the <absolute-or-relative-convergence-measure abs-limit=\"{}\" rel-limit=\"{}\" data=\"{}\" mesh=\"{}\" /> subtag "
                "in your <coupling-scheme ... /> in the preCICE configuration file.",
                absLimit, relLimit, dataName, meshName);
  ConvergenceMeasureDefintion convMeasureDef;
  convMeasureDef.data     = getData(dataName, meshName);
  convMeasureDef.suffices = suffices;
  convMeasureDef.strict   = strict;
  convMeasureDef.meshName = meshName;
  convMeasureDef.measure  = std::make_shared<impl::AbsoluteOrRelativeConvergenceMeasure>(absLimit, relLimit);
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
                "Relative convergence limit has to be in ]0;1]. "
                "Please check the <relative-convergence-measure limit=\"{}\" data=\"{}\" mesh=\"{}\" /> subtag "
                "in your <coupling-scheme ... /> in the preCICE configuration file.",
                limit, dataName, meshName);
  PRECICE_WARN_IF(
      limit < 10 * math::NUMERICAL_ZERO_DIFFERENCE,
      "The relative convergence limit=\"{}\" is close to the hard-coded numerical resolution=\"{}\" of preCICE. "
      "This may lead to instabilities. The minimum relative convergence limit should be > \"{}\"  ",
      limit, math::NUMERICAL_ZERO_DIFFERENCE, 10 * math::NUMERICAL_ZERO_DIFFERENCE);

  ConvergenceMeasureDefintion convMeasureDef;
  convMeasureDef.data     = getData(dataName, meshName);
  convMeasureDef.suffices = suffices;
  convMeasureDef.strict   = strict;
  convMeasureDef.meshName = meshName;
  convMeasureDef.measure  = std::make_shared<impl::RelativeConvergenceMeasure>(limit);
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
                "Relative convergence limit has to be in ]0;1]. "
                "Please check the <residul-relative-convergence-measure limit=\"{}\" data=\"{}\" mesh=\"{}\" /> subtag "
                "in your <coupling-scheme ... /> in the preCICE configuration file.",
                limit, dataName, meshName);
  PRECICE_WARN_IF(
      limit < 10 * math::NUMERICAL_ZERO_DIFFERENCE,
      "The relative convergence limit=\"{}\" is close to the hard-coded numerical resolution=\"{}\" of preCICE. "
      "This may lead to instabilities. The minimum relative convergence limit should be > \"{}\"  ",
      limit, math::NUMERICAL_ZERO_DIFFERENCE, 10 * math::NUMERICAL_ZERO_DIFFERENCE);

  ConvergenceMeasureDefintion convMeasureDef;
  convMeasureDef.data     = getData(dataName, meshName);
  convMeasureDef.suffices = suffices;
  convMeasureDef.strict   = strict;
  convMeasureDef.meshName = meshName;
  convMeasureDef.measure  = std::make_shared<impl::ResidualRelativeConvergenceMeasure>(limit);
  _config.convergenceMeasureDefinitions.push_back(convMeasureDef);
}

mesh::PtrData CouplingSchemeConfiguration::getData(
    const std::string &dataName,
    const std::string &meshName) const
{
  PRECICE_CHECK(_meshConfig->hasMeshName(meshName) && _meshConfig->getMesh(meshName)->data(dataName),
                "Data \"{}\" used by mesh \"{}\" is not configured.", dataName, meshName);
  const mesh::PtrMesh &mesh = _meshConfig->getMesh(meshName);
  return mesh->data(dataName);
}

mesh::PtrData CouplingSchemeConfiguration::findDataByID(
    int ID) const
{
  for (const mesh::PtrMesh &mesh : _meshConfig->meshes()) {
    if (mesh->hasDataID(ID)) {
      return mesh->data(ID);
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
  SerialCouplingScheme *scheme = new SerialCouplingScheme(_config.maxTime, _config.maxTimeWindows, _config.timeWindowSize, _config.participants[0], _config.participants[1], accessor, m2n, _config.dtMethod, BaseCouplingScheme::Explicit);

  for (const auto &exchange : _config.exchanges) {
    if ((exchange.from == _config.participants[1]) && exchange.exchangeSubsteps) {
      PRECICE_WARN(
          "Exchange of substeps is activated in the serial-explicit coupling between the second participant \"{}\" and first participant \"{}\". This is inefficient as these substeps will never be used. You can turn this off in your preCICE configuration setting substeps=\"False\" in <exchange data=\"{}\" mesh=\"{}\" from=\"{}\" to=\"{}\" substeps=\"False\" />", exchange.from, exchange.to, exchange.data->getName(), exchange.mesh->getName(), exchange.from, exchange.to);
    }
  }

  addDataToBeExchanged(*scheme, accessor);

  return PtrCouplingScheme(scheme);
}

PtrCouplingScheme CouplingSchemeConfiguration::createParallelExplicitCouplingScheme(
    const std::string &accessor) const
{
  PRECICE_TRACE(accessor);
  PRECICE_ASSERT(_config.dtMethod == constants::TimesteppingMethod::FIXED_TIME_WINDOW_SIZE);
  m2n::PtrM2N m2n = _m2nConfig->getM2N(
      _config.participants[0], _config.participants[1]);
  ParallelCouplingScheme *scheme = new ParallelCouplingScheme(_config.maxTime, _config.maxTimeWindows, _config.timeWindowSize, _config.participants[0], _config.participants[1], accessor, m2n, BaseCouplingScheme::Explicit);

  for (const auto &exchange : _config.exchanges) {
    if (exchange.exchangeSubsteps) {
      PRECICE_WARN(
          "Exchange of substeps is activated in the parallel-explicit coupling between \"{}\" and \"{}\". This is inefficient as these substeps will never be used. You can turn this off in your preCICE configuration setting substeps=\"False\" in <exchange data=\"{}\" mesh=\"{}\" from=\"{}\" to=\"{}\" substeps=\"False\" />", exchange.from, exchange.to, exchange.data->getName(), exchange.mesh->getName(), exchange.from, exchange.to);
    }
  }

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
  SerialCouplingScheme *scheme = new SerialCouplingScheme(_config.maxTime, _config.maxTimeWindows, _config.timeWindowSize, first, second, accessor, m2n, _config.dtMethod, BaseCouplingScheme::Implicit, _config.minIterations, _config.maxIterations);

  addDataToBeExchanged(*scheme, accessor);
  PRECICE_CHECK(scheme->hasAnySendData(),
                "No send data configured. "
                "Use explicit scheme for one-way coupling. "
                "Please check your <coupling-scheme ... /> and make sure that you provide at least one <exchange .../> subtag, "
                "where from=\"{}\".",
                accessor);

  // Add convergence measures
  checkIterationLimits();
  addConvergenceMeasures(scheme, second, _config.convergenceMeasureDefinitions);

  // Set acceleration
  setSerialAcceleration(scheme, first, second);

  if (scheme->doesFirstStep() && _accelerationConfig->getAcceleration() && not _accelerationConfig->getAcceleration()->getPrimaryDataIDs().empty()) {
    DataID dataID = *(_accelerationConfig->getAcceleration()->getPrimaryDataIDs().begin());
    PRECICE_CHECK(not scheme->hasSendData(dataID),
                  "In case of serial coupling, acceleration can be defined for data of second participant only!");
  }

  return PtrCouplingScheme(scheme);
}

PtrCouplingScheme CouplingSchemeConfiguration::createParallelImplicitCouplingScheme(
    const std::string &accessor) const
{
  PRECICE_TRACE(accessor);
  PRECICE_ASSERT(_config.dtMethod == constants::TimesteppingMethod::FIXED_TIME_WINDOW_SIZE);
  m2n::PtrM2N m2n = _m2nConfig->getM2N(
      _config.participants[0], _config.participants[1]);
  ParallelCouplingScheme *scheme = new ParallelCouplingScheme(_config.maxTime, _config.maxTimeWindows, _config.timeWindowSize, _config.participants[0], _config.participants[1], accessor, m2n, BaseCouplingScheme::Implicit, _config.minIterations, _config.maxIterations);

  addDataToBeExchanged(*scheme, accessor);
  PRECICE_CHECK(scheme->hasAnySendData(),
                "No send data configured. Use explicit scheme for one-way coupling. "
                "Please check your <coupling-scheme ... /> and make sure that you provide at least one <exchange .../> subtag, "
                "where from=\"{}\".",
                accessor);

  // Add convergence measures
  checkIterationLimits();
  addConvergenceMeasures(scheme, _config.participants[1], _config.convergenceMeasureDefinitions);

  // Set acceleration
  setParallelAcceleration(scheme, _config.participants[1]);

  return PtrCouplingScheme(scheme);
}

PtrCouplingScheme CouplingSchemeConfiguration::createMultiCouplingScheme(
    const std::string &accessor) const
{
  PRECICE_TRACE(accessor);
  PRECICE_ASSERT(_config.dtMethod == constants::TimesteppingMethod::FIXED_TIME_WINDOW_SIZE);

  BaseCouplingScheme *scheme;

  std::map<std::string, m2n::PtrM2N> m2ns;
  for (const std::string &participant : _config.participants) {
    if (_m2nConfig->isM2NConfigured(accessor, participant)) {
      m2ns[participant] = _m2nConfig->getM2N(accessor, participant);
    }
  }

  scheme = new MultiCouplingScheme(
      _config.maxTime, _config.maxTimeWindows, _config.timeWindowSize, accessor, m2ns, _config.controller, _config.minIterations, _config.maxIterations);

  MultiCouplingScheme *castedScheme = dynamic_cast<MultiCouplingScheme *>(scheme);
  PRECICE_ASSERT(castedScheme, "The dynamic cast of CouplingScheme failed.");
  addMultiDataToBeExchanged(*castedScheme, accessor);

  PRECICE_CHECK(scheme->hasAnySendData(),
                "No send data configured. Use explicit scheme for one-way coupling. "
                "Please check your <coupling-scheme ... /> and make sure that you provide at least one "
                "<exchange .../> subtag, where from=\"{}\".",
                accessor);

  // Add convergence measures
  checkIterationLimits();
  if (accessor == _config.controller) {
    addConvergenceMeasures(scheme, _config.controller, _config.convergenceMeasureDefinitions);
  }

  // Set acceleration
  setParallelAcceleration(scheme, _config.controller);

  PRECICE_WARN_IF(
      not scheme->doesFirstStep() && _accelerationConfig->getAcceleration() && _accelerationConfig->getAcceleration()->getPrimaryDataIDs().size() < 3,
      "Due to numerical reasons, for multi coupling, the number of coupling data vectors should be at least 3, not: {}. "
      "Please check the <data .../> subtags in your <acceleration:.../> and make sure that you have at least 3.",
      _accelerationConfig->getAcceleration()->getPrimaryDataIDs().size());
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
    // We should never reach this point.
    PRECICE_UNREACHABLE("Unknown timestepping method '{}'.", method);
  }
}

void CouplingSchemeConfiguration::updateConfigForImplicitCoupling()
{
  if (_config.maxIterations == CouplingScheme::UNDEFINED_MAX_ITERATIONS) {
    _config.maxIterations = CouplingScheme::INFINITE_MAX_ITERATIONS;
  }
}

void CouplingSchemeConfiguration::checkIterationLimits() const
{
  if (_config.convergenceMeasureDefinitions.empty()) {
    PRECICE_CHECK(_config.maxIterations != -1,
                  "Not defining convergence measures without providing a maximum iteration limit is forbidden. "
                  "Please define a convergence measure or set a maximum iteration limit using <max-iterations value=\"...\" />.");

    PRECICE_INFO("No convergence measures were defined for an implicit coupling scheme. "
                 "It will always iterate the maximum amount iterations, which is {}. "
                 "You may want to add a convergence measure in your <coupling-scheme:.../> in your configuration.",
                 _config.maxIterations);
  }
}

void CouplingSchemeConfiguration::checkSubstepExchangeWaveformDegree(const Config::Exchange &exchange) const
{
  const auto &participant = _participantConfig->getParticipant(exchange.to);

  const auto &meshPtr = participant->findMesh(exchange.data->getName()); // related to https://github.com/precice/precice/issues/1694

  if (meshPtr == nullptr) {
    // Only warn, because might be valid configuration, if summation action is used. See Integration/Serial/SummationActionTwoSources.
    PRECICE_WARN("You defined <exchange data=\"{}\" ... to=\"{}\" /> in the <coupling-scheme:... />, but <participant name=\"{}\"> has no corresponding <read-data name=\"{}\" ... />. Usually this means that there is an error in your configuration.",
                 exchange.data->getName(), exchange.to, exchange.to, exchange.data->getName());
    return; // skip checks below
  }

  const auto &readDataContext = participant->readDataContext(meshPtr->getName(), exchange.data->getName());
  if (readDataContext.getWaveformDegree() == 0) {
    PRECICE_CHECK(!exchange.exchangeSubsteps,
                  "You configured <data:scalar/vector name=\"{}\" waveform-degree=\"{}\" />. Please deactivate exchange of substeps by setting substeps=\"false\" in the following exchange tag of your coupling scheme: <exchange data=\"{}\" mesh=\"{}\" from=\"{}\" to=\"{}\" />. Reason: For constant interpolation no exchange of data for substeps is needed. Please consider using waveform-degree=\"1\" or higher, if you want to use subcycling.",
                  readDataContext.getDataName(), readDataContext.getWaveformDegree(), exchange.data->getName(), exchange.mesh->getName(), exchange.from, exchange.to);
  } else if (readDataContext.getWaveformDegree() >= 2) {
    PRECICE_CHECK(exchange.exchangeSubsteps,
                  "You configured <data:scalar/vector name=\"{}\" waveform-degree=\"{}\" />. Please activate exchange of substeps by setting substeps=\"true\" in the following exchange tag of your coupling scheme: <exchange data=\"{}\" mesh=\"{}\" from=\"{}\" to=\"{}\" />. Reason: For higher-order interpolation exchange of data for substeps is required. If you don't want to activate exchange of additional data, please consider using waveform-degree=\"1\". Note that deactivating exchange of substep data might lead to worse results, if you use subcycling.",
                  readDataContext.getDataName(), readDataContext.getWaveformDegree(), exchange.data->getName(), exchange.mesh->getName(), exchange.from, exchange.to);
  } else { // For first degree there is no restriction for exchange of substeps
    PRECICE_ASSERT(readDataContext.getWaveformDegree() == 1);
  }
}

void CouplingSchemeConfiguration::addDataToBeExchanged(
    BiCouplingScheme  &scheme,
    const std::string &accessor) const
{
  PRECICE_TRACE();
  for (const Config::Exchange &exchange : _config.exchanges) {
    const std::string &from     = exchange.from;
    const std::string &to       = exchange.to;
    const std::string &dataName = exchange.data->getName();
    const std::string &meshName = exchange.mesh->getName();

    PRECICE_CHECK(to != from,
                  "You cannot define an exchange from and to the same participant. "
                  "Please check the <exchange data=\"{}\" mesh=\"{}\" from=\"{}\" to=\"{}\" /> tag in the <coupling-scheme:... /> of your precice-config.xml.",
                  dataName, meshName, from, to);

    PRECICE_CHECK((utils::contained(from, _config.participants) || from == _config.controller),
                  "Participant \"{}\" is not configured for coupling scheme. "
                  "Please check the <exchange data=\"{}\" mesh=\"{}\" from=\"{}\" to=\"{}\" /> tag in the <coupling-scheme:... /> of your precice-config.xml.",
                  from, dataName, meshName, from, to);

    PRECICE_CHECK((utils::contained(to, _config.participants) || to == _config.controller),
                  "Participant \"{}\" is not configured for coupling scheme. "
                  "Please check the <exchange data=\"{}\" mesh=\"{}\" from=\"{}\" to=\"{}\" /> tag in the <coupling-scheme:... /> of your precice-config.xml.",
                  to, dataName, meshName, from, to);

    const bool requiresInitialization = exchange.requiresInitialization;

    PRECICE_CHECK(
        !(requiresInitialization && _participantConfig->getParticipant(from)->isDirectAccessAllowed(exchange.mesh->getName())),
        "Participant \"{}\" cannot initialize data of the directly-accessed mesh \"{}\" from the participant\"{}\". "
        "Either disable the initialization in the <exchange /> tag or use a locally provided mesh instead.",
        from, meshName, to);

    const bool exchangeSubsteps = exchange.exchangeSubsteps;
    // Additional flag for direct access with write data:
    // For the  direct access is enabled
    const auto &fromParticipant  = _participantConfig->getParticipant(from);
    const bool  directAccessData = fromParticipant->isDirectAccessAllowed(exchange.mesh->getName()) &&
                                  fromParticipant->isDataWrite(exchange.mesh->getName(), exchange.data->getName());

    if (from == accessor) {
      scheme.addDataToSend(exchange.data, exchange.mesh, requiresInitialization, exchangeSubsteps, directAccessData);
    } else if (to == accessor) {
      checkSubstepExchangeWaveformDegree(exchange);
      scheme.addDataToReceive(exchange.data, exchange.mesh, requiresInitialization, exchangeSubsteps, directAccessData);
    } else {
      PRECICE_ASSERT(_config.type == VALUE_MULTI);
    }
  }
  scheme.determineInitialDataExchange();
}

void CouplingSchemeConfiguration::addMultiDataToBeExchanged(
    MultiCouplingScheme &scheme,
    const std::string   &accessor) const
{
  PRECICE_TRACE();
  for (const Config::Exchange &exchange : _config.exchanges) {
    const std::string &from     = exchange.from;
    const std::string &to       = exchange.to;
    const std::string &dataName = exchange.data->getName();
    const std::string &meshName = exchange.mesh->getName();

    PRECICE_CHECK(to != from,
                  "You cannot define an exchange from and to the same participant. "
                  "Please check the <exchange data=\"{}\" mesh=\"{}\" from=\"{}\" to=\"{}\" /> tag in the <coupling-scheme:... /> of your precice-config.xml.",
                  dataName, meshName, from, to);

    PRECICE_CHECK((utils::contained(from, _config.participants) || from == _config.controller),
                  "Participant \"{}\" is not configured for coupling scheme",
                  from);

    PRECICE_CHECK((utils::contained(to, _config.participants) || to == _config.controller),
                  "Participant \"{}\" is not configured for coupling scheme", to);

    const bool initialize       = exchange.requiresInitialization;
    const bool exchangeSubsteps = exchange.exchangeSubsteps;

    if (from == accessor) {
      scheme.addDataToSend(exchange.data, exchange.mesh, initialize, exchangeSubsteps, to);
    } else if (to == accessor) {
      scheme.addDataToReceive(exchange.data, exchange.mesh, initialize, exchangeSubsteps, from);
    }
  }
  scheme.determineInitialDataExchange();
}

void CouplingSchemeConfiguration::checkIfDataIsExchanged(
    DataID dataID, std::string_view participant) const
{
  const auto match = std::find_if(_config.exchanges.begin(),
                                  _config.exchanges.end(),
                                  [dataID, participant](const Config::Exchange &exchange) {
                                    // handle multi coupling
                                    if (exchange.from != participant && exchange.to != participant) {
                                      return false;
                                    } else {
                                      return exchange.data->getID() == dataID;
                                    }
                                  });
  if (match != _config.exchanges.end()) {
    return;
  }

  // Data is not being exchanged
  std::string dataName = "";
  auto        dataptr  = findDataByID(dataID);
  if (dataptr) {
    dataName = dataptr->getName();
  }

  PRECICE_ERROR("You need to exchange every data that you use for convergence measures and/or the iteration acceleration. "
                "Data \"{}\" is currently not exchanged over the respective mesh of participant \"{}\" on which it is used for convergence measures and/or iteration acceleration. "
                "Please check the <exchange ... /> and <...-convergence-measure ... /> tags in the <coupling-scheme:... /> of your precice-config.xml.",
                dataName, participant);
}

void CouplingSchemeConfiguration::checkSerialImplicitAccelerationData(
    int                dataID,
    const std::string &first,
    const std::string &second) const
{
  checkIfDataIsExchanged(dataID, second);
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

  PRECICE_ERROR(
      "You configured acceleration data \"{}\" in the serial implicit coupling scheme between participants \"{}\" and \"{}\". "
      "For serial implicit coupling schemes, only data exchanged from the second to the first participant can be used for acceleration. "
      "Here, from \"{}\" to \"{}\". "
      "However, you configured data \"{}\" for acceleration, which is exchanged from \"{}\" to \"{}\". "
      "Please remove this acceleration data tag or switch to a parallel implicit coupling scheme.",
      dataName, first, second, second, first, dataName, first, second);
}

void CouplingSchemeConfiguration::addConvergenceMeasures(
    BaseCouplingScheme                             *scheme,
    const std::string                              &participant,
    const std::vector<ConvergenceMeasureDefintion> &convergenceMeasureDefinitions) const
{
  for (auto &elem : convergenceMeasureDefinitions) {
    _meshConfig->addNeededMesh(participant, elem.meshName);
    checkIfDataIsExchanged(elem.data->getID(), participant);
    scheme->addConvergenceMeasure(elem.data->getID(), elem.suffices, elem.strict, elem.measure);
  }
}

void CouplingSchemeConfiguration::setSerialAcceleration(
    BaseCouplingScheme *scheme,
    const std::string  &first,
    const std::string  &second) const
{
  if (_accelerationConfig->getAcceleration().get() != nullptr) {
    for (std::string &neededMesh : _accelerationConfig->getNeededMeshes()) {
      _meshConfig->addNeededMesh(second, neededMesh);
    }
    for (const DataID dataID : _accelerationConfig->getAcceleration()->getPrimaryDataIDs()) {
      checkSerialImplicitAccelerationData(dataID, first, second);
    }
    scheme->setAcceleration(_accelerationConfig->getAcceleration());
  }
}

void CouplingSchemeConfiguration::setParallelAcceleration(
    BaseCouplingScheme *scheme,
    const std::string  &participant) const
{
  if (_accelerationConfig->getAcceleration().get() != nullptr) {
    for (std::string &neededMesh : _accelerationConfig->getNeededMeshes()) {
      _meshConfig->addNeededMesh(participant, neededMesh);
    }
    for (const DataID dataID : _accelerationConfig->getAcceleration()->getPrimaryDataIDs()) {
      checkIfDataIsExchanged(dataID, participant);
    }
    scheme->setAcceleration(_accelerationConfig->getAcceleration());

    PRECICE_WARN_IF(
        dynamic_cast<acceleration::AitkenAcceleration *>(_accelerationConfig->getAcceleration().get()) != nullptr,
        "You configured participant \"{}\" in a parallel-implicit coupling scheme with \"Aitken\" "
        "acceleration, which is known to perform bad in parallel coupling schemes. "
        "See https://precice.org/configuration-acceleration.html#dynamic-aitken-under-relaxation for details. "
        "Consider switching to a serial-implicit coupling scheme or changing the acceleration method.",
        participant);
  }
}

void CouplingSchemeConfiguration::setRemeshing(
    bool allowed)
{
  _allowRemeshing = allowed;
}

} // namespace precice::cplscheme
