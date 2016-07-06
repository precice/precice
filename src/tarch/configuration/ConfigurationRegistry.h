#ifndef _TARCH_CONFIGURATION_CONFIGURATION_REGISTRY_H_
#define _TARCH_CONFIGURATION_CONFIGURATION_REGISTRY_H_

#ifdef Parallel
#include <mpi.h>
#endif
#include "tarch/logging/Log.h"
#include "tarch/configuration/BaseTopLevelConfigurationFactory.h"
#include <string>
#include <map>
#include <list>
#include <vector>


namespace tarch {
  namespace configuration {
    class ConfigurationRegistry;
    class TopLevelConfiguration;
  }
  namespace irr {
    namespace io {
      class IXMLBase;

      template<class char_type, class super_class>
        class IIrrXMLReader;
    }
  }
}

/**
 * Register top level configuration classes.
 *
 * Use the macro registerTopLevelConfiguration in TopLevelConfiguration.cpp to register
 * top level configurations. It initializes a TopLevelConfigurationFactory,
 * which is added to this class. This happens before the main method of the application is called.
 *
 * Important: the method initTopLevelConfigurationFactories has to be called before reading a file.
 * This creates instances of the TopLevelConfiguration classes. This solves segmentation faults because
 * of uninitialized static variables.
 *
 * @author Tobias Weinzierl
 * @version $Revision: 1.21 $
 */
class tarch::configuration::ConfigurationRegistry {
  private:
    /**
     * Log device for the configuration component.
     */
    static tarch::logging::Log _log;

    typedef std::map<std::string,TopLevelConfiguration*> TopLevelConfigurationContainer;

    TopLevelConfigurationContainer _topLevelTags;

    std::vector<BaseTopLevelConfigurationFactory*> _topLevelTagFactories;

    /**
     * Standard constructor.
     */
    ConfigurationRegistry();
    
    /**
     * Called internally by readFile() and readString().
     */
    std::list<TopLevelConfiguration*> readConfiguration(
      tarch::irr::io::IIrrXMLReader<char, tarch::irr::io::IXMLBase>* xmlReader,
      const std::string& topLevelTag
    );
  public:
    virtual ~ConfigurationRegistry();

    static ConfigurationRegistry& getInstance();

    void addTopLevelConfiguration( TopLevelConfiguration* configuration );

    /**
     * Register a TopLevelConfiguration factory.
     */
    void addTopLevelConfigurationFactory(BaseTopLevelConfigurationFactory* configuration);

    /**
     * Initializes the registered TopLevelConfiguration factories. This method
     * should be called once in the main() method.
     */
    void initTopLevelConfigurationFactories();

    void enlistAvailableTopLevelTags() const;

    void writeDummyConfigFile(std::ostream& out) const;

    /**
     * Reads a configuration file and returns the runners. If the file is
     * invalid, the operation returns an empty set.
     *
     * @param filename Qualified name of the file to open.
     */
    std::list<TopLevelConfiguration*> readFile(
      const std::string& filename,
      const std::string& topLevelTag
    );
    
    /**
     * Reads a configuration string and returns the runners. If the string is
     * invalid, the operation returns an empty set.
     *
     * @param configString String containing the XML configuration data.
     */
    std::list<TopLevelConfiguration*> readString(
      const std::string& configString,
      const std::string& topLevelTag
    );

    void freeConfigurations(std::list<TopLevelConfiguration*>& configurations);
};


#endif

