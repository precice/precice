/*
 * BaseTopLevelConfigurationFactory.h
 *
 *  Created on: 20 mars 2013
 *      Author: tr
 */

#ifndef BASETOPLEVELCONFIGURATIONFACTORY_H_
#define BASETOPLEVELCONFIGURATIONFACTORY_H_

namespace tarch {
  namespace configuration {
    class BaseTopLevelConfigurationFactory;
  }
}

/**
 * This class is used to register the TopLevelConfiguration in the ConfigurationRegistry.
 *
 * It provides a init method, which is used to create an instance of the TopLevelConfiguration.
 * For further details see TopLevelConfigurationFactory and methods addTopLevelConfigurationFactory
 * and initTopLevelConfigurationFactory of ConfigurationRegistry.
 *
 * @author Thomas Rebele
 */
class tarch::configuration::BaseTopLevelConfigurationFactory {
  public:
	BaseTopLevelConfigurationFactory();
    virtual ~BaseTopLevelConfigurationFactory();

    virtual void init();
};



#endif /* BASETOPLEVELCONFIGURATIONFACTORY_H_ */
