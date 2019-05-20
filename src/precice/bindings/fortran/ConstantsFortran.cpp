extern"C" {
#include "ConstantsFortran.hpp"
}
#include "precice/SolverInterface.hpp"
#include "utils/assertion.hpp"

void precicef_action_write_iter_checkp_
(
  char*  nameAction,
  int lengthNameAction )
{
  const std::string& name = precice::constants::actionWriteIterationCheckpoint();
  assertion(name.size() < (size_t) lengthNameAction, name.size(), lengthNameAction);
  for (size_t i=0; i < name.size(); i++){
    nameAction[i] = name[i];
  }
}

void precicef_action_write_initial_data_(
  char*  nameAction,
  int lengthNameAction )
{
  const std::string& name = precice::constants::actionWriteInitialData();
  assertion(name.size() < (size_t) lengthNameAction, name.size(), lengthNameAction);
  for (size_t i=0; i < name.size(); i++){
    nameAction[i] = name[i];
  }
}

void precicef_action_read_iter_checkp_
(
  char*  nameAction,
  int lengthNameAction )
{
  const std::string& name = precice::constants::actionReadIterationCheckpoint();
  assertion(name.size() < (size_t) lengthNameAction, name.size(), lengthNameAction);
  for (size_t i=0; i < name.size(); i++){
    nameAction[i] = name[i];
  }
}

