#pragma once

#ifdef __cplusplus
extern"C" {
#endif

/**
 * @brief Name of action for writing iteration checkpoint.
 *
 * Fortran syntax:
 * precicef_action_write_iter_checkpoint( CHARACTER nameAction(*) )
 */
void precicef_action_write_iter_checkp_(
  char* nameAction,
  int   lengthNameAction );

/**
 * @brief Name of action for writing initial data.
 *
 * FortranSyntax:
 * precicef_action_write_initial_data( CHARACTER nameAction(*) )
 */
void precicef_action_write_initial_data_(
  char* nameAction,
  int   lengthNameAction );

/**
 * @brief Name of action for reading iteration checkpoint.
 *
 * Fortran syntax:
 * precicef_action_read_iter_checkpoint( CHARACTER nameAction(*) )
 */
void precicef_action_read_iter_checkp_(
  char* nameAction,
  int   lengthNameAction );

#ifdef __cplusplus
}
#endif

