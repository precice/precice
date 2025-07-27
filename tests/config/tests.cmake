# precice_test_config_invalid( <PATH> <EXPECTED EXPRESSION> [Solver] [Ranks])
# precice_test_config_valid( <PATH> [Solver] [Ranks])

precice_test_config_invalid(missing.xml "unable to open configuration file")
precice_test_config_invalid(empty.xml "Tag <mesh> was not found")

precice_test_config_valid(solverdummies.xml)
precice_test_config_valid(solverdummies.xml SolverTwo)
precice_test_config_valid(solverdummies.xml SolverTwo 3)
precice_test_config_valid(solverdummies.xml SolverOne)
