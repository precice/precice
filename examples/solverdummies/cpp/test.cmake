message("Running Two Solverdummies")
message("Using executable: ${DUMMY_EXE}")
message("Using config: ${DUMMY_CONFIG}")

execute_process(
  COMMAND ${DUMMY_EXE} ${DUMMY_CONFIG} SolverOne MeshOne
  COMMAND ${DUMMY_EXE} ${DUMMY_CONFIG} SolverTwo MeshTwo
)

message("Success!")
