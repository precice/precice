#
# Golden-output tests for the mocked preCICE library (extras/mock).
#
# The driver exercises the Participant API against libpreciceMocked and logs
# observable behavior as "T> " lines; check.sh compares them with the files in
# expected/. The expected outputs were validated against the real library with
# the differential harness (run.sh) — see extras/mock/tests/README.md.
#

add_executable(preciceMockedTestDriver ${CMAKE_CURRENT_LIST_DIR}/driver.cpp)
target_link_libraries(preciceMockedTestDriver PRIVATE preciceMocked)
set_target_properties(preciceMockedTestDriver PROPERTIES
  CXX_STANDARD ${PRECICE_CXX_STANDARD}
  CXX_STANDARD_REQUIRED Yes
  CXX_EXTENSIONS No
  )

set(PRECICE_MOCK_TEST_DIR ${CMAKE_CURRENT_LIST_DIR})

function(precice_add_mock_test name participant cfgdir)
  add_test(NAME precice.mock.${name}
    COMMAND bash ${PRECICE_MOCK_TEST_DIR}/check.sh
      $<TARGET_FILE:preciceMockedTestDriver>
      ${PRECICE_MOCK_TEST_DIR}
      ${name} ${participant} ${cfgdir} ${ARGN}
    )
  set_tests_properties(precice.mock.${name} PROPERTIES
    LABELS "mock"
    TIMEOUT 60
    )
endfunction()

# Lifecycle scenarios (control flow: windows, checkpoints, subcycling, max-time)
precice_add_mock_test(explicit-one       SolverOne explicit  lifecycle)
precice_add_mock_test(explicit-two       SolverTwo explicit  lifecycle)
precice_add_mock_test(explicit-sub-one   SolverOne explicit  lifecycle-sub)
precice_add_mock_test(explicit-sub-two   SolverTwo explicit  lifecycle-sub)
precice_add_mock_test(implicit-one       SolverOne implicit  lifecycle)
precice_add_mock_test(implicit-two       SolverTwo implicit  lifecycle)
precice_add_mock_test(implicit-sub-one   SolverOne implicit  lifecycle-sub)
precice_add_mock_test(implicit-sub-two   SolverTwo implicit  lifecycle-sub)
precice_add_mock_test(maxtime-one        SolverOne maxtime   lifecycle)
precice_add_mock_test(maxtime-two        SolverTwo maxtime   lifecycle)
precice_add_mock_test(firstpart-one      SolverOne firstpart lifecycle-fixed)
precice_add_mock_test(firstpart-two      SolverTwo firstpart lifecycle-fixed)

# Error scenarios (preconditions and messages)
precice_add_mock_test(err-ctor           SolverOne explicit errors ctor-badrank)
precice_add_mock_test(err-preinit-one    SolverOne explicit errors preinit)
precice_add_mock_test(err-preinit-two    SolverTwo explicit errors preinit)
precice_add_mock_test(err-emptymesh      SolverOne explicit errors emptymesh)
precice_add_mock_test(err-postinit-one   SolverOne explicit errors postinit)
precice_add_mock_test(err-postinit-two   SolverTwo explicit errors postinit)
precice_add_mock_test(err-afterloop-one  SolverOne explicit errors afterloop)
precice_add_mock_test(err-afterloop-two  SolverTwo explicit errors afterloop)
precice_add_mock_test(err-advzero        SolverOne explicit errors advzero)
precice_add_mock_test(err-advneg         SolverOne explicit errors advneg)
precice_add_mock_test(err-advinf         SolverOne explicit errors advinf)
precice_add_mock_test(err-advbig         SolverOne explicit errors advbig)
precice_add_mock_test(err-nocheckpoint   SolverOne implicit errors nocheckpoint)
