#
# This file lists all integration test sources and test suites
#
target_sources(testprecice
    PRIVATE
    tests/fundamental/DifferentConfigs.cpp
    tests/fundamental/initial-data/InterleavedCreation.cpp
    tests/fundamental/initial-data/InterleavedCreationWithGradients.cpp
    tests/fundamental/profiling/NotAllStopped.cpp
    tests/fundamental/profiling/NotStoppedAtFinalize.cpp
    tests/fundamental/profiling/UserProfiling.cpp
    tests/geometric-multiscale/AxialGeoMultiscale.cpp
    tests/geometric-multiscale/RadialGeoMultiscale.cpp
    tests/parallel/CouplingOnLine.cpp
    tests/parallel/ExportTimeseries.cpp
    tests/parallel/GlobalRBFPartitioning.cpp
    tests/parallel/GlobalRBFPartitioningPETSc.cpp
    tests/parallel/LocalRBFPartitioning.cpp
    tests/parallel/LocalRBFPartitioningPETSc.cpp
    tests/parallel/MappingTypeRestriction.cpp
    tests/parallel/NearestProjectionRePartitioning.cpp
    tests/parallel/PrimaryRankSockets.cpp
    tests/parallel/TestBoundingBoxInitialization.cpp
    tests/parallel/TestBoundingBoxInitializationEmpty.cpp
    tests/parallel/TestBoundingBoxInitializationPUM.cpp
    tests/parallel/TestBoundingBoxInitializationTwoWay.cpp
    tests/parallel/TestFinalize.cpp
    tests/parallel/UserDefinedMPICommunicator.cpp
    tests/parallel/UserDefinedMPICommunicatorPetRBF.cpp
    tests/parallel/direct-mesh-access/AccessReceivedMeshAndMapping.cpp
    tests/parallel/direct-mesh-access/AccessReceivedMeshEmptyPartition.cpp
    tests/parallel/direct-mesh-access/AccessReceivedMeshEmptyPartitionTwoLevelInit.cpp
    tests/parallel/direct-mesh-access/AccessReceivedMeshNoOverlap.cpp
    tests/parallel/direct-mesh-access/AccessReceivedMeshNoOverlapTwoLevelInit.cpp
    tests/parallel/direct-mesh-access/AccessReceivedMeshOverlap.cpp
    tests/parallel/direct-mesh-access/AccessReceivedMeshOverlapNoWrite.cpp
    tests/parallel/direct-mesh-access/AccessReceivedMeshOverlapNoWriteTwoLevelInit.cpp
    tests/parallel/direct-mesh-access/AccessReceivedMeshOverlapTwoLevelInit.cpp
    tests/parallel/direct-mesh-access/helpers.cpp
    tests/parallel/direct-mesh-access/helpers.hpp
    tests/parallel/distributed-communication/TestDistributedCommunicationGatherScatterMPI.cpp
    tests/parallel/distributed-communication/TestDistributedCommunicationP2PMPI.cpp
    tests/parallel/distributed-communication/TestDistributedCommunicationP2PSockets.cpp
    tests/parallel/distributed-communication/helpers.cpp
    tests/parallel/distributed-communication/helpers.hpp
    tests/parallel/gather-scatter/EnforceGatherScatterEmptyPrimaryRank.cpp
    tests/parallel/gather-scatter/EnforceGatherScatterEmptyReceivedPrimaryRank.cpp
    tests/parallel/gather-scatter/helpers.cpp
    tests/parallel/gather-scatter/helpers.hpp
    tests/parallel/just-in-time-mapping/EmptyPartition.cpp
    tests/parallel/just-in-time-mapping/ExplicitRead.cpp
    tests/parallel/just-in-time-mapping/ExplicitWrite.cpp
    tests/parallel/just-in-time-mapping/ImplicitPUM.cpp
    tests/parallel/lifecycle/ConstructAndExplicitFinalize.cpp
    tests/parallel/lifecycle/ConstructOnly.cpp
    tests/parallel/lifecycle/Full.cpp
    tests/parallel/lifecycle/ImplicitFinalize.cpp
    tests/parallel/map-initial-data/NonZeroData.cpp
    tests/parallel/map-initial-data/ZeroData.cpp
    tests/parallel/map-initial-data/helper.hpp
    tests/parallel/mapping-nearest-neighbor-gradient/GradientTestParallelScalar.cpp
    tests/parallel/mapping-nearest-neighbor-gradient/GradientTestParallelVector.cpp
    tests/parallel/mapping-nearest-neighbor-gradient/GradientTestParallelWriteVector.cpp
    tests/parallel/mapping-volume/ParallelCube1To3.cpp
    tests/parallel/mapping-volume/ParallelCube3To1.cpp
    tests/parallel/mapping-volume/ParallelCubeConservative1To3.cpp
    tests/parallel/mapping-volume/ParallelCubeConservative3To1.cpp
    tests/parallel/mapping-volume/ParallelSquare1To2.cpp
    tests/parallel/mapping-volume/ParallelSquare2To1.cpp
    tests/parallel/mapping-volume/ParallelSquareConservative1To2.cpp
    tests/parallel/mapping-volume/ParallelTriangleConservative2To1.cpp
    tests/quasi-newton/helpers.cpp
    tests/quasi-newton/helpers.hpp
    tests/quasi-newton/parallel/TestQN1.cpp
    tests/quasi-newton/parallel/TestQN10.cpp
    tests/quasi-newton/parallel/TestQN10EmptyPartition.cpp
    tests/quasi-newton/parallel/TestQN11.cpp
    tests/quasi-newton/parallel/TestQN1EmptyPartition.cpp
    tests/quasi-newton/parallel/TestQN2.cpp
    tests/quasi-newton/parallel/TestQN2EmptyPartition.cpp
    tests/quasi-newton/parallel/TestQN3.cpp
    tests/quasi-newton/parallel/TestQN3EmptyPartition.cpp
    tests/quasi-newton/parallel/TestQN4.cpp
    tests/quasi-newton/parallel/TestQN4EmptyPartition.cpp
    tests/quasi-newton/parallel/TestQN5.cpp
    tests/quasi-newton/parallel/TestQN5EmptyPartition.cpp
    tests/quasi-newton/parallel/TestQN6.cpp
    tests/quasi-newton/parallel/TestQN6EmptyPartition.cpp
    tests/quasi-newton/parallel/TestQN7.cpp
    tests/quasi-newton/parallel/TestQN7EmptyPartition.cpp
    tests/quasi-newton/parallel/TestQN8.cpp
    tests/quasi-newton/parallel/TestQN8EmptyPartition.cpp
    tests/quasi-newton/parallel/TestQN9.cpp
    tests/quasi-newton/parallel/TestQN9EmptyPartition.cpp
    tests/quasi-newton/serial/DefaultConfig.cpp
    tests/quasi-newton/serial/TestQN1.cpp
    tests/quasi-newton/serial/TestQN10.cpp
    tests/quasi-newton/serial/TestQN11.cpp
    tests/quasi-newton/serial/TestQN12.cpp
    tests/quasi-newton/serial/TestQN13.cpp
    tests/quasi-newton/serial/TestQN14.cpp
    tests/quasi-newton/serial/TestQN2.cpp
    tests/quasi-newton/serial/TestQN3.cpp
    tests/quasi-newton/serial/TestQN4.cpp
    tests/quasi-newton/serial/TestQN5.cpp
    tests/quasi-newton/serial/TestQN6.cpp
    tests/quasi-newton/serial/TestQN7.cpp
    tests/quasi-newton/serial/TestQN8.cpp
    tests/quasi-newton/serial/TestQN9.cpp
    tests/remeshing/ReadAfterReset.cpp
    tests/remeshing/ResetAfterSubcycling.cpp
    tests/remeshing/ResetWhileIterating.cpp
    tests/remeshing/SubcylingAfterReset.cpp
    tests/remeshing/parallel-explicit/changed-mapping/RemeshBoth.cpp
    tests/remeshing/parallel-explicit/changed-mapping/RemeshBoth2LI.cpp
    tests/remeshing/parallel-explicit/changed-mapping/RemeshBothSerial.cpp
    tests/remeshing/parallel-explicit/changed-mapping/RemeshInput.cpp
    tests/remeshing/parallel-explicit/changed-mapping/RemeshInput2LI.cpp
    tests/remeshing/parallel-explicit/changed-mapping/RemeshInputSerial.cpp
    tests/remeshing/parallel-explicit/changed-mapping/RemeshOutput.cpp
    tests/remeshing/parallel-explicit/changed-mapping/RemeshOutput2LI.cpp
    tests/remeshing/parallel-explicit/changed-mapping/RemeshOutputSerial.cpp
    tests/remeshing/parallel-explicit/changed-partition/OverlapBoth.cpp
    tests/remeshing/parallel-explicit/changed-partition/OverlapBoth2LI.cpp
    tests/remeshing/parallel-explicit/changed-partition/ScatterOutputs.cpp
    tests/remeshing/parallel-explicit/changed-partition/ScatterOutputs2LI.cpp
    tests/remeshing/parallel-explicit/changed-partition/SwapOutputs.cpp
    tests/remeshing/parallel-explicit/changed-partition/SwapOutputs2LI.cpp
    tests/remeshing/parallel-explicit/helper.hpp
    tests/remeshing/parallel-explicit/noop/RemeshBoth.cpp
    tests/remeshing/parallel-explicit/noop/RemeshBoth2LI.cpp
    tests/remeshing/parallel-explicit/noop/RemeshBothSerial.cpp
    tests/remeshing/parallel-explicit/noop/RemeshInput.cpp
    tests/remeshing/parallel-explicit/noop/RemeshInput2LI.cpp
    tests/remeshing/parallel-explicit/noop/RemeshInputSerial.cpp
    tests/remeshing/parallel-explicit/noop/RemeshOutput.cpp
    tests/remeshing/parallel-explicit/noop/RemeshOutput2LI.cpp
    tests/remeshing/parallel-explicit/noop/RemeshOutputSerial.cpp
    tests/remeshing/parallel-explicit/two-meshes/LocalMapping.cpp
    tests/remeshing/parallel-explicit/two-meshes/RemoteMapping.cpp
    tests/remeshing/parallel-implicit/acceleration/IQN-ILS/RemeshBothSerial.cpp
    tests/remeshing/parallel-implicit/acceleration/IQN-ILS/RemeshFirstSerial.cpp
    tests/remeshing/parallel-implicit/acceleration/IQN-ILS/RemeshSecondSerial.cpp
    tests/remeshing/parallel-implicit/acceleration/IQN-IMVJ/RemeshBothSerial.cpp
    tests/remeshing/parallel-implicit/acceleration/IQN-IMVJ/RemeshFirstSerial.cpp
    tests/remeshing/parallel-implicit/acceleration/IQN-IMVJ/RemeshSecondSerial.cpp
    tests/remeshing/parallel-implicit/acceleration/aitken/RemeshBothSerial.cpp
    tests/remeshing/parallel-implicit/acceleration/aitken/RemeshFirstSerial.cpp
    tests/remeshing/parallel-implicit/acceleration/aitken/RemeshSecondSerial.cpp
    tests/remeshing/parallel-implicit/acceleration/constant/RemeshBothSerial.cpp
    tests/remeshing/parallel-implicit/acceleration/constant/RemeshFirstSerial.cpp
    tests/remeshing/parallel-implicit/acceleration/constant/RemeshSecondSerial.cpp
    tests/remeshing/parallel-implicit/acceleration/helper.hpp
    tests/remeshing/parallel-implicit/changed-mapping/RemeshBoth.cpp
    tests/remeshing/parallel-implicit/changed-mapping/RemeshBoth2LI.cpp
    tests/remeshing/parallel-implicit/changed-mapping/RemeshBothSerial.cpp
    tests/remeshing/parallel-implicit/changed-mapping/RemeshFirst.cpp
    tests/remeshing/parallel-implicit/changed-mapping/RemeshFirst2LI.cpp
    tests/remeshing/parallel-implicit/changed-mapping/RemeshFirstSerial.cpp
    tests/remeshing/parallel-implicit/changed-mapping/RemeshSecond.cpp
    tests/remeshing/parallel-implicit/changed-mapping/RemeshSecond2LI.cpp
    tests/remeshing/parallel-implicit/changed-mapping/RemeshSecondSerial.cpp
    tests/remeshing/parallel-implicit/changed-partition/OverlapBoth.cpp
    tests/remeshing/parallel-implicit/changed-partition/OverlapBoth2LI.cpp
    tests/remeshing/parallel-implicit/changed-partition/ScatterFirst.cpp
    tests/remeshing/parallel-implicit/changed-partition/ScatterFirst2LI.cpp
    tests/remeshing/parallel-implicit/changed-partition/ScatterSecond.cpp
    tests/remeshing/parallel-implicit/changed-partition/ScatterSecond2LI.cpp
    tests/remeshing/parallel-implicit/changed-partition/SwapFirst.cpp
    tests/remeshing/parallel-implicit/changed-partition/SwapFirst2LI.cpp
    tests/remeshing/parallel-implicit/changed-partition/SwapSecond.cpp
    tests/remeshing/parallel-implicit/changed-partition/SwapSecond2LI.cpp
    tests/remeshing/parallel-implicit/convergence/RemeshBothSerial.cpp
    tests/remeshing/parallel-implicit/convergence/RemeshFirstSerial.cpp
    tests/remeshing/parallel-implicit/convergence/RemeshSecondSerial.cpp
    tests/remeshing/parallel-implicit/helper.hpp
    tests/remeshing/parallel-implicit/noop/RemeshBoth.cpp
    tests/remeshing/parallel-implicit/noop/RemeshBoth2LI.cpp
    tests/remeshing/parallel-implicit/noop/RemeshBothSerial.cpp
    tests/remeshing/parallel-implicit/noop/RemeshFirst.cpp
    tests/remeshing/parallel-implicit/noop/RemeshFirst2LI.cpp
    tests/remeshing/parallel-implicit/noop/RemeshFirstSerial.cpp
    tests/remeshing/parallel-implicit/noop/RemeshSecond.cpp
    tests/remeshing/parallel-implicit/noop/RemeshSecond2LI.cpp
    tests/remeshing/parallel-implicit/noop/RemeshSecondSerial.cpp
    tests/serial/ImplicitCheckpointing.cpp
    tests/serial/PreconditionerBug.cpp
    tests/serial/SendMeshToMultipleParticipants.cpp
    tests/serial/SummationActionTwoSources.cpp
    tests/serial/TestExplicitWithDataMultipleReadWrite.cpp
    tests/serial/TestExplicitWithSolverGeometry.cpp
    tests/serial/TestImplicit.cpp
    tests/serial/TestReadAPI.cpp
    tests/serial/TwoAccelerations.cpp
    tests/serial/action-timings/ActionTimingsParallelExplicit.cpp
    tests/serial/action-timings/ActionTimingsParallelImplicit.cpp
    tests/serial/action-timings/ActionTimingsSerialExplicit.cpp
    tests/serial/action-timings/ActionTimingsSerialImplicit.cpp
    tests/serial/aitken/DefaultInitialRelaxation.cpp
    tests/serial/aitken/DefinedInitialRelaxation.cpp
    tests/serial/aitken/WithSubsteps.cpp
    tests/serial/circular/Explicit.cpp
    tests/serial/circular/helper.hpp
    tests/serial/compositional/FiveParticipants.cpp
    tests/serial/compositional/OneActivatedMuscle.cpp
    tests/serial/compositional/TwoActivatedMuscles.cpp
    tests/serial/compositional/TwoActivatedMusclesWithFeedback.cpp
    tests/serial/compositional/data/implicit-first/AllParallel.cpp
    tests/serial/compositional/data/implicit-first/Parallel.cpp
    tests/serial/compositional/data/implicit-first/SerialFirst.cpp
    tests/serial/compositional/data/implicit-first/SerialSecond.cpp
    tests/serial/compositional/data/implicit-second/AllParallel.cpp
    tests/serial/compositional/data/implicit-second/Parallel.cpp
    tests/serial/compositional/data/implicit-second/SerialFirst.cpp
    tests/serial/compositional/data/implicit-second/SerialSecond.cpp
    tests/serial/convergence-measures/helpers.cpp
    tests/serial/convergence-measures/helpers.hpp
    tests/serial/convergence-measures/testConvergenceMeasures1.cpp
    tests/serial/convergence-measures/testConvergenceMeasures2.cpp
    tests/serial/convergence-measures/testConvergenceMeasures3.cpp
    tests/serial/convergence-measures/testConvergenceMeasures4.cpp
    tests/serial/direct-mesh-access/DirectAccessReadWrite.cpp
    tests/serial/direct-mesh-access/DirectAccessWithDataInitialization.cpp
    tests/serial/direct-mesh-access/DirectAccessWithWaveform.cpp
    tests/serial/direct-mesh-access/Explicit.cpp
    tests/serial/direct-mesh-access/ExplicitAndMapping.cpp
    tests/serial/direct-mesh-access/ExplicitRead.cpp
    tests/serial/direct-mesh-access/ExplicitReadNoFlag.cpp
    tests/serial/direct-mesh-access/Implicit.cpp
    tests/serial/explicit/TestExplicitMPI.cpp
    tests/serial/explicit/TestExplicitMPISingle.cpp
    tests/serial/explicit/TestExplicitSockets.cpp
    tests/serial/explicit/helpers.cpp
    tests/serial/explicit/helpers.hpp
    tests/serial/initialize-data/Explicit.cpp
    tests/serial/initialize-data/Implicit.cpp
    tests/serial/initialize-data/ImplicitBoth.cpp
    tests/serial/initialize-data/ReadMapping.cpp
    tests/serial/initialize-data/WriteMapping.cpp
    tests/serial/initialize-data/helpers.cpp
    tests/serial/initialize-data/helpers.hpp
    tests/serial/just-in-time-mapping/ExplicitMultipleReadWrite.cpp
    tests/serial/just-in-time-mapping/ExplicitNoReadWrite.cpp
    tests/serial/just-in-time-mapping/ExplicitRead.cpp
    tests/serial/just-in-time-mapping/ExplicitReadPUM.cpp
    tests/serial/just-in-time-mapping/ExplicitWrite.cpp
    tests/serial/just-in-time-mapping/ExplicitWritePUM.cpp
    tests/serial/just-in-time-mapping/Implicit.cpp
    tests/serial/just-in-time-mapping/ImplicitDataInitialization.cpp
    tests/serial/just-in-time-mapping/ImplicitNoSubsteps.cpp
    tests/serial/just-in-time-mapping/ImplicitWithWaveform.cpp
    tests/serial/just-in-time-mapping/MappingWithDirectAccessRead.cpp
    tests/serial/just-in-time-mapping/MappingWithDirectAccessWrite.cpp
    tests/serial/just-in-time-mapping/SingleWrite.cpp
    tests/serial/just-in-time-mapping/UnsupportedConstraint.cpp
    tests/serial/just-in-time-mapping/UnsupportedMapping.cpp
    tests/serial/lifecycle/ConstructAndExplicitFinalize.cpp
    tests/serial/lifecycle/ConstructNoConfig.cpp
    tests/serial/lifecycle/ConstructNullComm.cpp
    tests/serial/lifecycle/ConstructOnly.cpp
    tests/serial/lifecycle/ConstructOnlyWait.cpp
    tests/serial/lifecycle/ConstructWrongCommSize.cpp
    tests/serial/lifecycle/Full.cpp
    tests/serial/lifecycle/FullWait.cpp
    tests/serial/lifecycle/ImplicitFinalize.cpp
    tests/serial/lifecycle/reconstruction/ConstructOnly.cpp
    tests/serial/lifecycle/reconstruction/Full.cpp
    tests/serial/lifecycle/reconstruction/ImplicitFinalize.cpp
    tests/serial/map-if-necessary/three-solvers/helper.hpp
    tests/serial/map-if-necessary/three-solvers/mixed-substeps/Multi.cpp
    tests/serial/map-if-necessary/three-solvers/mixed-substeps/ParallelExplicit.cpp
    tests/serial/map-if-necessary/three-solvers/mixed-substeps/ParallelImplicit.cpp
    tests/serial/map-if-necessary/three-solvers/mixed-substeps/SerialExplicit.cpp
    tests/serial/map-if-necessary/three-solvers/mixed-substeps/SerialImplicit.cpp
    tests/serial/map-if-necessary/three-solvers/with-substeps/Multi.cpp
    tests/serial/map-if-necessary/three-solvers/with-substeps/ParallelExplicit.cpp
    tests/serial/map-if-necessary/three-solvers/with-substeps/ParallelImplicit.cpp
    tests/serial/map-if-necessary/three-solvers/with-substeps/SerialExplicit.cpp
    tests/serial/map-if-necessary/three-solvers/with-substeps/SerialImplicit.cpp
    tests/serial/map-if-necessary/three-solvers/without-substeps/Multi.cpp
    tests/serial/map-if-necessary/three-solvers/without-substeps/ParallelExplicit.cpp
    tests/serial/map-if-necessary/three-solvers/without-substeps/ParallelImplicit.cpp
    tests/serial/map-if-necessary/three-solvers/without-substeps/SerialExplicit.cpp
    tests/serial/map-if-necessary/three-solvers/without-substeps/SerialImplicit.cpp
    tests/serial/map-if-necessary/two-solvers/helper.hpp
    tests/serial/map-if-necessary/two-solvers/mixed-substeps/ParallelExplicit.cpp
    tests/serial/map-if-necessary/two-solvers/mixed-substeps/ParallelImplicit.cpp
    tests/serial/map-if-necessary/two-solvers/mixed-substeps/SerialExplicit.cpp
    tests/serial/map-if-necessary/two-solvers/mixed-substeps/SerialImplicit.cpp
    tests/serial/map-if-necessary/two-solvers/with-substeps/ParallelExplicit.cpp
    tests/serial/map-if-necessary/two-solvers/with-substeps/ParallelImplicit.cpp
    tests/serial/map-if-necessary/two-solvers/with-substeps/SerialExplicit.cpp
    tests/serial/map-if-necessary/two-solvers/with-substeps/SerialImplicit.cpp
    tests/serial/map-if-necessary/two-solvers/without-substeps/ParallelExplicit.cpp
    tests/serial/map-if-necessary/two-solvers/without-substeps/ParallelImplicit.cpp
    tests/serial/map-if-necessary/two-solvers/without-substeps/SerialExplicit.cpp
    tests/serial/map-if-necessary/two-solvers/without-substeps/SerialImplicit.cpp
    tests/serial/map-initial-data/helper.hpp
    tests/serial/map-initial-data/non-zero-data/ParallelRead.cpp
    tests/serial/map-initial-data/non-zero-data/ParallelWrite.cpp
    tests/serial/map-initial-data/non-zero-data/SerialRead.cpp
    tests/serial/map-initial-data/non-zero-data/SerialWrite.cpp
    tests/serial/map-initial-data/zero-data/ParallelRead.cpp
    tests/serial/map-initial-data/zero-data/ParallelWrite.cpp
    tests/serial/map-initial-data/zero-data/SerialRead.cpp
    tests/serial/map-initial-data/zero-data/SerialWrite.cpp
    tests/serial/mapping-nearest-neighbor-gradient/GradientTestBidirectionalReadScalar.cpp
    tests/serial/mapping-nearest-neighbor-gradient/GradientTestBidirectionalReadVector.cpp
    tests/serial/mapping-nearest-neighbor-gradient/GradientTestBidirectionalWriteScalar.cpp
    tests/serial/mapping-nearest-neighbor-gradient/GradientTestBidirectionalWriteVector.cpp
    tests/serial/mapping-nearest-neighbor-gradient/GradientTestUnidirectionalReadBlockVector.cpp
    tests/serial/mapping-nearest-neighbor-gradient/GradientTestUnidirectionalReadScalar.cpp
    tests/serial/mapping-nearest-neighbor-gradient/GradientTestUnidirectionalReadVector.cpp
    tests/serial/mapping-nearest-neighbor-gradient/helpers.cpp
    tests/serial/mapping-nearest-neighbor-gradient/helpers.hpp
    tests/serial/mapping-nearest-projection/MappingNearestProjectionEdges.cpp
    tests/serial/mapping-nearest-projection/QuadMappingDiagonalNearestProjectionEdgesTallKite.cpp
    tests/serial/mapping-nearest-projection/QuadMappingDiagonalNearestProjectionEdgesWideKite.cpp
    tests/serial/mapping-nearest-projection/QuadMappingNearestProjectionEdges.cpp
    tests/serial/mapping-nearest-projection/helpers.cpp
    tests/serial/mapping-nearest-projection/helpers.hpp
    tests/serial/mapping-rbf-gaussian/GaussianShapeParameter.cpp
    tests/serial/mapping-rbf-gaussian/GaussianSupportRadius.cpp
    tests/serial/mapping-rbf-gaussian/helpers.cpp
    tests/serial/mapping-rbf-gaussian/helpers.hpp
    tests/serial/mapping-scaled-consistent/helpers.cpp
    tests/serial/mapping-scaled-consistent/helpers.hpp
    tests/serial/mapping-scaled-consistent/testQuadMappingScaledConsistentOnA.cpp
    tests/serial/mapping-scaled-consistent/testQuadMappingScaledConsistentOnB.cpp
    tests/serial/mapping-scaled-consistent/testTetraOnA.cpp
    tests/serial/mapping-scaled-consistent/testTetraOnB.cpp
    tests/serial/mapping-scaled-consistent/testVolumetricOnA2D.cpp
    tests/serial/mapping-scaled-consistent/testVolumetricOnB2D.cpp
    tests/serial/mapping-volume/OneTetraConservativeRead.cpp
    tests/serial/mapping-volume/OneTetraConservativeWrite.cpp
    tests/serial/mapping-volume/OneTetraRead.cpp
    tests/serial/mapping-volume/OneTetraWrite.cpp
    tests/serial/mapping-volume/OneTriangleConservativeRead.cpp
    tests/serial/mapping-volume/OneTriangleConservativeWrite.cpp
    tests/serial/mapping-volume/OneTriangleRead.cpp
    tests/serial/mapping-volume/OneTriangleWrite.cpp
    tests/serial/mapping-volume/helpers.cpp
    tests/serial/mapping-volume/helpers.hpp
    tests/serial/mesh-requirements/NearestNeighborA.cpp
    tests/serial/mesh-requirements/NearestNeighborB.cpp
    tests/serial/mesh-requirements/NearestProjection2DA.cpp
    tests/serial/mesh-requirements/NearestProjection2DB.cpp
    tests/serial/mixed-time-window-sizes/explicit/ParallelParallel.cpp
    tests/serial/mixed-time-window-sizes/explicit/ParallelSerial.cpp
    tests/serial/mixed-time-window-sizes/explicit/SerialParallel.cpp
    tests/serial/mixed-time-window-sizes/explicit/SerialSerial.cpp
    tests/serial/mixed-time-window-sizes/helper.cpp
    tests/serial/mixed-time-window-sizes/helper.hpp
    tests/serial/mixed-time-window-sizes/implicit/ParallelParallel.cpp
    tests/serial/mixed-time-window-sizes/implicit/ParallelSerial.cpp
    tests/serial/mixed-time-window-sizes/implicit/SerialParallel.cpp
    tests/serial/mixed-time-window-sizes/implicit/SerialSerial.cpp
    tests/serial/multi-coupling/MultiCoupling.cpp
    tests/serial/multi-coupling/MultiCouplingFourSolvers1.cpp
    tests/serial/multi-coupling/MultiCouplingFourSolvers2.cpp
    tests/serial/multi-coupling/MultiCouplingThreeSolvers1.cpp
    tests/serial/multi-coupling/MultiCouplingThreeSolvers10.cpp
    tests/serial/multi-coupling/MultiCouplingThreeSolvers2.cpp
    tests/serial/multi-coupling/MultiCouplingThreeSolvers3.cpp
    tests/serial/multi-coupling/MultiCouplingThreeSolvers4.cpp
    tests/serial/multi-coupling/MultiCouplingThreeSolvers5.cpp
    tests/serial/multi-coupling/MultiCouplingThreeSolvers6.cpp
    tests/serial/multi-coupling/MultiCouplingThreeSolvers7.cpp
    tests/serial/multi-coupling/MultiCouplingThreeSolvers8.cpp
    tests/serial/multi-coupling/MultiCouplingThreeSolvers9.cpp
    tests/serial/multi-coupling/MultiCouplingTwoSolvers1.cpp
    tests/serial/multi-coupling/MultiCouplingTwoSolvers2.cpp
    tests/serial/multi-coupling/helpers.cpp
    tests/serial/multi-coupling/helpers.hpp
    tests/serial/multiple-mappings/MultipleReadFromMappings.cpp
    tests/serial/multiple-mappings/MultipleReadToMappings.cpp
    tests/serial/multiple-mappings/MultipleWriteFromMappings.cpp
    tests/serial/multiple-mappings/MultipleWriteFromMappingsAndData.cpp
    tests/serial/multiple-mappings/MultipleWriteToMappings.cpp
    tests/serial/parallel-coupling/SolverAFirst.cpp
    tests/serial/parallel-coupling/SolverAFirstSubsteps.cpp
    tests/serial/parallel-coupling/SolverBFirst.cpp
    tests/serial/parallel-coupling/SolverBFirstSubsteps.cpp
    tests/serial/parallel-coupling/helpers.cpp
    tests/serial/parallel-coupling/helpers.hpp
    tests/serial/three-solvers/ThreeSolversExplicitExplicit.cpp
    tests/serial/three-solvers/ThreeSolversExplicitImplicit.cpp
    tests/serial/three-solvers/ThreeSolversFirstParticipant.cpp
    tests/serial/three-solvers/ThreeSolversImplicitExplicit.cpp
    tests/serial/three-solvers/ThreeSolversParallel.cpp
    tests/serial/three-solvers/ThreeSolversParallelExplicitImplicit.cpp
    tests/serial/three-solvers/helpers.cpp
    tests/serial/three-solvers/helpers.hpp
    tests/serial/time/explicit/compositional/DoNothingWithSubcycling.cpp
    tests/serial/time/explicit/compositional/ReadWriteScalarDataWithSubcycling.cpp
    tests/serial/time/explicit/compositional/ReadWriteScalarDataWithSubcyclingNoSubsteps.cpp
    tests/serial/time/explicit/parallel-coupling/DoManySmallSteps.cpp
    tests/serial/time/explicit/parallel-coupling/DoNonfittingWindows.cpp
    tests/serial/time/explicit/parallel-coupling/ReadWriteScalarDataWithSubcycling.cpp
    tests/serial/time/explicit/parallel-coupling/ReadWriteScalarDataWithSubcycling6400Steps.cpp
    tests/serial/time/explicit/parallel-coupling/ReadWriteScalarDataWithSubcycling640Steps.cpp
    tests/serial/time/explicit/parallel-coupling/ReadWriteScalarDataWithSubcycling64Steps.cpp
    tests/serial/time/explicit/parallel-coupling/ReadWriteScalarDataWithSubcyclingNoSubsteps.cpp
    tests/serial/time/explicit/parallel-coupling/helpers.cpp
    tests/serial/time/explicit/parallel-coupling/helpers.hpp
    tests/serial/time/explicit/serial-coupling/DoNothingWithSmallSteps.cpp
    tests/serial/time/explicit/serial-coupling/DoNothingWithSubcycling.cpp
    tests/serial/time/explicit/serial-coupling/ReadWriteScalarDataFirstParticipant.cpp
    tests/serial/time/explicit/serial-coupling/ReadWriteScalarDataFirstParticipantInitData.cpp
    tests/serial/time/explicit/serial-coupling/ReadWriteScalarDataWithSubcycling.cpp
    tests/serial/time/explicit/serial-coupling/ReadWriteScalarDataWithSubcyclingNoSubsteps.cpp
    tests/serial/time/explicit/serial-coupling/SecondTimeBug.cpp
    tests/serial/time/implicit/compositional/DoNothingWithSubcycling.cpp
    tests/serial/time/implicit/compositional/ReadWriteScalarDataWithSubcycling.cpp
    tests/serial/time/implicit/helpers.cpp
    tests/serial/time/implicit/helpers.hpp
    tests/serial/time/implicit/multi-coupling/DoNothingWithSubcycling.cpp
    tests/serial/time/implicit/multi-coupling/ReadWriteScalarDataWithSubcycling.cpp
    tests/serial/time/implicit/multi-coupling/ReadWriteScalarDataWithSubcyclingNoSubsteps.cpp
    tests/serial/time/implicit/multi-coupling/ReadWriteScalarDataWithWaveformSamplingFirst.cpp
    tests/serial/time/implicit/multi-coupling/ReadWriteScalarDataWithWaveformSubcyclingDifferentDts.cpp
    tests/serial/time/implicit/parallel-coupling/InitialDataAitken1.cpp
    tests/serial/time/implicit/parallel-coupling/InitialDataAitken2.cpp
    tests/serial/time/implicit/parallel-coupling/InitialDataConstantAcceleration1.cpp
    tests/serial/time/implicit/parallel-coupling/InitialDataConstantAcceleration2.cpp
    tests/serial/time/implicit/parallel-coupling/InitialDataIQNILSAcceleration1.cpp
    tests/serial/time/implicit/parallel-coupling/InitialDataIQNILSAcceleration2.cpp
    tests/serial/time/implicit/parallel-coupling/InitialDataIQNIMVJAcceleration1.cpp
    tests/serial/time/implicit/parallel-coupling/InitialDataIQNIMVJAcceleration2.cpp
    tests/serial/time/implicit/parallel-coupling/ReadWriteScalarDataWithSubcycling.cpp
    tests/serial/time/implicit/parallel-coupling/ReadWriteScalarDataWithSubcyclingNoSubsteps.cpp
    tests/serial/time/implicit/parallel-coupling/ReadWriteScalarDataWithSubcyclingNoSubstepsUseInitFirst.cpp
    tests/serial/time/implicit/parallel-coupling/ReadWriteScalarDataWithWaveformSamplingFirst.cpp
    tests/serial/time/implicit/parallel-coupling/ReadWriteScalarDataWithWaveformSamplingFirstNoInit.cpp
    tests/serial/time/implicit/parallel-coupling/ReadWriteScalarDataWithWaveformSubcyclingDifferentDts.cpp
    tests/serial/time/implicit/parallel-coupling/ReadWriteScalarDataWithWaveformSubcyclingFirst.cpp
    tests/serial/time/implicit/parallel-coupling/ReadWriteScalarDataWithWaveformSubcyclingMixed.cpp
    tests/serial/time/implicit/parallel-coupling/ReadWriteScalarDataWithWaveformSubcyclingThird.cpp
    tests/serial/time/implicit/parallel-coupling/WaveformSubcyclingWithDifferentNumberOfDts.cpp
    tests/serial/time/implicit/serial-coupling/InitialDataAitken1.cpp
    tests/serial/time/implicit/serial-coupling/InitialDataAitken2.cpp
    tests/serial/time/implicit/serial-coupling/InitialDataConstantAcceleration1.cpp
    tests/serial/time/implicit/serial-coupling/InitialDataConstantAcceleration2.cpp
    tests/serial/time/implicit/serial-coupling/InitialDataIQNILSAcceleration1.cpp
    tests/serial/time/implicit/serial-coupling/InitialDataIQNILSAcceleration2.cpp
    tests/serial/time/implicit/serial-coupling/InitialDataIQNIMVJAcceleration1.cpp
    tests/serial/time/implicit/serial-coupling/InitialDataIQNIMVJAcceleration2.cpp
    tests/serial/time/implicit/serial-coupling/ReadWriteScalarDataFirstParticipant.cpp
    tests/serial/time/implicit/serial-coupling/ReadWriteScalarDataWithSubcycling.cpp
    tests/serial/time/implicit/serial-coupling/ReadWriteScalarDataWithSubcyclingNoSubsteps.cpp
    tests/serial/time/implicit/serial-coupling/ReadWriteScalarDataWithSubcyclingNoSubstepsUseInitFirst.cpp
    tests/serial/time/implicit/serial-coupling/ReadWriteScalarDataWithWaveformSamplingFirst.cpp
    tests/serial/time/implicit/serial-coupling/ReadWriteScalarDataWithWaveformSamplingFirstNoInit.cpp
    tests/serial/time/implicit/serial-coupling/ReadWriteScalarDataWithWaveformSubcyclingFirst.cpp
    tests/serial/time/implicit/serial-coupling/ReadWriteScalarDataWithWaveformSubcyclingMixed.cpp
    tests/serial/time/implicit/serial-coupling/ReadWriteScalarDataWithWaveformSubcyclingThird.cpp
    tests/serial/time/implicit/serial-coupling/WaveformSubcyclingWithConstantAcceleration.cpp
    tests/serial/time/implicit/serial-coupling/WaveformSubcyclingWithConstantAccelerationNoInit.cpp
    tests/serial/watch/WatchIntegralScaleAndNoScaleParallel.cpp
    tests/serial/watch/WatchIntegralScaleAndNoScaleSerial.cpp
    tests/serial/watch/WatchPointParallel.cpp
    tests/serial/watch/WatchPointSerial.cpp
    tests/serial/watch/helpers.cpp
    tests/serial/watch/helpers.hpp
    tests/serial/whitebox/TestConfigurationComsol.cpp
    tests/serial/whitebox/TestConfigurationPeano.cpp
    tests/serial/whitebox/TestExplicitWithDataScaling.cpp
    )
