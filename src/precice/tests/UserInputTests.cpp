#include "precice/Participant.hpp"
#include "testing/Testing.hpp"

#include <array>
#include <limits>

BOOST_AUTO_TEST_SUITE(PreciceTests)
BOOST_AUTO_TEST_SUITE(UserInput)

BOOST_AUTO_TEST_SUITE(Constructor)

PRECICE_TEST_SETUP(1_rank)
BOOST_AUTO_TEST_CASE(NoSolverName)
{
  PRECICE_TEST();

  auto config = precice::testing::getPathToSources() + "/precice/tests/config-checker.xml";
  BOOST_CHECK_EXCEPTION(precice::Participant("", config, context.rank, context.size),
                        ::precice::Error,
                        precice::testing::errorContains("name is an empty string"));
}

PRECICE_TEST_SETUP(1_rank)
BOOST_AUTO_TEST_CASE(WrongSolverName)
{
  PRECICE_TEST();

  auto config = precice::testing::getPathToSources() + "/precice/tests/config-checker.xml";
  BOOST_CHECK_EXCEPTION(precice::Participant("IDontExist", config, context.rank, context.size),
                        ::precice::Error,
                        precice::testing::errorContains("is not defined"));
}

PRECICE_TEST_SETUP(1_rank)
BOOST_AUTO_TEST_CASE(NoConfig)
{
  PRECICE_TEST();

  BOOST_CHECK_EXCEPTION(precice::Participant("SolverOne", "", context.rank, context.size), ::precice::Error, ::precice::testing::errorContains("unable to open"));
}

PRECICE_TEST_SETUP(1_rank)
BOOST_AUTO_TEST_CASE(MissingConfig)
{
  PRECICE_TEST();

  BOOST_CHECK_EXCEPTION(precice::Participant("SolverOne", "this/file/is/missing.xml", context.rank, context.size), ::precice::Error, ::precice::testing::errorContains("unable to open"));
}

PRECICE_TEST_SETUP(1_rank)
BOOST_AUTO_TEST_CASE(WrongConfigType)
{
  PRECICE_TEST();

  auto notaconfig = precice::testing::getPathToSources() + "/precice/tests/notaconfig.txt";
  // The error depends on the libxml2 version
  BOOST_CHECK_THROW(precice::Participant("SolverOne", notaconfig, context.rank, context.size), ::precice::Error);
}

PRECICE_TEST_SETUP(1_rank)
BOOST_AUTO_TEST_CASE(EmptyConfigFile)
{
  PRECICE_TEST();

  auto notaconfig = precice::testing::getPathToSources() + "/precice/tests/emptyconfig.xml";
  BOOST_CHECK_EXCEPTION(precice::Participant("SolverOne", notaconfig, context.rank, context.size), ::precice::Error, ::precice::testing::errorContains(" is empty."));
}

PRECICE_TEST_SETUP(1_rank)
BOOST_AUTO_TEST_CASE(RankEqualsSize)
{
  PRECICE_TEST();

  auto config = precice::testing::getPathToSources() + "/precice/tests/config-checker.xml";
  BOOST_CHECK_EXCEPTION(precice::Participant("SolverOne", config, context.size, context.size), ::precice::Error, ::precice::testing::errorContains("needs to be smaller than"));
}

PRECICE_TEST_SETUP(1_rank)
BOOST_AUTO_TEST_CASE(RankLargerThanSize)
{
  PRECICE_TEST();

  auto config = precice::testing::getPathToSources() + "/precice/tests/config-checker.xml";
  BOOST_CHECK_EXCEPTION(precice::Participant("SolverOne", config, context.size + 1, context.size), ::precice::Error, ::precice::testing::errorContains("needs to be smaller than"));
}

PRECICE_TEST_SETUP(1_rank)
BOOST_AUTO_TEST_CASE(NegativeRank)
{
  PRECICE_TEST();

  auto config = precice::testing::getPathToSources() + "/precice/tests/config-checker.xml";
  BOOST_CHECK_EXCEPTION(precice::Participant("SolverOne", config, -3, context.size), ::precice::Error, ::precice::testing::errorContains("non-negative number"));
}

PRECICE_TEST_SETUP(1_rank)
BOOST_AUTO_TEST_CASE(NegativeSize)
{
  PRECICE_TEST();

  auto config = precice::testing::getPathToSources() + "/precice/tests/config-checker.xml";
  BOOST_CHECK_EXCEPTION(precice::Participant("SolverOne", config, context.rank, -3), ::precice::Error, ::precice::testing::errorContains("positive number"));
}

PRECICE_TEST_SETUP(1_rank)
BOOST_AUTO_TEST_CASE(ZeroSize)
{
  PRECICE_TEST();

  auto config = precice::testing::getPathToSources() + "/precice/tests/config-checker.xml";
  BOOST_CHECK_EXCEPTION(precice::Participant("SolverOne", config, 0, 0), ::precice::Error, ::precice::testing::errorContains("positive number"));
}

PRECICE_TEST_SETUP(1_rank)
BOOST_AUTO_TEST_CASE(NullptrAsCommunicator)
{
  PRECICE_TEST();

  auto config = precice::testing::getPathToSources() + "/precice/tests/config-checker.xml";
  BOOST_CHECK_EXCEPTION(precice::Participant("SolverOne", config, context.rank, context.size, nullptr), ::precice::Error, ::precice::testing::errorContains("nullptr"));
}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE(MeshDimenstions)

PRECICE_TEST_SETUP(1_rank)
BOOST_AUTO_TEST_CASE(NoName)
{
  PRECICE_TEST();

  auto                 config = precice::testing::getPathToSources() + "/precice/tests/config-checker.xml";
  precice::Participant p("SolverTwo", config, context.rank, context.size);
  BOOST_CHECK_EXCEPTION(p.getMeshDimensions(""),
                        ::precice::Error,
                        precice::testing::errorContains("Available meshes are"));
}

PRECICE_TEST_SETUP(1_rank)
BOOST_AUTO_TEST_CASE(WrongName)
{
  PRECICE_TEST();

  auto                 config = precice::testing::getPathToSources() + "/precice/tests/config-checker.xml";
  precice::Participant p("SolverTwo", config, context.rank, context.size);
  BOOST_CHECK_EXCEPTION(p.getMeshDimensions("FaceCenters"),
                        ::precice::Error,
                        precice::testing::errorContains("Available meshes are"));
}

PRECICE_TEST_SETUP(1_rank)
BOOST_AUTO_TEST_CASE(TypoNameChanged)
{
  PRECICE_TEST();

  auto                 config = precice::testing::getPathToSources() + "/precice/tests/config-checker.xml";
  precice::Participant p("SolverTwo", config, context.rank, context.size);
  BOOST_CHECK_EXCEPTION(p.getMeshDimensions("MeshOno"),
                        ::precice::Error,
                        precice::testing::errorContains("Did you mean mesh"));
}

PRECICE_TEST_SETUP(1_rank)
BOOST_AUTO_TEST_CASE(TypoNameMissing)
{
  PRECICE_TEST();

  auto                 config = precice::testing::getPathToSources() + "/precice/tests/config-checker.xml";
  precice::Participant p("SolverTwo", config, context.rank, context.size);
  BOOST_CHECK_EXCEPTION(p.getMeshDimensions("MeshOn"),
                        ::precice::Error,
                        precice::testing::errorContains("Did you mean mesh"));
}

PRECICE_TEST_SETUP(1_rank)
BOOST_AUTO_TEST_CASE(TypoNameExtra)
{
  PRECICE_TEST();

  auto                  config = precice::testing::getPathToSources() + "/precice/tests/config-checker.xml";
  precice::Participant  p("SolverTwo", config, context.rank, context.size);
  std::array<double, 3> pos{1, 2, 3};
  BOOST_CHECK_EXCEPTION(p.setMeshVertex("MeshOnee", pos),
                        ::precice::Error,
                        precice::testing::errorContains("Did you mean mesh"));
}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE(DataDimensions)

PRECICE_TEST_SETUP(1_rank)
BOOST_AUTO_TEST_CASE(WrongMesh)
{
  PRECICE_TEST();

  auto                 config = precice::testing::getPathToSources() + "/precice/tests/config-checker.xml";
  precice::Participant p("SolverTwo", config, context.rank, context.size);
  BOOST_CHECK_EXCEPTION(p.getDataDimensions("CellCenters", "DataTwo"),
                        ::precice::Error,
                        precice::testing::errorContains("Available meshes are"));
}

PRECICE_TEST_SETUP(1_rank)
BOOST_AUTO_TEST_CASE(TypoMesh)
{
  PRECICE_TEST();

  auto                 config = precice::testing::getPathToSources() + "/precice/tests/config-checker.xml";
  precice::Participant p("SolverTwo", config, context.rank, context.size);
  BOOST_CHECK_EXCEPTION(p.getDataDimensions("MeshToo", "DataTwo"),
                        ::precice::Error,
                        precice::testing::errorContains("Did you mean mesh"));
}

PRECICE_TEST_SETUP(1_rank)
BOOST_AUTO_TEST_CASE(WrongData)
{
  PRECICE_TEST();

  auto                 config = precice::testing::getPathToSources() + "/precice/tests/config-checker.xml";
  precice::Participant p("SolverTwo", config, context.rank, context.size);
  BOOST_CHECK_EXCEPTION(p.getDataDimensions("MeshTwo", "Temperature"),
                        ::precice::Error,
                        precice::testing::errorContains("Available data are"));
}

PRECICE_TEST_SETUP(1_rank)
BOOST_AUTO_TEST_CASE(TypoData)
{
  PRECICE_TEST();

  auto                 config = precice::testing::getPathToSources() + "/precice/tests/config-checker.xml";
  precice::Participant p("SolverTwo", config, context.rank, context.size);
  BOOST_CHECK_EXCEPTION(p.getDataDimensions("MeshTwo", "DataTwi"),
                        ::precice::Error,
                        precice::testing::errorContains("Did you mean data"));
}

PRECICE_TEST_SETUP(1_rank)
BOOST_AUTO_TEST_CASE(WrongMeshAndData)
{
  PRECICE_TEST();

  auto                 config = precice::testing::getPathToSources() + "/precice/tests/config-checker.xml";
  precice::Participant p("SolverTwo", config, context.rank, context.size);
  BOOST_CHECK_EXCEPTION(p.getDataDimensions("CellCenters", "Temperature"),
                        ::precice::Error,
                        precice::testing::errorContains("Available meshes are"));
}

PRECICE_TEST_SETUP(1_rank)
BOOST_AUTO_TEST_CASE(TypoMeshAndData)
{
  PRECICE_TEST();

  auto                 config = precice::testing::getPathToSources() + "/precice/tests/config-checker.xml";
  precice::Participant p("SolverTwo", config, context.rank, context.size);
  BOOST_CHECK_EXCEPTION(p.getDataDimensions("MeshUne", "DataToo"),
                        ::precice::Error,
                        precice::testing::errorContains("Did you mean mesh"));
}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE(SetVertex)

PRECICE_TEST_SETUP(1_rank)
BOOST_AUTO_TEST_CASE(NoName)
{
  PRECICE_TEST();

  auto                  config = precice::testing::getPathToSources() + "/precice/tests/config-checker.xml";
  precice::Participant  p("SolverTwo", config, context.rank, context.size);
  std::array<double, 3> pos{1, 2, 3};
  BOOST_CHECK_EXCEPTION(p.setMeshVertex("", pos),
                        ::precice::Error,
                        precice::testing::errorContains("Available meshes are"));
}

PRECICE_TEST_SETUP(1_rank)
BOOST_AUTO_TEST_CASE(WrongName)
{
  PRECICE_TEST();

  auto                  config = precice::testing::getPathToSources() + "/precice/tests/config-checker.xml";
  precice::Participant  p("SolverTwo", config, context.rank, context.size);
  std::array<double, 3> pos{1, 2, 3};
  BOOST_CHECK_EXCEPTION(p.setMeshVertex("FaceCenters", pos),
                        ::precice::Error,
                        precice::testing::errorContains("Available meshes are"));
}

PRECICE_TEST_SETUP(1_rank)
BOOST_AUTO_TEST_CASE(TypoNameChanged)
{
  PRECICE_TEST();

  auto                  config = precice::testing::getPathToSources() + "/precice/tests/config-checker.xml";
  precice::Participant  p("SolverTwo", config, context.rank, context.size);
  std::array<double, 3> pos{1, 2, 3};
  BOOST_CHECK_EXCEPTION(p.setMeshVertex("MeshOno", pos),
                        ::precice::Error,
                        precice::testing::errorContains("Did you mean mesh"));
}

PRECICE_TEST_SETUP(1_rank)
BOOST_AUTO_TEST_CASE(TypoNameMissing)
{
  PRECICE_TEST();

  auto                  config = precice::testing::getPathToSources() + "/precice/tests/config-checker.xml";
  precice::Participant  p("SolverTwo", config, context.rank, context.size);
  std::array<double, 3> pos{1, 2, 3};
  BOOST_CHECK_EXCEPTION(p.setMeshVertex("MeshOn", pos),
                        ::precice::Error,
                        precice::testing::errorContains("Did you mean mesh"));
}

PRECICE_TEST_SETUP(1_rank)
BOOST_AUTO_TEST_CASE(TypoNameExtra)
{
  PRECICE_TEST();

  auto                  config = precice::testing::getPathToSources() + "/precice/tests/config-checker.xml";
  precice::Participant  p("SolverTwo", config, context.rank, context.size);
  std::array<double, 3> pos{1, 2, 3};
  BOOST_CHECK_EXCEPTION(p.setMeshVertex("MeshOnee", pos),
                        ::precice::Error,
                        precice::testing::errorContains("Did you mean mesh"));
}

PRECICE_TEST_SETUP(1_rank)
BOOST_AUTO_TEST_CASE(EmptyPositions)
{
  PRECICE_TEST();

  auto                  config = precice::testing::getPathToSources() + "/precice/tests/config-checker.xml";
  precice::Participant  p("SolverOne", config, context.rank, context.size);
  std::array<double, 0> pos{};
  BOOST_CHECK_EXCEPTION(p.setMeshVertex("MeshOne", pos),
                        ::precice::Error,
                        precice::testing::errorContains("but found 0"));
}

PRECICE_TEST_SETUP(1_rank)
BOOST_AUTO_TEST_CASE(WrongDimensionality)
{
  PRECICE_TEST();

  auto                  config = precice::testing::getPathToSources() + "/precice/tests/config-checker.xml";
  precice::Participant  p("SolverOne", config, context.rank, context.size);
  std::array<double, 3> pos{1, 2, 3};
  BOOST_CHECK_EXCEPTION(p.setMeshVertex("MeshOne", pos),
                        ::precice::Error,
                        precice::testing::errorContains("but found 3"));
}

/// TODO check this
PRECICE_TEST_SETUP(1_rank)
BOOST_AUTO_TEST_CASE(InfInPosition)
{
  PRECICE_TEST();

  auto                  config = precice::testing::getPathToSources() + "/precice/tests/config-checker.xml";
  precice::Participant  p("SolverOne", config, context.rank, context.size);
  std::array<double, 2> pos{std::numeric_limits<double>::infinity(), 2};
  // Currently allowed
  BOOST_CHECK_NO_THROW(p.setMeshVertex("MeshOne", pos));
}

PRECICE_TEST_SETUP(1_rank)
BOOST_AUTO_TEST_CASE(NanInPosition)
{
  PRECICE_TEST();

  auto                  config = precice::testing::getPathToSources() + "/precice/tests/config-checker.xml";
  precice::Participant  p("SolverOne", config, context.rank, context.size);
  std::array<double, 2> pos{std::numeric_limits<double>::quiet_NaN(), 2};
  // Currently allowed
  BOOST_CHECK_NO_THROW(p.setMeshVertex("MeshOne", pos));
}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE(SetVertices)

PRECICE_TEST_SETUP(1_rank)
BOOST_AUTO_TEST_CASE(NoName)
{
  PRECICE_TEST();

  auto                             config = precice::testing::getPathToSources() + "/precice/tests/config-checker.xml";
  precice::Participant             p("SolverTwo", config, context.rank, context.size);
  std::array<double, 3>            pos{1, 2, 3};
  std::array<precice::VertexID, 1> vids;
  BOOST_CHECK_EXCEPTION(p.setMeshVertices("", pos, vids),
                        ::precice::Error,
                        precice::testing::errorContains("Available meshes are"));
}

PRECICE_TEST_SETUP(1_rank)
BOOST_AUTO_TEST_CASE(WrongName)
{
  PRECICE_TEST();

  auto                             config = precice::testing::getPathToSources() + "/precice/tests/config-checker.xml";
  precice::Participant             p("SolverTwo", config, context.rank, context.size);
  std::array<double, 3>            pos{1, 2, 3};
  std::array<precice::VertexID, 1> vids;
  BOOST_CHECK_EXCEPTION(p.setMeshVertices("FaceCenters", pos, vids),
                        ::precice::Error,
                        precice::testing::errorContains("Available meshes are"));
}

PRECICE_TEST_SETUP(1_rank)
BOOST_AUTO_TEST_CASE(TypoNameChanged)
{
  PRECICE_TEST();

  auto                             config = precice::testing::getPathToSources() + "/precice/tests/config-checker.xml";
  precice::Participant             p("SolverTwo", config, context.rank, context.size);
  std::array<double, 3>            pos{1, 2, 3};
  std::array<precice::VertexID, 1> vids;
  BOOST_CHECK_EXCEPTION(p.setMeshVertices("MeshOno", pos, vids),
                        ::precice::Error,
                        precice::testing::errorContains("Did you mean mesh"));
}

PRECICE_TEST_SETUP(1_rank)
BOOST_AUTO_TEST_CASE(TypoNameMissing)
{
  PRECICE_TEST();

  auto                             config = precice::testing::getPathToSources() + "/precice/tests/config-checker.xml";
  precice::Participant             p("SolverTwo", config, context.rank, context.size);
  std::array<double, 3>            pos{1, 2, 3};
  std::array<precice::VertexID, 1> vids;
  BOOST_CHECK_EXCEPTION(p.setMeshVertices("MeshOn", pos, vids),
                        ::precice::Error,
                        precice::testing::errorContains("Did you mean mesh"));
}

PRECICE_TEST_SETUP(1_rank)
BOOST_AUTO_TEST_CASE(TypoNameExtra)
{
  PRECICE_TEST();

  auto                             config = precice::testing::getPathToSources() + "/precice/tests/config-checker.xml";
  precice::Participant             p("SolverTwo", config, context.rank, context.size);
  std::array<double, 3>            pos{1, 2, 3};
  std::array<precice::VertexID, 1> vids;
  BOOST_CHECK_EXCEPTION(p.setMeshVertices("MeshOnee", pos, vids),
                        ::precice::Error,
                        precice::testing::errorContains("Did you mean mesh"));
}

PRECICE_TEST_SETUP(1_rank)
BOOST_AUTO_TEST_CASE(EmptyPositions)
{
  PRECICE_TEST();

  auto                             config = precice::testing::getPathToSources() + "/precice/tests/config-checker.xml";
  precice::Participant             p("SolverOne", config, context.rank, context.size);
  std::array<double, 0>            pos{};
  std::array<precice::VertexID, 1> vids;
  BOOST_CHECK_EXCEPTION(p.setMeshVertices("MeshOne", pos, vids),
                        ::precice::Error,
                        precice::testing::errorContains("1 vertex indices and 0 position components"));
}

PRECICE_TEST_SETUP(1_rank)
BOOST_AUTO_TEST_CASE(NotEnoughPositions)
{
  PRECICE_TEST();

  auto                             config = precice::testing::getPathToSources() + "/precice/tests/config-checker.xml";
  precice::Participant             p("SolverOne", config, context.rank, context.size);
  std::array<double, 3>            pos{1, 2, 3};
  std::array<precice::VertexID, 2> vids;
  BOOST_CHECK_EXCEPTION(p.setMeshVertices("MeshOne", pos, vids),
                        ::precice::Error,
                        precice::testing::errorContains("2 vertex indices and 3 position components"));
}

PRECICE_TEST_SETUP(1_rank)
BOOST_AUTO_TEST_CASE(TooManyPositions)
{
  PRECICE_TEST();

  auto                             config = precice::testing::getPathToSources() + "/precice/tests/config-checker.xml";
  precice::Participant             p("SolverOne", config, context.rank, context.size);
  std::array<double, 7>            pos{1, 2, 3, 4, 5, 6, 7};
  std::array<precice::VertexID, 3> vids;
  BOOST_CHECK_EXCEPTION(p.setMeshVertices("MeshOne", pos, vids),
                        ::precice::Error,
                        precice::testing::errorContains("3 vertex indices and 7 position components"));
}

/// TODO check this
PRECICE_TEST_SETUP(1_rank)
BOOST_AUTO_TEST_CASE(InfInPosition)
{
  PRECICE_TEST();

  auto                             config = precice::testing::getPathToSources() + "/precice/tests/config-checker.xml";
  precice::Participant             p("SolverOne", config, context.rank, context.size);
  std::array<double, 2>            pos{std::numeric_limits<double>::infinity(), 2};
  std::array<precice::VertexID, 1> vids;
  // Currently allowed
  BOOST_CHECK_NO_THROW(p.setMeshVertices("MeshOne", pos, vids));
}

PRECICE_TEST_SETUP(1_rank)
BOOST_AUTO_TEST_CASE(NanInPosition)
{
  PRECICE_TEST();

  auto                             config = precice::testing::getPathToSources() + "/precice/tests/config-checker.xml";
  precice::Participant             p("SolverOne", config, context.rank, context.size);
  std::array<double, 2>            pos{std::numeric_limits<double>::quiet_NaN(), 2};
  std::array<precice::VertexID, 1> vids;
  // Currently allowed
  BOOST_CHECK_NO_THROW(p.setMeshVertices("MeshOne", pos, vids));
}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE(WriteData)

PRECICE_TEST_SETUP(1_rank)
BOOST_AUTO_TEST_CASE(WrongMesh)
{
  PRECICE_TEST();

  auto                             config = precice::testing::getPathToSources() + "/precice/tests/config-checker.xml";
  precice::Participant             p("SolverTwo", config, context.rank, context.size);
  std::array<double, 6>            pos{1, 2, 3, 4, 5, 6};
  std::array<precice::VertexID, 3> vids;
  p.setMeshVertices("MeshTwo", pos, vids);

  std::array<double, 3> data{1.0, 2.0, 3.0};
  BOOST_CHECK_EXCEPTION(p.writeData("CellCenters", "DataTwo", vids, data),
                        ::precice::Error,
                        precice::testing::errorContains("Available meshes are"));
}

PRECICE_TEST_SETUP(1_rank)
BOOST_AUTO_TEST_CASE(TypoMesh)
{
  PRECICE_TEST();

  auto                             config = precice::testing::getPathToSources() + "/precice/tests/config-checker.xml";
  precice::Participant             p("SolverTwo", config, context.rank, context.size);
  std::array<double, 6>            pos{1, 2, 3, 4, 5, 6};
  std::array<precice::VertexID, 3> vids;
  p.setMeshVertices("MeshTwo", pos, vids);

  std::array<double, 3> data{1.0, 2.0, 3.0};
  BOOST_CHECK_EXCEPTION(p.writeData("MeshToo", "DataTwo", vids, data),
                        ::precice::Error,
                        precice::testing::errorContains("Did you mean mesh"));
}

PRECICE_TEST_SETUP(1_rank)
BOOST_AUTO_TEST_CASE(WrongData)
{
  PRECICE_TEST();

  auto                             config = precice::testing::getPathToSources() + "/precice/tests/config-checker.xml";
  precice::Participant             p("SolverTwo", config, context.rank, context.size);
  std::array<double, 6>            pos{1, 2, 3, 4, 5, 6};
  std::array<precice::VertexID, 3> vids;
  p.setMeshVertices("MeshTwo", pos, vids);

  std::array<double, 3> data{1.0, 2.0, 3.0};
  BOOST_CHECK_EXCEPTION(p.writeData("MeshTwo", "Temperature", vids, data),
                        ::precice::Error,
                        precice::testing::errorContains("Available data are"));
}

PRECICE_TEST_SETUP(1_rank)
BOOST_AUTO_TEST_CASE(TypoData)
{
  PRECICE_TEST();

  auto                             config = precice::testing::getPathToSources() + "/precice/tests/config-checker.xml";
  precice::Participant             p("SolverTwo", config, context.rank, context.size);
  std::array<double, 6>            pos{1, 2, 3, 4, 5, 6};
  std::array<precice::VertexID, 3> vids;
  p.setMeshVertices("MeshTwo", pos, vids);

  std::array<double, 3> data{1.0, 2.0, 3.0};
  BOOST_CHECK_EXCEPTION(p.writeData("MeshTwo", "DataTwi", vids, data),
                        ::precice::Error,
                        precice::testing::errorContains("Did you mean data"));
}

PRECICE_TEST_SETUP(1_rank)
BOOST_AUTO_TEST_CASE(WrongMeshAndData)
{
  PRECICE_TEST();

  auto                             config = precice::testing::getPathToSources() + "/precice/tests/config-checker.xml";
  precice::Participant             p("SolverTwo", config, context.rank, context.size);
  std::array<double, 6>            pos{1, 2, 3, 4, 5, 6};
  std::array<precice::VertexID, 3> vids;
  p.setMeshVertices("MeshTwo", pos, vids);

  std::array<double, 3> data{1.0, 2.0, 3.0};
  BOOST_CHECK_EXCEPTION(p.writeData("CellCenters", "Temperature", vids, data),
                        ::precice::Error,
                        precice::testing::errorContains("Available meshes are"));
}

PRECICE_TEST_SETUP(1_rank)
BOOST_AUTO_TEST_CASE(TypoMeshAndData)
{
  PRECICE_TEST();

  auto                             config = precice::testing::getPathToSources() + "/precice/tests/config-checker.xml";
  precice::Participant             p("SolverTwo", config, context.rank, context.size);
  std::array<double, 6>            pos{1, 2, 3, 4, 5, 6};
  std::array<precice::VertexID, 3> vids;
  p.setMeshVertices("MeshTwo", pos, vids);

  std::array<double, 3> data{1.0, 2.0, 3.0};
  BOOST_CHECK_EXCEPTION(p.writeData("MeshUne", "DataToo", vids, data),
                        ::precice::Error,
                        precice::testing::errorContains("Did you mean mesh"));
}

PRECICE_TEST_SETUP(1_rank)
BOOST_AUTO_TEST_CASE(NoVertices)
{
  PRECICE_TEST();

  auto                             config = precice::testing::getPathToSources() + "/precice/tests/config-checker.xml";
  precice::Participant             p("SolverTwo", config, context.rank, context.size);
  std::array<double, 6>            pos{1, 2, 3, 4, 5, 6};
  std::array<precice::VertexID, 3> vids;
  p.setMeshVertices("MeshTwo", pos, vids);

  std::array<double, 3> data{1.0, 2.0, 3.0};
  // Call is ignored
  BOOST_CHECK_NO_THROW(p.writeData("MeshTwo", "DataTwo", vids, data));
}

PRECICE_TEST_SETUP(1_rank)
BOOST_AUTO_TEST_CASE(DataTooSmall)
{
  PRECICE_TEST();

  auto                             config = precice::testing::getPathToSources() + "/precice/tests/config-checker.xml";
  precice::Participant             p("SolverTwo", config, context.rank, context.size);
  std::array<double, 6>            pos{1, 2, 3, 4, 5, 6};
  std::array<precice::VertexID, 3> vids;
  p.setMeshVertices("MeshTwo", pos, vids);

  std::array<double, 1> data{1.0};
  BOOST_CHECK_EXCEPTION(p.writeData("MeshTwo", "DataTwo", vids, data),
                        ::precice::Error,
                        precice::testing::errorContains("3 vertex indices and 1 data components"));
}

PRECICE_TEST_SETUP(1_rank)
BOOST_AUTO_TEST_CASE(DataTooLarge)
{
  PRECICE_TEST();

  auto                             config = precice::testing::getPathToSources() + "/precice/tests/config-checker.xml";
  precice::Participant             p("SolverTwo", config, context.rank, context.size);
  std::array<double, 6>            pos{1, 2, 3, 4, 5, 6};
  std::array<precice::VertexID, 3> vids;
  p.setMeshVertices("MeshTwo", pos, vids);

  std::array<double, 4> data{1.0, 2.0, 3.0, 4.0};
  BOOST_CHECK_EXCEPTION(p.writeData("MeshTwo", "DataTwo", vids, data),
                        ::precice::Error,
                        precice::testing::errorContains("3 vertex indices and 4 data components"));
}

PRECICE_TEST_SETUP(1_rank)
BOOST_AUTO_TEST_CASE(RepeatVertex)
{
  PRECICE_TEST();

  auto                             config = precice::testing::getPathToSources() + "/precice/tests/config-checker.xml";
  precice::Participant             p("SolverTwo", config, context.rank, context.size);
  std::array<double, 6>            pos{1, 2, 3, 4, 5, 6};
  std::array<precice::VertexID, 3> vids;
  p.setMeshVertices("MeshTwo", pos, vids);
  vids.back() = vids.front(); // last == first

  std::array<double, 3> data{1.0, 2.0, 3.0};
  // This is currently allowed but is in discussion of being forbidden
  BOOST_CHECK_NO_THROW(p.writeData("MeshTwo", "DataTwo", vids, data));
}

PRECICE_TEST_SETUP(1_rank)
BOOST_AUTO_TEST_CASE(WrongVertexIDSingle)
{
  PRECICE_TEST();

  auto                             config = precice::testing::getPathToSources() + "/precice/tests/config-checker.xml";
  precice::Participant             p("SolverTwo", config, context.rank, context.size);
  std::array<double, 6>            pos{1, 2, 3, 4, 5, 6};
  std::array<precice::VertexID, 3> vids;
  p.setMeshVertices("MeshTwo", pos, vids);
  vids[1] = 666;

  std::array<double, 3> data{1.0, 2.0, 3.0};
  // This is currently allowed but is in discussion of being forbidden
  BOOST_CHECK_EXCEPTION(p.writeData("MeshTwo", "DataTwo", vids, data),
                        ::precice::Error,
                        precice::testing::errorContains("invalid Vertex ID at"));
}

PRECICE_TEST_SETUP(1_rank)
BOOST_AUTO_TEST_CASE(WrongVertexIDRange)
{
  PRECICE_TEST();

  auto                             config = precice::testing::getPathToSources() + "/precice/tests/config-checker.xml";
  precice::Participant             p("SolverTwo", config, context.rank, context.size);
  std::array<double, 10>           pos{1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
  std::array<precice::VertexID, 5> vids;
  p.setMeshVertices("MeshTwo", pos, vids);
  std::fill_n(&vids[1], 3, 666);

  std::array<double, 5> data{1.0, 2.0, 3.0, 4.0, 5.0};
  // This is currently allowed but is in discussion of being forbidden
  BOOST_CHECK_EXCEPTION(p.writeData("MeshTwo", "DataTwo", vids, data),
                        ::precice::Error,
                        precice::testing::errorContains("invalid Vertex ID at"));
}

PRECICE_TEST_SETUP(1_rank)
BOOST_AUTO_TEST_CASE(NaNInData)
{
  PRECICE_TEST();

  auto                             config = precice::testing::getPathToSources() + "/precice/tests/config-checker.xml";
  precice::Participant             p("SolverTwo", config, context.rank, context.size);
  std::array<double, 6>            pos{1, 2, 3, 4, 5, 6};
  std::array<precice::VertexID, 3> vids;
  p.setMeshVertices("MeshTwo", pos, vids);

  std::array<double, 3> data{std::numeric_limits<double>::quiet_NaN(), 2.0, 3.0};
#ifdef NDEBUG
  BOOST_CHECK_NO_THROW(p.writeData("MeshTwo", "DataTwo", vids, data));
#else
  BOOST_CHECK_EXCEPTION(p.writeData("MeshTwo", "DataTwo", vids, data),
                        ::precice::Error,
                        precice::testing::errorContains("NaN"));
#endif
}

PRECICE_TEST_SETUP(1_rank)
BOOST_AUTO_TEST_CASE(InfInData)
{
  PRECICE_TEST();

  auto                             config = precice::testing::getPathToSources() + "/precice/tests/config-checker.xml";
  precice::Participant             p("SolverTwo", config, context.rank, context.size);
  std::array<double, 6>            pos{1, 2, 3, 4, 5, 6};
  std::array<precice::VertexID, 3> vids;
  p.setMeshVertices("MeshTwo", pos, vids);

  std::array<double, 3> data{std::numeric_limits<double>::infinity(), 2.0, 3.0};
#ifdef NDEBUG
  BOOST_CHECK_NO_THROW(p.writeData("MeshTwo", "DataTwo", vids, data));
#else
  BOOST_CHECK_EXCEPTION(p.writeData("MeshTwo", "DataTwo", vids, data),
                        ::precice::Error,
                        precice::testing::errorContains("infinity"));
#endif
}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE(NotBeforeInitialize)

PRECICE_TEST_SETUP(1_rank)
BOOST_AUTO_TEST_CASE(ReadData)
{
  PRECICE_TEST();

  auto                             config = precice::testing::getPathToSources() + "/precice/tests/config-checker.xml";
  precice::Participant             p("SolverTwo", config, context.rank, context.size);
  std::array<precice::VertexID, 0> vids;
  std::array<double, 0>            data;
  BOOST_CHECK_EXCEPTION(p.readData("MeshTwo", "DataTwo", vids, 0.0, data),
                        ::precice::Error,
                        precice::testing::errorContains("cannot be called"));
}

PRECICE_TEST_SETUP(1_rank)
BOOST_AUTO_TEST_CASE(GetMaxTimeStepSize)
{
  PRECICE_TEST();

  auto                 config = precice::testing::getPathToSources() + "/precice/tests/config-checker.xml";
  precice::Participant p("SolverTwo", config, context.rank, context.size);
  BOOST_CHECK_EXCEPTION(p.getMaxTimeStepSize(),
                        ::precice::Error,
                        precice::testing::errorContains("has to be called before"));
}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE(MeshConnectivity)

PRECICE_TEST_SETUP(1_rank)
BOOST_AUTO_TEST_CASE(UnknownMesh)
{
  PRECICE_TEST();

  auto                 config = precice::testing::getPathToSources() + "/precice/tests/config-checker.xml";
  precice::Participant p("SolverTwo", config, context.rank, context.size);
  BOOST_CHECK_EXCEPTION(p.requiresMeshConnectivityFor("CellCenters"),
                        ::precice::Error,
                        precice::testing::errorContains("Available meshes are"));
}

PRECICE_TEST_SETUP(1_rank)
BOOST_AUTO_TEST_CASE(TypoedMesh)
{
  PRECICE_TEST();

  auto                 config = precice::testing::getPathToSources() + "/precice/tests/config-checker.xml";
  precice::Participant p("SolverTwo", config, context.rank, context.size);
  BOOST_CHECK_EXCEPTION(p.requiresMeshConnectivityFor("MoshOne"),
                        ::precice::Error,
                        precice::testing::errorContains("Did you mean mesh"));
}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE(ResetMesh)

PRECICE_TEST_SETUP(1_rank)
BOOST_AUTO_TEST_CASE(NotEnabled)
{
  PRECICE_TEST();

  auto                 config = precice::testing::getPathToSources() + "/precice/tests/config-checker.xml";
  precice::Participant p("SolverTwo", config, context.rank, context.size);
  BOOST_CHECK_EXCEPTION(p.resetMesh("MeshOne"),
                        ::precice::Error,
                        precice::testing::errorContains("unlock the full API"));
}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE_END()
BOOST_AUTO_TEST_SUITE_END()
