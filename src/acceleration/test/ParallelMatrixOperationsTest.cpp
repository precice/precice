#include <boost/test/unit_test_log.hpp>
#ifndef PRECICE_NO_MPI

#include <Eigen/Core>
#include <algorithm>
#include <math.h>
#include <memory>
#include <ostream>
#include <stdlib.h>
#include <string>
#include <vector>
#include "acceleration/impl/ParallelMatrixOperations.hpp"
#include "com/Communication.hpp"
#include "com/MPIPortsCommunication.hpp"
#include "com/SharedPointer.hpp"
#include "testing/TestContext.hpp"
#include "testing/Testing.hpp"
#include "utils/MasterSlave.hpp"

BOOST_AUTO_TEST_SUITE(AccelerationTests)

using namespace precice;
using namespace precice::acceleration;
using namespace precice::acceleration::impl;

BOOST_AUTO_TEST_SUITE(ParallelMatrixOperationsTests)

void validate_result_equals_reference(
    Eigen::MatrixXd &result_local,
    Eigen::MatrixXd &reference_global,
    int              offset,
    bool             partitionedRowWise)
{
  for (int i = 0; i < result_local.rows(); i++) {
    for (int j = 0; j < result_local.cols(); j++) {
      if (partitionedRowWise) {
        BOOST_TEST(testing::equals(result_local(i, j), reference_global(i + offset, j)),
                   "Failed: (" << i + offset << ", " << j << ") result, reference:" << result_local(i, j) << ", " << reference_global(i + offset, j));
      } else {
        BOOST_TEST(testing::equals(result_local(i, j), reference_global(i, j + offset)),
                   "Failed: (" << i << ", " << j + offset << ") result, reference:" << result_local(i, j) << ", " << reference_global(i, j + offset));
      }
    }
  }
}

BOOST_AUTO_TEST_CASE(ParVectorOperations)
{
  PRECICE_TEST(""_on(4_ranks).setupMasterSlaves());
  int              n_global = 10;
  int              n_local;
  double           a = 0;
  std::vector<int> vertexOffsets{0, 3, 7, 7, 10};

  // definition of vectors
  Eigen::VectorXd vec1(n_global);
  Eigen::VectorXd vec2(n_global);

  // l2norm: 1.502540907218387
  vec1 << 0.422885689100085,
      0.094229338887735,
      0.598523668756741,
      0.470924256358334,
      0.695949313301608,
      0.699887849928292,
      0.638530758271838,
      0.033603836066429,
      0.068806099118051,
      0.319599735180496;

  // l2norm: 6.076423472407709
  vec2 << 2.104516338543882,
      2.345060145175945,
      1.506380184916943,
      2.326775951409934,
      1.518391207198873,
      1.512172623678050,
      1.836391544384709,
      2.439901436460795,
      1.803781441584700,
      1.462976489192458;

  // <vec1, vec2> = 7.069617899295469

  if (context.isMaster()) {
    n_local = 3;
    a       = 1;
  } else if (context.isRank(1)) {
    n_local = 4;
    a       = 2;
  } else if (context.isRank(2)) {
    n_local = 0;
    a       = 3;
  } else {
    BOOST_REQUIRE(context.isRank(3));
    n_local = 3;
    a       = 4;
  }

  // ------ test Allreduce/ Reduce ---------------------------
  std::vector<double> aa   = {a, a};
  std::vector<double> res2 = {0, 0};
  std::vector<double> res3 = {0, 0};

  double res1  = 0;
  int    iaa   = (int) a;
  int    ires1 = 0, ires2 = 0;

  utils::MasterSlave::allreduceSum(a, res1, 1);
  utils::MasterSlave::allreduceSum(iaa, ires2, 1);
  utils::MasterSlave::allreduceSum(aa.data(), res2.data(), 2);

  utils::MasterSlave::reduceSum(aa.data(), res3.data(), 2);
  utils::MasterSlave::reduceSum(iaa, ires1, 1);

  BOOST_TEST(testing::equals(res1, 10.));
  BOOST_TEST(testing::equals(ires2, 10));
  BOOST_TEST(testing::equals(res2.at(0), 10.));
  BOOST_TEST(testing::equals(res2.at(1), 10.));

  if (utils::MasterSlave::isMaster()) {
    BOOST_TEST(testing::equals(res3.at(0), 10.));
    BOOST_TEST(testing::equals(res3.at(1), 10.));
    BOOST_TEST(testing::equals(ires1, 10));
  }
  // ---------------------------------------------------------

  Eigen::VectorXd vec1_local(n_local);
  Eigen::VectorXd vec2_local(n_local);

  // partition and distribute
  int off = vertexOffsets.at(context.rank);
  for (int i = 0; i < n_local; i++) {
    vec1_local(i) = vec1(i + off);
    vec2_local(i) = vec2(i + off);
  }

  double normVec1   = utils::MasterSlave::l2norm(vec1_local);
  double normVec2   = utils::MasterSlave::l2norm(vec2_local);
  double dotproduct = utils::MasterSlave::dot(vec1_local, vec2_local);

  //  std::cout<<"l2norm vec1: "<<normVec1<<'\n';
  //  std::cout<<"l2norm vec2: "<<normVec2<<'\n';
  //  std::cout<<"dotproduct: "<<dotproduct<<'\n';

  // validate
  BOOST_TEST(testing::equals(normVec1, 1.502540907218387));
  BOOST_TEST(testing::equals(normVec2, 6.076423472407709));
  BOOST_TEST(testing::equals(dotproduct, 7.069617899295469));
}

BOOST_AUTO_TEST_CASE(ParallelMatrixMatrixOp)
{
  PRECICE_TEST(""_on(4_ranks).setupMasterSlaves());

  int              n_global = 10, m_global = 5;
  int              n_local;
  std::vector<int> vertexOffsets{0, 3, 7, 7, 10};

  // Definition of Matrices
  Eigen::MatrixXd W_global(n_global, m_global);
  Eigen::MatrixXd J_global(n_global, n_global);
  Eigen::MatrixXd Z_global(m_global, n_global);
  Eigen::MatrixXd WZ_global(n_global, n_global);
  Eigen::MatrixXd JW_global(n_global, m_global);
  Eigen::MatrixXd res_global(n_global, 1);
  Eigen::MatrixXd Jres_global(n_global, 1);

  J_global << 0.162182308193243, 0.450541598502498, 0.106652770180584, 0.431413827463545, 0.853031117721894, 0.417267069084370, 0.780252068321138, 0.234779913372406, 0.547008892286345, 0.929385970968730,
      0.794284540683907, 0.083821377996933, 0.961898080855054, 0.910647594429523, 0.622055131485066, 0.049654430325742, 0.389738836961253, 0.353158571222071, 0.296320805607773, 0.775712678608402,
      0.311215042044805, 0.228976968716819, 0.004634224134067, 0.181847028302852, 0.350952380892271, 0.902716109915281, 0.241691285913833, 0.821194040197959, 0.744692807074156, 0.486791632403172,
      0.528533135506213, 0.913337361501670, 0.774910464711502, 0.263802916521990, 0.513249539867053, 0.944787189721646, 0.403912145588115, 0.015403437651555, 0.188955015032545, 0.435858588580919,
      0.165648729499781, 0.152378018969223, 0.817303220653433, 0.145538980384717, 0.401808033751942, 0.490864092468080, 0.096454525168389, 0.043023801657808, 0.686775433365315, 0.446783749429806,
      0.601981941401637, 0.825816977489547, 0.868694705363510, 0.136068558708664, 0.075966691690842, 0.489252638400019, 0.131973292606335, 0.168990029462704, 0.183511155737270, 0.306349472016557,
      0.262971284540144, 0.538342435260057, 0.084435845510910, 0.869292207640089, 0.239916153553658, 0.337719409821377, 0.942050590775485, 0.649115474956452, 0.368484596490336, 0.508508655381127,
      0.654079098476782, 0.996134716626885, 0.399782649098896, 0.579704587365570, 0.123318934835166, 0.900053846417662, 0.956134540229802, 0.731722385658670, 0.625618560729690, 0.510771564172110,
      0.689214503140008, 0.078175528753184, 0.259870402850654, 0.549860201836332, 0.183907788282417, 0.369246781120215, 0.575208595078466, 0.647745963136307, 0.780227435151377, 0.817627708322262,
      0.748151592823709, 0.442678269775446, 0.800068480224308, 0.144954798223727, 0.239952525664903, 0.111202755293787, 0.059779542947156, 0.450923706430945, 0.081125768865785, 0.794831416883453;

  Z_global << 0.644318130193692, 0.939001561999887, 0.207742292733028, 0.194764289567049, 0.311102286650413, 0.979748378356085, 0.594896074008614, 0.117417650855806, 0.085515797090044, 0.730330862855453,
      0.378609382660268, 0.875942811492984, 0.301246330279491, 0.225921780972399, 0.923379642103244, 0.438869973126103, 0.262211747780845, 0.296675873218327, 0.262482234698333, 0.488608973803579,
      0.811580458282477, 0.550156342898422, 0.470923348517591, 0.170708047147859, 0.430207391329584, 0.111119223440599, 0.602843089382083, 0.318778301925882, 0.801014622769739, 0.578525061023439,
      0.532825588799455, 0.622475086001227, 0.230488160211558, 0.227664297816554, 0.184816320124136, 0.258064695912067, 0.711215780433683, 0.424166759713807, 0.029220277562146, 0.237283579771521,
      0.350727103576883, 0.587044704531417, 0.844308792695389, 0.435698684103899, 0.904880968679893, 0.408719846112552, 0.221746734017240, 0.507858284661118, 0.928854139478045, 0.458848828179931;

  W_global << 0.963088539286913, 0.037738866239552, 0.106761861607241, 0.030540946304637, 0.182922469414914,
      0.546805718738968, 0.885168008202475, 0.653757348668560, 0.744074260367462, 0.239932010568717,
      0.521135830804001, 0.913286827639239, 0.494173936639270, 0.500022435590201, 0.886511933076101,
      0.231594386708524, 0.796183873585212, 0.779051723231275, 0.479922141146060, 0.028674152464106,
      0.488897743920167, 0.098712278655574, 0.715037078400694, 0.904722238067363, 0.489901388512224,
      0.624060088173690, 0.261871183870716, 0.903720560556316, 0.609866648422558, 0.167927145682257,
      0.679135540865748, 0.335356839962797, 0.890922504330789, 0.617666389588455, 0.978680649641159,
      0.395515215668593, 0.679727951377338, 0.334163052737496, 0.859442305646212, 0.712694471678914,
      0.367436648544477, 0.136553137355370, 0.698745832334794, 0.805489424529686, 0.500471624154843,
      0.987982003161633, 0.721227498581740, 0.197809826685929, 0.576721515614685, 0.471088374541939;

  WZ_global << 0.801898401838153, 1.122529091861272, 0.423201945415435, 0.300978558234156, 0.551563616222783, 1.054655768517680, 0.709477478632269, 0.264166315674375, 0.348583586200086, 0.874757870918319,
      1.698638905064191, 2.252495235306823, 1.062194901100270, 0.692015805803183, 1.623336848442643, 1.286934957259726, 1.533908621128989, 0.972679319556613, 1.047374497419335, 1.496714260553646,
      1.659966647284183, 2.392880960790855, 1.479843373796984, 0.892278843543251, 2.112634376079842, 1.457681533373469, 1.399610510316321, 1.151987966448499, 1.518178535932093, 1.638155802609954,
      1.348697901510644, 1.659051863989714, 0.789659275602735, 0.479726420963902, 1.257027472012560, 0.798463699122210, 1.163875903906441, 0.729876025630416, 0.893478486537412, 1.135898804351492,
      1.586570049108119, 1.789685309000292, 1.090184919236369, 0.659006002219424, 1.161370228579310, 1.035482281126420, 1.499868794445182, 0.947182660718915, 1.121957022950641, 1.258422095146877,
      1.618531220549951, 1.790772713269325, 0.916463926569826, 0.546989894699598, 1.089407693306137, 1.052790194022541, 1.455702372306257, 0.783021389674577, 1.019797210461514, 1.328312461070188,
      1.959962169636355, 2.380620645499870, 1.630339868211612, 0.927153862633654, 1.903968068499247, 1.470962702535233, 1.685349366051057, 1.222266173248617, 1.786843939598552, 1.770901564449963,
      1.491283328084721, 2.103999079781337, 1.244121457180890, 0.793826284013736, 1.698194683860993, 1.236033609893653, 1.384257591822518, 1.081117931366570, 1.167011155147404, 1.345250416067541,
      1.460249196277507, 1.644252093240875, 1.054732358743229, 0.623131414603751, 1.142741229794316, 0.909989695162785, 1.359481290557063, 0.902231079036841, 1.115371804173721, 1.160083620581694,
      1.542692246994852, 2.303841728541361, 1.046337587563584, 0.725683826909555, 1.591295952328042, 1.647853963194485, 1.410744976161187, 0.876907034839234, 0.886670387283173, 1.541394813254221;

  JW_global << 2.977456624461038, 2.205533591974792, 3.027360939693343, 3.287118956132616, 2.375186418688407,
      3.137723846199826, 2.752792907538314, 2.639787055824843, 2.828009251522860, 2.504187495011176,
      2.447896971257707, 1.726501702097287, 2.500010777571874, 2.873160045582897, 1.868530866136191,
      3.094354526680215, 2.530311421424636, 3.046051246285873, 2.916644286334997, 2.126621716340908,
      1.981542330537719, 1.649251376583201, 2.034057335063531, 2.167256602621555, 1.754106473495030,
      2.384569365251443, 2.196160002079869, 1.998784819449047, 2.050878308167098, 1.687398579448636,
      2.655317993318163, 2.542013168246354, 2.989992698710950, 3.020860435534452, 2.259852349870721,
      3.812464252001572, 3.252800808239158, 3.906358255932346, 3.917717391998067, 2.952210341980682,
      3.030961610967875, 2.214613005088035, 2.582511677568232, 2.876621280799403, 2.343324946211924,
      2.633853068910827, 2.229873948587862, 1.567502413943023, 2.054972261427936, 1.887633072097274;

  res_global << 0.422885689100085,
      0.094229338887735,
      0.598523668756741,
      0.470924256358334,
      0.695949313301608,
      0.699887849928292,
      0.638530758271838,
      0.033603836066429,
      0.068806099118051,
      0.319599735180496;

  Jres_global << 2.104516338543882,
      2.345060145175945,
      1.506380184916943,
      2.326775951409934,
      1.518391207198873,
      1.512172623678050,
      1.836391544384709,
      2.439901436460795,
      1.803781441584700,
      1.462976489192458;

  if (context.isMaster()) {
    n_local = 3;
  } else if (context.isRank(1)) {
    n_local = 4;
  } else if (context.isRank(2)) {
    n_local = 0;
  } else {
    BOOST_REQUIRE(context.isRank(3));
    n_local = 3;
  }

  Eigen::MatrixXd W_local(n_local, m_global);
  Eigen::MatrixXd J_local(n_global, n_local);
  Eigen::MatrixXd Z_local(m_global, n_local);
  Eigen::MatrixXd WZ_local(n_global, n_local);
  Eigen::MatrixXd JW_local(n_local, m_global);
  Eigen::MatrixXd res_local(n_local, 1);
  Eigen::VectorXd res_local_vec(n_local);
  Eigen::MatrixXd Jres_local(n_local, 1);

  // partition and distribute matrices

  int off = vertexOffsets.at(context.rank);
  for (int i = 0; i < n_global; i++)
    for (int j = 0; j < n_local; j++) {
      J_local(i, j)  = J_global(i, j + off);
      WZ_local(i, j) = WZ_global(i, j + off);
    }

  for (int i = 0; i < n_local; i++)
    for (int j = 0; j < m_global; j++) {
      W_local(i, j)  = W_global(i + off, j);
      JW_local(i, j) = JW_global(i + off, j);
    }

  for (int i = 0; i < m_global; i++)
    for (int j = 0; j < n_local; j++) {
      Z_local(i, j) = Z_global(i, j + off);
    }

  for (int i = 0; i < n_local; i++) {
    res_local(i)     = res_global(i + off);
    res_local_vec(i) = res_global(i + off);
    Jres_local(i)    = Jres_global(i + off);
  }

  // initialize ParallelMatrixOperations object
  ParallelMatrixOperations parMatrixOps{};
  parMatrixOps.initialize(true);

  /*
   * test parallel multiplications
   */
  BOOST_TEST_MESSAGE("Test 1");
  // 1.) multiply JW = J * W (n x m), parallel: (n_local x m)
  Eigen::MatrixXd resJW_local(n_local, m_global);
  parMatrixOps.multiply(J_local, W_local, resJW_local, vertexOffsets, n_global, n_global, m_global);
  validate_result_equals_reference(resJW_local, JW_global, vertexOffsets.at(context.rank), true);

  BOOST_TEST_MESSAGE("Test 2");
  // 2.) multiply WZ = W * Z (n x n), parallel: (n_global x n_local)
  Eigen::MatrixXd resWZ_local(n_global, n_local);
  parMatrixOps.multiply(W_local, Z_local, resWZ_local, vertexOffsets, n_global, m_global, n_global);
  validate_result_equals_reference(resWZ_local, WZ_global, vertexOffsets.at(context.rank), false);

  BOOST_TEST_MESSAGE("Test 3");
  // 3.) multiply Jres = J * res (n x 1), parallel: (n_local x 1)
  Eigen::MatrixXd resJres_local(n_local, 1);
  parMatrixOps.multiply(J_local, res_local, resJres_local, vertexOffsets, n_global, n_global, 1);
  validate_result_equals_reference(resJres_local, Jres_global, vertexOffsets.at(context.rank), true);

  BOOST_TEST_MESSAGE("Test 4");
  // 4.) multiply JW = J * W (n x m), parallel: (n_local x m) with block-wise multiplication
  Eigen::MatrixXd resJW_local2(n_local, m_global);
  parMatrixOps.multiply(J_local, W_local, resJW_local2, vertexOffsets, n_global, n_global, m_global, false);
  validate_result_equals_reference(resJW_local2, JW_global, vertexOffsets.at(context.rank), true);

  BOOST_TEST_MESSAGE("Test 5");
  // 5.) multiply Jres = J * res (n x 1), parallel: (n_local x 1) with block-wise multiplication
  Eigen::VectorXd resJres_local2(n_local); // use the function with parameter of type Eigen::VectorXd
  parMatrixOps.multiply(J_local, res_local_vec, resJres_local2, vertexOffsets, n_global, n_global, 1, false);
  Eigen::MatrixXd matrix_cast = resJres_local2;
  validate_result_equals_reference(matrix_cast, Jres_global, vertexOffsets.at(context.rank), true);
}

BOOST_AUTO_TEST_SUITE_END()
BOOST_AUTO_TEST_SUITE_END()

#endif // PRECICE_NO_MPI
