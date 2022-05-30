// tests.cpp
#include "../../wasm/p528json.cpp"
#include <gtest/gtest.h>

TEST(P528Test, Debug) {

  const auto expectedJson =
      "{\"A__db\":[126.46085394936013,163.0815203486053,236.28197179677775,268."
      "19292392138044,295.0633534171327,323.82517691745284],\"A_a__db\":0,\"A_"
      "fs__db\":[126.425056055803,157.56187253649026,163.55576213633907,167."
      "04435177821983,169.5660332409185,171.5603915860845],"
      "\"ModeOfPropagation\":[1,1,3,3,3,3],\"d__km\":[0.0,360.0008519649057,"
      "720.0,1080.0,1440.0,1800.0],\"theta_h1__rad\":0,\"warnings\":0}";

  std::string outputJson =
      calcP528("{\"h_1__user_units\":5,\"h_2__user_units\":10000,\"time\":50,"
               "\"Polarization\":0,\"f__mhz\":5000,\"nbr_of_points\":5}");

  EXPECT_EQ(expectedJson, outputJson);
}

int main(int argc, char **argv) {
  testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}