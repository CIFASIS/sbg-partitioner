/*****************************************************************************

 This file is part of QSSModelInstance Solver.

 QSSModelInstance Solver is free software: you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation, either version 3 of the License, or
 (at your option) any later version.

 QSSModelInstance Solver is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.

 You should have received a copy of the GNU General Public License
 along with QSSModelInstance Solver.  If not, see
 <http://www.gnu.org/licenses/>.

 ******************************************************************************/

#include <gtest/gtest.h>

#include <cstdlib>
#include <fstream>

class ISBGPartTests : public testing::TestWithParam<const char*> {
};

TEST_P(ISBGPartTests, GenerateGraph)
{
  const std::string NAME = GetParam();
  std::cout << "Testing model: " << NAME << std::endl;
  const std::string MODEL = " ./system/gt_data/" + NAME + "/" + NAME + ".json";
  const std::string FOLDER_CMD = "mkdir ./system/test_data/" + NAME;
  const std::string SBG_PART = "../../bin/sbg-partitioner";
  const std::string ARGS = " -p 2 -f " + MODEL;
  const std::string TEST_CMD = "./system/test_results.sh " + NAME;
  const std::string RESULT_FILE = "./system/test_data/" + NAME + "/" + NAME + ".passed";
  const std::string COMP_CMD = SBG_PART + ARGS + MODEL + "; mv SBG.log ./system/test_data/" + NAME + "/";

  std::cout << "Setup data folders for " << NAME << std::endl;
  std::system(FOLDER_CMD.c_str());
  std::cout << "Testing model: " << NAME << std::endl;
  std::system(COMP_CMD.c_str());
  std::cout << "Check results for model: " << NAME << std::endl;
  std::system(TEST_CMD.c_str());

  std::ifstream result(RESULT_FILE.c_str());
  EXPECT_TRUE(result.good());

  EXPECT_TRUE(true);
}

const char* models[] = {"air_conditioners", "toy_example", "advection2D"};

INSTANTIATE_TEST_SUITE_P(Graphs, ISBGPartTests, testing::ValuesIn(models));
