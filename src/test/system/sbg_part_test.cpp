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
  std::cout << "Testing graph: " << NAME << std::endl;
  EXPECT_TRUE(true);
}

const char* models[] = {"dummy"};

INSTANTIATE_TEST_SUITE_P(Graphs, ISBGPartTests, testing::ValuesIn(models));
