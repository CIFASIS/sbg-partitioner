/*****************************************************************************

 This file is part of SBG Partitioner.

 SBG Partitioner is free software: you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation, either version 3 of the License, or
 (at your option) any later version.

 SBG Partitioner is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.

 You should have received a copy of the GNU General Public License
 along with SBG Partitioner.  If not, see <http://www.gnu.org/licenses/>.

 ******************************************************************************/

#include <getopt.h>
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <list>
#include <optional>
#include <string>

#include "build_sb_graph.hpp"
#include "partition_graph.hpp"

#include "kernighan_lin_sbg.hpp"

using namespace std;

using namespace sbg_partitioner;

void usage()
{
  cout << "Usage sbg-partitioner" << endl;
  cout << endl;
  cout << "-h, --help      Display this information and exit" << endl;
  cout << "-v, --version   Display version information and exit" << endl;
  cout << endl;
  cout << "SBG Partitioner home page: https://github.com/CIFASIS/sbg-partitioner " << endl;
}

void version()
{
  cout << "SBG Partitioner 1.0.0" << endl;
  cout << "License GPLv3+: GNU GPL version 3 or later <http://gnu.org/licenses/gpl.html>" << endl;
  cout << "This is free software: you are free to change and redistribute it." << endl;
  cout << "There is NO WARRANTY, to the extent permitted by law." << endl;
}


PartitionGraph create_initial_partition(SBG::LIB::CanonSBG& graph, unsigned number_of_partitions,
               sbg_partitioner::PartitionAlgorithm algorithm, bool pre_order)
{
  PartitionGraph partition_graph(graph, number_of_partitions, algorithm, pre_order);

  return partition_graph;
}


int main(int argc, char** argv)
{
  int ret = 0;
  int opt;
  optional<string> filename = nullopt;
  optional<unsigned> number_of_partitions = nullopt;

  while (true) {

    static struct option long_options[] = {
      {"filename", required_argument, 0, 'f'},
      {"partitions", required_argument, 0, 'p'},
      {"version", no_argument, 0, 'v'},
      {"help", no_argument, 0, 'h'}
    };

    int option_index = 0;
    opt = getopt_long(argc, argv, "f:p:vh:", long_options, &option_index);
    if (opt == EOF) break;

    switch (opt) {
    case 'f':
      if (optarg) {
        filename = string(optarg);
      }
      break;

    case 'p':
      if (optarg) {
        number_of_partitions = atoi(optarg);
      }
      break;

    case 'v':
      version();
      exit(0);

    case 'h':
      usage();
      exit(0);

    case '?':
      usage();
      exit(-1);
      break;

    default:
      abort();
    }
  }

  if (!filename or !number_of_partitions) {
    usage();
    exit(1);
  }

  cout << "filename is " << *filename << endl;
  cout << "number of partitions is " << *number_of_partitions << endl;

  auto sb_graph = build_sb_graph(filename->c_str());

  PartitionGraph pgraph = create_initial_partition(sb_graph, *number_of_partitions, sbg_partitioner::GREEDY, true);
  cout << pgraph << endl;

  // PartitionGraph pgraph2 = create_initial_partition(sb_graph, *number_of_partitions, sbg_partitioner::GREEDY, false);
  // cout << pgraph2 << endl;

  // auto adj = get_adjacents(sb_graph, sb_graph.V()[0]);
  // cout << "Adjacents of " << sb_graph.V()[0] << " are " << adj << endl;

  auto p = pgraph.partitions()[0];
  // unsigned ic = pgraph.internal_cost(p[0], 1);
  // cout << "Internal cost for " << p[0] << " is " << ic << endl;

  unsigned ec = pgraph.internal_cost(p[0], 0);
  cout << "External cost for " << p[0] << " is " << ec << endl;

  int g = pgraph.gain({1, pgraph.partitions()[1][0]}, {0, pgraph.partitions()[0][0]});
  cout << "Gain between " << pgraph.partitions()[1] << " and " << pgraph.partitions()[0] << " is " << g << endl;

  cout << "Exit code: " << ret << endl;

  auto s =  generate_gain_matrix(pgraph.partitions()[1], pgraph.partitions()[0], sb_graph);
  cout << "let's see...\n";
  for (const auto &x : s) {
    cout << x << ", ";
  }
  std::cout << std::endl;

  return ret;
}
