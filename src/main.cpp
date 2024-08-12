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
#include <vector>

#include "build_sb_graph.hpp"
#include "kernighan_lin_sbg.hpp"
#include "partition_graph.hpp"


#include "partition_strategy.hpp"

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

  cout << sb_graph << endl;
  cout << "sb graph created!" << endl;

  auto partitions = best_initial_partition(sb_graph, *number_of_partitions);

  cout << partitions << endl;

  kl_sbg_partitioner(sb_graph, partitions);

  // just for debugging
  cout << endl;
  for (unsigned i = 0; i < partitions.size(); i++) {
    cout << i << " " << partitions[i] << endl;
  }
  cout << endl;

  sanity_check(sb_graph, partitions, *number_of_partitions);

  cout << "Exit code: " << ret << endl;

  return ret;
}
