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
#include <iostream>
#include <optional>
#include <string>

#include "kernighan_lin_partitioner.hpp"


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
  int opt;
  optional<string> filename = nullopt;
  optional<unsigned> number_of_partitions = nullopt;
  optional<string> output_sb_graph = nullopt;
  optional<float> epsilon = nullopt;

  while (true) {

    static struct option long_options[] = {
      {"filename", required_argument, 0, 'f'},
      {"partitions", required_argument, 0, 'p'},
      {"output", required_argument, 0, 'o'},
      {"version", no_argument, 0, 'v'},
      {"help", no_argument, 0, 'h'}
    };

    int option_index = 0;
    opt = getopt_long(argc, argv, "f:p:e:gvh:", long_options, &option_index);
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

    case 'g':
      output_sb_graph = "";
      break;

    case 'e':
    if (optarg) {
      epsilon = atof(optarg);
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

  if (not filename or not number_of_partitions) {
    usage();
    exit(1);
  }

  if (not epsilon) {
    epsilon = 0.0;
  }

  if (*epsilon < 0 or *epsilon > 1) {
    usage();
    exit(1);
  }

  cout << "filename is " << *filename << endl;
  cout << "number of partitions is " << *number_of_partitions << endl;

  auto partition_str = partitionate_nodes(*filename, *number_of_partitions, *epsilon, output_sb_graph);

  cout << "final results " << partition_str << endl;

  if (output_sb_graph) {
    cout << "output_sb_graph " << *output_sb_graph << endl;
  }

  return 0;
}
