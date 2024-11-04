/**
 This file is part of Set--Based Graph Library.

 SBG Library is free software: you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation, either version 3 of the License, or
 (at your option) any later version.

 SBG Library is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.

 You should have received a copy of the GNU General Public License
 along with SBG Library.  If not, see <http://www.gnu.org/licenses/>.

 ******************************************************************************/

#include <getopt.h>
#include <iostream>
#include <optional>
#include <string>

#include "kernighan_lin_partitioner.hpp"
#include "partition_metrics_api.hpp"


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
    optional<string> directory = nullopt;
    optional<unsigned> number_of_partitions = nullopt;
    optional<string> output_sb_graph = nullopt;
    optional<float> epsilon = 0.0;

    while (true) {

        static struct option long_options[] = {
            {"directory", required_argument, 0, 'd'},
            {"partitions", required_argument, 0, 'p'},
            // {"output", required_argument, 0, 'o'},
            {"version", no_argument, 0, 'v'},
            {"help", no_argument, 0, 'h'}
        };

        int option_index = 0;
        opt = getopt_long(argc, argv, "d:p:e:gvh:", long_options, &option_index);
        if (opt == EOF) break;

        switch (opt) {
        case 'd':
            if (optarg) {
                directory = string(optarg);
            }
            break;

        case 'p':
          if (optarg) {
            number_of_partitions = atoi(optarg);
          }
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

        default:
          abort();
        }
    }

    cout << *directory << ", " << *number_of_partitions << ", " << *epsilon << endl;

    const auto [wg, pm] = partitionate_nodes_for_metrics(*directory, *number_of_partitions, *epsilon);

    int edge_cut = metrics::edge_cut(pm, wg);

    auto [comm_volume, max_comm_volume] = metrics::communication_volume(pm, wg);

    auto max_imb = metrics::maximum_imbalance(pm, wg);

    cout << "Edge cut: " << edge_cut << " communication volume " << comm_volume << ", " << max_comm_volume << " maximum imbalance " << max_imb << endl;

    return 0;
}