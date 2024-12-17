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

#include <filesystem>
#include <getopt.h>
#include <iostream>
#include <fstream>
#include <optional>
#include <string>

#include "build_sb_graph.hpp"
#include "kernighan_lin_partitioner.hpp"
#include "partition_metrics_api.hpp"


using namespace std;

using namespace sbg_partitioner;


typedef std::vector<std::string> stringvec;
 
struct path_leaf_string
{
    std::string operator()(const std::filesystem::directory_entry& entry) const
    {
        return entry.path().string();
    }
};

void read_directory(const std::string& name, stringvec& v)
{
    std::filesystem::path p(name);
    std::filesystem::directory_iterator start(p);
    std::filesystem::directory_iterator end;
    std::transform(start, end, std::back_inserter(v), path_leaf_string());
}

void usage()
{
    cout << "Usage sbg-partitioner" << endl;
    cout << endl;
    cout << "-h, --help       Display this information and exit" << endl;
    cout << "-v, --version    Display version information and exit" << endl;
    cout << "-f, --filename   Path to the input file, a json file that represents "
            "the model we want to partitionate." << endl;
    cout << "-p, --partitions Number of partitions." << endl;
    cout << "-d,              path to a directory with txt files that indicate the "
            "model partitioning using other algorithms." << endl;
    cout << "-e               Imbalance epsilon, a value between 0 and 1." << endl;
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
    unsigned iterations = 5;
    optional<unsigned> number_of_partitions = nullopt;
    optional<string> output_sb_graph = nullopt;
    optional<float> epsilon = 0.0;

    while (true) {

        static struct option long_options[] = {
            {"filename", required_argument, 0, 'f'},
            // {"directory", required_argument, 0, 'd'},
            {"partitions", required_argument, 0, 'p'},
            // {"output", required_argument, 0, 'o'},
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

    string partition_str;
    long int total_time = 0;
    for (unsigned i = 0; i < iterations; i++) {
        auto start = chrono::high_resolution_clock::now();
        partition_str = partitionate_nodes(*filename, *number_of_partitions, *epsilon, output_sb_graph);
        auto end = chrono::high_resolution_clock::now();
        total_time += chrono::duration_cast<chrono::milliseconds>(end - start).count();
    }

    long int avg_time = total_time / iterations;

    cout << "Average time after " << iterations << " executions is: " << avg_time << " milliseconds" << endl;

    cout << "Results: " << partition_str << endl;

    filesystem::path f_path = filesystem::path(*filename);

    auto output_file = filesystem::path(*filename).replace_extension(filesystem::path("")).string();
    output_file += "_";
    output_file += to_string(*number_of_partitions);
    output_file += "_exec_time.txt";
    ofstream output_stream;
    output_stream.open(output_file);
    output_stream << "Average time after " << iterations << " executions is: " << avg_time << " milliseconds" << endl;
    output_stream << "Results: " << partition_str << endl;
    output_stream.close();

    return 0;
}