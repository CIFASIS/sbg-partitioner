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

#include "build_sb_graph.hpp"
#include "partition_graph.hpp"


using namespace std;

using namespace SBG::LIB;
using namespace SBG::Util;


PartitionGraph::PartitionGraph(SBG::LIB::CanonSBG& graph, int number_of_partitions)
    : _graph(graph),
    _partitions({})
{
    make_initial_partition(number_of_partitions);
}


void PartitionGraph::print_partition()
{
    cout << endl;
    for (const auto& [interval, proc] : _partitions) {
        cout << interval << " -- " << proc << endl;
    }
    cout << endl;
}


const SBG::LIB::CanonSBG& PartitionGraph::graph() const { return _graph; }


void PartitionGraph::set_partition(size_t interval_index, size_t partition_number)
{
    assert(interval_index < _graph.V().size() and "Invalid index!");
    assert(partition_number < _partitions.size() and "Invalid partition!");

    _partitions[_graph.V()[interval_index]] = partition_number;
}


void PartitionGraph::make_initial_partition(int number_of_partitions)
{
    size_t begin, end = 0;
    size_t min_size = _graph.V().size() / number_of_partitions;

    for (int i = 0; i < number_of_partitions; i++) {

        begin = i * min_size;
        end = size_t(i + 1) * min_size;

        for (size_t j = begin; j < end; j++) {
            _partitions[_graph.V()[j]] = i;
        }
    }

    for (auto j = end; j < _graph.V().size(); j++) {

        size_t partition = rand() % number_of_partitions;
        _partitions[_graph.V()[j]] = partition;
    }
}