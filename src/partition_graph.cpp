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

namespace sbg_partitioner {


PartitionGraph::PartitionGraph(::CanonSBG& graph, int number_of_partitions)
    : _graph(graph),
    _partitions({})
{
    make_initial_partition(number_of_partitions);
}


const ::CanonSBG& PartitionGraph::graph() const { return _graph; }


void PartitionGraph::set_partition(size_t interval_index, size_t partition_number)
{
    assert(interval_index < _graph.V().size() and "Invalid index!");
    assert(partition_number < _partitions.size() and "Invalid partition!");

    _partitions[_graph.V()[interval_index]] = partition_number;
}


map<SetPiece, int> PartitionGraph::partitions() const { return _partitions; }


unordered_set<size_t> PartitionGraph::get_connectivity_set(size_t edge_index) const
{
    unordered_set<size_t> connectivity_set = {};

    const size_t n = _graph.map1().size();
    assert (edge_index < n and "Invalid index!");

    auto e1 = _graph.map1()[edge_index];
    auto im1 = image(e1);

    auto e2 = _graph.map2()[edge_index];
    auto im2 = image(e2);

    for (const auto& [interval, proc] : _partitions) {
        if (not isEmpty(intersection(im1, interval))) {
            connectivity_set.insert(proc);
        }

        if (not isEmpty(intersection(im2, interval))) {
            connectivity_set.insert(proc);
        }
    }

    return connectivity_set;
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


size_t connectivity_set_cardinality(const PartitionGraph& pgraph, size_t edge_index)
{
    unordered_set<size_t> conn_set = pgraph.get_connectivity_set(edge_index);
    return conn_set.size();
}


std::ostream& operator<<(std::ostream& os, const PartitionGraph pgraph)
{
    os << "\n";
    for (const auto& [interval, proc] : pgraph.partitions()) {
        os << interval << " -- " << proc << "\n";
    }
    os << "\n";

    return os;
}

}