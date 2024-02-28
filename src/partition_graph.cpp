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

#include <set>

#include "build_sb_graph.hpp"
#include "dfs_on_sbg.hpp"
#include "partition_graph.hpp"


using namespace std;

using namespace SBG::LIB;

using namespace sbg_partitioner::search;

namespace sbg_partitioner {


PartitionGraph::PartitionGraph(CanonSBG& graph, unsigned number_of_partitions, PartitionAlgorithm algorithm, bool pre_order)
    : _graph(graph),
    _partitions({})
{
    make_initial_partition(number_of_partitions, algorithm, pre_order);
}


const ::CanonSBG& PartitionGraph::graph() const { return _graph; }


map<unsigned, OrdSet> PartitionGraph::partitions() const { return _partitions; }


unordered_set<size_t> PartitionGraph::get_connectivity_set(size_t edge_index) const
{
    unordered_set<size_t> connectivity_set = {};

    const size_t n = _graph.map1().size();
    assert (edge_index < n and "Invalid index!");

    auto e1 = _graph.map1()[edge_index];
    auto im1 = image(e1);

    auto e2 = _graph.map2()[edge_index];
    auto im2 = image(e2);

    for (const auto& [proc, set_piece] : _partitions) {
        if (not isEmpty(intersection(im1, set_piece))) {
            connectivity_set.insert(proc);
        }

        if (not isEmpty(intersection(im2, set_piece))) {
            connectivity_set.insert(proc);
        }
    }

    return connectivity_set;
}


// I wish this was a separate function, not part of PartitionGraph but there were a lot of
// compile problems if partitions map object is created locally and OrdSet objects are added.
void PartitionGraph::make_initial_partition(unsigned number_of_partitions, PartitionAlgorithm algorithm, bool pre_order)
{
    switch (algorithm) {
        case PartitionAlgorithm::GREEDY:
        default:
            DFS dfs(_graph, number_of_partitions, pre_order);
            dfs.start();
            dfs.iterate();

            map<unsigned, set<SetPiece>> partitions = dfs.partitions();

            for (auto& [id, set] : partitions) {
                OrdSet set_piece;
                for (auto& s : set) {
                    set_piece.emplaceBack(s.intervals().front());
                }
                _partitions[id] = set_piece;
            }
    }
}


size_t connectivity_set_cardinality(const PartitionGraph& pgraph, size_t edge_index)
{
    unordered_set<size_t> conn_set = pgraph.get_connectivity_set(edge_index);
    return conn_set.size();
}


std::ostream& operator<<(std::ostream& os, const PartitionGraph& pgraph)
{
    os << "\n";
    for (const auto& [interval, proc] : pgraph.partitions()) {
        os << interval << " -- " << proc << "\n";
    }
    os << "\n";

    return os;
}

}