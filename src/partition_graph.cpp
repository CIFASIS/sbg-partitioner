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


map<unsigned, UnordSet> make_initial_partition(BaseSBG& graph, unsigned number_of_partitions, PartitionAlgorithm algorithm, bool pre_order)
{
    map<unsigned, UnordSet> partitions_set;
    switch (algorithm) {
        case PartitionAlgorithm::GREEDY:
        default:
            DFS dfs(graph, number_of_partitions, make_unique<PartitionStrategyGreedy>(number_of_partitions, graph), pre_order);
            dfs.start();
            dfs.iterate();

            map<unsigned, set<SetPiece>> partitions = dfs.partitions();

            for (auto& [id, set] : partitions) {
                SBG::LIB::UnordSet set_piece;
                for (auto& s : set) {
                    set_piece.emplace(s.intervals().front());
                }
                partitions_set[id] = set_piece;
            }
    }

    return partitions_set;
}


unordered_set<size_t> get_connectivity_set(
    BaseSBG& graph,
    map<unsigned, UnordSet>& partitions,
    size_t edge_index)
{
    unordered_set<size_t> connectivity_set = {};

    const size_t n = graph.map1().size();
    assert (edge_index < n and "Invalid index!");

    auto e1 = graph.map1()[edge_index];
    auto im1 = image(e1);

    auto e2 = graph.map2()[edge_index];
    auto im2 = image(e2);

    for (const auto& [proc, set_piece] : partitions) {
        if (not isEmpty(intersection(im1, set_piece))) {
            connectivity_set.insert(proc);
        }

        if (not isEmpty(intersection(im2, set_piece))) {
            connectivity_set.insert(proc);
        }
    }

    return connectivity_set;
}


size_t connectivity_set_cardinality(BaseSBG& graph,
    map<unsigned, UnordSet>& partitions,
    size_t edge_index)
{
    unordered_set<size_t> conn_set = get_connectivity_set(graph, partitions, edge_index);
    return conn_set.size();
}

}