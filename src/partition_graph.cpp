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


/// This is an ad hoc function which gets the addition of the size of each interval
static unsigned get_multidim_interval_size(OrdPWMDInter& final_edges)
{
    unsigned acc = 0;
    for(size_t i = 0; i < final_edges.size(); i++) {
        for(size_t j = 0; j < final_edges[i].size(); j++) {
            acc += final_edges[i][j].end() - final_edges[i][j].begin() + 1;
        }
    }

    return acc;
}


/// Takes edges and maps and gets edges which connect "node" from departure_map to arrival_map
static OrdPWMDInter get_comunication_edges(const SetPiece& node, const OrdSet& partition_set, const CanonPWMap& departure_map, const CanonPWMap& arrival_map)
{
    auto edges_map_1 = preImage(OrdSet(node), departure_map);
    auto arrival_nodes_map_1 = image(edges_map_1, arrival_map);
    auto arrival_nodes_partition = intersection(arrival_nodes_map_1, partition_set);
    auto involved_edges_map_1 = preImage(arrival_nodes_partition, arrival_map);
    auto final_edges = intersection(involved_edges_map_1, edges_map_1);

    return final_edges;
}


/// Takes edges and maps and gets the number of edges.
static unsigned get_comunication_cost(const SetPiece& node, const OrdSet& partition_set, const CanonPWMap& departure_map, const CanonPWMap& arrival_map)
{
    auto final_edges_1 = get_comunication_edges(node, partition_set, departure_map, arrival_map);
    auto final_edges_2 = get_comunication_edges(node, partition_set, arrival_map, departure_map);
    std::cout << "final_edges_1 " << final_edges_1 << "\t final_edges_2 " << final_edges_2 << std::endl;

    auto final_edges = cup(final_edges_1, final_edges_2);

    std::cout << "final_edges " << final_edges << std::endl;

    unsigned acc = get_multidim_interval_size(final_edges);

    return acc;
}


unsigned int PartitionGraph::internal_cost(const SetPiece& node, unsigned partition)
{
    OrdSet partition_set = _partitions[partition];

    // We expect that node is part of the partition
    auto partition_intersection = intersection(node, partition_set);
    assert(partition_intersection == node);

    unsigned acc = get_comunication_cost(node, partition_set, _graph.map1(), _graph.map2());

    return acc;
}


unsigned int PartitionGraph::external_cost(const SBG::LIB::SetPiece& node, unsigned partition)
{
    OrdSet partition_set = _partitions[partition];

    // We expect that node is not part of the partition
    cout << "intersection between " << node << ", " << partition_set << endl;
    auto partition_intersection = intersection(node, partition_set);
    assert(isEmpty(partition_intersection));

    unsigned acc = get_comunication_cost(node, partition_set, _graph.map1(), _graph.map2());

    return acc;
}


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
            DFS dfs(_graph, number_of_partitions, make_unique<PartitionStrategyGreedy>(number_of_partitions, _graph), pre_order);
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


int PartitionGraph::gain(const pair<unsigned, SetPiece> node_a, const pair<unsigned, SetPiece> node_b)
{
    cout << "PartitionGraph::gain between " << node_a.second << ", " << node_b.second << endl;
    unsigned c_ab = get_comunication_cost(node_a.second, node_b.second, _graph.map1(), _graph.map2());

    unsigned external_a = external_cost(node_a.second, node_b.first);
    unsigned external_b = external_cost(node_b.second, node_a.first);

    unsigned internal_a = internal_cost(node_a.second, node_a.first);
    unsigned internal_b = internal_cost(node_b.second, node_b.first);

    int d_a = external_a - internal_a;
    int d_b = external_b - internal_b;
    cout << "Costs a " << external_a << ", " << internal_a << endl;
    cout << "Costs b " << external_b << ", " << internal_b << endl;
    cout << "comunnication between a & b is " << c_ab << endl;
    int gain = d_a + d_b - 2*c_ab;

    return gain;
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