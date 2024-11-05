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

#include <iostream>

#include "dfs_on_sbg.hpp"
#include "weighted_sb_graph.hpp"

using namespace std;

using namespace SBG::LIB;
using namespace SBG::Util;


namespace sbg_partitioner {

namespace search {

static bool initialized = false;
static DFS sort_object;

void initialize_partitioning(
    WeightedSBGraph& graph,
    unsigned number_of_partitions)
{
    initialized = true;
    sort_object = DFS(graph, number_of_partitions);
}


void add_strategy(PartitionStrategy& strategy, bool pre_order)
{
    sort_object.add_partition_strategy(strategy, pre_order);
}


vector<map<unsigned, set<SetPiece>>> partitionate()
{
    assert(initialized and "Partioning was not inialized");

    sort_object.start();
    sort_object.iterate();
    vector<map<unsigned, set<SetPiece>>> partitions = sort_object.partitions();

    return partitions;
}


DFS::DFS(WeightedSBGraph & graph, unsigned number_of_partitions)
    : _number_of_partitions(number_of_partitions),
    _graph(graph)
{
    initialize_adjacents();
}


void DFS::initialize_adjacents()
{
    // Fill adjacents
    for (node_identifier i = 0; i < _graph.V().size(); i++) {
        const auto incoming_node = _graph.V()[i];

        for (node_identifier edge_counter = 0; edge_counter < _graph.E().size() - 1; edge_counter++) {
            const auto edge = _graph.E()[edge_counter];

            for (node_identifier map_counter = 0; map_counter < _graph.map1().size(); map_counter++) {
                const auto map1 = _graph.map1()[map_counter];
                add_adjacent_nodes(i, map1, edge);
                const auto map2 = _graph.map2()[map_counter];
                add_adjacent_nodes(i, map2, edge);
            }
        }
    }

    // Choosing root node
    _root_node_idx = 0;
    node_identifier current_max_adj_values = _adjacent[_root_node_idx].size();
    for (node_identifier i = 1; i < _graph.V().size(); i++) {
        if (_adjacent[i].size() > current_max_adj_values) {
            _root_node_idx = i;
            current_max_adj_values = _adjacent[_root_node_idx].size();
        }
    }

    cout << "Root node is " << _root_node_idx << endl;
}


void DFS::add_adjacent_nodes(const node_identifier id, const CanonMap& map, const SetPiece& edge)
{
    const auto edge_map_intersection = intersection(edge, map.dom());
    if (not isEmpty(edge_map_intersection)) {
        const auto map_image = image(*edge_map_intersection.pieces().begin(), map.exp());

        for (node_identifier node_idx = 0; node_idx < _graph.V().size(); node_idx++) {
            if (node_idx == id) {
                continue;
            }

            const auto potential_arriving_node = _graph.V()[node_idx];
            if (intersection(potential_arriving_node, map_image) == potential_arriving_node) {
                _adjacent[id].insert(node_idx);
            }
        }
    }
}


void DFS::start()
{
    // just in case, clear visited arrays
    _visited.clear();
    _partially_visited.clear();
    _stack.clear();

    // let's start with the root node
    add_it_partially(_root_node_idx);
    fill_current_node_stack();
}


void DFS::iterate()
{
    bool finished = false;
    while (not finished) {
        while (not _partially_visited.empty()) {
            node_identifier node_candidate = _partially_visited.back();
            while (not _stack.empty()) {
                const node_identifier new_node_candidate = _stack.back();
                _stack.pop_back();
                node_candidate = new_node_candidate;
                add_it_partially(node_candidate);
                fill_current_node_stack();
            }

            add_it_definitely(node_candidate);
        }

        // At this point, _partially_visited is empty, that means all nodes connected to
        // the root node were visited. Let's see if there is any isolated node.
        finished = _visited.size() == _graph.V().size();
        if (not finished) {
            for (node_identifier i = 0; i < _graph.V().size(); i++) {
                if (find(_visited.begin(), _visited.end(), i) == _visited.end()) {
                    add_it_partially(i);
                    break;
                }
            }
        }
    }
}


void DFS::add_partition_strategy(PartitionStrategy& strategy, bool pre_order)
{
    if (pre_order) {
        _partition_strategy_pre_order.push_back(&strategy);
    } else {
        _partition_strategy_post_order.push_back(&strategy);
    }
}


void DFS::fill_current_node_stack()
{
    node_identifier id = _partially_visited.back();
    for (const node_identifier adj_node : _adjacent[id]) {
        if (not was_partially_visited(adj_node) and not was_visited(adj_node)
            and not already_added(adj_node)) {
            _stack.push_back(adj_node);
        }
    }
}


vector<map<unsigned, set<SetPiece>>> DFS::partitions() const
{
    vector<map<unsigned, set<SetPiece>>> partitions;

    for (size_t i = 0; i < _partition_strategy_pre_order.size(); i++) {
        partitions.push_back(_partition_strategy_pre_order[i]->partitions());
    }

    for (size_t i = 0; i < _partition_strategy_post_order.size(); i++) {
        partitions.push_back(_partition_strategy_post_order[i]->partitions());
    }

    return partitions;
}


bool DFS::was_visited(node_identifier id) { return find(_visited.begin(), _visited.end(), id) != _visited.end(); }


bool DFS::was_partially_visited(node_identifier id)
{
    return find(_partially_visited.begin(), _partially_visited.end(), id) != _partially_visited.end();
}


bool DFS::already_added(node_identifier id)
{
    return (find(_stack.cbegin(), _stack.cend(), id) != _stack.cend());
}


void DFS::add_it_partially(node_identifier id)
{
    add_it_to_a_partition(id, true);

    _partially_visited.push_back(id);
}


void DFS::add_it_definitely(node_identifier id)
{
    add_it_to_a_partition(id, false);

    _partially_visited.pop_back();
    _visited.push_back(id);
}

void DFS::add_it_to_a_partition(node_identifier id, bool pre_order)
{
    auto v = _graph.V()[id];
    if (pre_order) {
        for (size_t i = 0; i < _partition_strategy_pre_order.size(); i++) {
            (*_partition_strategy_pre_order[i])(v);
        }
    } else {
        for (size_t i = 0; i < _partition_strategy_post_order.size(); i++) {
            (*_partition_strategy_post_order[i])(v);
        }
    }
}

}
}