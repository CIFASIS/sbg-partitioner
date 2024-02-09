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

using namespace std;

using namespace SBG::LIB;
using namespace SBG::Util;

using namespace sbg_partitioner::search;


DFS::DFS(CanonSBG graph) : _graph(graph)
{
    initialize_adjacents();
}

void DFS::initialize_adjacents()
{
    // Fill adjacents
    for (size_t i = 0; i < _graph.V().size(); i++) {
        _visited[i] = false;
        const auto incoming_node = _graph.V()[i];

        for (size_t edge_counter = 0; edge_counter < _graph.E().size(); edge_counter++) {
            const auto edge = _graph.E()[edge_counter];

            for (size_t map_counter = 0; map_counter < _graph.map1().size(); map_counter++) {
                const auto map = _graph.map1()[map_counter];
                const auto edge_map_intersection = intersection(edge.intervals()[0], map.dom());
                if (not isEmpty(edge_map_intersection)) {
                    const auto map_image = image(edge_map_intersection.pieces().begin()->intervals()[0], map.exp());

                    for (size_t node_idx = 0; node_idx < _graph.V().size(); node_idx++) {
                        if (node_idx == i) {
                            continue;
                        }

                        const auto potential_arriving_node = _graph.V()[node_idx];
                        if (intersection(potential_arriving_node, map_image) == potential_arriving_node) {
                            _adjacent[i].insert(node_idx);
                        }
                    }
                }
            }
        }
    }

    // Choosing root node
    _root_node_idx = 0;
    size_t current_max_adj_values = _adjacent[_root_node_idx].size();
    for (size_t i = 1; i < _graph.V().size(); i++) {
        if (_adjacent[i].size() > current_max_adj_values) {
            _root_node_idx = i;
            current_max_adj_values = _adjacent[_root_node_idx].size();
        }
    }

    cout << "Root node is " << _root_node_idx << endl;
}


void DFS::start()
{
    _stack = stack<size_t>();

    _current_node_idx = _root_node_idx;
    update_stack();

    for (size_t i = 0; i < _graph.V().size(); i++) {
        _visited[i] = false;
    }

    _visited[_current_node_idx] = true;
}


SetPiece DFS::current() const { return _graph.V()[_current_node_idx]; }


bool DFS::next()
{
    if (_stack.empty()) {
        for (size_t i = 0; i < _graph.V().size(); i++) {
            if (not _visited[i]) {
                _stack.push(i);
                _stack_mirror.insert(i);
                break;
            }
        }
    }

    if (_stack.empty()) {
        return false;
    }

    _current_node_idx = _stack.top();


    assert(not _visited[_current_node_idx]);
    _stack.pop();

    _visited[_current_node_idx] = true;

    update_stack();
    return true;
}

void DFS::update_stack()
{
    auto adjacents_node_idxs = _adjacent[_current_node_idx];
    for (size_t idx : adjacents_node_idxs) {
        if (not _visited[idx] and _stack_mirror.find(idx) == _stack_mirror.end()) {
            _stack.push(idx);
            _stack_mirror.insert(idx);
        }
    }
}

// static void actual_dfs(const size_t i, map<size_t, bool>& visited, const map<size_t, set<size_t> >& adjacent)
// {
//     cout << "Visiting " << i << std::endl;

//     visited[i] = true;

//     for (auto it = adjacent.at(i).cbegin(); it != adjacent.at(i).cend(); ++it) {
//         if (!visited[*it]) {
//             actual_dfs(*it, visited, adjacent);
//         }
//     }
// }

// static void dfs(const CanonSBG& graph)
// {
//     map<size_t, bool> visited;
//     map<size_t, set<size_t> > adjacent;



//     actual_dfs(4, visited, adjacent);
// }