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

#pragma once

#include <map>
#include <set>
#include <stack>
#include <vector>

#include <sbg/sbg.hpp>


namespace sbg_partitioner {

namespace search {

class DFS {

public:
    typedef size_t node_identifier;

    DFS(SBG::LIB::CanonSBG graph);

    void start();

    void iterate();
private:
    std::vector<node_identifier> _visited;
    std::vector<node_identifier> _partially_visited;
    std::map<node_identifier, std::set<node_identifier> > _adjacent;

    std::map<node_identifier, std::vector<node_identifier>> _stack;

    size_t _root_node_idx;

    SBG::LIB::CanonSBG _graph;

    void initialize_adjacents();

    void fill_current_node_stack();

    void add_adjacent_nodes(const node_identifier id, const SBG::LIB::CanonMap& map, const SBG::LIB::SetPiece& edge);

    bool was_visited(node_identifier id);
    bool was_partially_visited(node_identifier id);
    bool already_added(node_identifier id);

    void add_it_partially(node_identifier id);
    void add_it_definitely(node_identifier id);
};

}
}