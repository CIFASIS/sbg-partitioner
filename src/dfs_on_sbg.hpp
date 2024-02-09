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

#include <sbg/sbg.hpp>


namespace sbg_partitioner {

namespace search {

class DFS {

public:
    DFS(SBG::LIB::CanonSBG graph);

    void start();

    SBG::LIB::SetPiece current() const;

    bool next();

private:
    std::map<size_t, bool> _visited;
    std::map<size_t, std::set<size_t> > _adjacent;

    std::stack<size_t> _stack;
    std::set<size_t> _stack_mirror;

    size_t _root_node_idx;
    size_t _current_node_idx;

    SBG::LIB::CanonSBG _graph;

    void initialize_adjacents();

    void update_stack();
};

}
}