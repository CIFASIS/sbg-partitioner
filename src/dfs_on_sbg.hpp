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
#include <memory>
#include <set>
#include <stack>
#include <vector>

#include <sbg/sbg.hpp>

#include "partition_strategy.hpp"


namespace sbg_partitioner {

namespace search {

class DFS {

public:
    typedef size_t node_identifier;

    /// pre_order: True means pre-order, False means post-order. In-order is not taken
    /// into account since the graph is not a binary tree.
    DFS(SBG::LIB::CanonSBG& graph, unsigned number_of_partitions, std::unique_ptr<PartitionStrategy> partition_strategy, bool pre_order);

    void start();

    void iterate();

    std::map<node_identifier, std::set<node_identifier> > adjacents() const;

    std::map<unsigned, std::set<SBG::LIB::SetPiece>> partitions() const;
private:
    unsigned _number_of_partitions;
    bool _pre_order;

    std::vector<node_identifier> _visited;
    std::vector<node_identifier> _partially_visited;
    std::map<node_identifier, std::set<node_identifier> > _adjacent;

    std::vector<node_identifier> _stack;

    size_t _root_node_idx;

    SBG::LIB::CanonSBG _graph;

    std::unique_ptr<PartitionStrategy> _partition_strategy;

    void initialize_adjacents();

    void fill_current_node_stack();

    void add_adjacent_nodes(const node_identifier id, const SBG::LIB::CanonMap& map, const SBG::LIB::SetPiece& edge);

    bool was_visited(node_identifier id);
    bool was_partially_visited(node_identifier id);
    bool already_added(node_identifier id);

    void add_it_partially(node_identifier id);
    void add_it_definitely(node_identifier id);

    void add_it_to_a_partition(node_identifier id);
};

}
}