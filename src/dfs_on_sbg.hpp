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
#include "weighted_sb_graph.hpp"


namespace sbg_partitioner {

namespace search {

void initialize_partitioning(WeightedSBGraph& graph, unsigned number_of_partitions);

void add_strategy(PartitionStrategy& strategy, bool pre_order);

std::vector<std::map<unsigned, std::set<SBG::LIB::SetPiece>>> partitionate();

class DFS {

public:
    DFS() = default;

    /// pre_order: True means pre-order, False means post-order. In-order is not taken
    /// into account since the graph is not a binary tree.
    DFS(WeightedSBGraph& graph, unsigned number_of_partitions);

    DFS& operator= (const DFS&) = delete;   //deleted copy-assignment operator
    DFS(DFS&&) = default;
    DFS& operator= (DFS&&) = default;   //added move assignment operator

    ~DFS() = default;

    void start();

    void iterate();

    std::vector<std::map<unsigned, std::set<SBG::LIB::SetPiece>>> partitions() const;

    /// @note strategy object should live while this class does
    void add_partition_strategy(PartitionStrategy& strategy, bool pre_order);

private:
    typedef size_t node_identifier;

    unsigned _number_of_partitions;

    std::vector<node_identifier> _visited;
    std::vector<node_identifier> _partially_visited;
    std::map<node_identifier, std::set<node_identifier> > _adjacent;

    std::vector<node_identifier> _stack;

    size_t _root_node_idx;

    WeightedSBGraph _graph;

    std::vector<PartitionStrategy*> _partition_strategy_pre_order;
    std::vector<PartitionStrategy*> _partition_strategy_post_order;

    void initialize_adjacents();

    void fill_current_node_stack();

    void add_adjacent_nodes(const node_identifier id, const SBG::LIB::CanonPWMap& map, const SBG::LIB::OrdSet& edge);

    bool was_visited(node_identifier id);
    bool was_partially_visited(node_identifier id);
    bool already_added(node_identifier id);

    void add_it_partially(node_identifier id);
    void add_it_definitely(node_identifier id);

    void add_it_to_a_partition(node_identifier id, bool pre_order);
};

}
}