/**
 This file is part of Set--Based Graph Library.

 SBG Library is free software: you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation, either version 3 of the License, or
 (at your option) any later version.

 SBG Library is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.

 You should have received a copy of the GNU General Public License
 along with SBG Library.  If not, see <http://www.gnu.org/licenses/>.

 ******************************************************************************/

#pragma once

#include <map>
#include <iostream>

#include <sbg/sbg.hpp>

namespace sbg_partitioner {

using EdgeCost = std::map<SBG::LIB::OrdSet, unsigned>;

using NodeWeight = std::map<SBG::LIB::OrdSet, int>;

struct WeightedSBGraph : public SBG::LIB::CanonSBG
{
public:
    WeightedSBGraph() = default;
    WeightedSBGraph(SBG::LIB::CanonSBG& graph) : SBG::LIB::CanonSBG(graph) {}
    WeightedSBGraph(SBG::LIB::CanonSBG&& graph) : SBG::LIB::CanonSBG(graph) {}

    void set_node_weights(NodeWeight& node_weights) { _node_weights = std::move(node_weights); }

    NodeWeight get_node_weights() const { return _node_weights; }

    void set_node_weight(const SBG::LIB::OrdSet& node_set, int weight) { _node_weights[node_set] = weight; }

    int get_node_weight(const SBG::LIB::OrdSet& node_set) const { return _node_weights.at(node_set); }


    void set_edge_costs(EdgeCost& edge_costs) { _edge_costs = std::move(edge_costs); }

    EdgeCost get_edge_costs() const { return _edge_costs; }

    void set_edge_cost(const SBG::LIB::OrdSet& edge_set, unsigned cost) { _edge_costs[edge_set] = cost; }

    unsigned get_edge_cost(const SBG::LIB::OrdSet& edge_set) const { return _edge_costs.at(edge_set); }

private:
    NodeWeight _node_weights;

    EdgeCost _edge_costs;
};

WeightedSBGraph addSVW(SBG::LIB::OrdSet nodes, NodeWeight weights, WeightedSBGraph g);

WeightedSBGraph addSEW(SBG::LIB::CanonPWMap pw1, SBG::LIB::CanonPWMap pw2, EdgeCost costs, WeightedSBGraph g);

std::ostream& operator<<(std::ostream& os, const WeightedSBGraph& graph);

}