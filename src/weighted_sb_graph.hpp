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


struct WeightedSBGraph : public SBG::LIB::BaseSBG
{
public:
    WeightedSBGraph() = default;
    WeightedSBGraph(SBG::LIB::BaseSBG& graph) : SBG::LIB::BaseSBG(graph) {}
    WeightedSBGraph(SBG::LIB::BaseSBG&& graph) : SBG::LIB::BaseSBG(graph) {}

    void set_weights(std::map<SBG::LIB::UnordSet, unsigned>& weights) { _weights = std::move(weights); }

    std::map<SBG::LIB::UnordSet, unsigned> get_weights() const { return _weights; }

    void set_weight(const SBG::LIB::UnordSet& edge_set, unsigned weight)
    {
        _weights[edge_set] = weight;
    }

    unsigned get_weight(const SBG::LIB::UnordSet& edge_set) const { return _weights.at(edge_set); }

private:
    std::map<SBG::LIB::UnordSet, unsigned> _weights;
};

WeightedSBGraph addSEW(SBG::LIB::BasePWMap pw1, SBG::LIB::BasePWMap pw2, std::map<SBG::LIB::UnordSet, unsigned> weights, WeightedSBGraph g);

std::ostream& operator<<(std::ostream& os, const WeightedSBGraph& graph);

}