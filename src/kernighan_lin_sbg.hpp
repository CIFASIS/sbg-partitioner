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

#include <set>

#include "partition_graph.hpp"

namespace sbg_partitioner
{
// void kl_sbg();

// void kl_sbg_bipart();

struct GainObject {
    size_t i;
    size_t j;
    int gain;
    size_t size;
};

struct GainObjectComparator {
    bool operator()(const GainObject& gain_1, const GainObject& gain_2) const
    {
        return gain_1.gain >= gain_2.gain;
    }
};

using ec_ic = std::pair<SBG::LIB::OrdPWMDInter , SBG::LIB::OrdPWMDInter>;

using CostMatrix = std::set<GainObject, GainObjectComparator>;

std::ostream& operator<<(std::ostream& os, const GainObject& gain);

ec_ic compute_EC_IC(
    const SBG::LIB::OrdSet& partition,
    const SBG::LIB::OrdSet& nodes,
    const SBG::LIB::CanonPWMap& departure_map,
    const SBG::LIB::CanonPWMap& arrival_map);


// en principio es un vector pero quien sabe
std::vector<ec_ic> compute_diff(
    const SBG::LIB::OrdSet& partition,
    const SBG::LIB::CanonSBG& graph);


CostMatrix generate_gain_matrix(
    SBG::LIB::OrdSet& partition_a,
    SBG::LIB::OrdSet& partition_b,
    const SBG::LIB::CanonSBG& graph);


GainObject max_diff(
    CostMatrix& cost_matrix,
    SBG::LIB::OrdSet& partition_a,
    SBG::LIB::OrdSet& partition_b,
    const SBG::LIB::CanonSBG& graph);

}; // namespace sbg_partitioner
