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

#include "partition_graph.hpp"

namespace sbg_partitioner
{
// void kl_sbg();

// void kl_sbg_bipart();

using ec_ic = std::pair<SBG::LIB::OrdPWMDInter , SBG::LIB::OrdPWMDInter>;

ec_ic compute_EC_IC(
    const SBG::LIB::OrdSet& partition,
    const SBG::LIB::OrdSet& nodes,
    const SBG::LIB::CanonPWMap& departure_map,
    const SBG::LIB::CanonPWMap& arrival_map);


// en principio es un vector pero quien sabe
std::vector<ec_ic> compute_diff(
    const SBG::LIB::OrdPWMDInter& partitions,
    const SBG::LIB::CanonSBG& graph);

}; // namespace sbg_partitioner
