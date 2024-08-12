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

struct GainObject {
    size_t i;
    size_t j;
    int gain;
    size_t size;
};

struct KLBipartResult {
    SBG::LIB::UnordSet A;
    SBG::LIB::UnordSet B;
    int gain;
};


using ec_ic = std::pair<SBG::LIB::UnordPWMDInter , SBG::LIB::UnordPWMDInter>;

struct GainObjectComparator {
    bool operator()(const GainObject& gain_1, const GainObject& gain_2) const
    {
        return gain_1.gain >= gain_2.gain;
    }
};

using CostMatrix = std::set<GainObject, GainObjectComparator>;


void kl_sbg_partitioner(const SBG::LIB::BaseSBG& graph, PartitionMap& partitions);

KLBipartResult kl_sbg_bipart(const SBG::LIB::BaseSBG& graph, SBG::LIB::UnordSet& partition_a, SBG::LIB::UnordSet& partition_b);

int kl_sbg(const SBG::LIB::BaseSBG& graph, SBG::LIB::UnordSet& partition_a, SBG::LIB::UnordSet& partition_b);

ec_ic compute_EC_IC(
    const SBG::LIB::UnordSet& partition,
    const SBG::LIB::UnordSet& nodes,
    const SBG::LIB::UnordSet& partition_2,
    const SBG::LIB::BasePWMap& departure_map,
    const SBG::LIB::BasePWMap& arrival_map);


std::ostream& operator<<(std::ostream& os, const GainObject& gain);

std::ostream& operator<<(std::ostream& os, const CostMatrix& cost_matrix);

std::ostream& operator<<(std::ostream& os, const KLBipartResult& result);

}; // namespace sbg_partitioner
