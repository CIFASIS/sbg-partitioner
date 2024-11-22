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

#include "partition_graph.hpp"

namespace sbg_partitioner {

namespace metrics {


struct communication_metrics {
    int edge_cut;
    int comm_volume;
    int max_comm_volume;
    float maximum_imbalance;
};


int edge_cut(const PartitionMap& partitions, const WeightedSBGraph& sb_graph);

std::pair<int, int> communication_volume(const PartitionMap& partitions, const WeightedSBGraph& sb_graph);

float maximum_imbalance(const PartitionMap& partitions, const WeightedSBGraph& sb_graph);

PartitionMap read_partition_from_file(const std::string& name, const WeightedSBGraph& sb_graph);

std::ostream& operator<<(std::ostream& os, const communication_metrics& comm_metrics);
}

}