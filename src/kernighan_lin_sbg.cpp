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

#include <sbg/sbg.hpp>
#include <vector>
#include <utility>

#include "kernighan_lin_sbg.hpp"


using namespace std;

using namespace SBG::LIB;

namespace sbg_partitioner {

ec_ic compute_EC_IC(
    const OrdSet& partition,
    const OrdSet& nodes,
    const SBG::LIB::CanonPWMap& departure_map,
    const SBG::LIB::CanonPWMap& arrival_map)
{
    auto d = preImage(nodes, departure_map);
    auto i = image(d, arrival_map);
    auto ic = SBG::LIB::intersection(partition, i);
    auto ec = SBG::LIB::difference(i, ic);

    return make_pair(ec, ic);
}

vector<ec_ic> compute_diff(
    const OrdPWMDInter& partitions,
    const CanonSBG& graph)
{
    vector<ec_ic> d;
    d.reserve(partitions.pieces().size());

    for (const auto& p : partitions.pieces()) {
        ec_ic comunication_nodes = compute_EC_IC(partitions, p, graph.map1(), graph.map2());
        ec_ic comunication_nodes_2 = compute_EC_IC(partitions, p, graph.map2(), graph.map1());
        comunication_nodes.first = cup(comunication_nodes.first, comunication_nodes_2.first);
        comunication_nodes.second = cup(comunication_nodes.second, comunication_nodes_2.second);
        d.push_back(comunication_nodes);
    }

    return d;
}

}