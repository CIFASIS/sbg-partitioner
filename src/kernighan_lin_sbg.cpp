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


ostream& operator<<(ostream& os, const GainObject& gain)
{
    os << "< Partitions: ("
       << gain.i
       << ", "
       << gain.j
       << "), gain: "
       << gain.gain
       << ", size: "
       << gain.size
       << ">";

    return os;
}


/// This is an ad hoc function which gets the addition of the size of each interval
static unsigned get_multidim_interval_size(OrdPWMDInter& final_edges)
{
    unsigned acc = 0;
    for(size_t i = 0; i < final_edges.size(); i++) {
        for(size_t j = 0; j < final_edges[i].size(); j++) {
            acc += final_edges[i][j].end() - final_edges[i][j].begin() + 1;
        }
    }

    return acc;
}

static size_t get_c_ab(
    const OrdSet& a, const OrdSet& b,
    const CanonPWMap& map_1,
    const CanonPWMap& map_2)
{
    auto f = [](auto& a, auto& b, auto& departure_map, auto& arrival_map) {
        auto d = preImage(a, departure_map);
        auto i = image(d, arrival_map);
        auto inters = intersection(b, i);
        return inters;
    };

    auto intersection1 = f(a, b, map_1, map_2);
    auto intersection2 = f(a, b, map_2, map_1);

    auto comm_nodes = cup(intersection1, intersection2);

    size_t comm_size = get_multidim_interval_size(comm_nodes);

    return comm_size;
}

ec_ic compute_EC_IC(
    const OrdSet& partition,
    const OrdSet& nodes,
    const CanonPWMap& departure_map,
    const CanonPWMap& arrival_map)
{
    auto d = preImage(nodes, departure_map);
    auto i = image(d, arrival_map);
    auto ic = intersection(partition, i);
    auto ec = difference(i, ic);

    return make_pair(ec, ic);
}


vector<ec_ic> compute_diff(
    const OrdSet& partition,
    const CanonSBG& graph)
{
    vector<ec_ic> d;
    d.reserve(partition.pieces().size());

    for (const auto& p : partition.pieces()) {
        ec_ic comunication_nodes = compute_EC_IC(partition, p, graph.map1(), graph.map2());
        ec_ic comunication_nodes_2 = compute_EC_IC(partition, p, graph.map2(), graph.map1());
        comunication_nodes.first = cup(comunication_nodes.first, comunication_nodes_2.first);
        comunication_nodes.second = cup(comunication_nodes.second, comunication_nodes_2.second);
        d.push_back(comunication_nodes);
    }

    return d;
}

CostMatrix generate_gain_matrix(
    OrdSet& partition_a,
    OrdSet& partition_b,
    const CanonSBG& graph)
{
    CostMatrix cost_matrix;

    for (size_t i = 0; i < partition_a.pieces().size(); i++) {
        for (size_t j = 0; j < partition_b.pieces().size(); j++) {
            auto a = partition_a[i];
            auto b = partition_b[j];
            size_t min_size = min(a[0].end() - a[0].begin(), b[0].end() - b[0].begin());
            // no problem here, a is just a copy of partition_a[i], same for b
            a[0].set_end(a[0].begin() + min_size);
            b[0].set_end(b[0].begin() + min_size);
            std::cout << "new ones: " << a << ", " << b << endl;
            std::cout << "original ones: " << partition_a[i] << ", " << partition_b[j] << endl;

            OrdSet ec_nodes_a_1, ic_nodes_a_1;
            std::tie(ec_nodes_a_1, ic_nodes_a_1) = compute_EC_IC(partition_a, a, graph.map1(), graph.map2());

            OrdSet ec_nodes_a_2, ic_nodes_a_2;
            std::tie(ec_nodes_a_2, ic_nodes_a_2) = compute_EC_IC(partition_a, a, graph.map2(), graph.map1());

            OrdSet ec_nodes_a, ic_nodes_a;
            ec_nodes_a = cup(ec_nodes_a_1, ec_nodes_a_2);
            ic_nodes_a = cup(ic_nodes_a_1, ic_nodes_a_2);

            size_t ec_a = get_multidim_interval_size(ec_nodes_a);
            size_t ic_a = get_multidim_interval_size(ic_nodes_a);
            int d_a = ec_a - ic_a;

            OrdSet ec_nodes_b_1, ic_nodes_b_1;
            std::tie(ec_nodes_b_1, ec_nodes_b_1) = compute_EC_IC(partition_b, b, graph.map1(), graph.map2());

            OrdSet ec_nodes_b_2, ic_nodes_b_2;
            std::tie(ec_nodes_b_2, ic_nodes_b_2) = compute_EC_IC(partition_b, b, graph.map2(), graph.map1());

            OrdSet ec_nodes_b, ic_nodes_b;
            ec_nodes_b = cup(ec_nodes_b_1, ec_nodes_b_2);
            ic_nodes_b = cup(ic_nodes_b_1, ic_nodes_b_2);

            size_t ec_b = get_multidim_interval_size(ec_nodes_b);
            size_t ic_b = get_multidim_interval_size(ic_nodes_b);
            int d_b = ec_b - ic_b;

            size_t c_ab = get_c_ab(a, b, graph.map1(), graph.map2());

            int gain = d_a + d_b - 2 * c_ab;

            const auto gain_obj = GainObject{i, j, gain, min_size};;
            cout << gain_obj << endl;
            cost_matrix.insert(gain_obj);
        }
    }

    return cost_matrix;
}

void compute_gain_matrix(const OrdSet& partition_a, const OrdSet& partition_b)
{

}

}