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


#include "weighted_sb_graph.hpp"


using namespace std;

using namespace SBG::LIB;


namespace sbg_partitioner {

WeightedSBGraph addSEW(BasePWMap pw1, BasePWMap pw2, std::map<UnordSet, unsigned> weights, WeightedSBGraph g)
{
    WeightedSBGraph graph = addSE(pw1, pw2, g);
    graph.set_weights(weights);

    return graph;
}


ostream& operator<<(ostream& os, const WeightedSBGraph& graph)
{
    os << BaseSBG(graph);

    os << "weights: ";
    for (const auto& [set, weight] : graph.get_weights()) {
        os << "{" << set << " : " << weight << "} ";
    }

    os << endl;

    return os;
}

}