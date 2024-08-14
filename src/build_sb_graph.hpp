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

#include <sbg/sbg.hpp>

#include "weighted_sb_graph.hpp"

namespace sbg_partitioner {

/// Takes a path to a json file, reads it and then builds a sb graph with
/// a node for each access to a variable and an edge for each connection
/// between variables.
/// If a variable appears on the left and on the right side, an edge is created.
WeightedSBGraph build_sb_graph(const std::string& filename);


/// Ad hoc function to get pre image of an expression from its image.
/// @param image_interval  Image we want to get the pre image from
/// @param expression  Expression to get the pre image
/// @return pre image as an interval
SBG::LIB::Interval get_pre_image(const SBG::LIB::Interval& image_interval, const SBG::LIB::LExp& expression);


/// Takes a graph and a set of nodes and return a set of nodes connected with the mentioned set.
/// If a node is connected at least by one edge with one node in the set, it will be in the returned
/// set.
/// @param graph the graph where we are looking for connections.
/// @param node set of nodes we want to know its connections.
/// @return a set of nodes connected to the function parameter.
SBG::LIB::OrdSet get_adjacents(const SBG::LIB::CanonSBG& graph, const SBG::LIB::SetPiece& node);


unsigned get_node_size(SBG::LIB::SetPiece node);

unsigned get_node_size(const SBG::LIB::UnordSet& node);

std::pair<SBG::LIB::UnordSet, SBG::LIB::UnordSet> cut_interval_by_dimension(SBG::LIB::UnordSet& set_piece, size_t size);
}