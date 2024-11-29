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


#include <cassert>
#include <fstream>
#include <iostream>
#include <map>
#include <optional>
#include <rapidjson/document.h>
#include <rapidjson/pointer.h>
#include <rapidjson/istreamwrapper.h>
#include <vector>
#include <util/defs.hpp>
#include <util/logger.hpp>

#include "build_sb_graph.hpp"
#include "weighted_sb_graph.hpp"


using namespace rapidjson;
using namespace std;

using namespace SBG::LIB;
using namespace SBG::Util;


namespace sbg_partitioner {

// Using an unnamed namespace to define functions with internal linkage
namespace {

struct Var {
    string id;
    vector<pair<INT, INT>> exps;
    vector<int> defs;
    unsigned cost = 1;
};


std::ostream& operator<<(std::ostream& os, const Var& var)
{
  os << "id: \"" << var.id << "\"";

  os << ", exps: ";
  for (const auto& exp : var.exps) {
    os << "[" << exp.first << ", " << exp.second << "]";
  }

  os << ", defs:";
  for (const auto& v : var.defs) {
    os << " " << v;
  }

  os << " cost: " << var.cost;

  return os;
}


struct Node {
    int id;
    int weight;
    vector<pair<int, int>> intervals;
    vector<Var> rhs;
    vector<Var> lhs;

    Node(int id, int weight, vector<pair<int, int>>&& intervals, vector<Var>&& rhs, vector<Var>&& lhs)
      : id(id),
      weight(weight),
      intervals(intervals),
      rhs(rhs),
      lhs(lhs)
    {}
};


std::ostream& operator<<(std::ostream& os, const Node& node)
{
  os << "id: \"" << node.id << "\" " << endl;

  os << "weight: " << node.weight << endl;

  os << "intervals: ";
  for (const auto& interval : node.intervals) {
    os << "[" << interval.first << ", " << interval.second << "] ";
  }

  os << endl << "lhs: ";
  for (const auto& l : node.lhs) {
    os << l << "; ";
  }

  os << "\b " << endl << "rhs: ";
  for (const auto& r : node.rhs) {
    os << r << "; ";
  }
  os << "\b " << endl;

  return os;
}


/// This funcion takes a json array and returns a list of parsed variable objects (Var)
vector<Var> read_var_object(const rapidjson::Value& var_array)
{
    vector<Var> vars;
    // vars.reserve(var_array.GetArray().Size());
    for (const auto &value : var_array.GetArray()) {
        assert(value.HasMember("id") and value["id"].IsString());
        string id = value["id"].GetString();

        assert(value.HasMember("exp") and value["exp"].IsArray());
        auto expressions = value["exp"].GetArray();

        vector<pair<INT, INT>> exps;
        // assert(expressions.Size() == 1 and "this is just for testing air conditioner example!!");
        for (const auto& exp_array : expressions) {
          assert(exp_array.IsArray());
          const auto expression = exp_array.GetArray();

          assert(expression.Size() == 2 and "Size of expression object is not as expected");

          int exp_a = expression[0].GetInt();
          int exp_b = expression[1].GetInt();

          exps.push_back(make_pair(exp_a, exp_b));
        }

        auto def_object = value["defs"].GetArray();
        vector<int> defs;
        defs.reserve(def_object.Size());
        for(const auto& def : def_object) {
            defs.push_back(def.GetInt());
        }

        Var var = Var{id, exps, defs};

        if (value.HasMember("cost")) {
          unsigned weight = value["cost"].GetUint();
          var.cost = weight;
        }

        vars.push_back(var);
    }

    return vars;
}


/// This funcion takes a json object and returns a list of parsed node objects (Node)
/// and its id as key/value.
map<int, Node> create_node_objects_from_json(const Document& document)
{
  auto nodes_array = document["nodes"].GetArray();
  map<int, Node> nodes;
  for (const auto& node : nodes_array) {

    assert(node.HasMember("id") and node["id"].IsInt());
    unsigned id = node["id"].GetInt();

    int node_weight = 1; // default value
    if (node.HasMember("weight")) {
      assert(node["weight"].IsInt());
      node_weight = node["weight"].GetInt();
    }

    string error_msg = "Interval " + to_string(id) + " format is wrong";
    assert(node.HasMember("interval") and node["interval"].IsArray() and error_msg.c_str());

    auto node_intervals = node["interval"].GetArray();

    vector<pair<int, int>> intervals;
    intervals.reserve(node_intervals.Size());
    for (const auto& interval : node_intervals) {
      assert(interval.IsArray() and interval.Size() == 2);

      int interval_start = int(interval[0].GetInt());
      int interval_end = int(interval[1].GetInt());

      intervals.push_back(make_pair(interval_start, interval_end));
    }

    vector<Var> rhs = read_var_object(node["rhs"]);
    vector<Var> lhs = read_var_object(node["lhs"]);

    Node node_element = Node{int(id), node_weight, std::move(intervals), std::move(rhs), std::move(lhs)};

    nodes.insert({node_element.id, node_element});
  }

  for (const auto& [i, n]: nodes) {
    cout << n << endl;
  }
  cout << endl;

  return nodes;
}


/// Creates a set of nodes, taking into accout the offset of each one to avoid collisions.
tuple<OrdSet, NodeWeight> create_set_of_nodes(const map<int, Node>& nodes, map<int, int>& node_offsets, int& max_value)
{
  // We start to build out set of intervals from 0
  int current_max = 0;
  OrdSet node_set;
  NodeWeight weights;

  for (const auto& [id, node] : nodes) {

    cout << "Defining interval for node " << id << " ";

    // Define the interval and add it to the node set taking into account the current offset
    // Create an offset for each equation node. We want that each equation has its
    // own domain.
    SetPiece array_of_nodes;
    for (size_t i = 0; i < node.intervals.size(); i++) {
      auto node_interval = node.intervals[i];

      int interval_begin, interval_end;
      if ( i == 0) {
        interval_begin = current_max;
        interval_end = (node_interval.second - node_interval.first) + current_max;
      } else {
        interval_begin = node_interval.first;
        interval_end = node_interval.second;
      }
      Interval interval = Interval(interval_begin, 1, interval_end);

      array_of_nodes.emplaceBack(interval);
      cout << interval << endl;

      // set this node offset, the difference between the interval and the original one
      node_offsets[id] =  interval_begin - node_interval.first;

      // udpate max value
      if (i == 0) {
        current_max = interval_end + 1;
      }
    }
    node_set.emplace_hint(node_set.end(), array_of_nodes);
    weights.insert({array_of_nodes, node.weight});
  }

  // We save max value so edge domain will not collide with node domain
  max_value = current_max;

  return {node_set, weights};
}


vector<pair<Var, Exp>> read_left_vars(const Node& node, const string& var_id = "")
{
  const auto& l_nodes = node.lhs;

  assert(l_nodes.size() > 0);

  vector<pair<Var, Exp>> exps;
  for (const auto& l_node : l_nodes) {
    if (not var_id.empty() and l_node.id != var_id) {
      continue;
    }

    Var var = l_node;
    Exp exp;

    for (const auto& values : var.exps) {
      exp.emplaceBack(LExp(RAT(values.first, 1), RAT(values.second, 1)));
    }

    exps.push_back({var, exp});
  }

  return exps;
}


/// This function creates a map to connect to different variables
/// @param pre_image  Subset of the domain of the expression we want to connect
/// @param edge_domain  Domain of the map
/// @param var_exp  Original expression of the variable
CanonMap create_set_edge_map(const OrdSet& pre_image, const OrdSet& edge_domain, const Exp& var_exps, int set_offset)
{
  CanonMap map;

  map.set_dom(edge_domain);

  Exp map_exps;
  int i =0;
  for (const auto& var_exp : var_exps.exps()) {
    // If the slope is 0, we just return the expression.
    if (var_exp.slope() == 0) {
      cout << "Creating constant interval" << endl;
      LExp map_exp = var_exp;
      INT offset = var_exp.offset().numerator();
      if (i == 0) {
        offset += set_offset;
      }
      map_exp.set_offset(RAT(offset, 1));
      map_exps.emplaceBack(map_exp);
      i++;
      continue;
    }

    // Slope remains the same
    LExp map_exp;
    map_exp.set_slope(var_exp.slope());

    // We start from the original offset and substract the domain offset
    INT offset = var_exp.offset().numerator();
    if (i == 0) {
      offset += set_offset;
    }
    offset = offset - minElem(edge_domain)[i];

    // Then we take the minimum element of the pre-image and add the offset
    INT min_elem = minElem(pre_image)[i];
    offset += min_elem;

    // We set the offset
    map_exp.set_offset(RAT(offset, 1));
    map_exps.emplaceBack(map_exp);

    i++;
  }
  cout << "created " << i << " maps out of " << var_exps << endl;

  // And create the map
  map.set_exp(map_exps);

  // useful for debugging
  // cout << "edge_domain " << edge_domain << ": " << image(map) << endl;

  return map;
}


template<typename Set>
Set get_node_domain(Node node)
{
  // Domain of the node candidate
  Set node_intervals;
  SetPiece node_set_piece;
  for (const auto& node_interval : node.intervals) {
    Interval interval(node_interval.first, 1, node_interval.second);
    node_set_piece.emplaceBack(interval);
  }
  node_intervals.emplace(node_set_piece);

  return node_intervals;
}


template<typename S, typename T>
S get_edge_domain(SetPiece image_intersection_set, T& edge_set, int& max_value)
{
  cout << "get_edge_domain " << image_intersection_set << ", " << edge_set << endl;
  SetPiece edge_domain_set;
  for (size_t i = 0; i < image_intersection_set.size(); i++) {
    // Edge domain will be from the current max value and will have the quantity as the intersection
    auto image_int = image_intersection_set[i];

    auto domain_offset = image_int.end() - image_int.begin();

    Interval edge_domain;
    if (i == 0) {
      // domain for both edges
      edge_domain = Interval(max_value, 1, max_value + domain_offset);

      // Update maximum value
      max_value = maxElem(edge_domain) + 1;
    } else {
      edge_domain = Interval(image_int.begin(), 1, image_int.end());
    }

    edge_domain_set.emplaceBack(edge_domain);
  }
  edge_set.emplace_hint(edge_set.end(), edge_domain_set);

  return edge_domain_set;
}


tuple<OrdSet, CanonPWMap, CanonPWMap, EdgeCost> create_graph_edges(
        const std::map<int, Node>& nodes,
        const map<int, int>& node_offsets,
        int& max_value)
{
  OrdSet edge_set;  // Our set of edges
  CanonPWMap rhs_maps;  // Map object of one of the sides
  CanonPWMap lhs_maps; // Map object of one of the other side
  EdgeCost costs; // Weight of edges

  for (const auto& [id, node] : nodes) {
    cout << "Looking for connections with " << id << endl;

    // Define the equation intervals (without offsets)
    OrdSet current_node_domain;
    SetPiece interval_set_piece;
    // assert(node.intervals.size() == 1);
    for (const auto& node_interval : node.intervals) {
      Interval interval(node_interval.first, 1, node_interval.second);
      interval_set_piece.emplaceBack(interval);
    }
    current_node_domain.emplace(interval_set_piece);

    vector<pair<Var, Exp>> this_node_exps = read_left_vars(node);

    // Now, iterate the right hand side expresions to connect them to their definitions.
    for (const Var &right_var : node.rhs) {
      Exp right_exps;
      for (const auto& exp_values : right_var.exps) {
        LExp right_exp = LExp(RAT(exp_values.first, 1), RAT(exp_values.second, 1));
        right_exps.emplaceBack(right_exp);
      }

      // Now, calculate the image of the rhs expression
      auto rhs_map = CanonMap(current_node_domain, right_exps);
      // This is the image *used* by this expression, wwe want to
      // check where is defined.
      auto used_node_image = image(CanonMap(current_node_domain, right_exps));

      // Definitions of this variable are on defs field. We want to check if
      // intersects with any node
      for (int i : right_var.defs) {
        cout << "Is it connected to " << i << "?" << endl;
        auto node_candidate = nodes.at(i);

        // Domain of the node candidate
        OrdSet node_candidate_domain = get_node_domain<OrdSet>(node_candidate);

        // look for definitions of the same variable
        auto exps_and_var_names = read_left_vars(node_candidate, right_var.id);
        for (const auto& [var, node_candidate_exps] : exps_and_var_names) {
          if (var.id != right_var.id) {
            continue;
          }

          // Now, get the image.
          auto node_candidate_CanonMap = CanonMap(node_candidate_domain, node_candidate_exps);
          auto node_candidate_image = image(node_candidate_CanonMap);

          // we want to see if the intersection of the images is not empty
          auto candidate_image_intersection = intersection(node_candidate_image, used_node_image);
          if (isEmpty(candidate_image_intersection)) {
            cout << "No, it is not" << endl;
            continue;
          }
          cout << "Yes, it is: " << candidate_image_intersection << endl;

          // Now we need to create both maps, let's create their domain.
          auto image_intersection_set = candidate_image_intersection[0];

          // we need to create an edge for each left hand side variable, that means a couple of maps for each one
          const auto& [_, exp] = *this_node_exps.begin();

          if (node_candidate_exps.exps()[0].slope() == 0) {
            cout << "This should be 1-N " << node_candidate_domain << endl;
            auto node_size = node_candidate_domain[0][0].end() - node_candidate_domain[0][0].begin();
            image_intersection_set[0] = Interval(image_intersection_set[0].begin(), 1, image_intersection_set[0].begin() + node_size);

            OrdSet edge_domain_set = get_edge_domain<OrdSet>(image_intersection_set, edge_set, max_value);

            int offset = node_candidate_domain[0][0].begin() + node_offsets.at(i) - edge_domain_set[0][0].begin();
            cout << "node offset " << node_offsets.at(i) << ", " << edge_domain_set << " so offset is " << offset << endl;
            CanonMap to_node_candidate = CanonMap(edge_domain_set, LExp(1, RAT(offset, 1)));
            cout << "to_node_candidate " << to_node_candidate << endl;

            cout << to_node_candidate << endl;

            auto im = Interval(node_candidate_exps.exps()[0].offset().numerator(), 1, node_candidate_exps.exps()[0].offset().numerator());
            CanonMap to_current_node = create_set_edge_map(OrdSet(im), edge_domain_set, Exp(LExp(0, node_candidate_exps.exps()[0].offset())), node_offsets.at(id));
            cout << "to_current_node " << to_current_node << endl;

            lhs_maps.emplace(to_current_node);
            rhs_maps.emplace(to_node_candidate);

            continue;
          } else if (exp.exps()[0].slope() == 0) {
            cout << "This should be N-1" << endl;
            auto node_size = current_node_domain[0][0].end() - current_node_domain[0][0].begin();
            image_intersection_set[0] = Interval(image_intersection_set[0].begin(), 1, image_intersection_set[0].begin() + node_size);

            OrdSet edge_domain_set = get_edge_domain<OrdSet>(image_intersection_set, edge_set, max_value);

            int offset = current_node_domain[0][0].begin() + node_offsets.at(id) - edge_domain_set[0][0].begin();
            CanonMap to_current_node = CanonMap(edge_domain_set, LExp(1, RAT(offset, 1)));
            cout << "to_current_node " << to_current_node << endl;

            auto im = Interval(exp.exps()[0].offset().numerator(), 1, exp.exps()[0].offset().numerator());
            CanonMap to_node_candidate = create_set_edge_map(OrdSet(im), edge_domain_set, Exp(LExp(0, node_candidate_exps.exps()[0].offset())), node_offsets.at(i));
            cout << "to_node_candidate " << to_node_candidate << endl;

            lhs_maps.emplace(to_current_node);
            rhs_maps.emplace(to_node_candidate);

            continue;
          }

          OrdSet edge_domain_set = get_edge_domain<OrdSet>(image_intersection_set, edge_set, max_value);

          // Create map to node candidate
          auto node_candidate_map = create_set_edge_map(candidate_image_intersection, edge_domain_set, node_candidate_exps, node_offsets.at(i));
          auto node_candidate_map_image = image(node_candidate_map);
          cout << "map is " << node_candidate_map << endl;
          cout << "image: " << node_candidate_map_image << endl;

          // Create map to current node
          auto pre_image_current_node = preImage(OrdSet(image_intersection_set), rhs_map);
          auto im = image(CanonMap(OrdSet(pre_image_current_node), exp));
          auto current_node_map = create_set_edge_map(im, edge_domain_set, exp, node_offsets.at(id));
          auto current_node_map_image = image(current_node_map);
          cout << "map is " << current_node_map << endl;
          cout << "image: " << current_node_map_image << endl;

          if (not (current_node_map_image == node_candidate_map_image)) {
              lhs_maps.emplace(current_node_map);
              rhs_maps.emplace(node_candidate_map);
          }
          cout << "----" << endl;

          costs.insert({edge_domain_set, var.cost});
        }
      }
    }
  }

  return {edge_set, rhs_maps, lhs_maps, costs};
}

/// @brief  Add documentation
/// @param nodes 
/// @return 
WeightedSBGraph create_sb_graph(const std::map<int, Node>& nodes)
{
  int max_value = 0;  // We track the max value, so we avoid domain collision between edges and nodes
  map<int, int> node_offsets;

  // Now, we create our set of nodes.
  auto [node_set, weights] = create_set_of_nodes(nodes, node_offsets, max_value);
  cout << "node_set " << node_set << endl;

  // Now, let's build a graph!
  WeightedSBGraph graph; // This will be our graph

  // Firstly, add nodes to the graph.
  graph = addSVW(node_set, weights, graph);

  // Then, create edges and maps.
  auto [edge_set, left_maps, right_maps, costs] = create_graph_edges(nodes, node_offsets, max_value);

  // Now add those edges and maps to the graph
  graph = addSEW(left_maps, right_maps, costs, graph);

  return graph;
}


template<typename Set>
unsigned add_adjacent_nodes(const CanonPWMap& incoming_map, const CanonPWMap& arrival_map, const Set& node, OrdSet& adjacents)
{
  auto map_image = image(incoming_map);
  auto node_map_intersection = intersection(map_image, node);

  unsigned qty = 0;
  if (not isEmpty(node_map_intersection)) {

    auto pre_image = preImage(OrdSet(node_map_intersection), incoming_map);
    auto adjs = image(pre_image, arrival_map);
    qty += get_node_size(adjs, NodeWeight());
    for_each(adjs.begin(), adjs.end(), [&adjacents](auto& b) { adjacents.emplace(b); });
  }

  return qty;
}


pair<SetPiece, SetPiece> cut_interval(const SetPiece &interval, int cut_value)
{
    int interval_begin = interval.intervals().front().begin();
    int interval_end = interval.intervals().front().end();
    Interval interval_1(interval_begin, 1, cut_value);
    Interval interval_2(cut_value + 1, 1, interval_end);

    SetPiece set_1;
    set_1.emplaceBack(interval_1);

    SetPiece set_2;
    set_2.emplaceBack(interval_2);

    return make_pair(set_1, set_2);
}

}


/// Ad hoc function to get pre image of an expression from its image.
/// @param image_interval  Image we want to get the pre image from
/// @param expression  Expression to get the pre image
/// @return pre image as an interval
Interval get_pre_image(const Interval& image_interval, const LExp& expression)
{
  LExp inv_exp = inverse(expression);
  Interval pre_image = image(image_interval, inv_exp);
  return pre_image;
}


WeightedSBGraph build_sb_graph(const string& filename)
{
  cout << "Reading " << filename << "..." << endl;

  // Parse json document
  Document document;
  ifstream ifs(filename);
  IStreamWrapper isw(ifs);
  document.ParseStream(isw);

  // Now read the document and convert it into a known type
  auto nodes = create_node_objects_from_json(document);

  // Now, let's get our graph
  auto graph = create_sb_graph(nodes);

  SBG_LOG << graph;

  return graph;
}


OrdSet get_adjacents(const CanonSBG& graph, const OrdSet& node)
{
  OrdSet adjacents = {};

  add_adjacent_nodes(graph.map1(), graph.map2(), node, adjacents);
  add_adjacent_nodes(graph.map2(), graph.map1(), node, adjacents);

  return adjacents;
}

OrdSet get_adjacents(const CanonSBG& graph, const SetPiece& node)
{
  OrdSet adjacents = {};

  // Fill adjacents
  unsigned acc = 0;
  for (size_t map_counter = 0; map_counter < graph.map1().size(); map_counter++) {
      const auto map1 = graph.map1()[map_counter];
      const auto map2 = graph.map2()[map_counter];

      acc += add_adjacent_nodes(map1, map2, node, adjacents);

      auto map2_minus_map1_dom = difference(map2.dom(), map1.dom());
      if (not isEmpty(map2_minus_map1_dom)){
        CanonMap map2_ = CanonMap(map2_minus_map1_dom, map2.exp());

        acc += add_adjacent_nodes(map2_, map1, node, adjacents);
      }
  }

  return adjacents;
}


pair<OrdSet, OrdSet> cut_bidimensional_interval(const SetPiece &set_piece, size_t s)
{
    cout << "cutting interval " << set_piece << ", " << s << endl;

    auto size_node_2 = get_node_size(SetPiece(set_piece.intervals()[1]), NodeWeight());

    unsigned ammount_of_rows = s / size_node_2;

    unsigned rest = s % size_node_2;

    cout << "Ammount of rows " << ammount_of_rows << endl;

    OrdSet OrdSet_ret;
    SetPiece interval_2 = *set_piece.intervals().begin();

    if (ammount_of_rows > 0) {
        SetPiece interval_1;
        tie(interval_1, interval_2) = cut_interval(set_piece.intervals().front(), set_piece.intervals().front().begin() + ammount_of_rows - 1);
        cout << "Interval cut in " << ammount_of_rows << ": " << interval_1 << ", " << interval_2 << endl;
        if (interval_2.size() == 0 and rest > 0) {
            interval_2 = Interval(interval_1.intervals().front().end(), 1, interval_1.intervals().front().end());
        }
        interval_1.emplaceBack(set_piece.intervals()[1]);
        OrdSet_ret.emplace(interval_1);
    }

    if (rest > 0){
        SetPiece rest_set_piece;
        rest_set_piece.emplaceBack(Interval(interval_2.intervals().front().begin(), 1, interval_2.intervals().front().begin()));
        SetPiece interval_3, interval_4;
        tie(interval_3, interval_4) = cut_interval(set_piece.intervals()[1], set_piece.intervals()[1].begin() + rest - 1);
        rest_set_piece.emplaceBack(interval_3.intervals().front());

        OrdSet_ret.emplace(rest_set_piece);
    }

    OrdSet remaining = difference(OrdSet(set_piece), OrdSet_ret);

    cout << "original " << OrdSet(set_piece) << ", " << OrdSet_ret << ", " << remaining << endl;

    return make_pair(OrdSet_ret, remaining);
}


pair<OrdSet, OrdSet> cut_interval_by_dimension(OrdSet& set_piece, const NodeWeight& node_weight, size_t size)
{
    if (set_piece.pieces().size() == 0) {
        return make_pair(OrdSet(), OrdSet());
    }

    if (size == 0) {
        return make_pair(OrdSet(), set_piece);
    }

    size_t actual_size = size / unsigned(get_set_cost(*set_piece.pieces().begin(), node_weight));

    if (set_piece.pieces().begin()->intervals().size() == 1) {
        SetPiece p_1, p_2;
        tie(p_1, p_2) = cut_interval(set_piece.pieces().begin()->intervals().front(), set_piece.pieces().begin()->intervals().front().begin() + actual_size - 1);
        return make_pair(OrdSet(p_1), OrdSet(p_2));
    }

    if (set_piece.pieces().begin()->intervals().size() == 2) {
        OrdSet cut_node;
        OrdSet remaining_node = set_piece;
        for (const auto& piece : set_piece.pieces()) {
          size_t pice_size = get_node_size(piece, node_weight);
          size_t size_to_cut = min(pice_size, actual_size);
          OrdSet cut_node_piece, remaining_node_piece;
          tie(cut_node_piece, remaining_node_piece) = cut_bidimensional_interval(piece, size_to_cut);
          cut_node = cup(cut_node, cut_node_piece);
          remaining_node = difference(remaining_node, cut_node);

          if (remaining_node.pieces().empty()) {
            return make_pair(cut_node, remaining_node);
          }
        }

        return make_pair(cut_node, remaining_node);
    }

    cout << "Unexpected dimension: " << set_piece.pieces().begin()->intervals().size() << endl;

    throw 1;
}


unsigned get_node_size(const SetPiece& node, const NodeWeight& node_weight)
{
    if (node.size() == 0) {
        return 0;
    }

    int weight = get_set_cost(node, node_weight);

    unsigned acc = node.intervals().front().end() - node.intervals().front().begin() + 1;

    for (size_t i = 1; i < node.intervals().size(); i++) {
        auto interval = node.intervals()[i];
        acc = acc * (interval.end() - interval.begin() + 1);
    }

    acc *= weight;

    return acc;
}


unsigned get_node_size(const OrdSet& node, const NodeWeight& node_weight)
{
    if (node.pieces().empty()) {
      return 0;
  }

  unsigned size = 0;
  for (const auto& set_piece : node.pieces()) {
    size += get_node_size(set_piece, node_weight);
  }

  return size;
}


unsigned get_edge_set_cost(const SBG::LIB::SetPiece& node, const EdgeCost& edge_cost)
{
    if (node.size() == 0) {
        return 0;
    }

    int weight = get_set_cost(node, edge_cost);

    unsigned acc = node.intervals().front().end() - node.intervals().front().begin() + 1;

    for (size_t i = 1; i < node.intervals().size(); i++) {
        auto interval = node.intervals()[i];
        acc = acc * (interval.end() - interval.begin() + 1);
    }

    acc *= weight;

    return acc;
}


unsigned get_edge_set_cost(const SBG::LIB::OrdSet& node, const EdgeCost& edge_cost)
{
    if (node.pieces().empty()) {
        return 0;
    }

    unsigned size = 0;
    for (const auto& set_piece : node.pieces()) {
      size += get_edge_set_cost(set_piece,edge_cost);
    }

    return size;
}


void flatten_set(OrdSet &set, const CanonSBG& graph)
{
    if (set.pieces().size() > 0 and set.pieces().begin()->intervals().size() > 1) {
        cout << "flatten_set for sets with "
             << set.pieces().size()
             << " is not implemented"
             << endl;
        return;
    }

    OrdSet new_partition;
    for (const auto& v : graph.V()) {
        MDInterOrdSet set_piece_this_node_vector;
        for (auto& set_piece : set.pieces()) {

            if (not isEmpty(intersection(v, set_piece))) {
                set_piece_this_node_vector.emplace(set_piece);
            }
        }

        set_piece_this_node_vector = canonize(set_piece_this_node_vector);
        new_partition = cup(new_partition, set_piece_this_node_vector);
    }

    auto diff = difference(set, new_partition);
    assert(isEmpty(diff));

    set = new_partition;
}




WeightedSBGraph create_air_conditioners_graph()
{
    OrdSet nodes = OrdSet();
    nodes.emplaceBack(Interval(0, 1, 99)); //th
    nodes.emplaceBack(Interval(100, 1, 100)); //ierr
    nodes.emplaceBack(Interval(101, 1, 101)); //ptotal
    nodes.emplaceBack(Interval(102, 1, 103)); //ev_1
    nodes.emplaceBack(Interval(104, 1, 105)); //ev_2
    nodes.emplaceBack(Interval(106, 1, 106)); //ev_3
    nodes.emplaceBack(Interval(107, 1, 107)); //ev_4
    nodes.emplaceBack(Interval(108, 1, 108)); //ev_5
    nodes.emplaceBack(Interval(109, 1, 158)); //ev_6
    nodes.emplaceBack(Interval(159, 1, 208)); //ev_7
    nodes.emplaceBack(Interval(209, 1, 308)); //ev_8

    // maps
    CanonPWMap lhs_maps;
    CanonPWMap rhs_maps;

    lhs_maps.emplaceBack(CanonMap(Interval(319, 1, 368), Exp(LExp(1, RAT(-319, 1))))); // E1
    rhs_maps.emplaceBack(CanonMap(Interval(319, 1, 368), Exp(LExp(1, RAT(-210, 1)))));

    lhs_maps.emplaceBack(CanonMap(Interval(369, 1, 418), Exp(LExp(1, RAT(-319, 1))))); // E2
    rhs_maps.emplaceBack(CanonMap(Interval(369, 1, 418), Exp(LExp(1, RAT(-210, 1)))));

    lhs_maps.emplaceBack(CanonMap(Interval(419, 1, 518), Exp(LExp(1, RAT(-419, 1))))); // E3
    rhs_maps.emplaceBack(CanonMap(Interval(419, 1, 518), Exp(LExp(1, RAT(-210,1 )))));

    lhs_maps.emplaceBack(CanonMap(Interval(519, 1, 519), Exp(LExp(0, 100)))); // E4
    rhs_maps.emplaceBack(CanonMap(Interval(519, 1, 519), Exp(LExp(0, 106))));

    lhs_maps.emplaceBack(CanonMap(Interval(520, 1, 520), Exp(LExp(0, 100)))); // E5
    rhs_maps.emplaceBack(CanonMap(Interval(520, 1, 520), Exp(LExp(0, 107))));

    lhs_maps.emplaceBack(CanonMap(Interval(521, 1, 521), Exp(LExp(0, 100)))); // E6
    rhs_maps.emplaceBack(CanonMap(Interval(521, 1, 521), Exp(LExp(0, 108))));

    lhs_maps.emplaceBack(CanonMap(Interval(522, 1, 523), Exp(LExp(1, RAT(-420, 1))))); // E7
    rhs_maps.emplaceBack(CanonMap(Interval(522, 1, 523), Exp(LExp(1, RAT(-418, 1)))));

    lhs_maps.emplaceBack(CanonMap(Interval(524, 1, 525), Exp(LExp(1, RAT(-420, 1))))); // E8
    rhs_maps.emplaceBack(CanonMap(Interval(524, 1, 525), Exp(LExp(0, 101))));

    lhs_maps.emplaceBack(CanonMap(Interval(526, 1, 526), Exp(LExp(0, 106)))); // E9
    rhs_maps.emplaceBack(CanonMap(Interval(526, 1, 526), Exp(LExp(0, 108))));

    lhs_maps.emplaceBack(CanonMap(Interval(527, 1, 527), Exp(LExp(0, 107)))); // E10
    rhs_maps.emplaceBack(CanonMap(Interval(527, 1, 527), Exp(LExp(0, 108))));

    lhs_maps.emplaceBack(CanonMap(Interval(528, 1, 528), Exp(LExp(0, 101)))); // E11
    rhs_maps.emplaceBack(CanonMap(Interval(528, 1, 528), Exp(LExp(0, 108))));

    lhs_maps.emplaceBack(CanonMap(Interval(529, 1, 578), Exp(LExp(0, 102)))); // E12
    rhs_maps.emplaceBack(CanonMap(Interval(529, 1, 578), Exp(LExp(1, RAT(-420, 1)))));

    lhs_maps.emplaceBack(CanonMap(Interval(579, 1, 628), Exp(LExp(0, 103)))); // E13
    rhs_maps.emplaceBack(CanonMap(Interval(579, 1, 628), Exp(LExp(1, RAT(-420, 1)))));

    lhs_maps.emplaceBack(CanonMap(Interval(629, 1, 678), Exp(LExp(0, 108)))); // E14
    rhs_maps.emplaceBack(CanonMap(Interval(629, 1, 678), Exp(LExp(1, RAT(-520, 1)))));

    lhs_maps.emplaceBack(CanonMap(Interval(679, 1, 729), Exp(LExp(0, 108)))); // E15
    rhs_maps.emplaceBack(CanonMap(Interval(679, 1, 729), Exp(LExp(1, RAT(-520, 1)))));

      // Now, let's build a graph!
    WeightedSBGraph graph; // This will be our graph

    // Firstly, add nodes to the graph.
    graph = addSVW(nodes, NodeWeight(), graph);

    // Now add those edges and maps to the graph
    graph = addSEW(lhs_maps, rhs_maps, EdgeCost(), graph);

    return graph;
}

}