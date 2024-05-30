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
#include <rapidjson/document.h>
#include <rapidjson/pointer.h>
#include <rapidjson/istreamwrapper.h>
#include <vector>
#include <util/defs.hpp>

#include "build_sb_graph.hpp"
#include "logger.hpp"


using namespace rapidjson;
using namespace std;

using namespace SBG::LIB;
using namespace SBG::Util;


namespace sbg_partitioner {

struct Var {
    string id;
    vector<pair<INT, INT>> exps;
    vector<int> defs;
};


struct Node {
    int id;
    vector<pair<int, int>> intervals;
    vector<Var> rhs;
    vector<Var> lhs;

    Node(int id, vector<pair<int, int>>&& intervals, vector<Var>&& rhs, vector<Var>&& lhs)
      : id(id),
      intervals(intervals),
      rhs(rhs),
      lhs(lhs)
    {}
};


/// This funcion takes a json array and returns a list of parsed variable objects (Var)
static vector<Var> read_var_object(const rapidjson::Value& var_array)
{
    vector<Var> vars;
    // vars.reserve(var_array.GetArray().Size());
    for (const auto &value : var_array.GetArray()) {
        assert(value.HasMember("id") and value["id"].IsString());
        string id = value["id"].GetString();

        assert(value.HasMember("exp") and value["exp"].IsArray());
        auto expressions = value["exp"].GetArray();

        vector<pair<INT, INT>> exps;
        assert(expressions.Size() == 1 and "this is just for testing air conditioner example!!");
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
        vars.push_back(var);
    }

    return vars;
}


/// This funcion takes a json object and returns a list of parsed node objects (Node)
/// and its id as key/value.
static map<int, Node> create_node_objects_from_json(const Document& document)
{
  auto nodes_array = document["nodes"].GetArray();
  map<int, Node> nodes;
  for (const auto& node : nodes_array) {

    assert(node.HasMember("id") and node["id"].IsInt());
    unsigned id = node["id"].GetInt();

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

    Node node_element = Node{int(id), std::move(intervals), std::move(rhs), std::move(lhs)};

    nodes.insert({node_element.id, node_element});
  }

  return nodes;
}


/// Creates a set of nodes, taking into accout the offset of each one to avoid collisions.
static UnordSet create_set_of_nodes(const map<int, Node>& nodes, map<int, int>& node_offsets, int& max_value)
{
  // We start to build out set of intervals from 0
  int current_max = 0;
  UnordSet node_set;
  for (const auto& [id, node] : nodes) {

    cout << "Defining interval for node " << id << " ";

    // Define the interval and add it to the node set taking into account the current offset
    // Create an offset for each equation node. We want that each equation has its
    // own domain.
    for (const auto& node_interval : node.intervals) {
      int interval_begin = current_max;
      int interval_end = (node_interval.second - node_interval.first) + current_max;
      Interval interval(interval_begin, 1, interval_end);

      node_set.emplace_hint(node_set.end(), interval);
      cout << interval << endl;

      // set this node offset, the difference between the interval and the original one
      node_offsets[id] =  interval_begin - node_interval.first;

      // udpate max value
      current_max = interval_end + 1;
    }
  }

  // We save max value so edge domain will not collide with node domain
  max_value = current_max;

  return node_set;
}


static vector<Exp> read_left_vars(const Node& node)
{
  const auto& l_nodes = node.lhs;

  assert(l_nodes.size() > 0);

  vector<Exp> exps;
  for (const auto& l_node : l_nodes) {
    Var var = l_node;
    Exp exp;

    for (const auto& values : var.exps) {
      exp.emplaceBack(LExp(RAT(values.first, 1), RAT(values.second, 1)));
    }

    exps.push_back(exp);
  }

  return exps;
}


/// This function creates a map to connect to different variables
/// @param pre_image  Subset of the domain of the expression we want to connect
/// @param edge_domain  Domain of the map
/// @param var_exp  Original expression of the variable
static BaseMap create_map_interval(const UnordSet pre_image, const UnordSet& edge_domain, const Exp& var_exps, int set_offset)
{
  BaseMap map;

  map.set_dom(edge_domain);

  Exp map_exps;
  int i =0;
  for (const auto& var_exp : var_exps.exps()) {
    // If the slope is 0, we just return the expression.
    if (var_exp.slope() == 0) {
      cout << "Creating constant interval" << endl;
      LExp map_exp = var_exp;
      INT offset = var_exp.offset().numerator();
      offset += set_offset;
      map_exp.set_offset(RAT(offset, 1));
      map_exps.emplaceBack(map_exp);
      continue;
    }

    // Slope remains the same
    LExp map_exp;
    map_exp.set_slope(var_exp.slope());

    // We start from the original offset and substract the domain offset
    INT offset = var_exp.offset().numerator();
    offset += set_offset;
    offset = offset - minElem(edge_domain)[i];

    // Then we take the minimum element of the pre-image and add the offset
    INT min_elem = minElem(pre_image)[i];
    offset += min_elem;

    // We set the offset
    map_exp.set_offset(RAT(offset, 1));
    map_exps.emplaceBack(map_exp);

    i++;
  }

  // And create the map
  map.set_exp(map_exps);

  return map;
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


/// @brief Add documentation
/// @param nodes 
/// @param node_offsets 
/// @param max_value 
/// @return 
static tuple<UnordSet, BasePWMap, BasePWMap> create_graph_edges(
        const std::map<int, Node>& nodes,
        const map<int, int>& node_offsets,
        int& max_value)
{
  UnordSet edge_set;  // Our set of edges
  BasePWMap left_maps;  // Map object of one of the sides
  BasePWMap right_maps; // Map object of one of the other side

  // We iterate the set of nodes read from the input file
  for (const auto& [id, node] : nodes) {
    cout << "Looking for connections with " << id << endl;

    // Define the equation intervals (without offsets)
    UnordSet intervals;
    SetPiece interval_set_piece;
    assert(node.intervals.size() == 1);
    for (const auto& node_interval : node.intervals) {
      Interval interval(node_interval.first, 1, node_interval.second);
      interval_set_piece.emplaceBack(interval);
    }
    intervals.emplace(interval_set_piece);

    vector<Exp> this_node_exps = read_left_vars(node);

    // Now, iterate the right expresions to connect them to their definitions.
    // When we say left_sth we are talking about the variable defined on the left side,
    // Analogous for right_sth
    for (const Var &right_var : node.rhs) {
      Exp right_exps;
      for (const auto& exp_values : right_var.exps) {
        LExp right_exp = LExp(RAT(exp_values.first, 1), RAT(exp_values.second, 1));
        right_exps.emplaceBack(right_exp);
      }

      auto base_map = BaseMap(intervals, right_exps);
      auto right_node_image = image(BaseMap(intervals, right_exps));

      // Definitions of this variable are on defs field. We iterate it so we
      // can map them.
      for (int i : right_var.defs) {
        cout << "Is it connected to " << i << "?" << endl;
        auto left_node = nodes.at(i);

        UnordSet left_intervals;
        SetPiece left_intervals_set_piece;
        for (const auto& node_interval : left_node.intervals) {
          Interval left_interval(node_interval.first, 1, node_interval.second);
          left_intervals_set_piece.emplaceBack(left_interval);
        }
        left_intervals.emplace(left_intervals_set_piece);

        vector<Exp> left_exps = read_left_vars(left_node);
        for (const auto& left_exp : left_exps) {
          auto left_node_image = image(left_intervals, BaseMap(left_intervals, left_exp));

          // we want to see if the intersection of the images is not empty
          auto image_intersection = intersection(right_node_image, left_node_image);
          if (isEmpty(image_intersection)) {
            cout << "No, it is not" << endl;
            continue;
          }
          cout << "Yes, it is: " << image_intersection << endl;

          UnordSet edge_domain_set;
          auto image_intersection_set = image_intersection[0];
          for (size_t i = 0; i < image_intersection_set.size(); i++){
            // Edge domain will be from the current max value and will have the quantity as the intersection
            auto image_int = image_intersection_set[i];

            auto domain_offset = image_int.end() - image_int.begin();

            // domain for both edges
            Interval edge_domain(max_value, 1, max_value + domain_offset);
            edge_set.emplace_hint(edge_set.end(), edge_domain);

            // Update maximum value
            max_value = maxElem(edge_domain) + 1;

            edge_domain_set.emplace(edge_domain);
          }

          // Create and save lhs map
          auto pre_image_left = preImage(BaseMap(image_intersection, left_exp));
          auto left_map_interval = create_map_interval(pre_image_left, edge_domain_set, left_exp, node_offsets.at(left_node.id));
          left_maps.emplace(left_map_interval);

          // Create and save rhs map
          auto pre_image_right = preImage(BaseMap(image_intersection, this_node_exps[0]));
          auto right_map_interval = create_map_interval(pre_image_right, edge_domain_set, this_node_exps[0], node_offsets.at(id));
          right_maps.emplace(right_map_interval);
        }
      }
    }
  }

  return {edge_set, left_maps, right_maps};
}


/// @brief  Add documentation
/// @param nodes 
/// @return 
static BaseSBG create_sb_graph(const std::map<int, Node>& nodes)
{
  int max_value = 0;  // We track the max value, so we avoid domain collision between edges and nodes
  map<int, int> node_offsets;

  // Now, we create our set of nodes.
  UnordSet node_set = create_set_of_nodes(nodes, node_offsets, max_value);

  // Now, let's build a graph!
  BaseSBG graph; // This will be our graph

  /// Firstly, add each node piece to the graph.
  for (SetPiece mdi : node_set.pieces())
    graph = addSV(UnordSet(mdi), graph);

  // Then, create edges and maps.
  const auto [edge_set, left_maps, right_maps] = create_graph_edges(nodes, node_offsets, max_value);

  // now add those edges and maps to the graph
  graph = addSE(left_maps, right_maps, graph);

  return graph;
}


BaseSBG build_sb_graph(const string& filename)
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

  logger("build_sb_graph.txt", graph);

  return graph;
}


static unsigned add_adjacent_nodes(const CanonMap& incoming_map, const CanonMap& arrival_map, const SetPiece& node, OrdSet& adjacents)
{
  auto map_image = image(incoming_map.dom()[0], incoming_map.exp());
  auto node_map_intersection = intersection(map_image.intervals()[0], node);

  unsigned qty = 0;
  if (not isEmpty(node_map_intersection)) {

    auto pre_image = preImage(OrdSet(node_map_intersection), incoming_map);
    auto adjs = image(pre_image[0], arrival_map.exp());
    qty += adjs[0].end() - adjs[0].begin() + 1;
    adjacents.emplaceBack(adjs[0]);
  }

  return qty;
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
      cout << map2.dom() << ", " << map1.dom() << map2_minus_map1_dom << endl;
      if (not isEmpty(map2_minus_map1_dom)){
        CanonMap map2_ = CanonMap(map2_minus_map1_dom, map2.exp());

        acc += add_adjacent_nodes(map2_, map1, node, adjacents);
      }
  }

  cout << "Comunication is " << acc <<endl;
  return adjacents;
}

}