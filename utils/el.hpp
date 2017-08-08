#include <iostream>
#include <fstream>
#include <string>
#include <stdlib.h>
#include <boost/algorithm/string.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/lexical_cast.hpp>

using namespace std;

//global defnition


template <typename Graph> int read_el(
        const string& iFileName,
        Graph& g,
        map<string, typename boost::graph_traits<Graph>::vertex_descriptor> &vIdMap ){
  ifstream iFile(iFileName.c_str());
  if(!iFile.is_open()){
      cerr << "Can not open file: " << iFileName << endl;
      return -1;
  }

  //at the moment we treat file as undirected. Later it support directed file

  //property map
  typename boost::property_map <Graph, boost::vertex_name_t>::type vNameMap = boost::get(boost::vertex_name, g);
  typename boost::property_map <Graph, boost::vertex_index2_t>::type vIndex2Map = boost::get(boost::vertex_index2, g);
  typename boost::property_map <Graph, boost::vertex_potential_t>::type vPotentialMap = boost::get(boost::vertex_potential, g);
  typename boost::property_map <Graph, boost::edge_weight_t>::type eWeightMap = boost::get(boost::edge_weight, g);

  //node id
  //map<string, typename boost::graph_traits<Graph>::vertex_descriptor> vIdMap;

  for(string line; getline(iFile, line);){
      //skip empty lines and comment lines(start with #)
      boost::trim(line);
      if (boost::starts_with(line, "#")) continue;

      //split the body text into fields
      vector<string> fields;
      boost::split(fields, line, boost::is_any_of(" ,;\t"));

      if (fields.size()<2) {
          cout << "File format incorrect. " << endl;
          return 1;
      }

      string nodeId1 = fields[0];
      string nodeId2 = fields[1];
      typename boost::graph_traits<Graph>::vertex_descriptor v1, v2;

      if (!vIdMap.count(nodeId1)) {
          //add node to graph
          v1=boost::add_vertex(g);

          //insert the node into id map
          vIdMap[nodeId1] = v1;

          //Node label is same as node id
          //add vertex label
          string nodeLabel=nodeId1;
          boost::trim(nodeLabel);
          boost::put(vNameMap, v1, nodeLabel);

          boost::trim(nodeId1);
          //add vertex id
          boost::put(vIndex2Map, v1, nodeId1);

          //add vertex potential
          //todo: shall be set inital potential to be 0 or 1 or any other number
          boost::put(vPotentialMap, v1, 0);
      } else {
          v1 = vIdMap[nodeId1];
      }

      if (!vIdMap.count(nodeId2)) {
          //add node to graph
          v2=boost::add_vertex(g);

          //insert the node into id map
          vIdMap[nodeId2] = v2;

          //Node label is same as node id
          //add vertex label
          string nodeLabel=nodeId2;
          boost::trim(nodeLabel);
          boost::put(vNameMap, v2, nodeLabel);

          boost::trim(nodeId2);
          //add vertex id
          boost::put(vIndex2Map, v2, nodeId2);

          //add vertex potential
          //todo: shall be set inital potential to be 0 or 1 or any other number
          boost::put(vPotentialMap, v2, 0);
      } else {
          v2 = vIdMap[nodeId2];
      }

      //add edge
      float weight = 1;
      if (fields.size()>=3) weight = atof(fields[2].c_str());

      typename boost::graph_traits<Graph>::edge_descriptor e;
      bool isAdded=false;
      boost::tie(e, isAdded)=boost::add_edge(vIdMap[nodeId1], vIdMap[nodeId2], g);
      if (isAdded) {
          boost::put(eWeightMap, e, weight);
      }

      //add another edge to make it undirected graph
      boost::tie(e, isAdded)=boost::add_edge(vIdMap[nodeId2], vIdMap[nodeId1], g);
      if (isAdded) {
          boost::put(eWeightMap, e, weight);
      }
  }

  iFile.close();

  return 0;
}
