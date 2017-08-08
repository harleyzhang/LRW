#include <iostream>
#include <fstream>
#include <string>
#include <stdlib.h>
#include <boost/algorithm/string.hpp>
#include <boost/graph/adjacency_list.hpp>
//#include <boost/graph/graphviz.hpp>

using namespace std;

//global defnition
const char NodeDefHeader[]="nodedef>";
const char EdgeDefHeader[]="edgedef>";
const int HeaderSize=sizeof(NodeDefHeader)-1;


template <typename Graph> int read_gdf(
        const string& iFileName,
        Graph& g,
        map<string, typename boost::graph_traits<Graph>::vertex_descriptor> &vIdMap ){
  ifstream iFile(iFileName.c_str());
  if(!iFile.is_open()){
      cerr << "Can not open file: " << iFileName << endl;
      return -1;
  }


  //define if we are processing the node definition
  bool isNodePart=true;

  //define indexes variables
  int nodeNameIdx=0, nodeLabelIdx=0;
  int edgeNode1Idx=0, edgeNode2Idx=0, edgeWeigthIdx=0;

  //property map
  typename boost::property_map <Graph, boost::vertex_name_t>::type vNameMap = boost::get(boost::vertex_name, g);
  typename boost::property_map <Graph, boost::vertex_index2_t>::type vIndex2Map = boost::get(boost::vertex_index2, g);
  typename boost::property_map <Graph, boost::vertex_potential_t>::type vPotentialMap = boost::get(boost::vertex_potential, g);
  typename boost::property_map <Graph, boost::edge_weight_t>::type eWeightMap = boost::get(boost::edge_weight, g);

  //node id
  //map<string, typename boost::graph_traits<Graph>::vertex_descriptor> vIdMap;

  for(string line; getline(iFile, line);){
      //check "nodedef"
      string lead=line.substr(0,HeaderSize);
      boost::algorithm::to_lower(lead);

      if (lead == "nodedef>") {
          //check node definition. only support: name, label
          string body=line.substr(HeaderSize, line.length()-HeaderSize);

          //split the body text into fields
          vector<string> fields;
          boost::split(fields, body, boost::is_any_of(","));

          for(vector<string>::iterator it=fields.begin(); it!=fields.end(); ++it){
              string fieldString(*it);
              //boost::trim(fieldString); //trim white space is not necessary

              stringstream fieldStream(fieldString);
              string fieldName;

              //read the first word
              fieldStream >> fieldName;
              boost::algorithm::to_lower(fieldName);

              if(fieldName=="name") {
                  nodeNameIdx=it - fields.begin();
              } else if (fieldName=="label"){
                  nodeLabelIdx=it - fields.begin();
              }
          }
          isNodePart=true;

      }else if (lead=="edgedef>") {
          //check edge definition. only support: node1, node2, weight

          string body=line.substr(HeaderSize, line.length()-HeaderSize);
          //split the body text into fields
          vector<string> fields;
          boost::split(fields, body, boost::is_any_of(","));

          for(vector<string>::iterator it=fields.begin(); it!=fields.end(); ++it){
              string fieldString(*it);

              stringstream fieldStream(fieldString);
              string fieldName;

              //read the first word
              fieldStream >> fieldName;
              boost::algorithm::to_lower(fieldName);

              if(fieldName=="node1") {
                  edgeNode1Idx=it - fields.begin();
              } else if (fieldName=="node2"){
                  edgeNode2Idx=it - fields.begin();
              } else if (fieldName=="weight"){
                  edgeWeigthIdx=it - fields.begin();
              }
          }
          isNodePart=false;

      }else{
          if (isNodePart){
              //handle nodes
              //split the line
              vector<string> fields;
              boost::split(fields, line, boost::is_any_of(","));


              //add vertex label
              typename boost::graph_traits<Graph>::vertex_descriptor v;
              v=boost::add_vertex(g);
              string nodeLabel=fields[nodeLabelIdx];
              boost::trim(nodeLabel);
              boost::put(vNameMap, v, nodeLabel);

              //add vertex potential
              //todo: shall be set inital potential to be 0 or 1 or any other number
              boost::put(vPotentialMap, v, 0);

              //add to node mapping
              string nodeId = fields[nodeNameIdx];
              boost::trim(nodeId);
              vIdMap[nodeId] = v;

              //add vertex id
              boost::put(vIndex2Map, v, nodeId);

              //cout << fields[nodeLabelIdx] << "|" <<fields[nodeNameIdx] << endl;
          }else{
              //handle edges
              vector<string> fields;
              boost::split(fields, line, boost::is_any_of(","));
              //add edge
              typename boost::graph_traits<Graph>::edge_descriptor e;
              bool isAdded=false;
              string node1Id=fields[edgeNode1Idx], node2Id=fields[edgeNode2Idx];
              boost::trim(node1Id);
              boost::trim(node2Id);
              boost::tie(e, isAdded)=boost::add_edge(vIdMap[node1Id], vIdMap[node2Id], g);
              if (isAdded) {
                  boost::put(eWeightMap, e, atof(fields[edgeWeigthIdx].c_str()));
              }
              //cout << fields[edgeNode1Idx] << "|" <<fields[edgeNode2Idx] << "|" << fields[edgeWeigthIdx] << endl;
          }
      }
  }

  iFile.close();

  return 0;
}
