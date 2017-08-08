#ifndef GRAPH_CLUSTERING_GLOBAL
#define GRAPH_CLUSTERING_GLOBAL


#include <boost/algorithm/string.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/copy.hpp>

#include <boost/bimap.hpp>
#include <boost/bimap/multiset_of.hpp>
#include <boost/bimap/unordered_set_of.hpp>

using namespace std;

//define properties
//index2 : nodeId in the original format: string
//name: node label, string
//potential: node potential, double
typedef boost::property< boost::vertex_index2_t, string, boost::property<boost::vertex_name_t, string,
    boost::property<boost::vertex_potential_t, double> > >vertex_p;
typedef boost::property<boost::edge_weight_t, double> edge_p;

//Graph type definition
//The default definition is a directed graph. In case of undirected graph, define two edges between each pair of
//connected nodes of the same weights.
typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::directedS, vertex_p, edge_p> Graph_t;

//clustering result type
typedef boost::bimap<boost::bimaps::unordered_set_of<Graph_t::vertex_descriptor>, boost::bimaps::multiset_of<Graph_t::vertex_descriptor> > Clustering_Result_t;

#endif
