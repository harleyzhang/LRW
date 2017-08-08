#ifndef GRAPH_CLUSTERING_TOOLS
#define GRAPH_CLUSTERING_TOOLS

#include <iostream>
#include <boost/graph/adjacency_list.hpp>
//#include "global.hpp"

//boost uBlas performace is bad. We use CSparse. 
extern "C" {
#include "csparse/csparse.h"
}

//print graph on screen
template <typename Graph>
void print_graph(const Graph& g){
    std::cout << "Graph:" << std::endl;
    //print vertices
    typedef typename boost::graph_traits<Graph>::vertex_iterator vertex_iter;
    std::pair<vertex_iter, vertex_iter> vp;
    std::cout << "number of vertices: " << boost::num_vertices(g) << std::endl;

    //property type
    typename boost::property_map<Graph, boost::vertex_name_t>::const_type vNameMap;
    vNameMap = boost::get(boost::vertex_name, g);

    //color map
    //typename boost::property_map<Graph, boost::vertex_color_t>::const_type vColorMap;
    //vColorMap = boost::get(boost::vertex_color, g);

    //vertex index map
    typename boost::property_map<Graph, boost::vertex_index_t>::const_type vIndexMap;
    vIndexMap = boost::get(boost::vertex_index, g);

    //vertex potential map
    typename boost::property_map<Graph, boost::vertex_potential_t>::const_type vPotentialMap;
    vPotentialMap = boost::get(boost::vertex_potential, g);

    //edge weight map
    typename boost::property_map<Graph, boost::edge_weight_t>::const_type eWeightMap;
    eWeightMap = boost::get(boost::edge_weight, g);

    for(vp=boost::vertices(g); vp.first!=vp.second; ++vp.first){
        //std::string vColor=boost::get(vColorMap, *vp.first);
        //int vIndex=boost::get(vIndexMap, *vp.first);
        //std::cout << *vp.first << " index: " << vIndex << "  color:  " << vColor << std::endl;
        std::cout << *vp.first << " : " << boost::get(vNameMap, *vp.first) << "  : " << boost::get(vPotentialMap, *vp.first) << std::endl;
    }

    //print edges
    typedef typename boost::graph_traits<Graph>::edge_iterator edge_iter;
    std::pair<edge_iter, edge_iter> ep;
    std::cout << "number of edges: " << boost::num_edges(g) << std::endl;
    for(ep=boost::edges(g); ep.first!=ep.second; ++ep.first){
        //std::cout << *ep.first << std::endl;
        //std::cout << boost::source(*ep.first, g) << "--" << boost::target(*ep.first, g) << std::endl;
        double weight = boost::get(eWeightMap, *ep.first);
        std::cout << boost::source(*ep.first, g) << "--" << boost::target(*ep.first, g) <<
                     "  (" << weight << ")  " << std::endl;
    }
}

//print clustering result on screen
template< typename ClusteringResult, typename Graph>
void show_clusteirng_result(const ClusteringResult &cr, const Graph &g) {
    const typename ClusteringResult::right_map &right = cr.right;

    typedef typename boost::property_map<Graph, boost::vertex_name_t>::const_type const_vname_map_t;
    const_vname_map_t vname_map;
    vname_map = boost::get(boost::vertex_name, g);


    //iterate all clusters
    for(  typename ClusteringResult::right_const_iterator it = right.begin(), end = right.end(); it != end; it = right.upper_bound(it->first))
      {
          typename ClusteringResult::right_const_iterator it2 = it;
          cout << it->first << " :";
          string tmp = boost::get(vname_map, it->first);
          cout << tmp << " : ";
          for(; it2!=right.upper_bound(it->first); ++it2){
              cout << " " << it2->second;
          }
          cout << endl;
      }
}

//return the adjacency matrix of graph. The adjacency matrix is the transpose of
//normal adjacency matrix. P_ij = P(x_n+1=i | x_n=j), which means j->i
//
//The ownership of the return value is transfered. If failed, NULL will be returned.
template<typename Graph>
cs_sparse* get_adjacency_matrix(const Graph &g, bool inTriplet=false){
    //get graph size
    int n = boost::num_vertices(g);
    int m = boost::num_edges(g);
    
    //define adjancecy matrix
    //sparse matrix, allocate memory, and triplet format
    cs_sparse *adjTriplet = cs_spalloc(n, n, m, 1, 1);

    //edge weight map
    typename boost::property_map<Graph, boost::edge_weight_t>::const_type eWeightMap;
    eWeightMap = boost::get(boost::edge_weight, g);

    //iterate all edges
    typedef typename boost::graph_traits<Graph>::edge_iterator edge_iter;

    std::pair<edge_iter, edge_iter> ep;
    for(ep=boost::edges(g); ep.first!=ep.second; ++ep.first){
        double weight = boost::get(eWeightMap, *ep.first);
        int ret = cs_entry(adjTriplet, boost::target(*ep.first, g), boost::source(*ep.first, g), weight);

        if (ret == 0) {
           return cs_spfree(adjTriplet);
        }
    }

    if (inTriplet) return adjTriplet;

    cs_sparse *adj = cs_triplet(adjTriplet);
    cs_spfree(adjTriplet);
    return adj;
}

#endif
