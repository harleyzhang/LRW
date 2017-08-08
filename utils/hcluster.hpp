#ifndef GRAPH_CLUSTERING_HCLUSTER
#define GRAPH_CLUSTERING_HCLUSTER

#include "global.hpp"
#include "tools.hpp"


using namespace std;

//hierarchical graph data structure
template <typename Graph>
struct hierarchy_graph{
    typedef Graph graph_type;

    //vmap: lower_vertex -> upper_vertex
    typedef boost::bimap<boost::bimaps::unordered_set_of<typename Graph::vertex_descriptor>, boost::bimaps::multiset_of<typename Graph::vertex_descriptor> > vmap_t;
    vector<vmap_t> vmap_list;
    vector<Graph> graph_list;

    typedef typename boost::graph_traits<Graph>::vertex_descriptor v_desc_t;
    typedef typename boost::graph_traits<Graph>::vertex_iterator v_iter_t;
    typedef typename boost::graph_traits<Graph>::edge_descriptor e_desc_t;
    typedef typename boost::graph_traits<Graph>::edge_iterator e_iter_t;
};

template<typename Graph>
void show_hierarchy_graph(const hierarchy_graph<Graph> &hgraph, int index=-1){
    //-1: show the whole hierarchy graph
    if(index==-1) index=hgraph.graph_list.size()-1;

    //0: return
    if(index==0) return;

    //get reference to upper and lower graph
    const Graph &upper_g(hgraph.graph_list[index]);
    const Graph &lower_g(hgraph.graph_list[index-1]);

    //get reference to vertices map
    const typename hierarchy_graph<Graph>::vmap_t &vmap(hgraph.vmap_list[index-1]);
    typedef typename hierarchy_graph<Graph>::vmap_t vmap_t;

    typedef typename hierarchy_graph<Graph>::v_desc_t v_desc_t;
    typedef typename hierarchy_graph<Graph>::v_iter_t v_iter_t;
    typedef typename hierarchy_graph<Graph>::e_desc_t e_desc_t;
    typedef typename hierarchy_graph<Graph>::e_iter_t e_iter_t;

    //get vertex name map
    typename boost::property_map<Graph, boost::vertex_name_t>::const_type vUpperNameMap, vLowerNameMap;
    vUpperNameMap = boost::get(boost::vertex_name, upper_g);
    vLowerNameMap = boost::get(boost::vertex_name, lower_g);

    //iterate all vertices in upper graph
    for(pair<v_iter_t, v_iter_t> vp=boost::vertices(upper_g); vp.first!=vp.second; ++vp.first){
        //show vertex name of upper node
        cout << "(" << boost::get(vUpperNameMap, *vp.first) << ")" << " :";

        //show vertices name in lower graph
        const typename vmap_t::right_map &right = vmap.right;
        typename vmap_t::right_const_iterator r_it = right.find(*vp.first);
        typename vmap_t::right_const_iterator r_end = right.upper_bound(r_it->first);
        for(; r_it!=r_end; ++r_it){
            cout << " " << boost::get(vLowerNameMap, r_it->second);
        }
        cout << endl;
    }

    cout << "------------------" << endl;
    show_hierarchy_graph(hgraph, --index);
}

//hierarchical graph clustering
template<typename HClusteringResult, typename ClusteringConfig>
void hierarchical_clustering(HClusteringResult &hgraph, const ClusteringConfig config){
    //graph type definition
    typedef typename HClusteringResult::graph_type Graph;
    typedef typename HClusteringResult::v_desc_t v_desc_t;
    typedef typename HClusteringResult::v_iter_t v_iter_t;
    typedef typename HClusteringResult::e_desc_t e_desc_t;
    typedef typename HClusteringResult::e_iter_t e_iter_t;

    //get lower graph and do clustering
    const Graph &lower_g(hgraph.graph_list.back());

    //top clustering if all has been clustered
    if(boost::num_vertices(lower_g)==1) return;

    //clustering
    Clustering_Result_t clustering_rslt;
    //omp parallel support
    mcl(lower_g, clustering_rslt, config);
    if(config.verbose) show_clusteirng_result(clustering_rslt, lower_g);

    //Construct the upper layer graph and build hierarchical tree
    Graph upper_g;

    //mapping vertices of upper and lower graphs
    typename HClusteringResult::vmap_t vmap;

    //get reference to the upper and lower graph vertex property map
    typedef typename boost::property_map<Graph, boost::vertex_name_t>::type vname_map_t;
    typedef typename boost::property_map<Graph, boost::vertex_name_t>::const_type const_vname_map_t;
    vname_map_t upper_vname_map;
    const_vname_map_t lower_vname_map;
    upper_vname_map = boost::get(boost::vertex_name, upper_g);
    lower_vname_map = boost::get(boost::vertex_name, lower_g);

    //get lower and upper potential map
    typename boost::property_map<Graph, boost::vertex_potential_t>::const_type vLowerPotentialMap;
    typename boost::property_map<Graph, boost::vertex_potential_t>::type vUpperPotentialMap;
    vLowerPotentialMap = boost::get(boost::vertex_potential, lower_g);
    vUpperPotentialMap = boost::get(boost::vertex_potential, upper_g);

    //iterate all vertices in lower graph to construct upper graph
    const Clustering_Result_t::right_map &right = clustering_rslt.right;
    for(  typename Clustering_Result_t::right_const_iterator it = right.begin(), end = right.end();
          it != end;
          it = right.upper_bound(it->first))
    {
        typename Clustering_Result_t::right_const_iterator it2 = it;

        //add a vertex in the upper graph
        v_desc_t upper_v = boost::add_vertex(upper_g);

        //copy the destination vertex name to the upper layer
        //this is not a best solution, but a vertex has to have a name
        boost::put(upper_vname_map, upper_v, boost::get(lower_vname_map, it->first));

        //set intial potential to zero
        double upper_potential = 0;


        for(; it2!=right.upper_bound(it->first); ++it2){
            //add mapping from lower graph to upper graph
            //lower_id -> upper_id
            vmap.insert(typename HClusteringResult::vmap_t::value_type(it2->second,upper_v));
            upper_potential += boost::get(vLowerPotentialMap, it2->second);
        }

        //update upper graph vertex potential
        boost::put(vUpperPotentialMap, upper_v, upper_potential);
    }

    //collect edge weights for upper graph
    typedef map<pair<v_desc_t, v_desc_t>, double> upper_weight_map_t;
    upper_weight_map_t upper_weight_map;

    //collect vertex potentials for upper graph
    map<v_desc_t, double> upper_v_potential_map;


    //get lower weight map
    typename boost::property_map<Graph, boost::edge_weight_t>::const_type eWeightMap;
    eWeightMap = boost::get(boost::edge_weight, lower_g);

    //itterate all edges of the lower graph to calculate the accumulated weights
    std::pair<e_iter_t, e_iter_t> ep;
    for(ep=boost::edges(lower_g); ep.first!=ep.second; ++ep.first){
        double weight = boost::get(eWeightMap, *ep.first);
        v_desc_t upper_src, upper_dest;
        upper_src = vmap.left.at(boost::source(*ep.first, lower_g));
        upper_dest = vmap.left.at(boost::target(*ep.first, lower_g));
        //do not create self-loop in upper graph, instead we update the potentials
        if(upper_src == upper_dest) {
            //note: weights of each edges are added to the upper graph vertex potential.
            //if(upper_v_potential_map.find(upper_src) != upper_v_potential_map.end()){
                upper_v_potential_map[upper_src] += weight;
            //} else {
                //upper_v_potential_map[upper_src] = weight;
            //}
            continue;
        }
        pair<v_desc_t, v_desc_t> upper_e(upper_src, upper_dest);
        if(upper_weight_map.find(upper_e) != upper_weight_map.end()){
            //found the edge in the map
            upper_weight_map[upper_e] += weight;
        } else {
            upper_weight_map[upper_e] = weight;
        }
        //cout << boost::source(*ep.first, lower_g) << "-- " << weight << " -->" << boost::target(*ep.first, lower_g) << endl;
    }

    //add edges to upper graph
    typename boost::property_map<Graph, boost::edge_weight_t>::type eUpperWeightMap;
    eUpperWeightMap = boost::get(boost::edge_weight, upper_g);

    for(typename upper_weight_map_t::iterator uwm_iter=upper_weight_map.begin(); uwm_iter!=upper_weight_map.end(); ++uwm_iter){
        e_desc_t e;
        bool isAdded = false;
        boost::tie(e, isAdded)=boost::add_edge(uwm_iter->first.first, uwm_iter->first.second, upper_g);
        if (isAdded) {
            boost::put(eUpperWeightMap, e, uwm_iter->second);
        }
    }

    //update potentials of upper graph
    pair<v_iter_t, v_iter_t> vit;
    for(vit=boost::vertices(upper_g); vit.first != vit.second; ++vit.first){
        boost::put(vUpperPotentialMap, *vit.first, boost::get(vUpperPotentialMap, *vit.first)+upper_v_potential_map[*vit.first]);
    }

    if(config.verbose) { cout << "upper graph" << endl; print_graph(upper_g);}

    //note: this function will invalidate all iterators.
    hgraph.graph_list.push_back(upper_g);
    hgraph.vmap_list.push_back(vmap);

    //recursively clustering
    hierarchical_clustering(hgraph, config);
}

#endif
