#ifndef GRAPH_CLUSTERING_MCL
#define GRAPH_CLUSTERING_MCL

#include <iostream>
#include <fstream>
#include <string>
#include <map>
#include <stdlib.h>
#include <boost/algorithm/string.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/copy.hpp>

#include <boost/crc.hpp>

#ifdef OpenMP_SUPPORT
#include <omp.h> //OpenMP for shared memory parallelisation
#endif

using namespace std;



//define the probability data of random walks
template <typename VertexDescriptor>
struct vertex_randwalk {
    vertex_randwalk() : visited_vertices(0), rand_walk_prob(){}
    vertex_randwalk(const VertexDescriptor& v) : visited_vertices(1){rand_walk_prob[v]=1;}
    int visited_vertices;
    map<VertexDescriptor, double> rand_walk_prob;
};

template <typename rand_walk_data_t, typename v_desc_t>
void show_rand_walk_data(const rand_walk_data_t& randwalk_data, int nMaxVertices=200){
    //output randwalk_data
    typename rand_walk_data_t::const_iterator rwd_iter;
    int nVertices = 0;
    for(rwd_iter=randwalk_data.begin(); rwd_iter!=randwalk_data.end(); ++rwd_iter){
        //vertex_randwalk<v_desc_t> data;
        //cout << "vertex: " << rwd_iter->first << "  visited_ vertices: " <<  rwd_iter->second.visited_vertices;

        nVertices++;
        if(nVertices>nMaxVertices) return;

        cout << "vertex: " << rwd_iter->first;
        typename map<v_desc_t, double>::iterator vit;
        map<v_desc_t, double> rand_walk_cur = rwd_iter->second.rand_walk_prob;
        for(vit=rand_walk_cur.begin(); vit!=rand_walk_cur.end(); ++ vit) {
            cout << " (" << vit->first << "," << vit->second << ")";
        }
        cout << " ===> " << rand_walk_cur.size();
        cout << endl;
    }
}

template <typename Graph> void show_graph_weights(const Graph& g){
    typedef typename boost::graph_traits<Graph>::vertex_iterator v_iter_t;
    typedef typename boost::graph_traits<Graph>::out_edge_iterator e_iter_t;

    pair<v_iter_t, v_iter_t> vip;
    pair<e_iter_t, e_iter_t> eip;

    typename boost::property_map<Graph, boost::edge_weight_t>::const_type eWeightMap;
    eWeightMap = boost::get(boost::edge_weight, g);

    cout << "weight map" << endl;
    for(vip=boost::vertices(g); vip.first!=vip.second; ++vip.first){
        cout << *vip.first << " : ";
        for(eip = boost::out_edges(*vip.first, g); eip.first!=eip.second; ++eip.first){
            cout << "(" << boost::target(*eip.first, g) << "," << boost::get(eWeightMap, *eip.first) << ") ";
        }
        cout << endl;
    }
}

struct mcl_config {
    mcl_config() : max_iteration(100), inflation_power(2), verbose(0), epsilon(1E-5), min_mem(0) {}
    int max_iteration;
    double inflation_power; //inflation value
    bool min_mem; //minimum memory consumption
    int verbose; //verbose output
    double epsilon; //minimum value that is significantly small
};


//mcl algorithm
//input: graph, config
//output: clustering result (bidirectional bimap data structure)
//  id1->id1, id2->id1, id3->id1, ... id4->id4, ...
template <typename Graph, typename ClusteringResult> void mcl(const Graph& g1, ClusteringResult &clustering_result, const mcl_config& config){
    if(config.verbose) cout << "start clustering..." << endl;

    //Copy the original graph to a mutable one
    Graph g;
    boost::copy_graph(g1, g);

    typedef typename boost::graph_traits<Graph>::vertex_iterator v_iter_t;
    typedef typename boost::graph_traits<Graph>::out_edge_iterator e_iter_t;
    typedef typename boost::graph_traits<Graph>::vertex_descriptor v_desc_t;
    typedef typename boost::graph_traits<Graph>::edge_descriptor e_desc_t;

    typename boost::property_map<Graph, boost::edge_weight_t>::type eWeightMap;
    eWeightMap = boost::get(boost::edge_weight, g);

    //add selfloop and normalize the weight values
    pair<v_iter_t, v_iter_t> vip;
    pair<e_iter_t, e_iter_t> eip;
    for(vip=boost::vertices(g); vip.first!=vip.second; ++vip.first){
        //output weights of all edges
        double outWeights=0;
        double averageWeight=1;
        int out_degree=0;
        for(eip = boost::out_edges(*vip.first, g); eip.first!=eip.second; ++eip.first){
            outWeights+=boost::get(eWeightMap, *eip.first);
            out_degree++;
        }
        if(out_degree>0) averageWeight=outWeights/out_degree;
        //add self edge
        e_desc_t self_edge;
        bool isAdded;
        boost::tie(self_edge, isAdded) = boost::add_edge(*vip.first, *vip.first, g);
        if(isAdded) {
            boost::put(eWeightMap, self_edge, averageWeight);
            outWeights += averageWeight;
        }
        for(eip = boost::out_edges(*vip.first, g); eip.first!=eip.second; ++eip.first){
            double w = boost::get(eWeightMap, *eip.first);
            boost::put(eWeightMap, *eip.first, w/outWeights);
        }
    }

    //show weight map
    //show_graph_weights(g1);

    //define random walk data
    typedef map<v_desc_t, vertex_randwalk<v_desc_t> > randwalk_data_t;
    randwalk_data_t randwalk_data;
    randwalk_data_t *p_new_randwalk_data(NULL);


    //init random walk data
    for(vip=boost::vertices(g); vip.first!=vip.second; ++vip.first){
        vertex_randwalk<v_desc_t> data(*vip.first);
        randwalk_data[*vip.first] = data;
    }

    if(config.verbose) cout << "Initialization done" << endl;
    if(config.verbose>=100) show_rand_walk_data<randwalk_data_t, v_desc_t>(randwalk_data);

    //use boost CRC function to check if the algorithm is converged
    uint16_t prev_cksum = 0;

    //clustering core algorithm
    for(int ii=0; ii<config.max_iteration; ii++){
        //iterate all vertices

        //allocate memory for new random walk data
        if(config.min_mem){
            p_new_randwalk_data = &randwalk_data;
        } else {
            p_new_randwalk_data = new randwalk_data_t;
        }


        //collect vertex into a array so it can be parallized.
        //This is required because OpenMP does not support iterator parallization
        int num_vertices=boost::num_vertices(g);
        v_desc_t vertex_array[num_vertices];

        int jj=0;
        for(vip=boost::vertices(g); vip.first!=vip.second; ++vip.first){
            vertex_array[jj++] =  *vip.first;
            //init the new_randwalk_data, otherwise, the parallel code will panic
            (*p_new_randwalk_data)[*vip.first];
        }

#pragma omp parallel for private(jj)
        for(jj=0; jj<num_vertices; jj++){
            //iterate all vertices
            v_desc_t c_vertex = vertex_array[jj];

            vertex_randwalk<v_desc_t> &new_r_data((*p_new_randwalk_data)[c_vertex]);


            map<v_desc_t, double> t_prob;
            typename map<v_desc_t, double>::iterator vit;

            pair<e_iter_t, e_iter_t> eip;
            v_desc_t v_source = c_vertex;
            for(eip=boost::out_edges(v_source, g); eip.first!=eip.second; ++eip.first){
                //iterate all neighbors of v_source
                v_desc_t v_neighbor = boost::target(*eip.first, g);
                double weight_neighbor = boost::get(eWeightMap, *eip.first);

                //iterate transition matrix of neighbors
                vertex_randwalk<v_desc_t> &r_neightbor_data(randwalk_data[v_neighbor]);
                for(vit=r_neightbor_data.rand_walk_prob.begin(); vit!=r_neightbor_data.rand_walk_prob.end(); ++vit){
                    t_prob[vit->first] = t_prob[vit->first] + weight_neighbor*vit->second;
                }
            }

            //add small pertubation to prevent tie situation
            for(vit=t_prob.begin(); vit!=t_prob.end(); ++vit){
                vit->second += double(rand()) / RAND_MAX * 1E-6;
            }
            new_r_data.rand_walk_prob = t_prob;
        }


        //store the checksum of each vertex
        uint16_t vertices_crc[num_vertices];

        //do inflation and normalisation after each interation
#pragma omp parallel for private(jj)
        for(int jj=0; jj<num_vertices; jj++){
            //iterate all vertices
            v_desc_t c_vertex = vertex_array[jj];
            vertex_randwalk<v_desc_t> &r_data((*p_new_randwalk_data)[c_vertex]);
            typename map<v_desc_t, double>::iterator vit;

            //inflation
            for(vit=r_data.rand_walk_prob.begin(); vit!=r_data.rand_walk_prob.end(); ++vit){
                vit->second = pow(vit->second, config.inflation_power);
            }

            //remove small values
            vit=r_data.rand_walk_prob.begin();
            double sumNorm = 0;
            while(vit!=r_data.rand_walk_prob.end()){
                if(vit->second<config.epsilon){
                    r_data.rand_walk_prob.erase(vit++);
                }else{
                    sumNorm += vit->second;
                    ++vit;
                }
            }

            boost::crc_16_type crc_processor;
            //normalize
            for(vit=r_data.rand_walk_prob.begin(); vit!=r_data.rand_walk_prob.end(); ++vit){
                vit->second = vit->second / sumNorm;
                crc_processor.process_bytes(&(vit->first), sizeof(vit->first));
                uint round_prob = uint(vit->second*1E3);
                crc_processor.process_bytes(&(round_prob), sizeof(round_prob));
            }
            vertices_crc[jj]=crc_processor.checksum();
        }

        //update the random walk data
        if(!config.min_mem){
            randwalk_data = *p_new_randwalk_data;
            delete p_new_randwalk_data;
            p_new_randwalk_data = NULL;
        }

        boost::crc_16_type crc_processor;
        for(int jj=0; jj<num_vertices; jj++){
          crc_processor.process_bytes(&vertices_crc[jj], sizeof(vertices_crc[jj]));
        }

        uint16_t new_cksum = crc_processor.checksum();
        if(config.verbose) cout << "iteration " << ii << endl;
        //show intermediate result
        if(config.verbose>=100){
            show_rand_walk_data<randwalk_data_t, v_desc_t>(randwalk_data);
            cout << "checksum = " << new_cksum << endl;
        }


        if (prev_cksum == new_cksum){
            if(config.verbose) cout << "==> Converged in " << ii << " steps!" << endl;
            break;
        }

        prev_cksum = new_cksum;

    }

    //collect clustering result
    //iterate all vertices in the graph g
    for(vip=boost::vertices(g); vip.first!=vip.second; ++vip.first){
        v_desc_t src_v=*vip.first;
        //find the attractor: src_v -> dest_v
        const vertex_randwalk<v_desc_t> &v_data(randwalk_data[src_v]);
        typename map<v_desc_t, double>::const_iterator vit;
        v_desc_t dest_v(src_v);
        double max_value=0;
        for(vit=v_data.rand_walk_prob.begin(); vit!=v_data.rand_walk_prob.end(); ++vit){
            if(vit->second>max_value){
                max_value=vit->second;
                dest_v=vit->first;
            }
        }

        clustering_result.insert(typename ClusteringResult::value_type(src_v,dest_v));
    }

    if(config.verbose) cout << "clustering done" << endl;
}


#endif //GRAPH_CLUSTERING_MCL
