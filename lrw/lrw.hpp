#ifndef GRAPH_CLUSTERING_LRW
#define GRAPH_CLUSTERING_LRW

#include <iostream>
#include <fstream>
#include <string>
#include <iterator>
#include <map>
#include <stdlib.h>

#include <boost/timer/timer.hpp>


#include <boost/algorithm/string.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/copy.hpp>
#include <boost/iterator/counting_iterator.hpp>
#include <boost/filesystem/path.hpp>

#include <boost/serialization/vector.hpp>
#include <boost/serialization/map.hpp>
#include <boost/serialization/set.hpp>

#include <boost/archive/text_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>

//boost uBlas performance is bad. we use CSparse instead.
#include "csparse/csparse.h"
auto csfree_fun = [](cs_sparse *A){cs_spfree(A);};

#include "tools.hpp"

#ifdef OpenMP_SUPPORT
#include <omp.h> //OpenMP for shared memory parallelisation
#endif

#ifdef MPI_SUPPORT
//#include <mpi.h>
#include <boost/mpi.hpp>
#include <boost/mpi/communicator.hpp>
#include <boost/mpi/status.hpp>
#endif

using namespace std;
//using namespace boost::numeric;

//helping function
#define vcout if(config.verbose) cout

//definition for master slave communication
const int C_WALK_COMPLETED = -1;


//lrw algorithm
struct lrw_config {
    lrw_config() :
        max_iteration(100),
        inflation_power(2),
        inflation_delay(0),
        neighbor_following(false),
        verbose(0),
        epsilon(1E-5),
        sig_threshold(0.5),
        int_threshold(0.5), 
        n_slave_work_load(50),
        n_max_explore_vertices(0),
        stage(0),
        seed(-1){}

    int seed; //if not -1, this is a local clustering algorithm that only find the cluster contain the given seed.
    int max_iteration; //max number of iteration
    double inflation_power; //inflation value
    int inflation_delay; //delay the inflation in couple of steps. Defaut is 0.
    double sig_threshold; //the significant threshold relative to the max value(0-1)
    double int_threshold; //the threshold for cluster merge, intersection above this will be merged(0-1)
    bool neighbor_following; //after clusters are found, refine each nodes to follow its major neighbors
    int verbose; //verbose output
    double epsilon; //minimum value that is significantly small
    string fname; //store the clustering figure file name
    string output_base; //output file name base(together with path)

    //programming configuration
    int n_slave_work_load; //The number of vertices that each slave process handle.
    int n_max_explore_vertices; //the maximum number of vertices that a walker reports

    //for extreme large network, split the work in stages. 0 indicate do as one shot
    int stage;
};

//type definition for clusters
typedef pair<int, double> vertex_prob_t;
typedef vector<vertex_prob_t> vertex_walk_rslt_t;
typedef vector< pair<int, vertex_walk_rslt_t> > walk_collection_t;

typedef vector<int> cluster_list_t;
typedef set < int > cluster_feature_t;
typedef map < int, pair<cluster_list_t, cluster_feature_t > >cluster_collector_t;

const string Collector_File_Suffix = ".explored.txt";

//helper class to serialze vertex walk results
namespace boost{
    namespace serialization{
        template <typename Archive> void serialize(Archive &ar, vertex_prob_t &a_pair, const unsigned version){
            ar & a_pair.first & a_pair.second;
        }

    /*
    template <typename Archive> void serialize(Archive &ar, pair<int, vertex_walk_rslt_t> &a_walk_vector, const unsigned version){
        ar & a_walk_vector.first & a_walk_vector.second;
    }*/

    }
}


// Explore one vertex, return the sorted result
// return empty vector: out of memory
//         otherwise: success
vertex_walk_rslt_t explore_vertex(const cs_sparse *adj, int vertex_id, const lrw_config &config) {
    vertex_walk_rslt_t walk_rslt;

    //random walk result vector
    //define as triplet format
    cs_sparse *x_triplet = cs_spalloc(adj->m, 1, 1, 1, 1);

    //init x_triplet
    int ret = cs_entry(x_triplet, vertex_id, 0, 1);
    if (ret != 1){
        cs_spfree(x_triplet);
        return(walk_rslt);
    }

    cs_sparse *x = cs_triplet(x_triplet);
    cs_spfree(x_triplet);
    //keep x safe
    unique_ptr<cs_sparse, decltype(csfree_fun)>pX(x, csfree_fun);

    //cout << "init completed, nnz: " << cs_nnz(pX.get()) << endl;


    //start walking iteration
    for (int ii=0; ii<config.max_iteration; ++ii){

        //cs_print(pX.get(), 0);
        //cout << "before walk" << endl << endl;

        //one step walk
        cs_sparse *x2 = cs_multiply(adj, pX.get());

        if (x2==NULL){
            //out of memory
            vcout << "Out of memory at iteration " << ii << "  size: " << cs_nnz(pX.get()) << endl;
            return(walk_rslt);
        }


        //inflate x
        if (ii>config.inflation_delay){
          cs_element_power(x2, config.inflation_power); //return value ignored
        }

        //normalise x
        cs_normalise_c(x2, 1); //ignore return value

        //drop small values
        cs_droptol(x2, config.epsilon);


        //check convergence
        cs_sparse *deltaX = cs_add(x2, pX.get(), 1, -1);
        if (deltaX==NULL){
            vcout << "Out of memory at iteration" << ii << " size:" << cs_nnz(pX.get()) << endl;
            cs_spfree(x2);
            return(walk_rslt);
        }

        double delta = cs_norm(deltaX);
        cs_spfree(deltaX);
        if (delta<config.epsilon){
            vcout << "iteration completed : " << ii << endl;
            cs_spfree(x2);
            break;
        }

        //cs_print(x2, 0);
        //cout << "Iteration: " << ii << " delta=" << delta << " nnz=" << cs_nnz(x2) << endl;

        //prevent oscillation, there can be other method to adapt to the step size
        if (ii>1) {
            cs_sparse *x3 = cs_add(pX.get(), x2, 0.5, 0.5);

            if (x3 == NULL) {
                //out of memory
                vcout << "Out of memory at iteration " << ii << "  size: " << cs_nnz(pX.get()) << endl;
                return(walk_rslt);
            }
            cs_spfree(x2);
            x2 = x3;
        }

        //release x and reset to x2
        pX.reset(x2);
    }

    //cs_print(pX.get(), 0);
    //cout << "pX" << endl;
    //go through result x and find the maximum value
    double max_val = 0;
    for (int p = pX.get()->p[0] ; p < pX.get()->p[1] ; p++) max_val = max(max_val, pX.get()->x[p]);

    //collect value above max*threshold

    vcout << "walking completed, nnz: " << cs_nnz(pX.get()) << endl;
    for (int p = pX.get()->p[0] ; p < pX.get()->p[1] ; p++) {
        if (pX.get()->x[p] > max_val*config.sig_threshold)
            walk_rslt.push_back(pair<int, double>(pX.get()->i[p], pX.get()->x[p]));
    }

    // sort by descend order
    sort( walk_rslt.begin(),
          walk_rslt.end(),
          [](const vertex_prob_t &v1, const vertex_prob_t &v2){return v1.second>v2.second;}
    );

    //restrict number of vertices returned
    if (config.n_max_explore_vertices>0 && walk_rslt.size()>config.n_max_explore_vertices)
        walk_rslt.resize(config.n_max_explore_vertices);

    //output sort result
    vcout << "walk result: " << vertex_id << endl;
    for (vector<vertex_prob_t>::iterator it=walk_rslt.begin(); it!=walk_rslt.end(); ++it)
        vcout << "(" << (*it).first << "," << (*it).second << " )";
    vcout << endl;

    return(walk_rslt);
}

void show_collector(const cluster_collector_t &collector){
    //show collector
    cout << "collector content (" << collector.size() << "):" << endl;
    for (auto c : collector) {
        cout << c.first << "(" << c.second.first.size() << "): (";
        for (auto d : c.second.first) cout << d << " ";
        cout << ")(";
        for (auto d : c.second.second) cout << d << " ";
        cout << ")" << endl;
    }
    cout << endl;
}

//this function saves the collector into a file. The nodes id are the internal graph ids, 0-n_vertices-1
void save_collector(const cluster_collector_t &collector, const string &fname) {
    ofstream collector_file(fname);
    if(!collector_file.is_open()){
        cout << "Unable to open file " << fname << " to write" << endl;
        return;
    }
    boost::archive::text_oarchive oarch(collector_file);
    oarch << collector;


    /*
    for (auto c : collector) {
        collector_file << c.first << "(" << c.second.first.size() << "): (";
        for (auto d : c.second.first) collector_file << d << " ";
        collector_file << ")(";
        for (auto d : c.second.second) collector_file << d << " ";
        collector_file << ")" << endl;
    }*/

    collector_file.close();
}

//this function loads the collector from a file.
void load_collector(cluster_collector_t &collector, const string &fname) {
    ifstream collector_file(fname);
    if (!collector_file.is_open()){
        cout << "Unable to open file" << fname << " to read" << endl;
        return;
    }

    boost::archive::text_iarchive iarch(collector_file);
    iarch >> collector;

    collector_file.close();
}

//this function saves the cluster result into a file. The format is a white space separated nodes ids.
//Each line is a cluster
template <typename Graph> void save_cluster_rslt(const cluster_collector_t &collector, const string &fname, const Graph& g){
    ofstream collector_file(fname);
    if(!collector_file.is_open()){
        cout << "Unable to open file " << fname << " to write" << endl;
        return;
    }

    //property type
    typename boost::property_map <Graph, boost::vertex_index2_t>::const_type vIndex2Map = boost::get(boost::vertex_index2, g);

    for (auto c : collector) {
        for (auto d : c.second.first) collector_file << boost::get(vIndex2Map, d) << " ";
        collector_file << endl;
    }
    collector_file.close();
}

//this function saves the local cluster result into a file. The format is a white space separated nodes ids.
//The file contains one line. The line is white space separated ids. The first id is the seed.
//input: clst_vector: the cluster vector after the first random walk, the first vertex is always the seed
template <typename Graph> void save_local_cluster_rslt(const vector<int> &clst_vector, const cluster_collector_t &collector, const string &fname, const Graph& g){
    ofstream collector_file(fname);
    if(!collector_file.is_open()){
        cout << "Unable to open file " << fname << " to write" << endl;
        return;
    }

    //property type
    typename boost::property_map <Graph, boost::vertex_index2_t>::const_type vIndex2Map = boost::get(boost::vertex_index2, g);

    //output seed and important vertices first
    for (auto v : clst_vector) collector_file << boost::get(vIndex2Map, v) << " ";

    for (auto c : collector) {
        if (find(c.second.first.begin(), c.second.first.end(), clst_vector[0])!=c.second.first.end()) {
            for (auto d : c.second.first)
                if (find(clst_vector.begin(), clst_vector.end(), d)==clst_vector.end())
                    collector_file << boost::get(vIndex2Map, d) << " ";
            break;
        }
    }
    collector_file << endl;
    collector_file.close();
}


//handle walk result, collect clusters from the most dominant one
void handle_walk_rslt(int vertexId, const vertex_walk_rslt_t& walk_rslt, cluster_collector_t &collector, const lrw_config& config) {
    if (walk_rslt.size()==0){
        //something wrong here
        return;
    }

    vector<int> ftrvec;
    for (vertex_walk_rslt_t::const_iterator rslt_it = walk_rslt.begin(); rslt_it!=walk_rslt.end(); ++rslt_it)
        ftrvec.push_back((*rslt_it).first);

    //get the first vertex id, which is the dominant vertices
    int clstId = walk_rslt[0].first;
    cluster_collector_t::iterator it = collector.find( clstId );

    collector[clstId].first.push_back(vertexId);
    collector[clstId].second.insert(ftrvec.begin(), ftrvec.end());

    //show collector
    //show_collector();
}

//merge duplicate clusters
void merge_clsts(cluster_collector_t &collector, const cs_sparse *adj, const lrw_config &config){
    vcout << "before cluster merge" << endl;
    if(config.verbose) show_collector(collector);

    vcout << "start merging clustering ... " << endl;
    //first do a slow version
    for (cluster_collector_t::iterator it = collector.begin(); it!=collector.end(); ++it){
        for (cluster_collector_t::iterator it2 = next(it,1); it2!=collector.end();){
            //cout << (*it).first << "--" << (*it2).first << endl;
            set<int> inter_set;

            cluster_feature_t &clst1_ftr((*it).second.second);
            cluster_feature_t &clst2_ftr((*it2).second.second);

            set_intersection(clst1_ftr.begin(), clst1_ftr.end(),
                             clst2_ftr.begin(), clst2_ftr.end(),
                             inserter(inter_set, inter_set.begin()));

            if (inter_set.size() > config.int_threshold*min(clst1_ftr.size(), clst2_ftr.size())) {
                vcout << (*it).first << "--" << (*it2).first << ":" << inter_set.size() << endl;
                //merge the second cluster to the first one
                (*it).second.first.reserve((*it).second.first.size()+(*it2).second.first.size());
                (*it).second.first.insert((*it).second.first.end(), (*it2).second.first.begin(), (*it2).second.first.end());
                (*it).second.second.insert((*it2).second.second.begin(), (*it2).second.second.end());
                //remove the second one
                it2 = collector.erase(it2);
            } else {
                ++it2;
            }
        }
    }

    if (config.neighbor_following) {
        vcout << "Start neighbor following stage" << endl;
        //get number of nodes
        int n=adj->n;

        bool allDone=false;

        //node assign history, to prevent infinite loop
        map <int, set<int> > nodeAssignHistory;

        //repeat untill all done
        while (!allDone) {
            allDone=true;

        //collect nodes cluster Id vector
        vector<int> clstVec(n);
        for (cluster_collector_t::iterator it = collector.begin(); it!=collector.end(); ++it){
            const cluster_list_t &nodeList((*it).second.first);
            int clstId=(*it).first;
            for (cluster_list_t::const_iterator it2=nodeList.begin(); it2!=nodeList.end(); ++it2) {
                clstVec[*it2]=clstId;
            }
        }

        //stores the reassign of nodes
        //tuple: <nodeId, originalClusterId, newClusterId>
        vector<tuple<int, int, int>> nodeReassignList;

        int p, *Ap, *Ai;
        Ap = adj->p; Ai=adj->i;

        for (int ii=0; ii<n; ++ii) {
            map<int, int> nbrClstCount;
            int lgstClstr, lgstClstrCount=0;
            int secondLgstClstr, secondLgstClstrCount=0;
            //find all non-zero elements of each column
            for (p = Ap[ii]; p < Ap[ii+1]; p++) {
                int nbrClstrId = clstVec[Ai[p]];
                if (nbrClstCount.find(nbrClstrId)!=nbrClstCount.end())
                    nbrClstCount[nbrClstrId]++;
                else
                    nbrClstCount[nbrClstrId]=1;

                if (nbrClstCount[nbrClstrId] > lgstClstrCount) {
                    lgstClstrCount = nbrClstCount[nbrClstrId];
                    lgstClstr = nbrClstrId;
                } else if (nbrClstCount[nbrClstrId] > secondLgstClstrCount) {
                    secondLgstClstrCount = nbrClstCount[nbrClstrId];
                    secondLgstClstr = nbrClstrId;
                }
            }
            //check if we need to change the cluster association
            if ((clstVec[ii]!=lgstClstr && lgstClstrCount>secondLgstClstrCount) ||
                    (clstVec[ii]!=lgstClstr && clstVec[ii]!=secondLgstClstr)) {
                if (nodeAssignHistory[ii].find(lgstClstr)==nodeAssignHistory[ii].end()) {
                    nodeReassignList.push_back(tuple<int, int,int>(ii, clstVec[ii], lgstClstr));
                    //nodeAssignHistory[ii].insert(clstVec[ii]);
                    nodeAssignHistory[ii].insert(lgstClstr);
                }
            }
        }

        //reassign nodes
        vcout << "reassigning nodes" << endl;
        for (auto it=nodeReassignList.begin(); it!=nodeReassignList.end(); it++) {
            allDone=false;
            vcout << "node:" << get<0>(*it) << "  " << get<1>(*it) << "==>" << get<2>(*it) << endl;
            int nodeId=get<0>(*it);
            cluster_collector_t::iterator fromIt=collector.find(get<1>(*it));
            cluster_collector_t::iterator toIt=collector.find(get<2>(*it));
            if (fromIt != collector.end() && toIt!=collector.end()) {
                cluster_list_t &fromList(fromIt->second.first);
                cluster_list_t &toList(toIt->second.first);
                fromList.erase(std::remove(fromList.begin(), fromList.end(), nodeId), fromList.end());
                toList.push_back(nodeId);
                if (fromList.size()==0) {
                    vcout << "cluster removed : " << fromIt->first << endl;
                    collector.erase(fromIt);
                }
            }
        }
        }
    }

    vcout << "merging completed" << endl;

}

//lrw algorithm
//input: graph, config
//output: clustering result (bidirectional bimap data structure)
//  id1->id1, id2->id1, id3->id1, ... id4->id4, ...
template <typename Graph, typename ClusteringResult> void lrw(const Graph& g1, ClusteringResult &clustering_result, const lrw_config& config){
    vcout << "start clustering..." << endl;

    //get adjacency matrix
    cs_sparse *adj = get_adjacency_matrix(g1);


    //keep the adj object safe
    //define a lambda destruction function
    unique_ptr<cs_sparse, decltype(csfree_fun)>pAdj(adj, csfree_fun);

    //normalise adjacency matrix by columns
    //cs_sparse is CRC format, thus all operation along column are quite efficient
    if (cs_normalise_c(adj, 1)!=1) return;

    //cs_print(adj, 0);
    //cout << "adj" << endl <<endl;

    //check performance of vector wise loop or matrix level walk
    //loop vectors one by one is far more faster than working on a large sparse matrix

    //define probe vertices
    const int n_probe_vertices = boost::num_vertices(g1);

    cluster_collector_t collector;
    //local clustering
    if (config.seed!=-1){
        lrw_config tmpConfig = config;

        //walk start from a seed
        tmpConfig.sig_threshold = 0.005;
        const vertex_walk_rslt_t& walk_rslt=explore_vertex(adj, config.seed, tmpConfig);

        //find the important attraction vertices
        vector<int> clst_vector;

        clst_vector.push_back(config.seed);

        vertex_walk_rslt_t::const_iterator first_it = walk_rslt.begin();
        for (vertex_walk_rslt_t::const_iterator rslt_it = first_it; rslt_it!=walk_rslt.end(); ++rslt_it) {
            if((*rslt_it).second>config.sig_threshold*(*first_it).second)
                if ((*rslt_it).first!=config.seed)
                    clst_vector.push_back((*rslt_it).first);
            else
                break;
        }

        //loop all vertices and walk on them
        for (vector<vertex_prob_t>::const_iterator it=walk_rslt.begin(); it!=walk_rslt.end(); ++it) {
            //vcout << "(" << (*it).first << "," << (*it).second << " )";
            const vertex_walk_rslt_t& walk_rslt2=explore_vertex(adj, (*it).first, config);
            handle_walk_rslt((*it).first, walk_rslt2, collector, config);
        }
        vcout << "Local clustering completed" << endl;

        //merge collector
        merge_clsts(collector, adj, config);

        //show local clustering result
        show_collector(collector);

        //save local clustering result
        save_local_cluster_rslt(clst_vector, collector, config.output_base+".LRW.local"+".txt", g1);
        return;
    }

    //global clustering algorithm
    boost::timer::cpu_timer t;

    //iterate all vertices
    int ii;
//#pragma omp parallel for private(ii)
    for(ii=0; ii<n_probe_vertices; ++ii) {
        const vertex_walk_rslt_t& walk_rslt=explore_vertex(adj, ii, config);
        if (walk_rslt.size()==0){
            //out of memory
            vcout << "Blah... Something wrong, might be out of memory or data incorrect!";
            //sleep(2);
            //ii--; //try again
        } else {
            handle_walk_rslt(ii, walk_rslt, collector, config);
        }
    }

    cout << "explore completed: " << t.elapsed().wall/1.0E9 << endl;

    //save collector
    //show_collector(collector);
    save_collector(collector, config.output_base+".explored.txt");

    t.start();

    merge_clsts(collector, adj, config);
    vcout << "merge completed: " << t.elapsed().wall/1.0E9 << endl;

    //output result
    show_collector(collector);

    vcout << "cluster completed!" << endl;

    //save collector result
    save_collector(collector, config.output_base+".collector.LRW.txt");

    //save cluster result
    save_cluster_rslt(collector, config.output_base+".LRW.txt", g1);
}


template <typename Graph> void merge_save_rslt(const Graph &g1, const cs_sparse *adj, cluster_collector_t &collector, const lrw_config& config) {
    boost::timer::cpu_timer t;

    merge_clsts(collector, adj, config);
    vcout << "merge completed: " << t.elapsed().wall/1.0E9 << endl;

    //output result
    //show_collector(collector);

    //save collector result
    save_collector(collector, config.output_base+".collector.LRW.txt");

    //save cluster result
    save_cluster_rslt(collector, config.output_base+".LRW.txt", g1);

    vcout << "Merge completed: " << t.elapsed().wall/1.0E9 << "seconds " << endl;
}

template <typename Graph, typename ClusteringResult> void lrw_mpi(const Graph& g1, ClusteringResult &clustering_result, const lrw_config& config){
#ifdef MPI_SUPPORT
    vcout << "start clustering mpi..." << endl;

    //MPI_Status mpi_stat;
    boost::mpi::communicator world;

    //MPI return error handling is not mandatory. If a function failed, the default error handler will
    //abort current task before function returns
    int n_tasks = world.size();
    int rank = world.rank();

    vcout << "Number of tasks: " << n_tasks << endl;


    //sprase matrix of the graph
    cs_sparse *adj=NULL;
    int n_vertices = 0;

    if (config.stage==1){
        if (n_tasks > 1) {
            cout << "This is single thread operation. No need to waste CPU resource." << endl;
            exit(1);
        }
        //do preprocessing
        //get adjancency matrix in triplet format
        adj = get_adjacency_matrix(g1);

        //save
        FILE *f = fopen((config.output_base+".s1").c_str(), "wb");
        int err = cs_save_compressed(f, adj);
        if (err!=0) cout << "save compressed matrix failed!\n" << endl;
        fclose(f);

        cs_spfree(adj);
        cout << "Stage 1 intermediate file created." << endl;
        return;
    }

    if (config.stage == 3) {
        if (n_tasks > 1) {
            cout << "This is single thread operation. No need to waste CPU resource." << endl;
            exit(1);
        }
        cluster_collector_t collector;
        load_collector(collector, config.output_base+Collector_File_Suffix);
        merge_save_rslt(g1, adj, collector, config);
        return;
    }

    if (config.stage==2 || config.stage==0) {
        if(n_tasks<2){
            cout << "Number of tasks must be greater than 2. Use mpirun or srun to launch the application." << endl;
            exit(1);
        }

        if (config.stage==2){
            //load sparse matrix from file
            FILE *f = fopen((config.output_base+".s1").c_str(), "rb");
            if (f==NULL) {
                cout << "unable to open the sparse matrix file from stage 1" << endl;
                exit(1);
            }
            adj = cs_load_compressed(f);
            fclose(f);

            if (adj==NULL) {
                cout << "unable to open the sparse matrix file from the stage 1" << endl;
                exit(1);
            }

        } else {
            adj = get_adjacency_matrix(g1);
        }
        n_vertices = adj->n;
    }

    //keep the adj object safe
    //define a lambda destruction function
    unique_ptr<cs_sparse, decltype(csfree_fun)>pAdj(adj, csfree_fun);



    //the number of vertices that slave will handle
    //int n_slave_work_load = ceil(double(n_vertices) / n_slaves);
    int n_slave_work_load = config.n_slave_work_load;

    if (rank==0){
        vcout << "master running" << endl;

        //store the last pointer
        int last;

        //send the vertices to explore to each slave task
        for(int ii=1; ii<n_tasks; ii++) {
            int first=(ii-1)*n_slave_work_load;
            last = min(ii*n_slave_work_load, n_vertices);
            cout << "slave: " << ii << "  range: " << first << "," << last << endl;
            pair<int, int> slave_work(first, last);

            world.send(ii, 0, slave_work);
        }

        //collector for explore result
        int completed_vertices = 0;
        cluster_collector_t collector;

        boost::timer::cpu_timer t;

        while (completed_vertices<n_vertices){
            cout << "probing ..." << endl;

            boost::mpi::status status = world.probe();

            cout << "ready to receive a result" << endl;

            //receiving collection
            walk_collection_t rslt_collection;
            world.recv(status.source(), status.tag(), rslt_collection);

            vcout << "completed: " << "received results: " << rslt_collection.size() << endl;

            //iterate all vertices
            for (pair<int, vertex_walk_rslt_t> &rslt : rslt_collection ) {
                handle_walk_rslt(rslt.first, rslt.second, collector, config);
            }

            completed_vertices += rslt_collection.size();

            //send next work load
            int first=last;
            last = min(first+n_slave_work_load, n_vertices);
            pair<int, int> slave_work(first, last);

            world.send(status.source(), 0, slave_work);
            cout << "verttices till " << last << " sent to slave " << status.source() << endl;
        }

        cout << "explore completed: " << t.elapsed().wall/1.0E9 << " seconds " << endl;

        //save collector
        //show_collector(collector);
        save_collector(collector, config.output_base+Collector_File_Suffix);

        if (config.stage == 0 || config.stage == 3) {
            merge_save_rslt(g1, adj, collector, config);
        }
        cout << "all done" << endl;

    } else {
        vcout << "slave " << rank << " running" << endl;
        //get adjacency matrix
        //cs_sparse *adj = get_adjacency_matrix(g1);



        //normalise adjacency matrix by columns
        //cs_sparse is CRC format, thus all operation along column are quite efficient
        if (cs_normalise_c(adj, 1)!=1) return;

        //check performance of vector wise loop or matrix level walk
        //loop vectors one by one is far more faster than working on a large sparse matrix

        pair<int, int> work_load;
        //receiving vertices allocation
        world.recv(0, 0, work_load);

        while(work_load.first<work_load.second){
            vcout << "slave: " << work_load.first << "," << work_load.second << endl;

            //walk result collections
            walk_collection_t rslts;

            for (int vertex_id=work_load.first; vertex_id<work_load.second; vertex_id++) {
                const vertex_walk_rslt_t& walk_rslt=explore_vertex(adj, vertex_id, config);
                rslts.push_back(pair<int, vertex_walk_rslt_t>(vertex_id, walk_rslt));
            }
            cout << "Slave " << rank << " is ready to send a result till " << work_load.second << endl;

            //send results to master
            world.send(0, 1, rslts);

            vcout << rank << " has sent completion" << endl;

            //receiving vertices allocation
            world.recv(0, 0, work_load);
        }

        vcout << rank << " has completed all tasks" << endl;
    }

#endif
}


#endif //GRAPH_CLUSTERING_LRW
