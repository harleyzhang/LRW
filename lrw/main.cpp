#include <iostream>
#include <string>
#include <stdlib.h>
#include <boost/graph/adjacency_list.hpp>
#include <boost/filesystem.hpp>
#include <boost/algorithm/string.hpp>

#include "global.hpp"
#include "lrw.hpp"
#include "hcluster.hpp"

#include "tools.hpp"
#include "gdf.hpp"
#include "el.hpp"

#ifdef MPI_SUPPORT
#include <boost/mpi.hpp>
#endif

using namespace std;

void print_help(string basename){
  cout << "Usage: ./" << basename << " [options] [gdf_file]" << endl << endl;
  cout << "Options: " << endl;
  cout << "  -p <float>: inflation power" << endl;
  cout << "  -i <int>: max iterations" << endl;
  cout << "  -v: verbose mode" << endl;
  cout << "  -V: verbose verbose mode" << endl;
  cout << "  -m: parallel version" << endl;
  cout << "  -c <string>: local clustering method for a seed(node id)" << endl;
  cout << "  -o <directory>: output directory" << endl;
  cout << "  -t <float>: a number between 0 and 1, defines significance thresholding value" << endl;
  cout << "  -d <int>:  the inflation delay value. Inflation is done when iteration is above this value. Default 0" << endl;
  cout << "  -n: neighbor following stage. after the the cluster has been made, each nodes is assigned to its majority neighbor" << endl;
  cout << "  -u <float>: a number between 0 and 1, defines the cluster merge threshold" << endl;
  cout << "  -l <num>: number of vertices that each slave process handles(programming feature)" << endl;
  cout << "  -r <num>: number of vertices that each slave reports in the explore result" << endl;
  cout << "  -s <num>: do clusterig in stages. stage 1: preprocessing, stage 2: vertices exploration, stage 3: merging. Only stage2 can be processed in parellel.";
  cout << "  -h: print this page" << endl;
}


int main(int argc, char **argv) {

    //call MPI_INIT before handling the arguments
#ifdef MPI_SUPPORT
    boost::mpi::environment env(argc, argv);
#endif

    //parsing input arguments
    string fname("tg_bowtie_node.gdf");
    int c;
    lrw_config config;
    bool is_parallel = false;
    string seedId="";

    boost::filesystem::path exe_path(argv[0]);
    string basename = exe_path.filename().string();
    //string basename = "lrw";

    while ((c = getopt (argc, argv, "mp:i:vVl:o:r:s:c:t:u:d:n")) != -1)
        switch (c)
        {
        case 'p':
            config.inflation_power = atof(optarg);
            break;
        case 'i':
            config.max_iteration = atoi(optarg);
            break;
        case 'v':
            config.verbose = 1;
            break;
        case 'V':
            config.verbose = 100;
            break;
        case 'l':
            config.n_slave_work_load = atoi(optarg);
            break;
        case 'r':
            config.n_max_explore_vertices = atoi(optarg);
            break;
        case 'm':
            is_parallel = true;
            break;
        case 'o':
            config.output_base=optarg;
            break;
        case 's':
            config.stage = atoi(optarg);
            break;
        case 'c':
            seedId = optarg;
            break;
        case 't':
            config.sig_threshold = atof(optarg);
            break;
        case 'u':
            config.int_threshold = atof(optarg);
            break;
        case 'd':
            config.inflation_delay = atoi(optarg);
            break;
        case 'n':
            config.neighbor_following = true;
            break;
        case 'h':
        default:
            print_help(basename);
            return 0;
        }

    //non-option argument is fname
    if(optind<argc){
        fname=argv[optind];
    } else {
        print_help(basename);
        return 0;
    }

    //set output path
    /*
    if (config.output_base.length()>0){
        basename=boost::filesystem::path(fname).filename().string();
        cout << fname << " basename: " << basename << endl;
        config.output_base += '/' + basename;
    }else{
        config.output_base=fname;
    }*/
    if (config.output_base.length()==0) {
        config.output_base=fname;
    }

    cout << "Start clustering: " << fname << endl;

    config.fname = fname;

    //define a graph
    Graph_t g(0);

    //id map, convert id to vertex descriptor
    map<string, Graph_t::vertex_descriptor> vIdMap;

    //check if the file has extension of gdf
    boost::filesystem::path ipath(fname);
    string ext_str = boost::filesystem::extension(ipath);
    boost::algorithm::to_lower(ext_str);
    if(config.stage!=2){
        int err;
        if (ext_str.compare(string(".gdf"))==0) {
            err = read_gdf(fname, g, vIdMap);
        } else if (ext_str.compare(string(".el"))==0 || ext_str.compare(string(".edges"))==0) {
            err = read_el(fname, g, vIdMap);
        } else {
            cout << "input file must be in gdf or el(undirected edge list) format." << endl;
            return 1;
        }

        if (err!=0) {
            return 1;
        }
    }

    //set seed
    if (!seedId.empty()){
        //get seed vertex index from the Id
        config.seed = vIdMap[seedId];
    }


    //cout << "Print graph:" << endl;
    //print_graph(g);

    //clustering_tree_t clustering_tree;
    Clustering_Result_t clustering_rslt;

    //get adjacency matrix
    //const boost::numeric::ublas::compressed_matrix<double>& adj = get_adjacency_matrix(g);
    //cout << "adjacency matrix: " << endl << adj << endl;

    if(is_parallel) {
        if(!seedId.empty()) {
            cout << "Not yet implemented. Using non-parallel method for local clustering";
            return 1;
        }
        lrw_mpi(g, clustering_rslt, config);
    } else {
        lrw(g, clustering_rslt, config);
    }

/*
    //check execution time
    boost::progress_timer timer;
    //for(int ii=0; ii<10; ii++)
        mcl(g, clustering_rslt, config);
*/

    //cout << endl << "-------------" << endl;
    //show_clusteirng_result(clustering_rslt, g);


    /*
    hierarchy_graph<Graph_t> hgraph;
    hgraph.graph_list.push_back(g);
    hierarchical_clustering(hgraph, config);

    show_hierarchy_graph(hgraph);
    */

    cout << "Testing end!" << endl;

    return 0;
}
