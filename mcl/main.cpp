#include <iostream>
#include <string>
#include <stdlib.h>
#include <boost/graph/adjacency_list.hpp>
#include <boost/filesystem.hpp>

#include "global.hpp"
#include "mcl.hpp"
#include "hcluster.hpp"

#include "tools.hpp"
#include "gdf.hpp"

using namespace std;

void print_help(string basename){
  cout << "Usage: ./" << basename << " [options] [gdf_file]" << endl << endl;
  cout << "Options: " << endl;
  cout << "  -p <num>: inflation power" << endl;
  cout << "  -i <num>: max iterations" << endl;
  cout << "  -c: use minimum memory consumption algorithm, does not support parallel" << endl;
  cout << "  -v: verbose mode" << endl;
  cout << "  -V: verbose verbose mode" << endl;
  cout << "  -h: print this page" << endl;
}

int main(int argc, char **argv) {
    cout << "starting" << endl;
    //parsing input arguments
    string fname("tg_bowtie_node.gdf");
    int c;
    mcl_config config;
    //string exe_name = argv[0];
    boost::filesystem::path exe_path(argv[0]);

    //string basename = exe_path.filename().string();
    string basename="mcl";

    while ((c = getopt (argc, argv, "cp:i:vV")) != -1)
        switch (c)
        {
        case 'p':
            config.inflation_power = atof(optarg);
            break;
        case 'i':
            config.max_iteration = atoi(optarg);
            break;
        case 'c':
            config.min_mem = true;
            break;
        case 'v':
            config.verbose = 1;
            break;
        case 'V':
            config.verbose = 100;
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

    cout << "Start testing..." << endl;


    //define a graph
    Graph_t g(0);

    //id map, convert id to vertex descriptor
    map<string, Graph_t::vertex_descriptor> vIdMap;

    int err = read_gdf(fname, g, vIdMap);
    if (err!=0) {
        return 1;
    }


    //cout << "Print graph:" << endl;
    //print_graph(g);

/*
    //clustering_tree_t clustering_tree;
    Clustering_Result_t clustering_rslt;
 
    //check execution time
    boost::progress_timer timer;
    //for(int ii=0; ii<10; ii++)
        mcl(g, clustering_rslt, config);

    show_clusteirng_result(clustering_rslt, g);
*/


    hierarchy_graph<Graph_t> hgraph;
    hgraph.graph_list.push_back(g);
    hierarchical_clustering(hgraph, config);

    show_hierarchy_graph(hgraph);

    cout << "Testing end!" << endl;
    return 0;
}
