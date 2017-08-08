#!/bin/bash

# This script test basic lrw algorithm
../build/lrw/lrw -p 2 -n -t 0.1 -u 0.1 -v -o 'test1' ../data/tg_bowtie.gdf

# This script measure time spend on graph clustering
#time ../build/lrw/lrw ../data/com-amazon.el

