#!/bin/bash

for i in {1..50}
do
    echo "Iteration $i"
    Rscript test/test_MGM.R 
    diff ../graph1.txt ../testgraph.txt 
done
