## Compile

g++ main.cpp -o main -std=c++11 -O2

(on Linux: add *-lpthread* in the end)


## Run

./main cost_node_sub, cost_node_delete, cost_edge_sub, cost_edge_delete, graph_origin, graph_target

For example:

./main 2 3 4 5 data/c.gxl data/d.gxl

​    
