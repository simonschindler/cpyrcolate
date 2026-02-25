#ifndef PERCOLATE_CORE_H
#define PERCOLATE_CORE_H

void newman_ziff_bond(
    int num_nodes_main,
    int total_nodes,
    int num_edges,
    const int* edges_u,
    const int* edges_v,
    const int* order,
    int aux_0,
    int aux_1,
    int* out_max_size,
    int* out_span
);

#endif
