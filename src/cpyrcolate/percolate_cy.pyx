# cython: language_level=3, boundscheck=False, wraparound=False, initializedcheck=False
import numpy as np
cimport numpy as cnp

# Declare the C function signature
cdef extern from "percolate_core.h":
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
    )


def run_percolation(
    int num_nodes_main, 
    int total_nodes, 
    cnp.ndarray[cnp.int32_t, ndim=1] edges_u, 
    cnp.ndarray[cnp.int32_t, ndim=1] edges_v, 
    cnp.ndarray[cnp.int32_t, ndim=1] order,  # <-- Added order as an argument
    int aux_0=-1, 
    int aux_1=-1
):
    cdef int num_edges = edges_u.shape[0]
    
    # Prepare output arrays
    cdef cnp.ndarray[cnp.int32_t, ndim=1] out_max_size = np.zeros(num_edges + 1, dtype=np.int32)
    cdef cnp.ndarray[cnp.int32_t, ndim=1] out_span = np.zeros(num_edges + 1, dtype=np.int32)
    
    # Call the pure C function
    newman_ziff_bond(
        num_nodes_main,
        total_nodes,
        num_edges,
        <const int*> &edges_u[0],
        <const int*> &edges_v[0],
        <const int*> &order[0],
        aux_0,
        aux_1,
        <int*> &out_max_size[0],
        <int*> &out_span[0]
    )
    
    return out_max_size, out_span
