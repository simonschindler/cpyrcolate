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
    int aux_0=-1, 
    int aux_1=-1
):
    """
    Runs bond percolation on a given edge list.
    Edge permutation is handled internally in Python/NumPy before passing to C.
    """
    cdef int num_edges = edges_u.shape[0]
    
    # Generate random permutation array in Python/NumPy
    cdef cnp.ndarray[cnp.int32_t, ndim=1] order = np.random.permutation(num_edges).astype(np.int32)
    
    # Prepare output arrays
    cdef cnp.ndarray[cnp.int32_t, ndim=1] out_max_size = np.zeros(num_edges + 1, dtype=np.int32)
    cdef cnp.ndarray[cnp.int32_t, ndim=1] out_span = np.zeros(num_edges + 1, dtype=np.int32)
    
    # Call the pure C function
    newman_ziff_bond(
        num_nodes_main,
        total_nodes,
        num_edges,
        &edges_u[0],
        &edges_v[0],
        &order[0],
        aux_0,
        aux_1,
        &out_max_size[0],
        &out_span[0]
    )
    
    return out_max_size, out_span
