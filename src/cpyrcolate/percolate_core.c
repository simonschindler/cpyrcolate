#include <stdlib.h>
#include "percolate_core.h"

// Path compression findroot
static inline int findroot(int *ptr, int i) {
    if (ptr[i] < 0) return i;
    ptr[i] = findroot(ptr, ptr[i]);
    return ptr[i];
}

// Weighted union. Returns the size of the newly merged component.
static inline int merge(int *ptr, int r1, int r2) {
    if (ptr[r1] > ptr[r2]) { // r2 is larger (more negative)
        ptr[r2] += ptr[r1];
        ptr[r1] = r2;
        return -ptr[r2];
    } else {
        ptr[r1] += ptr[r2];
        ptr[r2] = r1;
        return -ptr[r1];
    }
}

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
) {
    // Allocate two Union-Find forests
    int* ptr_main = (int*)malloc(num_nodes_main * sizeof(int));
    int* ptr_span = (int*)malloc(total_nodes * sizeof(int));
    
    // Initialize roots to -1 (size 1)
    for(int i = 0; i < num_nodes_main; i++) ptr_main[i] = -1;
    for(int i = 0; i < total_nodes; i++) ptr_span[i] = -1;
    
    int max_size = (num_nodes_main > 0) ? 1 : 0;
    int has_span = 0;
    
    // Base case: 0 edges added
    out_max_size[0] = max_size;
    out_span[0] = 0;
    
    for (int i = 0; i < num_edges; i++) {
        int e = order[i];
        int u = edges_u[e];
        int v = edges_v[e];
        
        // 1. Spanning Forest Update (includes auxiliary nodes)
        if (!has_span && aux_0 >= 0 && aux_1 >= 0) {
            int r1_s = findroot(ptr_span, u);
            int r2_s = findroot(ptr_span, v);
            if (r1_s != r2_s) {
                merge(ptr_span, r1_s, r2_s);
                // Check if the two boundaries are now connected
                if (findroot(ptr_span, aux_0) == findroot(ptr_span, aux_1)) {
                    has_span = 1;
                }
            }
        }
        
        // 2. Main Forest Update (ignores auxiliary nodes)
        if (u < num_nodes_main && v < num_nodes_main) {
            int r1_m = findroot(ptr_main, u);
            int r2_m = findroot(ptr_main, v);
            if (r1_m != r2_m) {
                int new_size = merge(ptr_main, r1_m, r2_m);
                if (new_size > max_size) {
                    max_size = new_size;
                }
            }
        }
        
        // Record stats after adding edge 'i'
        out_max_size[i + 1] = max_size;
        out_span[i + 1] = has_span;
    }
    
    free(ptr_main);
    free(ptr_span);
}
