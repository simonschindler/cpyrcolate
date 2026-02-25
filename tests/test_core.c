#include <stdio.h>
#include <stdlib.h>
#include "percolate_core.h"

int main() {
    // 1. Define a simple graph: 4 main nodes (0, 1, 2, 3) 
    // connected in a line: 0-1, 1-2, 2-3
    int num_nodes_main = 4;
    
    // We will add 2 auxiliary nodes (4 and 5) attached to the ends.
    // So total nodes = 6
    int total_nodes = 6;
    int aux_0 = 4;
    int aux_1 = 5;

    // 2. Define the edges (including connections to auxiliary nodes)
    // Edge 0: 0-1
    // Edge 1: 1-2
    // Edge 2: 2-3
    // Edge 3: aux_0 - 0
    // Edge 4: aux_1 - 3
    int num_edges = 5;
    int edges_u[] = {0, 1, 2, aux_0, aux_1};
    int edges_v[] = {1, 2, 3, 0, 3};

    // 3. Define the order in which edges are added (the permutation)
    // We'll add them in a specific order to watch the clusters grow
    int order[] = {0, 2, 3, 4, 1}; 

    // 4. Allocate output arrays (size = num_edges + 1)
    int *out_max_size = (int *)malloc((num_edges + 1) * sizeof(int));
    int *out_span = (int *)malloc((num_edges + 1) * sizeof(int));

    // 5. Run the algorithm
    printf("Running Newman-Ziff Bond Percolation...\n");
    newman_ziff_bond(
        num_nodes_main, 
        total_nodes, 
        num_edges, 
        edges_u, 
        edges_v, 
        order, 
        aux_0, 
        aux_1, 
        out_max_size, 
        out_span
    );

    // 6. Print the results
    printf("Step | Max Cluster Size | Spanning?\n");
    printf("-----------------------------------\n");
    for (int i = 0; i <= num_edges; i++) {
        printf("%4d | %16d | %9d\n", i, out_max_size[i], out_span[i]);
    }

    // 7. Clean up
    free(out_max_size);
    free(out_span);

    return 0;
}
