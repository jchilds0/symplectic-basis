/*
 * log.c
 *
 * Log the program state to a file
 */

#include "symplectic_kernel.h"

#define FILENAME    "symplectic_basis.log"

void log_gluing(Triangulation *, CuspStructure **, OscillatingCurves *);
void log_train_lines(Triangulation *, CuspStructure **, OscillatingCurves *);
void log_cusp_regions(Triangulation *, CuspStructure **, OscillatingCurves *);
void log_homology(Triangulation *, CuspStructure **, OscillatingCurves *);
void log_edge_classes(Triangulation *, CuspStructure **, OscillatingCurves *);
void log_dual_curves(Triangulation *, CuspStructure **, OscillatingCurves *);
void log_inside_edge(Triangulation *, CuspStructure **, OscillatingCurves *);
void log_graph(Triangulation *, CuspStructure **, OscillatingCurves *);
void log_endpoints(Triangulation *, CuspStructure **, OscillatingCurves *);

/*
 * Types: gluing, train_lines, cusp_regions, homology, edge_indices,
 * dual_curves, inside_edge, graph, endpoints
 */

void log_structs(Triangulation *manifold, CuspStructure **cusps, OscillatingCurves *curves, char *type) {
    if (strcmp(type, "gluing") == 0) {
        log_gluing(manifold, cusps, curves);
    } else if (strcmp(type, "train_lines") == 0) {
        log_train_lines(manifold, cusps, curves);
    } else if (strcmp(type, "cusp_regions") == 0) {
        log_cusp_regions(manifold, cusps, curves);
    } else if (strcmp(type, "homology") == 0) {
        log_homology(manifold, cusps, curves);
    } else if (strcmp(type, "edge_indices") == 0) {
        log_edge_classes(manifold, cusps, curves);
    } else if (strcmp(type, "dual_curves") == 0) {
        log_dual_curves(manifold, cusps, curves);
    } else if (strcmp(type, "inside_edge") == 0) {
        log_inside_edge(manifold, cusps, curves);
    } else if (strcmp(type, "graph") == 0) {
        log_graph(manifold, cusps, curves);
    } else if (strcmp(type, "endpoints") == 0) {
        log_endpoints(manifold, cusps, curves);
    } else {
        printf("Unknown type: %s\n", type);
    }
    printf("-------------------------------\n");
}

void log_gluing(Triangulation *manifold, CuspStructure **cusps, OscillatingCurves *curves){
    int i, j, x_vertex1, x_vertex2, y_vertex1, y_vertex2;
    CuspTriangle *tri;
    CuspStructure *cusp;
    FILE *file = fopen(FILENAME, "w");
    if (file == NULL) {
        uFatalError("write_output_to_file", "symplectic_basis");
    }

    fprintf(file, "Triangle gluing info\n");
    for (i = 0; i < manifold->num_cusps; i++) {
        fprintf(file, "Boundary %d\n", i);
        cusp = cusps[i];

        for (tri = cusp->cusp_triangle_begin.next; tri != &cusp->cusp_triangle_end; tri = tri->next) {
            for (j = 0; j < 4; j++) {
                if (j == tri->tet_vertex)
                    continue;

                x_vertex1 = (int) remaining_face[tri->tet_vertex][j];
                x_vertex2 = (int) remaining_face[j][tri->tet_vertex];
                y_vertex1 = EVALUATE(tri->tet->gluing[j], x_vertex1);
                y_vertex2 = EVALUATE(tri->tet->gluing[j], x_vertex2);

                fprintf(file, "    (Tet Index: %d, Tet Vertex: %d) Cusp Edge %d glues to "
                       "(Tet Index: %d, Tet Vertex: %d) Cusp Edge %d. (%d -> %d, %d -> %d)\n",
                       tri->tet_index,               // Tet Index
                       tri->tet_vertex,                // Tet Vertex
                       j,      // Cusp Edge
                       tri->tet->neighbor[j]->index,                              // Tet Index
                       EVALUATE(tri->tet->gluing[j], tri->tet_vertex),             // Tet Vertex
                       EVALUATE(tri->tet->gluing[j], j),   // Cusp Edge
                       x_vertex1, y_vertex1,
                       x_vertex2, y_vertex2
                );
            }
        }
    }
}

void log_train_lines(Triangulation *manifold, CuspStructure **cusps, OscillatingCurves *curves){
    int i, j, k;
    PathNode *path_node;
    CuspStructure *cusp;
    PathEndPoint *endpoint;

    printf("Train Lines\n");
    for (i = 0; i < manifold->num_cusps; i++) {
        printf("Boundary %d\n", i);

        cusp = cusps[i];
        printf("    Train Line Path: \n");

        for (path_node = cusp->train_line_path_begin.next; path_node != &cusp->train_line_path_end; path_node = path_node->next) {
            printf("        Node %d: (Tet Index %d, Tet Vertex %d) Next Face: %d, Prev Face: %d, Inside Vertex: %d\n",
                   path_node->cusp_region_index, path_node->tri->tet_index, path_node->tri->tet_vertex,
                   path_node->next_face, path_node->prev_face, path_node->inside_vertex
            );
        }

        printf("    Train Line Endpoints\n");
        for (j = 0; j < cusp->num_edge_classes; j++) {
            for (k = 0; k < 2; k++) {
                if (cusp->train_line_endpoint[k][j].tri == NULL)
                    continue;

                endpoint = &cusp->train_line_endpoint[k][j];
                printf("        Region %d (Tet Index %d, Tet Vertex %d) Face %d Vertex %d Edge Class (%d, %d)\n",
                       endpoint->region_index, endpoint->tri->tet_index,
                       endpoint->tri->tet_vertex, endpoint->face, endpoint->vertex,
                       endpoint->tri->vertices[endpoint->vertex].edge_class,
                       endpoint->tri->vertices[endpoint->vertex].edge_index);
            }
        }
    }
}

void log_cusp_regions(Triangulation *manifold, CuspStructure **cusps, OscillatingCurves *curves){
    int i, j, v1, v2, v3;
    CuspRegion *region;
    CuspStructure *cusp;

    printf("Cusp Region info\n");

    for (i = 0; i < manifold->num_cusps; i++) {
        printf("Boundary %d\n", i);

        cusp = cusps[i];
        for (j = 0; j < 4 * cusp->manifold->num_tetrahedra; j++) {
            printf("    Cusp Triangle (Tet Index %d Tet Vertex %d)\n", j / 4, j % 4);
            for (region = cusp->cusp_region_begin[j].next;
                 region != &cusp->cusp_region_end[j]; region = region->next) {
                v1 = edgesThreeToFour[region->tet_vertex][0];
                v2 = edgesThreeToFour[region->tet_vertex][1];
                v3 = edgesThreeToFour[region->tet_vertex][2];

                printf("    Region %d (Tet Index: %d, Tet Vertex: %d) (Adj Tri: %d, %d, %d) (Adj Regions: %d, %d, %d) "
                       " (Curves: [%d %d] [%d %d] [%d %d]) (Adj Curves: [%d %d] [%d %d] [%d %d]) (Dive: [%d %d] [%d %d] [%d %d])\n",
                       region->index, region->tet_index, region->tet_vertex,
                       region->adj_cusp_triangle[v1], region->adj_cusp_triangle[v2], region->adj_cusp_triangle[v3],
                       region->adj_cusp_regions[v1] == NULL ? -1 : region->adj_cusp_regions[v1]->index,
                       region->adj_cusp_regions[v2] == NULL ? -1 : region->adj_cusp_regions[v2]->index,
                       region->adj_cusp_regions[v3] == NULL ? -1 : region->adj_cusp_regions[v3]->index,
                       region->curve[v2][v1], region->curve[v3][v1],
                       region->curve[v1][v2], region->curve[v3][v2],
                       region->curve[v1][v3], region->curve[v2][v3],
                       region->num_adj_curves[v2][v1], region->num_adj_curves[v3][v1],
                       region->num_adj_curves[v1][v2], region->num_adj_curves[v3][v2],
                       region->num_adj_curves[v1][v3], region->num_adj_curves[v2][v3],
                       region->dive[v2][v1], region->dive[v3][v1],
                       region->dive[v1][v2], region->dive[v3][v2],
                       region->dive[v1][v3], region->dive[v2][v3]
                );
            }
        }
    }

}

void log_homology(Triangulation *manifold, CuspStructure **cusps, OscillatingCurves *curves){
    int i;
    CuspTriangle *tri;
    CuspStructure *cusp;

    printf("Homology info\n");
    for (i = 0; i < manifold->num_cusps; i++) {
        cusp = cusps[i];

        printf("Boundary %d\n", i);
        printf("Intersect Tet Index %d, Intersect Tet Vertex %d\n", cusp->intersect_tet_index, cusp->intersect_tet_vertex);
        printf("    Meridian\n");

        for (tri = cusp->cusp_triangle_begin.next; tri != &cusp->cusp_triangle_end; tri = tri->next) {
            printf("        (Tet Index: %d, Tet Vertex: %d) %d %d %d %d\n",
                   tri->tet_index,
                   tri->tet_vertex,
                   tri->tet->curve[M][right_handed][tri->tet_vertex][0],
                   tri->tet->curve[M][right_handed][tri->tet_vertex][1],
                   tri->tet->curve[M][right_handed][tri->tet_vertex][2],
                   tri->tet->curve[M][right_handed][tri->tet_vertex][3]
            );
        }
        printf("    Longitude\n");
        for (tri = cusp->cusp_triangle_begin.next; tri != &cusp->cusp_triangle_end; tri = tri->next) {
            printf("        (Tet Index: %d, Tet Vertex: %d) %d %d %d %d\n",
                   tri->tet_index,
                   tri->tet_vertex,
                   tri->tet->curve[L][right_handed][tri->tet_vertex][0],
                   tri->tet->curve[L][right_handed][tri->tet_vertex][1],
                   tri->tet->curve[L][right_handed][tri->tet_vertex][2],
                   tri->tet->curve[L][right_handed][tri->tet_vertex][3]
            );
        }
    }
}

void log_edge_classes(Triangulation *manifold, CuspStructure **cusps, OscillatingCurves *curves){
    int i, v1, v2, v3;
    CuspTriangle *tri;
    CuspStructure *cusp;

    printf("Edge classes\n");

    for (i = 0; i < manifold->num_cusps; i++) {
        printf("Boundary %d\n", i);

        cusp = cusps[i];
        for (tri = cusp->cusp_triangle_begin.next; tri != &cusp->cusp_triangle_end; tri = tri->next) {
            v1 = edgesThreeToFour[tri->tet_vertex][0];
            v2 = edgesThreeToFour[tri->tet_vertex][1];
            v3 = edgesThreeToFour[tri->tet_vertex][2];

            printf("    (Tet Index: %d, Tet Vertex: %d) Vertex %d: (%d %d), "
                   "Vertex %d: (%d %d), Vertex %d: (%d %d)\n",
                   tri->tet_index, tri->tet_vertex,
                   v1, tri->vertices[v1].edge_class, tri->vertices[v1].edge_index,
                   v2, tri->vertices[v2].edge_class, tri->vertices[v2].edge_index,
                   v3, tri->vertices[v3].edge_class, tri->vertices[v3].edge_index
            );
        }
    }
}

void log_dual_curves(Triangulation *manifold, CuspStructure **cusps, OscillatingCurves *curves){
    int i, j;
    PathNode *path_node;
    CurveComponent *path;

    printf("Oscillating curve paths\n");

    // which dual curve
    for (i = 0; i < curves->num_curves; i++) {
        j = 0;

        printf("Dual Curve %d\n", i);
        // which curve component
        for (path = curves->curve_begin[i].next; path != &curves->curve_end[i]; path = path->next) {
            printf("    Part %d: \n", j);

            for (path_node = path->path_begin.next;
                 path_node != &path->path_end;
                 path_node = path_node->next)
                printf("        Node %d: (Tet Index %d, Tet Vertex %d) Next Face: %d, Prev Face: %d, Inside Vertex: %d\n",
                       path_node->cusp_region_index, path_node->tri->tet_index, path_node->tri->tet_vertex,
                       path_node->next_face, path_node->prev_face, path_node->inside_vertex
                );
            j++;
        }
    }
}

void log_inside_edge(Triangulation *manifold, CuspStructure **cusps, OscillatingCurves *curves){
    int i;
    CuspTriangle *tri;
    CuspStructure *cusp;

    printf("Inside edge info\n");

    for (i = 0; i < manifold->num_cusps; i++) {
        printf("Boundary %d\n", i);

        cusp = cusps[i];
        for (tri = cusp->cusp_triangle_begin.next; tri != &cusp->cusp_triangle_end; tri = tri->next) {
            printf("    (Tet Index: %d, Tet Vertex: %d) Edge label (%d, %d, %d)\n",
                   tri->tet_index,               // Tet Index
                   tri->tet_vertex,                // Tet Vertex
                   edge3_between_faces[edgesThreeToFour[tri->tet_vertex][1]][edgesThreeToFour[tri->tet_vertex][2]],
                   edge3_between_faces[edgesThreeToFour[tri->tet_vertex][0]][edgesThreeToFour[tri->tet_vertex][2]],
                   edge3_between_faces[edgesThreeToFour[tri->tet_vertex][0]][edgesThreeToFour[tri->tet_vertex][1]]
            );
        }
    }
}

void log_graph(Triangulation *manifold, CuspStructure **cusps, OscillatingCurves *curves){
    int i, j;
    EdgeNode *edge_node;
    Graph *g;
    CuspStructure *cusp;

    printf("Graph info\n");

    for (i = 0; i < manifold->num_cusps; i++) {
        cusp = cusps[i];

        printf("Boundary %d\n", i);
        g = cusp->dual_graph;
        for (j = 0; j < g->num_vertices; j++) {
            if (cusp->dual_graph_regions[j] == NULL)
                continue;

            printf("    Vertex %d (Tet Index: %d, Tet Vertex: %d): ", j,
                   cusp->dual_graph_regions[j]->tet_index,
                   cusp->dual_graph_regions[j]->tet_vertex
            );
            for (edge_node = g->edge_list_begin[j].next;
                 edge_node != &g->edge_list_end[j];
                 edge_node = edge_node->next)
                printf("%d ", edge_node->y);

            printf("\n");
        }
    }
}

void log_endpoints(Triangulation *manifold, CuspStructure **cusps, OscillatingCurves *curves){
    int i, j, k;
    CurveComponent *path;

    printf("EndPoint Info\n");

    // which curve
    for (i = 0; i < curves->num_curves; i++) {
        printf("Dual Curve %d\n", i);

        j = 0;
        // which component
        for (path = curves->curve_begin[i].next; path != &curves->curve_end[i]; path = path->next) {
            printf("    Part %d Cusp %d\n", j, path->endpoints[0].tri->tet->cusp[path->endpoints[0].tri->tet_vertex]->index);
            for (k = 0; k < 2; k++) {
                if (k == 0)
                    printf("        Start: ");
                else
                    printf("        End:   ");

                printf("Region %d (Tet Index %d, Tet Vertex %d) Face %d Vertex %d Edge Class (%d, %d) Adj Curves %d\n",
                       path->endpoints[k].region_index, path->endpoints[k].tri->tet_index,
                       path->endpoints[k].tri->tet_vertex, path->endpoints[k].face, path->endpoints[k].vertex,
                       path->endpoints[k].tri->vertices[path->endpoints[k].vertex].edge_class,
                       path->endpoints[k].tri->vertices[path->endpoints[k].vertex].edge_index,
                       path->endpoints[k].num_adj_curves);
            }

            j++;
        }
    }
}

