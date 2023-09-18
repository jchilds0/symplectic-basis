/*
 * log.c
 *
 * Log the program state to a file
 */

#include "symplectic_kernel.h"
#include <time.h>

#define MAX_BUF_SIZE    100

void time_string(char *);

void log_gluing(Triangulation *, CuspStructure **, OscillatingCurves *);
void log_train_lines(Triangulation *, CuspStructure **, OscillatingCurves *);
void log_cusp_regions(Triangulation *, CuspStructure **, OscillatingCurves *);
void log_homology(Triangulation *, CuspStructure **, OscillatingCurves *);
void log_edge_classes(Triangulation *, CuspStructure **, OscillatingCurves *);
void log_dual_curves(Triangulation *, CuspStructure **, OscillatingCurves *);
void log_inside_edge(Triangulation *, CuspStructure **, OscillatingCurves *);
void log_graph(Triangulation *, CuspStructure **, OscillatingCurves *);
void log_endpoints(Triangulation *, CuspStructure **, OscillatingCurves *);

static FILE *file = NULL;
static int console_log = 0;

void time_string(char *buf) {
    memset(buf, '\0', MAX_BUF_SIZE);

    time_t t = time(NULL);
    struct tm tm = *localtime(&t);
    sprintf(buf, "[%d-%d-%d, %02d:%02d:%02d]", tm.tm_mday, tm.tm_mon + 1, tm.tm_year + 1900, tm.tm_hour, tm.tm_min, tm.tm_sec);
}

void start_logging(Triangulation *manifold) {
    char buf[MAX_BUF_SIZE];

    memset(buf, '\0', sizeof buf);
    sprintf(buf, "manifold-%s.log", manifold->name);

    if (debug) {
        file = fopen(buf, "w");

        if (file == NULL) {
            uFatalError("start_logging", "symplectic_basis");
        }
    }

    time_string(buf);

    if (console_log) {
        fprintf(stdout, "%s\tManifold: %s, Num. of Tetrahedra: %d, Num. of Cusps: %d\n",
                buf, manifold->name, manifold->num_tetrahedra, manifold->num_cusps);
    }
}

void finish_logging() {
    fclose(file);
}

/*
 * Types: gluing, train_lines, cusp_regions, homology, edge_indices,
 * dual_curves, inside_edge, graph, endpoints
 */

void log_structs(Triangulation *manifold, CuspStructure **cusps, OscillatingCurves *curves, char *type) {
    char buf[MAX_BUF_SIZE];
    static int i = 0;

    if (strcmp(type, "oscillating curve") == 0) {
        if (!debug)
            return;

        log_gluing(manifold, cusps, curves);
        fprintf(file,"\n");
        fprintf(file,"Oscillating Curve %d\n", i);
        fprintf(file,"\n");
        fprintf(file,"-------------------------------\n");

        log_dual_curves(cusps[0]->manifold, cusps, curves);
        log_endpoints(cusps[0]->manifold, cusps, curves);
        log_cusp_regions(cusps[0]->manifold, cusps, curves);
        log_graph(cusps[0]->manifold, cusps, curves);
        i++;
    } else if (strcmp(type, "train line") == 0) {
        if (!debug)
            return;

        fprintf(file, "\n");
        fprintf(file, "Manifold Train Lines\n");
        fprintf(file, "\n");
        fprintf(file, "-------------------------------\n");

        log_train_lines(manifold, cusps, curves);
        log_cusp_regions(manifold, cusps, curves);
        log_graph(manifold, cusps, curves);

    } else if (strcmp(type, "cusp structure") == 0) {
        if (!debug)
            return;

        fprintf(file,"\n");
        fprintf(file,"Struct Initialisation\n");
        fprintf(file,"\n");

        log_gluing(manifold, cusps, curves);
        log_homology(manifold, cusps, curves);
        log_edge_classes(manifold, cusps, curves);
        log_inside_edge(manifold, cusps, curves);
        log_cusp_regions(manifold, cusps, curves);
    } else if (console_log) {
        time_string(buf);
        fprintf(stdout, "%s\t%s\n", buf, type);
        fflush(stdout);
        return;
    }

    if (debug) {
        fprintf(file, "-------------------------------\n");
    }
}

void log_gluing(Triangulation *manifold, CuspStructure **cusps, OscillatingCurves *curves){
    int i, j, x_vertex1, x_vertex2, y_vertex1, y_vertex2;
    CuspTriangle *tri;
    CuspStructure *cusp;

    fprintf(file,"Triangle gluing info\n");
    for (i = 0; i < manifold->num_cusps; i++) {
        fprintf(file,"\tBoundary %d\n", i);
        cusp = cusps[i];

        for (tri = cusp->cusp_triangle_begin.next; tri != &cusp->cusp_triangle_end; tri = tri->next) {
            for (j = 0; j < 4; j++) {
                if (j == tri->tet_vertex)
                    continue;

                x_vertex1 = (int) remaining_face[tri->tet_vertex][j];
                x_vertex2 = (int) remaining_face[j][tri->tet_vertex];
                y_vertex1 = EVALUATE(tri->tet->gluing[j], x_vertex1);
                y_vertex2 = EVALUATE(tri->tet->gluing[j], x_vertex2);

                fprintf(file,"\t\t(Tet Index: %d, Tet Vertex: %d) Cusp Edge %d glues to "
                       "(Tet Index: %d, Tet Vertex: %d) Cusp Edge %d. (%d -> %d, %d -> %d)\n",
                       tri->tet_index, tri->tet_vertex, j, tri->tet->neighbor[j]->index,
                       EVALUATE(tri->tet->gluing[j], tri->tet_vertex), EVALUATE(tri->tet->gluing[j], j),
                       x_vertex1, y_vertex1, x_vertex2, y_vertex2
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

    fprintf(file,"Train Lines\n");
    for (i = 0; i < manifold->num_cusps; i++) {
        fprintf(file,"Boundary %d\n", i);

        cusp = cusps[i];
        fprintf(file,"    Train Line Path: \n");

        for (path_node = cusp->train_line_path_begin.next; path_node != &cusp->train_line_path_end; path_node = path_node->next) {
            fprintf(file,"        Node %d: (Tet Index %d, Tet Vertex %d) Next Face: %d, Prev Face: %d, Inside Vertex: %d\n",
                   path_node->cusp_region_index, path_node->tri->tet_index, path_node->tri->tet_vertex,
                   path_node->next_face, path_node->prev_face, path_node->inside_vertex
            );
        }

        fprintf(file,"    Train Line Endpoints\n");
        for (j = 0; j < cusp->num_edge_classes; j++) {
            for (k = 0; k < 2; k++) {
                if (cusp->train_line_endpoint[k][j].tri == NULL)
                    continue;

                endpoint = &cusp->train_line_endpoint[k][j];
                fprintf(file,"        Region %d (Tet Index %d, Tet Vertex %d) Face %d Vertex %d Edge Class (%d, %d)\n",
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

    fprintf(file,"Cusp Region info\n");

    for (i = 0; i < manifold->num_cusps; i++) {
        fprintf(file,"Boundary %d\n", i);

        cusp = cusps[i];
        for (j = 0; j < 4 * cusp->manifold->num_tetrahedra; j++) {
            fprintf(file,"    Cusp Triangle (Tet Index %d Tet Vertex %d)\n", j / 4, j % 4);
            for (region = cusp->cusp_region_begin[j].next;
                 region != &cusp->cusp_region_end[j]; region = region->next) {
                v1 = edgesThreeToFour[region->tet_vertex][0];
                v2 = edgesThreeToFour[region->tet_vertex][1];
                v3 = edgesThreeToFour[region->tet_vertex][2];

                fprintf(file,"    Region %d (Tet Index: %d, Tet Vertex: %d) (Adj Tri: %d, %d, %d) (Adj Regions: %d, %d, %d) "
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

    fprintf(file,"Homology info\n");
    for (i = 0; i < manifold->num_cusps; i++) {
        cusp = cusps[i];

        fprintf(file,"Boundary %d\n", i);
        fprintf(file,"Intersect Tet Index %d, Intersect Tet Vertex %d\n", cusp->intersect_tet_index, cusp->intersect_tet_vertex);
        fprintf(file,"    Meridian\n");

        for (tri = cusp->cusp_triangle_begin.next; tri != &cusp->cusp_triangle_end; tri = tri->next) {
            fprintf(file,"        (Tet Index: %d, Tet Vertex: %d) %d %d %d %d\n",
                   tri->tet_index,
                   tri->tet_vertex,
                   tri->tet->curve[M][right_handed][tri->tet_vertex][0],
                   tri->tet->curve[M][right_handed][tri->tet_vertex][1],
                   tri->tet->curve[M][right_handed][tri->tet_vertex][2],
                   tri->tet->curve[M][right_handed][tri->tet_vertex][3]
            );
        }
        fprintf(file,"    Longitude\n");
        for (tri = cusp->cusp_triangle_begin.next; tri != &cusp->cusp_triangle_end; tri = tri->next) {
            fprintf(file,"        (Tet Index: %d, Tet Vertex: %d) %d %d %d %d\n",
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

    fprintf(file,"Edge classes\n");

    for (i = 0; i < manifold->num_cusps; i++) {
        fprintf(file,"Boundary %d\n", i);

        cusp = cusps[i];
        for (tri = cusp->cusp_triangle_begin.next; tri != &cusp->cusp_triangle_end; tri = tri->next) {
            v1 = edgesThreeToFour[tri->tet_vertex][0];
            v2 = edgesThreeToFour[tri->tet_vertex][1];
            v3 = edgesThreeToFour[tri->tet_vertex][2];

            fprintf(file,"    (Tet Index: %d, Tet Vertex: %d) Vertex %d: (%d %d), "
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

    fprintf(file,"Oscillating curve paths\n");

    // which dual curve
    for (i = 0; i < curves->num_curves; i++) {
        j = 0;

        fprintf(file,"Dual Curve %d\n", i);
        // which curve component
        for (path = curves->curve_begin[i].next; path != &curves->curve_end[i]; path = path->next) {
            fprintf(file,"    Part %d: \n", j);

            for (path_node = path->path_begin.next;
                 path_node != &path->path_end;
                 path_node = path_node->next)
                fprintf(file,"        Node %d: (Tet Index %d, Tet Vertex %d) Next Face: %d, Prev Face: %d, Inside Vertex: %d\n",
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

    fprintf(file,"Inside edge info\n");

    for (i = 0; i < manifold->num_cusps; i++) {
        fprintf(file,"Boundary %d\n", i);

        cusp = cusps[i];
        for (tri = cusp->cusp_triangle_begin.next; tri != &cusp->cusp_triangle_end; tri = tri->next) {
            fprintf(file,"    (Tet Index: %d, Tet Vertex: %d) Edge label (%d, %d, %d)\n",
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

    fprintf(file,"Graph info\n");

    for (i = 0; i < manifold->num_cusps; i++) {
        cusp = cusps[i];

        fprintf(file,"Boundary %d\n", i);
        g = cusp->dual_graph;
        for (j = 0; j < g->num_vertices; j++) {
            if (cusp->dual_graph_regions[j] == NULL)
                continue;

            fprintf(file,"    Vertex %d (Tet Index: %d, Tet Vertex: %d): ", j,
                   cusp->dual_graph_regions[j]->tet_index,
                   cusp->dual_graph_regions[j]->tet_vertex
            );
            for (edge_node = g->edge_list_begin[j].next;
                 edge_node != &g->edge_list_end[j];
                 edge_node = edge_node->next)
                fprintf(file,"%d ", edge_node->y);

            fprintf(file,"\n");
        }
    }
}

void log_endpoints(Triangulation *manifold, CuspStructure **cusps, OscillatingCurves *curves){
    int i, j, k;
    CurveComponent *path;

    fprintf(file,"EndPoint Info\n");

    // which curve
    for (i = 0; i < curves->num_curves; i++) {
        fprintf(file,"Dual Curve %d\n", i);

        j = 0;
        // which component
        for (path = curves->curve_begin[i].next; path != &curves->curve_end[i]; path = path->next) {
            fprintf(file,"    Part %d Cusp %d\n", j, path->endpoints[0].tri->tet->cusp[path->endpoints[0].tri->tet_vertex]->index);
            for (k = 0; k < 2; k++) {
                if (k == 0)
                    fprintf(file,"        Start: ");
                else
                    fprintf(file,"        End:   ");

                fprintf(file,"Region %d (Tet Index %d, Tet Vertex %d) Face %d Vertex %d Edge Class (%d, %d) Adj Curves %d\n",
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

