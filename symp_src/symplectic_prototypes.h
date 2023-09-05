//
// Created by joshu on 5/09/2023.
//

#ifndef SYMPLECTIC_PROTOTYPES_H
#define SYMPLECTIC_PROTOTYPES_H

#include "symplectic_typedefs.h"

/************************************************************************/
/*                                                                      */
/*                          cusp_regions.c                              */
/*                                                                      */
/************************************************************************/

void                    init_cusp_region(CuspStructure *);

void                    free_cusp_region(CuspStructure *);

void                    update_adj_region_data(CuspStructure *);

void                    copy_region(CuspRegion *, CuspRegion *);

/************************************************************************/
/*                                                                      */
/*                         cusp_structure.c                             */
/*                                                                      */
/************************************************************************/

CuspStructure           *init_cusp_structure(Triangulation *, Cusp *, EndMultiGraph *);

void                    free_cusp_structure(CuspStructure *);

void                    construct_cusp_region_dual_graph(CuspStructure *);

/************************************************************************/
/*                                                                      */
/*                        cusp_train_line.c                             */
/*                                                                      */
/************************************************************************/

void                    init_train_line(CuspStructure *);

void                    free_train_line(CuspStructure *);

void                    do_manifold_train_lines(Triangulation *, CuspStructure **, EndMultiGraph *);

CurveComponent          *setup_train_line_component(CuspStructure *, EndMultiGraph *, CurveComponent *, CurveComponent *, CuspEndPoint *, int);

void                    do_curve_component_on_train_line(CuspStructure *, CurveComponent *);

/************************************************************************/
/*                                                                      */
/*                         end_multi_graph.c                            */
/*                                                                      */
/************************************************************************/

EndMultiGraph           *init_end_multi_graph(Triangulation *);

void                    free_end_multi_graph(EndMultiGraph *);

void                    find_multi_graph_path(Triangulation *, EndMultiGraph *, CuspEndPoint *, CuspEndPoint *, int);

/************************************************************************/
/*                                                                      */
/*                              graph.c                                 */
/*                                                                      */
/************************************************************************/

/* Generic Graph */

Graph                   *init_graph(int, Boolean);

void                    free_graph(Graph *);

int                     insert_edge(Graph *, int, int, Boolean);

void                    delete_edge(Graph *, int, int, Boolean);

Boolean                 edge_exists(Graph *, int, int);

Graph                   *spanning_tree(Graph *, int, int *);

/* Breadth First Search */

void                    init_search(Graph *, Boolean *, Boolean *, int *);

void                    bfs(Graph *, int, Boolean *, Boolean *, int *);

void                    find_path(int, int, int *, EdgeNode *, EdgeNode *);

Boolean                 cycle_exists(Graph *, int, Boolean *, Boolean *, int *, int *, int *);

int                     **ford_fulkerson(Graph *, int, int);

int                     bfs_target_list(Graph *, int, int *, int, Boolean *, Boolean *, int *);

void                    free_edge_node(EdgeNode *, EdgeNode *);

/************************************************************************/
/*                                                                      */
/*                              log.c                                   */
/*                                                                      */
/************************************************************************/

void                    log_structs(Triangulation *, CuspStructure **, OscillatingCurves *, char *);

/************************************************************************/
/*                                                                      */
/*                       oscillating_curves.c                           */
/*                                                                      */
/************************************************************************/

void                    do_oscillating_curves(Triangulation *manifold, Boolean *edge_classes);

CurveComponent          *init_curve_component(int, int, int);

void                    update_adj_curve_on_cusp(CuspStructure *);

void                    graph_path_to_dual_curve(CuspStructure *, EdgeNode *, EdgeNode *, PathNode *, PathNode *, PathEndPoint *, PathEndPoint *);

void                    split_cusp_regions_along_path(CuspStructure *, PathNode *, PathNode *, PathEndPoint *, PathEndPoint *);

void                    endpoint_edge_node_to_path_node(CuspRegion *, PathNode *, EdgeNode *, PathEndPoint *, int);

void                    interior_edge_node_to_path_node(CuspRegion *, PathNode *, EdgeNode *);

void                    split_cusp_region_path_interior(CuspRegion *, CuspRegion *, PathNode *, int);

void                    split_cusp_region_path_endpoint(CuspRegion *, CuspRegion *, PathNode *, PathEndPoint *, int, int);

void                    update_cusp_triangle_path_interior(CuspRegion *, CuspRegion *, CuspRegion *, PathNode *);

void                    update_cusp_triangle_endpoints(CuspRegion *, CuspRegion *, CuspRegion *, PathEndPoint *, PathNode *, int);

/************************************************************************/
/*                                                                      */
/*                             queue.c                                  */
/*                                                                      */
/************************************************************************/

Queue                   *init_queue(int);

Queue                   *enqueue(Queue *, int);

int                     dequeue(Queue *);

Boolean                 empty_queue(Queue *);

void                    free_queue(Queue *);

#endif /* SYMPLECTIC_PROTOTYPES_H */
