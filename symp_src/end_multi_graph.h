//
// Created by joshu on 5/09/2023.
//

#ifndef END_MULTI_GRAPH_H
#define END_MULTI_GRAPH_H

#include "graph.h"

typedef struct CuspEndPoint {
    int                         cusp_index;
    int                         edge_class[2];
    struct CuspEndPoint         *next;
    struct CuspEndPoint         *prev;
} CuspEndPoint;

typedef struct EndMultiGraph {
    int                         e0;                      /** edge connecting vertices of the same color */
    int                         num_edge_classes;
    int                         num_cusps;
    int                         **edges;                 /** edge_class[u][v] is the edge class of the edge u->v */
    Boolean                     *edge_classes;           /** which edge classes are in the multigraph */
    Graph                       *multi_graph;            /** tree with extra edge of cusps */
} EndMultiGraph;

EndMultiGraph           *init_end_multi_graph(Triangulation *);
void                    free_end_multi_graph(EndMultiGraph *);
void                    find_multi_graph_path(Triangulation *, EndMultiGraph *, CuspEndPoint *, CuspEndPoint *, int);

#endif /* END_MULTI_GRAPH_H */
