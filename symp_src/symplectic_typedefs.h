/*
 * symplectic_typedefs.h
 *
 * Type definitions for symplectic_basis
 */

#ifndef SYMPLECTIC_TYPEDEFS_H
#define SYMPLECTIC_TYPEDEFS_H

#include "kernel_typedefs.h"

#define TRI_TO_INDEX(tet_index, tet_vertex)     (4 * (tet_index) + (tet_vertex))
#define ATLEAST_TWO(a, b, c)                    ((a) && (b)) || ((a) && (c)) || ((b) && (c))

#define COPY_PATH_ENDPOINT(new, old)    {                                                       \
                                            (new)->vertex = (old)->vertex;                      \
                                            (new)->face = (old)->face;                          \
                                            (new)->tri = (old)->tri;                            \
                                            (new)->region_index = (old)->region_index;          \
                                            (new)->region = (old)->region;                      \
                                            (new)->node = (old)->node;                          \
                                            (new)->num_adj_curves = (old)->num_adj_curves;      \
                                        }

#define COPY_PATH_NODE(new, old)        {                                                           \
                                            (new)->next = NULL;                                     \
                                            (new)->prev = NULL;                                     \
                                            (new)->next_face = (old)->next_face;                    \
                                            (new)->prev_face = (old)->prev_face;                    \
                                            (new)->inside_vertex = (old)->inside_vertex;            \
                                            (new)->cusp_region_index = (old)->cusp_region_index;    \
                                            (new)->tri = (old)->tri;                                \
                                        }

static int debug = 0;

enum pos {
    START,
    FINISH
};

static int edgesThreeToFour[4][3] = {{1, 2, 3},
                                     {0, 2, 3},
                                     {0, 1, 3},
                                     {0, 1, 2}};

/**
 * Oscillating Curves
 *
 * Each oscillating curve contributes combinatorial holonomy, we store this in
 * curve[4][4] in a similar way to the curve[4][4] attribute of a Tetrahedron.
 * An array of size num_edge_classes is attached to each Tetrahedron.
 * tet->extra[edge_class]->curve[v][f] is the intersection number of
 * the oscillating curve associated to edge_class with the face 'f' of the
 * cusp triangle at vertex 'v' of tet.
 */

struct extra {
    int                         curve[4][4];            /** oscillating curve holonomy for a cusp triangle */
};

/**
 * Queue
 */

typedef struct Queue {
    int                         front;
    int                         rear;
    int                         len;
    int                         size;
    int                         *array;
} Queue ;

/**
 * Graph
 */

typedef struct EdgeNode {
    int                         y;
    struct EdgeNode             *next;
    struct EdgeNode             *prev;
} EdgeNode;

typedef struct Graph {
    EdgeNode                    *edge_list_begin;        /** header node of doubly linked list */
    EdgeNode                    *edge_list_end;          /** tailer node ... */
    int                         *degree;                 /** degree of each vertex */
    int                         *color;                  /** color a tree bipartite */
    int                         num_vertices;            /** number of vertices in the graph */
    Boolean                     directed;                /** is the graph directed */
} Graph;

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

/**
 * Path End Points
 *
 * Path endpoints can have different states of initialisation. As a convention
 * if the pointers tri and region are NULL then endpoint is not initialised.
 * If tri is not NULL and region is NULL then the endpoint is initialised but
 * the region is not known either because it has not been choosen or we have
 * split along the curve. In either of the previous cases the tri pointer is
 * still valid. If tri is not NULL and region is not NULL then the region
 * pointer is valid and tri = region->tri.
 */

typedef struct PathEndPoint {
    FaceIndex                   face;                   /** face containg the short rectangle carrying the curve */
    VertexIndex                 vertex;                 /** vertex we dive through the manifold along */
    int                         region_index;           /** index of the region the endpoint lies in */
    int                         num_adj_curves;         /** where the curve dives into the manifold */
    struct PathNode             *node;                  /** pointer to the path node which connects to the endpoint */
    struct CuspRegion           *region;                /** pointer to the region the endpoint lies in */
    struct CuspTriangle         *tri;                   /** pointer to the cusp triangle the endpoint lies in */
} PathEndPoint;

typedef struct PathNode {
    int                         cusp_region_index;
    FaceIndex                   next_face;               /** face the path crosses to the next node */
    FaceIndex                   prev_face;               /** face the path crosses to the prev node */
    VertexIndex                 inside_vertex;           /** inside vertex of the path */
    struct CuspTriangle         *tri;                    /** cusp triangle the node lies in */
    struct PathNode             *next;                   /** next node in doubly linked list */
    struct PathNode             *prev;
} PathNode;

typedef struct CurveComponent {
    int                         edge_class[2];          /** edge classes at path end points */
    int                         cusp_index;             /** which cusp does the curve lie in */
    PathNode                    path_begin;             /** header node of doubbly linked list */
    PathNode                    path_end;               /** tailer node of ... */
    PathEndPoint                endpoints[2];           /** path end points */
    struct CurveComponent       *next;                  /** next curve component in doubly linked list */
    struct CurveComponent       *prev;                  /** prev ... */
} CurveComponent;

typedef struct OscillatingCurves {
    int                         num_curves;
    int                         *edge_class;
    CurveComponent              *curve_begin;          /** array of doubly linked lists of dual curves */
    CurveComponent              *curve_end;            /** array of doubly linkek lists of dual curves */
} OscillatingCurves;

/**
 * Cusp Triangulation
 *
 * CuspTriangle stores information about a triangle in the cusp triangulation.
 * The homology curves bound a fundamental domain, and cusp regions store the
 * information for intersection of this domain with each cusp triangle. When
 * we add oscillating curves, these regions are divided further.
*/

typedef struct CuspVertex {
    int                         edge_class;
    int                         edge_index;
    EdgeClass                   *edge;
    VertexIndex                 v1;
    VertexIndex                 v2;
} CuspVertex;

typedef struct CuspTriangle {
    Tetrahedron                 *tet;                   /** tetrahedron the triangle comes from */
    Cusp                        *cusp;                  /** cusp the triangle lies in */
    int                         tet_index;              /** tet->index */
    VertexIndex                 tet_vertex;             /** vertex the triangle comes from */
    CuspVertex                  vertices[4];            /** information about each vertex */
    struct CuspTriangle         *neighbours[4];         /** triangle neighbouring a face */
    struct CuspTriangle         *next;                  /** next cusp triangle on doubly linked list */
    struct CuspTriangle         *prev;                  /** prev cusp triangle on doubly linkled list */
} CuspTriangle;

typedef struct CuspRegion {
    CuspTriangle                *tri;                   /** cusp triangle the region lies on */
    int                         tet_index;              /** tri->tetIndex */
    VertexIndex                 tet_vertex;             /** tri->tet_vertex */
    int                         index;                  /** index of the cusp region */
    int                         curve[4][4];            /** looking at face, number of curves between the region and vertex */
    Boolean                     adj_cusp_triangle[4];   /** does the region meet this edge of the cusp triangle */
    Boolean                     dive[4][4];             /** can we dive along the face into this vertex */
    int                         num_adj_curves[4][4];   /** stores the number of curves between a region and a face */
    int                         temp_adj_curves[4][4];  /** store the adj curve until pathfinding is complete */
    struct CuspRegion           *adj_cusp_regions[4];   /** index of the adjacent regions */
    struct CuspRegion           *next;                  /** next cusp region in doubly linked list */
    struct CuspRegion           *prev;                  /** prev cusp region in doubly linked list */
} CuspRegion;

typedef struct CuspStructure {
    int                         intersect_tet_index;    /** index of the intersection triangle */
    VertexIndex                 intersect_tet_vertex;   /** vertex of the intersection triangle */
    int                         num_edge_classes;       /** number of edge classes in the cusp */
    int                         num_cusp_triangles;     /** number of cusp triangle in the cusp */
    int                         num_cusp_regions;       /** number of cusp regions in the cusp */
    Triangulation               *manifold;              /** manifold */
    Cusp                        *cusp;                  /** which manifold cusp does the struct lie in */
    Graph                       *dual_graph;            /** dual graph of the cusp region */
    CuspRegion                  **dual_graph_regions;
    CuspTriangle                cusp_triangle_begin;    /** header node of doubly linked list of cusp triangles */
    CuspTriangle                cusp_triangle_end;      /** tailer node of ... */
    CuspRegion                  *cusp_region_begin;     /** array of header nodes for cusp regions, index by cusp tri */
    CuspRegion                  *cusp_region_end;       /** array of tailer nodes for ...*/
    PathNode                    train_line_path_begin;  /** header node of doubly linked list of train line node */
    PathNode                    train_line_path_end;    /** tailer node of ... */
    PathEndPoint                *train_line_endpoint[2];/** train line endpoints for [edge_index][edge_class] */
} CuspStructure;

#endif /* SYMPLECTIC_TYPEDEFS_H */
