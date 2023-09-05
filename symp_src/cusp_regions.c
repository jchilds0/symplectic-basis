/*
 * cusp_regions.c
 *
 * Provides the functions
 *
 *      void init_cusp_region(CuspStructure *);
 *
 *      void free_cusp_region(CuspStructure *);
 *
 * which are used by cusp_structure.c to initialise and free the cusp regions.
 * It also provides the functions
 *
 *      void update_adj_region_data(CuspStructure *);
 *
 *      void copy_region(CuspRegion *, CuspRegion *);
 *
 * which are used in various places for cusp region operations.
 */


#include "symplectic_kernel.h"

int                     init_intersect_cusp_region(CuspStructure *, CuspTriangle *, int);
int                     init_intersect_vertex_two_zero_flows(CuspStructure *, CuspTriangle *, int);
int                     init_normal_cusp_region(CuspStructure *, CuspTriangle *, int);
void                    set_cusp_region_data(CuspStructure *, CuspTriangle *, const int [4], const Boolean [4], int);
CuspRegion              *find_adj_region(CuspRegion *, CuspRegion *, CuspRegion *, int);
int                     net_flow_around_vertex(CuspTriangle *, int);

/*
 * Initialise the cusp region doubly linked list to cotain the regions bounded
 * by the meridian and longitude curves.
 */

void init_cusp_region(CuspStructure *cusp) {
    int index;
    CuspTriangle *tri;

    // Header and tailer nodes.
    cusp->cusp_region_begin = NEW_ARRAY(4 * cusp->manifold->num_tetrahedra, CuspRegion);
    cusp->cusp_region_end   = NEW_ARRAY(4 * cusp->manifold->num_tetrahedra, CuspRegion);

    for (index = 0; index < 4 * cusp->manifold->num_tetrahedra; index++) {
        cusp->cusp_region_begin[index].next    = &cusp->cusp_region_end[index];
        cusp->cusp_region_begin[index].prev    = NULL;
        cusp->cusp_region_end[index].next      = NULL;
        cusp->cusp_region_end[index].prev      = &cusp->cusp_region_begin[index];
    }

    index = 0;
    for (tri = cusp->cusp_triangle_begin.next; tri != &cusp->cusp_triangle_end; tri = tri->next) {
        // Intersection vertex doesn't have a center
        if (tri->tet_index == cusp->intersect_tet_index && tri->tet_vertex == cusp->intersect_tet_vertex) {
            index = init_intersect_cusp_region(cusp, tri, index);
            continue;
        }

        index = init_normal_cusp_region(cusp, tri, index);
    }

    update_adj_region_data(cusp);
    cusp->num_cusp_regions = index;
}

void free_cusp_region(CuspStructure *cusp) {
    CuspRegion *region;

    for (int i = 0; i < 4 * cusp->manifold->num_tetrahedra; i++) {
        while (cusp->cusp_region_begin[i].next != &cusp->cusp_region_end[i]) {
            region = cusp->cusp_region_begin[i].next;
            REMOVE_NODE(region)
            my_free(region);
        }
    }

    my_free(cusp->cusp_region_begin);
    my_free(cusp->cusp_region_end);
}

/*
 * Assume peripheral_curves() has been called, and as a result the only curves
 * on the intersection triangle are those which intersect, and they give a
 * valid intersection.
 */

int init_intersect_cusp_region(CuspStructure *cusp, CuspTriangle *tri, int index) {
    int i, curve_index, vertex, v1, v2, v3;
    int distance[4];
    Boolean adj_triangle[4];

    // which vertex are we inside the flow of
    for (vertex = 0; vertex < 4; vertex++) {
        if (vertex == tri->tet_vertex) {
            continue;
        }

        v1 = (int) remaining_face[tri->tet_vertex][vertex];
        v2 = (int) remaining_face[vertex][tri->tet_vertex];

        for (i = 1; i < net_flow_around_vertex(tri, vertex); i++) {
            for (curve_index = 0; curve_index < 2; curve_index++) {
                distance[v1]                    = i;
                distance[v2]                    = MIN(distance[v1], 2 * net_flow_around_vertex(tri, vertex) - distance[v1])
                                                  + net_flow_around_vertex(tri, v2) + net_flow_around_vertex(tri, v1);
                distance[vertex]                = net_flow_around_vertex(tri, vertex)
                                                  - distance[v1] + net_flow_around_vertex(tri, v1);
                distance[tri->tet_vertex]       = -1;

                adj_triangle[v1]                = 1;
                adj_triangle[v2]                = 0;
                adj_triangle[vertex]            = 0;
                adj_triangle[tri->tet_vertex]   = -1;

                set_cusp_region_data(cusp, tri, distance, adj_triangle, index);
                index++;

                // Swap vertices
                v1 = (int) remaining_face[vertex][tri->tet_vertex];
                v2 = (int) remaining_face[tri->tet_vertex][vertex];
            }
        }

        // Region in the middle of face vertex
        if (net_flow_around_vertex(tri, v1) && net_flow_around_vertex(tri, v2)) {
            distance[v1]                    = net_flow_around_vertex(tri, v2);
            distance[v2]                    = net_flow_around_vertex(tri, v1);
            distance[vertex]                = MIN(net_flow_around_vertex(tri, v1) + distance[v1],
                                                  net_flow_around_vertex(tri, v2) + distance[v2])
                                              + net_flow_around_vertex(tri, vertex);
            distance[tri->tet_vertex]       = -1;

            adj_triangle[v1]                = 0;
            adj_triangle[v2]                = 0;
            adj_triangle[vertex]            = 1;
            adj_triangle[tri->tet_vertex]   = -1;

            set_cusp_region_data(cusp, tri, distance, adj_triangle, index);
            index++;
        }
    }

    // Region of distance 0 to vertex
    v1 = edgesThreeToFour[tri->tet_vertex][0];
    v2 = edgesThreeToFour[tri->tet_vertex][1];
    v3 = edgesThreeToFour[tri->tet_vertex][2];

    // Edge Case: Two vertices with 0 flow
    if (ATLEAST_TWO(!net_flow_around_vertex(tri, v1),
                    !net_flow_around_vertex(tri, v2),
                    !net_flow_around_vertex(tri, v3)))
        return init_intersect_vertex_two_zero_flows(cusp, tri, index);

    for (vertex = 0; vertex < 4; vertex++) {
        if (vertex == tri->tet_vertex)
            continue;

        v1 = (int) remaining_face[tri->tet_vertex][vertex];
        v2 = (int) remaining_face[vertex][tri->tet_vertex];

        distance[vertex]               = 0;
        distance[v1]                   = net_flow_around_vertex(tri, vertex) + net_flow_around_vertex(tri, v1);
        distance[v2]                   = net_flow_around_vertex(tri, vertex) + net_flow_around_vertex(tri, v2);
        distance[tri->tet_vertex]      = -1;

        adj_triangle[vertex]           = 0;
        adj_triangle[v1]               = 1;
        adj_triangle[v2]               = 1;
        adj_triangle[tri->tet_vertex]  = 0;

        set_cusp_region_data(cusp, tri, distance, adj_triangle, index);
        index++;
    }

    return index;
}

int init_intersect_vertex_two_zero_flows(CuspStructure *cusp, CuspTriangle *tri, int index) {
    int vertex, v1, v2, v3, distance[4];
    Boolean adj_triangle[4];

    v1 = (int) edgesThreeToFour[tri->tet_vertex][0];
    v2 = (int) edgesThreeToFour[tri->tet_vertex][1];
    v3 = (int) edgesThreeToFour[tri->tet_vertex][2];

    distance[v1]                   = net_flow_around_vertex(tri, v1);
    distance[v2]                   = net_flow_around_vertex(tri, v2);
    distance[v3]                   = net_flow_around_vertex(tri, v3);
    distance[tri->tet_vertex]      = -1;

    adj_triangle[v1]               = 1;
    adj_triangle[v2]               = 1;
    adj_triangle[v3]               = 1;
    adj_triangle[tri->tet_vertex]  = -1;

    set_cusp_region_data(cusp, tri, distance, adj_triangle, index);
    index++;

    // Find vertex with non-zero flow
    for (vertex = 0; vertex < 4; vertex++) {
        if (vertex == tri->tet_vertex)
            continue;

        if (net_flow_around_vertex(tri, vertex)) {
            v1 = vertex;
            v2 = (int) remaining_face[tri->tet_vertex][v1];
            v3 = (int) remaining_face[v1][tri->tet_vertex];
            break;
        }
    }
    distance[v1]                    = 0;
    distance[v2]                    = net_flow_around_vertex(tri, v1);
    distance[v3]                    = net_flow_around_vertex(tri, v1);
    distance[tri->tet_vertex]       = -1;

    adj_triangle[v1]                = 0;
    adj_triangle[v2]                = 1;
    adj_triangle[v3]                = 1;
    adj_triangle[tri->tet_vertex]   = 0;

    set_cusp_region_data(cusp, tri, distance, adj_triangle, index);

    return index + 1;
}

int init_normal_cusp_region(CuspStructure *cusp, CuspTriangle *tri, int index) {
    int i, vertex, v1, v2;
    int distance[4];
    Boolean adj_triangle[4];

    // which vertex are we inside the flow of
    for (vertex = 0; vertex < 4; vertex++) {
        if (vertex == tri->tet_vertex) {
            continue;
        }

        v1 = (int) remaining_face[tri->tet_vertex][vertex];
        v2 = (int) remaining_face[vertex][tri->tet_vertex];

        for (i = 0; i < net_flow_around_vertex(tri, vertex); i++) {
            distance[vertex]                = i;
            distance[v1]                    = net_flow_around_vertex(tri, v1)
                                              + (net_flow_around_vertex(tri, vertex) - distance[vertex]);
            distance[v2]                    = net_flow_around_vertex(tri, v2)
                                              + (net_flow_around_vertex(tri, vertex) - distance[vertex]);
            distance[tri->tet_vertex]       = -1;

            adj_triangle[vertex]            = 0;
            adj_triangle[v1]                = 1;
            adj_triangle[v2]                = 1;
            adj_triangle[tri->tet_vertex]   = 0;

            set_cusp_region_data(cusp, tri, distance, adj_triangle, index);
            index++;
        }

    }

    // center region
    for (vertex = 0; vertex < 4; vertex++) {
        if (vertex == tri->tet_vertex)
            continue;

        distance[vertex]        = net_flow_around_vertex(tri, vertex);
        adj_triangle[vertex]    = 1;
    }

    distance[tri->tet_vertex]       = -1;
    adj_triangle[tri->tet_vertex]   = 0;

    set_cusp_region_data(cusp, tri, distance, adj_triangle, index);
    index++;
    return index;
}

/*
 * Helper function to init_cusp_regions which allocates the attributes of the
 * cusp region
 */

void set_cusp_region_data(CuspStructure *cusp, CuspTriangle *tri, const int distance[4],
                          const Boolean adj_cusp_triangle[4], int index) {
    int i, j, v1, v2, v3;
    CuspRegion *region = NEW_STRUCT( CuspRegion );
    INSERT_BEFORE(region, &cusp->cusp_region_end[TRI_TO_INDEX(tri->tet_index, tri->tet_vertex)])

    region->tri             = tri;
    region->tet_index       = region->tri->tet_index;
    region->tet_vertex      = region->tri->tet_vertex;
    region->index           = index;

    // default values
    for (i = 0; i < 4; i++) {
        region->adj_cusp_triangle[i] = FALSE;
        region->adj_cusp_regions[i]  = NULL;

        for (j = 0; j < 4; j++) {
            region->curve[i][j]             = -1;
            region->dive[i][j]              = 0;
            region->num_adj_curves[i][j]    = 0;
            region->temp_adj_curves[i][j]   = 0;
        }
    }

    for (i = 0; i < 3; i++) {
        v1 = edgesThreeToFour[tri->tet_vertex][i];
        v2 = edgesThreeToFour[tri->tet_vertex][(i + 1) % 3];
        v3 = edgesThreeToFour[tri->tet_vertex][(i + 2) % 3];

        region->curve[v2][v1]   = distance[v1];
        region->curve[v3][v1]   = distance[v1];
        region->dive[v2][v1]    = distance[v1] ? FALSE : TRUE;
        region->dive[v3][v1]    = distance[v1] ? FALSE : TRUE;

        region->adj_cusp_triangle[v1] = adj_cusp_triangle[v1];
    }
}
/*
 * Calculate which regions are located across cusp edges and store the result
 * in the adj_cusp_regions attribute
 */

void update_adj_region_data(CuspStructure *cusp) {
    CuspTriangle *adj_triangle;
    CuspRegion *region;
    FaceIndex f;
    int i, adj_index;

    // Add adjacent region info
    for (i = 0; i < 4 * cusp->manifold->num_tetrahedra; i++) {
        for (region = cusp->cusp_region_begin[i].next; region != &cusp->cusp_region_end[i]; region = region->next) {
            for (f = 0; f < 4; f++) {
                if (!region->adj_cusp_triangle[f] || region->tet_vertex == f) {
                    region->adj_cusp_regions[f] = NULL;
                    continue;
                }

                adj_triangle = region->tri->neighbours[f];
                adj_index = TRI_TO_INDEX(adj_triangle->tet_index, adj_triangle->tet_vertex);
                region->adj_cusp_regions[f] = find_adj_region(&cusp->cusp_region_begin[adj_index],
                                                              &cusp->cusp_region_end[adj_index],
                                                              region, f);
            }
        }
    }
}

/*
 * Find the cusp region which is adjacent to x across face.
 */

CuspRegion *find_adj_region(CuspRegion *cusp_region_begin, CuspRegion *cusp_region_end,
                            CuspRegion *x, int face) {
    int v1, v2, y_vertex1, y_vertex2, y_face, distance_v1, distance_v2, tet_index, tet_vertex;
    Boolean adj_face;
    CuspTriangle *tri = x->tri;
    CuspRegion *region;

    v1 = (int) remaining_face[tri->tet_vertex][face];
    v2 = (int) remaining_face[face][tri->tet_vertex];

    y_vertex1    = EVALUATE(tri->tet->gluing[face], v1);
    y_vertex2    = EVALUATE(tri->tet->gluing[face], v2);
    y_face       = EVALUATE(tri->tet->gluing[face], face);

    // Check current adj region first
    if (x->adj_cusp_regions[face] != NULL) {
        distance_v1      = (x->curve[face][v1] == x->adj_cusp_regions[face]->curve[y_face][y_vertex1]);
        distance_v2      = (x->curve[face][v2] == x->adj_cusp_regions[face]->curve[y_face][y_vertex2]);
        adj_face         = x->adj_cusp_regions[face]->adj_cusp_triangle[y_face];

        if (distance_v1 && distance_v2 && adj_face)
            return x->adj_cusp_regions[face];
    }

    /*
     * We search through the regions in reverse as the new regions
     * are added to the end of the doubly linked list
     */
    for (region = cusp_region_end->prev; region != cusp_region_begin; region = region->prev) {
        tet_index    = (tri->neighbours[face]->tet_index == region->tet_index);
        tet_vertex   = (tri->neighbours[face]->tet_vertex == region->tet_vertex);

        if (!tet_index || !tet_vertex)
            continue;

        distance_v1      = (x->curve[face][v1] == region->curve[y_face][y_vertex1]);
        distance_v2      = (x->curve[face][v2] == region->curve[y_face][y_vertex2]);
        adj_face         = region->adj_cusp_triangle[y_face];

        // missing distance
        if (region->curve[y_face][y_vertex1] == -1 || region->curve[y_face][y_vertex2] == -1)
            uFatalError("find_adj_region", "symplectic_basis");

        if (distance_v1 && distance_v2 && adj_face)
            return region;
    }

    // We didn't find a cusp region
    //uFatalError("find_cusp_region", "symplectic_basis");
    return NULL;
}

/*
 * region1 splits into region1 and region2, set them up to be split
 */

void copy_region(CuspRegion *region1, CuspRegion *region2) {
    int i, j;

    if (region1 == NULL || region2 == NULL || region1->tri == NULL)
        uFatalError("copy_region", "symplectic_basis");

    region2->tri            = region1->tri;
    region2->tet_index      = region1->tet_index;
    region2->tet_vertex     = region1->tet_vertex;

    for (i = 0; i < 4; i++) {
        region2->adj_cusp_triangle[i]   = region1->adj_cusp_triangle[i];
        region2->adj_cusp_regions[i]    = NULL;

        for (j = 0; j < 4; j++) {
            region2->curve[i][j]            = region1->curve[i][j];
            region2->dive[i][j]             = FALSE;
            region2->num_adj_curves[i][j]   = region1->num_adj_curves[i][j];
            region2->temp_adj_curves[i][j]  = region1->temp_adj_curves[i][j];
        }
    }
}

int net_flow_around_vertex(CuspTriangle *tri, int vertex) {
    int mflow, lflow, retval;

    // Contribution from meridian curves
    mflow = FLOW(tri->tet->curve[M][right_handed][tri->tet_vertex][remaining_face[tri->tet_vertex][vertex]],
                 tri->tet->curve[M][right_handed][tri->tet_vertex][remaining_face[vertex][tri->tet_vertex]]);

    // Contribution from longitudinal curves
    lflow = FLOW(tri->tet->curve[L][right_handed][tri->tet_vertex][remaining_face[tri->tet_vertex][vertex]],
                 tri->tet->curve[L][right_handed][tri->tet_vertex][remaining_face[vertex][tri->tet_vertex]]);

    retval = ABS(mflow) + ABS(lflow);
    return retval;
}
