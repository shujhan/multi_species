#include "AMRStructure.hpp"

// =====================================================================
//  Shared-edge point lookup for directional (x / v / xv) refinement.
//
//  When a panel splits it creates points on the two edges perpendicular
//  to the split:
//
//      v-split  ->  new points on the LEFT and RIGHT edges
//      x-split  ->  new points on the BOTTOM and TOP edges
//      xv-split ->  new points on all four edges
//
//  Each of those points may already exist on the neighbour's side of the
//  edge.  In the isotropic scheme that was a single test, because one
//  tree level meant one level of refinement in both directions, so the
//  neighbour either had the point or it didn't.
//
//  With directional splits that is no longer true.  An x-split does not
//  subdivide a vertical edge at all, so a neighbour may be several tree
//  levels deep and still not have reached the v-level we are asking
//  about -- or it may have reached it below an x-split that we have to
//  step through first.  So the test becomes a short loop: descend toward
//  our edge, checking at each step whether that panel's facing edge
//  carries a point at the coordinate we want.
//
//  Everything here is decided from coordinates, not from stored levels,
//  so it stays correct no matter how large the level difference is.
//
//  Point ordering within a panel (column-major, index = 3*i + j):
//
//        2 ----- 5 ----- 8        left  edge : 0, 1, 2
//        |               |        right edge : 6, 7, 8
//        1 ----- 4 ----- 7        bottom edge: 0, 3, 6
//        |               |        top   edge : 2, 5, 8
//        0 ----- 3 ----- 6
//
//  Child ordering:
//      xv-split (4) : [0]=BL [1]=TL [2]=BR [3]=TR
//      v-split  (2) : [0]=bottom [1]=top
//      x-split  (2) : [0]=left   [1]=right
// =====================================================================


// ---------------------------------------------------------------------
//  Coordinate accessors that see both committed and staged points.
// ---------------------------------------------------------------------
double AMRStructure::point_x(int ind) const {
    if (ind < (int) xs.size()) { return xs[ind]; }
    return (*staged_xs)[ind - (int) xs.size()];
}
double AMRStructure::point_v(int ind) const {
    if (ind < (int) ps.size()) { return ps[ind]; }
    return (*staged_ps)[ind - (int) ps.size()];
}


// ---------------------------------------------------------------------
//  Read / write a neighbour slot by edge id.
// ---------------------------------------------------------------------
int AMRStructure::edge_nbr_ind(const Panel& P, int edge) const {
    switch (edge) {
        case EDGE_LEFT   : return P.left_nbr_ind;
        case EDGE_RIGHT  : return P.right_nbr_ind;
        case EDGE_BOTTOM : return P.bottom_nbr_ind;
        default          : return P.top_nbr_ind;
    }
}

void AMRStructure::set_edge_nbr(Panel& P, int edge, int val) {
    switch (edge) {
        case EDGE_LEFT   : P.left_nbr_ind   = val; break;
        case EDGE_RIGHT  : P.right_nbr_ind  = val; break;
        case EDGE_BOTTOM : P.bottom_nbr_ind = val; break;
        default          : P.top_nbr_ind    = val; break;
    }
}


// ---------------------------------------------------------------------
//  Does Q's facing edge carry a point at `target`?
//
//  "Facing" means the edge of Q that touches us:
//      our LEFT   edge  <->  Q's RIGHT  edge  (points 6, 7, 8)
//      our RIGHT  edge  <->  Q's LEFT   edge  (points 0, 1, 2)
//      our BOTTOM edge  <->  Q's TOP    edge  (points 2, 5, 8)
//      our TOP    edge  <->  Q's BOTTOM edge  (points 0, 3, 6)
//
//  `target` is a v-coordinate for the vertical edges and an x-coordinate
//  for the horizontal ones.  Returns the point index, or -1.
// ---------------------------------------------------------------------
int AMRStructure::facing_edge_hit(const Panel& Q, int edge, double target) const {
    int a, b, c;
    bool compare_v;
    switch (edge) {
        case EDGE_LEFT   : a = 6; b = 7; c = 8; compare_v = true;  break;
        case EDGE_RIGHT  : a = 0; b = 1; c = 2; compare_v = true;  break;
        case EDGE_BOTTOM : a = 2; b = 5; c = 8; compare_v = false; break;
        default          : a = 0; b = 3; c = 6; compare_v = false; break;
    }
    // Coordinates on both sides are produced by repeated halving from the
    // same origin, so they agree to a few ULP; this tolerance is many
    // orders of magnitude below the smallest panel at any usable
    // max_height, and many orders above the rounding.
    const double tol = 1e-12 * (compare_v ? Lp : Lx);

    const int ii[3] = {a, b, c};
    for (int k = 0; k < 3; ++k) {
        int pind = Q.point_inds[ii[k]];
        double coord = compare_v ? point_v(pind) : point_x(pind);
        if (fabs(coord - target) <= tol) { return pind; }
    }
    return -1;
}


// ---------------------------------------------------------------------
//  Step one level down the neighbour subtree, toward our edge.
//
//  Two independent moves, applied together:
//    - move toward us  : the direction normal to the shared edge.  A left
//                        neighbour hands back its right-hand child, a
//                        bottom neighbour its top child, and so on.
//    - move toward the target : the direction along the shared edge.  Pick
//                        the child whose span brackets `target`.
//
//  A split in only the normal direction does not change the neighbour's
//  refinement in the governing direction -- it is a step we take without
//  making progress on the question, which is exactly why this has to be a
//  loop.  Returns -1 when Q is a leaf.
// ---------------------------------------------------------------------
int AMRStructure::child_toward(const Panel& Q, int edge, double target) const {
    if (Q.is_leaf()) { return -1; }
    const int cs = Q.child_inds_start;
    const double xmid = point_x(Q.point_inds[4]);
    const double vmid = point_v(Q.point_inds[4]);

    switch (edge) {
        case EDGE_LEFT :   // want Q's right side, select by v
            if (Q.is_refined_xp) { return cs + 2 + (target < vmid ? 0 : 1); }
            if (Q.is_refined_p)  { return cs +     (target < vmid ? 0 : 1); }
            return cs + 1;                       // x-split: right child
        case EDGE_RIGHT :  // want Q's left side, select by v
            if (Q.is_refined_xp) { return cs +     (target < vmid ? 0 : 1); }
            if (Q.is_refined_p)  { return cs +     (target < vmid ? 0 : 1); }
            return cs + 0;                       // x-split: left child
        case EDGE_BOTTOM : // want Q's top side, select by x
            if (Q.is_refined_xp) { return cs + 1 + (target < xmid ? 0 : 2); }
            if (Q.is_refined_x)  { return cs +     (target < xmid ? 0 : 1); }
            return cs + 1;                       // v-split: top child
        default :          // EDGE_TOP: want Q's bottom side, select by x
            if (Q.is_refined_xp) { return cs +     (target < xmid ? 0 : 2); }
            if (Q.is_refined_x)  { return cs +     (target < xmid ? 0 : 1); }
            return cs + 0;                       // v-split: bottom child
    }
}


// ---------------------------------------------------------------------
//  The whole question, in one place.
//
//  Returns the index of an existing point at `target` on the given edge
//  of panel `panel_ind`, or -1 if the caller must create one.  Creating
//  one is always a legal answer: interpolation is per-panel over 9 points
//  and the quadrature is a sum of per-leaf panel integrals, so a
//  duplicated point costs memory and nothing else.  That is what makes
//  this routine safe to get conservatively wrong.
//
//  `found_panel_ind` receives the abutting neighbour panel, for use as
//  the new child's neighbour pointer, or -1.
// ---------------------------------------------------------------------
int AMRStructure::lookup_edge_point(int panel_ind, int edge, double target,
                                    int& found_panel_ind)
{
    found_panel_ind = -1;

    const Panel& P = panels[panel_ind];
    int q = edge_nbr_ind(P, edge);

    if (q < 0)          { return -1; }  // -2 boundary, -1 coarser: create
    if (q == panel_ind) { return -1; }  // periodic self-neighbour (1 panel in x)

    // Across the periodic seam the two sides sit a full domain apart in x,
    // so the points are distinct even though the panels touch.
    if (bcs == periodic_bcs) {
        if (edge == EDGE_LEFT  && P.is_left_bdry ) { return -1; }
        if (edge == EDGE_RIGHT && P.is_right_bdry) { return -1; }
    }

    // Descend.  The bound is generous; it exists only so that a corrupted
    // tree produces a missing point rather than an infinite loop.
    const int max_steps = 2 * max_height + 8;
    for (int step = 0; step < max_steps; ++step) {
        const Panel& Q = panels[q];

        int hit = facing_edge_hit(Q, edge, target);
        if (hit >= 0) { found_panel_ind = q; return hit; }

        int c = child_toward(Q, edge, target);
        if (c < 0) { return -1; }       // leaf, and it has no point there
        q = c;
    }
    return -1;
}