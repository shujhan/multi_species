#include "AMRStructure.hpp"


/* new features for refining in v
---
in generate_mesh:

if num_p_levels > 0 : 
    create_preerefined_mesh_plus_v
else : 
    Create_prerefined_mesh

refine_panels : no change

refine_panels_in_v :

*/

// #define DEBUG
// #define DEBUG_L2

int AMRStructure::create_prerefined_mesh_p_refinement() {
    // printf("setting initial mesh of height %i.\n", initial_height);
    if (initial_height + p_height < 1) {
        throw std::invalid_argument("height + p_height must be greater than 1");
    }
    double dx = (x_max - x_min) / 4;
    double dp = (p_max - p_min) / 4;
    std::vector<double> xs_init, ps_init;
    for(int ii = 0; ii < 5; ++ii) {
        xs_init.push_back(x_min + ii * dx);
        ps_init.push_back(p_min + ii * dp);
    }
    panels.clear();
    xs.clear();
    ps.clear();
    fs.clear();
    xs.reserve(15);
    ps.reserve(15);
    for (int ii = 0; ii < 5; ii += 2) {
        for (int jj = 0; jj < 5; jj+=2) { 
            xs.push_back(xs_init[ii]); 
            ps.push_back(ps_init[jj]);
        }
    }

    for (int ii = 0; ii < 5; ii += 2) {
        for (int jj = 1; jj < 5; jj += 2) {
            xs.push_back(xs_init[ii]); 
            ps.push_back(ps_init[jj]);
        }
    }
    // for (int ii = 0; ii < 2; ++ii) {
    //     xs.push_back(xs_init[2*ii]); xs.push_back(xs_init[2*ii]);
    //     ps.push_back(ps_init[1]); ps.push_back(ps_init[3]);
    //     for (int jj = 0; jj < 5; ++jj) {
    //         xs.push_back(xs_init[1 + 2*ii]);
    //         ps.push_back(ps_init[jj]);
    //     }
    // }
    // for (int jj = 1; jj < 5; jj+=2) {
    //     xs.push_back(xs_init[4]);
    //     ps.push_back(ps_init[jj]);
    // }
    fs = std::vector<double>(15, 1.0);
    //   2        5         8 (2 by periodic bcs)
    
    //   10      12 [2]    14 (10 by pbcs)
    
    //   1        4         7 (1 by periodic bcs)
    
    //   9       11 [1]    13 (9 by pbcs)
    
    //   0        3         6 (0 by periodic bcs)
    if (bcs == periodic_bcs) {
        panels.push_back(Panel{});
        panels[0].is_left_bdry = true;
        panels[0].is_right_bdry = true;
        panels.push_back(Panel{1, 1, 0, 0, 1, 2, 1, -2});
        panels[1].is_left_bdry = true;
        panels[1].is_right_bdry = true;
        panels.push_back(Panel{2, 1, 0, 1, 2, -2, 2, 1});
        panels[2].is_left_bdry = true;
        panels[2].is_right_bdry = true;
    } else if (bcs == open_bcs) {
        panels.push_back(Panel{0, 0, -1, -1, -2, -2, -2, -2});
        panels[0].set_point_inds(0,1,2,3,4,5,6,7,8);
        panels.push_back(Panel{1, 1, 0, 0, -2, 2, -2, -2});
        panels.push_back(Panel{2, 1, 0, 1, -2, -2, -2, 1});
    }
    panels[1].set_point_inds(0,9,1,3,11,4,6,13,7);
    panels[1].needs_refinement = true;
    panels[2].set_point_inds(1,10,2,4,12,5,7,14,8);
    panels[2].needs_refinement = true;

    bool is_refined_p = true;
    panels[0].set_child_inds_start(1,is_refined_p);
    minimum_unrefined_index = 1;

    // call refine
    for (int level = 1; level < p_height; ++level) {
        int num_panels_pre_refine = panels.size();

        // route the uniform v pre-refinement through the same dispatch the
        // adaptive passes use, so there is only one v-split implementation
        for (auto panel_it = panels.begin() + minimum_unrefined_index; panel_it != panels.end(); ++panel_it) {
            panel_it->needs_refinement = true;
            panel_it->needs_refine_x  = false;
            panel_it->needs_refine_v  = true;
        }
        bool do_adaptive_refine = false;
        refine_panels( [] (double x, double v) {return 1.0;} , do_adaptive_refine);
        minimum_unrefined_index = num_panels_pre_refine;

    }
    for (int level = 0; level < initial_height; ++level) {
        int num_panels_pre_refine = panels.size();

        for (auto panel_it = panels.begin() + minimum_unrefined_index; panel_it != panels.end(); ++panel_it) {
            panel_it->needs_refinement = true;
            panel_it->needs_refine_x  = true;
            panel_it->needs_refine_v  = true;
        }
        refine_panels( [] (double x, double v) {return 1.0;} , false);
        minimum_unrefined_index = num_panels_pre_refine;
    }
    is_initial_mesh_set = true;

    return 0;
}


int AMRStructure::create_prerefined_mesh() {
    // printf("setting initial mesh of height %i.\n", initial_height);
    if (initial_height < 1) {
        throw std::invalid_argument("height must be greater than 1");
    }
    double dx = (x_max - x_min) / 4;
    double dp = (p_max - p_min) / 4;
    std::vector<double> xs_init, ps_init;
    for(int ii = 0; ii < 5; ++ii) {
        xs_init.push_back(x_min + ii * dx);
        ps_init.push_back(p_min + ii * dp);
    }
    panels.clear();
    xs.clear();
    ps.clear();
    fs.clear();
    xs.reserve(25);
    ps.reserve(25);
    for (int ii = 0; ii < 5; ii += 2) {
        for (int jj = 0; jj < 5; jj+=2) { 
            xs.push_back(xs_init[ii]); 
            ps.push_back(ps_init[jj]);
        }
    }
    for (int ii = 0; ii < 2; ++ii) {
        xs.push_back(xs_init[2*ii]); xs.push_back(xs_init[2*ii]);
        ps.push_back(ps_init[1]); ps.push_back(ps_init[3]);
        for (int jj = 0; jj < 5; ++jj) {
            xs.push_back(xs_init[1 + 2*ii]);
            ps.push_back(ps_init[jj]);
        }
    }
    for (int jj = 1; jj < 5; jj+=2) {
        xs.push_back(xs_init[4]);
        ps.push_back(ps_init[jj]);
    }
    fs = std::vector<double>(25, 1.0);
    //   2   15      5    22    8 (2 by periodic bcs)
    
    //   10  14[2]   17   21[4] 24 (10 by pbcs)
    
    //   1   13      4    20     7 (1 by periodic bcs)
    
    //   9   12[1]   16   19[3] 23 (9 by pbcs)
    
    //   0   11     3    18     6 (0 by periodic bcs)
    if (bcs == periodic_bcs) {
        panels.push_back(Panel{});
        panels[0].is_left_bdry = true;
        panels[0].is_right_bdry = true;
        panels.push_back(Panel{1, 1, 0, 0, 3, 2, 3, -2});
        panels[1].is_left_bdry = true;
        panels.push_back(Panel{2, 1, 0, 1, 4, -2, 4, 1});
        panels[2].is_left_bdry = true;
        panels.push_back(Panel{3, 1, 0, 2, 1, 4, 1, -2});
        panels[3].is_right_bdry = true;
        panels.push_back(Panel{4, 1, 0, 3, 2,-2,2,3});
        panels[4].is_right_bdry = true;
    } else if (bcs == open_bcs) {
        panels.push_back(Panel{0, 0, -1, -1, -2, -2, -2, -2});
        panels[0].set_point_inds(0,1,2,3,4,5,6,7,8);
        panels.push_back(Panel{1, 1, 0, 0, -2, 2, 3, -2});
        panels.push_back(Panel{2, 1, 0, 1, -2, -2, 4, 1});
        panels.push_back(Panel{3, 1, 0, 2, 1, 4, -2, -2});
        panels.push_back(Panel{4, 1, 0, 3, 2,-2,-2,3});
    }
    panels[1].set_point_inds(0,9,1,11,12,13,3,16,4);
    panels[1].needs_refinement = true;
    panels[2].set_point_inds(1,10,2,13,14,15,4,17,5);
    panels[2].needs_refinement = true;
    panels[3].set_point_inds(3,16,4,18,19,20,6,23,7);
    panels[3].needs_refinement = true;
    panels[4].set_point_inds(4,17,5,20,21,22,7,24,8);
    panels[4].needs_refinement = true;

    panels[0].set_child_inds_start(1);
    minimum_unrefined_index = 1;

    // call refine
    for (int level = 1; level < initial_height; ++level) {
        int num_panels_pre_refine = panels.size();

        for (auto panel_it = panels.begin() + minimum_unrefined_index; panel_it != panels.end(); ++panel_it) {
            panel_it->needs_refinement = true;
            panel_it->needs_refine_x  = true;
            panel_it->needs_refine_v  = true;
        }
        refine_panels( [] (double x, double v) {return 1.0;} , false);
        minimum_unrefined_index = num_panels_pre_refine;
    }
    for (int level = 0; level < p_height; ++level) {
        int num_panels_pre_refine = panels.size();

        // route the uniform v pre-refinement through the same dispatch the
        // adaptive passes use, so there is only one v-split implementation
        for (auto panel_it = panels.begin() + minimum_unrefined_index; panel_it != panels.end(); ++panel_it) {
            panel_it->needs_refinement = true;
            panel_it->needs_refine_x  = false;
            panel_it->needs_refine_v  = true;
        }
        bool do_adaptive_refine = false;
        refine_panels( [] (double x, double v) {return 1.0;} , do_adaptive_refine);
        minimum_unrefined_index = num_panels_pre_refine;
    }


    is_initial_mesh_set = true;

    return 0;
}


void AMRStructure::refine_panels_refine_v(std::function<double (double,double)> f, bool do_adaptive_refine) {

    // Note: this assumes that we are refining in v uniformly before any xp refinement;
    // No compatibility with xp refined panels is guaranteed
    std::vector <double> new_xs;
    std::vector <double> new_ps;
    std::vector <double> new_fs;
    std::vector <int> prospective_leaf_inds;
    // int new_vert_ind = particles.size();
    int new_vert_ind = xs.size();
    int num_panels_before_this_iter = panels.size();
    


    for (int jj = minimum_unrefined_index; jj < num_panels_before_this_iter; ++jj) {
        Panel* panel= &(panels[jj]);
        
        if (panel->needs_refinement ) {
            std::vector<double> panel_xs;
            std::vector<double> panel_ps;
            double dx, dp;

            const int* panel_points = panel->point_inds;
            for (int ii = 0; ii < 9; ++ii) {
                int point_ind = panel_points[ii];
                panel_xs.push_back(xs[point_ind]);
                panel_ps.push_back(ps[point_ind]);
            }
            dx = panel_xs[3] - panel_xs[0];
            dp = panel_ps[1] - panel_ps[0];
            double sub_dp = 0.5 * dp;
            double sub_dx = 0.5 * dx;

            int num_new_panels = panels.size();
            double subpanel_xs[5], subpanel_ps[5];

            for (int ii = 0; ii < 5; ii ++) {
                subpanel_ps[ii] = panel_ps[0] + sub_dp * ii;
                subpanel_xs[ii] = panel_xs[0] + sub_dx * ii;
            }
            //   2        5         8 (2 by periodic bcs)
            
            //   10      12 [2]    14 (10 by pbcs)
            
            //   1        4         7 (1 by periodic bcs)
            
            //   9       11 [1]    13 (9 by pbcs)
            
            //   0        3         6 (0 by periodic bcs)

            int point_9_ind, point_10_ind, point_13_ind, point_14_ind;
            int child_0_bottom_nbr_ind = -1;
            int child_0_left_nbr_ind = num_new_panels; 
            int child_0_right_nbr_ind = child_0_left_nbr_ind;
            int child_1_left_nbr_ind = num_new_panels + 1;
            int child_1_top_nbr_ind = -1;
            int child_1_right_nbr_ind = child_1_left_nbr_ind;

            // generate new vertices

            Panel* panel_parent;
            // check left neighbor
            if (panel->left_nbr_ind == -2) {
                child_0_left_nbr_ind = -2;
                child_1_left_nbr_ind = -2;
                point_9_ind = new_vert_ind++;
                point_10_ind = new_vert_ind++;
                new_xs.push_back(subpanel_xs[0]); new_xs.push_back(subpanel_xs[0]);
                new_ps.push_back(subpanel_ps[1]); new_ps.push_back(subpanel_ps[3]);
            } else if (panel->left_nbr_ind == -1) {
                panel_parent = &(panels[panel->parent_ind]);
                Panel* parent_left = &(panels[panel_parent->left_nbr_ind]);
                if (! (parent_left->is_refined_xp || parent_left->is_refined_p) ) {
                    parent_left->needs_refinement = true;
                    need_further_refinement = true;
                    // cout << "refine: setting refinement flag in panel " << jj << endl;
                }
                point_9_ind = new_vert_ind;
                // point_10_ind = new_vert_ind++;
                point_10_ind = point_9_ind + 1;
                new_vert_ind += 2;
                new_xs.push_back(subpanel_xs[0]); new_xs.push_back(subpanel_xs[0]);
                new_ps.push_back(subpanel_ps[1]); new_ps.push_back(subpanel_ps[3]);
            } else {
                Panel* panel_left = &(panels[panel->left_nbr_ind]);
                if (! (panel_left->is_refined_xp || panel_left->is_refined_p) ) {
                    point_9_ind = new_vert_ind++;
                    point_10_ind = new_vert_ind++;
                    new_xs.push_back(subpanel_xs[0]); new_xs.push_back(subpanel_xs[0]);
                    new_ps.push_back(subpanel_ps[1]); new_ps.push_back(subpanel_ps[3]);
                }
                else {
                    if (panel_left->is_refined_xp) {
                        child_0_left_nbr_ind = panel_left->child_inds_start +2;
                        child_1_left_nbr_ind = panel_left->child_inds_start + 3;
                    } else { // panel_left is refined in v
                        child_0_left_nbr_ind = panel_left->child_inds_start +0;
                        child_1_left_nbr_ind = panel_left->child_inds_start + 1;
                    }
                    Panel* child_0_left_nbr = &(panels[child_0_left_nbr_ind]);
                    child_0_left_nbr->right_nbr_ind = num_new_panels;
                    Panel* child_1_left_nbr = &(panels[child_1_left_nbr_ind]);
                    child_1_left_nbr->right_nbr_ind = num_new_panels + 1;
                    if (panel->is_left_bdry && bcs==periodic_bcs) {
                        point_9_ind = new_vert_ind++;
                        point_10_ind = new_vert_ind++;
                        new_xs.push_back(subpanel_xs[0]); new_xs.push_back(subpanel_xs[0]);
                        new_ps.push_back(subpanel_ps[1]); new_ps.push_back(subpanel_ps[3]);
                    } else {
                        point_9_ind = child_0_left_nbr->point_inds[7];
                        point_10_ind = child_1_left_nbr->point_inds[7];
                    }
                }
            }
            
            // check bottom neighbor
            int bottom_nbr_ind = panel->bottom_nbr_ind;
            if (bottom_nbr_ind == -2) {
                child_0_bottom_nbr_ind = -2;
            } else if (bottom_nbr_ind == -1) {
                cout << "not allowed to refine in v if panel doesn't have bottom neighbor!" << endl;
            } else {
                Panel* panel_bottom = &(panels[bottom_nbr_ind]);
                if (panel_bottom->is_refined_xp ) {
                    cout << "Shouldn't be allowed to call refine in v if bottom neighbor is refined in x and v!" << endl;
                }
                else {
                    if (!panel_bottom->is_refined_p) {
                        child_0_bottom_nbr_ind = bottom_nbr_ind;
                        panel_bottom->top_nbr_ind = num_new_panels;
                    } else { //panel_bottom is refined in v
                        child_0_bottom_nbr_ind = panel_bottom->child_inds_start + 1;
                        Panel* child_0_bottom_nbr = &(panels[child_0_bottom_nbr_ind]);
                        child_0_bottom_nbr->top_nbr_ind = num_new_panels;
                    }
                }
            }

            // check top neighbor
            int top_nbr_ind = panel->top_nbr_ind;
            if (top_nbr_ind == -2) {
                child_1_top_nbr_ind = -2;
            } else if (top_nbr_ind == -1) {
                cout << "not allowed to refine in v if panel doesn't have top neighbor!" << endl;
            } else {
                Panel* panel_top = &(panels[top_nbr_ind]);
                if (panel_top->is_refined_xp ) {
                    cout << "Shouldn't be allowed to call refine in v if bottom neighbor is refined in x and v!" << endl;
                } else {
                    if (!panel_top->is_refined_p) {
                        child_1_top_nbr_ind = -1;
                        // panel_top->bottom_nbr_ind = num_new_panels + 1;
                    }
                    else {
                        child_1_top_nbr_ind = panel_top->child_inds_start;
                        Panel* child_1_top_nbr = &(panels[child_1_top_nbr_ind]);
                        child_1_top_nbr->bottom_nbr_ind = num_new_panels + 1;
                    }
                }
            }

            // check right neighbor
            if (panel->right_nbr_ind == -2) {
                child_0_right_nbr_ind = -2;
                child_1_right_nbr_ind = -2;
                point_13_ind = new_vert_ind++;
                point_14_ind = new_vert_ind++;
                new_xs.push_back(subpanel_xs[4]); new_xs.push_back(subpanel_xs[4]);
                new_ps.push_back(subpanel_ps[1]); new_ps.push_back(subpanel_ps[3]);
            } else if (panel->right_nbr_ind == -1) {
                panel_parent = &(panels[panel->parent_ind]);
                Panel* parent_right = &(panels[panel_parent->right_nbr_ind]);
                if (!(parent_right->is_refined_xp || parent_right->is_refined_p) ) {
                    parent_right->needs_refinement = true;
                    need_further_refinement = true;
                    // cout << "refine: setting refinement flag in panel " << jj << endl;
                }
                point_13_ind = new_vert_ind++;
                point_14_ind = new_vert_ind++;
                new_xs.push_back(subpanel_xs[4]); new_xs.push_back(subpanel_xs[4]);
                new_ps.push_back(subpanel_ps[1]); new_ps.push_back(subpanel_ps[3]);
            } else {
                Panel* panel_right = &(panels[panel->right_nbr_ind]);
                if (! (panel_right->is_refined_xp || panel_right->is_refined_p) ) {
                    point_13_ind = new_vert_ind++;
                    point_14_ind = new_vert_ind++;
                    new_xs.push_back(subpanel_xs[4]); new_xs.push_back(subpanel_xs[4]);
                    new_ps.push_back(subpanel_ps[1]); new_ps.push_back(subpanel_ps[3]);
                }
                else {
                    child_0_right_nbr_ind = panel_right->child_inds_start;
                    child_1_right_nbr_ind = panel_right->child_inds_start + 1;
                    Panel* child_0_right_nbr = &(panels[child_0_right_nbr_ind]);
                    child_0_right_nbr->left_nbr_ind = num_new_panels;
                    Panel* child_1_right_nbr = &(panels[child_1_right_nbr_ind]);
                    child_1_right_nbr->left_nbr_ind = num_new_panels + 1;
                    if (panel->is_right_bdry && bcs==periodic_bcs) {
                        point_13_ind = new_vert_ind++;
                        point_14_ind = new_vert_ind++;
                        new_xs.push_back(subpanel_xs[4]); new_xs.push_back(subpanel_xs[4]);
                        new_ps.push_back(subpanel_ps[1]); new_ps.push_back(subpanel_ps[3]);
                    } else {
                        point_13_ind = child_0_right_nbr->point_inds[1];
                        point_14_ind = child_1_right_nbr->point_inds[1];
                    }
                }
            } // end check right neighbor

            // add interior points
            int point_11_ind = new_vert_ind;
            for (int ii = 0; ii < 2; ii++) {
                new_xs.push_back(subpanel_xs[2]);
                new_ps.push_back(subpanel_ps[1+2*ii]);
            }
            new_vert_ind += 2;

            // generate new panels
            // add these to list of prospective_panel_indices
            if (do_adaptive_refine) {
                for (int ii = num_new_panels; ii < num_new_panels + 4; ++ii) {
                    prospective_leaf_inds.push_back(ii);
                }
            }
            // panel->child_inds_start = num_new_panels;
            bool refining_in_p = true;
            panel->set_child_inds_start(num_new_panels, refining_in_p);
            // printf("post refinement, panel looks like:\n");
            // panel->print_panel();
            int child_level = panel->level + 1;
            int panel_ind = panel->panel_ind;
            int* point_inds = panel->point_inds;
            // for (int ii = 0; ii < ; ii++) {
            //     panel_vertex_inds[ii] = panel->point_inds[ii];
            // }
            panels.push_back(Panel {num_new_panels, child_level, panel_ind, 0, 
                    point_inds[0], point_9_ind, point_inds[1],
                    point_inds[3], point_11_ind, point_inds[4],
                    point_inds[6], point_13_ind, point_inds[7],
                    child_0_left_nbr_ind, num_new_panels + 1, 
                    child_0_right_nbr_ind, child_0_bottom_nbr_ind,
                    panel->is_left_bdry, panel->is_right_bdry});
            panels.push_back(Panel {num_new_panels+1, child_level, panel_ind, 1, 
                    point_inds[1], point_10_ind, point_inds[2],
                    point_inds[4], point_11_ind+1, point_inds[5],
                    point_inds[7], point_14_ind, point_inds[8],
                    child_1_left_nbr_ind, child_1_top_nbr_ind,
                    child_1_right_nbr_ind, num_new_panels,
                    panel->is_left_bdry, panel->is_right_bdry});

        } // end if panel is flagged
        
    } //end for loop through panels

    // set fs
    new_fs.reserve(new_xs.size() );
    for (int ii = 0; ii < new_xs.size(); ++ii) {
        new_fs.push_back( f(new_xs.at(ii), new_ps.at(ii)) );
    }

    for (int ii = 0; ii < new_xs.size(); ++ii) {
        xs.push_back(new_xs[ii]); ps.push_back(new_ps[ii]); fs.push_back(new_fs[ii]);
        // particles.push_back(Particle(new_xs.at(ii), new_ps.at(ii), new_fs.at(ii), 0.0));
    }


    // test 
    // if (do_adaptive_refine) { 
    //     for (int ii = 0; ii < prospective_leaf_inds.size(); ++ii) {
    //         test_panel(prospective_leaf_inds.at(ii), false);
    //     }
    // }
}

// =====================================================================
//  Directional AMR
//
//  A panel may now be refined in x only (2 children, [0]=left,[1]=right),
//  in v only (2 children, [0]=bottom,[1]=top), or in both (4 children,
//  column-major [0]=BL [1]=TL [2]=BR [3]=TR).  Every child is still a 3x3
//  panel; only its aspect ratio changes.  Any split costs one level.
//
//  Which edges receive new points depends on the split direction:
//      refine in v   -> new points on the LEFT and RIGHT edges
//      refine in x   -> new points on the BOTTOM and TOP edges
//      refine in xv  -> new points on all four edges
//  Only those "live" edges need to be queried for existing points.  The
//  remaining "quiet" edges are not subdivided, so their pointers are
//  inherited from the parent.
//
//  Point sharing on a live edge depends on ONE bit of the neighbour:
//      left/right edges -> is the neighbour refined in v?
//      top/bottom edges -> is the neighbour refined in x?
//  The neighbour's refinement in the other direction only changes WHICH
//  of its children abut the edge, never whether the points exist.
// =====================================================================


// Walk down the subtree of `nbr_ind` along a vertical (left/right) shared
// edge until reaching the leaf that actually abuts that edge, and flag it
// for refinement in v.  Refining in x never subdivides a vertical edge, so
// we only ever descend through x-splits.
//   from_left == true : the neighbour lies to our left, so we want its
//                       right-hand child.
void AMRStructure::flag_v_refinement(int nbr_ind, bool from_left) {
    while (nbr_ind >= 0) {
        Panel* q = &(panels[nbr_ind]);
        if (q->refined_in_v()) { return; }   // the edge points already exist
        if (! q->is_refined_x) {             // leaf: flag it
            q->needs_refinement = true;
            q->needs_refine_v   = true;
            need_further_refinement = true;
            return;
        }
        nbr_ind = q->child_inds_start + (from_left ? 1 : 0);
    }
}

// Mirror image: descend along a horizontal (top/bottom) shared edge through
// v-splits and flag the abutting leaf for refinement in x.
//   from_below == true : the neighbour lies below us, so we want its top child.
void AMRStructure::flag_x_refinement(int nbr_ind, bool from_below) {
    while (nbr_ind >= 0) {
        Panel* q = &(panels[nbr_ind]);
        if (q->refined_in_x()) { return; }
        if (! q->is_refined_p) {
            q->needs_refinement = true;
            q->needs_refine_x   = true;
            need_further_refinement = true;
            return;
        }
        nbr_ind = q->child_inds_start + (from_below ? 1 : 0);
    }
}


// =====================================================================
//  Directional refinement.
//
//  All three routines have the same shape:
//
//    1. copy everything needed out of the parent, because push_back on
//       `panels` invalidates any reference into it;
//    2. for each point on a subdivided edge, ask lookup_edge_point
//       whether the neighbour already has it, and create it if not;
//    3. create the interior points unconditionally -- they are interior
//       to the parent, so nobody else can own them;
//    4. push the children and wire the neighbour pointers.
//
//  Edge bookkeeping follows one rule:
//
//    live edge  (subdivided by this split) -> the child's neighbour is
//        whatever lookup_edge_point landed on, or -1 if it found nothing
//        (which means everything across that edge is coarser than the
//        child);
//    quiet edge (untouched by this split)  -> the child's edge IS the
//        parent's edge, so inherit the parent's pointer verbatim.
//
//  That second line is the whole correction relative to the previous
//  version, which wrote -1 on quiet edges whenever the neighbour was
//  refined in the perpendicular direction.  That claimed "coarser" about
//  a neighbour that was actually finer, and the leaf search then walked
//  off the front of the panel array.
// =====================================================================


// ---------------------------------------------------------------------
//  Split in x and v -> 4 children, [0]=BL [1]=TL [2]=BR [3]=TR.
// ---------------------------------------------------------------------
void AMRStructure::refine_one_xv(int jj,
                                 std::vector<double>& new_xs,
                                 std::vector<double>& new_ps,
                                 int& new_vert_ind,
                                 std::vector<int>& prospective_leaf_inds,
                                 bool do_adaptive_refine)
{
    int pi[9];
    for (int ii = 0; ii < 9; ++ii) { pi[ii] = panels[jj].point_inds[ii]; }
    const double x0 = xs[pi[0]], xm = xs[pi[3]], x1 = xs[pi[6]];
    const double v0 = ps[pi[0]], vm = ps[pi[1]], v1 = ps[pi[2]];
    const double q1 = 0.5 * (v0 + vm), q3 = 0.5 * (vm + v1);
    const double m1 = 0.5 * (x0 + xm), m3 = 0.5 * (xm + x1);

    const int nbrL = panels[jj].left_nbr_ind;
    const int nbrR = panels[jj].right_nbr_ind;
    const int nbrB = panels[jj].bottom_nbr_ind;
    const int nbrT = panels[jj].top_nbr_ind;
    const bool p_left_bdry  = panels[jj].is_left_bdry;
    const bool p_right_bdry = panels[jj].is_right_bdry;
    const int child_level   = panels[jj].level + 1;
    const int panel_ind     = panels[jj].panel_ind;

    const int c0 = (int) panels.size(), c1 = c0 + 1, c2 = c0 + 2, c3 = c0 + 3;

    auto make_point = [&](double px, double pv) -> int {
        new_xs.push_back(px); new_ps.push_back(pv);
        return new_vert_ind++;
    };
    auto edge_point = [&](int edge, double target,
                          double px, double pv, int& found) -> int {
        int ind = lookup_edge_point(jj, edge, target, found);
        return (ind >= 0) ? ind : make_point(px, pv);
    };
    auto child_nbr = [](int parent_slot, int found) -> int {
        return (parent_slot == -2) ? -2 : found;
    };

    int fL0 = -1, fL1 = -1, fR0 = -1, fR1 = -1;
    int fB0 = -1, fB1 = -1, fT0 = -1, fT1 = -1;

    int pL0 = edge_point(EDGE_LEFT,   q1, x0, q1, fL0);
    int pL1 = edge_point(EDGE_LEFT,   q3, x0, q3, fL1);
    int pR0 = edge_point(EDGE_RIGHT,  q1, x1, q1, fR0);
    int pR1 = edge_point(EDGE_RIGHT,  q3, x1, q3, fR1);
    int pB0 = edge_point(EDGE_BOTTOM, m1, m1, v0, fB0);
    int pB1 = edge_point(EDGE_BOTTOM, m3, m3, v0, fB1);
    int pT0 = edge_point(EDGE_TOP,    m1, m1, v1, fT0);
    int pT1 = edge_point(EDGE_TOP,    m3, m3, v1, fT1);

    // interior: 8 points, never shared with anyone
    int i0 = make_point(m1, q1);
    int i1 = make_point(m1, vm);
    int i2 = make_point(m1, q3);
    int i3 = make_point(xm, q1);
    int i4 = make_point(xm, q3);
    int i5 = make_point(m3, q1);
    int i6 = make_point(m3, vm);
    int i7 = make_point(m3, q3);

    int c0L = child_nbr(nbrL, fL0), c1L = child_nbr(nbrL, fL1);
    int c2R = child_nbr(nbrR, fR0), c3R = child_nbr(nbrR, fR1);
    int c0B = child_nbr(nbrB, fB0), c2B = child_nbr(nbrB, fB1);
    int c1T = child_nbr(nbrT, fT0), c3T = child_nbr(nbrT, fT1);

    // single panel across x with periodic bcs: the panel is its own
    // left and right neighbour, so the children wrap onto each other
    if (nbrR == jj) { c2R = c0; c3R = c1; }
    if (nbrL == jj) { c0L = c2; c1L = c3; }

    panels.push_back(Panel {c0, child_level, panel_ind, 0,
            pi[0], pL0,  pi[1],
            pB0,   i0,   i1,
            pi[3], i3,   pi[4],
            c0L, c1, c2, c0B,
            p_left_bdry, false});
    panels.push_back(Panel {c1, child_level, panel_ind, 1,
            pi[1], pL1,  pi[2],
            i1,    i2,   pT0,
            pi[4], i4,   pi[5],
            c1L, c1T, c3, c0,
            p_left_bdry, false});
    panels.push_back(Panel {c2, child_level, panel_ind, 2,
            pi[3], i3,   pi[4],
            pB1,   i5,   i6,
            pi[6], pR0,  pi[7],
            c0, c3, c2R, c2B,
            false, p_right_bdry});
    panels.push_back(Panel {c3, child_level, panel_ind, 3,
            pi[4], i4,   pi[5],
            i6,    i7,   pT1,
            pi[7], pR1,  pi[8],
            c1, c3T, c3R, c2,
            false, p_right_bdry});

    if (fL0 >= 0) { panels[fL0].right_nbr_ind  = c0; }
    if (fL1 >= 0) { panels[fL1].right_nbr_ind  = c1; }
    if (fR0 >= 0) { panels[fR0].left_nbr_ind   = c2; }
    if (fR1 >= 0) { panels[fR1].left_nbr_ind   = c3; }
    if (fB0 >= 0) { panels[fB0].top_nbr_ind    = c0; }
    if (fB1 >= 0) { panels[fB1].top_nbr_ind    = c2; }
    if (fT0 >= 0) { panels[fT0].bottom_nbr_ind = c1; }
    if (fT1 >= 0) { panels[fT1].bottom_nbr_ind = c3; }

    if (do_adaptive_refine) {
        for (int ii = c0; ii <= c3; ++ii) { prospective_leaf_inds.push_back(ii); }
    }
    panels[jj].set_child_inds_start(c0);
}


// ---------------------------------------------------------------------
//  Split in v -> 2 children, [0]=bottom [1]=top.
//  Live edges: left, right.   Quiet edges: bottom, top.
// ---------------------------------------------------------------------
void AMRStructure::refine_one_v(int jj,
                                std::vector<double>& new_xs,
                                std::vector<double>& new_ps,
                                int& new_vert_ind,
                                std::vector<int>& prospective_leaf_inds,
                                bool do_adaptive_refine)
{
    int pi[9];
    for (int ii = 0; ii < 9; ++ii) { pi[ii] = panels[jj].point_inds[ii]; }
    const double x0 = xs[pi[0]], xm = xs[pi[3]], x1 = xs[pi[6]];
    const double v0 = ps[pi[0]], vm = ps[pi[1]], v1 = ps[pi[2]];
    const double q1 = 0.5 * (v0 + vm), q3 = 0.5 * (vm + v1);

    const int nbrL = panels[jj].left_nbr_ind;
    const int nbrR = panels[jj].right_nbr_ind;
    const int nbrB = panels[jj].bottom_nbr_ind;
    const int nbrT = panels[jj].top_nbr_ind;
    const bool p_left_bdry  = panels[jj].is_left_bdry;
    const bool p_right_bdry = panels[jj].is_right_bdry;
    const int child_level   = panels[jj].level + 1;
    const int panel_ind     = panels[jj].panel_ind;

    const int c0 = (int) panels.size(), c1 = c0 + 1;

    auto make_point = [&](double px, double pv) -> int {
        new_xs.push_back(px); new_ps.push_back(pv);
        return new_vert_ind++;
    };
    auto edge_point = [&](int edge, double target,
                          double px, double pv, int& found) -> int {
        int ind = lookup_edge_point(jj, edge, target, found);
        return (ind >= 0) ? ind : make_point(px, pv);
    };
    auto child_nbr = [](int parent_slot, int found) -> int {
        return (parent_slot == -2) ? -2 : found;
    };

    int fL0 = -1, fL1 = -1, fR0 = -1, fR1 = -1;
    int pL0 = edge_point(EDGE_LEFT,  q1, x0, q1, fL0);
    int pL1 = edge_point(EDGE_LEFT,  q3, x0, q3, fL1);
    int pR0 = edge_point(EDGE_RIGHT, q1, x1, q1, fR0);
    int pR1 = edge_point(EDGE_RIGHT, q3, x1, q3, fR1);

    int pC0 = make_point(xm, q1);   // interior: the children's centres
    int pC1 = make_point(xm, q3);

    int c0L = child_nbr(nbrL, fL0), c1L = child_nbr(nbrL, fL1);
    int c0R = child_nbr(nbrR, fR0), c1R = child_nbr(nbrR, fR1);
    if (nbrL == jj) { c0L = c0; c1L = c1; }
    if (nbrR == jj) { c0R = c0; c1R = c1; }

    // quiet edges: the child's bottom edge IS the parent's bottom edge,
    // same extent in x, so the parent's pointer is still exactly right --
    // including when the neighbour is refined further, which the leaf
    // search resolves by descending.
    const int c0B = nbrB;
    const int c1T = nbrT;

    panels.push_back(Panel {c0, child_level, panel_ind, 0,
            pi[0], pL0, pi[1],
            pi[3], pC0, pi[4],
            pi[6], pR0, pi[7],
            c0L, c1, c0R, c0B,
            p_left_bdry, p_right_bdry});
    panels.push_back(Panel {c1, child_level, panel_ind, 1,
            pi[1], pL1, pi[2],
            pi[4], pC1, pi[5],
            pi[7], pR1, pi[8],
            c1L, c1T, c1R, c0,
            p_left_bdry, p_right_bdry});

    if (fL0 >= 0) { panels[fL0].right_nbr_ind = c0; }
    if (fL1 >= 0) { panels[fL1].right_nbr_ind = c1; }
    if (fR0 >= 0) { panels[fR0].left_nbr_ind  = c0; }
    if (fR1 >= 0) { panels[fR1].left_nbr_ind  = c1; }
    if (nbrB >= 0 && nbrB != jj) { panels[nbrB].top_nbr_ind    = c0; }
    if (nbrT >= 0 && nbrT != jj) { panels[nbrT].bottom_nbr_ind = c1; }

    if (do_adaptive_refine) {
        prospective_leaf_inds.push_back(c0);
        prospective_leaf_inds.push_back(c1);
    }
    bool refining_in_p = true;
    panels[jj].set_child_inds_start(c0, refining_in_p);
}


// ---------------------------------------------------------------------
//  Split in x -> 2 children, [0]=left [1]=right.
//  Live edges: bottom, top.   Quiet edges: left, right.
// ---------------------------------------------------------------------
void AMRStructure::refine_one_x(int jj,
                                std::vector<double>& new_xs,
                                std::vector<double>& new_ps,
                                int& new_vert_ind,
                                std::vector<int>& prospective_leaf_inds,
                                bool do_adaptive_refine)
{
    int pi[9];
    for (int ii = 0; ii < 9; ++ii) { pi[ii] = panels[jj].point_inds[ii]; }
    const double x0 = xs[pi[0]], xm = xs[pi[3]], x1 = xs[pi[6]];
    const double v0 = ps[pi[0]], vm = ps[pi[1]], v1 = ps[pi[2]];
    const double m1 = 0.5 * (x0 + xm), m3 = 0.5 * (xm + x1);

    const int nbrL = panels[jj].left_nbr_ind;
    const int nbrR = panels[jj].right_nbr_ind;
    const int nbrB = panels[jj].bottom_nbr_ind;
    const int nbrT = panels[jj].top_nbr_ind;
    const bool p_left_bdry  = panels[jj].is_left_bdry;
    const bool p_right_bdry = panels[jj].is_right_bdry;
    const int child_level   = panels[jj].level + 1;
    const int panel_ind     = panels[jj].panel_ind;

    const int c0 = (int) panels.size(), c1 = c0 + 1;

    auto make_point = [&](double px, double pv) -> int {
        new_xs.push_back(px); new_ps.push_back(pv);
        return new_vert_ind++;
    };
    auto edge_point = [&](int edge, double target,
                          double px, double pv, int& found) -> int {
        int ind = lookup_edge_point(jj, edge, target, found);
        return (ind >= 0) ? ind : make_point(px, pv);
    };
    auto child_nbr = [](int parent_slot, int found) -> int {
        return (parent_slot == -2) ? -2 : found;
    };

    int fB0 = -1, fB1 = -1, fT0 = -1, fT1 = -1;
    int pB0 = edge_point(EDGE_BOTTOM, m1, m1, v0, fB0);
    int pB1 = edge_point(EDGE_BOTTOM, m3, m3, v0, fB1);
    int pT0 = edge_point(EDGE_TOP,    m1, m1, v1, fT0);
    int pT1 = edge_point(EDGE_TOP,    m3, m3, v1, fT1);

    int pM0 = make_point(m1, vm);   // interior: the children's centres
    int pM1 = make_point(m3, vm);

    int c0B = child_nbr(nbrB, fB0), c1B = child_nbr(nbrB, fB1);
    int c0T = child_nbr(nbrT, fT0), c1T = child_nbr(nbrT, fT1);

    // quiet edges: inherit verbatim
    int c0L = nbrL;
    int c1R = nbrR;
    if (nbrL == jj) { c0L = c1; }   // periodic, one panel across x
    if (nbrR == jj) { c1R = c0; }

    panels.push_back(Panel {c0, child_level, panel_ind, 0,
            pi[0], pi[1], pi[2],
            pB0,   pM0,   pT0,
            pi[3], pi[4], pi[5],
            c0L, c0T, c1, c0B,
            p_left_bdry, false});
    panels.push_back(Panel {c1, child_level, panel_ind, 1,
            pi[3], pi[4], pi[5],
            pB1,   pM1,   pT1,
            pi[6], pi[7], pi[8],
            c0, c1T, c1R, c1B,
            false, p_right_bdry});

    if (fB0 >= 0) { panels[fB0].top_nbr_ind    = c0; }
    if (fB1 >= 0) { panels[fB1].top_nbr_ind    = c1; }
    if (fT0 >= 0) { panels[fT0].bottom_nbr_ind = c0; }
    if (fT1 >= 0) { panels[fT1].bottom_nbr_ind = c1; }
    if (nbrL >= 0 && nbrL != jj) { panels[nbrL].right_nbr_ind = c0; }
    if (nbrR >= 0 && nbrR != jj) { panels[nbrR].left_nbr_ind  = c1; }

    if (do_adaptive_refine) {
        prospective_leaf_inds.push_back(c0);
        prospective_leaf_inds.push_back(c1);
    }
    panels[jj].set_child_inds_start_x(c0);
}

void AMRStructure::refine_panels(std::function<double (double,double)> f, bool do_adaptive_refine) {
    std::vector <double> new_xs;
    std::vector <double> new_ps;
    std::vector <double> new_fs;
    std::vector <int> prospective_leaf_inds;
    int new_vert_ind = xs.size();
    int num_panels_before_this_iter = panels.size();

    // Points created during this sweep live in new_xs/new_ps until the end,
    // so the edge lookup needs to be able to see them.
    staged_xs = &new_xs;
    staged_ps = &new_ps;

    for (int jj = minimum_unrefined_index; jj < num_panels_before_this_iter; ++jj) {
        if (! panels[jj].needs_refinement) { continue; }

        bool rx = panels[jj].needs_refine_x;
        bool rv = panels[jj].needs_refine_v;
        if (!rx && !rv) { rx = true; rv = true; }   // flagged with no direction: split both

        if (rx && rv) {
            refine_one_xv(jj, new_xs, new_ps, new_vert_ind,
                          prospective_leaf_inds, do_adaptive_refine);
        } else if (rv) {
            refine_one_v (jj, new_xs, new_ps, new_vert_ind,
                          prospective_leaf_inds, do_adaptive_refine);
        } else {
            refine_one_x (jj, new_xs, new_ps, new_vert_ind,
                          prospective_leaf_inds, do_adaptive_refine);
        }
    }

    // set fs at the newly created points
    new_fs.reserve(new_xs.size());
    for (int ii = 0; ii < new_xs.size(); ++ii) {
        new_fs.push_back( f(new_xs.at(ii), new_ps.at(ii)) );
    }
    for (int ii = 0; ii < new_xs.size(); ++ii) {
        xs.push_back(new_xs[ii]); ps.push_back(new_ps[ii]); fs.push_back(new_fs[ii]);
    }
    staged_xs = nullptr;
    staged_ps = nullptr;
}

void AMRStructure::generate_mesh(std::function<double (double,double)> f, 
                                 bool do_adaptive_refine, bool is_initial_step) 
{
    bool verbose=false;

    auto start = high_resolution_clock::now();
    // if (p_height > 0){
        // create_prerefined_mesh_p_refinement();
    // } else {
    create_prerefined_mesh();
    // }
    auto stop = high_resolution_clock::now();
    add_time(tree_build_time,  duration_cast<duration<double>>(stop - start) );



    start = high_resolution_clock::now();
    if (is_initial_step) {
        for (int ii = 0; ii < xs.size(); ii++) {
            fs[ii] = (*f0)(xs[ii],ps[ii]);
        }
    } else {
        int nx_points = 2*npanels_x + 1;
        int np_points = 2*npanels_p + 1;

        #ifdef DEBUG
        cout << "interpolating to grid " << endl;
        #endif
        #ifdef DEBUG_L2
        cout << "xs size " << xs.size() << endl;
        cout << "ps size " << ps.size() << endl;
        #endif

        interpolate_to_initial_xps(fs,xs,ps, nx_points, np_points,verbose);
        #ifdef DEBUG
        cout << "done interpolating to grid" << endl;
        #endif
    }
    stop = high_resolution_clock::now();
    add_time(interp_time, duration_cast<duration<double>>(stop - start) );

    int num_panels_pre_refine = panels.size();

// for debugging
#ifdef DEBUG
    // if (iter_num >= 236) {
    //     for (int ii = 0; ii < xs.size(); ++ii) {
    //         if (ps[ii] > 0.013) {
    //             if (xs[ii] > 0.005 && xs[ii] < 0.015) {
    //                 cout << "(x,v,f)_" << ii << "=(" << xs[ii] << ", " << ps[ii] << ", " << fs[ii]<< ")"<<endl;
    //             }
    //         }
    //     }
    // }
        // cout << "test initial grid for refinement" << endl;
        // for (int ii = minimum_unrefined_index; ii < num_panels_pre_refine; ++ii) {
        //     test_panel(ii, false);
        // }
#endif /* DEBUG */

    if (do_adaptive_refine) {
        start = high_resolution_clock::now();

        for (int ii = minimum_unrefined_index; ii < panels.size(); ++ii) {
            test_panel(ii, verbose);
        }
        // stop = high_resolution_clock::now();
        // add_time(panel_test_time,  duration_cast<duration<double>>(stop - start) );


        while (need_further_refinement) {
            #ifdef DEBUG
            cout << "making additional refinement passes" << endl;
            int current_max_height = 0;
            for (int ii = minimum_unrefined_index; ii < panels.size(); ++ii) {
                current_max_height = std::max(current_max_height, panels[ii].level);
            }
            cout << "highest current level is " << current_max_height << ", max allowed is " << max_height << endl;
            #endif /* DEBUG */
            need_further_refinement = false;
            #ifdef DEBUG
            cout << "refining panels" << endl;
            #endif
            auto amr_start = high_resolution_clock::now();
            refine_panels(f, do_adaptive_refine);
            auto amr_stop = high_resolution_clock::now();
            add_time(amr_refine_time, duration_cast<duration<double>>(amr_stop-amr_start) );

            amr_start = high_resolution_clock::now();
            // cout << "test initial grid for refinement" << endl;
            for (int ii = minimum_unrefined_index; ii < panels.size(); ++ii) {
                if (panels[ii].is_leaf()) {
                    test_panel(ii, verbose);
                }
            }
            amr_stop = high_resolution_clock::now();
            add_time(amr_test_time,  duration_cast<duration<double>>(amr_stop - amr_start) );
        }
        stop = high_resolution_clock::now();
        add_time(tree_build_time,  duration_cast<duration<double>>(stop - start) );
  
    }
    #ifdef DEBUG

    if (iter_num >= 236) {
        cout << "trying to debug" << endl;
        for (int ii = 0; ii < xs.size(); ++ii) {
            if (ps[ii] >= 0.013) {
                if (xs[ii] >= 0.005 && xs[ii] <= 0.015) {
                    cout << "(x,v,f)_" << ii << "=(" << xs[ii] << ", " << ps[ii] << ", " << fs[ii] <<")"<<endl;
                }
            }
        }
    }
    #endif

// #ifdef DEBUG
// cout << "fs at initialization" << endl;
// std::copy(fs.begin(), fs.end(), std::ostream_iterator<double>(cout, " "));
// cout << endl;
// #endif /* DEBUG */

    set_leaves_weights();

    { // Optional diagnostic -- delete freely.  Prints the split-type mix
      // and the point count each remesh.  Point-count growth is the signal
      // for how often the edge lookup is falling back to creating a point.
        int n_xv = 0, n_x = 0, n_v = 0, n_leaf = 0;
        for (size_t ii = 0; ii < panels.size(); ++ii) {
            if      (panels[ii].is_refined_xp) { ++n_xv; }
            else if (panels[ii].is_refined_x)  { ++n_x;  }
            else if (panels[ii].is_refined_p)  { ++n_v;  }
            else                               { ++n_leaf; }
        }
        std::cout << "[AMR] panels=" << panels.size()
                  << "  xv-split=" << n_xv
                  << "  x-split=" << n_x
                  << "  v-split=" << n_v
                  << "  leaves=" << n_leaf
                  << "  points=" << xs.size() << std::endl;
    }

    Q0 = 0;
    for (int ii = 0; ii < q_ws.size(); ii++) {
        Q0 += q_ws[ii];
    }
}

// ---------------------------------------------------------------------
//  Directional refinement criterion.
//
//  amr_epsilons[0] = eps_x   : threshold on the variation along x
//  amr_epsilons[1] = eps_v   : threshold on the variation along v
//  amr_epsilons[2] = eps_rel : relative threshold, applied to both
//                              directions; disabled when <= 0
//
//  var_x > eps_x only          -> split in x  (2 children, left/right)
//  var_v > eps_v only          -> split in v  (2 children, bottom/top)
//  both                        -> split in x and v (4 children)
// ---------------------------------------------------------------------
void AMRStructure::test_panel(int panel_ind, bool verbose) {

    // The refinement indicator is the background-subtracted distribution
    //     df = f - f_M(p)
    // not f itself.  Applied to full f, the equilibrium slope df_M/dp
    // dominates every panel in the bulk and the criterion saturates: at
    // t = 0 it either flags everything or, once normalised, flags the
    // tails preferentially -- neither of which tracks the perturbation
    // that the simulation is actually resolving.
    //
    // panel_fb keeps the background itself, because it is the natural
    // scale for the relative criterion below.
    double panel_fs[9];   // full f, kept for the interpolation-trouble check
    double panel_fb[9];   // f_M(p), the control variate
    double panel_df[9];   // f - f_M, the indicator
    auto panel_it = panels.begin() + panel_ind;
    for (int ii = 0; ii < 9; ++ii) {
        int pind = panel_it->point_inds[ii];
        double fb = (f0 != nullptr) ? f0->background(ps[pind]) : 0.0;
        panel_fs[ii] = fs[pind];
        panel_fb[ii] = fb;
        panel_df[ii] = fs[pind] - fb;
    }

    // Points are stored column-major, index = 3*i + j, i the x-column and
    // j the v-row:
    //      2 ----- 5 ----- 8
    //      1 ----- 4 ----- 7
    //      0 ----- 3 ----- 6
    // so v-row    j is {j, 3+j, 6+j}     (varies in x)
    // and x-column i is {3i, 3i+1, 3i+2} (varies in v).
    double var_x = 0.0;   // largest range along x, taken over the three rows
    double var_v = 0.0;   // largest range along v, taken over the three columns

    for (int jj = 0; jj < 3; ++jj) {
        double hi = panel_df[jj], lo = panel_df[jj];
        for (int ii = 1; ii < 3; ++ii) {
            double dij = panel_df[3*ii + jj];
            if (dij > hi) { hi = dij; }
            if (dij < lo) { lo = dij; }
        }
        if (hi - lo > var_x) { var_x = hi - lo; }
    }
    for (int ii = 0; ii < 3; ++ii) {
        double hi = panel_df[3*ii], lo = panel_df[3*ii];
        for (int jj = 1; jj < 3; ++jj) {
            double dij = panel_df[3*ii + jj];
            if (dij > hi) { hi = dij; }
            if (dij < lo) { lo = dij; }
        }
        if (hi - lo > var_v) { var_v = hi - lo; }
    }

    // full-panel range, kept only for the interpolation-trouble diagnostic
    double max_f = panel_fs[0], min_f = panel_fs[0];
    for (int ii = 1; ii < 9; ++ii) {
        if (panel_fs[ii] > max_f) { max_f = panel_fs[ii]; }
        if (panel_fs[ii] < min_f) { min_f = panel_fs[ii]; }
    }
    if (max_f - min_f >= 800) {
        cout << "interpolation trouble at panel " << panel_ind << endl;
        cout << "max f " << max_f << ", min f" << min_f << ", difference= " << max_f - min_f << endl;
        for (int ii = 0; ii < 9; ++ii) {
            int pind = panel_it->point_inds[ii];
            cout << "point " << pind << ": (x,v,f)=(" << xs[pind] << ", " << ps[pind] << ", " << panel_fs[ii] << ")" << endl;
        }
        cout << endl;
    }

    //   amr_epsilons = [eps_x, eps_v, eps_rel]
    //
    // Absolute thresholds on the range of df along each direction, then the
    // relative test in exactly the form used previously:
    //
    //     range > rtol * f_ref,     f_ref = max(|max df|, |min df|)
    //
    // written as a product so there is no division to guard.  It is also
    // self-guarding: f_ref == 0 only when every df on the panel is zero, in
    // which case var_x = var_v = 0 and nothing fires.
    bool refine_x = false;
    bool refine_v = false;
    if (amr_epsilons.size() > 0) { refine_x = (var_x > amr_epsilons[0]); }
    if (amr_epsilons.size() > 1) { refine_v = (var_v > amr_epsilons[1]); }
    if (amr_epsilons.size() > 2) {
        double max_df = panel_df[0], min_df = panel_df[0];
        for (int ii = 1; ii < 9; ++ii) {
            if (panel_df[ii] > max_df) { max_df = panel_df[ii]; }
            if (panel_df[ii] < min_df) { min_df = panel_df[ii]; }
        }
        const double f_ref = std::max(std::fabs(max_df), std::fabs(min_df));
        refine_x = refine_x || (var_x > amr_epsilons[2] * f_ref);
        refine_v = refine_v || (var_v > amr_epsilons[2] * f_ref);
    }

    if (panel_it->level < max_height && (refine_x || refine_v)) {
        panel_it->needs_refinement = true;
        panel_it->needs_refine_x   = refine_x;
        panel_it->needs_refine_v   = refine_v;
        need_further_refinement = true;
        if (verbose) {
            cout << "panel " << panel_ind << " is level " << panel_it->level
                 << ", max height " << max_height << ", flagged for refinement in ";
            if (refine_x && refine_v) { cout << "x and v"; }
            else if (refine_v)        { cout << "v only"; }
            else                      { cout << "x only"; }
            cout << " (var_x=" << var_x << ", var_v=" << var_v << ")" << endl;
        }
    }
    else if (verbose)
    {
        cout << "panel " << panel_ind << " is level " << panel_it->level
             << ", max height " << max_height;
        if (refine_x || refine_v) {
            cout << ", meets the criterion but is at max height" << endl;
        } else {
            cout << ", and is not flagged for refinement"
                 << " (var_x=" << var_x << ", var_v=" << var_v << ")" << endl;
        }
    }
}


void AMRStructure::set_leaves_weights() {
    leaf_inds = std::vector<int> ();
    // for (int ii = 0; ii < particles.size(); ++ii) {
    //     // particles[ii].q_w = 0;
    // }
    q_ws = std::vector<double> (xs.size());
    // rho_ws = std::vector<double> (xs.size());
    recursively_set_leaves_weights(0);

    // for (int ii = 0; ii < particles.size(); ++ii) {
    //     particles[ii].q_w *= particles[ii].f;
    // }
    for (int ii = 0; ii < xs.size(); ++ii) {
        q_ws[ii] *= fs[ii];
        // rho_ws[ii] *= fs[ii];
    }

}

void AMRStructure::recursively_set_leaves_weights(int panel_ind) {
    auto panel_it = panels.begin() + panel_ind;
    if (panel_it->is_refined_p || panel_it->is_refined_x) {
        // 2 children: bottom/top for a v-split, left/right for an x-split
        int child_start = panel_it->child_inds_start;
        for (int ii = 0; ii < 2; ii++) {
            recursively_set_leaves_weights(child_start + ii);
        }

    } else if (panel_it->is_refined_xp) {
        int child_start = panel_it->child_inds_start;
        for (int ii = 0; ii < 4; ii++) {
            recursively_set_leaves_weights(child_start + ii);
        }
    }
    else {
        leaf_inds.push_back(panel_ind);
        double dx = xs[panel_it->point_inds[3]] - xs[panel_it->point_inds[0]];
        // double v0 = particles[panel_it->vertex_inds[0]].v;
        // double v1 = particles[panel_it->vertex_inds[1]].v;
        double p0 = ps[panel_it->point_inds[0]];
        double p1 = ps[panel_it->point_inds[1]];
        double dp = p1 - p0;
        switch (quad) {
            case simpsons : {
                double qdxdp9 = q*dx * dp / 9.0;
                double weights[9] = {1.0,4.0,1.0, 4.0, 16.0,4.0,1.0,4.0,1.0};
                for (int ii = 0; ii < 9; ii++) {
                    q_ws[panel_it->point_inds[ii]] += qdxdp9 * weights[ii];
                }
                break;
            }
            default : {// trap 
                double qdxdp4 = q*dx * dp / 4.0;
                // cout << "area factor " << qdxdp4 << endl;
                double weights[9] = {1.0,2.0,1.0,2.0,4.0,2.0,1.0, 2.0, 1.0};
                for (int ii = 0; ii < 9; ii++) {
                    q_ws[panel_it->point_inds[ii]] += qdxdp4 * weights[ii];
                    // q_ws[panel_it->point_inds[ii]] += weights[ii];
                }
                break;
            }
        }
        // double qdp4 = q*dp/2.0;
        // double rho_weights[9] = {1.0,2.0,1.0,1.0,2.0,1.0,1.0, 2.0, 1.0};
        // for (int ii = 0; ii < 9; ii++) {
        //     rho_ws[panel_it->point_inds[ii]] += qdp4 * rho_weights[ii];
        // }     
    }
}

double mysqrt(double a) { return sqrt(a);};
double mysqr(double a) { return a*a;};

void AMRStructure::remesh() {

    if (sqrt_f) {
        // std::transform(fs.begin(), fs.end(), fs.begin(), mysqrt);
        for (int ii = 0; ii < fs.size(); ++ii) {
            fs[ii] = 2 * fs[ii];//sqrt(fs[ii]);
        }
    }
    // create copy of current panels and particle data

    // auto start = high_resolution_clock::now();
    old_panels = std::vector<Panel> (); old_panels.reserve(panels.size() );
    old_xs = std::vector<double> (xs); //old_xs.reserve(xs.size());
    old_ps = std::vector<double> (ps); //old_ps.reserve(xs.size());
    old_fs = std::vector<double> (fs); //old_fs.reserve(xs.size());


    for (const auto& panel : panels) {
        old_panels.push_back(Panel (panel));
    }
    
    // for (int ii = 0; ii < xs.size(); ++ii ) {
    //     old_xs.push_back(xs[ii]);
    //     old_ps.push_back(ps[ii]);
    //     old_fs.push_back(fs[ii]);
    // }

    // auto stop = high_resolution_clock::now();
    // auto duration = duration_cast<microseconds>(stop - start);

    // cout << "Old data copy time " << duration.count() << " microseconds." << endl << endl;

    bool is_initial_step = false;
    #ifdef DEBUG
    cout << "Generating mesh" << endl;
    #endif
    generate_mesh([&] (double x, double v) { return interpolate_from_mesh(x,v,false);} , do_adaptively_refine, is_initial_step);
    

    if (sqrt_f) {
        // std::transform(fs.begin(), fs.end(), fs.begin(), mysqr);
        for (int ii = 0; ii < fs.size(); ++ii) {
            fs[ii] = 0.5 * fs[ii];// * fs[ii];
        }
    }
    // init_e();
}