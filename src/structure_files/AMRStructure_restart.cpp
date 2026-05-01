/**
 * @file AMRStructure_restart.cpp
 * @brief Restart support: dump the full panel tree (not just leaves) and
 *        load (xs, ps, fs, panels) back from disk so a simulation can be
 *        continued from a previous diagnostic dump.
 *
 * Add this file to src/structure_files/ and to your build (Makefile / CMake).
 *
 * On-disk layout (per species, in sim_dir/simulation_output/<species_name>/):
 *   xs/xs_<iter>             : already produced, raw doubles, one per particle
 *   ps/ps_<iter>             : already produced, raw doubles, one per particle
 *   fs/fs_<iter>             : already produced, raw doubles, one per particle
 *   panels/leaf_point_inds_<iter>     : already produced (leaves only)
 *   panels/full_tree_<iter>           : NEW, full panel tree (this file)
 *
 * The full_tree file is a sequence of fixed-size records, one per Panel.
 * Each record is packed as POD ints/bools in the order matching the
 * struct layout (see PanelDumpRecord below). We do NOT rely on the C++
 * struct layout directly because Panel mixes ints and bools; we marshal
 * explicitly to be portable across compilers.
 */

#include "AMRStructure.hpp"
#include <fstream>
#include <iostream>
#include <cstdint>
#include <stdexcept>

namespace {

// Plain layout used for panel-tree serialization. Keep this stable.
// All ints are int32; bools are stored as int8.
struct PanelDumpRecord {
    int32_t panel_ind;
    int32_t level;
    int32_t parent_ind;
    int32_t which_child;
    int32_t point_inds[9];
    int32_t left_nbr_ind;
    int32_t top_nbr_ind;
    int32_t right_nbr_ind;
    int32_t bottom_nbr_ind;
    int32_t child_inds_start;
    int8_t  is_left_bdry;
    int8_t  is_right_bdry;
    int8_t  is_refined_xp;
    int8_t  is_refined_p;
    int8_t  needs_refinement;
    int8_t  _pad[3]; // keep 4-byte alignment
};

} // anonymous namespace


int AMRStructure::write_panel_tree_to_file(int iter_num) {
    std::ofstream tree_file;
    std::string fname = sim_dir + "simulation_output/" + species_name
                      + "/panels/full_tree_" + std::to_string(iter_num);
    tree_file.open(fname, std::ios::out | std::ios::binary);
    if (!tree_file) {
        std::cout << "Unable to open " << fname << " for writing panel tree."
                  << std::endl;
        return 1;
    }

    // First: write the number of panels as int64 so the loader knows.
    int64_t n_panels = static_cast<int64_t>(panels.size());
    tree_file.write(reinterpret_cast<const char*>(&n_panels), sizeof(int64_t));

    PanelDumpRecord rec;
    for (const Panel& p : panels) {
        rec.panel_ind   = p.panel_ind;
        rec.level       = p.level;
        rec.parent_ind  = p.parent_ind;
        rec.which_child = p.which_child;
        for (int i = 0; i < 9; ++i) rec.point_inds[i] = p.point_inds[i];
        rec.left_nbr_ind   = p.left_nbr_ind;
        rec.top_nbr_ind    = p.top_nbr_ind;
        rec.right_nbr_ind  = p.right_nbr_ind;
        rec.bottom_nbr_ind = p.bottom_nbr_ind;
        rec.child_inds_start = p.child_inds_start;
        rec.is_left_bdry      = p.is_left_bdry      ? 1 : 0;
        rec.is_right_bdry     = p.is_right_bdry     ? 1 : 0;
        rec.is_refined_xp     = p.is_refined_xp     ? 1 : 0;
        rec.is_refined_p      = p.is_refined_p      ? 1 : 0;
        rec.needs_refinement  = p.needs_refinement  ? 1 : 0;
        rec._pad[0] = rec._pad[1] = rec._pad[2] = 0;
        tree_file.write(reinterpret_cast<const char*>(&rec), sizeof(rec));
    }

    if (!tree_file.good()) {
        std::cout << "Error writing panel tree for step " << iter_num
                  << std::endl;
        return 1;
    }
    tree_file.close();
    return 0;
}


int AMRStructure::load_particles_from_file(int iter_num) {
    // Read xs, ps, fs binary files for this iter_num. Each is a stream of
    // doubles; all three must have the same length.
    auto read_doubles = [&](const std::string& sub, std::vector<double>& out) -> int {
        std::string fname = sim_dir + "simulation_output/" + species_name
                          + "/" + sub + "/" + sub + "_" + std::to_string(iter_num);
        std::ifstream fin(fname, std::ios::in | std::ios::binary | std::ios::ate);
        if (!fin) {
            std::cout << "Unable to open " << fname << " for restart load."
                      << std::endl;
            return 1;
        }
        std::streamsize bytes = fin.tellg();
        if (bytes < 0 || (bytes % sizeof(double)) != 0) {
            std::cout << "Bad file size for " << fname << std::endl;
            return 1;
        }
        size_t n = static_cast<size_t>(bytes) / sizeof(double);
        out.resize(n);
        fin.seekg(0, std::ios::beg);
        fin.read(reinterpret_cast<char*>(out.data()),
                 static_cast<std::streamsize>(n * sizeof(double)));
        if (!fin) {
            std::cout << "Error reading " << fname << std::endl;
            return 1;
        }
        return 0;
    };

    if (read_doubles("xs", xs)) return 1;
    if (read_doubles("ps", ps)) return 1;
    if (read_doubles("fs", fs)) return 1;

    if (xs.size() != ps.size() || ps.size() != fs.size()) {
        std::cout << "Mismatched xs/ps/fs sizes on restart load: "
                  << xs.size() << " / " << ps.size() << " / " << fs.size()
                  << std::endl;
        return 1;
    }
    std::cout << "Restart: loaded " << xs.size() << " particles for species "
              << species_name << " at iter " << iter_num << std::endl;
    return 0;
}


int AMRStructure::load_panel_tree_from_file(int iter_num) {
    std::string fname = sim_dir + "simulation_output/" + species_name
                      + "/panels/full_tree_" + std::to_string(iter_num);
    std::ifstream fin(fname, std::ios::in | std::ios::binary);
    if (!fin) {
        std::cout << "Unable to open " << fname << " for restart load. "
                  << "Did you re-run with the patched write_to_file?"
                  << std::endl;
        return 1;
    }

    int64_t n_panels = 0;
    fin.read(reinterpret_cast<char*>(&n_panels), sizeof(int64_t));
    if (!fin || n_panels <= 0) {
        std::cout << "Bad header in panel tree file " << fname << std::endl;
        return 1;
    }

    panels.clear();
    panels.reserve(static_cast<size_t>(n_panels));

    PanelDumpRecord rec;
    for (int64_t i = 0; i < n_panels; ++i) {
        fin.read(reinterpret_cast<char*>(&rec), sizeof(rec));
        if (!fin) {
            std::cout << "Error reading record " << i << " from "
                      << fname << std::endl;
            return 1;
        }
        Panel p; // default-constructed
        p.panel_ind   = rec.panel_ind;
        p.level       = rec.level;
        p.parent_ind  = rec.parent_ind;
        p.which_child = rec.which_child;
        for (int k = 0; k < 9; ++k) p.point_inds[k] = rec.point_inds[k];
        p.left_nbr_ind     = rec.left_nbr_ind;
        p.top_nbr_ind      = rec.top_nbr_ind;
        p.right_nbr_ind    = rec.right_nbr_ind;
        p.bottom_nbr_ind   = rec.bottom_nbr_ind;
        p.child_inds_start = rec.child_inds_start;
        p.is_left_bdry     = (rec.is_left_bdry     != 0);
        p.is_right_bdry    = (rec.is_right_bdry    != 0);
        p.is_refined_xp    = (rec.is_refined_xp    != 0);
        p.is_refined_p     = (rec.is_refined_p     != 0);
        p.needs_refinement = (rec.needs_refinement != 0);
        panels.push_back(p);
    }

    std::cout << "Restart: loaded " << panels.size() << " panels for species "
              << species_name << " at iter " << iter_num << std::endl;
    return 0;
}


int AMRStructure::load_restart(int iter_num) {
    // 1. Load particle arrays
    if (load_particles_from_file(iter_num)) return 1;
    // 2. Load full panel tree
    if (load_panel_tree_from_file(iter_num)) return 1;

    // 3. Rebuild leaf_inds and q_ws from the loaded tree + xs/ps/fs.
    //    set_leaves_weights walks the tree starting at root (panel 0),
    //    populates leaf_inds, computes the trapezoid/Simpson weights from
    //    the panel geometry, and multiplies by fs at each point.
    set_leaves_weights();

    // 4. Mark internal flags consistent with a freshly remeshed state.
    is_initial_mesh_set    = true;
    need_further_refinement = false;
    minimum_unrefined_index = 0;

    // 5. The 'height' member is informational; recompute as max over panels.
    int max_lvl = 0;
    for (const auto& p : panels) {
        if (p.level > max_lvl) max_lvl = p.level;
    }
    height = max_lvl;

    std::cout << "Species " << species_name
              << " restart complete: " << xs.size() << " particles, "
              << panels.size() << " panels, " << leaf_inds.size()
              << " leaves, max level " << height << std::endl;
    return 0;
}