#include "AMRSimulation.hpp"

AMRSimulation::AMRSimulation() {}

AMRSimulation::AMRSimulation(std::string sim_dir, std::string deck_address) 
    : need_scatter (false), outFile("debug.txt")
{
    this->sim_dir = sim_dir;
    this->deck_address = deck_address;

    // create ptree
    pt::ptree deck;
    load_deck(deck_address, deck);
    cout << "deck address: " << deck_address << endl;

    get_box_t_params(deck);

    // create e solver
    calculate_e = make_field_return_ptr(deck);
    make_external_field(deck);

    // load species
    try {
        pt::ptree &species_list_deck = deck.get_child("species_list");
        for (pt::ptree::value_type &sp : species_list_deck) {
            pt::ptree &sp_deck = sp.second;
            distribution* f0 = make_f0_return_ptr(sp_deck);
            ic_list.push_back(f0 );
            species_list.push_back(make_species_return_ptr(sp_deck, f0));
        }
    } catch(std::exception& e) {
        cout << "Invalid deck format.  Must have at least one species in a species list" << endl;
        return;
    }
    N_sp = species_list.size();
    get_qms();

#if TESTFLAG
    outFile << "Debugging is enabled!" << std::endl;
    outFile << "All species after loading deck: " << endl;
    for (int i = 0; i < N_sp; i++) {
        outFile << species_list[i]->species_name << "  , xs:  size = " << species_list[i]->xs.size() << endl;
        for (int j = 0; j < species_list[i]->xs.size(); j++) {
            outFile << species_list[i]->xs[j] << "  ";
        }
        outFile << endl;
        outFile << species_list[i]->species_name << "  , ps:  size = " << species_list[i]->ps.size() << endl;
        for (int j = 0; j < species_list[i]->ps.size(); j++) {
            outFile << species_list[i]->ps[j]  << "  ";
        }
        outFile << endl;
        outFile << species_list[i]->species_name << "  , fs: size = " << species_list[i]->fs.size() << endl;
        for (int j = 0; j < species_list[i]->fs.size(); j++) {
            outFile << species_list[i]->fs[j]  << "  ";
        }
        outFile << endl;
        outFile << species_list[i]->species_name << "  , q_ws: size = " << species_list[i]->q_ws.size() << endl;
        for (int j = 0; j < species_list[i]->q_ws.size(); j++) {
            outFile << species_list[i]->q_ws[j]  << "  ";
        }
        outFile << endl;
    }
#endif

   if (restart) {
        // We still build initial meshes from f0 above (in make_species_return_ptr),
        // but immediately overwrite that state from disk.
        iter_num = restart_iter;
        t = restart_iter * dt;
        if (load_restart_state() != 0) {
            cout << "ERROR: failed to load restart state at iter "
                 << restart_iter << ". Aborting." << endl;
            // leave the object in a defined-but-empty state; run() will
            // just immediately exit if iter_num >= num_steps anyway.
            return;
        }
        // Recompute E on the loaded state.
        evaluate_field_uniform_grid(t);
        // Do NOT call write_to_file() here -- that file already exists on
        // disk from the previous run.
        cout << "Restart complete. Resuming at iter " << iter_num
             << ", t = " << t
             << ". Will run to iter " << num_steps << "." << endl;
    } else {
        iter_num = 0;
        t = 0;
        evaluate_field_uniform_grid(t);
        write_to_file();
    }
 
    print_sim_setup();

}

//destructor
AMRSimulation::~AMRSimulation() {
    for (int ii = 0; ii < ic_list.size(); ++ii) {
        delete ic_list[ii];
        delete species_list[ii];
    }
    delete calculate_e;
}


int AMRSimulation::load_restart_state() {
    cout << "Loading restart state at iter " << restart_iter
         << " from " << sim_dir << " ..." << endl;
    for (auto& species : species_list) {
        if (species->load_restart(restart_iter) != 0) {
            cout << "Restart load failed for species "
                 << species->species_name << endl;
            return 1;
        }
    }
    // After loading, the per-species xs/ps/fs/q_ws are valid. We need to
    // mark that the simulation-level (gathered) buffers are stale so the
    // next step()/evaluate_field call regathers them.
    need_gather  = true;
    need_scatter = false;
    return 0;
}

