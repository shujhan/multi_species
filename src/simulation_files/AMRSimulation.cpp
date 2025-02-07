#include "AMRSimulation.hpp"

AMRSimulation::AMRSimulation() {}

AMRSimulation::AMRSimulation(std::string sim_dir, std::string deck_address) 
    : need_scatter (false)
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
    std::ofstream outFile("debug.txt");
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
        outFile << species_list[i]->species_name << "  , q_ws: size = " << species_list[i]->q_ws.size() << endl;
        for (int j = 0; j < species_list[i]->q_ws.size(); j++) {
            outFile << species_list[i]->q_ws[j]  << "  ";
        }
        outFile << endl;
    }
#endif

    // initialize e
    iter_num = 0;
    t = 0;
    evaluate_field_uniform_grid(t);
#if TESTFLAG
    outFile << "Reduced xs and es for each species before RK4: " << std::endl;
    for (int i = 0; i < N_sp; i++) {
        outFile << species_list[i]->species_name << "  , reduced_xs:  size = " << species_list[i]->reduced_xs.size() << endl;
        for (int j = 0; j < species_list[i]->reduced_xs.size(); j++) {
            outFile << species_list[i]->reduced_xs[j]  << "  ";
        }
        outFile << endl;
        outFile << species_list[i]->species_name << "  , reduced_ws:  size = " << species_list[i]->reduced_ws.size() << endl;
        for (int j = 0; j < species_list[i]->reduced_ws.size(); j++) {
            outFile << species_list[i]->reduced_ws[j]  << "  ";
        }
        outFile << endl;
         outFile << species_list[i]->species_name << "  , es:  size = " << species_list[i]->es.size() << endl;
        for (int j = 0; j < species_list[i]->es.size(); j++) {
            outFile << species_list[i]->es[j]  << "  ";
        }
        outFile << endl;
    }
    // TODO add the general simulation::reduced_es and reduced_xs, reduced_ws
#endif

    //write to file
    write_to_file();

    // print AMR description
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

