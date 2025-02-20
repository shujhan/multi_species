#include "AMRSimulation.hpp"


int AMRSimulation::evaluate_field_uniform_grid(double t) {
    // get reduced xs, reduced_ws from each species
    // combine reduced x, reduced_ws
    std::vector<double> reduced_xs, reduced_ws;
    for (auto &species : species_list) {
        species->get_reduced_xs_ws();
        reduced_xs.insert( reduced_xs.end(), species->reduced_xs.begin(), species->reduced_xs.end());
        reduced_ws.insert( reduced_ws.end(), species->reduced_ws.begin(), species->reduced_ws.end());
    }

    // duplicate reduced x
    std::vector<double> reduced_xs_cpy (reduced_xs);
    // evaluate reduced e
    std::vector<double> reduced_es(reduced_ws.size());
    (*calculate_e)(reduced_es.data(), reduced_xs.data(), reduced_xs.size(),
                 reduced_xs_cpy.data(), reduced_ws.data(), reduced_xs.size() );
    // if (use_external_field) {
    //     (*calculate_e_external)(reduced_es.data(), reduced_xs.data(), reduced_xs.size(), t);
    // }

   // distribute to each species
    size_t start_ind = 0;
    for (auto &species : species_list) {
        species->get_reduced_es(reduced_es.data() + start_ind);
        start_ind += species->reduced_xs.size();
    }
    need_gather = true;

#if TESTFLAG
    outFile << "Reduced xs and es for each species: " << std::endl;
    for (int i = 0; i < N_sp; i++) {
        outFile << species_list[i]->species_name << "  , reduced_xs:  size = " << species_list[i]->reduced_xs.size() << endl;
        for (int j = 0; j < species_list[i]->reduced_xs.size(); j++) {
            outFile << setprecision(15) << species_list[i]->reduced_xs[j]  << "  ";
        }
        outFile << endl;
        outFile << species_list[i]->species_name << "  , reduced_ws:  size = " << species_list[i]->reduced_ws.size() << endl;
        for (int j = 0; j < species_list[i]->reduced_ws.size(); j++) {
            outFile << setprecision(15) << species_list[i]->reduced_ws[j]  << "  ";
        }
        outFile << endl;
        outFile << species_list[i]->species_name << "  , reduced_es:  size = " << species_list[i]->sort_es.size() << endl;
        for (int j = 0; j < species_list[i]->sort_es.size(); j++) {
            outFile << setprecision(15) << species_list[i]->sort_es[j]  << "  ";
        }
        outFile << endl;

        outFile << species_list[i]->species_name << "  , general reduced_xs:  size = " << reduced_xs.size() << endl;
        for (int j = 0; j < reduced_xs.size(); j++) {
            outFile << setprecision(15) << reduced_xs[j]  << "  ";
        }
        outFile << endl;

        outFile << species_list[i]->species_name << "  , general reduced_ws:  size = " << reduced_ws.size() << endl;
        for (int j = 0; j < reduced_ws.size(); j++) {
            outFile << setprecision(15) << reduced_ws[j]  << "  ";
        }
        outFile << endl;

        outFile << species_list[i]->species_name << "  , e compare to farsight :  size = " << reduced_es.size() << endl;
        for (int j = 0; j < reduced_es.size(); j++) {
            outFile << setprecision(15) << reduced_es[j]  << "  ";
        }
        outFile << endl;
    }
#endif


    return 0;
}

int AMRSimulation::evaluate_field(std::vector<double>& es_local, std::vector<double>& xs_local, std::vector<double>& q_ws_local, double t) {

    std::vector<double> xtemp_cpy(xs_local), xtemp_cpy2(xs_local);
    (*calculate_e)(es_local.data(), xtemp_cpy.data(), es_local.size(),
                    xtemp_cpy2.data(), q_ws_local.data(), xtemp_cpy.size());
    if (use_external_field) {
        (*calculate_e_external)(es_local.data(), xs_local.data(), es.size(), t);
    }
    return 0;
}
