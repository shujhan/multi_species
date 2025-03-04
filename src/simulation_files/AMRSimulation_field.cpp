#include "AMRSimulation.hpp"
#include <unordered_map>
#include <unordered_set>

int AMRSimulation::evaluate_field_uniform_grid(double t) {
    // get reduced xs, reduced_ws from each species
    // combine reduced x, reduced_ws
    std::vector<double> reduced_xs, reduced_ws;
    for (auto &species : species_list) {
        species->get_reduced_xs_ws();
    }

// combine all species xs and ws without duplicates
    unordered_set<double> seen; // check unique elements
    unordered_map<double, size_t> indexMap; // Map values to indices
    for (size_t this_sp = 0; this_sp < species_list.size(); this_sp++) {
        for (double val : species_list[this_sp]->reduced_xs) {
            if (seen.insert(val).second) {
                indexMap[val] = reduced_xs.size();
                reduced_xs.push_back(val);
            }
        }
    }
    // store index, index_multi represents the index in the general reduced_xs
    for (size_t this_sp = 0; this_sp < species_list.size(); this_sp++) {
        for (double val : species_list[this_sp]->reduced_xs) {
            species_list[this_sp]->index_multi.push_back(indexMap[val]);
        }
    }
    // get reduced_ws 
    reduced_ws.resize(reduced_xs.size());
    for (size_t this_sp = 0; this_sp < species_list.size(); this_sp++) {
        for (size_t i = 0; i < species_list[this_sp]->reduced_ws.size(); i++) {
            reduced_ws[species_list[this_sp]->index_multi[i]] += species_list[this_sp]->reduced_ws[i];
        }
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


// For amr, distribute reduced es here to each species
// then in each species it will distribute again
    for (auto &species : species_list) {
        std::vector<double> this_es;
        // species_list[this_sp]->index_multi
        for (size_t index : species->index_multi) {
            this_es.push_back(reduced_es[index]);
        }
        species->get_reduced_es(this_es.data());
    }

    need_gather = true;

    for (size_t this_sp = 0; this_sp < species_list.size(); this_sp++) {
        species_list[this_sp]->index_multi.clear();
    }

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
