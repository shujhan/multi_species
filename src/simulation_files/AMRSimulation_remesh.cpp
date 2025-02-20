#include "AMRSimulation.hpp"

int AMRSimulation::remesh() {
    if (need_scatter) {
        //scatter to species
        bool send_e = false;
        scatter(send_e);
    }

    // remesh
    for (auto &species : species_list) {
        species->remesh();
    }


// after remesh print out 
#if TESTFLAG
    outFile << "post remesh!" << std::endl;
    for (int i = 0; i < N_sp; i++) {
        outFile << species_list[i]->species_name << "  , xs:  size = " << species_list[i]->xs.size() << endl;
        for (int j = 0; j < species_list[i]->xs.size(); j++) {
            outFile << setprecision(15) << species_list[i]->xs[j] << "  ";
        }
        outFile << endl;
        outFile << species_list[i]->species_name << "  , ps:  size = " << species_list[i]->ps.size() << endl;
        for (int j = 0; j < species_list[i]->ps.size(); j++) {
            outFile << setprecision(15) << species_list[i]->ps[j]  << "  ";
        }
        outFile << endl;
        outFile << species_list[i]->species_name << "  , fs: size = " << species_list[i]->fs.size() << endl;
        for (int j = 0; j < species_list[i]->fs.size(); j++) {
            outFile << setprecision(15) << species_list[i]->fs[j]  << "  ";
        }
        outFile << endl;
        outFile << species_list[i]->species_name << "  , q_ws: size = " << species_list[i]->q_ws.size() << endl;
        for (int j = 0; j < species_list[i]->q_ws.size(); j++) {
            outFile << setprecision(15) << species_list[i]->q_ws[j]  << "  ";
        }
        outFile << endl;
    }
#endif


    evaluate_field_uniform_grid(t);

    need_gather = true;
    return 0;
}
