#include "FieldStructure.hpp"
#include <cmath>
#include <cassert>
#include <cstring>
#include <sys/times.h>
// #include <openacc.h>
#include <iostream>
#include <cfloat> // dbl_min
#include <cstddef>
using namespace std;
#if OPENACC_ENABLED
#include <accelmath.h>
#endif

ExternalElectricField::~ExternalElectricField() = default;
// -----------------
ExternalTanh::ExternalTanh()
        : Am(0.2), k(0.26), omega(0.37) {
            eps = drive_fun(t0);
        }
ExternalTanh::ExternalTanh(double Am)
        : Am(Am), k(0.26), omega(0.37) {
            eps = drive_fun(t0);
        }    
ExternalTanh::ExternalTanh(double Am, double k, double omega)
        : Am(Am), k(k), omega(omega) {
            eps = drive_fun(t0);
        }   
double ExternalTanh::drive_fun(double t) {
    return  0.5 * (tanh((t - tL)/twL) - tanh((t-tR)/twR));
}
void ExternalTanh::operator() (double* es, double* targets, int nt, double t) {
    double a_t_k = (drive_fun(t) - eps) / (1-eps) * Am * k;
    double omt = omega * t;
    for (int ii = 0; ii < nt; ++ii) {
        es[ii] += a_t_k * sin(k*targets[ii] - omt);
    }
}
void ExternalTanh::print_field_obj() {
    cout << "External field : Tanh" << endl;
    cout << "Am " << Am << ", k " << k << ", omega " << omega << endl;
    cout << "t0 " << t0 << ", tL " << tL << ", tR " << tR << endl;
    cout << "twL " << twL << ", twR " << twR << endl;
    cout << "epsilon driver " << eps << endl;
}
ExternalTanh::~ExternalTanh() = default;
// -----------------
ExternalSine::ExternalSine()
        : Am(0.4), k(0.26), omega(0.37) {}
ExternalSine::ExternalSine(double Am)
        : Am(Am), k(0.26), omega(0.37) {}    
ExternalSine::ExternalSine(double Am, double k, double omega)
        : Am(Am), k(k), omega(omega) {}   
void ExternalSine::operator() (double* es, double* targets, int nt, double t) {
    double a0 = 0;
    if (t < t1) {
        a0 = Am * sin (t * M_PI / 100);
    } else if (t < t2) {
        a0 = Am;
    } else if (t < t3) {
        a0 = Am * cos((t-150)*M_PI / 100);
    }
    double omt = omega * t;
    for (int ii = 0; ii < nt; ++ii) {
        es[ii] += a0 * sin(k*targets[ii] - omt);
    }
}
void ExternalSine::print_field_obj() {
    cout << "External field : sine" << endl;
    cout << "Am " << Am << ", k " << k << ", omega " << omega << endl;
}
ExternalSine::~ExternalSine() = default;
//--------------------------
ExternalLogistic::ExternalLogistic()
        : Am(0.4), k(0.26), omega(0.37) {}
ExternalLogistic::ExternalLogistic(double Am)
        : Am(Am), k(0.26), omega(0.37) {}    
ExternalLogistic::ExternalLogistic(double Am, double k, double omega)
        : Am(Am), k(k), omega(omega) {}   
void ExternalLogistic::operator() (double* es, double* targets, int nt, double t) {
    double a0;
    if (t < t1) {
        a0 = Am / (1.0 + exp(-40.0*(t-10.0)));
    } else {
        a0 = Am * (1.0 - 1.0 / (1.0 + exp(-40.0 * (t - 110.0))));
    }
    double omt = omega * t;
    for (int ii = 0; ii < nt; ++ii) {
        es[ii] += a0 * sin(k*targets[ii] - omt);
    }
}
void ExternalLogistic::print_field_obj() {
    cout << "External field : logistic" << endl;
    cout << "Am " << Am << ", k " << k << ", omega " << omega << endl;
}
ExternalLogistic::~ExternalLogistic() = default;
//================== end external ======

ElectricField::~ElectricField() = default;


E_MQ_DirectSum::E_MQ_DirectSum() {}
E_MQ_DirectSum::E_MQ_DirectSum(double L, double epsilon) : L(L), epsilon(epsilon) {}
E_MQ_DirectSum::~E_MQ_DirectSum() = default;

void E_MQ_DirectSum::operator() (double* es, double* targets, int nt, 
                        double* sources, double* q_ws, int ns)
{    
    double epsLsq = epsilon * epsilon / L / L;
    double norm_epsL = sqrt(1 + 4 * epsLsq );

#ifdef OPENACC_ENABLED
std::cout << "Running with OpenACC" << std::endl;
#pragma acc parallel loop independent
#else
std::cout << "Running without OpenACC" << std::endl;
#pragma omp parallel for
#endif
    for (int ii = 0; ii < nt; ++ii) {
        double xi = targets[ii];
        double ei = 0.0;
#ifdef OPENACC_ENABLED
#pragma acc loop independent reduction(+:ei)
#endif
        for (int jj = 0; jj < ns; ++jj) {
            double z = (xi - sources[jj]) / L;
            z = z - round(z);
            // while (z < -0.5) { z += 1.0; }
            // while (z >= 0.5) { z -= 1.0; }
            ei += q_ws[jj] * (0.5 * z * norm_epsL / sqrt(z * z + epsLsq) - z);
        }
        es[ii] = ei;
    }
}
void E_MQ_DirectSum::print_field_obj() {
    cout << "-------------" << endl;
    cout << "Field object: " << endl;
    cout << "MQ kernel, direct sum" << endl;
    cout << "epsilon = " << epsilon << endl;
    cout << "Periodic boundary conditions, domain size = " << L << endl;
    cout << "-------------" << endl;
}

E_MQ_DirectSum_openbcs::E_MQ_DirectSum_openbcs() {}
E_MQ_DirectSum_openbcs::E_MQ_DirectSum_openbcs(double epsilon) : epsilon(epsilon) {}
E_MQ_DirectSum_openbcs::~E_MQ_DirectSum_openbcs() = default;

void E_MQ_DirectSum_openbcs::operator() (double* es, double* targets, int nt, 
                        double* sources, double* q_ws, int ns)
{    
    double eps_sq = epsilon * epsilon;

#ifdef OPENACC_ENABLED
#pragma acc parallel loop independent
#else
#pragma omp parallel for
#endif
    for (int ii = 0; ii < nt; ++ii) {
        double xi = targets[ii];
        double ei = 0.0;
#ifdef OPENACC_ENABLED
#pragma acc loop independent reduction(+:ei)
#endif
        for (int jj = 0; jj < ns; ++jj) {
            double z = xi - sources[jj];
            ei += q_ws[jj] * 0.5 * z / sqrt(z * z + eps_sq);
        }
        es[ii] = ei - xi;
    }
}
void E_MQ_DirectSum_openbcs::print_field_obj() {
    cout << "-------------" << endl;
    cout << "Field object: " << endl;
    cout << "MQ kernel, direct sum" << endl;
    cout << "epsilon = " << epsilon << endl;
    cout << "Open boundary conditions "<< endl;
    cout << "-------------" << endl;
}

E_MQ_Treecode::E_MQ_Treecode() {}
E_MQ_Treecode::E_MQ_Treecode(double L, double epsilon, double beta) : 
L(L), epsilon(epsilon), mac(mac), beta(beta) {}

// Definitions for E_MQ_Treecode
E_MQ_Treecode::E_MQ_Treecode(double L, double epsilon,
                             double mac, int degree, int max_source, int max_target,
                             int verbosity)
    : L(L), epsilon(epsilon), mac(mac), degree(degree),
      max_source(max_source), max_target(max_target), verbosity(verbosity),
      lambda(nullptr), particles_x(nullptr), keps_tc_reord(nullptr), keps_tc_noreord(nullptr), 
      tree_members{nullptr, nullptr}, leaf_members{nullptr, nullptr}, interaction_list_far(nullptr), 
      interaction_list_far_size(nullptr), max_far_size(0), interaction_list_near(nullptr), interaction_list_near_size(nullptr), max_near_size(0), 
      cluster_list_t1(nullptr), cluster_list_moments(nullptr) {
        #ifdef OPENACC_ENABLED
        #pragma acc enter data copyin(this)
        #endif
      }


    E_MQ_Treecode::~E_MQ_Treecode() 
    {
    #ifdef OPENACC_ENABLED
      #pragma acc exit data delete(this)
    #endif
    };


void E_MQ_Treecode::operator() (double* es, double* targets, int nt, 
                double* sources, double* q_ws, int ns)
{
    //  dont need these i think: check later
    // numpars_s = static_cast<size_t> (ns);
    numpars_s = ns;
    particles_x = new double[numpars_s];
    lambda = new double[numpars_s];


    for (size_t i = 0; i < numpars_s; i++) {
        // particles_x[i] = sources[i];
        particles_x[i] = fmod(sources[i],L);
        lambda[i] = q_ws[i];
    }

    P = degree;
    PP = P+1;
    Pflat = PP;
    sq_theta = mac * mac;
    epsLsq = epsilon * epsilon / (L * L);
    norm_epsL = sqrt(1 + 4 * epsLsq);
    N0 = max_source;

#if OPENACC_ENABLED
std::cout << "Running with OpenACC" << std::endl;
#pragma acc enter data copyin(lambda[0:numpars_s])
#pragma acc enter data copyin(particles_x[0:numpars_s])
#else
std::cout << "Running without OpenACC" << std::endl;
#endif


    //============  BLTC ============
    keps_tc_reord = new double[numpars_s];
    keps_tc_noreord = new double[numpars_s];

    #if OPENACC_ENABLED
    #pragma acc enter data create(keps_tc_reord[0:numpars_s])
    #endif

    compute_RHS_BLTC(); // e-field results saved in keps_tc_noreord

    // Copy results back to es
    for (size_t i = 0; i < numpars_s; i++) {
        es[i] = keps_tc_noreord[i];
    }

    #if OPENACC_ENABLED 
    #pragma acc exit data delete(lambda[0:numpars_s])
    #endif
    if (lambda != nullptr) {
        delete[] lambda;
        lambda = nullptr;
    }

    #if OPENACC_ENABLED
    #pragma acc exit data delete(particles_x[0:numpars_s])
    #endif
    if (particles_x != nullptr) {
        delete[] particles_x;
        particles_x = nullptr;
    }

    #if OPENACC_ENABLED
    #pragma acc exit data delete(keps_tc_reord[0:numpars_s])
    #endif

    if (keps_tc_reord != nullptr) {
        delete[] keps_tc_reord;
        keps_tc_reord = nullptr;
    }

    if (keps_tc_noreord != nullptr) {
        delete[] keps_tc_noreord;
        keps_tc_noreord = nullptr;
    }

    // Deallocate arrays of pointers, if necessary
    for (size_t i = 0; i < 2; ++i) {
        if (tree_members[i]) {
            delete[] tree_members[i];
            tree_members[i] = nullptr;
        }
        if (leaf_members[i]) {
            delete[] leaf_members[i];
            leaf_members[i] = nullptr;
        }
    }

    // Assuming max_far_size defines the number of far lists
    if (interaction_list_far) {
        for (size_t i = 0; i < max_far_size; ++i) {
            if (interaction_list_far[i]) {
                delete[] interaction_list_far[i];
                interaction_list_far[i] = nullptr;
            }
        }
        delete[] interaction_list_far;
        interaction_list_far = nullptr;
    }
    if (interaction_list_far_size) {
        delete[] interaction_list_far_size;
        interaction_list_far_size = nullptr;
    }

    // Assuming max_near_size defines the number of near lists
    if (interaction_list_near) {
        for (size_t i = 0; i < max_near_size; ++i) {
            if (interaction_list_near[i]) {
                delete[] interaction_list_near[i];
                interaction_list_near[i] = nullptr;
            }
        }
        delete[] interaction_list_near;
        interaction_list_near = nullptr;
    }

    if (interaction_list_near_size) {
        delete[] interaction_list_near_size;
        interaction_list_near_size = nullptr;
    }

// Deallocate 2D arrays or lists
if (cluster_list_t1) {
    for (size_t i = 0; i < leaf_count; ++i) {
        if (cluster_list_t1[i]) {
            delete[] cluster_list_t1[i];
            cluster_list_t1[i] = nullptr;
        }
    }
    delete[] cluster_list_t1;
    cluster_list_t1 = nullptr;
}

if (cluster_list_moments) {
    for (size_t i = 0; i < leaf_count; ++i) {
        if (cluster_list_moments[i]) {
            delete[] cluster_list_moments[i];
            cluster_list_moments[i] = nullptr;
        }
    }
    delete[] cluster_list_moments;
    cluster_list_moments = nullptr;
}
tree.clear();
leaf.clear();

}

void E_MQ_Treecode::compute_RHS_BLTC() {
    long build_tree_time = getTickCount();

    xminmax[0] = minval(particles_x, numpars_s);
    xminmax[1] = maxval(particles_x, numpars_s);

    // Reset rep_tc_reord_x
    for (size_t i = 0; i < numpars_s; i++)
        keps_tc_reord[i] = 0.0;



#if OPENACC_ENABLED
#pragma acc update device(keps_tc_reord[0:numpars_s])
#endif


    int* pt_temp_index = new int[numpars_s];
    int* pt_temp_old_index = new int[numpars_s];
    for (size_t i = 0; i < numpars_s; i++) {
        pt_temp_index[i] = -1;
        pt_temp_old_index[i] = i;
    }

    build_tree_init();

    build_tree_1D_Recursive(0, 0, particles_x, pt_temp_index, pt_temp_old_index);

#if OPENACC_ENABLED
#pragma acc update device(particles_x[0:numpars_s])
#endif

    assert(node_count == tree.size());
    cout << "Tree size = " << node_count << endl;

    assert(leaf_count == leaf.size());
    cout << "Leaf size = " << leaf_count << endl;

    //====================== Compute interaction list =================

    vector<vector<size_t>> Interaction_List_far;
    vector<vector<size_t>> Interaction_List_near;

    Interaction_List_far.resize(leaf_count);
    Interaction_List_near.resize(leaf_count);

    for (size_t leaf_index = 0; leaf_index < leaf_count; leaf_index++) {
        build_interaction_list(leaf_index, 0, Interaction_List_far, Interaction_List_near);
    }

    //=============== Copy interaction list and cluster list to the device ===============
    alloc_set_interaction_list(Interaction_List_far, Interaction_List_near);
    alloc_set_cluster_list(particles_x);

    tree_members[0] = new size_t[node_count];
    tree_members[1] = new size_t[node_count];

    for (size_t panel_index = 0; panel_index < node_count; panel_index++) {
        tree_members[0][panel_index] = tree[panel_index].members[0];
        tree_members[1][panel_index] = tree[panel_index].members[1];
    }

    leaf_members[0] = new size_t[leaf_count];
    leaf_members[1] = new size_t[leaf_count];

    for (size_t leaf_index = 0; leaf_index < leaf_count; leaf_index++) {
        leaf_members[0][leaf_index] = tree[leaf[leaf_index]].members[0];
        leaf_members[1][leaf_index] = tree[leaf[leaf_index]].members[1];
    }

#if OPENACC_ENABLED
#pragma acc enter data copyin(tree_members[0:2][0:node_count])
#pragma acc enter data copyin(leaf_members[0:2][0:leaf_count])
#endif

    Compute_SUM();

#if OPENACC_ENABLED
#pragma acc exit data delete(tree_members[0:2][0:node_count])
#pragma acc exit data delete(leaf_members[0:2][0:leaf_count])
#endif

    delete[] tree_members[0];
    delete[] tree_members[1];

    delete[] leaf_members[0];
    delete[] leaf_members[1];

//========= Change back to original order =====================
    double* lambda_temp;
    lambda_temp = new double[numpars_s];

#if OPENACC_ENABLED
#pragma acc update host(keps_tc_reord[0:numpars_s])
#endif

    for (size_t i = 0; i < numpars_s; i++) {
        keps_tc_noreord[pt_temp_old_index[i]] = keps_tc_reord[i];
        lambda_temp[i] = lambda[i];
    }

    for (size_t i = 0; i < numpars_s; i++)
        lambda[pt_temp_old_index[i]] = lambda_temp[i];

#if OPENACC_ENABLED
#pragma acc update device(lambda[0:numpars_s])
#endif

    delete[] pt_temp_old_index;
    delete[] pt_temp_index;

    // Free the information in tree and leaf vectors
    tree.clear();
    leaf.clear();
    free_cluster_list();
    free_interaction_list();
    leaf_count = 0;
    node_count = 0;
    delete[] lambda_temp;
}


void E_MQ_Treecode::build_tree_init() {
    panel temp_panel;

    // Indices of particles belonging to the panel
    temp_panel.members[0] = 0;
    temp_panel.members[1] = numpars_s - 1;

    // Interval defining the panel
    temp_panel.xinterval[0] = xminmax[0];
    temp_panel.xinterval[1] = xminmax[1];

    temp_panel.xc = 0.5 * (temp_panel.xinterval[0] + temp_panel.xinterval[1]);

    double xL = temp_panel.xinterval[1] - temp_panel.xinterval[0];

    
    double sq_r = 0.25 * (xL * xL); // r^2
    temp_panel.MAC = sq_r / sq_theta; // MAC = r^2 / theta^2

    tree.push_back(temp_panel);
    node_count = 1;
}


void E_MQ_Treecode::build_tree_1D_Recursive(size_t panel_index, int level, double* pt_x, int *pt_index, int *pt_old_index)
{
    size_t n = tree[panel_index].members[1] - tree[panel_index].members[0] + 1;
    if (n >= (size_t)N0) {
        split_tree_node(panel_index, pt_x, pt_index, pt_old_index);

        for (size_t i = 0; i < tree[panel_index].children.size(); i++) {
            size_t panel_index_new = tree[panel_index].children[i];
            build_tree_1D_Recursive(panel_index_new, level + 1, pt_x, pt_index, pt_old_index);
        }
    }
    else {
        leaf.push_back(panel_index);
        leaf_count++;
    }
}


void E_MQ_Treecode::split_tree_node(size_t panel_index, double* pt_x, int *pt_index, int *pt_old_index)
{
    split_2(panel_index, pt_x, pt_index, pt_old_index);
}



void E_MQ_Treecode::split_2(size_t panel_index, double* pt_x, int* pt_index, int* pt_old_index) {
    panel child[2];
    /*
     |Child 0 |Child 1 |
     -----------------------------> axis x
     start     mid      end
     */

    double tp_x0 = tree[panel_index].xinterval[0];
    double tp_x1 = tree[panel_index].xinterval[1];

    for (int i = 0; i < 2; i++) {
        child[i].xinterval[0] = tp_x0;
        child[i].xinterval[1] = tp_x1;
    }

    double midpoint = 0.5 * (tp_x0 + tp_x1);

    child[0].xinterval[1] = midpoint;
    child[1].xinterval[0] = midpoint;

    double xL = 0.5 * (tp_x1 - tp_x0);
    double sq_r = 0.25 * (xL * xL); // r^2
    double MAC = sq_r / sq_theta; // MAC = r^2 / theta^2

    for (int i = 0; i < 2; i++) {
        child[i].xc = 0.5 * (child[i].xinterval[0] + child[i].xinterval[1]);
        child[i].MAC = MAC;
    }

    vector<size_t> v[2];
    size_t start = tree[panel_index].members[0];
    size_t end = tree[panel_index].members[1];
    size_t* addr_table = new size_t[end - start + 1];

    size_t index;
    for (index = start; index <= end; index++) {
        pt_index[index] = index;
        addr_table[index - start] = index;
        
        if (pt_x[index] <= midpoint)
            v[0].push_back(index);
        else
            v[1].push_back(index);
    }

    size_t seq = start;
    for (size_t j = 0; j < 2; j++) {
        size_t size = v[j].size();

        if (size >= 1) {
            for (size_t k = 0; k < size; k++) {
                if (k == 0)
                    child[j].members[0] = seq;
                if (k == size - 1)
                    child[j].members[1] = seq;

                index = v[j][k];

                // This uses an address table
                size_t pos = addr_table[index - start];
                size_t out = pt_index[seq];
                Swap(pos, seq, pt_x, pt_index, pt_old_index);
                addr_table[index - start] = seq;
                addr_table[out - start] = pos;

                seq++;
            }

            node_count++;
            tree[panel_index].children.push_back(node_count - 1);
            tree.push_back(child[j]);
            v[j].clear();
        }
    }

    delete[] addr_table;
}


void E_MQ_Treecode::Swap(size_t i, size_t j, double* pt_x, int *pt_index, int *pt_old_index)
{
    if (i == j)
        return;

    double x = pt_x[i];
    size_t index = pt_index[i];
    size_t old_index = pt_old_index[i];

    double temp_lambda;
    temp_lambda = lambda[i];

    pt_x[i] = pt_x[j];
    pt_index[i] = pt_index[j];
    pt_old_index[i] = pt_old_index[j];

    lambda[i] = lambda[j];

    pt_x[j] = x;
    pt_index[j] = index;
    pt_old_index[j] = old_index;

    lambda[j] = temp_lambda;
}




void E_MQ_Treecode::build_interaction_list(size_t leaf_index, size_t panel_index,
                            vector<vector<size_t> >& Interaction_List_far,
                            vector<vector<size_t> >& Interaction_List_near)
{
    double xcb = tree[leaf[leaf_index]].xc;

    double xcp = tree[panel_index].xc;

    double x_vb = tree[leaf[leaf_index]].xinterval[1];

    double x_vp = tree[panel_index].xinterval[1];

    // double dfbpx = xcb - xcp;
    double r_dist = fabs(xcb - xcp);
    double dfbpx = fmin(r_dist, L - r_dist);

    double dfbx = x_vb - xcb;

    double dfpx = x_vp - xcp;

    double sq_R = dfbpx * dfbpx;
    double sq_rb = dfbx * dfbx;
    double sq_rp = dfpx * dfpx;

    if (tree[panel_index].children.size() == 0) // A leaf panel
        Interaction_List_near[leaf_index].push_back(panel_index);
    else {
        if ((sqrt(sq_rb) + sqrt(sq_rp)) / sqrt(sq_R) < sqrt(sq_theta))
            Interaction_List_far[leaf_index].push_back(panel_index);
        else {
            size_t size = tree[panel_index].children.size();
            if (size == 0) // It is a leaf
                Interaction_List_near[leaf_index].push_back(panel_index);
            else {
                for (size_t i_children = 0; i_children < size; i_children++)
                    build_interaction_list(leaf_index, tree[panel_index].children[i_children], Interaction_List_far, Interaction_List_near);
            }
        }
    }
}


void E_MQ_Treecode::alloc_set_interaction_list(const vector<vector<size_t> >& Interaction_List_far,
                                const vector<vector<size_t> >& Interaction_List_near)
{
    interaction_list_far = new size_t*[leaf_count];
    interaction_list_near = new size_t*[leaf_count];

    interaction_list_far_size = new size_t[leaf_count];
    interaction_list_near_size = new size_t[leaf_count];

    max_far_size = 0;
    max_near_size = 0;

    for (size_t leaf_index = 0; leaf_index < leaf_count; leaf_index++) {
        size_t far_size = Interaction_List_far[leaf_index].size();
        interaction_list_far_size[leaf_index] = far_size;
        if (far_size > max_far_size)
            max_far_size = far_size;

        size_t near_size = Interaction_List_near[leaf_index].size();
        interaction_list_near_size[leaf_index] = near_size;
        if (near_size > max_near_size)
            max_near_size = near_size;
    }

    for (size_t leaf_index = 0; leaf_index < leaf_count; leaf_index++) {
        interaction_list_far[leaf_index] = new size_t[max_far_size];
        std::copy(Interaction_List_far[leaf_index].begin(), Interaction_List_far[leaf_index].end(), interaction_list_far[leaf_index]);

        interaction_list_near[leaf_index] = new size_t[max_near_size];
        std::copy(Interaction_List_near[leaf_index].begin(), Interaction_List_near[leaf_index].end(), interaction_list_near[leaf_index]);
    }

#if OPENACC_ENABLED
    #pragma acc enter data create(interaction_list_far[0:leaf_count][0:max_far_size], \
                                  interaction_list_far_size[0:leaf_count], \
                                  interaction_list_near[0:leaf_count][0:max_near_size], \
                                  interaction_list_near_size[0:leaf_count])

    #pragma acc update device(interaction_list_far[0:leaf_count][0:max_far_size], \
                              interaction_list_far_size[0:leaf_count], \
                              interaction_list_near[0:leaf_count][0:max_near_size], \
                              interaction_list_near_size[0:leaf_count])
#endif
}


void E_MQ_Treecode::free_interaction_list() 
{
#if OPENACC_ENABLED
    #pragma acc exit data delete(interaction_list_far[0:leaf_count][0:max_far_size], \
                                  interaction_list_far_size[0:leaf_count], \
                                  interaction_list_near[0:leaf_count][0:max_near_size], \
                                  interaction_list_near_size[0:leaf_count])
#endif

    for (size_t leaf_index = 0; leaf_index < leaf_count; leaf_index++) {
        delete[] interaction_list_far[leaf_index];
        delete[] interaction_list_near[leaf_index];
    }

    delete[] interaction_list_far;
    delete[] interaction_list_near;

    delete[] interaction_list_far_size;
    delete[] interaction_list_near_size;
}




void E_MQ_Treecode::alloc_set_cluster_list(double* pt_x)
{
    cluster_list_t1 = new double*[node_count];
    cluster_list_moments = new double*[node_count];

    for (size_t tree_index = 0; tree_index < node_count; tree_index++) {
        cluster_list_t1[tree_index] = new double[PP];
        cluster_list_moments[tree_index] = new double[Pflat];
        memset(cluster_list_t1[tree_index], 0, sizeof(double) * PP);
        memset(cluster_list_moments[tree_index], 0, sizeof(double) * Pflat);
    }

#if OPENACC_ENABLED
    #pragma acc enter data create(cluster_list_t1[0:node_count][0:PP], \
                                  cluster_list_moments[0:node_count][0:Pflat])
#endif

    double h = pi / P;
    double t[PP] = {0.0};
    for (int i = 0; i < PP; i++)
        t[i] = cos(i * h); // Chebyshev interpolation points [-1, 1]

    double x1, x2;

    double w1i[PP];
    double dj[PP];
    int i, j;
    dj[0] = 0.5;
    dj[P] = 0.5;
    for (j = 1; j < P; j++)
        dj[j] = 1.0;

    for (j = 0; j < PP; j++)
        w1i[j] = ((j % 2 == 0)? 1 : -1) * dj[j];

    double a1i[PP];
    double x;
    double dx;
    double SumA1;
    double D;
    double s;
    int a1exactind;

    for (size_t tree_index = 1; tree_index < node_count; tree_index++) { // Skip tree root
        x1 = tree[tree_index].xinterval[0];
        x2 = tree[tree_index].xinterval[1];

        for (j = 0; j < PP; j++)
            cluster_list_t1[tree_index][j] = x1 + (t[j] + 1.0) * 0.5 * (x2 - x1);

        size_t tp0 = tree[tree_index].members[0];
        size_t tp1 = tree[tree_index].members[1];
        size_t tp_j;

        for (tp_j = tp0; tp_j <= tp1; tp_j++) {
            x = pt_x[tp_j];

            a1exactind = -1;

            SumA1 = 0.0;

            for (j = 0; j < PP; j++) {
                dx = x - cluster_list_t1[tree_index][j];
                if (fabs(dx) <= DBL_MIN)
                    a1exactind = j;
                else {
                    a1i[j] = w1i[j] / dx;
                    SumA1 += a1i[j];
                }
            }

            if (a1exactind > -1) {
                SumA1 = 1.0;
                for (j = 0; j < PP; j++)
                    a1i[j] = 0.0;
                a1i[a1exactind] = 1.0;
            }

            D = 1.0 / SumA1;

            for (i = 0; i < Pflat; i++) {
                s = a1i[i] * D;
                cluster_list_moments[tree_index][i] += s * lambda[tp_j];
            }
        }
    }

#if OPENACC_ENABLED
    #pragma acc update device(cluster_list_t1[0:node_count][0:PP], \
                              cluster_list_moments[0:node_count][0:Pflat])
#endif
}


void E_MQ_Treecode::free_cluster_list()
{
#if OPENACC_ENABLED
    #pragma acc exit data delete(cluster_list_t1[0:node_count][0:PP], \
                                 cluster_list_moments[0:node_count][0:Pflat])
#endif

    for (size_t tree_index = 0; tree_index < node_count; tree_index++) {
        delete[] cluster_list_t1[tree_index];
        delete[] cluster_list_moments[tree_index];
    }

    delete[] cluster_list_t1;
    delete[] cluster_list_moments;
}




void E_MQ_Treecode::Compute_SUM()
{
    long Call_BL_cpu_time = getTickCount();
    Call_BL();
    Call_BL_cpu_time = getTickCount() - Call_BL_cpu_time;

    // cout << "Call_BL_cpu_time " << Call_BL_cpu_time << endl;

    long Call_Ds_cpu_time = getTickCount();
    Call_Ds();
    Call_Ds_cpu_time = getTickCount() - Call_Ds_cpu_time;

    // cout << "Call_Ds_cpu_time " << Call_Ds_cpu_time << endl;
}


void E_MQ_Treecode::Call_BL()
{
#if OPENACC_ENABLED
#pragma acc kernels copyin(leaf_count) \
present(leaf_members[0:2][0:leaf_count], \
particles_x[0:numpars_s],  \
keps_tc_reord[0:numpars_s], \
interaction_list_far[0:leaf_count][0:max_far_size], \
interaction_list_far_size[0:leaf_count], \
cluster_list_t1[0:node_count][0:PP], \
cluster_list_moments)
{ // Begin acc kernels region
#endif
    
#if OPENACC_ENABLED
    #pragma acc loop independent
#endif
    for (size_t leaf_index = 0; leaf_index < leaf_count; leaf_index++) {
        size_t batch_limit_1 = leaf_members[0][leaf_index];
        size_t batch_limit_2 = leaf_members[1][leaf_index];
        size_t far_list_size = interaction_list_far_size[leaf_index];

#if OPENACC_ENABLED
        #pragma acc loop vector(128) independent
#endif
        for (size_t ii = batch_limit_1; ii <= batch_limit_2; ii++) {
            double tempx = 0.0;
            double p_x = particles_x[ii];

#if OPENACC_ENABLED
            #pragma acc loop seq
#endif
            for (size_t jj = 0; jj < far_list_size; jj++) {
                size_t far_index = interaction_list_far[leaf_index][jj];

#if OPENACC_ENABLED
                #pragma acc loop seq
#endif
                for (int kk = 0; kk < Pflat; kk++) {
                    double s = (p_x - cluster_list_t1[far_index][kk]) / L;
                    s = s - round(s);
                    // while (s < -0.5) { s += 1.0; }
                    // while (s >= 0.5) { s -= 1.0; }
                    tempx += cluster_list_moments[far_index][kk] * (0.5 * s * norm_epsL / sqrt(s * s + epsLsq) - s);
                } // kk
            } // jj

            keps_tc_reord[ii] += tempx;
        } // ii
    } // leaf_index
#if OPENACC_ENABLED
} // End acc kernels region
// cout << "done call_BL" << endl;
#endif
}

void E_MQ_Treecode::Call_Ds()
{
#if OPENACC_ENABLED
#pragma acc kernels copyin(leaf_count) \
present(leaf_members[0:2][0:leaf_count], \
tree_members[0:2][0:node_count], \
particles_x[0:numpars_s], \
lambda[0:numpars_s], \
keps_tc_reord[0:numpars_s], \
interaction_list_near[0:leaf_count][0:max_near_size], \
interaction_list_near_size[0:leaf_count])
    { // Begin acc kernels region
#endif
#if OPENACC_ENABLED
        #pragma acc loop independent
#endif
        for (size_t leaf_index = 0; leaf_index < leaf_count; leaf_index++) {
            size_t limit_1_b = leaf_members[0][leaf_index];
            size_t limit_2_b = leaf_members[1][leaf_index];

            size_t near_list_size = interaction_list_near_size[leaf_index];

#if OPENACC_ENABLED
            #pragma acc loop vector(128) independent
#endif
            for (size_t ii = limit_1_b; ii <= limit_2_b; ii++) {
                double tempx = 0.0;

#if OPENACC_ENABLED
                #pragma acc loop seq
#endif
                for (size_t kk = 0; kk < near_list_size; kk++) {
                    size_t limit_1_c = tree_members[0][interaction_list_near[leaf_index][kk]];
                    size_t limit_2_c = tree_members[1][interaction_list_near[leaf_index][kk]];

#if OPENACC_ENABLED
                    #pragma acc loop seq
#endif
                    for (size_t jj = limit_1_c; jj <= limit_2_c; jj++) {
                        double s = (particles_x[ii] - particles_x[jj])/L;
                        s = s - round(s);
                        // while (s < -0.5) { s += 1.0; }
                        // while (s >= 0.5) { s -= 1.0; }
                        tempx += lambda[jj] * (0.5 * s * norm_epsL / sqrt(s * s + epsLsq) - s);
                    } // jj
                } // kk

                keps_tc_reord[ii] += tempx;
            } // ii
        } // size_t leaf_index
#if OPENACC_ENABLED
    } // End acc kernels region
#endif
}





long E_MQ_Treecode::getTickCount() {
    tms tm;
    return times(&tm);
}


double E_MQ_Treecode::minval(double* x, size_t len) {
    double MinVal = x[0];
    for (size_t i = 1; i < len; i++) {
        if (MinVal > x[i]) {
            MinVal = x[i];
        }
    }
    return MinVal - 0.1234;
}

double E_MQ_Treecode::maxval(double* x, size_t len) {
    double MaxVal = x[0];
    for (size_t i = 1; i < len; i++) {
        if (MaxVal < x[i]) {
            MaxVal = x[i];
        }
    }
    return MaxVal + 0.001;
}


// void E_MQ_Treecode::cleanup() {

// }



void E_MQ_Treecode::print_field_obj() {
    cout << "-------------" << endl;
    cout << "Field object: " << endl;
    cout << "MQ kernel, treecode" << endl;
    cout << "mac = " << mac << endl;
    cout << "degree = " << degree << endl;
    cout << "max_source leaf = " << max_source << endl;
    cout << "max target leaf = " << max_target << endl;
    cout << "epsilon = " << epsilon << endl;
    cout << "Periodic boundary conditions, domain size = " << L << endl;
    cout << "-------------" << endl;
}


E_MQ_Treecode_openbcs::E_MQ_Treecode_openbcs() {}
// E_MQ_Treecode_openbcs::E_MQ_Treecode_openbcs(double epsilon, double beta) : 
//     kernel(MQ), 
//     singularity(SKIPPING), approximation(LAGRANGE), compute_type(PARTICLE_CLUSTER),
//     beta(beta), theta(-1.0), interpDegree(-1), maxPerSourceLeaf(1), maxPerTargetLeaf(1),
//     verbosity(0)
// {
//     kernelParams.push_back(-1.0); kernelParams.push_back(epsilon);
// }
// E_MQ_Treecode_openbcs::E_MQ_Treecode_openbcs(double epsilon,
//     double theta, int interpDegree, int maxPerSourceLeaf, int maxPerTargetLeaf,
//     int verbosity) :
//         kernel(MQ), 
//         singularity(SKIPPING), approximation(LAGRANGE), compute_type(PARTICLE_CLUSTER),
//         beta(-1.0), theta(theta), interpDegree(interpDegree), 
//         maxPerSourceLeaf(maxPerSourceLeaf), maxPerTargetLeaf(maxPerTargetLeaf),
//         verbosity(verbosity)
// {
//     kernelParams.push_back(-1.0); kernelParams.push_back(epsilon);
// }

// void E_MQ_Treecode_openbcs::operator() (double* es, double* targets, int nt, 
//                 double* sources, double* q_ws, int ns)
// {
//     std::vector <double> xS(ns);
//     std::vector <double> yS(ns);
//     std::vector <double> wS(ns, 1.0);

//     std::vector <double> xT(nt);
//     std::vector <double> yT(nt);
//     std::vector <double> qT(nt, 1.0);

//     BaryTreeInterface(nt, ns, xT.data(), yT.data(), 
//                     targets, qT.data(),
//                     xS.data(), yS.data(), 
//                     sources, q_ws, wS.data(),
//                     es,
//                     kernel, kernelParams.size(), kernelParams.data(),
//                     singularity, approximation, compute_type,
//                     theta, interpDegree, maxPerSourceLeaf, maxPerTargetLeaf,
//                     1.0, beta, verbosity);
//     for (int ii = 0; ii < nt; ++ii) {
//         es[ii] -= targets[ii];
//     }
// }
// void E_MQ_Treecode_openbcs::print_field_obj() {
//     cout << "-------------" << endl;
//     cout << "Field object: " << endl;
//     cout << "MQ kernel, treecode" << endl;
//     cout << "epsilon = " << kernelParams[0] << endl;
//     cout << "Open boundary conditions" << endl;
//     if (beta >= 0) {
//         cout << "treecode parameters set with beta = " << beta << endl;
//     } else {
//         cout << "mac = " << theta << ", n = " << interpDegree << endl;
//         cout << " max source = " << maxPerSourceLeaf << " , max target = " << maxPerTargetLeaf << endl;
//     }
//     cout << "-------------" << endl;
// }
    
E_MQ_Treecode_openbcs::~E_MQ_Treecode_openbcs() = default;