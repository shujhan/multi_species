#ifndef FIELD_STRUCTURE_HPP
#define FIELD_STRUCTURE_HPP

#include <algorithm>
#include <iterator>
#include <math.h>
#include <omp.h>
#include <iostream>
using std::cout;
using std::endl;
#include <vector>


#include <cstring>
#include <cfloat>
#include <cstddef>
#include <iomanip>
#include <cmath>
#include <sys/times.h>
#include <cassert>
#include <fstream>

using namespace std;

// Helper struct definition for panel
struct panel
{
    size_t members[2];
    double xinterval[2];
    double xc; // Panel center x coordinate
    double MAC; // r^2 / theta^2
    std::vector<size_t> children;

    // Initialization
    panel() : xc(0.0), MAC(0.0)
    {
        std::memset(members, 0, sizeof(members));
        std::memset(xinterval, 0, sizeof(xinterval));
    }
};

class ExternalElectricField {
    public:
        virtual void operator() (double* es, double* targets, int nt, double t) = 0;
        virtual void print_field_obj() = 0;
        virtual ~ExternalElectricField();
};

class ExternalTanh : public ExternalElectricField {
    public:
        double Am, k, omega;
        double t0 = 0, tL = 69, tR = 307, twL = 20, twR = 20;
        double eps;
        double drive_fun(double t);
        ExternalTanh();
        ExternalTanh(double Am);
        ExternalTanh(double Am, double k, double omega);
        void operator() (double* es, double* targets, int nt, double t);
        void print_field_obj();
        ~ExternalTanh();
};

class ExternalSine : public ExternalElectricField {
    public:
        double Am, k, omega;
        double t1 = 50, t2 = 150, t3 = 200;
        ExternalSine();
        ExternalSine(double Am);
        ExternalSine(double Am, double k, double omega);
        void operator() (double* es, double* targets, int nt, double t);
        void print_field_obj();
        ~ExternalSine();
};

class ExternalLogistic : public ExternalElectricField {
    public:
        double Am, k, omega;
        double t1 = 60;
        ExternalLogistic();
        ExternalLogistic(double Am);
        ExternalLogistic(double Am, double k, double omega);
        void operator() (double* es, double* targets, int nt, double t);
        void print_field_obj();
        ~ExternalLogistic();
};

class ElectricField {
    public: 
        virtual void operator()     (double* es, double* targets, int nt, 
                                    double* sources, double* q_ws, int ns) = 0;
        virtual void print_field_obj() = 0;
        virtual ~ElectricField();
};

class E_MQ_DirectSum : public ElectricField {
    public:
        E_MQ_DirectSum();
        E_MQ_DirectSum(double L, double epsilon);
        double epsilon;
        double L;
        void operator() (double* es, double* targets, int nt, 
                        double* sources, double* q_ws, int ns);
        void print_field_obj();
        ~E_MQ_DirectSum();
};

class E_MQ_DirectSum_openbcs : public ElectricField {
    public:
        E_MQ_DirectSum_openbcs();
        E_MQ_DirectSum_openbcs(double epsilon);
        double epsilon;
        void operator() (double* es, double* targets, int nt, 
                        double* sources, double* q_ws, int ns);
        void print_field_obj();
        ~E_MQ_DirectSum_openbcs();
};

class E_MQ_Treecode : public ElectricField {
    public:
        E_MQ_Treecode();
        E_MQ_Treecode(double L, double epsilon, double beta);
        E_MQ_Treecode(double L, double epsilon,
                  double mac, int degree, int max_source, int max_target,
                  int verbosity);
        void operator() (double* es, double* targets, int nt, 
                        double* sources, double* q_ws, int ns) override;
        void print_field_obj();
        ~E_MQ_Treecode();

    // private:
        void compute_RHS_BLTC();
        void cleanup();
        double minval(double* x, size_t len);
        double maxval(double* x, size_t len);

        void build_tree_init();
        void build_tree_1D_Recursive(size_t panel_index, int level, double* pt_x, int* pt_index, int* pt_old_index);
        void build_interaction_list(size_t leaf_index, size_t panel_index,
                                    std::vector<std::vector<size_t>>& Interaction_List_far,
                                    std::vector<std::vector<size_t>>& Interaction_List_near);
        void Swap(size_t i, size_t j, double* pt_x, int* pt_index, int* pt_old_index);
        void split_tree_node(size_t panel_index, double* pt_x, int* pt_index, int* pt_old_index);
        void split_2(size_t panel_index, double* pt_x, int* pt_index, int* pt_old_index);
        void alloc_set_interaction_list(const std::vector<std::vector<size_t>>& Interaction_List_far,
                                        const std::vector<std::vector<size_t>>& Interaction_List_near);
        void free_interaction_list();
        void alloc_set_cluster_list(double* pt_x);
        void free_cluster_list();
        long getTickCount();
        void Compute_SUM();
        void Call_BL();
        void Call_Ds();

        double beta;
        double L;
        double epsilon;
        double mac;
        int degree;
        int max_source;
        int max_target;
        int verbosity;

        // Other member variables...
        int P; // Order of Far-field approximation
        int PP;
        int Pflat;
        size_t numpars_s;
        size_t N0;
        double sq_theta;
        double norm_epsL;
        double epsLsq;
        double pi = 3.14159265358979323846;

        double* lambda;
        double* particles_x;
        double* keps_tc_reord;
        double* keps_tc_noreord;
        
        // Following the variables from the previous implementation...
        size_t node_count = 0;
        size_t leaf_count = 0;
        std::vector<panel> tree;
        std::vector<size_t> leaf;
        size_t* tree_members[2];
        size_t* leaf_members[2];
        size_t** interaction_list_far;
        size_t* interaction_list_far_size;
        size_t max_far_size;
        size_t** interaction_list_near;
        size_t* interaction_list_near_size;
        size_t max_near_size;
        double** cluster_list_t1;
        double** cluster_list_moments;

        double xminmax[2];
        // double* lambda_temp; 
};

class E_MQ_Treecode_openbcs : public ElectricField {
    public:
        // KERNEL kernel;
        // std::vector<double> kernelParams;
        // SINGULARITY singularity;
        // APPROXIMATION approximation;
        // COMPUTE_TYPE compute_type;
        // double beta, theta;
        // int interpDegree, maxPerSourceLeaf, maxPerTargetLeaf;
        // int verbosity;

        // E_MQ_Treecode_openbcs();
        // E_MQ_Treecode_openbcs(double epsilon, double beta);
        // E_MQ_Treecode_openbcs(double epsilon,
        //     double theta, int interpDegree, int maxPerSourceLeaf, int maxPerTargetLeaf,
        //     int verbosity);
        // void operator() (double* es, double* targets, int nt, 
        //                 double* sources, double* q_ws, int ns);
        // void print_field_obj();
        // ~E_MQ_Treecode_openbcs();
        E_MQ_Treecode_openbcs();
        ~E_MQ_Treecode_openbcs();


};

#endif /* FIELD_STRUCTURE_HPP */