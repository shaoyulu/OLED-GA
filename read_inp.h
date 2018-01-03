#ifndef READINP_H
#define READINP_H


#include "stringtools.h"
#include <string.h>
//todo: check connection error for xyzfile
//todo: generate a complete input file including all the default value


class readinp {
  
   private:

//   vector<string> vars;
//     string  rfolder;
//   vector<string> r_sym_vars;

   public:
 
   string rfolder;

   vector<vector<string> > r_sym;

   string fitness_f;            /*filename of fitness.inp*/

//   void gen_basisfile(string rfolder);

   string ref_f;
   
   int n_temp;
  
   int n_rfiles;                  /*number of files in R library*/

   vector<string> feature_list;                // indices of features we are interested in at the begining
  
   vector<string> temp_f;         /*filename of base cats*/

   void read_all();          /*read .inp & gen .basis files*/
   void read_ga(string filename);       /*read GA.inp*/

   void read_dft(string filename);      /*read DFT.inp*/

   void read_mopac(string filename);   /*read MOPAC.inp*/

//   void read_refxyzf(string inpfile);    /*read reference .xyz file*/

   void read_ref(string filename);

   void read_tempxyzf(vector<string> & inpfiles);   /*read template .xyz file(s)*/

   void get_basisfile(string filedir);

   void read_seed(string filename);
   
   int load_lib(int n, string libname);

//store variables in GA.inp
   struct INP_SGEN
   {
     int mode;
     int nlayers;
     int layer_min;
     int layer_max;
     string lib;
     int n_regs;
     int n_caps;
     int n_seeds;
     bool symm;                    
     string seed;            //filename of seed
     double HOST_HOMO, HOST_LUMO, ETL_HOMO, ETL_LUMO, HTL_HOMO, HTL_LUMO;
     int jobtype =1;
     bool bestRestart;
     bool isWriteFragScore;
    double fscore;
     int TransportType =1 ;
    } inp_sgen;

   struct SEED
   {
     int nhandles;
     int nbranches;    // different branches
     int natoms;
     vector<int> xyzw;
     vector<string> label;
     vector<string> name;
     vector<double> coords;
   } seed; 

   struct REF_XYZ
   {
     int natoms;
     vector<string> anames;
     vector<double> coords;
   } ref_xyz;

   struct INP_GA
   {
   int ox_state;
   int num_r;
   int use_dft;
   int use_mopac;
   int pop_size;
   int ga_stop;
   double cross_r;
   double muta_r;
   int out_molden;
   int out_excel;
   bool RWS;
   bool SUS; 
   bool B_Tournament;
   bool L_Rank;
   double L_Rank_rate;
   bool E_Rank;  
   double E_Rank_base;
   bool force_mutation; 
   
   INP_GA()
   {
        RWS = false;
        SUS = false;
        B_Tournament = false;
        L_Rank = false;
        E_Rank = false;
        force_mutation = true;
   } 
   
   } inp_ga;

   struct xyzfile
   {
          int natoms;
          int nlabels;
          int charge;
          int multi;
	  vector<string> element;
	  vector<double> coords;
	  vector<vector<string> > label;
	  vector<vector<int> > joint;
       void reset()
       {
       natoms = -9999;
       nlabels = -9999;
       charge = -99999;
       multi = -9999;
       vector<string>().swap(element);
       vector<double>().swap(coords);
       vector<vector<string> >().swap(label);
       vector<vector<int> >().swap(joint);
       }
    };
    xyzfile xyzf_ref;
    vector<xyzfile> xyzf_temp;    

    struct INP_MOPAC {
    
    string method;
    string walltime;
    bool use_aux;
    double gnorm;
    bool uhf;
    } mopac;

//store variables in DFT.inp
    struct singlepoint {

    string functional;
    string unres;
    string scf_algo;
    int scf_max_cycles;
    string basis;
    string wavefunction_analysis;

    int scf_conv;

    string sym_ignore;
    string sym;

    void initial() {

    functional = "B3LYP";      //set B3LYP as default functional
    unres = "TRUE";
    scf_algo = "rca_diis";
    scf_max_cycles = 150;
    basis = "LANL2DZ";        // set LAN2DZ as default basis set
    wavefunction_analysis = "FALSE";
    scf_conv = 6;

    sym_ignore = "TRUE";
    sym = "FALSE";
    }

    } sp;

//
//
  struct optimization {

   string functional;
   string unres;
   string scf_algo;
   int scf_max_cycles;
   string basis;
   string wavefunction_analysis;
   int opt_max_cycles;
   int opt_tol_displacement;
   int opt_tol_gradient;
   int opt_tol_energy;
   int scf_conv;

   string sym_ignore;
   string sym;

   void initial() {
       functional = "B3LYP";      //set B3LYP as default functional
       unres = "TRUE";
       scf_algo = "rca_diis";
       scf_max_cycles = 150;
       basis = "LANL2DZ";        // set LAN2DZ as default basis set
       wavefunction_analysis = "FALSE";
       scf_conv = 6;
       opt_max_cycles = 300;
       opt_tol_displacement = 2500;
       opt_tol_gradient = 800;
       opt_tol_energy = 5000;
       sym_ignore = "TRUE";
       sym = "FALSE";
       }

    } opt;

//
//
  struct transitionstate {

    string functional;
    string unres;
    string scf_algo;
    int scf_max_cycles;
    string basis;
    string wavefunction_analysis;
    int scf_conv;

    string sym_ignore;
    string sym;

    void initial() {
        functional = "B3LYP";      //set B3LYP as default functional
        unres = "TRUE";
        scf_algo = "rca_diis";
        scf_max_cycles = 150;
        basis = "LAN2DZ";        // set LAN2DZ as default basis set
        wavefunction_analysis = "FALSE";
        scf_conv = 6;
        sym_ignore = "TRUE";
        sym = "FALSE";
        }

     } ts;


	  
};


#endif
