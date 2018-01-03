#ifndef NBO_H
#define NBO_H

#include "stringtools.h"

class NBO {
  
  private:
 
   string label;    // used for geninp and read output 
   int natoms;
//   int* anumbers;
   vector<int> anumbers;
//   string* anames;
   vector<string> anames;

//   double* xyz;
   vector<double> xyz;

   int isQC; //QC or MOPAC
   int mheadersize;
//   string* mheader;
   vector<string> mheader;

   int nalpha; 
   int nbeta;
   int nelec;
   int nao; //# atomic orbitals
   int nmo; //# molecular orbitals
   int read_mo(string filename);
   int read_mopac_mo(string filename);
   void alloc_nbo();

  public:

   void alloc(int natoms);
//   void init(int natoms0, int* anumbers0, string* anames0, double* xyz0);

   void init(int natoms0, vector<string> & anames0, vector<double> & xyz0);
//   void reset(int natoms0, int* anumbers0, string* anames0, double* xyz0);
   void freemem();

   int read_nbo_file(string filename);
   void print_nbo();

//   double* q;
   vector<double> q;
//   double* MO;
   vector<double> MO;
//   int* wAO;
   vector<int> wAO;
//   int* tAO;
   vector<int> tAO;

   int hasNBO;
   int bmo; //occupied bonding orbitals
   int vmo; //unoccupied bonding orbitals
//   double* mo_occ;
   vector<double> mo_occ;
//   double* bmo_occ;
   vector<double> bmo_occ;
//   double* bmo_polar;
   vector<double> bmo_polar;
//   int* bmo_atoms;
   vector<int> bmo_atoms;
//   int* bmo_num;
   vector<int> bmo_num;
//   double* vmo_occ;
   vector<double> vmo_occ;    
//   int* vmo_atoms;
   vector<int> vmo_atoms;
//   int* vmo_num;
   vector<int> vmo_num;

  //MOPAC's COSMO quantities
   double volume;
   double area;
 
  //MOPAC's homo lumo
  double homo;
  double lumo;
  
   void print_molden_orbs(int norbs, int* olist, string filename);
   int compare_nbo(NBO nbo1, string filename);
   double get_pol(int a1, int a2, double& ev1);

};

#endif
