#ifndef FRAGMENT_H
#define FRAGMENT_H

#include "read_inp.h"
#include "mopac.h"
#include "nbo.h"

   struct F_feature{
     string descriptor;
     int idx;
     double val;
     };

class Fragment{    

    public:

    Fragment(readinp* pInp);             //constructor

    private:

    int PLG;
    int PLG_nbr;

    readinp* pInp;
    bool dft;
    bool mopac;

    void read_geo();

    void read_nbo_mopac(string dir);  //geninp->run->run
//    void read_geo_and_nbo(ifstream & infile);

//   void read_xyzfile(string filename);

    void opt_mopac(string dir);  //geninp->run->read

    void opt_dft();
    public:

//rgroup is saturated by H atom

     int natoms;
     string name;
     int index;
     vector<string> anames;
     vector<double> coords;
     vector<F_feature> f_props;      //list of properties; 

    void init(string path,int i);  //i: index of R group = i+1        

    void print_props_csv(ofstream & output);       //print features in csv format         
    void print_fragment_prop();  //print all properties
    void print_fragment_prop(int idx);   //print one property with feature index = idx
    double read_prop(int prop_idx);
    double read_prop(string descriptor);
    
    void get_all_props(string dir);    //read all desired properties based on inpfile

};

#endif
