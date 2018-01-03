#include "fragment.h"

#include <openbabel/obconversion.h>
#include <openbabel/mol.h>
#include <openbabel/builder.h>
#include <openbabel/obiter.h>
#include <memory>
#include "stringtools.h"

using namespace std;

Fragment::Fragment(readinp* pointer)
{
   this->pInp = pointer;
}

//initialize anames & coords & name from xyzfile
void Fragment::init(string path, int i) {

    string filename;
    ifstream infile;
    string line;
    vector<string> tok_line;
 
    
       this->index = i+1;

       this->name = "fragment" + StringTools::int2str(i+1,1,"0");
       filename = path + StringTools::int2str(i+1,1,"0") + ".xyz";

       infile.open(filename.c_str());
    
       
       if (!infile)
       {
          cout << "cannot open " << filename << endl;
          exit(-1);
       }
        
        getline(infile,line);
        natoms = atoi(line.c_str());
        getline(infile,line);
//        tok_line = StringTools::tokenize(line," ");
//        this->PLG = atoi(tok_line[0].c_str());
//        this->PLG_nbr = atoi(tok_line[1].c_str());

//        name = StringTools::newCleanString(line);
        
        for (int j=0; j < natoms; j++)
        {
             getline(infile, line);
             tok_line = StringTools::tokenize(line, " ");
             anames.push_back(tok_line[0]);
             coords.push_back(atof(tok_line[1].c_str()));
             coords.push_back(atof(tok_line[2].c_str()));
             coords.push_back(atof(tok_line[3].c_str()));
             if (tok_line[4].find("PLG") != string::npos)
            {
                this->PLG = j+1;
            }
        } 

      infile.close();
//      cout << "initialize R_group based on " << pInp->rlist[i] << endl;  
      unique_ptr<OpenBabel::OBMol> mol(new OpenBabel::OBMol);

      OpenBabel::OBConversion conv;
      conv.SetInFormat("xyz");
      conv.ReadFile(mol.get(),filename);

      OpenBabel::OBAtom* plg = mol->GetAtom(this->PLG);
      FOR_NBORS_OF_ATOM(nbr,plg)
      {
         this->PLG_nbr = nbr->GetIdx();   
      }
     return;
}

void Fragment::get_all_props(string dir)
{
    opt_mopac(dir);
    read_nbo_mopac(dir);       
}

//hard-coded!!!!!
void Fragment::read_nbo_mopac(string dir)
{
   NBO nbo;
   nbo.init(natoms, anames, coords);
   
    string s_index = StringTools::int2str(this->index, 4, "0");
    string nbofile = dir + "mopac/" + s_index + "-mopt.out";

    nbo.read_nbo_file(nbofile.c_str());
    
    F_feature f1;

    f1.descriptor = "polar";
    double ev1 = -1000.0;
    f1.val = nbo.get_pol(this->PLG,this->PLG_nbr,ev1);
    f1.idx = 0;
    f_props.push_back(f1);

#if 0  
    F_feature f2;
    
    f2.descriptor = "volume";
    f2.val = nbo.volume;
    f2.idx = 1; 
    f_props.push_back(f2);
#endif    
    F_feature f2;
    f2.descriptor = "HOMO";
    f2.val = nbo.homo;
    f2.idx = 1;
    f_props.push_back(f2);

    F_feature f3;
    f3.descriptor = "LUMO";
    f3.val = nbo.lumo;
    f3.idx = 2;
    f_props.push_back(f3);              
}

void Fragment::opt_mopac(string dir)
{

    string s_index = StringTools::int2str(this->index, 4, "0");
    string inpfile = dir + "mopac/" + s_index + "-mopt";      
    
    if (!Mopac::opt_fragment(inpfile, pInp, anames, coords))
    {
       cout << "There is a problem to optimize Fragment: " << name << " through MOPAC" << endl;
       exit(-1);
    }

    return;     
}

double Fragment::read_prop(int prop_idx) {

     double prop = f_props.at(prop_idx-1).val;
     return prop;
}

double Fragment::read_prop(string descriptor) {

     vector<F_feature>::iterator it;
     for (it = f_props.begin(); it != f_props.end(); it++)
     {
          if ((*it).descriptor == descriptor)
          return (*it).val;
     }  
     return 0.;
}

void Fragment::print_fragment_prop() {

     printf("Fragment[%2i] %10s: \n", index, name.c_str());   
     vector<F_feature>::iterator it;
      
     for (it = f_props.begin(); it != f_props.end(); it++)
     {
         printf("%10s   %8.5f \n", (*it).descriptor.c_str(), (*it).val);
     }    
}

void Fragment::print_fragment_prop(int idx) {

      printf("Fragment[%2i] %10s: ", index, name.c_str());
      F_feature & tmp = f_props[idx-1]; 

          printf("%10s:   %8.5f \n", tmp.descriptor.c_str(), tmp.val);
      }
//one comma before newline
void Fragment::print_props_csv(ofstream & output) {

    vector<F_feature>::iterator it;
    for (it = f_props.begin(); it != f_props.end(); it++)
    {
        output << (*it).val << "," ;
    }
//    output << endl;
}
