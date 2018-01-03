#include <iostream>
#include <fstream>
#include <stdio.h>
#include <sys/types.h>
#include <sys/stat.h>
using namespace std;
#include "read_inp.h"
#include "icoord.h"
#include "genetic.h"  
#include "simple_genetic.h"
#include "mol2xyz.h"

int main() {

  cout.setf(ios::fixed, ios::floatfield);
/* 
  string molfile = "1.mol";

  MOL2XYZ  mol2xyz;
 
  unique_ptr<OpenBabel::OBMol> mol(mol2xyz.readMOL(molfile));

  cout << "finish reading and manipulating" << endl;

  mol2xyz.title = "1";

  string dir = "/export/zimmerman/yzhao/S-Generator/runs/MOL2XYZ/";

  mol2xyz.writeXYZ(mol,dir);

  cout << "finish writing xyz file" << endl;
  exit(-1);
*/   

  readinp* pInp = new readinp();
//  LinkedList* pList = new LinkedList();

  pInp->read_all();

#if 0
  string ga_inp;
 
  ga_inp="ga.inp";

  cout << "*****************************" << endl;
  cout << "**   Reading Input Files   **" << endl;
  cout << "*****************************" << endl;

  pInp->read_ga(ga_inp);

  cout << "\n" << endl;

  pInp->read_refxyzf(pInp->ref_f);

  cout << "\n" << endl;

  pInp->read_tempxyzf(pInp->temp_f);

  cout << "\n" << endl;

  pInp->get_basisfile(pInp->rfolder);
 
  cout << "totally " << pInp->n_rfiles << " R basis" << endl;
#endif
#if 0
  xyznode* pRef = new xyznode();
  pRef->charge = pInp->xyzf_ref.charge;
  pRef->multi = pInp->xyzf_ref.multi;
  pRef->natoms = unsigned(pInp->xyzf_ref.natoms);
  pRef->element = pInp->xyzf_ref.element;
  pRef->coords = pInp->xyzf_ref.coords;
  pRef->ntemp = 0;
  pRef->next = NULL;
  pRef->comment = "Ref";

  cout << "*******All about reference cat*******" << endl; 
  cout << "Charge: " << pRef->charge << endl;
  cout << "Multiplicity: " << pRef->multi << endl;
  cout << "Number of Atoms: " << pRef->natoms << endl;
  cout << "Cartesian coordinate: " << endl;
  
  for (size_t i = 0; i < pRef->natoms; i++)
  {
     cout << pRef->element[i] << " " << pRef->coords[3*i+0] << " " << pRef->coords[3*i+1] << " " << pRef->coords[3*i+2] << endl;
  }

  cout << endl;
 
//  pList->addnode(pRef);
#endif

//  string filename;
//  filename = "DFT.inp";
//  pInp->read_dft(filename);

//  filename = "MOPAC.inp";
//  pInp->read_mopac(filename);

#if 0
  Fitness fitness(pList,pInp); 
  fitness.read_fitness(pInp->fitness_f);
  cout << "**************************" << endl; 
  cout << "Done reading fitness_f.inp" << endl;
  cout << "**************************" << endl;

  fitness.geninp_Binding(pRef);
  fitness.geninp_Haffinity(pRef);
  fitness.geninp_Hydricity(pRef);
#endif


    S_Genetic sgenetic(pInp);
    sgenetic.run_ga();

#if 0
//test case : 3 branches on seed, 1 socket on each adding block, 2 layers
  int i;
  int j;
  int m;
  int n;
  int x;
  int y;
  int nfrags = pInp->inp_sgen.nblocks;
  int ncaps = pInp->inp_sgen.nblocks_endp;
  string code;
  int count=0;

//  gTree gtree(pInp);
  for (i=0; i< nfrags; i++) 
  {
   string istr = StringTools::int2str(i+1,1,"0");
   for (j=0; j < ncaps; j++)
   {
    string jstr = StringTools::int2str(j+1+nfrags,1,"0");
//    for (m=0; m < nfrags; m++)
//    {
//     string mstr = StringTools::int2str(m+1, 1, "0");
//     for (n=0; n < ncaps ; n++)
//     {
//      string nstr = StringTools::int2str(n+1, 1, "0");
//      for (x=0; x < nfrags; x++)
//      {
//       string xstr =StringTools::int2str(x+1, 1, "0");
//       for (y=0; y < ncaps; y++)
//       {
//         gTree gtree(pInp);

//         string ystr = StringTools::int2str(y+1, 1, "0");
           code = "0(" + istr + "-1(" + jstr + "-1))";
//         code = "0("+ istr + "-1(" + jstr + "-1))(" + mstr + "-2(" + nstr + "-1))(" + xstr + "-3(" +ystr + "-1))";
         count++;   
        gTree gtree(pInp);      
        gtree.genecode = code;
        string cwd = StringTools::getcwd_str();
        string mmopt_folder = cwd + "/Scratch/MMOPT/";
        mkdir(mmopt_folder.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
        string nstr = StringTools::int2str(count, 4, "0");
        string cname = mmopt_folder + "mmff94.xyz" + nstr;

        ofstream xyzfile(cname.c_str());
        gtree.gene2struc(gtree.genecode,xyzfile);
        xyzfile.close();

        ifstream output(cname.c_str());
//        gtree.get_xyz(output);
        GANode

//        gtree.print_xyz("sgentest.xyz");
//       } 
//      }
//     }
//    }
   }
  } 
//  gtree.installation();

//  gtree.gene2struc(gtree.genecode,cname);
#endif
    

//  xyznode* pNode = new xyznode();
//  pNode->natoms = 163;
//  pNode->comment = "0001";
//  pNode->element.resize(pNode->natoms);
//  pNode->coords.resize(3*pNode->natoms);  
//  string filedir = "./";
//  DFT::save_opted_xyz_energy(filedir,pNode); 

  
  cout << "Done here" << endl;

//  pList->print_xyzfile();

//  pList->clear();

  delete pInp;
//  delete pNode;  
//  delete pList;
//  delete pRef;
  

  return 0;
  
}

