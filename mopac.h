#ifndef MOPAC_H
#define MOPAC_H

#include "stringtools.h"
#include "read_inp.h"
#include "node.h"
#include "linkedlist.h"
//#include "icoord.h"
//#include "CatsGen.h"

#define MAX_TIME_WAIT 25000

namespace Mopac {

  void opt_header(ofstream& inpfile, readinp* pInp, Node<GANode>* pNode);
    int opt_fragment(string mopac_folder, readinp* pInp,
                        vector<string> & anames, vector<double> & coords);
//  void write_ic_input(ofstream & inpfile, xyznode* pNode);
   
//  double opt();
//  double opt(string filename);
//  static int opt(string mopac_folder, readinp* pointer, xyznode* pNode, int* frezlist);
    void mopac_genqsh(string filedir, vector<int> & list);
   void mopac_genslurm(string filedir, vector<int> & list); 
   int gen_inp(string mopac_folder, readinp* pInp, Node<GANode>* pNode);
//  void opt_write();
//  void opt_write(string filename);
//  void opt_write(string filename, readinp* pointer, xyznode* pNode);

//  double read_output(string filename);
    void readout(string mopac_folder,Node<GANode>* pNode);
    int save_all_into_list(string filedir, LinkedList<Node<GANode> >* pList);
//    int hard_code_save_alltolist(string filedir, LinkedList<Node<GANode> >* pList);
    void xyz_read(string filename, vector<string> & anames, vector<double> & coords);
//    void xyz_read(string mopac_folder, Node<GANode>* pNode);
    int check_success(string filename);
    bool run_opt_jobs(string filedir, vector<int> & list, int jobtype);

 //  void xyz_read_aux(string filename);
//  void xyz_save(string filename);

};

#endif
