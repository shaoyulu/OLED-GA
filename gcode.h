#ifndef gCode_H
#define gCode_H

#include "stringtools.h"
//#include "icoord.h"
#include "read_inp.h"
//#include "OBabel.h"
#include <deque>
#include <openbabel/builder.h>
#include <openbabel/obiter.h>
#include <memory>
#include <openbabel/forcefield.h>
#include <openbabel/mol.h>
#include <openbabel/obconversion.h>

struct TNode {

   int label;           //label in the lib                
   int c_index;
   //only save index in tree (so none of two nodes has same val )
   int tag;
   TNode* firstchild;
   TNode* nextsibling;
   };

struct data {

   int label;         // label in the lib
   int c_index;         // index in the cArray (for use of gen new code)
   int socket;       // ->pos on parent
   int plug;         //pos->parent  on this node
   int acceptor;        //which it docks to
   int layer;           //belongs to i-th layer
//   int nhandles;
   int tag;                //belongs to X/Y/Z/W... branches
   vector<int> xyzw;         // only for seed < x y z w > = <3> only 3 handles labeled as X
   int natoms;
//   vector<string> name;
//   vector<double> coords;
   vector<string> label_h;       //label of chemical handle
   };

class gTree {                   //binary tree representation of genetic code
    
//   int nnode;
  
//   int nlayers;
 
   readinp* pInp;
   
   void cpynode(data &newnode, const data oldnode);
//   vector<data> cArray;        //conversion array for tree building

/*for high symmetry ONLY now*/     
   void parser (string code);          //parse genetic code  
   void save_node(vector<string>::iterator & itv, int layer);
   TNode* Tnode(int index);
   TNode* root;               //point to root of tree

   TNode* create(int & i);       //add i-th node into tree and return its pointer         
// int readlib(data & bblock,int libn);    //return position of PLG
//   int readlib_lite(int libn);   //just return n different branches

   int readelib(data & bblock,int elibn);   //read end-point lib
//   int getpos(data & dockto);    //return position of SOKT on parent
   int getpos(data & seed, int index); //resturn position of SOKT with index i 

//  void destory(TNode* p);         
   void getParent(TNode* pNode, bool & find, int index, TNode* & parent);
   void tree2code(TNode* pNode, string & code);
   vector<int> serialize_tree(TNode* root, int height);
   void getNode(TNode* pNode, bool & find, int index, TNode* & target);
//   void findallKids();  //get all kids for a Parent in tree structure (not "binary" tree)
//   void install(TNode* pNode, ICoord & ic);      //follow order of  preOrder traversal
//   void postOrder();
//   void inOrder();
//   void merge1(int t1, int t2, ICoord & ic1, ICoord & ic2);
//   void merge_ic(ICoord & ic1, int socket, ICoord & ic2, int plug, ICoord & ic3);
//   void dock(int tracking, vector<data>::iterator it, ICoord & ic,
//             vector<int> & addorder, vector<vector<int> > & newpos);     
   // dock one block to another
//   void newbranch(data & seed, int tag, int socket);
     void newbranch(data & parent, int  i, int nlayers);
     int max_depth(TNode* root);
//   void newbranch_symm(data & seed, int socket);
//   void copybranch(vector<data> & tmp);
//   void install_B(ICoord & ic, vector<data>::iterator & itc);  //install a mainbranch with "tag" from cArray
//  unique_ptr<ICoord> create_icptr(int natoms);   //allocating memory for ICoord class as well
   
public:

   int nlayers;

   int nnodes;

   vector<data> cArray;
   vector<data> tempC;

   string genecode;

   void get_xyz(ifstream &instream);  //get coords from xyz file  

   vector<string> anames;
   vector<double> coords;

   gTree(readinp* pointer);               //constructor
   ~gTree(void);                          //destructor    

   vector<data> getcArray();
   void print_xyz(string filename);
   int readcode(string code);
   int readcode(string code, vector<data> tempC);    //generate cArray from genetic code
   int newStruct();         //generate a cArray for a new structure randomly
   string printCode(); 
   string printCode(vector<data> cArray);      //generate genetic code from cArray
   void destory(TNode* p);
   void destoryTree();
   int maxDepth();       //return hight of the tree
   TNode* get_root();
   string tree_to_code();
   vector<int> serialize(int height);   //serialize binary tree in a complete-binary tree form     

   TNode* getNodeTree(int index);  // return a ptr which pointer to node with c_index = index

   void createTree();
// void installTree();   
   unique_ptr<OpenBabel::OBMol> install();         //install the molecule
  
//   void install_All(ICoord & ic);        //install all branches to seed based on cArray   
//   void postOrderTraverse();
//   void inOrderTraverse();
   TNode* getParentTree(int index);   //return a p which point to parent
   void add_branch(TNode* & pNode, int tag, int nbranches, int layer); 
   int readlib_lite(int libn);   //just return n different branches
   
    int readlib(data & bblock,int libn);  
   void installation();
   void gene2struc(string code, ofstream& ostream);
//   static unique_ptr<ICoord> create_icptr(int natoms); //allocating memory for ICoord class as well
   void gene2struc2(string code, ofstream& ostram, vector<data> cArray);
   void align_to_x(int natoms, int t1, int t2, vector<double> & xyz); 
//   static void merge1(int t1, int t2, const unique_ptr<ICoord> & ic1, const unique_ptr<ICoord> & ic2);
//   static void merge1(int t1, int t2, ICoord & ic1, ICoord & ic2);
//   static unique_ptr<ICoord> merge_ic(const unique_ptr<ICoord> & ic1, int socket, const unique_ptr<ICoord> & ic2, int plug);
   unique_ptr<OpenBabel::OBMol> GetMol(string filename);

   void merge_mol(unique_ptr<OpenBabel::OBMol> & mol1, int socket, unique_ptr<OpenBabel::OBMol> & mol2, int plug);

   int minimize(unique_ptr<OpenBabel::OBMol> & mol);

//  static void merge_ic(ICoord & ic1, int socket, ICoord & ic2, int plug, ICoord & ic3);
//   static void dock(vector<data>::iterator it, ICoord & ic,
//             vector<int> & addorder, vector<vector<int> > & newpos);
   // dock one block to another   
  
};
 



#endif
