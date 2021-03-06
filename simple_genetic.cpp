#include <fstream>
#include <stdio.h>
#include <stdlib.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <dirent.h>
#include <math.h>
#include <string>

using namespace std;
#include "pTable.h"
#include "simple_genetic.h"
#include "mopac.h"
#include "gcode.h"
#include "fragment.h"
//#include "icoord.h"
//#include "xyzlist.h"
//#include "CatsGen.h"
//#include "DFT.h"
#include "utils.h"
int gcount=0;
bool isFirstrun=true;
double bestsofar, bestsofarHTL, bestsofarHOST =0;
double bestHOMOETL, bestHOMOHTL, bestLUMOETL, bestLUMOHTL=0;
string bestgene, bestgeneHTL, bestsofarHOSTstr="";
vector<string> testRandom, databanks, tempdata;
double refHOST_HOMO, refHOST_LUMO, refHOST_score=0;
vector<string> topoled, thisoleds;

vector<string> tempstring_ETL,tempstring_HTL;
string ptr1string,ptr2string;
string strFrag[] = {"quinazoline","pyridine","triazine","benzimidazole_1",
                    "1,3-xylene","dibenzofuran","1,4-xylene","naphthalene",
                   "diphenylfluorene","triphenylamine","N-phenylcarbazole_1","N-phenylcarbazole_2",
                   "dihydroacridine","quinazoline_1","1,5-naphthyridine",
                   "benzimidazole_2","acridine","quinazoline_2","benzimidazole","N-phenylcarbazole",
                  "phenoxazine","dimethylfluorene","biphenyl_1","biphenyl_2","nitrile",
                  "triphenylamine","pyridine","naphthalene_1","naphthalene_2",
                  "haloalkane","phenyl","methyl"};


vector<Node<GANode>* >  ptr_array_HTL_temp;
vector<Node<GANode>* >  ptr_array_ETL_temp;
vector<Node<GANode>* > pTmp2;
#define NUM_ALL_PERM 9450 

void S_Genetic::fragment_prop()
{
    string path = "./Scratch/features/mopac";
    mkdir(path.c_str(),S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
//  string dir = pInp->inp_sgen.lib;
    string dir = "./Scratch/features/";

    int nreg = pInp->inp_sgen.n_regs;
    int ncap = pInp->inp_sgen.n_caps;

    for (int i=0; i < nreg + ncap; i++)
    {   
  

       Fragment* p_frag = new Fragment(pInp);

       p_frag->init(dir,i);
       p_frag->get_all_props(dir);

       pFrags.push_back(p_frag);
    }

    ofstream csv_reg;
    string csv_reg_fname = "fragment_reg_prop.csv";
    csv_reg.open(csv_reg_fname.c_str());
    csv_reg << setprecision(6);

    ofstream csv_cap;
    string csv_cap_fname = "fragment_cap_prop.csv";
    csv_cap.open(csv_cap_fname.c_str());
    csv_cap << setprecision(6);

    csv_cap << "frag_index,";
    csv_reg << "frag_index,";

    vector<string>::iterator it_flist;
    for (it_flist = pInp->feature_list.begin(); it_flist < pInp->feature_list.end(); it_flist++)
    {
        csv_cap << *it_flist << ",";
        csv_reg << *it_flist << ",";
    }    
    csv_cap << endl;
    csv_reg << endl;    

   cout << "Calculated Properties of Fragments:" << endl;

    vector<Fragment*>::iterator it;
   int i = 0;
   for (it = pFrags.begin(); it != pFrags.end(); it++)
   { 
       i++;
       if (i == 1)  
         cout << "regular fragments:" << endl;
       if (i == nreg+1)
         cout << "cap fragments:" << endl;
       
       (*it)->print_fragment_prop();
       if (i <= nreg)
        {
        csv_reg << i << ",";
       (*it)->print_props_csv(csv_reg);
        csv_reg << endl;
        }
       else
        {
        csv_cap << i << ",";
       (*it)->print_props_csv(csv_cap);
        csv_cap << endl;
        }
   }

    cout << endl << endl;
    csv_reg.close();
    csv_cap.close();
#if 0
vector<int>::iterator it_i;
for (it_i = pInp->feature_list.begin(); it_i < pInp->feature_list.end(); i++)
{
   vector<Fragment*> tmp(pFrags);
   auto comp = [](const Fragment* x, const Fragment* y)
   {return x->f_props[*it_i].val < y->f_props[*it_i].val ; };

   sort(tmp.begin(), tmp.end(), comp);

   for (it = tmp.begin(); it != tmp.end(); it++)
   {
      (*it)->print_frag_prop(*it_i);    // sort by "pol" and print (which feature index is 1     )
   }
}
#endif
   return;   
}

//cluster all candidates in the given feature space the for first sampling 
void S_Genetic::init_feature_space()
{
  int nfrags = pInp->inp_sgen.n_regs;
  int ncaps = pInp->inp_sgen.n_caps;
  string code;
  Fragment* pfrag;
 // int count = 0;

    ofstream summary;
    summary.open("entire_//space.csv");
    summary << setprecision(6);
    
   for (int i=0; i < nfrags; i++)
 //    for (int i=0; i < 1; i++)
     {
       string istr = StringTools::int2str(i+1,1,"0");
 //        for (int j=0; j < 2; j++)
     for (int k=0; k < nfrags; k++)
     {
       string kstr = StringTools::int2str(k+1,1,"0");
       for (int j=0; j < ncaps; j++)
       {   
         string jstr = StringTools::int2str(j+1+nfrags, 1, "0");
         code = "0(" + istr + "-1(" + kstr + "-1(" + jstr + "-1)))";
//       cout << "********genecode: " << code << "  *********"<< endl;
         summary << i+1 << "-" << k+1 << "-" << j+nfrags+1 << ",";
         
         pfrag = pFrags[i];
         pfrag->print_props_csv(summary);
         pfrag = pFrags[k];
         pfrag->print_props_csv(summary);
         pfrag = pFrags[j+nfrags];
         pfrag->print_props_csv(summary);      
         summary << endl;           
       }
     }
    }
}
/*
void S_Genetic::get_fit4test(Node<GANode>* pNode)
{
   double MW=0;
   double tmp;
   int anumber;
   for (size_t i=0; i< pNode->data.natoms; i++)
   {
     anumber = PTable::atom_number(pNode->data.anames[i]);
     tmp = PTable::atom_mass(anumber);
     MW = MW + tmp;
   }

//     pNode->score = exp((MW-1000)/100);
     pNode->data.MW = MW;
//     pNode->data.score = MW;

   return;
}
*/
void S_Genetic::hashinit()
{
   int ngen = pInp->inp_ga.ga_stop + 1;
   int popsize = pInp->inp_ga.pop_size;

   int size = 2*(ngen*popsize) + 1;
   
   cout << "size of hash_all:" << size << endl; 
 
     genehash_all.resize(size);
  
   size = 2*popsize + 1;
 
   cout << "size of hash_pop:" << size << endl;

     genehash_pop.resize(size);

   hash_record.resize(2*NUM_ALL_PERM+1);
}

long S_Genetic::BKDRHash(string str)
{
   long seed = 131; // 31 131 1313 13131 131313 etc..  
   long hash = 0;  
   for(int i = 0; i < str.length(); i++)  
   {  
      hash = (hash * seed) + str.at(i);  
   }  
//   cout << "hash: " << hash << endl;
   return (hash & 0x7FFFFFFF);  
}
  
int S_Genetic::addtohash(vector<unique_ptr<LinkedList<Node<string> > > > & hashtable, string & str)
{

//cout << "hash_here" << hashtable.size() << endl;
  long size = hashtable.size();
  long key = BKDRHash(str) % size;
//  cout << key << endl;

  if (hashtable[key].get() == NULL)
  {
    unique_ptr<LinkedList<Node<string> > >  ptr(new LinkedList<Node<string> >());
    hashtable[key] = move(ptr);
    Node<string>* pNode = new Node<string>(str, NULL);
    hashtable[key]->addnode(pNode);
  }
  else
  {
     Node<string>* pNode = new Node<string>(str,NULL);
     hashtable[key]->addnode(pNode);
  }
//  cout << "Successfully add to hash table" << endl;
  return 1;
  
}
  
int S_Genetic::addtohash_record(vector<unique_ptr<LinkedList<Node<record> > > > & hash_record, Node<record>* pNode)
{
   long size = hash_record.size();
   long key = BKDRHash(pNode->data.gene) % size;

   if(hash_record[key].get() == NULL)
   {
   unique_ptr<LinkedList<Node<record> > >  ptr(new LinkedList<Node<record> >());
   hash_record[key] = move(ptr);
   hash_record[key]->addnode(pNode);
//   cout << "add new key" << endl;
   }
   else
   {
     hash_record[key]->addnode(pNode);
//     cout << "exist key" << endl;
   }
   return 1;
}

  
int S_Genetic::searchhash(vector<unique_ptr<LinkedList<Node<string> > > > & hashtable, string & genecode)
{
// cout << "search hash" << hashtable.size() << endl;
   long size = hashtable.size();
   long key = BKDRHash(genecode) % size;
   
//   cout << hashtable[key].get() << endl;
   if (hashtable[key].get() == NULL)  return 0;  //new entry
   else
   {
     if(hashtable[key]->findnode(genecode)) return 1; //not a new entry
     else return 0;   //new entry
   }
}  
  

void S_Genetic::searchhash_record(vector<unique_ptr<LinkedList<Node<record> > > > & hash_record, Node<GANode>* pNode)
{
// cout << "run to here" << endl;
   string genecode = pNode->data.gene;

   cout << genecode << endl;

   long size = hash_record.size();
   cout << size << endl;
   long key = BKDRHash(genecode) % size;
   cout << key << endl;

   Node<record>* ptr;

   cout << hash_record[key].get() << endl;
   
   ptr = hash_record[key]->getnode(genecode);

   cout << ptr << endl;

   pNode->data.HOMO = ptr->data.HOMO;
   pNode->data.LUMO = ptr->data.LUMO;
   pNode->data.MW = ptr->data.MW;
 
   cout << ptr->data.HOMO << " " << ptr->data.LUMO << endl;
     
   pNode->data.score_ETL = oled_fit(ptr->data.HOMO,ptr->data.LUMO)[0];
   pNode->data.score_HTL = oled_fit(ptr->data.HOMO,ptr->data.LUMO)[1];
   pNode->data.score = oled_fit(ptr->data.HOMO,ptr->data.LUMO)[2]; 

// cout << "score ETL =  " << pNode->data.score_ETL << "  score HTL = "<< pNode->data.score_HTL << endl;
   return;
}

void S_Genetic::initialHOSTref(Node<GANode>* pNode){
   
double diff;
double score_homo, score_lumo;
double thisscore;   

  
   diff = pNode->data.HOMO - (pInp->inp_sgen.HOST_HOMO);
   cout << "homo  " << pNode->data.HOMO<< endl; 
cout << "lumo  " << pNode->data.LUMO<< endl;
   score_homo = 1000 * exp(-(diff*diff)/0.125);
   diff = pNode->data.LUMO - (pInp->inp_sgen.HOST_LUMO);      /*hard-coded*/
   score_lumo = 1000 * exp(-(diff*diff)/0.125);  //2sigma = 0.5 (0.95 chance);
   thisscore=sqrt(score_lumo*score_homo);
   cout << "this score =  " <<thisscore << endl;
   if(thisscore>=refHOST_score) {
   refHOST_score=thisscore; 
   refHOST_HOMO=pNode->data.HOMO;
   refHOST_LUMO=pNode->data.LUMO;
}
 cout << "final score " << refHOST_score << endl;
}

void S_Genetic::oled_score(Node<GANode>* pNode, char OLEDtype)
{
   double score_homo;
   double score_lumo;
   double score;
   int type= pInp->inp_sgen.TransportType;
   double refHOMO;
   double refLUMO;     
   double diff;   
   
// cout << " Reference HOST updated " << refHOST_HOMO <<"/" << refHOST_LUMO <<endl;

switch(OLEDtype){
  case '1':
  cout << "calculating ETL score ...." << endl;
  refHOMO = refHOST_HOMO;
  refLUMO = refHOST_LUMO;
  score_homo = 1/(1+exp(-5.0* (refHOMO - pNode->data.HOMO)));
  diff = pNode->data.LUMO - (refLUMO);      /*hard-coded*/
  score_lumo = 1000 * exp(-(diff*diff)/0.125);  //2sigma = 0.5 (0.95 chance);
  score= score_lumo*score_homo;
  break;
//for ETL fitness function : refHOMO -8.60 refLUMO -0.667 // 

   case '2':
   cout << "calculating HTL score ...." << endl;
     refHOMO = refHOST_HOMO;
  refLUMO = refHOST_LUMO;
   diff = pNode->data.HOMO - (refHOMO);
   score_homo = 1000 * exp(-(diff*diff)/0.125);
   score_lumo = 1/(1+exp(5.0*((refLUMO)-pNode->data.LUMO)));
   score= score_lumo*score_homo;
   break;
}
//for HTL fitness function : refHOMO -8.02 refLUMO -0.85 //
   cout << score << endl;
   pNode->data.score = score;
   pNode->data.score_ETL = oled_fit(pNode->data.HOMO,pNode->data.LUMO)[0];
   pNode->data.score_HTL = oled_fit(pNode->data.HOMO,pNode->data.LUMO)[1];
    pNode->data.score = oled_fit(pNode->data.HOMO,pNode->data.LUMO)[2];

   return;
}

double * S_Genetic::oled_fit(double homo, double lumo)
{
// cout << "calculating score .... " << endl;
   double score_homo;
   double score_lumo;
   double *finalscore=new double[2];
   int type= pInp->inp_sgen.TransportType;  
   double refHOMO; 
   double refLUMO;
   double diff;


// cout << "ETL ref energy : " << pInp->inp_sgen.ETL_HOMO << "  " << pInp->inp_sgen.ETL_LUMO << endl;
// cout << "HTL ref energy : " << pInp->inp_sgen.HTL_HOMO << "  " << pInp->inp_sgen.HTL_LUMO << endl;
  if(homo==0&&lumo==0){
   finalscore[1]=0;
   finalscore[0]=0;
   finalscore[2]=0;
} 

else{ 

    refHOMO = refHOST_HOMO;
  refLUMO = refHOST_LUMO;
   score_homo = 1/(1+exp(-5.0* (refHOMO - homo)));
  diff = lumo - (refLUMO);      /*hard-coded*/
  score_lumo = 1000 * exp(-(diff*diff)/0.125);  //2sigma = 0.5 (0.95 chance);
  finalscore[0]=score_lumo*score_homo;

//for ETL fitness function : refHOMO -8.60 refLUMO -0.667 // 

   diff = homo - (refHOMO);
   score_homo = 1000 * exp(-(diff*diff)/0.125);
   score_lumo = 1/(1+exp(5.0*((refLUMO)-lumo)));
   finalscore[1]=score_lumo*score_homo;
  
   diff = homo - (refHOMO);
   score_homo = 1000 * exp(-(diff*diff)/0.125);
    diff = lumo - (refLUMO);      /*hard-coded*/ 
   score_lumo = 1000 * exp(-(diff*diff)/0.125);  //2sigma = 0.5 (0.95 chance);
    finalscore[2]=sqrt(score_lumo*score_homo);

  } return finalscore;
}

void S_Genetic::crossover(string & genecode1, string & genecode2, vector<data> tempC ){

   gTree gtree1(pInp);
   gTree gtree2(pInp);
  
   cout << genecode1 << endl;
   cout << genecode2 << endl;   

   gtree1.readcode(genecode1,tempC);
   gtree2.readcode(genecode2,tempC);
   
   vector<data> cArrayptr1, cArrayptr2,tempc1,tempc2;
   cArrayptr1=gtree1.getcArray();
   tempc1=cArrayptr1;
   cArrayptr2=gtree2.getcArray();
   tempc2=cArrayptr2; 
     int splice_layer = Utils::randomi2(gtree1.nlayers);
//   int splice_layer =1 ; 
//   cout << "Parent1 before cross over " << endl;

/* for (int i=0;i<cArrayptr1.size();i++){

         cout << cArrayptr1[i].label << "\t" << cArrayptr1[i].tag << "\t" << cArrayptr1[i].plug << "\t" << cArrayptr1[i].acceptor << "\t" << cArrayptr1[i].layer << "\t" << cArrayptr1[i].xyzw.size()<< endl;

   }
*/

//cout << "Parent2 before cross over " << endl;

/* for (int i=0;i<cArrayptr2.size();i++){

         cout << cArrayptr2[i].label << "\t" << cArrayptr2[i].tag << "\t" << cArrayptr2[i].plug << "\t" << cArrayptr2[i].acceptor << "\t" << cArrayptr2[i].layer << "\t" << cArrayptr2[i].xyzw.size()<< endl;

   }
*/

   cout << "determine splice point: " << splice_layer << endl;	


   for (int i=0; i<cArrayptr1.size();i++){
   // cout << "cArrayptr1.layer " << cArrayptr1[i].layer <<" splice " << splice_layer << endl;
      
     if(cArrayptr1[i].layer >= splice_layer ){
      
      cArrayptr1[i].label= tempc2[i].label;
      cArrayptr1[i].plug= tempc2[i].plug;
      
      cArrayptr2[i].label= tempc1[i].label;
      cArrayptr2[i].plug= tempc1[i].plug;

    }
  }
   
// cout << "Parent1 after cross over " << endl;
/*
   for (int i=0;i<cArrayptr1.size();i++){

       cout << cArrayptr1[i].label << "\t" << cArrayptr1[i].tag << "\t" << cArrayptr1[i].plug << "\t" << cArrayptr1[i].acceptor << "\t" << cArrayptr1[i].layer << "\t" << cArrayptr1[i].xyzw.size()<< endl;

   }

*/
//cout << "Parent2 after cross over " << endl;
/*
   for (int i=0;i<cArrayptr2.size();i++){

       cout << cArrayptr2[i].label << "\t" << cArrayptr2[i].tag << "\t" << cArrayptr2[i].plug << "\t" << cArrayptr2[i].acceptor << "\t" << cArrayptr2[i].layer << "\t" << cArrayptr2[i].xyzw.size()<< endl;

   }
*/
   newcode1 = gtree1.printCode(cArrayptr1);
   newcode2 = gtree2.printCode(cArrayptr2);
   
// ptr1->data.gene=newcode1;
// ptr2->data.gene=newcode2;

   cout << "parent1(crossover) : " << newcode1 << endl;
   cout << "parent2(crossover) : " << newcode2 << endl;

}


void S_Genetic::crossover(string & genecode1, string & genecode2)
{

  gTree gtree1(pInp);
  gTree gtree2(pInp);
 
  gtree1.readcode(genecode1);
  gtree2.readcode(genecode2);

  
  
  
  int count = 0;

  int splice_layer = Utils::randomi2(gtree1.nlayers);
 
  cout << "determine splice point: " << splice_layer << endl;
  vector<int> splice_1;
  vector<int> splice_2;

  for (int i=0; i< gtree1.nnodes; i++)
  {
    if (gtree1.cArray[i].layer == splice_layer)
    splice_1.push_back(i);
  }
  
  for (int i=0; i< gtree2.nnodes; i++)
  {
    if (gtree2.cArray[i].layer == splice_layer)
     splice_2.push_back(i);
  }

  int splice1 = Utils::randomi(splice_1.size()); 
//  cout << "splice 1: " << splice1 << endl;
  int splice2 = Utils::randomi(splice_2.size());
//  cout << "splice 2: " << splice2 << endl;

  int node1_idx = gtree1.cArray[splice_1.at(splice1)].c_index;
//  cout << "parent 1: " << parent1_idx << endl;
  int node2_idx = gtree2.cArray[splice_2.at(splice2)].c_index;

  TNode* parent1 = gtree1.getParentTree(node1_idx);
//  cout << "parent 1: "  << parent1 << " " << parent1->c_index << endl;

  TNode* parent2 = gtree2.getParentTree(node2_idx);

  if ((parent1->firstchild->c_index == splice_1.at(splice1)) && (parent2->firstchild->c_index == splice_2.at(splice2)))
  {
    swap_ptr(parent1->firstchild, parent2->firstchild);
    swap_ptr(parent1->firstchild->nextsibling, parent2->firstchild->nextsibling);
  }
  
  else if ((parent1->firstchild->c_index == splice_1.at(splice1)) && (parent2->nextsibling->c_index == splice_2.at(splice2)))
  {
     swap_ptr(parent1->firstchild, parent2->nextsibling);
     swap_ptr(parent1->firstchild->nextsibling, parent2->nextsibling->nextsibling);
  }
  else if ((parent1->nextsibling->c_index == splice_1.at(splice1)) && (parent2->firstchild->c_index == splice_2.at(splice2)))
  {
     swap_ptr(parent1->nextsibling, parent2->firstchild);
     swap_ptr(parent1->nextsibling->nextsibling, parent2->firstchild->nextsibling);
  }
  else if ((parent1->nextsibling->c_index == splice_1.at(splice1)) && (parent2->nextsibling->c_index == splice_2.at(splice2)))
  {
     swap_ptr(parent1->nextsibling, parent2->nextsibling);
     swap_ptr(parent1->nextsibling->nextsibling, parent2->nextsibling->nextsibling);
  }

  newcode1 = gtree1.tree_to_code(); 
  newcode2 = gtree2.tree_to_code();
}

void S_Genetic::swap_ptr(TNode* & ptr1, TNode* & ptr2)
{
   TNode* tmp;
   tmp = ptr1;
   ptr1 = ptr2;
   ptr2 = tmp;
 
   return;
}




string S_Genetic::mutation_asym(string genecode, vector<data> tempC){

   string newcode = genecode;
   gTree gtree(pInp);
   gtree.readcode(newcode,tempC);
   vector<data> cArray;
   cArray=gtree.getcArray();
   struct data seed;  
  
   int splice = Utils::randomi(gtree.nlayers); //splice = [0, nnodes-1];
   int layer = gtree.cArray[splice].layer;

   cout << "mutation layer: " << layer << endl;
  
   pTmp2 = ptr_array_ETL;   

   for (int i=0; i< cArray.size(); i++)
   {
    if(cArray[i].layer  >=layer){

       if (cArray[i].layer == 0){
   cArray[i].label = Utils::randomi2(pInp->inp_sgen.n_seeds);
   cArray[i].plug = gtree.readlib(seed, cArray[i].label);}
   else if (cArray[i].layer < gtree.nlayers){
   cArray[i].label = Utils::randomir(pInp->inp_sgen.n_seeds+1, pInp->inp_sgen.n_seeds +pInp->inp_sgen.n_regs);
   cArray[i].plug = gtree.readlib(seed, cArray[i].label);}
   else if (cArray[i].layer == gtree.nlayers){
   cArray[i].label = Utils::randomir(pInp->inp_sgen.n_seeds+pInp->inp_sgen.n_regs+1, pInp->inp_sgen.n_seeds + pInp->inp_sgen.n_regs + pInp->inp_sgen.n_caps); 
   cArray[i].plug = gtree.readlib(seed, cArray[i].label);
}
   }
}
   ptr_array_ETL = pTmp2;
   newcode = gtree.printCode(cArray);

   return newcode;

}

string S_Genetic::mutation(string genecode,vector<data> tempC)
{
   string newcode = genecode;
   
     int exist = 1;
     while(exist)
     {   
   gTree gtree(pInp);
   gtree.readcode(newcode,tempC);

//   int splice = Utils::randomi(gtree.nnodes); //splice = [0, nnodes-1];
   int splice = 2;
    cout << "check point 1 "<< splice << endl;

   TNode* ptr = gtree.getNodeTree(splice);
   
   gtree.destory(ptr->firstchild);
   ptr->firstchild = NULL;
 
   int layer = gtree.cArray[splice].layer;
   
   cout << "mutation layer: " << layer << endl;

   if (layer == 0)
   ptr->label = Utils::randomi2(pInp->inp_sgen.n_seeds); 
   else if (layer < gtree.nlayers)
   ptr->label = Utils::randomir(pInp->inp_sgen.n_seeds+1, pInp->inp_sgen.n_seeds +pInp->inp_sgen.n_regs);
   else if (layer == gtree.nlayers)
   ptr->label = Utils::randomir(pInp->inp_sgen.n_seeds+pInp->inp_sgen.n_regs+1, pInp->inp_sgen.n_seeds + pInp->inp_sgen.n_regs + pInp->inp_sgen.n_caps);
   
   int nbranches = gtree.readlib_lite(ptr->label);
   gtree.add_branch(ptr->firstchild, 1, nbranches, layer+1);

// cout << "nbranches : " << nbranches << endl; 

// gtree.add_branch(ptr->firstchild, 2, nbranches, layer+1);

   newcode = gtree.tree_to_code();
   
   cout << "mutated gene:" << newcode << endl;
 
     exist = searchhash(genehash_all, newcode); 

     if (exist) cout << "in the list already, mutate again" << endl;
     }
   cout << endl;

   return newcode;
}



void S_Genetic::gen_node_random(Node<GANode>* pNode, string mmopt_temp, string mopac_temp,LinkedList<Node<GANode> >* pNext )
{
  int exist = 1;
  ncat++;
  string genecode;
  vector <data> cArray;

    while(exist == 1)
    {
    gTree gtree(pInp);
  if(gtree.newStruct())
  genecode = gtree.printCode();
   cArray=gtree.getcArray();    
    exist = searchhash(genehash_all, genecode);
    }

  cout << "******unique entry ******" << endl;
    addtohash(genehash_all,genecode);
    addtohash(genehash_pop,genecode);

  pNode->data.gene = genecode;

  pNode->data.comment = StringTools::int2str(ncat, 4, "0");
  
  cout << "randomly create a molecule: " << genecode << endl;

     
    if(!search_from_databank(genecode)){cout << "A new gene generated "<< endl;}


#if 1
  string cname = mmopt_temp + "ff" + pNode->data.comment+".xyz";
  
  ofstream xyzfile(cname.c_str());
  gTree gtree2(pInp);
  gtree2.gene2struc2(pNode->data.gene, xyzfile,cArray);
  xyzfile.close();

   ifstream infile(cname.c_str());
   pNode->data.get_xyz(infile);

// cout << "cwd : " << mopac_folder << endl;
   Mopac::gen_inp(mopac_temp, pInp, pNode);

   infile.close();
#endif

  pNode->data.flag = 0;
  pNode->next = NULL;
  pNext->addnode(pNode);
  optr_list.push_back(ncat);
  
  return;
}

void S_Genetic::gen_node(Node<GANode>* pNode, vector<data> tempC, string mmopt_folder, string mopac_folder,LinkedList<Node<GANode> >* pNext)
 {
   ncat++;


    int old_label;
//   if (ncat % 1000 ==0){
     addtohash(genehash_all,pNode->data.gene);
     addtohash(genehash_pop, pNode->data.gene);


  pNode->data.comment = StringTools::int2str(ncat, 4, "0");
/*hard-coded*/
#if 0
   old_label = ncat + off_set; 

   string tmp = StringTools::int2str(old_label,4,"0");
#endif

#if 1
 string cname = mmopt_folder + "ff" + pNode->data.comment+".xyz";

ofstream xyzfile(cname.c_str());
gTree gtree(pInp);
gtree.gene2struc2(pNode->data.gene, xyzfile ,tempC);
xyzfile.close();

ifstream infile(cname.c_str());
pNode->data.get_xyz(infile);


#endif

pNode->data.flag = 0;
pNode->next = NULL;


pNext->addnode(pNode);
optr_list.push_back(ncat);

 Mopac::gen_inp(mopac_folder, pInp, pNode);


   cout << "Done gen_node()" << endl << endl;
//}
   return;
 }


void S_Genetic::read_databank(){

  string cwd = StringTools::getcwd_str();
  string fname = cwd + "/data";

  ifstream databank(fname.c_str(),ios::in);
  string line;
//vector<string> tok_
//
   getline(databank,line);
  
  while(getline(databank,line)){

  databanks.push_back(line); 
 
}
  databank.close();

//for(int i=0;i<databanks.size();i++){

// cout << "print data: " << databanks.at(i) << endl;
// }

  return;

}
    
vector<string> S_Genetic::dataSort(vector<string> datas){
    double datascore[datas.size()];
    int index[datas.size()];
    double scores[datas.size()];
    double temp;
    vector<string> reverseData;
    vector<string> sortData;    

    vector<string> tok_line;
    for(int i=0;i<datas.size();i++){
    tok_line = StringTools::tokenize(datas.at(i)," ");
    scores[i]=atof(tok_line[3].c_str());
    }
    vector<pair<double,int> >V;
    for(int i=0;i<datas.size();i++){
        pair<double,int>P=make_pair(scores[i],i);
        V.push_back(P);
    }
    
    sort(V.begin(),V.end());
    for(int i=0;i<datas.size();i++){
    string test=StringTools::tokenize(datas.at(V[i].second)," ")[0].c_str();
    string test1=StringTools::tokenize(datas.at(V[i].second)," ")[1].c_str();
    string test2=StringTools::tokenize(datas.at(V[i].second)," ")[2].c_str();
    string newString = string_trim(test)+"  "+string_trim(test1)+"  "+string_trim(test2)+"  "+to_string(V[i].first);
    sortData.push_back(newString);
}
   
   for(int i=datas.size()-1;i>=0;i--){
     reverseData.push_back(sortData.at(i));
  }
return  reverseData;
}

vector<string> S_Genetic::dataSort_EML(vector<string> datas){
    double datascore[datas.size()];
    int index[datas.size()];
    double scores[datas.size()];
    double temp;
    vector<string> reverseData;
    vector<string> sortData;
    vector<string> tok_line;
    for(int i=0;i<datas.size();i++){
    tok_line = StringTools::tokenize(datas.at(i)," ");
    scores[i]=atof(tok_line[5].c_str());
    }
    vector<pair<double,int> >V;
    for(int i=0;i<datas.size();i++){
        pair<double,int>P=make_pair(scores[i],i);
        V.push_back(P);
    }
    sort(V.begin(),V.end());
    for(int i=0;i<datas.size();i++){
    string test=StringTools::tokenize(datas.at(V[i].second)," ")[0].c_str();
    string test1=StringTools::tokenize(datas.at(V[i].second)," ")[1].c_str();
    string test2=StringTools::tokenize(datas.at(V[i].second)," ")[2].c_str();
    string test3=StringTools::tokenize(datas.at(V[i].second)," ")[3].c_str();
    string test4=StringTools::tokenize(datas.at(V[i].second)," ")[4].c_str();
    string newString = string_trim(test)+"  "+string_trim(test1)+"  "+string_trim(test2)+"  "+string_trim(test3)+"  "+string_trim(test4)+"  "+to_string(V[i].first);
    sortData.push_back(newString);
}
   for(int i=datas.size()-1;i>=0;i--){
     reverseData.push_back(sortData.at(i));
  }

    for(int i=0; i<3;i++){

// cout << "Top " << i+1 << " targets : " << reverseData.at(i)<< endl;
   }

return  reverseData;
}



vector<string> S_Genetic::dataSort_ETL(vector<string> datas){
    double datascore[datas.size()];
    int index[datas.size()];
    double scores[datas.size()];
    double temp;
    vector<string> reverseData;
    vector<string> sortData;
    vector<string> tok_line;
    for(int i=0;i<datas.size();i++){
    tok_line = StringTools::tokenize(datas.at(i)," ");
    scores[i]=atof(tok_line[3].c_str());
    }
    vector<pair<double,int> >V;
    for(int i=0;i<datas.size();i++){
        pair<double,int>P=make_pair(scores[i],i);
        V.push_back(P);
    }
    sort(V.begin(),V.end());
    for(int i=0;i<datas.size();i++){
    string test=StringTools::tokenize(datas.at(V[i].second)," ")[0].c_str();
    string test1=StringTools::tokenize(datas.at(V[i].second)," ")[1].c_str();
    string test2=StringTools::tokenize(datas.at(V[i].second)," ")[2].c_str();
    string test3=StringTools::tokenize(datas.at(V[i].second)," ")[4].c_str();
    string test4=StringTools::tokenize(datas.at(V[i].second)," ")[5].c_str();
    string newString = string_trim(test)+"  "+string_trim(test1)+"  "+string_trim(test2)+"  "+to_string(V[i].first)+"  "+string_trim(test3)+"  "+string_trim(test4);
    sortData.push_back(newString);
}
   for(int i=datas.size()-1;i>=0;i--){
     reverseData.push_back(sortData.at(i));
  }

    for(int i=0; i<3;i++){

// cout << "Top " << i+1 << " targets : " << reverseData.at(i)<< endl;
   }

return  reverseData;
}

vector<string> S_Genetic::dataSort_HTL(vector<string> datas){
    double datascore[datas.size()];
    int index[datas.size()];
    double scores[datas.size()];
    double temp;
    vector<string> reverseData;
    vector<string> sortData;
    vector<string> tok_line;
    for(int i=0;i<datas.size();i++){
    tok_line = StringTools::tokenize(datas.at(i)," ");
    scores[i]=atof(tok_line[4].c_str());
    }
    vector<pair<double,int> >V;
    for(int i=0;i<datas.size();i++){
        pair<double,int>P=make_pair(scores[i],i);
        V.push_back(P);
    }
    sort(V.begin(),V.end());
    for(int i=0;i<datas.size();i++){
    string test=StringTools::tokenize(datas.at(V[i].second)," ")[0].c_str();
    string test1=StringTools::tokenize(datas.at(V[i].second)," ")[1].c_str();
    string test2=StringTools::tokenize(datas.at(V[i].second)," ")[2].c_str();
    string test3=StringTools::tokenize(datas.at(V[i].second)," ")[3].c_str();
    string test4=StringTools::tokenize(datas.at(V[i].second)," ")[5].c_str();
    string newString = string_trim(test)+"  "+string_trim(test1)+"  "+string_trim(test2)+"  "+string_trim(test3)+"  "+to_string(V[i].first)+"  "+string_trim(test4);
    sortData.push_back(newString);
}
   for(int i=datas.size()-1;i>=0;i--){
     reverseData.push_back(sortData.at(i));
  }

  for(int i=0; i<3;i++){

// cout << "Top " << i+1 << " targets : " << reverseData.at(i)<< endl;
   }

return  reverseData;
}

string S_Genetic::string_trim(string in) {

    stringstream ss;
    string out;
    ss << in;
    ss >> out;
    return out;

}



string S_Genetic::dataParser(string line){
  
   string parseLine="";
   vector<string> tok_line;
   int count =0;
   double score=0;
   double score_ETL=0;
   double score_HTL=0;
   tok_line = StringTools::tokenize(line," ");
   
   refHOST_HOMO=pInp->inp_sgen.HOST_HOMO;
   refHOST_LUMO=pInp->inp_sgen.HOST_LUMO; 

  for(int i=0;i<tok_line.size();i++){
//   }
//     
   score_ETL=oled_fit(atof(tok_line[1].c_str()),atof(tok_line[2].c_str()))[0];
   score_HTL=oled_fit(atof(tok_line[1].c_str()),atof(tok_line[2].c_str()))[1];
   score=oled_fit(atof(tok_line[1].c_str()),atof(tok_line[2].c_str()))[2];
   parseLine=tok_line[0]+"    "+tok_line[1]+"    "+tok_line[2]+"    "+to_string(score_ETL)+"    "+to_string(score_HTL)+"    "+to_string(score);
}
// cout << "after parser : " << parseLine << endl;
   return parseLine;
}

bool S_Genetic::search_from_databank(string gene){

   bool isInDatabank=false;
   vector<string> tok_line;
   int count=0;
   for(int i=0;i<databanks.size();i++){
   tok_line = StringTools::tokenize(databanks.at(i)," ");
// cout << "Data check : " << tok_line[0] << "compared with "<< gene << endl; 
   if(tok_line[0]==gene){
    count ++;
    isInDatabank=true;
    cout << "Data exist ! " << endl;
   }
 }

  return isInDatabank;

}

string S_Genetic::parseGene(string genestring){

  vector<string> tok_line;
  string seed_num;
  string frag_num;

  tok_line  = StringTools::tokenize(genestring,"(");
  vector<string> geneparse;
  string temp;
  string check;

  for(int i=0;i<genestring.length();i++){
   if(genestring.at(i)!='('){
       temp=temp+genestring.at(i);
   }
   else{
    geneparse.push_back(temp);
    temp="";
   }
}
  geneparse.push_back(tok_line.at(tok_line.size()-1));
  
    seed_num = geneparse.at(0);
    for(int i=1;i<geneparse.size();i++){
   check=geneparse.at(i);
   if(check.at(check.length()-1)='-'){
    tok_line  = StringTools::tokenize(geneparse.at(i),"-");
    frag_num=frag_num+tok_line[0]+";";
   }
} 
   string total = seed_num+";"+frag_num;
   return total;
 }

vector<string> S_Genetic::screenData(vector<string> datas){
vector<string> datatemp;
vector<string> tok_line;
int count=0;
  
for(int i=0;i<datas.size();i++){
datatemp.push_back(datas.at(i));
for(int j=0;j<datatemp.size();j++){
if(datas.at(i)==datatemp.at(j)){
count++;
}
}
  if(count >1){
  datatemp.pop_back();
  cout << "i = " << i << " count = "<<count<< endl;   
}
count=0;
}
return datatemp;
}

void S_Genetic::calculateFragScore(vector<string> datas, char oledtype){

  vector <string> tok_line, tok_line2;
  int temp_number=0;
  int totallib = pInp->inp_sgen.n_seeds + pInp->inp_sgen.n_regs+ pInp->inp_sgen.n_caps;
  double fragScore[totallib],avgfragScore[totallib];
  int fragCount[totallib];  
//cout << "check" << endl;
//cout << pInp->inp_sgen.fragscorethreshold << endl;
  double threshold=pInp->inp_sgen.fscore;                                 
  bool isCount=false;

  for (int i=0;i<totallib;i++){
    fragScore[i]=0;
    avgfragScore[i]=0;
    fragCount[i]=0;
    }

  double totalscore_seed, totalscore_frag, totalscore_cap; 
   for (int i=0;i<datas.size();i++){
   
// cout << "string data " << datas.at(i) << endl;
   tok_line  = StringTools::tokenize(datas.at(i)," ");
   tok_line2 = StringTools::tokenize( parseGene(datas.at(i)) ,";");
   for(int i=0;i<tok_line2.size();i++){
    for(int j=0;j<totallib;j++){
   if((j+1)==atoi(tok_line2.at(i).c_str())){
   switch(oledtype){
   
   case '1':
   if(atof(tok_line[tok_line.size()-3].c_str())>threshold){
   fragScore[j]=fragScore[j]+atof(tok_line[tok_line.size()-3].c_str());
   fragCount[j]+=1;}
   oledtype='1';
   break;
   case '2':
   if(atof(tok_line[tok_line.size()-2].c_str())>threshold){
   fragScore[j]=fragScore[j]+atof(tok_line[tok_line.size()-2].c_str()); 
   fragCount[j]+=1;}
   oledtype='2';
   break;
   case '3':
   if(atof(tok_line[tok_line.size()-1].c_str())>threshold){
   fragScore[j]=fragScore[j]+atof(tok_line[tok_line.size()-2].c_str());
   fragCount[j]+=1;}
   oledtype='3';
   }

 }

}
  }
   }
   
totalscore_seed=0;
totalscore_frag=0;
totalscore_cap=0;


if(isCount){

for(int i=0;i<totallib;i++){
if(fragCount[i]>0){
fragScore[i]=fragScore[i]/fragCount[i];
}else{fragScore[i]=0; cout<<"come"<<endl;}
}}



  for(int g=0;g<totallib;g++){
   if(g<pInp->inp_sgen.n_seeds){
    totalscore_seed+=fragScore[g];
   }else if(g>=pInp->inp_sgen.n_seeds && g<pInp->inp_sgen.n_seeds+pInp->inp_sgen.n_regs){totalscore_frag+=fragScore[g];}
   else if(g>=pInp->inp_sgen.n_seeds+pInp->inp_sgen.n_regs){totalscore_cap+=fragScore[g];}
   }   

  for(int g=0;g<totallib;g++){
     if(g<pInp->inp_sgen.n_seeds){
         avgfragScore[g]=100*fragScore[g]/totalscore_seed;
           }else if(g>=pInp->inp_sgen.n_seeds && g<pInp->inp_sgen.n_seeds+pInp->inp_sgen.n_regs){avgfragScore[g]=100*fragScore[g]/totalscore_frag;}
               else if(g>=pInp->inp_sgen.n_seeds+pInp->inp_sgen.n_regs){avgfragScore[g]=100*fragScore[g]/totalscore_cap;}
                  }

   cout << "" << endl;

   switch(oledtype){

   case '1':
   cout << "Fragment's score prediction for Electron Transporting Layer" << endl;
   break;
   case '2':
   cout << "Fragment's score prediction for Hole Transporting Layer" << endl;
   break;
   case '3':
   cout << "Fragment's score prediction for Host Layer" << endl;
   break;
}

   cout << "  " << endl;
   for(int g=0;g<totallib;g++){
 
if(g<pInp->inp_sgen.n_seeds){
cout << "seed(" <<g+1<<") : "<< " count "<< fragCount[g]<<"  " << strFrag[g] <<"'s score : " << fragScore[g] << " avg "<<avgfragScore[g] << " per " << fragScore[g]/fragCount[g] <<endl;
}
else if(g>=pInp->inp_sgen.n_seeds && g<pInp->inp_sgen.n_seeds+pInp->inp_sgen.n_regs){
cout << "regs(" <<g+1<<") : "<< " count "<< fragCount[g]<<"  "<< strFrag[g] <<"'s score : " << fragScore[g] << " avg "<<avgfragScore[g] << " per " << fragScore[g]/fragCount[g] <<endl;}
else if(g>=pInp->inp_sgen.n_seeds+pInp->inp_sgen.n_regs){
cout << "caps(" <<g+1<<") : "<< " count "<< fragCount[g]<<"  "<< strFrag[g] <<"'s score : " << fragScore[g] << " avg "<<avgfragScore[g] << " per " << fragScore[g]/fragCount[g] <<endl;
}

}

}
void S_Genetic::updatedHOST(){

  string cwd = StringTools::getcwd_str();
  string fname = cwd + "/data";
  ifstream record0(fname.c_str(),ios::in);
  string total="";
  string gene="";
  string vecstr="";
  string line;
  vector<string> tok_line;
  vector<string> databanks2;
  double max;

while(getline(record0,line))
  {
  databanks.push_back(dataParser(line));
}
databanks=dataSort_EML(databanks);

  cout <<""<<endl;
  cout << "========Reference HOST updated=======" <<endl;
//cout  << databanks.at(0) << endl;
  tok_line=StringTools::tokenize(databanks.at(0)," ");
  refHOST_HOMO=atof(tok_line[1].c_str());
  refHOST_LUMO=atof(tok_line[2].c_str());
  total = to_string(refHOST_HOMO)+"/"+to_string(refHOST_LUMO);
  cout << "The Reference HOST(HOMO/LUMO) is set to " <<endl;
  cout << refHOST_HOMO << "/" << refHOST_LUMO << " based on current genes "<< endl;
  cout << "=====================================" << endl;
  cout << "  "<< endl;
  vecstr = tok_line[0] + "   " + total +"  " +tok_line[5].c_str();
  
  max = atof(tok_line[5].c_str());
  if(max >= bestsofarHOST) {
   bestsofarHOSTstr = vecstr;
   bestsofarHOST = max;
   topoled.push_back(vecstr);
  }else{
   topoled.push_back(vecstr);

  }
  thisoleds.push_back(vecstr); 
  record0.close();
  return;


}


void S_Genetic::read_record(){


  string cwd = StringTools::getcwd_str();
  string fname = cwd + "/data";
  ifstream record0(fname.c_str(),ios::in);
  
  string line;
  vector<string> tok_line;
  vector<string> databanks2;
 
//getline(record0,line);

  while(getline(record0,line))
  {
  databanks.push_back(dataParser(line));  
}
//databanks=screenData(databanks);
  databanks=dataSort_EML(databanks);
  cout << "Read best HOST from record ....  " << endl; 
//cout << databanks.at(0) << endl;
  tok_line=StringTools::tokenize(databanks.at(0)," ");
  refHOST_HOMO=atof(tok_line[1].c_str());
  refHOST_LUMO=atof(tok_line[2].c_str());
  cout << "The Reference HOST(HOMO/LUMO) is set to " << endl;
  cout << refHOST_HOMO << "/" << refHOST_LUMO << " based on HOST database "<< endl;
  databanks=dataSort_ETL(databanks);
  databanks2=dataSort_HTL(databanks);
  //screenData(databanks);
  if(pInp->inp_sgen.isWriteFragScore){
  calculateFragScore(databanks,'1');
  calculateFragScore(databanks,'2');
  calculateFragScore(databanks,'3'); 
  }
  for(int i=0;i<databanks.size();i++){
   Node<record>* pNode = new Node<record>();
   tok_line=StringTools::tokenize(databanks.at(i)," ");
   pNode->data.gene = tok_line[0];
   parseGene(tok_line[0]);
   pNode->data.HOMO = atof(tok_line[1].c_str());
   pNode->data.LUMO = atof(tok_line[2].c_str());
   
   
   addtohash_record(hash_record, pNode); 
  }
  cout << "Finish reading record file" << endl;
  cout << "******************************" << endl;
  cout << "*** ETL/HTL GA Search Start **" << endl;
  cout << "******************************"<< endl; 
  record0.close();
  return;
}
#if 0
double Genetic::getfit_Co(xyznode* pNode)
{
  int index;
  double alpha = 2.0;
  index = 484*(pNode->ntemp - 1) + 22*(pNode->gene[0] - 1) + (pNode->gene[1] - 1);
  pNode->binding = Eb[index];
  pNode->h_affinity = Ha[index];
  pNode->hydricity = Hyd[index];

  double ref_Eb = 22.015367;
  double ref_Ha = 338.684573;
  double ref_Hyd = 181.023427;
  double score;
//need to screen out some extremely weak binding first; 
  if (pNode->binding > 1000.0)
  {
    pNode->score = 0.0000;
  }
  else
  {
  pNode->score = exp((ref_Ha - pNode->h_affinity)/alpha);
//  pNode->score = exp((ref_Ha - pNode->h_affinity)/alpha) + exp((ref_Hyd - pNode->hydricity)/50.0);
//   double diff = pNode->binding - 20.;
//  pNode->score = exp(-(diff*diff)/50.0);  //standard dev = 5;
  }
  score = pNode->score;

  return score;
}
#endif


void S_Genetic::spawn(string mmopt_temp,string mopac_temp,LinkedList<Node<GANode> >* pNtrp, char oledtype)
{
   int pool = pInp->inp_ga.pop_size;
   int npairs = pool/2;
   int nextra = pool%2;    // if odd pop_size
   double random;
   bool ptrs_identic;
   vector<data> tempC;
   gTree gtree3(pInp);
   gtree3.newStruct();
   tempC=gtree3.getcArray();
   bool isBestRestart=pInp->inp_sgen.bestRestart;
   double ratio=1;

//      for (int m=0;m<pool;m++){
//   ptr_array_ETL[m]->data.gene="";
 //  }

   cout << "is Best Restart mode on ? " << isBestRestart << endl;
   vector<int> tmp;
   vector<int>::iterator it;
   if (ngen == 0)
   {
    #if 1
     if(isBestRestart){
   cout << "The gene of 1st round will be generated from exist record" << endl;   
     for (int i=0; i<pool;i++){
       Node<GANode>* pNode0 = new Node<GANode>();
       vector<string> tok_line;
       switch(oledtype){
       case '1':
       tok_line = StringTools::tokenize(dataSort_ETL(databanks).at(rand()%(int)databanks.size()*0.005)," ");
       break;
       case '2':
       tok_line = StringTools::tokenize(dataSort_HTL(databanks).at(rand()%(int)databanks.size()*0.005)," ");
       break;
       }
       pNode0->data.gene = tok_line[0];
       gen_node(pNode0,tempC,mmopt_temp,mopac_temp,pNtrp);
       cout << "Generate Gene from Database : " << pNode0->data.gene << endl;
       }
   }
  else{
     cout << "The gene of 1st round will be generated randomly" << endl;  
     for (int i=0; i< pool; i++)
     {
       Node<GANode>* pNode = new Node<GANode>();
       gen_node_random(pNode,mmopt_temp,mopac_temp,pNtrp);     
    }
   
  }
    #endif
/*enhanced-initiation*/
   #if 0
    enhanced_init_gen();
   #endif 
   ncat =0;
  }
   else
   {
     
     pNtrp->clear();
     cout << "Current length of pNext List: " << pNtrp->listlength << endl;
//   gen_pairs(pool);
     for (int i=0; i < npairs; i++)
     {
     switch(oledtype){
     case '1':
     cout << oledtype << endl;
     ptrs_identic = get_ptrs_ETL(i);
     cout <<"ptr1 and ptr2  => " << ptr1string << "  " << ptr2string << endl;
     break;
     case '2':
     ptrs_identic = get_ptrs_HTL(i);
     cout <<"ptr1 and ptr2  => " << ptr1string << "  " << ptr2string << endl;
     break;
     }     
     cout << "Applying genetic operator to " << i+1 << "-th pair " << endl;
     cout << "parent1 : " << ptr1string << endl;
     cout << "parent2 : " << ptr2string << endl;
     cout << endl;
    
     Node<GANode>* pNode1 = new Node<GANode>();
     Node<GANode>* pNode2 = new Node<GANode>();
/*
     if(ptrs_identic)
     {
       cout << "***parent1 and parent2 are identical!***" << endl;
       cout << "test => " << ptr1->data.gene << "  " << ptr2->data.gene << endl;
       int exist = searchhash(genehash_pop, ptr1->data.gene);
      if (!exist)
       {
        cout << "They are not in the current population" << endl;
        pNode1->data.gene = mutation_asym(ptr1->data.gene,tempC);
        cout << "new mutation ptr1 : " << pNode1->data.gene << endl;
        gen_node(pNode1,tempC,mmopt_temp,mopac_temp,pNtrp);
        
        pNode2->data.gene = mutation_asym(pNode1->data.gene,tempC);
        ptr2=pNode2;   
        cout << "new mutation ptr2 : " << pNode2->data.gene << endl;
        gen_node(pNode2,tempC,mmopt_temp,mopac_temp,pNtrp);

       }

       else 
       {
         cout << "They are already in the Current population" << endl;
         pNode1->data.gene = mutation_asym(ptr1->data.gene,tempC);
         ptr1=pNode1;
         cout << "new mutation ptr1 : " << pNode1->data.gene << endl;
         gen_node(pNode1,tempC,mmopt_temp,mopac_temp,pNtrp);
         pNode2->data.gene = mutation_asym(ptr2->data.gene,tempC);
         ptr2=pNode2;
         cout << "new mutation ptr1 : " << pNode2->data.gene << endl;
         gen_node(pNode2,tempC,mmopt_temp,mopac_temp,pNtrp);
        
       }       
     }
     else
     {
   
*/
     random = Utils::randomf(1.0);     
     cout << "Get random number: " << random << "; Compared with crossover rate " << pInp->inp_ga.cross_r << endl;
 

     if (random < pInp->inp_ga.cross_r)
     {
      cout << "********Crossover and Mutation:********** " << endl;
      
      bool done = false;
      int count = 0;
      crossover(ptr1string,ptr2string,tempC);   
      if ((!search_from_databank(newcode1)) && (!search_from_databank(newcode2))) 
      {
       cout << "find unique offsprings" << endl; 
       done = true; 
    }

      if (!done)  
      {
        newcode1 = ptr1string;
        newcode2 = ptr2string;

      }
      double random;
      random = Utils::randomf(1.0);
      if (random < pInp->inp_ga.muta_r)   
      {
        cout << "ptr1 is under mutation" << endl;
        newcode1 = mutation_asym(newcode1,tempC);    
     	cout << "new mutation ptr1 : " << newcode1 << endl;
     	search_from_databank(newcode1);

    }  
      random = Utils::randomf(1.0);
      if (random < pInp->inp_ga.muta_r)   
      {
        cout << "ptr2 is under mutation" << endl;
        newcode2 = mutation_asym(newcode2,tempC);
        cout << "new mutation ptr2 : " << newcode2 << endl;
        search_from_databank(newcode2);
      }

      if (search_from_databank(newcode1))  
      {
      cout << "newcode 1 is already in the list" << endl;
      newcode1 = mutation_asym(newcode1,tempC);
      cout << "new mutation ptr1 : " << newcode1 << endl;
      search_from_databank(newcode1);
      }
      if (search_from_databank(newcode2))  
      {   
       cout << "newcode 2 is already in the list" << endl;
       newcode2 = mutation_asym(newcode2,tempC);
       cout << "new mutation ptr2 : " << newcode2 << endl;
       search_from_databank(newcode2);
      }
       pNode1->data.gene = newcode1;
       pNode2->data.gene = newcode2;
       cout << "save new molecules to the list" << endl; 
       gen_node(pNode1,tempC,mmopt_temp,mopac_temp,pNtrp);
       gen_node(pNode2,tempC,mmopt_temp,mopac_temp,pNtrp);
  }  //end if
    else 
     {
       cout << "*********Clone Parent********" << endl;
       double random;
       random = Utils::randomf(1.0);
       if (random < pInp->inp_ga.muta_r)  
       {
       cout << "prt1 is under mutation" << endl;
       newcode1 = mutation_asym(ptr1string,tempC);
        pNode1->data.gene = newcode1;   
       cout << "new mutation ptr1 : " << newcode1 << endl; 
       gen_node(pNode1,tempC,mmopt_temp,mopac_temp,pNtrp);
       }
       else
       {
 //        pNode1 = cpy_node(ptr1);
 //        pNtrp->addnode(pNode1);   //clone one of them
      pNode1->data.gene = ptr1string;
      gen_node(pNode1,tempC,mmopt_temp,mopac_temp,pNtrp);
    }

       random = Utils::randomf(1.0);
       if (random < pInp->inp_ga.muta_r)  
       {
         cout << "prt2 is under mutation" << endl;
         newcode2 = mutation_asym(ptr2string, tempC);
         pNode2->data.gene = newcode2;
         cout << "new mutation ptr2 : " << newcode2 << endl;
         gen_node(pNode2,tempC,mmopt_temp,mopac_temp,pNtrp);
       }
       else
       {
 //tt   if (search_from_databank(ptr2->data.gene))
       if(search_from_databank(ptr2string))
   {
//      pNode2 = cpy_node(ptr2);
//      pNtrp->addnode(pNode2);   //clone one of them
       pNode2->data.gene = ptr2string;
      gen_node(pNode2,tempC,mmopt_temp,mopac_temp,pNtrp);
      }
        else
        {
        cout << "prt2 is already in the Current population" << endl;
        pNode2->data.gene = mutation_asym(ptr2string,tempC); 
       gen_node(pNode2,tempC,mmopt_temp,mopac_temp,pNtrp);
        }
       }
     } //else for ****Clone****
 //   }  // end of (prts_identic) else
    
    }   //end of loop
     if (nextra)
     { 
      cout << "nextra" << endl;
      Node<GANode>* pOdd = new Node<GANode>();
      ncat++;
      int flag;
      Node<GANode>* pointer;
//    cout << "oledtypeall " << oledtypeall << endl;
      if(oledtype='1'){
      permutation_ETL[pool-1];
      pointer = ptr_array_ETL[select_list_ETL[flag]];}
      else{
      flag=permutation_HTL[pool-1];
      pointer = ptr_array_HTL[select_list_HTL[flag]];}
      pOdd->data.gene = mutation(pointer->data.gene,tempC);
     }
//
   } //end of else
}

//stochastic acceptance O(1) version

void S_Genetic::select_parents()
{
/*
  int pool = pInp->inp_ga.pop_size;
  int i = 0;

  while(i < pool)
  {
//    int random = Utils::randomi(pool);
    double random2 = Utils::randomf(1.0);
    int index = int(random2 * pool);
    if (random2 < prob_array[index]/max_prob)
       {
       select_list[i] = index;              //accepted
       i++;
       }  
  }
*/
    if (pInp->inp_ga.RWS == true)
    {
    cout << "run RWS" << endl;
    select_RWS();
    }
    else if (pInp->inp_ga.SUS == true)
    {    
    select_RWS_SUS();
    }
    else if (pInp->inp_ga.B_Tournament == true)
    select_binary_tournament();
    else if (pInp->inp_ga.L_Rank == true)
    select_linear_rank_based(pInp->inp_ga.L_Rank_rate);
    else if (pInp->inp_ga.E_Rank == true)
    {
    cout << "run E Rank" << endl;
    select_exp_rank_based(pInp->inp_ga.E_Rank_base);
    }

    cout << "Done selection operation" << endl;
        
}

void S_Genetic::select_RWS()    //roulette Wheel selection
{
    int pool = pInp->inp_ga.pop_size;

    vector<double> wheel(pool+1);
    int i = 0;
    wheel[0] = 0.0;

    while(i < pool)
    {
    //  cout << "oledtypeall_RWS " << oledtypeall << endl;
        switch(oledtypeall){
        case '1':
        wheel[i+1] = prob_array_ETL[i] + wheel[i];
        break;
        case '2':
        wheel[i+1] = prob_array_HTL[i] + wheel[i];
        break;
        }
        cout << "wheel " << i+1 << " " << wheel[i+1] << endl;
        i++;
    }
    
    i = 0;

    while(i<pool)
    {
       double random = Utils::randomf(wheel.back());
       cout << random << endl;
    // cout << "oledtypeall_RWS2 " << oledtypeall << endl;
       switch(oledtypeall){
       case '1':
       select_list_ETL[i] = select_binary_search(wheel, random);
//     cout << "select " << select_list_ETL[i] << endl;
       break;
       case '2':
       select_list_HTL[i] = select_binary_search(wheel, random);
//     cout << "select " << select_list_HTL[i] << endl;
       break; 
       }
       i++;          
    }
}

int S_Genetic::select_binary_search(vector<double> & wheel, double random)
{
    int left = 0;
    int right = wheel.size()-1;
    int mid;
    
    while (left <= right)
    {
        mid = left + (right-left)/2;
        if (wheel[mid] == random) return mid;
        else if (wheel[mid] < random)
        left = mid+1;
        else 
        right = mid-1; 
    }

    return (left-1);
}

void S_Genetic::select_RWS_SUS()   //roulette Wheel + stocastic universal sampling
{
   int  pool = pInp->inp_ga.pop_size;
   
   vector<double> wheel(pool+1);
   int i = 0;
   wheel[0] = 0.0;
// cout << "oledtypeall (selectRWSSUS) " << oledtypeall << endl;
   while(i < pool)
   {
      switch(oledtypeall){
      case '1':
      wheel[i+1] = prob_array_ETL[i] + wheel[i];
      break;
      case '2':
      wheel[i+1] = prob_array_HTL[i] + wheel[i];
      break;
      i++;
   }}
 
   i = 0;

   double random = Utils::randomf(wheel.back());
   double stepsize = wheel.back()/pool;

   while (i < pool)
   {
//    cout << "oledtypeall (RWSSUS2) " << oledtypeall << endl;
      switch(oledtypeall){
      case '1':
      select_list_ETL[i] = select_binary_search(wheel, random);
      break;
      case '2':
      select_list_HTL[i] = select_binary_search(wheel, random);
      break;
      } random += stepsize;
      if (random > wheel.back())
            random = fmod(random, wheel.back());
      i++;         
   }        
}

void S_Genetic::select_linear_rank_based(double rate) //the reproduction rate of worst individual
{
   int pool = pInp->inp_ga.pop_size;
   vector<double> wheel(pool+1);
   int i = 0;
   wheel[0] = 0.0;
   while(i < pool)
   {
      wheel[i+1]= (rate + 2 * (1 - rate) * (pool-i-1)/(pool-1))/double(pool) + wheel[i]; 
      cout << "wheel " << i+1 << " " << wheel[i+1] << endl;
      i++;      
   }
 
   i = 0;

   double stepsize = wheel.back()/pool;
   double random = Utils::randomf(wheel.back());
//cout << "oledtypeall (select LE) " << oledtypeall << endl;
   while(i< pool)
   {
     cout << random << endl;
//   cout << "oledtypeall " << oledtypeall << endl;
     switch(oledtypeall){
     case '1':
     select_list_ETL[i] = select_binary_search(wheel,random);
//   cout << "select " << select_list_ETL[i] << endl;
     break;
     case '2':
     select_list_HTL[i] = select_binary_search(wheel,random);
//   cout << "select " << select_list_HTL[i] << endl;
     break;}   
     random += stepsize;
     if (random > wheel.back())
        random = fmod(random,wheel.back());
     i++;
   }    
}

void S_Genetic::select_exp_rank_based(double base) // 0 < base < 1
{
   int pool = pInp->inp_ga.pop_size;
   vector<double> wheel(pool+1);
   int i = 0;
   wheel[0] = 0.0;
   while(i < pool)
   {
     wheel[i+1] = (base-1)/(pow(base,double(pool))-1)*pow(base,double(i)) + wheel[i];
     cout << "wheel " << i+1 << " " << wheel[i+1] << endl;
     i++;
   } 

   i = 0;

// cout << "oledtypeall (select ER) " << oledtypeall << endl;
   double stepsize = wheel.back()/pool;
   double random = Utils::randomf(wheel.back());
   while ( i < pool)
   {
     cout << random << endl;
//  cout << "oledtypeall(ER) " << oledtypeall << endl; 
     switch(oledtypeall){
     case '1':
     select_list_ETL[i] = select_binary_search(wheel,random);
//   cout << "select " << select_list_ETL[i] << endl;
     break;
     case '2': 
     select_list_HTL[i] = select_binary_search(wheel,random);
//   cout << "select " << select_list_HTL[i] << endl;
     break;   
     }
     random += stepsize;
     if (random > wheel.back())
        random = fmod(random, wheel.back());
     i++;
   }
}

void S_Genetic::select_binary_tournament()
{
    int pool = pInp->inp_ga.pop_size;
    int i = 0;
    
    while (i < pool)
   {
      double random = Utils::randomf(1.0);
      int index1 = int(floor(pool * random + 0.5));  //rounding half up only work for positive values
      random = Utils::randomf(1.0);
      int index2 = int(floor(pool * random + 0.5));
      while(index1 == index2)
        {
           random = Utils::randomf(1.0);
           index2 = int(floor(pool * random + 0.5));
        }         
//   cout << "oledtypeall (select BR) " << oledtypeall << endl;
      switch(oledtypeall){
      case '1':
      select_list_ETL[i] = (prob_array_ETL[index1]>prob_array_ETL[index2])?index1:index2;
      break;
      case '2':
      select_list_HTL[i] = (prob_array_HTL[index1]>prob_array_HTL[index2])?index1:index2;
      break;
      }i++;   
   }    
}

void S_Genetic::gen_pairs(int pool)
{
   int flag;
   for (int i= pool; i>1; i--){
     cout << "permutation : " << permutation_ETL[i] << endl;
     flag = Utils::randomi(i);
     swap(permutation_ETL[flag],permutation_ETL[i-1]);        //shuffle     
   }

   cout << "Selected parents:" << endl;
   for (int i=0; i< pool; i++){
     cout << select_list_ETL[permutation_ETL[i]] << endl;
   }
   cout << "******************" << endl;
}


//int n: the number of n-pairs
//save pointers for this pair to ptr1 & ptr2
//if ptr1 & ptrs2 are identical return true, otherwise return false
// cout << "parent 2: " << ptr2 << endl;

bool S_Genetic::get_ptrs_ETL(int n)
{
   int flag;
   int index;
   
// cout << "oledtypeall_getETL " << oledtypeall << endl;
   for(int i=0;i<pInp->inp_ga.pop_size;i++){
// cout <<"select : " << i <<"  "<<select_list_ETL[i]<< endl;}
    }
// flag = permutation_ETL[2*n+0];
// index = select_list_ETL[flag];
   ptr1string=tempstring_ETL.at(select_list_ETL[2*n+0]);

// flag = permutation_ETL[2*n+1];
// index = select_list_ETL[flag];
   ptr2string=tempstring_ETL.at(select_list_ETL[2*n+1]);

    if(ptr1string == ptr2string)
   return true;
   else
   return false;
}
bool S_Genetic::get_ptrs_HTL(int n)
{
// cout << "oledtypeall_getHTL " << oledtypeall << endl;
   int flag;
   int index;
   for(int i=0;i<pInp->inp_ga.pop_size;i++){
// cout <<"select : " << i <<"  "<<select_list_HTL[i]<< endl;}
   }
//   flag = permutation_HTL[2*n+0];
//   index = select_list_HTL[flag];
   ptr1string = tempstring_HTL.at(select_list_HTL[2*n+0]);
//   flag = permutation_HTL[2*n+1];
//   index = select_list_HTL[flag];
   ptr2string = tempstring_HTL.at(select_list_HTL[2*n+1]);
   if(ptr1string == ptr2string)
   return true;
   else
   return false;
}


#if 0
bool Genetic::mutate(vector<int> & offspring)
{
   double random;
   int bits = int(offspring.size());
   bool success=false;

   for (int i=0; i< bits; i++)
   {
     if (i == 0)
     {
       random = Utils::randomf(1.0);
       if (random < pInp->inp_ga.muta_r)
       {
         offspring[i] = Utils::randomi2(pInp->n_temp);
         success = true;
       }
     }
     else
     {
        random = Utils::randomf(1.0);
        if (random < pInp->inp_ga.muta_r)
        {
          offspring[i] = Utils::randomi2(pInp->n_rfiles);
          success = true;
        }
      }
   }
   return success;
}

bool Genetic::mutate2(xyznode* prt, xyznode* pKid)
{
  vector<int> string;
  vector<int>::iterator it;

  string.push_back(prt->ntemp);
  for (it = prt->gene.begin(); it != prt->gene.end(); ++it)
  {
    string.push_back(*it);
  } 

  if(mutate(string) && check_unique2(genelist,string))
  {
    cout << "Going through a mutation: " << endl;
    pKid->ntemp = string[0];
    vector<int> tmp;
    cout << "The new gene is: " << pKid->ntemp << " ";
    for (it = string.begin()+1; it != string.end(); ++it)
    {
       tmp.push_back(*it);
       cout << *it << " ";
    }
    cout << endl;
    pKid->gene = tmp;
    
    return true;
  }
  else
  return false;
}


void Genetic::force_mutate(xyznode* prt, xyznode* pKid)
{
   vector<int> string;
   vector<int>::iterator it;
   int random;
   bool done = false;
  
   string.push_back(prt->ntemp);
   for (it = prt->gene.begin(); it != prt->gene.end(); ++it)
   {
     string.push_back(*it);
   }

   vector<int> offspring;

   int bits = int(string.size());

   cout << "Force mutation to occur: " << endl;

   while (!done)
   {
     random = Utils::randomi(bits);
     offspring = string;
     if (random == 0) 
     {
       offspring[random] = Utils::randomi2(pInp->n_temp);
       done = check_unique2(genelist, offspring);
     }      
     else
     {
       offspring[random] = Utils::randomi2(pInp->n_rfiles);
       done = check_unique2(genelist, offspring);
     }
   }
   
//   cout << "Force mutation to occur: " << endl;
   pKid->ntemp = offspring[0];
   vector<int> tmp;
   cout << "The new gene is: " << pKid->ntemp << " ";
   for (it = offspring.begin()+1; it != offspring.end(); ++it)
   {
     tmp.push_back(*it);
     cout << *it << " ";
   }
   cout << endl;
   pKid->gene = tmp;

   return;
}

//int n: the number of n-pairs
bool Genetic::splice_mutate()
{
   vector<int> gene;
   vector<int>::iterator it;
   vector<int> offspring1;
   vector<int> offspring2;

//     get_prts(n);
     
     gene.push_back(prt1->ntemp);
     for (it = prt1->gene.begin(); it != prt1->gene.end(); ++it)
     {
        gene.push_back(*it);
     }
     int gene_bits = int(gene.size());

     parents[0] = gene; 

     cout << "gene for parent1: ";
     for (it = gene.begin(); it != gene.end(); ++it)
     cout << *it << " " ;
     cout << endl;
    
     gene.clear();
  
     gene.push_back(prt2->ntemp);
     for (it = prt2->gene.begin(); it != prt2->gene.end(); ++it)
     {
        gene.push_back(*it);
     }
     
     parents[1] = gene;

     cout << "gene for parent2: ";
     for (it = gene.begin(); it != gene.end(); ++it)
     cout << *it << " " ;
     cout << endl;
    
     int count = 1;
     bool unique1;
     bool unique2;
     bool success = false;

     while (count < gene_bits)
     {
     count++;
     vector<int>().swap(offspring1);
     vector<int>().swap(offspring2);

//e.g. splice_point = 2, between bit 1 & 2 (gene starting from bit 0)
     int splice_point = Utils::randomi2(gene_bits-1);
     cout << "Determine splice_point as " << splice_point << endl;
     
     for (int i=0; i != splice_point; i++)
     {
     offspring1.push_back(parents[0][i]);
     offspring2.push_back(parents[1][i]);
     }

     for (int i= splice_point; i < gene_bits; i++)
     {
     offspring1.push_back(parents[1][i]);
     offspring2.push_back(parents[0][i]);
     }

   cout << "gene for offspring1: ";
   for (it = offspring1.begin(); it != offspring1.end(); ++it)
   cout << *it << " " ;
   cout << endl;
  
   cout << "gene for offspring2: ";
   for (it = offspring2.begin(); it != offspring2.end(); ++it)
   cout << *it << " " ;
   cout << endl;

   if(mutate(offspring1))
   cout << "offspring 1 went through a mutation:" << endl;
   cout << "The new gene is: ";
   for (it = offspring1.begin(); it != offspring1.end(); ++it)
   cout << *it << " " ;
   cout << endl<<endl;

   if(mutate(offspring2))
   cout << "offspring 2 went through a mutation:" << endl;
   cout << "The new gene is: ";
   for (it = offspring2.begin(); it != offspring2.end(); ++it)
   cout << *it << " " ;
   cout << endl<<endl;

   unique1 = check_unique2(genelist,offspring1);
   unique2 = check_unique2(genelist,offspring2);

   if (unique1 && unique2)
   {
   success = true;
   cout << "They are unique combinations" << endl;
   offsprings[0] = offspring1;
   offsprings[1] = offspring2;
   break;
   }
   else
   cout << "Generated already, Try again!" << endl;
  } 

     return success;
}   
#endif

S_Genetic::S_Genetic(readinp *p1)
{
  pInp = p1;  
  parents.resize(2);
  offsprings.resize(2);
  ptr_array_ETL.resize(pInp->inp_ga.pop_size);
  ptr_array_HTL.resize(pInp->inp_ga.pop_size);
  prob_array_ETL.resize(pInp->inp_ga.pop_size);
  prob_array_HTL.resize(pInp->inp_ga.pop_size);
  select_list_ETL.resize(pInp->inp_ga.pop_size);
  select_list_HTL.resize(pInp->inp_ga.pop_size);
  permutation_ETL.resize(pInp->inp_ga.pop_size); 
  permutation_HTL.resize(pInp->inp_ga.pop_size);
  hashinit();
}

S_Genetic::~S_Genetic()
{
    vector<Fragment*>::iterator it;
    for (it = pFrags.begin(); it != pFrags.end(); it++)
    {
       delete *it;
    }
}
#if 0
int Genetic::Tempcode()
{
  int nt = pInp->n_temp;
  int ntemp = Utils::randomi2(nt);
  return ntemp;
}

int Genetic::Rcode()
{
  int nrs = pInp->n_rfiles;
  int rcode = Utils::randomi2(nrs);
  return rcode;
}

bool Genetic::check_unique2(vector<vector<int> > & target_list, vector<int> & offspring)
{
  size_t ncats = target_list.size();
  size_t gene_bits = offspring.size();
  int status=1;

  if (ncats == 0)
  {
   return true;
  } 
 
  for (size_t i=0; i < ncats; i++)
  {
    if (target_list[i][0] == offspring[0])
    {
      status = 0;
      for (size_t j=1; (j<gene_bits) && (status == 0); j++)
      {
        if (target_list[i][j] == offspring[j])
        status = 0;
        else 
        {
          status = 1;
        }
      }
    } 
      if (status == 0)
      {
      cout << "This is not unique combination" << endl;
      break;
      }
  }

  if (status == 0)
  {
//    vector<int>().swap(offspring);
    return false;
  }
  else
    return true;
}

bool Genetic::check_unique(vector<vector<int> > & target_list, int tempcode, vector<int>& rcodes)
{
  size_t ncats = target_list.size();
  size_t nvars = rcodes.size();
  int status = 1;

  if (ncats == 0)
  {
    return true;
  }
  
  for (size_t i=0; i < ncats; i++)
  {
    if (target_list[i][0] == tempcode)
    {  
      status = 0;     
      for (size_t j=0; (j < nvars) && (status == 0) ; j++)
      {
        if (target_list[i][j+1] == rcodes[j])
        status = 0;
        else 
        {
         status = 1;
        }
      }
     }
      if (status == 0)
      break; 
   }

   if (status == 0) 
     { 
//      vector<int>().swap(rcodes);
      return false;
     }
   else 
      return true;
}

void Genetic::save_genelist(int tempcode, vector<int> rcodes)
{
   vector<int>::iterator it;

   vector<int> tmp;

   tmp.push_back(tempcode);
   for (it = rcodes.begin(); it != rcodes.end(); it++)
   {
     tmp.push_back(*it);
   }
   
   genelist.push_back(tmp); 
}

void Genetic::save_gene_in_pop(int tempcode, vector<int> rcodes)
{
   vector<int>::iterator it;
   vector<int> tmp;

   tmp.push_back(tempcode);
   
   for (it = rcodes.begin(); it != rcodes.end(); it++)
   {
     tmp.push_back(*it);
   }
  
   gene_in_pop.push_back(tmp); 
}

void Genetic::gencode_random(xyznode* pNode)
{

   bool unique = false;
   int tempcode;
   int nvars = pInp->inp_ga.num_r;
   vector<int> rcodes;
   
 
   while(unique == false)

   {
     tempcode = Tempcode();
     cout << tempcode << " ";
     vector<int>().swap(rcodes); 
     
     for (int i=0; i < nvars; i++)
     {
       int r = Rcode();
       cout << r << " ";
       rcodes.push_back(r);
     }   

     cout << endl;

     unique = check_unique(genelist,tempcode, rcodes);

   }
       cout << "This is a unique combination" << endl;
//       save_genelist(tempcode, rcodes);
       pNode->ntemp = tempcode;
       pNode->gene = rcodes;
}
#endif
/*
//for enhanced initiation and enhanced mutation
void Genetic::all_candidates()
{
  int nfrags = pInp->inp_sgen.n_regs;
  int ncaps = pInp->inp_sgen.n_caps;
  string code;
  vector<gTree*> all_genes_tree;
  int MaxHeight = 0;
 // int count = 0;

   for (int i=0; i < nfrags; i++)
 //    for (int i=0; i < 1; i++)
     {
       string istr = StringTools::int2str(i+1,1,"0");
 //        for (int j=0; j < 2; j++)
     for (int k=0; k < nfrags; k++)
     {
       string kstr = StringTools::int2str(k+1,1,"0");
       for (int j=0; j < ncaps; j++)
       {
         string jstr = StringTools::int2str(j+1+nfrags, 1, "0");
         code = "0(" + istr + "-1(" + kstr + "-1(" + jstr + "-1)))";
         cout << "********genecode: " << code << "  *********"<< endl;
         all_genes_origin.push_back(code);
         gTree pTree = new gTree(pInp);           
         pTree->readcode(code);
         all_genes_tree(pTree);
         int height = gtree.maxDepth();
         if(height > MaxHeight) MaxHeight = height;      //update the MaxHeight of all tree        
       }
     }
    }
    vector<gTree*>::iterator it; 
    for(it = all_genes_tree.begin(); it != all_genes_tree.end(); it++)
    {
       serial_tree.push_back(serialze(MaxHeight));
       delete *it;
    }    
 
     return;
}
*/
//for benchmark case



void S_Genetic::gen_all_perm()
{
 int nseeds = pInp->inp_sgen.n_seeds;
 int nfrags = pInp->inp_sgen.n_regs;
 int ncaps = pInp->inp_sgen.n_caps;
 string code;
 int ntotal = nseeds*nfrags*nfrags*ncaps;
 string codes[ntotal];
 int count=0;
// int count = 0;

      cout << "========Generation All Param=========" << endl;

  for (int x=0; x < nseeds; x++)
  {
    string xstr = StringTools::int2str(x+1,1,"0");    
  for (int i=0; i < nfrags; i++)
//    for (int i=0; i < 1; i++)
    {
      string istr = StringTools::int2str(nseeds+i+1,1,"0");
//        for (int j=0; j < 2; j++)
    for (int k=0; k < nfrags; k++)
    {
      string kstr = StringTools::int2str(nseeds+k+1,1,"0");     
      for (int j=0; j < ncaps; j++)
      {
        string jstr = StringTools::int2str(j+1+nfrags+nseeds, 1, "0");
        code = xstr + "(" + istr + "-1(" + kstr + "-1(" + jstr + "-1)))";
        cout << "********genecode: " << code << "  *********"<< endl;       
        codes[count]=code;
        Node<GANode>* pNode = new Node<GANode>();
        pNode->data.gene = code;
        count++;
     // gen_node(pNode);

      }
    }
   }
  }
    cout << "=======================================" << endl;
    
    for (int f=0 ; f< ntotal; f++){
     cout << f+1 << "   " << codes[f] << endl;
     }

    return;
}
#if 0
int* Genetic::gen_frezlist(xyznode* pNode)
{
   int* frezlist = new int[pNode->natoms];

   int ntemp = pNode->ntemp-1;
   int temp_natoms = pInp->xyzf_temp[ntemp].natoms;
   vector<vector<int> > joint = pInp->xyzf_temp[ntemp].joint;

   for (size_t i=0; i < pNode->natoms; i++) frezlist[i]=0;
   for (int i=0; i < temp_natoms; i++)  frezlist[i] = 1;
   for (size_t i =0; i< joint.size(); i++)
   {
      for(size_t j=0; j< joint[i].size(); j++)
         frezlist[joint[i][j]] = 0;
   }   

   return frezlist;
}

//gen_node after gencode
void Genetic::gen_node(xyznode* pNode)
{

  save_gene_in_pop(pNode->ntemp, pNode->gene);
  save_genelist(pNode->ntemp, pNode->gene);

  ncat++;
  pNode->comment = StringTools::int2str(ncat, 4, "0"); 
  pNode->flag = 0;
  pNode->next = NULL;

  pNext->addnode(pNode);

//for test ONLY! 
#if 1
  double score;
  score = getfit_Co(pNode);
  cout << "The score for the " << ncat << "-th catalyst is " << score << endl;
#endif
}
#endif

void S_Genetic::array_reset_ETL()
{
   for (int i=0; i< pInp->inp_ga.pop_size; i++)
   {
     ptr_array_ETL[i] = NULL;
     prob_array_ETL[i] = -1.0;
     select_list_ETL[i] = -9999;
     permutation_ETL[i] = i;
   }
}

void S_Genetic::array_reset_HTL()
{
   for (int i=0; i< pInp->inp_ga.pop_size; i++)
   {
     ptr_array_HTL[i] = NULL;
     prob_array_HTL[i] = -1.0;
     select_list_HTL[i] = -9999;
     permutation_HTL[i]=i;
   }
}


void S_Genetic::sort_and_select_ETL(LinkedList<Node<GANode> >* pNtrp)
{
   cout << "" << endl;
   cout << "Sorting by ETL score..." << endl;
   
  int pool = pInp->inp_ga.pop_size;
  int i=0;
  double max, max_HOMO, max_LUMO;
  double sum = 0.0;
   string thisoledstr="";
   string topoledstr="";
   string max_gene="";
   ofstream average;
   string av_fname = "average";

   average.open(av_fname.c_str(), ios::app);

   ofstream datafile;
   string data_fname = "data";
   datafile.open(data_fname.c_str(), ios::app);
   datafile.setf(ios::fixed);
   datafile.setf(ios::left);
   datafile << setprecision(6);

    ofstream datafile2;
   string data_fname2 = "dataParsers";
   datafile2.open(data_fname2.c_str(), ios::app);
   datafile2.setf(ios::fixed);
   datafile2.setf(ios::left);
   datafile2 << setprecision(6);

   ofstream best_file;
   string best_fname = "best_ETL";

   best_file.open(best_fname.c_str(), ios::app);

   ofstream bestsofar_file;
   string bestsofar_fname="bestsofar_ETL";
   bestsofar_file.open(bestsofar_fname.c_str(), ios::app);
   
  Node<GANode>* p1 = pNtrp->gethead();
  Node<GANode>* newhead;
  newhead = pNtrp->sortlist_ETL(p1);
  pNtrp->newhead(newhead);
  p1 = pNtrp->gethead();  
  max=p1->data.score_ETL;
  max_gene = p1->data.gene;  
  max_HOMO = p1->data.HOMO;
  max_LUMO = p1->data.LUMO;
   if(p1->data.score_ETL >= bestsofar){
   if(pInp->inp_sgen.bestRestart){
    if(ngen > 0){
   bestsofar = p1->data.score_ETL;
   bestgene = p1->data.gene;
   bestHOMOETL = p1->data.HOMO;
   bestLUMOETL = p1->data.LUMO;
   bestsofar_file << bestsofar << "     " << bestgene << "   " << bestHOMOETL  << "   "<< bestLUMOETL <<endl;
     }
   }else{
   bestsofar = p1->data.score_ETL;
   bestgene = p1->data.gene;
   bestHOMOETL= p1->data.HOMO;
   bestLUMOETL=p1->data.LUMO;
   bestsofar_file << bestsofar << "     " << bestgene << "   " << bestHOMOETL << "   "<< bestLUMOETL <<endl;
 }
 }else{
  bestsofar_file << bestsofar << "     " << bestgene << "   " << bestHOMOETL << "   "<< bestLUMOETL <<endl; 
}

//array_reset_ETL();
  string cname="";
  while (i != pool)
  {
   ptr_array_ETL[i] = p1;
    
    if(i<1 && p1->data.score_ETL > 900){
    
      cname = "OLED_XYZ/ETL_" + p1->data.gene +"_"+to_string((int)(p1->data.score_ETL))+ ".xyz";
      ofstream outfile;
      outfile.open(cname.c_str());
      outfile.setf(ios::fixed);
      outfile.setf(ios::left);
      outfile << setprecision(6);
      outfile <<  p1->data.natoms << endl;
      outfile <<  p1->data.gene << endl;
      for (size_t i=0; i< p1->data.natoms; i++)
      {
      outfile.width(2);
      outfile << std::left << p1->data.anames[i] << "\t ";
      outfile.width(9);
      outfile << std::right << p1->data.coords[3*i+0] << "\t";
      outfile.width(9);
      outfile << std::right << p1->data.coords[3*i+1] << "\t";
      outfile.width(9);
      outfile << std::right << p1->data.coords[3*i+2] << "\t" << endl;
      }
}




    if(!search_from_databank(p1->data.gene) && p1->data.HOMO!=0){

     string sdata="";
     datafile.width(20);
     datafile << std::left << p1->data.gene << "     ";
     sdata=sdata+p1->data.gene+"     "+to_string(p1->data.HOMO)+"     "+to_string(p1->data.HOMO)+"     ";
     datafile.width(9);
     datafile << std::right << p1->data.HOMO << "     ";
     datafile.width(9);
     datafile << std::right << p1->data.LUMO << "     ";
     datafile.width(9);

     datafile2.width(20);
     datafile2 << std::left << p1->data.gene << "     ";
     sdata=sdata+p1->data.gene+"     "+to_string(p1->data.HOMO)+"     "+to_string(p1->data.HOMO)+"     ";
     datafile2.width(9);
     datafile2 << std::right << p1->data.HOMO << "     ";
     datafile2.width(9);
     datafile2 << std::right << p1->data.LUMO << "     ";
     datafile2.width(9);

     if(p1!=NULL){
      
     sdata=sdata+to_string(p1->data.score_ETL);
 
       
    datafile << std::right << p1->data.score_HTL << "   " << p1->data.score << "   "<< p1->data.score_ETL << endl;
     datafile2 << std::right << p1->data.score_HTL << "   " << p1->data.score << "   "<< p1->data.score_ETL << endl;     
}else{
    sdata=sdata+"99.999";
     datafile2 << std::right << "99.999" << endl;
     datafile << std::right << "99.999" << endl;}
   databanks.push_back(sdata);
  }
     score_list.push_back(ptr_array_ETL[i]->data.score_ETL);
     sum=ptr_array_ETL[i]->data.score_ETL + sum;
     p1 = p1->next;
     i++;
   }
   p1 = pNtrp->gethead();
   i = 0;
   while (i != pool)
   {
      cout << p1->data.gene << "HTL = " << p1->data.score_HTL << " HOST = " << p1->data.score << " ETL = " << p1->data.score_ETL << endl; 
    prob_array_ETL[i] = (p1->data.score_ETL)/sum;
      p1 = p1->next;
     i++;
   }
   cout << "AverageScore " << sum/pool << " ; MaxScore "<< max << endl;
   average <<  sum/pool << "   " << endl;
   best_file << max << "    " << max_gene << "   "<< max_HOMO  << "   " << max_LUMO << endl;
   
   thisoledstr = "  " + max_gene + "   " + to_string(max_HOMO)+"/"+to_string(max_LUMO)+"  "+to_string(max);
   thisoleds.push_back(thisoledstr); 
   
   topoledstr = "  " + bestgene + "   " + to_string(bestHOMOETL)+"/"+to_string(bestLUMOETL)+"  "+to_string(bestsofar);
   topoled.push_back(topoledstr);
   select_parents();
   best_file.close();
   bestsofar_file.close();
   average.close();

}

void S_Genetic::sort_and_select_HTL(LinkedList<Node<GANode> >* pNtrp)
{
  
 cout << "" << endl;
 
 cout << "Sorting by HTL score..." << endl;
  int pool = pInp->inp_ga.pop_size;
  int i=0;
  double max, max_HOMO, max_LUMO;
 string thisoledstr="";
 string topoledstr="";
   double sum = 0.0;
   string max_gene="";

   ofstream average;
   string av_fname = "average";

   average.open(av_fname.c_str(), ios::app);

   ofstream datafile;
   string data_fname = "data";
   datafile.open(data_fname.c_str(), ios::app);
   datafile.setf(ios::fixed);
   datafile.setf(ios::left);
   datafile << setprecision(6);

    ofstream datafile2;
   string data_fname2 = "dataParsers";
   datafile2.open(data_fname2.c_str(), ios::app);
   datafile2.setf(ios::fixed);
   datafile2.setf(ios::left);
   datafile2 << setprecision(6);

   ofstream best_file;
   string best_fname = "best_HTL";

   best_file.open(best_fname.c_str(), ios::app);
     ofstream bestsofar_file;
   string bestsofar_fname="bestsofar_HTL";
   bestsofar_file.open(bestsofar_fname.c_str(), ios::app);

  Node<GANode>* p1 = pNtrp->gethead();
  Node<GANode>* newhead;
  newhead = pNtrp->sortlist_HTL(p1);
  pNtrp->newhead(newhead);
  p1 = pNtrp->gethead();
  max=p1->data.score_HTL;
  max_gene = p1->data.gene;
  max_HOMO = p1->data.HOMO;
  max_LUMO = p1->data.LUMO;
  if(p1->data.score_HTL >= bestsofarHTL){
   if(pInp->inp_sgen.bestRestart){
    if(ngen > 0){
   bestsofarHTL = p1->data.score_HTL;
   bestgeneHTL = p1->data.gene;
   bestHOMOHTL=p1->data.HOMO;
   bestLUMOHTL=p1->data.LUMO;
   bestsofar_file << bestsofarHTL << "     " << bestgeneHTL << "   " << bestHOMOHTL << "   "<< bestLUMOHTL <<endl;  
   }
   }else{
   bestsofarHTL = p1->data.score_HTL;
   bestgeneHTL = p1->data.gene;
   bestHOMOHTL=p1->data.HOMO;
   bestLUMOHTL=p1->data.LUMO;
   bestsofar_file << bestsofarHTL << "     " << bestgeneHTL << "   " << bestHOMOHTL << "   "<< bestLUMOHTL <<endl;
 }
 }else{
   bestsofar_file << bestsofarHTL << "     " << bestgeneHTL << "   " << bestHOMOHTL << "   "<< bestLUMOHTL <<endl; 
}

//array_reset_HTL();
  string cname="";
  while (i != pool)
  {
   ptr_array_HTL[i] = p1;
     if(i<1 && p1->data.score_HTL > 900){

      cname = "OLED_XYZ/HTL_" + p1->data.gene +"_"+to_string((int)(p1->data.score_HTL))+ ".xyz";
      ofstream outfile;
      outfile.open(cname.c_str());
      outfile.setf(ios::fixed);
      outfile.setf(ios::left);
      outfile << setprecision(6);
      outfile <<  p1->data.natoms << endl;
      outfile <<  p1->data.gene << endl;
      for (size_t i=0; i< p1->data.natoms; i++)
      {
      outfile.width(2);
      outfile << std::left << p1->data.anames[i] << "\t ";
      outfile.width(9);
      outfile << std::right << p1->data.coords[3*i+0] << "\t";
      outfile.width(9);
      outfile << std::right << p1->data.coords[3*i+1] << "\t";
      outfile.width(9);
      outfile << std::right << p1->data.coords[3*i+2] << "\t" << endl;
      }
}
 

   if(!search_from_databank(p1->data.gene) && p1->data.HOMO!=0){

     string sdata="";
     datafile.width(20);
     datafile << std::left << p1->data.gene << "     ";
     sdata=sdata+p1->data.gene+"     "+to_string(p1->data.HOMO)+"     "+to_string(p1->data.HOMO)+"     ";
     datafile.width(9);
     datafile << std::right << p1->data.HOMO << "     ";
     datafile.width(9);
     datafile << std::right << p1->data.LUMO << "     ";
     datafile.width(9);

     datafile2.width(20);
     datafile2 << std::left << p1->data.gene << "     ";
     sdata=sdata+p1->data.gene+"     "+to_string(p1->data.HOMO)+"     "+to_string(p1->data.HOMO)+"     ";
     datafile2.width(9);
     datafile2 << std::right << p1->data.HOMO << "     ";
     datafile2.width(9);
     datafile2 << std::right << p1->data.LUMO << "     ";
     datafile2.width(9);

     if(p1!=NULL){

     sdata=sdata+to_string(p1->data.score_HTL);
     datafile << std::right << p1->data.score_HTL << "   " << p1->data.score << "   "<< p1->data.score_ETL << endl;  
     datafile2 << std::right << p1->data.score_HTL << "   " << p1->data.score << "   "<< p1->data.score_ETL << endl;     

}else{
    sdata=sdata+"99.999";
     datafile2 << std::right << "99.999" << endl;
     datafile << std::right << "99.999" << endl;}
  databanks.push_back(sdata);
  }
     score_list.push_back(ptr_array_HTL[i]->data.score_HTL);
     sum=ptr_array_HTL[i]->data.score_HTL + sum;
     p1 = p1->next;
     i++;
   }
   p1 = pNtrp->gethead();
   i = 0;
   while (i != pool)
   {
     cout << p1->data.gene << "HTL = " << p1->data.score_HTL << " HOST = " << p1->data.score << " ETL = " << p1->data.score_ETL << endl;

    prob_array_HTL[i] = (p1->data.score_HTL)/sum;
      p1 = p1->next;
     i++;
   }
   cout << "AverageScore " << sum/pool << " ; MaxScore "<< max << endl;
   average <<  sum/pool << "   " << endl;
   
    topoledstr = "  " + bestgeneHTL + "   " + to_string(bestHOMOHTL)+"/"+to_string(bestLUMOHTL)+"  "+to_string(bestsofarHTL);
   topoled.push_back(topoledstr);


   best_file << max << "    " << max_gene << "   "<< max_HOMO  << "   " << max_LUMO << endl;
   thisoledstr = "  " + max_gene + "   " + to_string(max_HOMO)+"/"+to_string(max_LUMO)+"  "+to_string(max);
   thisoleds.push_back(thisoledstr);

   select_parents();
   best_file.close();
   bestsofar_file.close();
   average.close();

}

//  xyznode* newhead = pNext->sortlist(head);
//  pNext->newhead(newhead);
//  string filename = StringTools::int2str(ngen,4,"0");
//  filename = "Work/catslist.xyz"+filename;
//  pNext->print_xyzfile(filename);
//}

#if 0
void Genetic::print_all_score()
{
   int i = 0;
   size_t s = score_list.size();
   ofstream outfile("Work/scorelist.dat");
   for(size_t j=0; j<s; j++)
   {
     i++;
     outfile << i << "\t" ;
     outfile.width(9);
     outfile << score_list[j] << endl; 
//     outfile.width(9);
//     outfile << MW_list[j] << endl;
   }
   return;
}
#endif

LinkedList<Node<GANode> >* S_Genetic::after_sort(int ngen,LinkedList<Node<GANode> >* pNext_ETL, LinkedList<Node<GANode> >* pCtrp)
{
//    Node<GANode>* head = pNext->gethead();
//    Node<GANode>* newhead = pNext->sortlist(head);
//    pNext->newhead(newhead);
#if 0
    string filename = StringTools::int2str(ngen,4,"0");
    filename = "./Scratch/opted_list.xyz"+filename;
    pNext->print_xyzfile(filename);
#endif

//  cout << "Rename New Gen as Current Gen" << endl;
//  cout << endl << endl;

  pTmp = pNext_ETL;
  pNext_ETL = pCtrp;
  pCtrp = pTmp;     //swap pNext and pCurrent
  pTmp = NULL;

  return pCtrp;
#if 0
//  int i = 0;
//  double sum = 0.0;
  
//  xyznode* p1;

//  array_reset();

//  p1 = pCurrent->gethead();

//  while (i != pInp->inp_ga.pop_size)
//  {
//    ptr_array[i] = p1;
    cout << p1 << endl;
    sum = p1->score + sum;
    p1 = p1->next; 
    i++;
  }
  
  p1 = pCurrent->gethead();
  i = 0;
  while (i != pInp->inp_ga.pop_size)
  {
//    if (i==0)
    prob_array[i] = (p1->score)/sum;
    p1 = p1->next;
    i++;
//    else
//    prob_array[i] = (p1->score)/sum + prob_array[i-1];       
  } 
#endif
}

Node<GANode>* S_Genetic::cpy_node(Node<GANode>* pNode)
{
   Node<GANode>* pNew = new Node<GANode>();
   ncat++;
//   pNew->charge = pNode->charge;
//   pNew->multi = pNode->multi;
   pNew->data.natoms = pNode->data.natoms;
//   pNew->ntemp = pNode->ntemp;
   pNew->data.coords = pNode->data.coords;
   pNew->data.anames = pNode->data.anames;
   pNew->data.gene = pNode->data.gene;
   pNew->data.comment = pNode->data.comment;
   pNew->data.score = pNode->data.score;
//   pNew->energy = pNode->energy;
//   pNew->binding = pNode->binding;
//   pNew->h_affinity = pNode->h_affinity;
//   pNew->hydricity = pNode->hydricity;
   pNew->data.HOMO = pNode->data.HOMO;
   pNew->data.LUMO = pNode->data.LUMO;   
   pNew->data.flag = -1;              /*set flag to avoid opt in next round*/
   pNew->data.MW = pNode->data.MW; 
   pNew->next = NULL;

//   addtohash(genehash_pop, pNew->data.gene);
  
   return pNew;
}

int S_Genetic::mopac_calc(string filedir, vector<int> & joblist,  LinkedList<Node<GANode> >* pList,char OLEDtype, int jobtype)
{
   bool success = true;
   Node<GANode>* pNode;
   pNode = pList->gethead();

   cout << "mopac cal" << endl;

     if(jobtype==1){
   Mopac::mopac_genslurm(filedir, joblist);   //hard-coded for 0-gen case (joblist.size() = pop_size)
 }
  else{
   Mopac::mopac_genqsh(filedir, joblist);
}
   success = Mopac::run_opt_jobs(filedir, joblist, jobtype);   

 
   if(success)
   { 
//     Mopac::hard_code_save_alltolist(filedir,pList);
   Mopac::save_all_into_list(filedir,pList);
//   return 1;
   } 
   
   pNode = pList->gethead();
  if(isFirstrun){
    refHOST_HOMO=pInp->inp_sgen.HOST_HOMO;
    refHOST_LUMO=pInp->inp_sgen.HOST_LUMO;
  }   

for(int i=0; i< pList->listlength;i++){
     oled_score(pNode,OLEDtype); 
     pNode = pNode->next;
   }
// }


     return 1;
}
#if 0
int Genetic::dft_opt(string filedir, vector<int> & joblist, LinkedList* pList)
{
   bool success = true;

   DFT::opt_genqsh(filedir, joblist);
   success = DFT::run_opt_jobs(filedir, joblist);

   if(success)
   {
   DFT::save_all_into_list(filedir, pList);
   check_connect(pList);
   return 1;
   }
   else
   return 0;
}

/*check if TM and ligands disconnect after DFT opt*/
void Genetic::check_connect(LinkedList* pList)
{
   xyznode* pNode;
   pNode = pList->gethead();
   int count = 0;

   cout << "Checking bonding changes..." << endl;
  
   while(count < pList->listlength)
   { 
      int ntemp = pNode->ntemp - 1;
      vector<string> element = pInp->xyzf_temp[ntemp].element;
      vector<double> coords = pInp->xyzf_temp[ntemp].coords;
      int natoms = pInp->xyzf_temp[ntemp].natoms;

      ICoord ict;
      ICoord icopt;
  
      ict.init(natoms, element, coords);
      icopt.init(pNode->natoms, pNode->element, pNode->coords);

      for (int i=0; i < ict.nbonds; i++)
      {
        if ((ict.bonds[i][0] != icopt.bonds[i][0]) || (ict.bonds[i][1] != icopt.bonds[i][1])) 
        {
          pNode->flag = -1;
          pNode->score = 0; 
          cout << "****Catalyst " << count << ": " << "Connection between TM and ligand(s) is changed after DFT opt!" << endl; 
          break;
        }
      }
      count++;
      pNode = pNode->next;
      ict.freemem();
      icopt.freemem();
   }

   cout << "Done checking." << endl;
   return;
}


int Genetic::get_fit(string filedir, vector<int> &list, Fitness & fitness)
{  

   xyznode* pNode;
   pNode = pNext->gethead();
   int count=0;

   while(count < pNext->listlength)
   {
     if(pNode->flag == -1)
     {
       pNode = pNode->next;
       count++;
     }
     else if (pNode->flag == 0)
     {
       fitness.geninp(filedir, pNode);
       pNode = pNode->next;
       count++;
     }
   }

  fitness.genqsh(filedir, list);
  fitness.run_all(filedir, list);

   pNode = pNext->gethead();
   count=0;

   while (count < pNext->listlength)
   {
     if (pNode->flag == -1)
     {
        pNode = pNode->next;
        count++;
     }
     else if (pNode->flag == 0)
     {
        fitness.get_fit(filedir, pNode);
        pNode = pNode->next;
        count++;
     }
   }
   return 1;
}
#endif
#if 0
/*write in binary format maybe?*/
bool Genetic::save_restart(xyznode* pRef)
{
   ofstream restart;
   string restart_fname = "ga.restart";
   restart.open(restart_fname.c_str(),ios::trunc);  
   restart.setf(ios::fixed);
   restart.setf(ios::left);
   restart << setprecision(6);

   restart << "Current Generation   " << ngen << endl;
   
   restart << "Reference Catalyst Info. " << endl;
   restart << "Charge" << "\t" << pRef->charge << endl;
   restart << "Multi" << "\t" << pRef->multi << endl;
   restart << "Score" << "\t" << pRef->score << endl;
   restart << "Energy-I" << "\t" << pRef->energy << endl;

/*hard-coded*/
   restart << "Energy-II" << "\t" << pRef->energy2 << endl;
   restart << "Energy-III" << "\t" << pRef->energy3 << endl;   

 
   restart << "E_binding" << "\t" << pRef->binding << endl;
   restart << "E_haffinity" << "\t" << pRef->h_affinity << endl;
   restart << "E_hydricity" << "\t" << pRef->hydricity << endl;
   restart << "Natoms" << "\t" << pRef->natoms << endl;
   for(size_t i = 0; i < pRef->natoms; i++)
   {
     restart.width(2);
     restart << std::left << pRef->element[i] << "\t";
     restart.width(9);
     restart << std::right << pRef->coords[3*i+0] << "\t";
     restart.width(9);
     restart << std::right << pRef->coords[3*i+1] << "\t";
     restart.width(9);
     restart << std::right << pRef->coords[3*i+2] << endl;
   }

//   restart << "Current Generation   " << ngen << endl;
   
   xyznode* pNode = pNext->gethead();
   int listlength = pNext->listlength;
   int i = 0;
   vector<int>::iterator it;

   while(i != listlength)
   {
   string nstr = StringTools::int2str(i, 1, "0");
   restart << "$node" << nstr << endl;
   restart << "Charge" << "\t" << pNode->charge << endl;
   restart << "Multi" << "\t" << pNode->multi<< endl;
   restart << "Template" << "\t" << pNode->ntemp << endl;
   restart << "Gene"; 
   for(it=pNode->gene.begin(); it!=pNode->gene.end(); ++it)
   {
     restart << "\t" << *it; 
   }
   restart << endl; 
   restart << "Score" << "\t" << pNode->score << endl;
   restart << "Energy-I" << "\t" << pNode->energy << endl;

/*hard-coded*/
   restart << "Energy-II" << "\t" << pNode->energy2 << endl;
   restart << "Energy-III" << "\t" << pNode->energy3 << endl;

   restart << "E_binding" << "\t" << pNode->binding << endl;
   restart << "E_haffinity" << "\t" << pNode->h_affinity << endl;
   restart << "E_hydricity" << "\t" << pNode->hydricity << endl;
   restart << "Natoms" << "\t" << pNode->natoms << endl;

    for (size_t j=0; j < pNode->natoms; j++)
    {
     restart.width(2);
     restart << std::left << pNode->element[j] << "\t";
     restart.width(9);
     restart << std::right << pNode->coords[3*j+0] << "\t";
     restart.width(9);
     restart << std::right << pNode->coords[3*j+1] << "\t";
     restart.width(9);
     restart << std::right << pNode->coords[3*j+2] << endl;
    }
     i++; 
     pNode = pNode->next;
   }

   restart.close();

   ofstream glist;
   glist.open("genelist.restart", ios::trunc);
   glist.setf(ios::fixed);
   glist.setf(ios::left);
  
//   vector<int>::iterator it;
   vector<vector<int> >::iterator it2;

   for (it2 = genelist.begin(); it2 != genelist.end(); ++it2)
   {
      for (it = it2->begin(); it != it2->end(); ++it)
      {
         glist.width(3);
         glist << std::right << *it << "\t"; 
      }
      glist << endl;
   }
  
   glist.close();
   }
#endif

bool S_Genetic::clear_folder_restart(int ngen)
{
   string ngen_str = StringTools::int2str(ngen, 3, "0");
   string cwd = StringTools::getcwd_str();
   string path = cwd + "/Work/MM/" + ngen_str; 
   struct stat statbuf;
   if (stat(path.c_str(), &statbuf) != -1)
   {
     if(S_ISDIR(statbuf.st_mode))
     {
       string cmd = "rm -rf " + path;
       system(cmd.c_str());
     }
   }
   path = cwd + "/Work/mopac/" + ngen_str; 
   if (stat(path.c_str(), &statbuf) != -1)
   {
     if(S_ISDIR(statbuf.st_mode))
     {
        string cmd = "rm -rf " + path;
        system(cmd.c_str());
     }
   }
   path = cwd + "/Work/opt/" + ngen_str;
   if (stat(path.c_str(), &statbuf) != -1)
   {
     if(S_ISDIR(statbuf.st_mode))
     {
        string cmd = "rm -rf " + path;
        system(cmd.c_str());
     }
   }
   path = cwd + "/rgroups/data";
   if (stat(path.c_str(), &statbuf) != -1)
   {
     if(S_ISDIR(statbuf.st_mode))
     {
        string cmd = "rm -rf " + path;
        system(cmd.c_str());
     }
   }

    cout << "Clear incompleted generation!" << endl;
}
#if 0
int Genetic::eead_restart()
{
   int ngen;
//   pRef = new xyznode();
   string line;
   vector<string> tok_line;

   string cwd = StringTools::getcwd_str();
   string fname = cwd + "/ga.restart";
   ifstream restart(fname.c_str(), ios::in);
   vector<string>::iterator it;   

   string listname = cwd + "/genelist.restart";
   ifstream glist(listname.c_str(),ios::in);
 
   if((!restart)||(!glist))
   {
      cout << "Cannot find restart file " << endl;
      return 0;
   }
   else
   {
      cout << "****** read restart file ******" << endl << endl;
      pRef = new xyznode();
      getline(restart, line);
      cout << line << endl;
      tok_line = StringTools::tokenize(line, " ");
      ngen = atoi(tok_line[1].c_str());

      while(getline(restart,line))
      {
       if (line.find("Reference") != string::npos)
       {
          getline(restart, line);
          if (line.find("Charge") != string::npos)
          {
            tok_line = StringTools::tokenize(line," ");
            pRef->charge = atoi(tok_line[1].c_str());
            cout << line << endl;
          } 
           
          getline(restart, line);
          if (line.find("Multi") != string::npos)
          {
             tok_line = StringTools::tokenize(line," ");
             pRef->multi = atoi(tok_line[1].c_str());
             cout << line << endl;
          }
    
          getline(restart, line);
          if (line.find("Score") != string::npos)
          {
              tok_line = StringTools::tokenize(line, " ");
              pRef->score = atof(tok_line[1].c_str());
              cout << line << endl;
          }

          getline(restart, line);
          if (line.find("Energy-I") != string::npos)
          {
              tok_line = StringTools::tokenize(line, " ");
              pRef->energy = atof(tok_line[1].c_str());
              cout << line << endl;
          }


          getline(restart, line);
          if (line.find("Energy-II") != string::npos)
          {
             tok_line = StringTools::tokenize(line, " ");
             pRef->energy2 =  atof(tok_line[1].c_str());
             cout << line << endl;
          }

          getline(restart, line);
          if (line.find("Energy-III") != string::npos)
          {
             tok_line = StringTools::tokenize(line," ");
             pRef->energy3 = atof(tok_line[1].c_str());
             cout << line << endl;
          }

          getline(restart, line);
          if (line.find("E_binding") != string::npos)
          {
              tok_line = StringTools::tokenize(line, " ");
              pRef->binding = atof(tok_line[1].c_str());
              cout << line << endl;
          }
 
          getline(restart, line);
          if (line.find("E_haffinity") != string::npos)
          {
              tok_line = StringTools::tokenize(line," ");
              pRef->h_affinity = atof(tok_line[1].c_str());
              cout << line << endl;
          }
         
          getline(restart, line);
          if (line.find("E_hydricity") != string::npos)
          {
              tok_line = StringTools::tokenize(line, " ");
              pRef->hydricity = atof(tok_line[1].c_str());
              cout << line << endl;
          }

          getline(restart, line);
          if (line.find("Natoms") != string::npos)
          {
             tok_line = StringTools::tokenize(line, " ");
             pRef->natoms = atoi(tok_line[1].c_str());
             cout << line << endl;
             pRef->init(int(pRef->natoms));
             for (size_t i=0; i < pRef->natoms; i++)
             {
                getline(restart, line);
                tok_line = StringTools::tokenize(line, " ");
                pRef->element[i] = tok_line[0];
                pRef->coords[3*i+0] = atof(tok_line[1].c_str());
                pRef->coords[3*i+1] = atof(tok_line[2].c_str());
                pRef->coords[3*i+2] = atof(tok_line[3].c_str());
                cout << line << endl;
             }
          }
       }
       else if (line.find("$node") != string::npos)
       {
            xyznode* pNode = new xyznode();
            pNext->addnode(pNode);
            cout << line << endl;

           getline(restart, line);
           if (line.find("Charge") != string::npos)
           {
            tok_line = StringTools::tokenize(line," ");
            pNode->charge = atoi(tok_line[1].c_str());
            cout << line << endl;
           }

           getline(restart, line);
           if (line.find("Multi") != string::npos)
           {
             tok_line = StringTools::tokenize(line, " ");
             pNode->multi = atoi(tok_line[1].c_str());
             cout << line << endl;
           }

           getline(restart, line); 
           if (line.find("Template") != string::npos)
           {
             tok_line = StringTools::tokenize(line," ");
             pNode->ntemp = atoi(tok_line[1].c_str());
             cout << line << endl;
           }

           getline(restart, line);
           if (line.find("Gene")!= string::npos)
           {
             tok_line = StringTools::tokenize(line," ");
           for (it=tok_line.begin()+1; it!=tok_line.end(); ++it)
           {
             int tmp = atoi((*it).c_str());
             pNode->gene.push_back(tmp);
           } 
           cout << line << endl;
           }

           getline(restart, line);          
           if (line.find("Score")!= string::npos)
           {
             tok_line = StringTools::tokenize(line, " ");
             pNode->score = atof(tok_line[1].c_str());
             cout << line << endl;
           }
            
           getline(restart, line);
           if (line.find("Energy-I") != string::npos)
           {
           tok_line = StringTools::tokenize(line," ");
           pNode->energy = atof(tok_line[1].c_str());
           cout << line << endl;
           }

          getline(restart, line);
          if (line.find("Energy-II") != string::npos)
          {
             tok_line = StringTools::tokenize(line, " ");
             pNode->energy2 =  atof(tok_line[1].c_str());
             cout << line << endl;
          }

          getline(restart, line);
          if (line.find("Energy-III") != string::npos)
          {
             tok_line = StringTools::tokenize(line," ");
             pNode->energy3 = atof(tok_line[1].c_str());
             cout << line << endl;
          }

           getline(restart, line);
           if (line.find("E_binding") != string::npos)
           {
             tok_line = StringTools::tokenize(line, " ");
             pNode->binding = atof(tok_line[1].c_str());
             cout << line << endl;
           }
    
           getline(restart, line);
           if (line.find("E_haffinity") != string::npos)
           {
              tok_line = StringTools::tokenize(line, " ");
              pNode->h_affinity = atof(tok_line[1].c_str());
              cout << line << endl;
           }
           
           getline(restart, line);
           if (line.find("E_hydricity") != string::npos)
           { 
              tok_line = StringTools::tokenize(line, " ");
              pNode->hydricity = atof(tok_line[1].c_str());
              cout << line << endl;
           }

           getline(restart, line);
           if (line.find("Natoms") != string::npos)
          {
          tok_line = StringTools::tokenize(line," ");
          cout << line << endl;
          pNode->natoms = atoi(tok_line[1].c_str());
          pNode->init(int(pNode->natoms));
          for (size_t i=0; i < pNode->natoms; i++)
          {
             getline(restart,line);
             tok_line = StringTools::tokenize(line," ");
             pNode->element[i] = tok_line[0];
             pNode->coords[3*i+0] = atof(tok_line[1].c_str());
             pNode->coords[3*i+1] = atof(tok_line[2].c_str());
             pNode->coords[3*i+2] = atof(tok_line[3].c_str());
             cout << line << endl;
          } 
//          cout << line << endl;
          }
       }//end of else if
     } // end of while loop
     
     vector<int> tmp;
     int code;

     cout << "read saved genelist for restart file" << endl;

     while(getline(glist, line))
     {
        tok_line = StringTools::tokenize(line, " ");
        size_t size = tok_line.size();
        for (size_t i=0; i < size; i++)
        {
           code = atoi(tok_line[i].c_str());
           tmp.push_back(code);   
        }
        genelist.push_back(tmp);
        tmp.clear();
        cout << line << endl;
     } 
   }
   pre_sort();  
   return ngen+1;
}
#endif
#if 0
void Genetic::read_ref(string folder, Node<GANode>* pRef)
{
//  pRef->charge = pInp->xyzf_ref.charge;
//  pRef->multi = pInp->xyzf_ref.multi;
  pRef->data.natoms = pInp->ref_xyz.natoms;
  pRef->data.anames = pInp->ref_xyz.anames;
  pRef->data.coords = pInp->ref_xyz.coords;
//  pRef->ntemp = 0;
  pRef->next = NULL;
  pRef->data.comment = "0000";
 
  pRef->data.flag = 0;

  cout << "****   Referene Molecule:  ****" << endl;
//  cout << "Charge: " << pRef->charge << endl;
//  cout << "Multiplicity: " << pRef->multi << endl;
  cout << "Number of Atoms: " << pRef->data.natoms << endl;
  cout << "Cartesian Coordinate: " << endl;

  for (size_t i = 0; i < pRef->data.natoms; i++)
  {
    cout.width(2);
    cout << std::left << pRef->data.anames[i] << "\t";
    cout.width(9); 
    cout << std::right << pRef->data.coords[3*i+0]<< "\t";
    cout.width(9);
    cout << std::right << pRef->data.coords[3*i+1] << "\t";
    cout.width(9);
    cout << std::right << pRef->data.coords[3*i+2] << endl; 
  }

//    Mopac::gen_inp(folder, pInp, pRef);

//    string cmd = "/tmp/MOPAC2012.exe " + folder+ "m0000";
    
//    system(cmd.c_str());
   
    cout << "reading mopac result ... " << endl;

    Mopac::readout(folder, pRef);  

    cout << "optimized coordinates:" << endl;
    for (size_t i = 0; i < pRef->data.natoms; i++)
   {
     cout.width(2);
     cout << std::left << pRef->data.anames[i] << "\t";
     cout.width(9);
     cout << std::right << pRef->data.coords[3*i+0]<< "\t";
     cout.width(9);
     cout << std::right << pRef->data.coords[3*i+1] << "\t";
     cout.width(9);
     cout << std::right << pRef->data.coords[3*i+2] << endl;
   }
    cout << "HOMO: " << pRef->data.HOMO << " eV" << endl;
    cout << "LUMO: " << pRef->data.LUMO << " eV" << endl;

//    DFT::opt_geninp(folder, pInp, pRef);
//    vector<int> list(1,0);
//    DFT::opt_genqsh(folder,list);
//    DFT::run_opt_jobs(folder,list);
//    DFT::save_opted_xyz_energy(folder, pRef);
}
#endif

void S_Genetic::run_ga()
{ 

  cout << "******************************" << endl;
  cout << "*****  Starting GA runs  *****" << endl;
  cout << "*******************************" << endl;
  
  LinkedList<Node<GANode> > List0, List00;
  LinkedList<Node<GANode> > List1;
  LinkedList<Node<GANode> > List2;   
  vector<string> tok_line;
  bool noRecord=false;

  fstream frecord;
  frecord.open("data");
  if(frecord.is_open()==true){
   cout << "Start reading record.... "<<endl;
   cout << "" << endl;
   read_record();
   }else{
   noRecord=true; 
   frecord.width(20);
   frecord << std::left << "gene" << "     ";
   frecord.width(9);
   frecord << std::right << "HOMO (eV)" << "     ";
   frecord.width(9);
   frecord << std::right << "LUMO (eV)" << "     ";
   frecord.width(9);
   frecord << std::right << "score" << endl;
   cout << "run newstart " << endl;
   }

   ofstream datafile2;
   string data_fname2 = "dataParsers";
   datafile2.open(data_fname2.c_str(), ios::app);
   datafile2.setf(ios::fixed);
   datafile2.setf(ios::left);
   datafile2 << setprecision(6);


   for(int i=0;i<databanks.size();i++){
   datafile2 << databanks.at(i) << endl;
   }
   datafile2.close();
   
  ngen = 0;
  
  if(ngen == 0)
  {
   string cwd = StringTools::getcwd_str();
   string folder = cwd + "/Scratch/ref/";
   pRef = new Node<GANode>();
   int status;
   string path = "./Scratch/mm";
   status = mkdir(path.c_str(),S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
  }

  else
  {
  cout << "clear_folder" << endl;
  clear_folder_restart(ngen);
  }

  pCurrent_ETL = &List0;
  pCurrent_HTL = &List00;
  pNext_ETL = &List1;
  pNext_HTL = &List2;


  string ngen_str;

  while(ngen <= pInp->inp_ga.ga_stop){
   

   ncat = 0;
   ngen_str = StringTools::int2str(ngen, 3, "0"); 
   string cwd = StringTools::getcwd_str();  
// cout << "cwd : " << cwd <<endl;
   mopac_folder_ETL = cwd + "/Scratch/mopac/ETL_" + ngen_str + "/";
   mmopt_folder_ETL = cwd + "/Scratch/mm/ETL_" + ngen_str + "/";
   
   mopac_folder_HTL = cwd + "/Scratch/mopac/HTL_" + ngen_str + "/";
   mmopt_folder_HTL = cwd + "/Scratch/mm/HTL_" + ngen_str + "/";

   int status_ETL;
   int status_HTL;
    
   status_ETL = mkdir(mmopt_folder_ETL.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
   status_ETL = mkdir(mopac_folder_ETL.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
   status_HTL = mkdir(mmopt_folder_HTL.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
   status_HTL = mkdir(mopac_folder_HTL.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
   
   cout << "***************************************"<<endl;    
   oledtypeall='1';

  
   cout << "The " << ngen << "-th Generation for ETL" << endl;   
   spawn(mmopt_folder_ETL,mopac_folder_ETL,pNext_ETL, '1');   
   cout << "****************************************"<< endl << endl;
   mopac_calc(mopac_folder_ETL, optr_list,  pNext_ETL, '1', pInp->inp_sgen.jobtype);
// updatedHOST();
   sort_and_select_ETL(pNext_ETL);
   pCurrent_ETL=after_sort(ngen,pNext_ETL,pCurrent_ETL);
   pCurrent_ETL->print_record(ngen);

   tempstring_ETL.clear();

  for(int i=0;i<pInp->inp_ga.pop_size;i++){
   ptr_array_ETL_temp.push_back(ptr_array_ETL[i]);
   tempstring_ETL.push_back(ptr_array_ETL[i]->data.gene);
}
   genehash_pop.clear();
   genehash_pop.resize(2*pInp->inp_ga.pop_size+1);
   optr_list.clear();

    
   ncat=0; 
  
   cout << "***************************************"<<endl;
   cout << "The " << ngen << "-th Generation for HTL" << endl;
   oledtypeall='2';
   spawn(mmopt_folder_HTL,mopac_folder_HTL,pNext_HTL, '2');
   cout << "****************************************"<< endl << endl;
   mopac_calc(mopac_folder_HTL, optr_list, pNext_HTL,'2', pInp->inp_sgen.jobtype);
   sort_and_select_HTL(pNext_HTL);
   pCurrent_HTL=after_sort(ngen,pNext_HTL,pCurrent_HTL);
   pCurrent_HTL->print_record(ngen);

   tempstring_HTL.clear();   

   for(int i=0;i<pInp->inp_ga.pop_size;i++){
   ptr_array_HTL_temp.push_back(ptr_array_HTL[i]);
   tempstring_HTL.push_back(ptr_array_HTL[i]->data.gene);
}

  
   cout << "The length of pNext List" << pNext_HTL->listlength << endl;
   cout << "The length of pCurrent List" << pCurrent_HTL->listlength << endl;
   
   updatedHOST();

   cout << "" << endl;
   cout << "======This Generation("<< ngen << "-th)====="<<endl;
   cout << "Best HTL : " << thisoleds.at(1) << endl;
   cout << "Best HOST : " << thisoleds.at(2) << endl;
   cout << "Best ETL : " << thisoleds.at(0) << endl;
   double totalscore =0;
   for(int i=0;i<3;i++){
   tok_line = StringTools::tokenize(thisoleds.at(i), " ");
   totalscore+=atof(tok_line[2].c_str());  
   }
   cout << "Total OLED score = " << totalscore << endl;
   cout << "================================"<<endl;
   cout << "" << endl;

   thisoleds.clear();  

    cout << "" << endl;
   cout << "======BestSofar OLED device====="<<endl;
   cout << "Best HTL : " << topoled.at(1) << endl;
   cout << "Best HOST : " << topoled.at(2) << endl;
   cout << "Best ETL : " << topoled.at(0) << endl;
   double totaltopscore =0;
   for(int i=0;i<3;i++){
   tok_line = StringTools::tokenize(topoled.at(i), " ");
   totaltopscore+=atof(tok_line[2].c_str());
   }
   cout << "Total OLED score = " << totaltopscore << endl;
   cout << "================================"<<endl;
   cout << "" << endl;

   topoled.clear();



   ngen++;
   genehash_pop.clear();
   genehash_pop.resize(2*pInp->inp_ga.pop_size+1);
   optr_list.clear();

} //end of while
   
   cout << "Total size of databanks : " << databanks.size() << endl;

   List0.clear();
   List00.clear();
   List1.clear();
   List2.clear();
 }
