#include <iostream>
#include <fstream>
#include <stdio.h>
#include <stdlib.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <dirent.h>
#include <math.h>

using namespace std;
#include "pTable.h"
#include "genetic.h"
#include "mopac.h"
#include "gcode.h"
#include "fragment.h"
//#include "icoord.h"
//#include "xyzlist.h"
//#include "CatsGen.h"
//#include "DFT.h"
#include "utils.h"

#define NUM_ALL_PERM 9450

void Genetic::fragment_prop()
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
void Genetic::init_feature_space()
{
  int nfrags = pInp->inp_sgen.n_regs;
  int ncaps = pInp->inp_sgen.n_caps;
  string code;
  Fragment* pfrag;
 // int count = 0;

    ofstream summary;
    summary.open("entire_space.csv");
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

void Genetic::get_fit4test(Node<GANode>* pNode)
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

void Genetic::hashinit()
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

long Genetic::BKDRHash(string str)
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

int Genetic::addtohash(vector<unique_ptr<LinkedList<Node<string> > > > & hashtable, string & str)
{
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

int Genetic::addtohash_record(vector<unique_ptr<LinkedList<Node<record> > > > & hash_record, Node<record>* pNode)
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

int Genetic::searchhash(vector<unique_ptr<LinkedList<Node<string> > > > & hashtable, string & genecode)
{
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

void Genetic::searchhash_record(vector<unique_ptr<LinkedList<Node<record> > > > & hash_record, Node<GANode>* pNode)
{
   cout << "run to here" << endl;
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
     
   double score = oled_fit(ptr);
 
   pNode->data.score = score;

   cout << "score = " << score << endl;
   return;
}

void Genetic::oled_score(Node<GANode>* pNode)
{
    cout << "calculating score .... " << endl;
   double score_homo;

    score_homo = 1/(1+exp(-5.0* (-9.2 - pNode->data.HOMO)));
//   if(pNode->data.HOMO < -9.5)  score_homo = 1.0;
//   else if (pNode->data.HOMO > -8.5) score_homo = 0.0;
//   else score_homo = -8.5 - pNode->data.HOMO;

//or try sigmoid function for homo, score = homo*lumo

   double score_lumo;
  double diff = pNode->data.LUMO - (-1.0);      /*hard-coded*/
   score_lumo = 1000 * exp(-(diff*diff)/0.125);  //2sigma = 0.5 (0.95 chance);

   double score = score_lumo*score_homo;

   cout << score << endl;

   pNode->data.score = score;
   return;
}

double Genetic::oled_fit(Node<record>* pNode)
{
   cout << "calculating score .... " << endl;
   double score_homo;

   score_homo = 1/(1+exp(-5.0*(-9.2 - pNode->data.HOMO))); 
//   if(pNode->data.HOMO < -9.5)  score_homo = 1.0;
//   else if (pNode->data.HOMO > -8.5) score_homo = 0.0;
//   else score_homo = -8.5 - pNode->data.HOMO;

//or try sigmoid function for homo, score = homo*lumo

   double score_lumo; 
  double diff = pNode->data.LUMO - (-1.0);
   score_lumo = 1000 * exp(-(diff*diff)/0.125);  //2sigma = 0.5 (0.95 chance);   
   
   double score = score_lumo*score_homo;

   return score;
}

void Genetic::crossover(string & genecode1, string & genecode2)
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

void Genetic::swap_ptr(TNode* & ptr1, TNode* & ptr2)
{
   TNode* tmp;
   tmp = ptr1;
   ptr1 = ptr2;
   ptr2 = tmp;
 
   return;
}

string Genetic::mutation(string genecode)
{
   string newcode = genecode;
   
   int exist = 1;
   while(exist)
   {   
   gTree gtree(pInp);
   gtree.readcode(newcode);
   int splice = Utils::randomi(gtree.nnodes); //splice = [0, nnodes-1];
  
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

   newcode = gtree.tree_to_code();
   
   cout << "mutated gene:" << newcode << endl;
 
   exist = searchhash(genehash_all, newcode); 

   if (exist) cout << "in the list already, mutate again" << endl;
   }
   cout << endl;

   return newcode;
}



void Genetic::gen_node_random(Node<GANode>* pNode)
{
  int exist = 1;
  ncat++;
  string genecode;
  
  while(exist == 1)
  {
  gTree gtree(pInp);
  if(gtree.newStruct())
  genecode = gtree.printCode();
  exist = searchhash(genehash_all, genecode);
//  cout << exist << endl;  
  }
  
  cout << "******unique entry ******" << endl;
  addtohash(genehash_all,genecode);
  addtohash(genehash_pop,genecode);

  pNode->data.gene = genecode;

  pNode->data.comment = StringTools::int2str(ncat, 4, "0");
  
  cout << "come here check:" << pNode->data.comment<< endl;

  cout << "1randomly create a molecule: " << genecode << endl;

#if 0
  string cname = mmopt_folder + "ff.xyz" + pNode->data.comment;
  
  ofstream xyzfile(cname.c_str());
  gTree gtree2(pInp);
  gtree2.gene2struc(pNode->data.gene, xyzfile);
  xyzfile.close();

   ifstream infile(cname.c_str());
   pNode->data.get_xyz(infile);
   infile.close();
#endif

  pNode->data.flag = 0;
  pNode->next = NULL;

  pNext->addnode(pNode);
  optr_list.push_back(ncat);

//   get_fit4test(pNode);    //get MW info
//  Mopac::gen_inp(mopac_folder, pInp, pNode);
    searchhash_record(hash_record, pNode);
  
  return;
}

//if pNode->data.gene is unique
void Genetic::gen_node(Node<GANode>* pNode, int off_set)
 {
   ncat++;

    int old_label;

   addtohash(genehash_all,pNode->data.gene);
   addtohash(genehash_pop, pNode->data.gene);

  pNode->data.comment = StringTools::int2str(ncat, 4, "0");
/*hard-coded*/
#if 0
   old_label = ncat + off_set; 

   string tmp = StringTools::int2str(old_label,4,"0");
#endif

#if 1
   string cname = mmopt_folder + "ff.xyz" + pNode->data.comment;
//    string cname = mmopt_folder + "ff.xyz" + tmp;

#if 0
   ofstream xyzfile(cname.c_str());
   gTree gtree(pInp);
   gtree.gene2struc(pNode->data.gene, xyzfile);
   xyzfile.close();
#endif

   ifstream infile(cname.c_str());
   pNode->data.get_xyz(infile);
   
//   cout << ncat << " " << old_label << " " << pNode->data.natoms << endl;
   infile.close();
#endif
   pNode->data.flag = 0;
   pNode->next = NULL;

   pNext->addnode(pNode);
   optr_list.push_back(ncat);
//   Mopac::gen_inp(mopac_folder, pInp, pNode);

//   get_fit4test(pNode);    //get MW info
//    oled_score(pNode);
//   searchhash_record(hash_record, pNode);

   cout << "Done gen_node()" << endl << endl;
   return;
 }

void Genetic::read_record()
{
  string cwd = StringTools::getcwd_str();
  string fname = cwd + "/record0";

  ifstream record0(fname.c_str(),ios::in);
  string line;
  vector<string> tok_line;

  getline(record0,line);

  while(getline(record0,line))
  {
    Node<record>* pNode = new Node<record>();
    
    tok_line = StringTools::tokenize(line," ");
    
    pNode->data.gene = tok_line[0];
    pNode->data.HOMO = atof(tok_line[1].c_str());
    pNode->data.LUMO = atof(tok_line[2].c_str());
    pNode->data.MW = atof(tok_line[3].c_str());
   
    cout << pNode->data.gene << endl;
//    pNode->data.score = atof(tok_line[3].c_str());
//    cout << pNode->data.HOMO << " " << pNode->data.LUMO << endl;
    addtohash_record(hash_record, pNode);
  }
  
  cout << "Finish reading record file" << endl;
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


void Genetic::spawn()
{
//   int count = 0;
   int pool = pInp->inp_ga.pop_size;
//   vector<vector<int> > offsprings;
   int npairs = pool/2;
   int nextra = pool%2;    // if odd pop_size
   double random;
   bool ptrs_identic;

   vector<int> tmp;
   vector<int>::iterator it;
   if (ngen == 0)
   {
    #if 1
     for (int i=0; i< pool; i++)
     {
//       ncat++;
       Node<GANode>* pNode = new Node<GANode>();
       gen_node_random(pNode);     
    }
    #endif
/*enhanced-initiation*/
   #if 0
    enhanced_init_gen();
   #endif 

   }
   else
   {
     
     pNext->clear();
     cout << "Current length of pNext List: " << pNext->listlength << endl;
//     select_parents();
     gen_pairs(pool);
     for (int i=0; i < npairs; i++)
     {
//     offsprings.clear();
     ptrs_identic = get_ptrs(i);     
     cout << "Applying genetic operator to " << i+1 << "-th pair " << endl;
     cout << "parent1 : " << ptr1->data.gene << endl;
     cout << "parent2 : " << ptr2->data.gene << endl;
     cout << endl;

     Node<GANode>* pNode1 = new Node<GANode>();
     Node<GANode>* pNode2 = new Node<GANode>();
 
     if(ptrs_identic)
     {
       cout << "***parent1 and parent2 are identical!***" << endl;
       int exist = searchhash(genehash_pop, ptr1->data.gene);
       if (!exist)
       {
         cout << "They are not in the current population" << endl;
         pNode1 = cpy_node(ptr1);
         pNext->addnode(pNode1);   //clone one of them
          
         pNode2->data.gene = mutation(ptr2->data.gene); 
//       pNode2->data.gene = enhanced_mutation(ptr2->data.gene);    

         gen_node(pNode2);
       }

       else 
       {
         cout << "They are already in the Current population" << endl;
         pNode1->data.gene = mutation(ptr1->data.gene);
//       pNode1->data.gene = enhanced_mutation(ptr1->data.gene);
         gen_node(pNode1);
   
         pNode2->data.gene = mutation(ptr2->data.gene);
//       pNode2->data.gene = enhanced_mutation(ptr2->data.gene);    
         gen_node(pNode2);        
       }       
     }
     else
     {
     random = Utils::randomf(1.0);     
     cout << "Get random number: " << random << "; Compared with crossover rate " << pInp->inp_ga.cross_r << endl;
 
     if (random < pInp->inp_ga.cross_r)
     {
      cout << "********Crossover and Mutation:********** " << endl;
      bool done = false;
      int count = 0;
      while (count < 5)    //hard-coded!!!!!!!!
      {
      crossover(ptr1->data.gene, ptr2->data.gene);
      int exist1 = searchhash(genehash_all,newcode1);
      int exist2 = searchhash(genehash_all,newcode2);
      if ((!exist1) && (!exist2)) 
      {
       cout << "find unique offsprings" << endl; 
       done = true; 
       break;
      }
      count++;
      }

//      cout << "offspring 1: " << newcode1 << endl;
//      cout << "offspring 2: " << newcode2 << endl;

      if (!done)  
      {
         newcode1 = ptr1->data.gene;
         newcode2 = ptr2->data.gene;
      }
 
      cout << "offspring 1: " << newcode1 << endl;
      cout << "offspring 2: " << newcode2 << endl;


      double random;
      random = Utils::randomf(1.0);
      if (random < pInp->inp_ga.muta_r)   
      {
        cout << "ptr1 is under mutation" << endl;
        newcode1 = mutation(newcode1);
      }  
      random = Utils::randomf(1.0);
      if (random < pInp->inp_ga.muta_r)   
      {
         cout << "ptr2 is under mutation" << endl;
         newcode2 = mutation(newcode2);
      }


      int exist1 = searchhash(genehash_all, newcode1);
      int exist2 = searchhash(genehash_all, newcode2);
      
      if (exist1)  
      {
      cout << "newcode 1 is already in the list" << endl;
      newcode1 = mutation(newcode1);

//    newcode1 = enhanced_mutation(newcode1);

      }
      if (exist2)  
      {   
       cout << "newcode 2 is already in the list" << endl;
       newcode2 = mutation(newcode2);
//     newcode2 = enhanced_mutation(newcode2);
//
      }

      pNode1->data.gene = newcode1;
      pNode2->data.gene = newcode2;

       cout << "save new molecules to the list" << endl; 

       gen_node(pNode1);
       gen_node(pNode2);
      }  //end if
    else 
     {
       cout << "*********Clone********" << endl;
       double random;
       random = Utils::randomf(1.0);
       if (random < pInp->inp_ga.muta_r)  
       {
       cout << "prt1 is under mutation" << endl;
       newcode1 = mutation(ptr1->data.gene);
       pNode1->data.gene = newcode1;
       gen_node(pNode1);
       }
       else
       {
       int exist = searchhash(genehash_pop, ptr1->data.gene);
       if (!exist)
       {
         cout << "prt 1 not in the current population" << endl;
         pNode1 = cpy_node(ptr1);
         pNext->addnode(pNode1);   //clone one of them
       }
       else
       {
         cout << "prt1 already in the Current population" << endl;
         pNode1->data.gene = mutation(ptr1->data.gene);
//         pNode1->data.gene = enhanced_mutation(ptr1->data.gene);
         gen_node(pNode1);
       }
       }

       random = Utils::randomf(1.0);
       if (random < pInp->inp_ga.muta_r)  
       {
         cout << "prt2 is under mutation" << endl;
         newcode2 = mutation(ptr2->data.gene);
         pNode2->data.gene = newcode2;
         gen_node(pNode2);
       }
       else
       {
      int exist = searchhash(genehash_pop, ptr2->data.gene);
      if (!exist)
      {
        cout << "prt2 is not in the current population" << endl;
        pNode2 = cpy_node(ptr2);
        pNext->addnode(pNode2);   //clone one of them
      }
      else
      {
        cout << "prt2 is already in the Current population" << endl;
        pNode2->data.gene = mutation(ptr2->data.gene);
//      pNode2->data.gene = enhanced_mutation(ptr2->data.gene);
        gen_node(pNode2);
      }
       }
     } //else for ****Clone****
    }  // end of (prts_identic) else
   }   //end of loop
     if (nextra)
     { 
      Node<GANode>* pOdd = new Node<GANode>();
      ncat++;
      int flag = permutation[pool-1];
      Node<GANode>* pointer = ptr_array[select_list[flag]];
      pOdd->data.gene = mutation(pointer->data.gene);
     }
   } //end of else
}

//stochastic acceptance O(1) version

void Genetic::select_parents()
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

void Genetic::select_RWS()    //roulette Wheel selection
{
    int pool = pInp->inp_ga.pop_size;

    vector<double> wheel(pool+1);
    int i = 0;
    wheel[0] = 0.0;
    
    while(i < pool)
    {
        wheel[i+1] = prob_array[i] + wheel[i];
        cout << "wheel " << i+1 << " " << wheel[i+1] << endl;
        i++;
    }
    
    i = 0;

    while(i<pool)
    {
       double random = Utils::randomf(wheel.back());
       cout << random << endl;
       select_list[i] = select_binary_search(wheel, random);
       cout << "select " << select_list[i] << endl;
       i++;          
    }
}

int Genetic::select_binary_search(vector<double> & wheel, double random)
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

void Genetic::select_RWS_SUS()   //roulette Wheel + stocastic universal sampling
{
   int  pool = pInp->inp_ga.pop_size;
   
   vector<double> wheel(pool+1);
   int i = 0;
   wheel[0] = 0.0;
   while(i < pool)
   {
       wheel[i+1] = prob_array[i] + wheel[i];
       i++;
   }
 
   i = 0;

   double random = Utils::randomf(wheel.back());
   double stepsize = wheel.back()/pool;

   while (i < pool)
   {
      select_list[i] = select_binary_search(wheel, random);
      random += stepsize;
      if (random > wheel.back())
            random = fmod(random, wheel.back());
      i++;         
   }        
}

void Genetic::select_linear_rank_based(double rate) //the reproduction rate of worst individual
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

   while(i< pool)
   {
     cout << random << endl;
     select_list[i] = select_binary_search(wheel,random);
     cout << "select " << select_list[i] << endl;   
     random += stepsize;
     if (random > wheel.back())
        random = fmod(random,wheel.back());
     i++;
   }    
}

void Genetic::select_exp_rank_based(double base) // 0 < base < 1
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

   double stepsize = wheel.back()/pool;
   double random = Utils::randomf(wheel.back());
   while ( i < pool)
   {
     cout << random << endl;
     select_list[i] = select_binary_search(wheel,random);
     cout << "select " << select_list[i] << endl;   
     random += stepsize;
     if (random > wheel.back())
        random = fmod(random, wheel.back());
     i++;
   }
}

void Genetic::select_binary_tournament()
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
      select_list[i] = (prob_array[index1]>prob_array[index2])?index1:index2;
      i++;    
   }    
}

void Genetic::gen_pairs(int pool)
{ 
//  int pool = pInp->inpga.pop_size();  

//generate a random permutation in a finite set
  
   int flag;
//  int index;
   
   for (int i= pool; i>1; i--)
   {
     flag = Utils::randomi(i);
     swap(permutation[flag],permutation[i-1]);        //shuffle     
   }

   cout << "Selected parents:" << endl;
   for (int i=0; i< pool; i++)
   {
     cout << select_list[permutation[i]] << endl;
   }
   
   cout << "******************" << endl;
    
//permuation array save the one random permutation of array index of select_list
//2i and 2i+1 element in select_list would be one pair
}

//int n: the number of n-pairs
//save pointers for this pair to ptr1 & ptr2
//if ptr1 & ptrs2 are identical return true, otherwise return false
bool Genetic::get_ptrs(int n)
{
   int flag;
   int index;

   flag = permutation[2*n+0];
   index = select_list[flag];
   ptr1 = ptr_array[index];
   cout << "parent 1: " << ptr1 << endl; 

   flag = permutation[2*n+1];
   index = select_list[flag]; 
   ptr2 = ptr_array[index];
   cout << "parent 2: " << ptr2 << endl;

   if(ptr1 == ptr2)
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

Genetic::Genetic(readinp *p1)
{
  pInp = p1;  
  parents.resize(2);
  offsprings.resize(2);
  ptr_array.resize(pInp->inp_ga.pop_size);
  prob_array.resize(pInp->inp_ga.pop_size);
  select_list.resize(pInp->inp_ga.pop_size);
  permutation.resize(pInp->inp_ga.pop_size); 
  hashinit();
}

Genetic::~Genetic()
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
void Genetic::gen_all_perm()
{
 int nseeds = pInp->inp_sgen.n_seeds;
 int nfrags = pInp->inp_sgen.n_regs;
 int ncaps = pInp->inp_sgen.n_caps;
 string code;
// int count = 0;


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
        Node<GANode>* pNode = new Node<GANode>();
        pNode->data.gene = code;

        gen_node(pNode);

      }
    }
   }
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

void Genetic::array_reset()
{
   for (int i=0; i< pInp->inp_ga.pop_size; i++)
   {
     ptr_array[i] = NULL;
     prob_array[i] = -1.0;
     select_list[i] = -9999;
     permutation[i] = i;
   }
}

void Genetic::sort_and_select()
{
  int pool = pInp->inp_ga.pop_size;
//  xyznode* head = pNext->gethead();
//  pNext->save_score(score_list, MW_list);
  int i=0;
  double max;
  double sum = 0.0;
  Node<GANode>* p1 = pNext->gethead();
  max = p1->data.score;

  Node<GANode>* newhead = pNext->sortlist(p1);
  pNext->newhead(newhead);
  p1 = pNext->gethead();  

  array_reset();

  while (i != pool)
  {
     ptr_array[i] = p1;
     cout << p1 << " " << p1->data.score << endl;
     score_list.push_back(p1->data.score);
     sum = p1->data.score + sum;
     p1 = p1->next;
     i++;
   }
 
   p1 = pNext->gethead();
   i = 0;
   while (i != pool)
   {
 //    if (i==0)
     prob_array[i] = (p1->data.score)/sum;
     p1 = p1->next;
     i++;
 //    else
 //    prob_array[i] = (p1->score)/sum + prob_array[i-1];
   }
   select_parents();
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

void Genetic::after_sort(int ngen)
{
//    Node<GANode>* head = pNext->gethead();
//    Node<GANode>* newhead = pNext->sortlist(head);
//    pNext->newhead(newhead);
#if 1
    string filename = StringTools::int2str(ngen,4,"0");
    filename = "./Scratch/opted_list.xyz"+filename;
    pNext->print_xyzfile(filename);
#endif

//  cout << "Rename New Gen as Current Gen" << endl;
//  cout << endl << endl;

  pTmp = pNext;
  pNext = pCurrent;
  pCurrent = pTmp;     //swap pNext and pCurrent
  pTmp = NULL;
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

Node<GANode>* Genetic::cpy_node(Node<GANode>* pNode)
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

   addtohash(genehash_pop, pNew->data.gene);
  
   return pNew;
}

int Genetic::mopac_calc(string filedir, vector<int> & joblist,  LinkedList<Node<GANode> >* pList)
{
   bool success = true;
   Node<GANode>* pNode;
   pNode = pList->gethead();

#if 0
   for(int i = 0; i < pList->listlength; i++)
   {
    Mopac::gen_inp(filedir, pInp, pNode);
    pNode = pNode->next;
   }
   
   Mopac::mopac_genqsh(filedir, joblist);   //hard-coded for 0-gen case (joblist.size() = pop_size)
   success = Mopac::run_opt_jobs(filedir, joblist);
#endif
 
   if(success)
   { 
//     Mopac::hard_code_save_alltolist(filedir,pList);
   Mopac::save_all_into_list(filedir,pList);
//   return 1;
   } 
   
   pNode = pList->gethead();
   for(int i=0; i< pList->listlength;i++)
   {
      oled_score(pNode);
      pNode = pNode->next;
   } 

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

bool Genetic::clear_folder_restart(int ngen)
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
int Genetic::read_restart()
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

void Genetic::run_ga()
{ 

  cout << "******************************" << endl;
  cout << "*****  Starting GA runs  *****" << endl;
  cout << "************test2*************" << endl;

#if 0
   fragment_prop();
   init_feature_space();
   cout << "finish GA part" << endl;
   exit(-1); 
#endif


  LinkedList<Node<GANode> > List0;
  LinkedList<Node<GANode> > List1;

//  Fitness fitness(pInp);
//  fitness.readinp_f(pInp->fitness_f);

//   read_clusters(); 
   read_record();
   ngen = 0;

//  ngen = read_restart();
  
  if(ngen == 0)
  {
  string cwd = StringTools::getcwd_str();

//  string folder = "/export/zimmerman/yzhao/S-Generator/runs/all_perm/Scratch/ref/";
//  string folder = cwd + "/Scratch/ref/";
// int status = mkdir(folder.c_str(), S_IRWXU | S_IROTH | S_IXOTH);
//  pRef = new Node<GANode>();
//  read_ref(folder,pRef);  
//  fitness.calc_ref(folder, pRef);

   int status;
   string path = "./Scratch/mm";
//  status = mkdir(path.c_str(),S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
  }

  else
  {
  clear_folder_restart(ngen);
  }

  pCurrent = &List0;
  pNext = &List1;

  string ngen_str;

while(ngen <= pInp->inp_ga.ga_stop)
{
  ncat = 0;
   ngen_str = StringTools::int2str(ngen, 3, "0"); 
   string cwd = StringTools::getcwd_str();  

//   mopac_folder = "/export/zimmerman/yzhao/S-Generator/runs/all_perm/Scratch/mopac/000/";
   mopac_folder = cwd + "/Scratch/mopac/" + ngen_str + "/";


//   dftopt_folder = cwd + "/Work/opt/" + ngen_str + "/";
//    mmopt_folder = "/export/zimmerman/yzhao/S-Generator/runs/all_perm/Scratch/mm/000/";
   mmopt_folder = cwd + "/Scratch/mm/" + ngen_str + "/";
//   fitness_folder = cwd + "/Work/fitness/" + ngen_str + "/";
   int status;
    
//   status = mkdir(mmopt_folder.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
//   status = mkdir(mopac_folder.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);

//   cout << endl << endl;
   cout << "***************************************"<<endl;    
   cout << "The " << ngen << "-th Generation" << endl;   
   gen_all_perm();
//   spawn();   
   cout << "****************************************"<< endl << endl;
//   exit(-1);
   mopac_calc(mopac_folder, optr_list,  pNext);

//   pNext->print_record(ngen);

//   save_restart(pRef);
   sort_and_select();
   after_sort(ngen);

   pCurrent->print_record(ngen);

   cout << "The length of pNext List" << pNext->listlength << endl;
   cout << "The length of pCurrent List" << pCurrent->listlength << endl;
   ngen++;

    genehash_pop.clear();

    genehash_pop.resize(2*pInp->inp_ga.pop_size+1);
    optr_list.clear();
  } //end of while

//    print_all_score();


    List0.clear();
    List1.clear();

//    delete pRef;
}

