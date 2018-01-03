#ifndef GENETIC_H
#define GENETIC_H

#include "read_inp.h"
#include "node.h"
#include "linkedlist.h"
#include "gcode.h"
#include <memory>
#include "fragment.h"
//#include "fitness.h"

class Genetic {

private:

//void all_candidates() // for enhanced initiation and enhanced mutation

vector<string> entire_space;         // all possible solutions
vector<vector<string>* > cluster_info;

/*hard-coded enhanced initiation and mutation*/
void read_clusters();
void enhanced_init_gen();
string enhanced_mutation(string genecode);

/*hard-coded enhance initiation and mutation*/

vector<string> all_genes_origin;
vector<vector<int> > all_serial_tree;   // serilazed binary tree in complete binary tree form

vector<Fragment*> pFrags;

void fragment_prop();
vector<data> returncArray();

void init_feature_space();

//double max_prob;

string mmopt_folder;

string mopac_folder;

string dftopt_folder;

string fitness_folder;

vector<int> optr_list;

vector<double> score_list;

LinkedList<Node<GANode> >* pCurrent;

LinkedList<Node<GANode> >* pNext;

LinkedList<Node<GANode> >* pTmp;   //use it when sawp pCurrent and pNext;

readinp* pInp;

Node<GANode>* pRef;

//pointers to a pair
Node<GANode>* ptr1;
Node<GANode>* ptr2;

vector<vector<int> > parents;
vector<vector<int> > offsprings;

//newcode after crossover
string newcode1;
string newcode2;

int ncat;

int ngen;

vector<unique_ptr<LinkedList<Node<record> > > > hash_record;

vector<string> genelist;
vector<unique_ptr<LinkedList<Node<string> > > > genehash_all;     //all previous gene strings
vector<unique_ptr<LinkedList<Node<string> > > > genehash_pop;     //gene strings in current pop

void hashinit();                  //initiate size of hash
//int savegene_hash(string str);

int addtohash(vector<unique_ptr<LinkedList<Node<string> > > > & hashtable, string & str);

int addtohash_record(vector<unique_ptr<LinkedList<Node<record> > > > & hash_record, Node<record>* pNode);

int searchhash(vector<unique_ptr<LinkedList<Node<string> > > > & hashtable, string & str);

void searchhash_record(vector<unique_ptr<LinkedList<Node<record> > > > & hash_record, Node<GANode>* pNode);

long BKDRHash(string str);                //Hash function

//vector<vector<int> > genelist;
vector<vector<int> > gene_in_pop; 

//vector of pointers to each node in the LinkedList
vector<Node<GANode>* >  ptr_array;   

//vector of probilities of selecting each node in the LinkedList
vector<double> prob_array;

//save selected parents
vector<int> select_list;

//save a random permutation, 2i and 2i+1 would be one pair
vector<int> permutation;   

void gen_pairs(int pool);

//encode new cats, save elements, coords and etc. into *pNode and add to xyzlist
//int* gen_frezlist(xyznode* pNode);
//void gen_node(xyznode* pNode);

void gen_node_random(Node<GANode>* pNode);

//with genecode in pNode->data.gene
void gen_node(Node<GANode>* pNode, int off_set = 0);
//make a replica
Node<GANode>* cpy_node(Node<GANode>* pNode);

void after_sort(int ngen);
//save index_array and prob_array
void sort_and_select();

//select a Template
int Tempcode();

//select a R group for current variable
int Rcode();

bool check_unique(vector<vector<int> > & target_list, int tempcode, vector<int> & rcodes);
bool check_unique2(vector<vector<int> > & target_list, vector<int> & offspring); 

//void gencode_random(xyznode* pNode);

//for test only, use molecular weight as score
void get_fit4test(Node<GANode>* pNode);
vector<double> MW_list;

void oled_score(Node<GANode>* pNode);
double oled_fit(Node<record>* pNode);

void spawn();

void gen_all_perm();

void select_parents();
void select_RWS();    //roulette Wheel selection
int select_binary_search(vector<double> & wheel, double prob);
void select_RWS_SUS();   // roulette Wheel + stocastic universal sampling
void select_binary_tournament();
void select_linear_rank_based(double rate);
void select_exp_rank_based(double base);

bool get_ptrs(int n);

void swap_ptr(TNode* & ptr1, TNode* & ptr2);

void crossover(string & gene1, string & gene2);

string mutation(string genecode);

//bool mutate(vector<int> & offspring);

//bool mutate2(xyznode* pCurrent, xyznode* pKid);

//void force_mutate(xyznode* pCurrent, xyznode* pKid);

//bool splice_mutate();

void save_genelist(int tempcode, vector<int> rcodes);

void save_gene_in_pop(int tempcode, vector<int> rcodes);

void array_reset();

void print_all_score();

int mopac_calc(string filedir, vector<int> & joblist,  LinkedList<Node<GANode> >* pList);

//int dft_opt(string filedir, vector<int> & joblist, LinkedList* pList);

bool save_restart(Node<GANode>* pRef);

int read_restart();  /*return Current Generation +1*/

void read_record();

bool clear_folder_restart(int ngen);

//void check_connect(LinkedList* pList);

//int get_fit(string filedir, vector<int> & list, Fitness & fitness);

void read_ref(string filedir, Node<GANode>* pRef);

public:

Genetic(readinp* p1);
~Genetic();

void run_ga();
};

#endif
