#include <fstream>
#include <stdio.h>
#include <stdlib.h>

using namespace std;
#include "genetic.h"
#include "utils.h"

#define NUM_CLUSTERS 20    //5,10,20,40
#define CLUSTERING_METHOD "kmeans_20"

void Genetic::read_clusters() {

  int nfrags = pInp->inp_sgen.n_regs;
  int ncaps = pInp->inp_sgen.n_caps;

  for (int i=0; i < nfrags; i++)
   {
    string istr = StringTools::int2str(i+1,1,"0");
    for (int k=0; k < nfrags; k++)
     {
       string kstr = StringTools::int2str(k+1,1,"0"); 
       for (int j=0; j < ncaps; j++)
        {
            string jstr = StringTools::int2str(j+1+nfrags, 1, "0");
            string code = "0(" + istr + "-1(" + kstr + "-1(" + jstr + "-1)))";
            cout << code << endl;
            this->entire_space.push_back(code);
        }
     }
    }

    ifstream infile;
    infile.open("./inp/clustering.csv");
    
    if (!infile)
    {
        cout << "cannot open clustering result." << endl;
        exit(-1);
    }
    
    string line;
    vector<string> tok_line;
    
    getline(infile,line);
    tok_line = StringTools::tokenize(line,",");
    vector<string>::iterator it;    
    int col_num = 0;
    for (it = tok_line.begin(); it != tok_line.end(); it++)
    {
        if (*it == CLUSTERING_METHOD) break;  
        col_num++;
    }

    col_num++;    // take first "index" column into count
    cout << "read " << col_num << "-th column" << endl;
 
    for (int i = 0; i < NUM_CLUSTERS; i++)
    {
        vector<string>* ptr_v = new vector<string>;
        this->cluster_info.push_back(ptr_v);         
    }   

    int count = 0;
    while(getline(infile,line))
    {
        tok_line = StringTools::tokenize(line,",");
        int cluster_idx = StringTools::str2int(tok_line[col_num]);
        cluster_info[cluster_idx]->push_back(entire_space[count]);
//        cout << cluster_idx << " : " << entire_space[count] << endl; 
        count++;
    }

    for (int i=0; i < NUM_CLUSTERS; i++)
    {
        cout << "cluster " << i << " : " << cluster_info[i]->size() << endl;
    }    
    cout << "finish clustering info reading." << endl;
    return;    
}

void Genetic::enhanced_init_gen() {

    int pool = pInp->inp_ga.pop_size;
    int num_pick_each_cluster = pool/cluster_info.size();
    int exist;
    string genecode;

    for (int i=0; i< cluster_info.size(); i++)
    {
      for (int j =0 ; j < num_pick_each_cluster; j++)
      { 
         exist = 1;
         Node<GANode>* pNode = new Node<GANode>();
         this->ncat++;
         while (exist == 1)
         {   
         int random_idx = Utils::randomi(cluster_info[i]->size());
         cout << random_idx << endl;   
         genecode = cluster_info[i]->at(random_idx);
         exist = searchhash(genehash_all,genecode);
         }

        addtohash(genehash_all,genecode);
        addtohash(genehash_pop,genecode);
        
        cout << "cluster " << i << " : " << genecode << endl;

        pNode->data.gene = genecode;
        pNode->data.comment = StringTools::int2str(ncat, 4 ,"0");  
        pNode->data.flag = 0;
        pNode->next = NULL;
     
        this->pNext->addnode(pNode);
        searchhash_record(hash_record, pNode);
      }  
    }
}

string Genetic::enhanced_mutation(string genecode){

    string newcode = genecode;
    int found = 0;
    int cluster;

    for (int i =0; i < cluster_info.size() && found != 1; i++)
    { 
       int size_cluster = cluster_info[i]->size();
       for (int j =0; j < size_cluster; j++)
       {  
       if (cluster_info[i]->at(j) == newcode)
            { 
                found = 1;
                cluster = i;
                cout << "found " << newcode << " in cluster " << cluster << endl;    
                break;
            }      
       }   
    }
    
    if (found == 0)
    {
        cout << "Something wrong with genecode lookup module." << endl;
        exit(-1);
    }
    int exist = 1;
    int count = 0;
    found = 0;
    int cluster_size = cluster_info[cluster]->size();
    while (exist && count < cluster_size)    
    {
    count++;
    int random_idx = Utils::randomi(cluster_size);
    newcode = cluster_info[cluster]->at(random_idx);
    exist = searchhash(genehash_all, newcode);

    if (exist) 
    cout << "in the list already, mutate again" << endl;
    else
    found = 1;
    }
    
    if (found == 0)   
//  newcode = mutation(genecode);

    cout << endl;
    return newcode;
}

