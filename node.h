#ifndef NODE_H
#define NODE_H
#include <iostream>
#include <vector>
#include <fstream>
#include "stringtools.h"
using namespace std;

template<class T>
class Node
{
    public:
        Node();
        Node(const T& item, Node<T>* ptrnext = NULL);
        T data;
// access to the next node
        Node<T>* NextNode();
// list modification methods
        void InsertAfter(Node<T>* p);
        Node<T>* DeleteAfter();
        Node<T> * GetNode(const T& item, Node<T>* nextptr = NULL);
//    private:

        Node<T> * next;
                     };
struct record
{
    string gene;
    double HOMO;
    double LUMO;
    double MW;
};

struct GANode
{
//       int charge;
//       int multi;
       size_t natoms;
//       int ntemp;                  /*based on which template, ref=0*/
       vector<double>  coords;
       vector<string> anames;
       string    gene;      /*R code:  save as |V1_code|V2_code|...|*/
       string comment;            /*name of the molecule? use as filename?*/
       double  score;
       double score_ETL;
       double score_HTL;
       double HOMO;       /*Energy read after OPT*/
       double LUMO;       /*Energy for specie II*/    //hard-coded
       double MW;
       int flag;            /*flag = 0: ready to update; flag = -1: no need for update*/

       int multi = 1;
       int charge = 0;

       void init(int natoms)
       { coords.resize(3*natoms);
         anames.resize(natoms);
       }
    
 void get_xyz(ifstream &ostream)
 {
    if (!ostream)
    {
     cout << "Cannot read MM output file " << endl;
     exit(-1);
    }

    string line;
    vector<string> tok_line;

    getline(ostream,line);
    line = StringTools::newCleanString(line);
    natoms = atoi(line.c_str());
    getline(ostream,line);

    anames.resize(natoms);
    coords.resize(3*natoms);

    for (int i=0; i<natoms; i++)
    {
      getline(ostream, line);
      tok_line = StringTools::tokenize(line," ");
      anames.at(i) = tok_line[0];
      coords.at(3*i+0) = atof(tok_line[1].c_str());
      coords.at(3*i+1) = atof(tok_line[2].c_str());
      coords.at(3*i+2) = atof(tok_line[3].c_str());
    }
    return;
 }

};

template<class T>
Node<T>::Node()
{
// default constructor
// this is to allow us to create an object without any initialization

}

//  This constructor is just to set next pointer of a node and the data contained.
template<class T>
Node<T>::Node(const T& item,Node<T>* ptrnext)
{
   this->data = item;
   this->next = ptrnext;
}

template<class T>
Node<T>*Node<T>::NextNode()
{
   return this->next;
}

//  This methods inserts a node just after the node that the method belongs to
//  TO-DO: Consider a better implementatio template<class T>
template<class T>
void Node<T>::InsertAfter(Node<T> *p)
{
    p->next = this->next;
    this->next = p;
}

// Deletes the node from the list and returns the deleted node
template<class T>
Node<T>* Node<T>::DeleteAfter()
{
   Node<T>* tempNode = next;
   if(next != NULL)
   next = next->next;

   return tempNode;
}

template<class T>
Node<T> * GetNode(const T& item, Node<T>* nextptr = NULL)
{
    Node<T>* newnode; // Local ptr for new node
    newnode = new Node<T>(item,nextptr);
    if ( newnode == NULL)
    {
       cerr << "Memory allocation failed." << endl;
       exit(1);
    }
    return newnode;
}
#endif // NODE_H 
