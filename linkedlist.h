#ifndef LINKEDLIST_H
#define LINKEDLIST_H

 #include "linkedlist.h"
 #include "stringtools.h"
 #include "node.h"
 #include <iostream>
 #include <fstream>
 #include <stdio.h>

using namespace std;

template<class T>
class LinkedList{

private:
    T* tail;
   
    T* head;
   
public:

   LinkedList();

  ~LinkedList();
   
   void clear(); 

   int listlength;  

//return the head of linkedlist 
   T* gethead();

//assign newhead (after sort) to class 
   void newhead(T* pNode);

// add node at tail
   void addnode(T* pNode);      

// delete node at any position, with the pointer which point to the node before the deleting position.
   void delnode(T* pPrev);

   int findnode(string str);
//keep node if score is in the range, otherwise delete it.
   void select(double target);

//merge two sorted list (a & b) 
  T* mergelist(T* a, T* b );
  T* mergelist_ETL(T* a, T* b );
  T* mergelist_HTL(T* a, T* b );

//divide list
  T* split(T* source);
 
//sort selected list, return head pointer.
  T* sortlist(T* a);
  T* sortlist_ETL(T* a);
  T* sortlist_HTL(T* a);

 void print_xyzfile(string fname);

 void print_record(int ngen);

 T* getnode(string str);  
};

template<class T>
T* LinkedList<T>::getnode(string str)
{
  T* ptr = head;
  while (ptr != NULL)
  {
     cout << ptr << "  " << ptr->data.gene << endl;
     if (ptr->data.gene == str)
      {
       break;
      }
     else
      ptr = ptr->next;
  }
  return ptr;
}

 template<class T>
 LinkedList<T>::LinkedList(void)
 {
       head = NULL;
       tail = NULL;    //This is the constructor
       listlength = 0;
 }

 template<class T>
 LinkedList<T>::~LinkedList(void)
 {
   clear();
 }

 template<class T>
 //clear all nodes
 void LinkedList<T>::clear()
 {
    T* pDel = head;

    while(listlength > 0)
    {
     head = head->next;
     delete pDel;

     pDel = head;
     listlength--;
    }

 //   delete pDel;
 //   listlength--;
    head = tail = NULL;
 }

 template<class T>
 //return the head of linkedlist
 T* LinkedList<T>::gethead()
 {
       T* pNode;
       pNode = head;

       return pNode;
 }

 template<class T>
 void LinkedList<T>::newhead(T* pNode)
 {
    head = pNode;
 }
 template<class T>
 void LinkedList<T>::print_record(int ngen)
 {
   T* pNode;
   pNode = head;
   int count = 0;

   ofstream record;
   string record_fname = "record";

   record.open(record_fname.c_str(), ios::app);
   record.setf(ios::fixed);
   record.setf(ios::left);
   record << setprecision(6);

   if (ngen == 0)
   {
   record.width(20);
   record << std::left << "gene" << "     ";
//   record.width(9);
//   record << std::right << "score" << "\t";
 //  record.width(4);
 //  record << std::right << "No." <<  "  ";
   record.width(9);
   record << std::right << "HOMO (eV)" << "     ";
   record.width(9);
   record << std::right << "LUMO (eV)" << "     ";
   record.width(9); 
//   record << std::right << "MW " << "     "; 
//   record.width(9);
   record << std::right << "score" << endl;
 //  record.width(12);
 //  record << std::right << "Energy-I" << "  ";
 //  record.width(12);
 //  record << std::right << "Energy-II" << "  ";
 //  record.width(12);
 //  record << std::right << "Energy-III" << endl;
   }

   while(count < listlength)
   {
 //  string tmp;
 //  tmp = StringTools::int2str(int(pNode->ntemp),1," ");
  // record << std::left << pNode->ntemp;
 //  for (vector<int>::iterator it = pNode->gene.begin(); it != pNode->gene.end(); ++it)
 //  {
 //    tmp = tmp + "-";
 //    tmp = tmp + StringTools::int2str(*it, 1, " ");
 //    record << "-" << *it;
 //  }

 //  record.width(10);
 //  record << std::left << tmp;
 //  record << "\t";
 //  record << "\t";

   record.width(20);
   record << std::left << pNode->data.gene << "     ";
//   record.width(9);
//   record << std::right << pNode->data.score << "\t ";
   record.width(9);
   record << std::right << pNode->data.HOMO << "     ";
   record.width(9);
   record << std::right << pNode->data.LUMO << "     ";
//   record.width(9);
//   record << std::right << pNode->data.MW << "     ";  
   record.width(9);
   if(pNode!=NULL){
   record << std::right << pNode->data.score << endl;
   }else{
   record << std::right << "99.999" << endl;}
 //  record.width(12);
 //  record << std::right << pNode->energy << "  ";
 //  record.width(12);
 //  record << std::right << pNode->energy2 << "  ";
 //  record.width (12);
 //  record << std::right << pNode->energy3 << endl;

   pNode = pNode->next;
   count++;
   }

   record.close();
 }

 template <class T>
 void LinkedList<T>::print_xyzfile(string fname)
 {
    T* pNode;
    pNode = head;
    int count = 0;

    ofstream outfile;
    outfile.open(fname.c_str());
    outfile.setf(ios::fixed);
    outfile.setf(ios::left);
    outfile << setprecision(6);

    while (count < listlength)
    {
      outfile <<  pNode->data.natoms << endl;
      outfile <<  pNode->data.gene << endl;
 //     outfile <<  "Genecode: " <<pNode->ntemp;
 //     size_t s = pNode->gene.size();
 //     for (size_t i=0; i < s; i++)
 //     {
 //     outfile << " " << pNode->gene[i];
 //     }
 //     outfile << endl;
      for (size_t i=0; i< pNode->data.natoms; i++)
      {
      outfile.width(2);
      outfile << std::left << pNode->data.anames[i] << "\t ";
      outfile.width(9);
      outfile << std::right << pNode->data.coords[3*i+0] << "\t";
      outfile.width(9);
      outfile << std::right << pNode->data.coords[3*i+1] << "\t";
      outfile.width(9);
      outfile << std::right << pNode->data.coords[3*i+2] << "\t" << endl;
      }

      pNode = pNode->next;
      count++;
    }

 }

 template<class T>
 // add node at tail
 void LinkedList<T>::addnode(T* pNode)
        {
 //         xyznode* n = pNode;
 //        n->charge = val->charge;
 //       n->multi  = val->multi;
 //        n->natoms  = val->natoms;
 //        n->ntemp = val->ntemp;
 //        n->coords  = val->coords;
 //        n->element = val->element;
 //        n->comment = val->comment;
 //        n->gene = val->gene;
 //        n->energy = val->energy;
 //        n->score = val->score;

         if (head == NULL)
          {  head = pNode;
             tail = pNode;
             tail->next = NULL;
          }
         else
         {
            tail->next = pNode;  //make it point to the next node
            tail= tail->next;     //reset tail point at the new node
            tail->next = NULL;
          }
      listlength++;
   }
 template<class T>
 int LinkedList<T>::findnode(string str)
 //find a node in a list, return 1; otherwise return 0;
 {
    T* ptr = head;
    int flag = 0;
    while (ptr != NULL)
    {
       if (ptr->data == str)
        {
         flag =1;
         break;
        }
       else
        ptr = ptr->next;
    }
    return flag;
 }

 template<class T>
 // delete node at any position, with the pointer which point to the node before the deleti    ng position.
 void LinkedList<T>::delnode(T* pPrev)
    {
     T* pDel;
         pDel = pPrev->next;

    if (pPrev->next == head)
      {
         head = head->next;
         pPrev->next = pDel->next;
         delete pDel;
         pDel = NULL;
      }
    else if (pPrev->next == tail)
      {
         tail = pPrev;
         delete pDel;
         pDel = NULL;
      }
    else
      {
        pPrev->next  = pDel->next;
        delete pDel;
        pDel = NULL;
      }

     --listlength;

    return;
    }
 template<class T>
 T* LinkedList<T>::mergelist(T* a, T* b )
  {
    if (a == NULL)
      return b;
    else if (b == NULL)
      return a;
    T* result;
    if (a->data.score >= b->data.score)
    {
      result = a;
      result->next = mergelist(a->next, b);
    }
    else
    {
      result = b;
      result->next = mergelist(a, b->next);
    }
    return result;
 }
 template<class T>
 T* LinkedList<T>::mergelist_ETL(T* a, T* b )
  {
    if (a == NULL)
      return b;
    else if (b == NULL)
      return a;
    T* result;
    if (a->data.score_ETL >= b->data.score_ETL)
    {
      result = a;
      result->next = mergelist_ETL(a->next, b);
    }
    else
    {
      result = b;
      result->next = mergelist_ETL(a, b->next);
    }
    return result;
 }
template<class T>
 T* LinkedList<T>::mergelist_HTL(T* a, T* b )
  {
    if (a == NULL)
      return b;
    else if (b == NULL)
      return a;
    T* result;
    if (a->data.score_HTL >= b->data.score_HTL)
    {
      result = a;
      result->next = mergelist_HTL(a->next, b);
    }
    else
    {
      result = b;
      result->next = mergelist_HTL(a, b->next);
    }
    return result;
 }

 
 template<class T>
 //divide list
  T* LinkedList<T>::split(T* source)

 {
    T* fast;
    T* slow;

  if ((source == NULL) || (source->next == NULL))
  {
    return source;
  }
  else
  {
    slow = source;
    fast = slow->next;

    while(fast != NULL)
    {
      fast = fast->next;
      if(fast != NULL)
       {
         slow = slow->next;
         fast = fast->next;
       }
    }
  }
  return slow;
 }

   template<class T>
 T* LinkedList<T>::sortlist(T* a)
   {
    if ((a == NULL) || (a->next ==NULL))
      return a;

     T* mid;
     T* b = NULL;

    mid = split(a);

    if (mid != NULL)
    {
       b = mid->next;
       mid->next = NULL;
    }

    a = sortlist(a);
    b = sortlist(b);

    a = mergelist(a,b);
    return a;
   }
   template<class T>
 T* LinkedList<T>::sortlist_ETL(T* a)
   {
    if ((a == NULL) || (a->next ==NULL))
      return a;

     T* mid;
     T* b = NULL;

    mid = split(a);

    if (mid != NULL)
    {
       b = mid->next;
       mid->next = NULL;
    }

    a = sortlist_ETL(a);
    b = sortlist_ETL(b);

    a = mergelist_ETL(a,b);
    return a;
   }

  template<class T>
 T* LinkedList<T>::sortlist_HTL(T* a)
   {
    if ((a == NULL) || (a->next ==NULL))
      return a;

     T* mid;
     T* b = NULL;

    mid = split(a);

    if (mid != NULL)
    {
       b = mid->next;
       mid->next = NULL;
    }

    a = sortlist_HTL(a);
    b = sortlist_HTL(b);

    a = mergelist_HTL(a,b);
    return a;
   }
#endif   
   


  
