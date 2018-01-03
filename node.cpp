#include "node.h"
using namespace std;

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
