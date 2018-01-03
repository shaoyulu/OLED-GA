#include "utils.h"
#include "gcode.h"
#include <ctype.h>
#include <algorithm>
#include <math.h>

using namespace std;

gTree::gTree(readinp* pointer)
{
  pInp = pointer;  
  root = NULL;
}

gTree::~gTree(void)
{
  destoryTree();
}

//unique_ptr<OBMol> gTree::create_molptr()
//{
//  return unique_ptr<IOBMol>(new OBMol());
//}

int gTree::readlib_lite(int libn)
{
   string filename = StringTools::int2str(libn, 1, "0");
   filename += ".xyz";
//   cout << "Loading building block " << filename << endl;
   ifstream infile;
   
   filename = "./Scratch/lib/" + filename;
//   filename = pInp->inp_sgen.lib + "/basis/" + filename;
   infile.open(filename.c_str());

 //  bblock.nhandles = 0;

   if(!infile)
   {
     cout << "Cannot read " << filename << endl;
     exit(-1);
   }  

   string line;
   vector<string> tok_line;

   getline(infile, line);
   getline(infile,line);
   line = StringTools::newCleanString(line);
   tok_line = StringTools::tokenize(line, " ");
   
   int nbranches = tok_line.size();
 
   infile.close();
   return nbranches;
}

int gTree::readlib(data & bblock,int libn)
{
  int posP = 0;

  bblock.label = libn;
//  bblock.label = Utils::randomi2(pInp->inp_sgen.nblocks);
//  bblock.layer = 0;

  string filename = StringTools::int2str(bblock.label, 1, "0");
  filename += ".xyz";
//  cout << "Loading building block " << filename << endl;
  ifstream infile;
  
  filename = "./Scratch/lib/" + filename;
//  filename = pInp->inp_sgen.lib + "/basis/" + filename;
  infile.open(filename.c_str());

//  bblock.nhandles = 0;

  if(!infile)
  {
    cout << "Cannot read " << filename << endl;
    exit(-1);
  }

  string line;
  vector<string> tok_line;

  getline(infile, line);
  line = StringTools::newCleanString(line);
  bblock.natoms = atoi(line.c_str());
  getline(infile,line);
  line = StringTools::newCleanString(line);
  tok_line = StringTools::tokenize(line, " ");
  vector<string>::iterator it;
  for (it = tok_line.begin(); it != tok_line.end(); ++it)
  {
   int tmp = atoi((*it).c_str());
//   bblock.nhandles += tmp;
   bblock.xyzw.push_back(tmp);
  }
  
//  bblock.nhandles += 1;    //number of sockets and 1 plug
//  bblock.nbranches = int(bblock.xyzw.size());
//  bblock.name.resize(bblock.natoms);
//  bblock.coords.resize(3*bblock.natoms);
    bblock.label_h.resize(bblock.natoms);

//  bool find = false;

  for(int i =0; i < bblock.natoms; i++)
  {
    getline(infile,line);
    tok_line = StringTools::tokenize(line, " ");
    bblock.label_h[i] = tok_line[4];
//    bblock.name[i] = tok_line[0];
//    bblock.coords[3*i+0] = atof(tok_line[1].c_str());
//    bblock.coords[3*i+1] = atof(tok_line[2].c_str());
//    bblock.coords[3*i+2] = atof(tok_line[3].c_str()); 

    if (tok_line[4].find("PLG") != string::npos)
    {
      posP = i;
//      find == true;
//      cout << "find PLG position at " << i+1 << endl;
      bblock.label_h[i] = "X";
      
    }
  }
   infile.close();
   return posP;

}


int gTree::readelib(data & bblock, int elibn)
{
  int posP = 0;
  bblock.label = elibn;
//  bblock.label = Utils::randomi2(pInp->inp_sgen.nblocks_endp);
//  bblock.layer = 0;
  
  string filename = StringTools::int2str(bblock.label, 1, "0");
  filename += ".xyz";
//  cout << "Loading building block " << filename << endl;
  ifstream infile;
  
  filename = "./Scratch/lib/" + filename;
//  filename = pInp->inp_sgen.lib + "/endpoint/basis/" + filename;
  infile.open(filename.c_str());

//  bblock.nhandles = 0;

  if(!infile)
  {
    cout << "Cannot read " << filename << endl;
    exit(-1);
  }

  string line;
  vector<string> tok_line;
 
  getline(infile, line);
  line = StringTools::newCleanString(line);
  bblock.natoms = atoi(line.c_str());
  getline(infile,line);
  line = StringTools::newCleanString(line);
//  tok_line = StringTools::tokenize(line, " ");
//  vector<string>::iterator it;
//  for (it = tok_line.begin(); it != tok_line.end(); ++it)
//  {
//   int tmp = atoi((*it).c_str());
//   bblock.nhandles = 1;          //only 1 plug
   bblock.xyzw.push_back(0);
//  }

//  bblock.name.resize(bblock.natoms);
//  bblock.coords.resize(3*bblock.natoms);
  bblock.label_h.resize(bblock.natoms);

//  bool find = false;

  for(int i =0; i < bblock.natoms; i++)
  {

    getline(infile,line);
    tok_line = StringTools::tokenize(line, " ");
//    bblock.name[i] = tok_line[0];
    bblock.label_h[i] = tok_line[4];
//    bblock.coords[3*i+0] = atof(tok_line[1].c_str());
//    bblock.coords[3*i+1] = atof(tok_line[2].c_str());
//    bblock.coords[3*i+2] = atof(tok_line[3].c_str());

    if (tok_line[4].find("PLG") != string::npos)
    {
//       bblock.name[i] = "TBD";
       posP = i;
//       find = true;
//       cout << "find PLG position at " << i+1 << endl;
       bblock.label_h[i] = "X";      
    }
//    bblock.label_h[i] = "X";
  } 
   infile.close();
   return posP;
}
/*
//return docking position (socket) from parent
int gTree::getpos(data & dockto)
{
  int posS = 0;
  for (int i=0; i < dockto.natoms; i++)
  {
  if (dockto.label_h[i].find("X") == string::npos)
  {
    dockto.label_h[i] = "USED";
    posS = i;
    break;
  } 
  } 
  return posS; 
}
*/

//return docking position (socket) with index i from parent
int gTree::getpos(data & seed,int index)
{
  int pos=0;
  string tag = StringTools::int2str(index,0," ");
  for(int i = 0; i < seed.natoms; i++)
  {
    if (seed.label_h[i].find(tag) != string::npos)
    {
       seed.label_h[i] = "X";
       pos = i;
       break;
    }  
  }
  return pos; 
}
/*
void gTree::cpynode(data &newnode, const data & oldnode)
{
    newnode.label = oldnode.label;
    newnode.plug = oldnode.plug;
    newnode.layer = oldnode.layer;
    newnode.nhandles = oldnode.nhandles;
    newnode.natoms = oldnode.natoms;
    newnode.coords = oldnode.coords;
    newnode.label_h = oldnode.label_h;
    return;
}
*/

vector<data>gTree::getcArray(){

  return cArray;
}

int gTree::newStruct()
{

   struct data seed;

   int libn = Utils::randomi2(pInp->inp_sgen.n_seeds);     
   readlib(seed, libn); 
   seed.layer = 0;
//   seed.label = 0;
//   seed.xyzw = pInp->seed.xyzw;
//   seed.natoms = pInp->seed.natoms;
//   seed.name = pInp->seed.name;
//   seed.coords = pInp->seed.coords;
//   seed.label_h = pInp->seed.label;
   seed.c_index = 0;
   seed.acceptor = -1;
   seed.socket = 0;
   seed.plug = 0;
   seed.tag = 0;
//   seed.nhandles = pInp->seed.nhandles;

   cArray.push_back(seed);

   int nlayers;

  if(pInp->inp_sgen.nlayers == 0)
    {
    nlayers = Utils::randomir(pInp->inp_sgen.layer_min, pInp->inp_sgen.layer_max);
    }
  else
    {
    nlayers = pInp->inp_sgen.nlayers;
    }

//cout << "user defined number of layers:" << nlayers << endl;

    newbranch(seed, 1, nlayers);

   cout << "successfully generate cArray" << endl;
   
   for (int i=0; i< cArray.size(); i++)
   {
// cout << cArray[i].label << "\t" << cArray[i].tag << "\t" << cArray[i].plug << "\t" << cArray[i].acceptor << "\t" << cArray[i].layer << "\t" << cArray[i].xyzw.size()<< endl;  
   }
   
   return 1;      
}

void gTree::newbranch(data & parent, int i, int nlayers)
{
  while (parent.xyzw.size()!= 0)
 {  
  struct data bblock;
    
  if ((i < nlayers) && (parent.xyzw.size() !=0))
  {
   int libn = Utils::randomir(pInp->inp_sgen.n_seeds+1, pInp->inp_sgen.n_regs + pInp->inp_sgen.n_seeds);
///cout << pInp->inp_sgen.n_seeds << " &&  " << pInp->inp_sgen.n_regs << endl;
   bblock.plug = readlib(bblock, libn)+1;
   int tag = cArray[parent.c_index].xyzw.size() - parent.xyzw.size() + 1;
//  bblock.socket = getpos(parent,tag);
   bblock.tag = tag;
   bblock.acceptor = parent.c_index;
   bblock.layer = i;
   bblock.c_index = int(cArray.size());
   
   cArray.push_back(bblock);  
   parent.xyzw.pop_back();
   
   newbranch(bblock, i+1, nlayers);
   
//   newbranch(parent, i, nlayers);

  }
  else if ((i == nlayers) && (parent.xyzw.size() != 0))
  {
    int elibn = Utils::randomir(pInp->inp_sgen.n_seeds+pInp->inp_sgen.n_regs+1, pInp->inp_sgen.n_seeds + pInp->inp_sgen.n_regs + pInp->inp_sgen.n_caps);
    bblock.plug = readelib(bblock, elibn)+1;
    int tag = cArray[parent.c_index].xyzw.size() - parent.xyzw.size() + 1;
//    bblock.socket = getpos(parent, tag);
    bblock.tag = tag;
    bblock.acceptor = parent.c_index;
    bblock.c_index = int(cArray.size());
    bblock.layer = i;

    cArray.push_back(bblock);
    parent.xyzw.pop_back();
    
//    newbranch(parent, i, nlayers);
  }
 } 
  return;
}

#if 0
//turnon symm
void gTree::newbranch_symm(data & seed, int socket)
{
  int nlayers;

  if(pInp->inp_sgen.nlayers == 0)
    {
    nlayers = Utils::randomir(pInp->inp_sgen.layer_min, pInp->inp_sgen.layer_max);
    }
  else
    {
    nlayers = pInp->inp_sgen.nlayers;
    }

  cout << "user defined number of layers:" << nlayers << endl;

  vector<data> tmp;
  vector<int> l_exist(nlayers,-1);

  tmp.push_back(seed);
  
  int nblocks = pInp->inp_sgen.nblocks;
  int nblocks_endp = pInp->inp_sgen.nblocks_endp;

  int round = 0;
  while(!tmp.empty())
  {
    round++;
    cout << "round: " << round << endl;

    cout << tmp.size() << endl;

    cout << "the last node in the stack: " << tmp.back().label << endl;

    struct data node;
    node.layer = tmp.back().layer + 1;
    cout << "working on layer: " << node.layer << endl; 
    cout << l_exist[node.layer-1] << endl;

    if (node.layer == nlayers && l_exist[node.layer-1] == -1)
    {
      int elibn = Utils::randomir(nblocks+1, nblocks+nblocks_endp);
      cout << "-1 cap the structure with " << elibn << endl;
      l_exist[node.layer-1] = elibn;
      node.plug = readelib(node,elibn);
      node.nhandles--;
//      node.tag = tmp.back().nbranches - int(tmp.back().xyzw.size()) +1;
 //     tmp.back().xyzw.pop_back();
    }
    else if (node.layer == nlayers && l_exist[node.layer-1] != -1)
    {
      if (tmp.back().xyzw.size() != 1)
      {
         int elibn = Utils::randomir(nblocks+1, nblocks+nblocks_endp);
         node.plug = readelib(node,elibn);
         node.nhandles--;
//         node.tag = tmp.back().nbranches - int(tmp.back().xyzw.size()) +1;
         tmp.back().xyzw.pop_back();
      }
      else
      {      
      int elibn = l_exist[node.layer-1];
      cout << "0 cap the structure with " << elibn << endl;
      node.plug = readelib(node,elibn);
       node.nhandles--;
//      node.tag = tmp.back().nbranches - int(tmp.back().xyzw.size())+1;
      }  
    }
    else if (node.layer != nlayers && l_exist[node.layer-1] == -1)
    {
      int libn = Utils::randomi2(nblocks);
      l_exist[node.layer-1] = libn;
      node.plug = readlib(node, libn);
      node.nhandles--;
//      node.tag = tmp.back().nbranches - int(tmp.back().xyzw.size()) + 1;
//      tmp.back().xyzw.pop_back();
    }
    else if (node.layer != nlayers && l_exist[node.layer-1] != -1)
    { 
       if (tmp.back().xyzw.size() != 1)
       {
          int libn = Utils::randomi2(nblocks);
          node.plug = readlib(node,libn);
          node.nhandles--;
//          node.tag = tmp.back().nbranches - int(tmp.back().xyzw.size()) + 1;
          tmp.back().xyzw.pop_back();
       }
       else
       {
       int libn = l_exist[node.layer-1];
       node.plug = readlib(node, libn);
       node.nhandles--;
//       node.tag = tmp.back().nbranches - int(tmp.back().xyzw.size()) + 1;
       }
    }
    
    if (node.layer == 1)
    {
       node.socket = socket;
    }
    else node.socket = getpos(tmp.back());
  
    node.acceptor = tmp.back().c_index;

    tmp.back().nhandles--;
 
    cout << "current node: "<< node.label << "with " << node.nhandles << " rest" << endl;
    cout << "its parent: " << tmp.back().label << "with " << tmp.back().nhandles << " rest" << endl;

   if (tmp.back().nhandles == 0)
   {
      tmp.pop_back();
   }
   
   node.c_index = int(cArray.size());
   cArray.push_back(node);
  
   if (node.nhandles != 0)
   {
     tmp.push_back(node);
   }

   }  

  cout << "Done installing one branch" << endl;

  return;
}

#endif

#if 0
//turnoff symm
void gTree::newbranch(data & seed, int tag, int socket)

{
   int nlayers;

   if(pInp->inp_sgen.nlayers == 0)
    {
    nlayers = Utils::randomir(pInp->inp_sgen.layer_min, pInp->inp_sgen.layer_max);
    }
    else
    {
    nlayers = pInp->inp_sgen.nlayers;
    }

//    string lib = pInp->inp_sgen.lib;

    vector<data> tmp;

    tmp.push_back(seed);

    while(!tmp.empty())
    {
      struct data node;
 //     node.nhandles--;
      node.layer = tmp.back().layer+1;
      node.tag = tag;   /*belongs to i-th branch*/

      if (node.layer == nlayers)
      {
//          string dir = lib + "/endpoint/";
      int nblocks = pInp->inp_sgen.nblocks;
      int nblocks_endp = pInp->inp_sgen.nblocks_endp;
      int elibn = Utils::randomir(nblocks+1, nblocks+nblocks_endp+1);
      node.plug=readelib(node, elibn);
      node.nhandles--;
 // could only select a block from end-point lib
      }
      else
      {
 //   random select a block from lib and also a pos
//      string dir = pInp->seed.lib;
      int libn = Utils::randomi2(pInp->inp_sgen.nblocks);
      node.plug=readlib(node, libn);
      node.nhandles--;
      }
      if (node.layer == 1)
      {
         node.socket = socket;
      }
      else node.socket = getpos(tmp.back());

      node.acceptor = tmp.back().c_index;

      if (tmp.back().nhandles == 0)
      {
        tmp.pop_back();
      }
      
      cArray.push_back(node);
      node.c_index = int(cArray.size())-1;

      if (node.nhandles != 0)
      {
      tmp.push_back(node);
      }
   }
   cout << "Done installing one branch" << endl;
   return;
}
#endif

string gTree::printCode(vector<data> cArray){

//  cout << "generating gene based on cArray" << endl;

  vector<data>::iterator itd;
  string code;
  deque<int> stack;
  code = StringTools::int2str(cArray[0].label,1,"0");

  for (itd = cArray.begin()+1; itd !=cArray.end(); ++itd)
  {
    if(stack.empty())
    {
      stack.push_back((*itd).layer);
    }
    else
    {
       while ((!stack.empty()) && ((*itd).layer <= stack.back()))
         {
            code = code + ")";
            stack.pop_back();
         }
       stack.push_back((*itd).layer);
    }
    code = code + "(";
    code = code + StringTools::int2str((*itd).label,1,"0");
    code = code + "-";
    code = code + StringTools::int2str((*itd).tag,1,"0");
//  cout << code << endl;
  }
  while (!stack.empty())
  {
     code = code + ")";
     stack.pop_back();
  }
  return code;

}


string gTree::printCode()
{
  cout << "generating gene based on cArray" << endl;

  vector<data>::iterator itd;
  string code;
  deque<int> stack;
  code = StringTools::int2str(cArray[0].label,1,"0");  
 
  for (itd = cArray.begin()+1; itd !=cArray.end(); ++itd)
  { 
    if(stack.empty())
    {
      stack.push_back((*itd).layer);
    }
    else
    {
       while ((!stack.empty()) && ((*itd).layer <= stack.back()))
         {
            code = code + ")";
            stack.pop_back();   
         }
       stack.push_back((*itd).layer);
    }   
    code = code + "(";
    code = code + StringTools::int2str((*itd).label,1,"0");
    code = code + "-";
    code = code + StringTools::int2str((*itd).tag,1,"0");
//  cout << code << endl;   
  }
  while (!stack.empty())
  {
     code = code + ")";
     stack.pop_back();
  }   
  return code;
}


int gTree::readcode(string code, vector<data> tempC)
{

   cArray.clear();
   parser(code);
   
// cout << "cArray size "<< cArray.size()<<endl;
  for (size_t i=0; i< cArray.size(); i++)
   {
//     int j=1;
//     if (cArray[i].acceptor == 0)
//     {
//        cArray[i].tag = j;
//        j++;
//     }
      
     cArray[i].tag=tempC[i].tag;
     cArray[i].acceptor=tempC[i].acceptor;

//   cout << cArray[i].label << "\t" << cArray[i].tag << "\t" << cArray[i].plug << "\t" << cArray[i].acceptor << "\t" << cArray[i].layer << "\t" << cArray[i].xyzw.size() << endl;
   } 
 
   createTree();
   return 1;
}

int gTree::readcode(string code)
{

   parser(code);
   createTree();
   return 1;

}


void gTree::parser(string code)
{
   deque<char> q;
   vector<string> newcode;
   int layer = 0;
   string tmp;
   int max_layer=0;

//   cout << "reading given genetic representation:" << endl;
   string::iterator its;
   for (its=code.begin(); its != code.end(); ++its)
   {
//     cout << "parser : " << *its << endl;
     q.push_back(*its);
   }
//    cout << "q_size : " << q.size() << endl;
//   cout << endl;
   while(!q.empty())
   {
      if(isdigit(q.front()))
      {
       tmp.clear();
       do
       {
        tmp = tmp + string(1,q.front());
        q.pop_front();
//      cout << "tmp : " << tmp << endl;
       } while (isdigit(q.front()));
 
        newcode.push_back(tmp);
      }
      else
      {
        tmp = string(1,q.front());
        newcode.push_back(tmp);
        q.pop_front();
      }
   }
   
//genetic code is like A(B(E(J))(F)(G))(C(H(K)))(D(I))
//high sym code is A(B(E(J))(F)(G))  all branches are same
//B= XX-XX (lib#-tag#)

  vector<string>::iterator itv; 
  itv = newcode.begin();
  struct data seed;
    
  readlib(seed, atoi((*itv).c_str()));  
  seed.layer = 0;
  seed.socket = 0;
  seed.plug = 0;
//  seed.xyzw = pInp->seed.xyzw;
//  seed.natoms = pInp->seed.natoms;
//  seed.name = pInp->seed.name;
//  seed.coords = pInp->seed.coords;
//  seed.label_h = pInp->seed.label;
  seed.c_index = 0;
  seed.acceptor = -1;
  seed.tag = -1;
//  seed.nhandles = pInp->seed.nhandles;
  cArray.push_back(seed);

  nnodes = 1;

  for (itv= itv+1; itv!=newcode.end(); ++itv)
  {
//     cout << *itv << endl;

    if((*itv).find("(") != string::npos)
     {
       itv++;
       layer++;
       save_node(itv,layer); 
       if (layer > max_layer) max_layer = layer; 
       nnodes++;  
     }
     else if((*itv).find(")") != string::npos)
     {
       layer--;    
     }
  }
   nlayers = max_layer;

   return;
}
//root is like Axx (none socket and plug info)
//B: XX-XX-XX
//load all info into node (coords as well)
void gTree::save_node(vector<string>::iterator & itv, int layer)
{
//   cout << "saving node to cArray:" << endl;
   struct data node;
   node.label = atoi((*itv).c_str());
   itv = itv+2;          //skip "-"
   node.tag = atoi((*itv).c_str());
//   itv = itv+2; 
//   node.plug = atoi((*itv).c_str());
   node.layer = layer;
   if (node.label > pInp->inp_sgen.n_regs)
   node.plug = readelib(node, node.label)+1;
   else
   node.plug = readlib(node, node.label)+1;
   node.c_index = nnodes;
//   node.tag = -1;
   cArray.push_back(node);
}

void gTree::add_branch(TNode* & pNode, int tag, int nbranches, int layer)
{
  
  if (layer > nlayers) {pNode =  NULL; return;}
  if (tag > nbranches) {pNode =  NULL; return;}

  int n_branches;

  pNode = new struct TNode;
  pNode->tag = tag;
 
  if (layer == 0)
  {
    pNode->label = Utils::randomi2(pInp->inp_sgen.n_seeds);
    n_branches = readlib_lite(pNode->label);
  }
  else if (layer < nlayers) 
  {
  pNode->label = Utils::randomir(pInp->inp_sgen.n_seeds+1, pInp->inp_sgen.n_seeds+pInp->inp_sgen.n_regs);
  n_branches = readlib_lite(pNode->label);
  }
  else if (layer == nlayers)
  {
  pNode->label = Utils::randomir(pInp->inp_sgen.n_seeds+pInp->inp_sgen.n_regs+1, pInp->inp_sgen.n_seeds + pInp->inp_sgen.n_regs + pInp->inp_sgen.n_caps);
  n_branches = 0;
  }

//cout << "in_add_branch (layer) " << layer  << endl;
//cout << "in_add_branch (n_branches) " << n_branches  << endl;
//cout << "in_add_branch (nbranches) " << nbranches  << endl;
  add_branch(pNode->firstchild, 1, n_branches, layer+1);
  add_branch(pNode->nextsibling, tag+1, nbranches, layer);
  
  return;
}

TNode* gTree::Tnode(int i)
{
   struct TNode* pNode = new struct TNode;
   pNode->label = cArray[i].label;
   pNode->tag = cArray[i].tag;
   pNode->c_index = i;
   pNode->firstchild = NULL;
   pNode->nextsibling = NULL;
  
   return pNode;
}

TNode* gTree::create(int & i)
{
   struct TNode* pNode = Tnode(i);
  
   if((i!= nnodes-1) && (cArray[i+1].layer > cArray[i].layer)) 
   {
   cArray[i+1].acceptor = i;
   pNode->firstchild = create(++i);
   }
   if((i!= nnodes-1) && (cArray[i+1].layer == cArray[i].layer))
   {
   cArray[i+1].acceptor = cArray[i].acceptor; 
   pNode->nextsibling = create(++i);
   }
   return pNode;
}

void gTree::createTree()
{
    int i = 0;
    root = create(i);
 
//    cout << "done tree creation" << endl;
    return;
}

TNode* gTree::get_root()
{
   TNode* pNode;
   pNode = root;
   
   return pNode;
}

void gTree::tree2code(TNode* pNode, string & code)
{
  if(pNode == NULL) return;
  code = code + "(" + StringTools::int2str(pNode->label, 1, "0") + "-" + StringTools::int2str(pNode->tag, 1, "0");
  tree2code(pNode->firstchild, code);
  code = code + ")";
  tree2code(pNode->nextsibling, code);

  return;
}

string gTree::tree_to_code()
{
   string code;
   code = StringTools::int2str(root->label,1,"0");
   tree2code(root->firstchild, code);
   return code;
}

vector<int> gTree::serialize_tree(TNode* root, int height)
{
   if (root == NULL)  
   { 
      cout << "the tree representation is not exist" << endl;
      exit(-1);
   }

   vector<int> result;
   queue<TNode*> q;
   
   int total = pow(2.0, height-1) +1;
   result.push_back(root->label);

   q.push(root->firstchild);
   while (!q.empty() && result.size()!=total)
   {
      if (q.front() != NULL)
      result.push_back(q.front()->label);
      else result.push_back(-1);

      if (q.front() == NULL)
      {
        q.push(NULL);
        q.push(NULL);
      }
      else
      {
      if(q.front()->firstchild == NULL)
      {
         q.push(NULL);
      }
      else  
      q.push(q.front()->firstchild);
       
      if(q.front()->nextsibling == NULL)
      {
          q.push(NULL);
      }
      else
      q.push(q.front()->nextsibling);
      }
      q.pop();     
   }

   return result;   
}

vector<int> gTree::serialize(int height)
{
  return serialize_tree(this->root, height);
}

void gTree::destory(TNode* pNode)
{
   if(pNode == NULL)   return;
   destory(pNode->firstchild);
   destory(pNode->nextsibling); 
   delete pNode; 

   return;
}

void gTree::destoryTree()
{
   destory(root);
   root = NULL;

   return;
}

void gTree::getNode(TNode* pNode, bool & find, int index, TNode* & target)
{
   if (pNode == NULL || find == true)
   return;
   
   else
   {
     if (pNode->c_index == index)
     {
        find = true;
        target = pNode;
        return;
     }
     else
     {
       getNode(pNode->firstchild, find, index, target);
       getNode(pNode->nextsibling,find, index, target);
     }
   }
   return;
}

TNode* gTree::getNodeTree(int index)
{
  bool find = false;
  TNode* target = NULL;
  getNode(root, find, index, target);
  return target;
}


void gTree::getParent(TNode* pNode, bool & find, int index, TNode* & parent)
 {
//   cout << pNode << endl;

   if (pNode == NULL || find == true)
   return;

   else
   {
     if ((pNode->firstchild != NULL && pNode->firstchild->c_index == index) || (pNode->nextsibling != NULL) && (pNode->nextsibling->c_index == index))
     {
       find = true;
       parent = pNode;
       return;
     }
     else
     {
     getParent(pNode->firstchild, find, index, parent);
     getParent(pNode->nextsibling, find, index, parent);
     }
   }
 }

 TNode* gTree::getParentTree(int index)
 {
    if (index == root->c_index)
    return NULL;

    else
    {
      bool find = false;
      TNode* parent = NULL;
      getParent(root,find,index,parent);
//      cout << parent << endl;
      return parent;
    }
 }


int gTree::maxDepth()
{
   return max_depth(this->root);
}  

int gTree::max_depth(TNode* pnode)
{
   if (pnode == NULL) return 0;
   int left = max_depth(pnode->firstchild);
   int right = max_depth(pnode->nextsibling);
   
   return max(left,right)+1; 
}

//installation from tail to head
unique_ptr<OpenBabel::OBMol> gTree::install()
{
    vector<unique_ptr<OpenBabel::OBMol> > stack;
    vector<int> stack_index;
  
// cout << "mol installation check : " << endl; 
//r (size_t i=0; i< cArray.size(); i++)
// {
//     cout << cArray[i].label << "\t" << cArray[i].tag << "\t" << cArray[i].plug << "\t" << cArray[i].acceptor << "\t" << cArray[i].layer << "\t" << cArray[i].xyzw.size() << endl;
//

    int index = int (cArray.size())-1;

    int tag;

//   bblock.label = libn;

   string filename = StringTools::int2str(cArray[index].label, 1, "0");
   filename += ".xyz";
// cout << "Loading building block " << filename << endl;
// ifstream infile;

    string dir = "./Scratch/lib/";
 
//   if (cArray[index].label <= pInp->inp_sgen.n_regs)
//   dir = pInp->inp_sgen.lib;
//   else
//   dir = pInp->inp_sgen.lib + "/endpoint";

//   filename = dir + "/basis/" + filename;
   filename = dir + filename;  

   unique_ptr<OpenBabel::OBMol> mol = GetMol(filename);

//   ptr->init(cArray[index].name, cArray[index].coords);

//    ptr->print_xyz();
    
    stack.push_back(move(mol));   //now ptr->NULL

//    stack.back()->print_xyz();
  
    stack_index.push_back(index);

    int plug;
    int socket;
//    unique_ptr<OpenBabel::OBMol> molnew;
//    unique_ptr<OpenBabel::OBMol> moladd;

    while (!stack.empty())
    {
      index--;

//     cout << index << endl << endl;
     filename = StringTools::int2str(cArray[index].label, 1, "0");
     filename += ".xyz";

//    if (cArray[index].label <= pInp->inp_sgen.n_regs)
//    dir = pInp->inp_sgen.lib;
//    else
//    dir = pInp->inp_sgen.lib + "/endpoint";
      filename = "./Scratch/lib/" + filename; 
//    filename = dir + "/basis/" + filename;

    mol = GetMol(filename);

//    icptr = create_icptr(cArray[index].natoms);
//      icptr->init(cArray[index].name, cArray[index].coords);
 
      while ((!stack_index.empty()) && (cArray[stack_index.back()].acceptor == index))
      {
            tag = cArray[stack_index.back()].tag;
                 
            for (int k=0; k < cArray[index].natoms; k++)
            {
//            cout << k << ":" << cArray[index].label_h[k] << endl;
            if (cArray[index].label_h[k].find("X") == string::npos)
            {
             int flag = atoi(cArray[index].label_h[k].c_str());
//             cout << "flag = " << flag << " tag = " << tag << endl;
             if (flag == tag)
             {
              socket = k+1;  
              plug = cArray[stack_index.back()].plug;
//              cout << "add " << stack_index.back() << " to " << index << " at " << socket << endl << endl;
              merge_mol(mol, socket, stack.back(), plug);
              OpenBabel::OBConversion conv;
              conv.SetOutFormat("xyz");
//              conv.Write(mol.get(), &cout);
//              moladd.reset();
//              moladd = move(molnew);
//              icptr->print_xyz();
              }  //end of if
             } 
            }   // end of k
          stack_index.pop_back();
          stack.back().reset();     
          stack.pop_back();
//          stack.push_back(move(icptr));
//          stack.back()->print_xyz();
       }    // end of while stack_index

        stack.push_back(move(mol));
        stack_index.push_back(index);

     if (index == 0)
     {
       mol = move(stack.back());
       stack.pop_back();
     }
    }   // end of while stack

//    ptr->update_bonds();
//    ptr->ic_create_nobonds();
    cout << "Done installation" << endl;
    return mol;
}

//install from a tree
/*
void gTree::installT(TNode* pNode, ICoord & ic)
{
  if (pNode != NULL)
  {
     docking(ICoord & ic);
     install(pNode->firstchild, ICoord & ic);
     install(pNode->nextslibling, ICoord & ic);
  }
  return;
}
*/

#if 0
//install all branches to seed 
//fixit!!!!
void gTree::install_All(ICoord & icseed)
{

  ICoord icbr;
  ICoord icnew;
  vector<data>::iterator itc;
  itc = cArray.begin()+1;
  vector<data>::iterator itc2;
  itc2 = cArray.begin();

  icseed.init((*itc).natoms,(*itc).name, (*itc).coords);
  
  icseed.print_xyz();
  cout << "starting 3D structure assembly" << endl;

  cout << "number of different branches need to install: " << pInp->seed.nhandles << endl; 
  for (int i=0; i < pInp->seed.nbranches; i++)
  { 
    cout << "installing the " << i+1 << "-th unique branch to seed" << endl;
 
    if (i == 0)
    {
       cout << "first branch:" << endl;
       install_B(icbr, itc2);
       icbr.print_xyz();
    }   
    else
    {
      icbr.freemem(); 
      install_B(icbr, itc2);
    }
   
    int tag = i+1;
    int socket;
    int plug;
    
    for (int k=0; k<pInp->seed.xyzw[i];k++)
    {
     if (k == 0)
     {
     for (itc; itc !=cArray.end(); ++itc)
     {
       if ((*itc).acceptor == 0)
       {
          plug = (*itc).plug;
          socket = (*itc).socket;
          cArray[0].label_h[socket] = "X";
//          count << "found plug:" << plug+1 << endl;
          break;
       }       
      }
     }
     else
     {
       for (int j=0; j < cArray[0].natoms; j++)
       {
          if (cArray[0].label_h[j].find("X") == string::npos)
          {
            int flag = atoi(cArray[0].label_h[j].c_str());
            if (flag == tag) 
             {
                socket =j;
                cArray[0].label_h[j] = "X";
                break;
             }
          }
       }
     }
    cout << "found plug: " << plug+1 << endl;
    cout << "add to socket: " << socket+1 << " on the seed " << endl; 
          
         merge_ic(icseed, socket, icbr, plug, icnew);
         icseed.freemem();
//         icbr.freemem(); 
         icseed.duplicate(icnew);   
//         (*itc).label_h[j] = "X"; 
         icnew.freemem();      
     } 
       
   }
//     icbr.freemem();
  
   icseed.update_bonds();
   icseed.ic_create_nobonds();
   icbr.freemem();
//   icnew.freemem();
   cout << "Done installation" << endl;

   return;
}

//install a branch from cArray
void gTree::install_B(ICoord & ic, vector<data>::iterator & it )
{
 vector<int> addorder(cArray.size());  
      //store the cArray index follow adding order;
//installation follows order of cArray (preOrder Traversal of gTree)

  vector<vector<int> > newpos;       
// store new pos of each atom after installation by adding order 

//  vector<data>::iterator it = cArray.begin()+1;
  while (it!=cArray.end())
  {
       if ((*it).layer == 1)
       {
        ic.init((*it).natoms,(*it).name, (*it).coords);
        vector<int> tmp;
        for (int i=0; i< (*it).natoms; i++)
        {
        tmp.push_back(i);
        }
        newpos.push_back(tmp);
        addorder[(*it).c_index] = int(newpos.size())-1;

//       cout << "Done here" << endl;
       
         for(it=it+1; it!=cArray.end();++it)
         {
         if ((*it).layer > 1)
         {
         cout << "dock " << (*it).c_index << " to " << (*it).acceptor << endl;
         dock(it,ic,addorder,newpos);               
         }
         else break;
         }
         break;
       } 
    ++it;
  }
  cout << "assembly a branch " << endl;
}
#endif

#if 0
void gTree::merge1(int t1, int t2, const unique_ptr<ICoord> & ic1, const unique_ptr<ICoord> & ic2)
 {
//   printf(" single atom substitution of %s at %i  \n",ic2.anames[t2].c_str(),t1+1);

   int natoms1 = ic1->natoms;

   //icnew.reset(natomsn,ic1.anames,ic1.anumbers,ic1.coords);
   ic1->anames[t1] = ic2->anames[t2];

   return;
 }
#endif

#if 0
//merge ic1 and ic2 and save into ic3
unique_ptr<ICoord> gTree::merge_ic(const unique_ptr<ICoord> & ic1, int socket, const unique_ptr<ICoord> & ic2, int plug)
 {

//    cout << "ic1: " << ic1.get() << endl;
//    cout << "ic2: " << ic2.get() << endl;
    vector<string> l1nam;
    vector<double> l1xyz;
    vector<double> l2xyz;

    int t1;
    int t2;

   for (int i=0; i < ic1->natoms; i++)
   {
     l1nam.push_back(ic1->anames[i]);
     l1xyz.push_back(ic1->coords[3*i+0]);
     l1xyz.push_back(ic1->coords[3*i+1]);
     l1xyz.push_back(ic1->coords[3*i+2]);
   }


   for(int i=0; i < ic2->natoms; i++)
   {
     l2xyz.push_back(ic2->coords[3*i+0]);
     l2xyz.push_back(ic2->coords[3*i+1]);
     l2xyz.push_back(ic2->coords[3*i+2]);
   }
    if ((ic1->natoms == 2) || (ic2->natoms == 2))
    {
      if (ic1->natoms == 2)
      {
         merge1(plug, socket, ic2, ic1);
         l1nam.clear();
         l1xyz.clear();

         for (int i=0; i < ic2->natoms; i++)
         {
         l1nam.push_back(ic2->anames[i]);
         l1xyz.push_back(ic2->coords[3*i+0]);
         l1xyz.push_back(ic2->coords[3*i+1]);
         l1xyz.push_back(ic2->coords[3*i+2]);
         }
 //        return;
      }
      else if (ic2->natoms == 2)
      {
        merge1(socket, plug, ic1, ic2);
        l1nam.clear();
        l1xyz.clear();

        for (int i=0; i< ic1->natoms; i++)
        {
        l1nam.push_back(ic1->anames[i]);
        l1xyz.push_back(ic1->coords[3*i+0]);
        l1xyz.push_back(ic1->coords[3*i+1]);
        l1xyz.push_back(ic1->coords[3*i+2]);
        }
 //       return;
      }
    }    /*single atom installation*/
    else
    {
 //     vector<double> xyza(3*ic2.natoms);
 //     for (int i=0; i < 3*ic2.natoms; i++) xyza[i] = ic2.coords[i];

      t1 = -1;     /*position of the atom next to socket*/
      int found = 0;
      for (int i = 0; i < ic1->natoms; i++)
      {
        if (ic1->bond_exists(socket,i))
        {
        cout << " found atom bond to socket: " << i+1 << endl;
        t1 = i;
        break;
         }
      }
      if (t1 == -1)
        {
          cout << " cannot find atom bond to socket" << endl;
          exit(1);
        }


      t2 = -1; /*position of the atom next to plug*/
      found = 0;

      for (int i = 0; i < ic2->natoms; i++)
      {
      if (ic2->bond_exists(plug,i))
      {
        cout << " found atom bond to plug: " << i+1 << endl;
        t2 = i;
        break;
      }
      }

      double d0 = (ic1->getR(t1) + ic2->getR(t2))/2.4;
      cout << "set distance between two building blocks as " << d0 << endl;

 //     for (int i=0; i < 3*l1nat; i++) xyznew.push_back(ic1.coords[i]);

      align_to_x(ic1->natoms, t1, socket, l1xyz);
      align_to_x(ic2->natoms, plug, t2, l2xyz);

 //     for (int i=0; i < l1nat; i++) namenew.push_back(ic1.anames[i]);

      l1nam.at(socket) = ic2->anames[t2];
      for (int i=0; i < ic2->natoms; i++)
       {
          if ((i != t2) && (i!=plug))
        {
          l1nam.push_back(ic2->anames[i]);
        }
       }


      for (int i=0; i < ic2->natoms; i++)
      {
        if ((i != t2) && (i!=plug))
        {
          l2xyz.at(3*i+0) = l2xyz.at(3*i+0) - l2xyz.at(3*t2+0);
          l2xyz.at(3*i+1) = l2xyz.at(3*i+1) - l2xyz.at(3*t2+1);
          l2xyz.at(3*i+2) = l2xyz.at(3*i+2) - l2xyz.at(3*t2+2);
        }
      }

      l1xyz.at(3*socket+0) = d0;
      l1xyz.at(3*socket+1) = 0.0;
      l1xyz.at(3*socket+2) = 0.0;

      for (int i=0; i < ic2->natoms; i++)
      {
        if ((i != t2) && (i!= plug))
        {
          l1xyz.push_back(l2xyz.at(3*i+0)+d0);
          l1xyz.push_back(l2xyz.at(3*i+1));
          l1xyz.push_back(l2xyz.at(3*i+2));
        }
      }


 //     return;
    }   /*merge two ligands*/
    int n = int(l1nam.size());
   
//    cout << n << endl;

    unique_ptr<ICoord> ic3(create_icptr(n));
//    cout << ic3.get() << endl;
    ic3->initxyz(l1nam, l1xyz);
    ic3->bonds_union(socket, plug, t2, ic1, ic2);
 //   cout << "finish merging ligands" << endl;
 //   ic1.freemem();
 //   ic2.freemem();
    return ic3;
 }

#endif

void gTree::align_to_x(int natoms, int t1, int t2, vector<double> & xyz)
 {

   //printf(" aligning structure to x axis \n");

   //t1 --> 0,0,0
   //t2 --> a,0,0


   double* x1 = new double[3];
   x1[0] = xyz[3*t2+0] - xyz[3*t1+0];
   x1[1] = xyz[3*t2+1] - xyz[3*t1+1];
   x1[2] = xyz[3*t2+2] - xyz[3*t1+2];

   double n1 = sqrt(x1[0]*x1[0]+x1[1]*x1[1]+x1[2]*x1[2]);
   for (int i=0;i<3;i++) x1[i] = x1[i]/n1;

   //printf(" vector from t1 to t2: %6.4f %6.4f %6.4f \n",x1[0],x1[1],x1[2]);

   double* xyzn = new double[3*natoms];
   for (int i=0;i<3*natoms;i++) xyzn[i] = 0.;
   double* xyz1 = new double[3*natoms];
   for (int i=0;i<natoms;i++)
   {
     xyz1[3*i+0] = xyz[3*i+0] - xyz[3*t1+0];
     xyz1[3*i+1] = xyz[3*i+1] - xyz[3*t1+1];
     xyz1[3*i+2] = xyz[3*i+2] - xyz[3*t1+2];
   }
   //print_xyz_gen(natoms,anames,xyz1);
   //printf("b  %s %4.3f %4.3f %4.3f \n",anames[t1].c_str(),xyz1[3*t1+0],xyz1[3*t1+1],xyz1[    3*t1+2]);
   //printf("b  %s %4.3f %4.3f %4.3f \n",anames[t2].c_str(),xyz1[3*t2+0],xyz1[3*t2+1],xyz1[    3*t2+2]);

   double* angles = new double[3];
   double** rotm = new double*[3];
   rotm[0] = new double[3];
   rotm[1] = new double[3];
   rotm[2] = new double[3];

   double* u1 = new double[3];
  //align second atom to x axis

   x1[0] = xyz1[3*t2+0] - xyz1[3*t1+0];
   x1[1] = xyz1[3*t2+1] - xyz1[3*t1+1];
   x1[2] = xyz1[3*t2+2] - xyz1[3*t1+2];

 //  n1 = sqrt(x1[0]*x1[0]+x1[1]*x1[1]+x1[2]*x1[2]);
   n1 = sqrt(x1[0]*x1[0]+x1[2]*x1[2]);
   for (int i=0;i<3;i++) x1[i] = x1[i]/n1;

   u1[0] = x1[0];
   u1[1] = x1[1];
   u1[2] = x1[2];
   //printf(" u1: %4.3f %4.3f %4.3f \n",u1[0],u1[1],u1[2]);

   angles[0] = angles[1] = angles[2] = 0.;
   angles[1] = acos(u1[0]);
   if (u1[2]<0.) angles[1] = angles[1] * -1;
   //angles[2] = acos(u1[1]);
   //printf(" angle of %i %i to x axis: %4.3f \n",t1,t2,angles[1]/3.14159);

   Utils::get_rotation_matrix(rotm,angles);

   for (int i=0;i<3*natoms;i++) xyzn[i] = 0.;
   for (int i=0;i<natoms;i++)
   for (int j=0;j<3;j++)
   for (int k=0;k<3;k++)
     xyzn[3*i+j] += rotm[j][k]*xyz1[3*i+k];
   for (int i=0;i<3*natoms;i++)
     xyz1[i] = xyzn[i];

   //print_xyz_gen(natoms,anames,xyzn);
   //printf("1  %s %4.3f %4.3f %4.3f \n",anames[t1].c_str(),xyz1[3*t1+0],xyz1[3*t1+1],xyz1[    3*t1+2]);
   //printf("1  %s %4.3f %4.3f %4.3f \n",anames[t2].c_str(),xyz1[3*t2+0],xyz1[3*t2+1],xyz1[    3*t2+2]);


  //start second rotation
   x1[0] = xyz1[3*t2+0] - xyz1[3*t1+0];
   x1[1] = xyz1[3*t2+1] - xyz1[3*t1+1];
   x1[2] = xyz1[3*t2+2] - xyz1[3*t1+2];

 //  n1 = sqrt(x1[0]*x1[0]+x1[1]*x1[1]+x1[2]*x1[2]);
   n1 = sqrt(x1[0]*x1[0]+x1[1]*x1[1]);
   for (int i=0;i<3;i++) x1[i] = x1[i]/n1;

   u1[0] = x1[0];
   u1[1] = x1[1];
   u1[2] = x1[2];
   //printf(" u1: %4.3f %4.3f %4.3f \n",u1[0],u1[1],u1[2]);

   angles[0] = angles[1] = angles[2] = 0.;
   angles[2] = acos(u1[0]);
   if (u1[1]>0.) angles[2] = angles[2] * -1;
   //printf(" angle of %i %i to x axis: %4.3f \n",t1,t2,angles[2]/3.14159);

   Utils::get_rotation_matrix(rotm,angles);

   for (int i=0;i<3*natoms;i++) xyzn[i] = 0.;
   for (int i=0;i<natoms;i++)
   for (int j=0;j<3;j++)
   for (int k=0;k<3;k++)
     xyzn[3*i+j] += rotm[j][k]*xyz1[3*i+k];
   for (int i=0;i<3*natoms;i++)
     xyz1[i] = xyzn[i];

   //print_xyz_gen(natoms,anames,xyzn);
   //printf("2  %s %4.3f %4.3f %4.3f \n",anames[t1].c_str(),xyz1[3*t1+0],xyz1[3*t1+1],xyz1[    3*t1+2]);
   //printf("2  %s %4.3f %4.3f %4.3f \n",anames[t2].c_str(),xyz1[3*t2+0],xyz1[3*t2+1],xyz1[    3*t2+2]);

   //checking final geom
   double THRESH = 0.1;
   if (abs(xyzn[3*t1+0]) > THRESH
    || abs(xyzn[3*t1+1]) > THRESH
    || abs(xyzn[3*t1+2]) > THRESH
    || abs(xyzn[3*t2+1]) > THRESH
    || abs(xyzn[3*t2+2]) > THRESH)
   {
     printf("\n WARNING: align did not set geom to 0,0,0/a,0,0 \n");
     //align_to_x(natoms,t1,t2,xyzn,anames);
   }

   for (int i=0;i<3*natoms;i++)
     xyz[i] = xyzn[i];

   delete [] angles;
   delete [] rotm[0];
   delete [] rotm[1];
   delete [] rotm[2];
   delete [] rotm;

   delete [] xyzn;
   delete [] xyz1;
   delete [] x1;
   delete [] u1;

   return;
 }

#if 0
void gTree::dock(vector<data>::iterator it, ICoord & ic, 
                    vector<int> & addorder, vector<vector<int> > & newpos)
{

  if((*it).acceptor == 0)   //if the block dock to seed
    { 
      vector<int> tmp;
      ic.init((*it).natoms,(*it).name, (*it).coords);
      for (int i=0; i< (*it).natoms; i++)
      {
      tmp.push_back(i);
      }
      newpos.push_back(tmp);
      addorder[(*it).c_index] = int(newpos.size())-1;
    }
  else
  {
    ICoord icadd;
    ICoord icnew;
    icadd.init((*it).natoms, (*it).name, (*it).coords);
    icadd.print_xyz();
    vector<double> tmp_coords(3*(*it).natoms);
    
    int t1 = -1;
    for (int i=0; i< icadd.natoms; i++)
    {
       if (icadd.bond_exists(i, (*it).plug))
       {
          cout << " found atom " << i << " bond to plug " << (*it).plug << " in block: " << (*it).label << endl;
          t1 = i;
          break;  
       }
     }
    if (t1 == -1)
    {
      cout << " cannot find atom bond to given position " << endl;
      exit(1); 
    }

    align_to_x((*it).natoms, (*it).plug, t1, (*it).coords);

    for (int i=0; i< (*it).natoms; i++)
    {
       if ((i!=(*it).plug) && (i!=t1))
        { 
         tmp_coords[3*i+0] = (*it).coords[3*i+0] - (*it).coords[3*t1+0];
         tmp_coords[3*i+1] = (*it).coords[3*i+1] - (*it).coords[3*t1+1];
         tmp_coords[3*i+2] = (*it).coords[3*i+2] - (*it).coords[3*t1+2];
        }
    }

    int t2 = -1;
    int ns_pos = newpos[addorder[(*it).acceptor]][(*it).socket];  //new sokt position

    cout << "new socket is atom:" << ns_pos << endl;

//    ic.print_xyz();

    for (int i=0; i < ic.natoms; i++)
    { 
       if(ic.bond_exists(i,ns_pos))
        {
           cout << " found atom bond to socket " << ns_pos << ": "<< i << endl;
           t2 = i;
           break;
        }     
    }
     
    if (t2 == -1)
    {

       cout << " Cannot find atom bond to given position " << endl;
       exit(1);
    }
    
    double d0 = (ic.getR(t2) + icadd.getR(t1))/2.4;
    cout << "set bond length as " << d0 << endl;

    
    vector<string> Tmpname(ic.natoms);
    vector<double> Tmpxyz(3*(ic.natoms));
 
 //   ic.print_xyz();
 

    for (int i=0; i < ic.natoms; i++)
    {
    Tmpname[i] = ic.anames[i];
    Tmpxyz[3*i+0] = ic.coords[3*i+0];
    Tmpxyz[3*i+1] = ic.coords[3*i+1];
    Tmpxyz[3*i+2] = ic.coords[3*i+2];
    }


    align_to_x(ic.natoms, t2, ns_pos, Tmpxyz);

//    cout << "Done alignment" << endl;
 
    Tmpname[ns_pos] = (*it).name[t1];

    Tmpxyz[3*ns_pos+0] = d0;
    Tmpxyz[3*ns_pos+1] = 0;
    Tmpxyz[3*ns_pos+2] = 0;

    int count=ic.natoms;
    vector<int> tmppos;
    for (int i=0; i < icadd.natoms; i++)
    {
       if (i == (*it).plug)
          tmppos.push_back(-1);
       else if (i == t1)
       { 
        tmppos.push_back(ns_pos);          
       }
       else 
       {
          tmppos.push_back(count);
          Tmpname.push_back((*it).name[i]);
          Tmpxyz.push_back(tmp_coords[3*i+0]+d0);
          Tmpxyz.push_back(tmp_coords[3*i+1]);
          Tmpxyz.push_back(tmp_coords[3*i+2]);
          count++; 
       }                 
     }
     newpos.push_back(tmppos);
     addorder[(*it).c_index] = int(newpos.size())-1;
     cout << count << endl;

     icnew.initxyz(count, Tmpname, Tmpxyz);
     icnew.print_xyz();
     cout << endl;
     icnew.bonds_union(ns_pos, (*it).plug, t1, ic, icadd);  // revise this function!!!!
     
     ic.freemem();
     
     ic.initxyz(count, Tmpname, Tmpxyz);
     ic.copy_bonds(icnew);

     icnew.freemem();  
     icadd.freemem();
     }

     cout << "Current natoms = " << ic.natoms << endl;
     return;    
  }
#endif

//for test ONLY
void gTree::installation()
{
   if(newStruct())
   genecode = printCode();

   cout << "randomly create a molecule: " << genecode << endl;
   
   unique_ptr<OpenBabel::OBMol> molptr(install());

//   icptr->print_xyz();
   
//   ic.freemem();

   return;
}

void gTree::gene2struc(string code, ofstream &ostream)
{

    readcode(code);

    OpenBabel::OBConversion conv;

    conv.SetOutFormat("xyz");
    unique_ptr<OpenBabel::OBMol> mol(install());

    if(minimize(mol))
    {
    cout << "finish 1st optimization by MM:" << endl;
    conv.Write(mol.get());
    conv.Write(mol.get(), &ostream);
    }
    return;
}

void gTree::gene2struc2(string code, ofstream &ostream, vector<data> cArray2)
{

    cArray.clear();
    for (int i=0; i<(cArray2.size());i++){
    cArray.push_back(cArray2.at(i));}
    
//  cout << "in gene2struc2 : " << cArray2.size() << "  " << cArray.size() << endl;
    readcode(code,cArray);

    OpenBabel::OBConversion conv;

    conv.SetOutFormat("xyz");
    unique_ptr<OpenBabel::OBMol> mol(install());

    if(minimize(mol))
    {
    cout << "finish 1st optimization by MM:" << endl;
    conv.Write(mol.get());
    conv.Write(mol.get(), &ostream);
    }
    return;
}


//get coords & anames from a .xyz file
void gTree::get_xyz(ifstream &ostream)
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
   int natoms = atoi(line.c_str());
   getline(ostream,line);   

   cout << "get_xyz natoms =  " << natoms << endl;
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

void gTree::print_xyz(string filename)
{
   ofstream xyzfile;
   
   xyzfile.open(filename.c_str(), ios::app);
   xyzfile.setf(ios::fixed);
   xyzfile.setf(ios::left);
   xyzfile << setprecision(6);

   xyzfile << anames.size() << endl;
   xyzfile << genecode << endl;
  
   for (size_t i=0; i < anames.size(); i++)
   {
     xyzfile.width(2);
     xyzfile << std::left << anames[i] << "\t";
     xyzfile.width(9);
     xyzfile << std::right << coords[3*i+0] << "\t";
     xyzfile.width(9);
     xyzfile << std::right << coords[3*i+1] << "\t";
     xyzfile.width(9);
     xyzfile << std::right << coords[3*i+2] << endl;
   }
  xyzfile.close();
  return;
}

unique_ptr<OpenBabel::OBMol> gTree::GetMol(string filename)
{
//   // Create the OBMol object.
     unique_ptr<OpenBabel::OBMol> mol(new OpenBabel::OBMol);
//
//       // Create the OBConversion object.
     OpenBabel::OBConversion conv;
     OpenBabel::OBFormat *format = conv.FormatFromExt(filename.c_str());
     conv.SetOutFormat("xyz");
     if (!format || !conv.SetInFormat(format))
     {
     cout << "Could not find input format for file " << filename << endl;
     return mol;
     }
// Open the file.
//     ifstream ifs(filename.c_str());
//     if (!ifs)
//     {
//     cout << "Could not open " << filename << " for reading." << endl;
//     return mol;
//     }

//     cout << mol->NumAtoms() << endl;
// Read the molecule.
     if(!conv.ReadFile(mol.get(), filename))
     {
     cout << "Could not read molecule from file " << filename << endl;
     return mol;
     }
//     cout << mol->NumAtoms() << endl;
//     cout << "reading .xyz file: " << endl;

//     conv.Write(mol.get(), &cout);

     return mol;
}

//merge two fragment
 void gTree::merge_mol(unique_ptr<OpenBabel::OBMol> & mol1, int socket, unique_ptr<OpenBabel::OBMol> & mol2, int plug)
 {
   int n = mol1->NumAtoms();

 // Add second molecule
    (*mol1) += (*mol2);

    int i = plug + n;

 // Connect the fragments
    OpenBabel::OBBuilder::Connect(*mol1, socket, i);

   return;

 }

int gTree::minimize(unique_ptr<OpenBabel::OBMol> & mol)
{
//      OpenBabel::OBConversion conv;
//      conv.SetOutFormat("xyz");
      string ff = "MMFF94s";
//      unique_ptr<OpenBabel::OBMol> ptr(new OpenBabel::OBMol());
      OpenBabel::OBForceField* pFF = OpenBabel::OBForceField::FindForceField(ff.c_str());

//Make sure we have a valid pointer
     if (!pFF)  
        {
          cout << "cannot find ff: " << ff << endl;
          exit(-1);
        }

        pFF->SetLogFile(&cout);
        pFF->SetLogLevel(OBFF_LOGLVL_NONE);

     if (!pFF->Setup(*mol)) 
        {
        cout << "ERROR: could not setup force field." << endl;
        exit(-1);
        }

      pFF->SteepestDescent(50);

      pFF->ConjugateGradients(1000);

      pFF->GetCoordinates(*mol);

//      conv.Write(mol.get(), &cout);

      return 1;
}     
