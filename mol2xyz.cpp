#include "mol2xyz.h"
#include <iostream>
#include <fstream>
#include <stdio.h>
#include <stdlib.h>
#include <memory>

using namespace std;

unique_ptr<OpenBabel::OBMol> MOL2XYZ::readXYZ(string filename)
{
//cout << "reading user provided XYZ file..." << endl;
  ifstream infile;
  infile.open(filename.c_str());

  if(!infile)
  {
  cout << "can't open the file" << endl;
  exit(-1);
  }

  string line;
  vector<string> tok_line;
  
//  getline(infile, line);
//  getline(infile, line);
 // getline(infile, line);

  getline(infile, line);
//  line = StringTools::newCleanString(line);

//  tok_line = StringTools::tokenize(line," ");

//  int natoms = atoi(tok_line[0].c_str());
//  vector<int> add_H(natoms,1);

//  int nbonds = atoi(tok_line[1].c_str());

  getline(infile, line);

  int idx = 0;
 
  while(getline(infile, line))
  {
   line = StringTools::newCleanString(line);
    tok_line = StringTools::tokenize(line, " ");
   idx++;

   if(tok_line[4].find("PLG") != string::npos)
   {
   index.push_back(idx);
   label.push_back("PLG"); 
   }
   else if (tok_line[4].find("1") != string::npos)
    {
      index.push_back(idx);
      label.push_back("1");
    }
    else if (tok_line[4].find("2") != string::npos)
    {
        index.push_back(idx);
        label.push_back("2");
    }
    else if (tok_line[4].find("3") != string::npos)
    {
        index.push_back(idx);
        label.push_back("3");
    }
    else if (tok_line[4].find("4") != string::npos)
    {
        index.push_back(idx);
        label.push_back("4");
    }
  }
 
  infile.close();

//cout << "Done here" << endl;

  unique_ptr<OpenBabel::OBMol> mol(new OpenBabel::OBMol);
  

  OpenBabel::OBConversion conv;
  conv.SetInFormat("xyz");
  conv.SetOutFormat("xyz");
  conv.ReadFile(mol.get(), filename); 
  
  conv.Write(mol.get());

//cout << "stophere  " << index.size() << endl;

  for(int i=0; i<index.size(); i++)
  {
    OpenBabel::OBAtom* atom = mol->GetAtom(index[i]);

    FOR_NBORS_OF_ATOM(nbr, atom)
    {
       handles.push_back(&*nbr);
    }   
  }

  vector<OpenBabel::OBAtom*> tmp_atom;
  
  for(int i=0; i<index.size();i++)
  {
   tmp_atom.push_back(mol->GetAtom(index[i]));
  }

  for(int i=0; i<tmp_atom.size(); i++)
  {
    mol->DeleteAtom(tmp_atom[i]);
  }

//  mol->AddHydrogens();

//  minimize(mol);
/*
  vector<OpenBabel::OBAtom*>::iterator it;

  for (it = handles.begin(); it!=handles.end(); it++)
  {
     FOR_NBORS_OF_ATOM(nbr, *it)
     {
       if(nbr->IsHydrogen())
       {
         mol->DeleteAtom(&*nbr);
         break;
       }
     }   
  } 
*/
//cout << "Done here" << endl;

  return mol;
}

void MOL2XYZ::writeXYZ(unique_ptr<OpenBabel::OBMol> & mol, string dir)
{
   int natoms = mol->NumAtoms();

   vector<string> column(natoms,"X");
   vector<int> xyzw(4,0);
   int idx;
  
   int size = label.size();
// cout << "size of label:" << size << endl;

   for(int i =0; i < size; i++)
   {
 //  cout << i << endl;
     if (label.at(i).find("PLG") != string::npos)   
     {
     idx = handles.at(i)->GetIdx();
//   cout << "new PLG:" << idx << endl;
     column.at(idx-1) = "PLG";
     } 
     else if (label.at(i).find("1") != string::npos)
     {
        idx = handles.at(i)->GetIdx();
        column.at(idx-1) = "1";
        if (xyzw.size() != 1) 
        xyzw[0]=1;
        else xyzw[0]++;        
     }
     else if (label.at(i).find("2") != string::npos)
     {
       idx = handles.at(i)->GetIdx();
       column.at(idx-1) = "2";
       if (xyzw.size() != 2)
       xyzw[1]=1;
       else xyzw[1]++;
     }
     else if (label.at(i).find("3") != string::npos)
     {
        idx = handles.at(i)->GetIdx();
        column.at(idx-1) = "3";
        if (xyzw.size() != 3)
        xyzw[2]=1;
        else xyzw[2]++;
     }   
     else if (label.at(i).find("4") != string::npos)
     {
        idx = handles.at(i)->GetIdx();
        column.at(idx-1) = "4";
        if (xyzw.size() != 4)
        xyzw[3]=1;
        else xyzw[3]++;
     }    
   }

   string filename = dir + "lib/" + title + ".xyz";
// cout << filename << endl;

   ofstream outfile;
   outfile.open(filename.c_str());
   outfile.setf(ios::fixed);
   outfile << setprecision(6);

   natoms = mol->NumAtoms();   
   outfile << natoms << endl;

   vector<int>::iterator it_xyzw;
   for (it_xyzw = xyzw.begin(); it_xyzw!= xyzw.end(); it_xyzw++)
   {
   if(*it_xyzw!=0)
   outfile << *it_xyzw << "   ";
   }  
   outfile << endl;

   int i=0;

   FOR_ATOMS_OF_MOL(atom,*mol)
   {
     outfile << PTable::atom_name(atom->GetAtomicNum()) << "   ";
     outfile << atom->x() << "    ";
     outfile << atom->y() << "    ";
     outfile << atom->z() << "    ";
     outfile << column[i] << endl;
     i++;
   }

  outfile.close();
  return;
}

int MOL2XYZ::minimize(unique_ptr<OpenBabel::OBMol> & mol)
{
//      OpenBabel::OBConversion conv;
//      conv.SetOutFormat("xyz");
//      string ff = "MMFF94s";
        string ff = "UFF";
//      unique_ptr<OpenBabel::OBMol> ptr(new OpenBabel::OBMol());
      OpenBabel::OBForceField* pFF = OpenBabel::OBForceField::FindForceField(ff.c_str());

//Make sure we have a valid pointer
     if (!pFF)
        {
          cout << "cannot find ff: " << ff << endl;
          exit(-1);
        }

        pFF->SetLogFile(&cout);
//        pFF->SetLogLevel(OBFF_LOGLVL_NONE);
        pFF->SetLogLevel(OBFF_LOGLVL_LOW);

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
