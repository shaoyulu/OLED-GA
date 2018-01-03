using namespace std;
#include "OBabel.h"

// Helper function to read molecule from file
unique_ptr<OpenBabel::OBMol> GetMol(string filename)
{
//   // Create the OBMol object.
     unique_ptr<OpenBabel::OBMol> mol(new OpenBabel::OBMol);
//
//       // Create the OBConversion object.
     OpenBabel::OBConversion conv;
     OpenBabel::OBFormat *format = conv.FormatFromExt(filename.c_str());
     if (!format || !conv.SetInFormat(format)) 
     {
     cout << "Could not find input format for file " << filename << endl;
     return mol;
     }
// Open the file.
     ifstream ifs(filename.c_str());
     if (!ifs) 
     {
     cout << "Could not open " << filename << " for reading." << endl;
     return mol;
     }
// Read the molecule.
     if (!conv.Read(mol.get(), &ifs)) {
     cout << "Could not read molecule from file " << filename << endl;
     return mol;
     }

     return mol;
}

//merge two fragment 
void merge_mol(unique_ptr<OpenBabel::OBMol> & mol1, int socket, unique_ptr<OpenBabel::OBMol> & mol2, int plug)
{
   OpenBabel::OBAtom* atom = mol1->GetAtom(socket);

   int i;
   int j;

   FOR_BONDS_OF_ATOM(a, *atom)
   {
     if (socket == a->GetBeginAtomIdx())
       i = a->GetEndAtomIdx();
    else if (socket == a->GetEndAtomIdx())
       i = a->GetBeginAtomIdx();
   }   

  mol1->DeleteAtom(atom);

  int n = mol1->NumAtoms();

  atom = mol2->GetAtom(plug);

  FOR_BONDS_OF_ATOM(b, *atom)
  {
    if (plug == b->GetBeginAtomIdx())
      j = b->GetEndAtomIdx();
    else if (plug == b->GetEndAtomIdx())
      j = b->GetBeginAtomIdx();
  }
  
  mol2->DeleteAtom(atom);

// Add second molecule
   (*mol1.get()) += (*mol2.get());

   j = j + n;    

// Connect the fragments
   OpenBabel::OBBuilder::Connect(*mol1.get(), i, j);

  return;

}


unique_ptr<OpenBabel::OBMol> minimize(unique_ptr<OpenBabel::OBMol> & mol)
{
      unique_ptr<OpenBabel::OBMol> ptr(new OpenBabel::OBMol);
 //     OpenBabel::OBForceField* pFF = OpenBabel::OBForceField::FindForceField("UFF");
        OpenBabel::OBForceField* pFF = OpenBabel::OBForceField::FindForceField("MMFF94");

//Make sure we have a valid pointer
     if (!pFF)    exit(-1);

        pFF->SetLogFile(&cout);
        pFF->SetLogLevel(OBFF_LOGLVL_LOW);
     if (!pFF->Setup(*mol.get())) {
        cerr << "ERROR: could not setup force field." << endl;                                                   }

      pFF->SteepestDescent(1000);
      
      pFF->GetCoordinates(*ptr.get());

      return ptr; 
}
