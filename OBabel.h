#ifndef OB_H
#define OB_H

#include <openbabel/obconversion.h>
#include <openbabel/mol.h>
#include <openbabel/builder.h>
#include <openbabel/obiter.h>
#include <memory> 
#include <iostream>
#include <openbabel/forcefield.h>

namespace OBabel
{
     unique_ptr<OpenBabel::OBMol> GetMol(string filename);

     void merge_mol(unique_ptr<OpenBabel::OBMol> & mol1, int socket, unique_ptr<OpenBabel::OBMol> & mol2, int plug);

     unique_ptr<OpenBabel::OBMol> minimize(unique_ptr<OpenBabel::OBMol> & mol);
};

#endif
