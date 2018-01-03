#include <iostream>
#include <fstream>
#include <stdio.h>
#include <stdlib.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <dirent.h>
using namespace std;

#include "read_inp.h"
#include "mol2xyz.h"

void readinp::read_all()
{
  string ga_inp = "ga.inp";
  
  cout << "************************************" << endl;
  cout << "****    Reading Input Files     ****" << endl;
  cout << "************************************" << endl;


//  string path = "./Scratch/";
//  mkdir(path.c_str(),S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);  
  read_ga(ga_inp);

  string lib_path = "./Scratch/lib/";
  mkdir(lib_path.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
  
  string path = "./Scratch/features/";
  mkdir(path.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
/*hard coded !!!!!!*/
//  int hard_code[] = {0,1,2,3}
//  feature_list(hard_code, hard_code + sizeof(hard_code) / sizeof(int))
/*hard coded !!!!!!*/

  cout << endl << endl;

  string dir = inp_sgen.lib + "/seed";

  inp_sgen.n_seeds = load_lib(0, dir);    

//  read_seed(inp_sgen.seed);
  
  cout << endl << "find " << inp_sgen.n_seeds << " core building blocks in the library" << endl;

  dir = inp_sgen.lib + "/reg";  

  inp_sgen.n_regs = load_lib(inp_sgen.n_seeds, dir);
  
  cout << "find " << inp_sgen.n_regs << " multiple handles building block in the library" << endl;
  
  dir = inp_sgen.lib + "/cap";
  
  inp_sgen.n_caps = load_lib(inp_sgen.n_regs + inp_sgen.n_seeds, dir);
 
  cout << "find " << inp_sgen.n_caps << " single handle building block in the library" << endl;

  cout << endl << endl;
 
  read_ref(ref_f);

//  cout << endl << endl;

//  read_tempxyzf(temp_f);
 
//  cout << endl << endl;

//  get_basisfile(rfolder);

//  cout << "total" << n_rfiles << " R basis " << endl;     
 
  string fname = "dft.inp";
 
  if(inp_ga.use_dft == 1) 
  read_dft(fname);
  
  fname = "mopac.inp";
  
  if(inp_ga.use_mopac == 1)
  read_mopac(fname);

  cout << "********* All inputs are read ********" << endl;
}

/*read GA.inp*/
void readinp::read_ga(string inpfile)
{
  
  cout << "Reading GA.inp: " << inpfile << '\n';
  ifstream infile;
  string filename = "inp/";
  filename += inpfile;
  infile.open(filename.c_str());
  if (!infile)
  {
    cout << "Cannot read inpfile" << '\n';
    exit(-1);
  }

  string line;
  string block;
  string subblock;
  vector<string>::iterator it;       /* iterator for vector<string> */
  vector<string> tok_line; 
  vector<string> tok_line2;
  vector<string> tok_line3;
  int flag=0;
  inp_ga.num_r = 0;          /*number of independent R groups*/
//  bool success=true; 

  while(infile.good() && !infile.eof())
  {
    getline(infile, line);       
//    if (StringTools::iscomment(line)) continue;
//    cout << line << endl;
    line = StringTools::newCleanString(line);   /*discard comments and first a few whitespace*/

   if (line.length() == 0) continue;         /*skip blank lines*/
   else if (line[0] == '$')                  /*start of a new block*/
     { 
      line = StringTools::lowerCase(line);      /*convert to lowerCase*/
      if (line.find("end") == string::npos)           /*not the end of a block*/
        {
          flag++;                                   /*flag of block*/
          if (flag == 1)  
             {
               block = line.substr(1);
               cout << '\n';
               cout << "Start reading $" << block << '\n';
             } 
          else if (flag == 2)  
             {  
                cout << '\n';
                subblock = line.substr(1);
                cout << "Start reading subblock $" << subblock << " in $" << block << '\n';
             } 
        }   /*  end of "if  */
      else  
        { 
          flag--;
//          if (flag == 1)    
//          {
//             cout << "Finish reading subblock $" << subblock << '\n';
//             cout << '\n';
//             subblock.clear();
//          }
   
//          else if (flag == 0 ) 
//          { 
//             cout << "Finish reading $" << block << '\n';
//             cout << "*********************************" << '\n';
//             block.clear();
//          }         
        } /*end of else*/      
    }    /*end of else if*/            
   else 
    {
      tok_line = StringTools::tokenize(line, " ");
      
      line = StringTools::lowerCase(line);

      if (block.find("features") != string::npos)
      {
        if (line.find("polarization") != string::npos)
            {
            if (line.find("on") != string::npos)
            feature_list.push_back("polar");
            }
        if (line.find("homo") != string::npos)
            {
            if (line.find("on") != string::npos)
            feature_list.push_back("HOMO");
            }
        if (line.find("lumo") != string::npos)
            {
            if (line.find("on") != string::npos)
            feature_list.push_back("LUMO");
            }
      }     
      if (block.find("s-generator") != string::npos)
      {
         if (line.find("mode") != string::npos)
         {
            inp_sgen.mode = atoi(tok_line[1].c_str());
            cout << "Building Mode is set to " << inp_sgen.mode << endl;
         }
         else if (line.find("layers") !=string::npos)
         {
            if (tok_line[1].find("-")!= string::npos)
              {
                tok_line2 = StringTools::tokenize(tok_line[1],"-");
                inp_sgen.layer_min = atoi(tok_line2[0].c_str());
                inp_sgen.layer_max = atoi(tok_line2[1].c_str());
                cout << "Will build " << inp_sgen.layer_min << " to " << inp_sgen.layer_max << " layers" << endl; 
              }
            else 
               {
                 inp_sgen.nlayers = atoi(tok_line[1].c_str());
                 cout << "Will  build " << inp_sgen.nlayers << " layers" << endl;
               }
         }
         else if (line.find("symm") != string::npos)
         {
            string onoff = StringTools::lowerCase(tok_line[1]);
            if (onoff.find("on") !=string::npos)
                inp_sgen.symm = true;
            else  inp_sgen.symm = false;
         }
         else if (line.find("seed") != string::npos)
         {
            inp_sgen.seed = tok_line[1];
            cout << "The seed is given by user: " << inp_sgen.seed << endl;            
         }
         else if (line.find("library") != string::npos)
         {
            inp_sgen.lib = tok_line[1];
            cout << "The library is given by user: " << inp_sgen.lib << endl;
         }

        else if (line.find("reference_electrode") != string::npos)
         {
          cout<< "Start reading Reference Electrode..." << endl;
          inp_sgen.HOST_HOMO = atof(tok_line[1].c_str());
          inp_sgen.HOST_LUMO = atof(tok_line[2].c_str());
        
          inp_sgen.TransportType = 1;
           cout << "The ref Anode(eV) = " << inp_sgen.HOST_HOMO << " and Cathode(eV) = " <<
           inp_sgen.HOST_LUMO << endl;
         }

           else if (line.find("jobtype") != string::npos)
         {
            string onoff = StringTools::lowerCase(tok_line[1]);
            if (onoff.find("slurm") !=string::npos)
                inp_sgen.jobtype =1;
            else  inp_sgen.jobtype =2;
         }
         
          else if (line.find("writefragscore")!= string::npos){
            string onoff = StringTools::lowerCase(tok_line[1]);
            if (onoff.find("on") !=string::npos)
                inp_sgen.isWriteFragScore = true;
            else  inp_sgen.isWriteFragScore = false;
          }
         else if (line.find("fragscorethreshold") != string::npos)
         {
          inp_sgen.fscore = atof(tok_line[1].c_str());
         }


         else if (line.find("best_restart")!= string::npos){
            string onoff = StringTools::lowerCase(tok_line[1]);
            if (onoff.find("on") !=string::npos)
                inp_sgen.bestRestart = true;
            else  inp_sgen.bestRestart = false; 
          } 

      }
#if 0     
      if (block.find("templete_cats")!=string::npos)
      {
        if (line.find("ox_state")!=string::npos) 
          {
             inp_ga.ox_state = atoi(tok_line[1].c_str());
             cout << "Oxidative State of TM is " << inp_ga.ox_state << endl;  
          }
//        else if (line.find("charge")!=string::npos)
//          {           
//             charge = atoi(tok_line[1].c_str());
//             cout << "Charge of reference cat is " << charge << endl;
//          }
        else if (line.find("reference")!=string::npos)
          {
             ref_f = tok_line[1];
             cout << "Filename of reference cat is " << ref_f << endl;
          } 
        else if (line.find("template")!=string::npos)
          {
          tok_line2 = StringTools::tokenize(tok_line[1],",");
          for (size_t i=0; i < tok_line2.size(); i++)
            {
            if (tok_line2[i].find("-")!=string::npos)
              {
               tok_line3 = StringTools::tokenize(tok_line2[i],"-");
               for (int j= atoi(tok_line3[0].c_str()); j < atoi(tok_line3[1].c_str())+1; j++)
                {
                temp_f.push_back(StringTools::genfilename("Template_cats_","xyz",4,j));
                cout << "Filename of template cat is " << temp_f.back() << '\n'; 
                }
//                tok_line3.clear();
//    line.clear();
               }
            else 
              {
              temp_f.push_back(StringTools::genfilename("Template_cats_","xyz",4,atoi(tok_line2[i].c_str())));
              cout << "Filename of template cat is " << temp_f.back() << '\n';
              }          
            }   /*end of for loop*/
/*write a tools for converting string of range to int*/
//          for (int i) 
//          StringTools:genfilename("Templete_cats_","xyz",4,nxyz)
           }  /*end of else if*/ 
//catalyst named as Templete_cats_*.xyz
        }  /*end of if*/

     else if (block.find("initiation")!=string::npos)
     { 
       if (subblock.find("r_library")!=string::npos) 
       {
        rfolder = tok_line[0];
//        for (it = tok_line.begin(); it!= tok_line.end(); it++)
//        {
//        tmp.push_back(*it);
//        }   
//        rfiles.push_back(tmp);
       cout << "Library of R basis is " << rfolder << '\n'; 
       }
       else if (subblock.find("r_symmetry")!=string::npos)
       {
        vector<string> tmp;
        for (it = tok_line.begin(); it != tok_line.end(); it++)
         {
           tmp.push_back(*it);
           if(it != tok_line.begin()) 
           cout << *it << " ";
         }
         cout << "are identical." << '\n';
         r_sym.push_back(tmp);
       
       inp_ga.num_r++;
       }
     }     /*end of else if*/
#endif
   
     else if (block.find("fitness")!=string::npos)
     {
       if (subblock.find("qc_method")!=string::npos)
       {
         if (line.find("dft")!=string::npos) 
         {
           /*enable dft calculation & read DFT.inp*/
           if (line.find("on")!=string::npos)  
              {
               int status;
               string path = "./Scratch/dft";
               status = mkdir(path.c_str(),S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
               inp_ga.use_dft = 1;
               if(status)
               cout << "DFT calculation is ON." << '\n';
              }   
           else if (line.find("off")!=string::npos) inp_ga.use_dft = 0;
           else cout << "unrecognized setting for DFT calculation" << endl; 
         }  
         else if (line.find("semi_empirical")!=string::npos)
         {
           /*enable semi-empirical calculation & read MOPAC.inp*/
           if (line.find("on")!=string::npos) 
              {
                int status;
                string path = "./Scratch/mopac";
                status = mkdir(path.c_str(),S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
                inp_ga.use_mopac =1;
                if(status)
                cout << "Semi-empirical calculation is ON." << '\n';
              }
           else if (line.find("off")!=string::npos) inp_ga.use_mopac = 0;
           else cout << "unrecognized setting for semi-empirical calculation" << endl;
         }
        } 
       else if  (line.find("def_file")!=string::npos)
          {
          fitness_f = tok_line[1]; 
          cout << "Using file " << fitness_f << " as fitness function" << '\n'; 
          }
       else if (line.find("reference") != string::npos)
         {
         ref_f = tok_line[1];
         cout << "Reference molecule is given by user:" << ref_f << endl;
         }
     }
   
     else if (block.find("ga")!=string::npos)
     {
       if (subblock.find("selection")!=string::npos) 
       {
            if (line.find("method") != string::npos)
            {
                if (tok_line[1].find("RWS") != string::npos)
                    inp_ga.RWS = true;
                if (tok_line[1].find("SUS") != string::npos)
                    inp_ga.SUS = true;           
                if (tok_line[1].find("B_Tournament")!= string::npos)
                    inp_ga.B_Tournament = true;
                if (tok_line[1].find("L_Rank") != string::npos)
                    inp_ga.L_Rank = true;
                if (tok_line[1].find("E_Rank") != string::npos)
                    inp_ga.E_Rank = true;
            }
            else if (line.find("parameter") != string::npos)
            {
                if (inp_ga.L_Rank == true) 
                inp_ga.L_Rank_rate = atof(tok_line[1].c_str());
                else if (inp_ga.E_Rank == true)
                    inp_ga.E_Rank_base = atof(tok_line[1].c_str());
            }
       }
       else if (line.find("force_mutation") != string::npos)
         {
       if (tok_line[1].find("true") != string::npos){
          inp_ga.force_mutation = true;
          cout << "Force mutation will be applied" << endl;
        }else
          inp_ga.force_mutation = false;    
      
         }     
       else if (line.find("population_size")!=string::npos)
          inp_ga.pop_size = atoi(tok_line[1].c_str());
       else if (line.find("termination")!=string::npos)
          inp_ga.ga_stop = atoi(tok_line[1].c_str());
       else if (line.find("crossoverrate") != string::npos)
          inp_ga.cross_r = atof(tok_line[1].c_str());
       else if (line.find("mutationrate")!= string::npos)
          inp_ga.muta_r = atof(tok_line[1].c_str());
     } 

     else if (block.find("output")!=string::npos)
     {
       if (line.find("molden_format")!=string::npos)
       {
         if (line.find("on")!=string::npos)  
            {
              inp_ga.out_molden = 1;
              cout << "Output as Molden Format" << '\n';
            }
                
         else if (line.find("off")!=string::npos) inp_ga.out_molden = 0;
         else cout << "unrecognized setting for output format" << endl;
       }      
       else if (line.find("excel_format")!=string::npos)
       {
         if (line.find("on")!=string::npos) 
            {
               inp_ga.out_excel = 1;
               cout << "Output as Excel Format" << '\n';
            }
         else if (line.find("off")!=string::npos) inp_ga.out_excel = 0;
         else cout << "unrecognized setting for output format" << endl;
       }    
     }    

    }     /*end of else*/
   }                /*end of while loop*/

   infile.close();

   cout << "******* Finish Reading GA.inp File *******" << endl;

  return;
}

void readinp::read_ref(string filename)
{
  cout << "Reading reference provided by user: " << endl;
  ifstream infile;
  string fname = "./inp/";
  fname += filename;
  infile.open(fname.c_str());

  if(!infile)
  {
    cout << "Cannot read " << fname << endl;
    exit(-1); 
  }

  string line;
  vector<string> tok_line;

  getline(infile,line);
  line = StringTools::newCleanString(line);
  ref_xyz.natoms = atoi(line.c_str());
  getline(infile,line);

  ref_xyz.anames.resize(ref_xyz.natoms);
  ref_xyz.coords.resize(3*ref_xyz.natoms);

  for(int i =0; i< ref_xyz.natoms; i++)
  {
     getline(infile, line);
   tok_line = StringTools::tokenize(line, " ");
   ref_xyz.anames[i] = tok_line[0];
   ref_xyz.coords[3*i+0] = atof(tok_line[1].c_str());
   ref_xyz.coords[3*i+1] = atof(tok_line[2].c_str());
   ref_xyz.coords[3*i+2] = atof(tok_line[3].c_str());
  }     
  
  return;
}
void readinp::read_seed(string filename)
{
    string lib_path = "./Scratch/lib/";
    mkdir(lib_path.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);

  cout << "Reading seed block" << endl;
//  ifstream infile;
  string fname = "./inp/";
  fname += filename;
//  infile.open(fname.c_str());
/* 
  MOL2XYZ mol2xyz;
  unique_ptr<OpenBabel::OBMol> mol(mol2xyz.readMOL(fname));
  mol2xyz.title = "0";   //or = "seed" depend on gcode class

  mol2xyz.writeXYZ(mol,lib_path);
*/

  string cmd = "cp " + fname + " " + lib_path + "0.xyz";
  system(cmd.c_str());
   
  ifstream seed_xyz;
//  fname = lib_path + mol2xyz.title + ".xyz";  
  seed_xyz.open(fname.c_str());     
      
  seed.nhandles = 0;

  if(!seed_xyz)
  {
  cout << "Cannot read " << fname << endl;
  exit(-1);
  }
   
  string line;
  vector<string> tok_line;

  getline(seed_xyz,line);
  line = StringTools::newCleanString(line);
  seed.natoms = atoi(line.c_str());
  getline(seed_xyz,line);
  line = StringTools::newCleanString(line);
  tok_line = StringTools::tokenize(line, " ");
  vector<string>::iterator it;
  for (it = tok_line.begin(); it != tok_line.end(); ++it)
  {
  int tmp = atoi((*it).c_str());
//  seed.nhandles += tmp;
  seed.xyzw.push_back(tmp);
  }  

  seed.name.resize(seed.natoms);
  seed.label.resize(seed.natoms);
  seed.coords.resize(3*seed.natoms);

  for(int i=0; i<seed.natoms; i++)
  {
    getline(seed_xyz, line);
    tok_line = StringTools::tokenize(line, " ");
    seed.name[i] = tok_line[0];
    seed.label[i] = tok_line[4];
    seed.coords[3*i+0] = atof(tok_line[1].c_str());
    seed.coords[3*i+1] = atof(tok_line[2].c_str());
    seed.coords[3*i+2] = atof(tok_line[3].c_str());
  }

  seed.nbranches = int(seed.xyzw.size());

  cout << "seed has " << seed.nhandles << " handles" << endl;

  for (int i=0; i<seed.natoms; i++)
  {
  cout << seed.name[i] << "\t" << seed.coords[3*i+0] << "\t" << seed.coords[3*i+1] << "\t" << seed.coords[3*i+2] << endl;
  }
/*
  if (inp_sgen.mode == 1)
  seed.nhandles = 1; 
*/           // for high symmetric building

  string path = "./Scratch/features/";
  mkdir(path.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
}

int readinp::load_lib(int n, string libname)
{

  struct dirent *pDirent;
  DIR *pDir;
  int count = 0;
  string str;
  vector<string> basislist;
  string filename;
  string cmd;
  string newfilename;

  pDir = opendir(libname.c_str());

  if (pDir != NULL) 
  {
    while ((pDirent = readdir(pDir)) != NULL) 
      {
        str = string(pDirent->d_name);

        if (str.find(".xyz") != string::npos) 
          {
               basislist.push_back(str);
               count++;
               cout << "new label for " <<  str << " in library would be " << n+count << endl; 
          }
        }
    closedir (pDir);
   }
  
   string path = "./Scratch/";
//   mkdir(path.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);

   for (int i=0;  i < count; i++)
   {
      MOL2XYZ mol2xyz;
      filename = libname + "/" + basislist[i];

      unique_ptr<OpenBabel::OBMol> mol(mol2xyz.readXYZ(filename));

      newfilename = StringTools::int2str(i+n+1, 1, "0");

      mol2xyz.title = newfilename;

      newfilename = path + "features/" + newfilename + ".xyz";
      cout << newfilename << endl;

      mol2xyz.writeXYZ(mol,path);
      cmd = "cp " + filename + " " + newfilename;
      system(cmd.c_str());
   }

   return count;
}

void readinp::read_dft(string inpfile)
{
  cout << "Reading QChem-DFT Inpfile" << endl;
  ifstream infile;
  string filename = "./inp/";
  filename += inpfile;
  infile.open(filename.c_str());

  if(!infile)
  {
  cout << "Cannot read " << inpfile << endl;
  exit(1);
  }

  sp.initial();
  opt.initial();
//  ts.initial();

  string block;
  string line;
  vector<string> tok_line;

  while(!infile.eof())
  {
    getline(infile,line);
    line = StringTools::newCleanString(line);

    if (line.length()==0) continue;
    else if (line[0] == '$')
    {
      line = StringTools::lowerCase(line);
      if (line.find("end") == string::npos)
      {
       block = line.substr(1);

//      cout << "find block:" << block << endl;
      }
    }

    else
    {
      line = StringTools::lowerCase(line);
      tok_line = StringTools::tokenize(line," ");

      if (block.find("sp") != string::npos)
      {

       if(line.find("functional") != string::npos)
          sp.functional = tok_line[1];

       else if (line.find("unrestricted") != string::npos)
          sp.unres = tok_line[1];

       else if (line.find("scf_algorithm") != string::npos)
          sp.scf_algo = tok_line[1];

       else if (line.find("scf_max_cycles") != string::npos)
         sp.scf_max_cycles = atoi(tok_line[1].c_str());

       else if (line.find("basis") != string::npos)
         sp.basis = tok_line[1];

       else if (line.find("wavefunction_analysis") != string::npos)
         sp.wavefunction_analysis = tok_line[1];

       else if (line.find("scf_convergence") != string::npos)
         sp.scf_conv = atoi(tok_line[1].c_str());

       else if (line.find("sym_ignore") != string::npos)
         sp.sym_ignore = tok_line[1];

       else if (line.find("symmetry") != string::npos)
         sp.sym = tok_line[1];

      }    /*end of $SP block*/

      else if (block.find("opt") != string::npos)
      {
       if(line.find("functional") != string::npos)
          opt.functional = tok_line[1];

       else if (line.find("unrestricted") != string::npos)
          opt.unres = tok_line[1];

       else if (line.find("scf_algorithm") != string::npos)
          opt.scf_algo = tok_line[1];

       else if (line.find("scf_max_cycles") != string::npos)
         opt.scf_max_cycles = atoi(tok_line[1].c_str());

       else if (line.find("basis") != string::npos)
         opt.basis = tok_line[1];

       else if (line.find("wavefunction_analysis") != string::npos)
         opt.wavefunction_analysis = tok_line[1];

       else if (line.find("scf_convergence") != string::npos)
         opt.scf_conv = atoi(tok_line[1].c_str());

       else if (line.find("sym_ignore") != string::npos)
         opt.sym_ignore = tok_line[1];

       else if (line.find("symmetry") != string::npos)
         opt.sym = tok_line[1];

       else if (line.find("geom_opt_max_cycles")!= string::npos)
         opt.opt_max_cycles = atoi(tok_line[1].c_str());

       else if (line.find("geom_opt_tol_displacement")!= string::npos)
         opt.opt_tol_displacement = atoi(tok_line[1].c_str());

       else if (line.find("geom_opt_tol_gradient")!= string::npos)
         opt.opt_tol_gradient = atoi(tok_line[1].c_str());

       else if (line.find("geom_opt_tol_energy") != string::npos)
         opt.opt_tol_energy = atoi(tok_line[1].c_str());
      }  /*end of $OPT block*/
    }

  }    /* end of while*/
  infile.close();

  return;
}

#if 0
void readinp::read_refxyzf(string inpfile) 
{
  
   cout << "Reading inpfile for catalyst:\n " << inpfile << endl;
   string filename = "scratch/";
   filename += inpfile;

//   cout << filename << endl;

   ifstream infile;
   infile.open(filename.c_str());

   if (!infile)
   {
     cout << "cannot read inpfile" << endl;
     exit(-1);
   }

  xyzf_ref.natoms=0;
  xyzf_ref.nlabels=0;
  string line;
  vector<string> tok_line;
  vector<string> tok_line2; 
  vector<string> tok_line3; 

  getline(infile, line);
  xyzf_ref.natoms = atoi(line.c_str());
  
  getline(infile,line);
  line = StringTools::newCleanString(line);
  tok_line = StringTools::tokenize(line,",");
  xyzf_ref.charge = atoi(tok_line[0].c_str());
  xyzf_ref.multi = atoi(tok_line[1].c_str()); 

  int i=0;
  size_t natoms = unsigned(xyzf_ref.natoms);
  xyzf_ref.element.reserve(natoms);
  xyzf_ref.coords.reserve(3*natoms);

  for(;;)
  { 
//    line.clear();
    getline(infile, line);
    line = StringTools::newCleanString(line);

    if (line.length()!=0 && line.find("$") == string::npos )
    {  
       tok_line = StringTools::tokenize(line, " ");
       i++;  
       xyzf_ref.element.push_back(tok_line[0]);
       for (int j=1; j < 4; j++)
       {
        xyzf_ref.coords.push_back(atof(tok_line[j].c_str()));       
       }
//       tmp.clear();
//       tok_line.clear();
     }
     else break;
  } 
                               /*store .xyz of reference in string format*/
  
  if (i != xyzf_ref.natoms)  
  {
    cout << "The number of atoms in ref.xyz does't match" << endl;
    exit(-1);
  }                       /*check the number of atoms in ref_cat*/

  string block;
  vector<string> tmp2;
  while (!infile.eof())
  {
//    line.clear();

    getline(infile,line);

//    tok_line.clear();
//    tok_line2.clear(); 
    tmp2.clear();

 //    if (StringTools::iscomment(line)) continue;
    line = StringTools::newCleanString(line);
    if (line.empty()) continue;               /*skip blank lines*/
    else if (line[0] == '$')                  /*start of a new block*/
    {
      line = StringTools::lowerCase(line);      /*convert to lowerCase*/
      if (line.find("end") == string::npos)           /*not the end of a block*/
         block = line.substr(1);    
//      else  block.clear();
     }

    else 
    {
      tok_line = StringTools::tokenize(line," ");
      line = StringTools::lowerCase(line);

      if (block.find("label")!=string::npos)
      {  
	 xyzf_ref.nlabels++;        /*need this info?*/
         tmp2.push_back(tok_line[0]);
         tok_line2 = StringTools::tokenize(tok_line[1],",");
         for (size_t i=0; i < tok_line2.size(); i++)
           {
           if (tok_line2[i].find("-")!=string::npos)
             {
             tok_line3 = StringTools::tokenize(tok_line2[i],"-");
             for (int j = atoi(tok_line3[0].c_str());j< atoi(tok_line3[1].c_str())+1; j++)
                tmp2.push_back(StringTools::int2str(j,4," ")); 
//             tok_line3.clear();	 
             }    
           else
              tmp2.push_back(tok_line2[i]);
           }
	 xyzf_ref.label.push_back(tmp2);      /*size == nlabels*/
        } 
      #if 0       
      else if (block.find("connect")!=string::npos) 
      {	
	 tmp2.push_back(tok_line[0]);
	 tok_line2 = StringTools::tokenize(tok_line[1],":");
         tmp2.push_back(tok_line2[0]);     /*order may vary, scaff:R or R:scaff*/
	 tmp2.push_back(tok_line2[1]);
	 xyzf_ref.nodes.push_back(tmp2);
      }
      #endif
    }
  }    /*end of while*/
  infile.close();

  return;
}  /*end of function*/

#endif

void readinp::read_tempxyzf(vector<string> & inpfiles)
{

 using namespace std;

 xyzfile tmp3;

 n_temp = 0;

 for(size_t i=0; i < inpfiles.size() ; i++)
 { 

  cout << "Reading inpfile for template:\n " << inpfiles[i] << endl; 

  string filename;
  filename = "scratch/";
  filename += inpfiles[i];  
  ifstream infile;
  infile.open(filename.c_str());

  if (!infile)
  {
    cout << "cannot read inpfile" << endl;
    exit(-1);
  } 

  int natoms =0;
  int nlabels=0;
  string line;
  vector<string> tok_line;
  vector<string> tok_line2;
  vector<string> tok_line3;

  getline(infile, line);
  natoms = atoi(line.c_str());
  getline(infile,line);   /*read the comment line of xyz file*/
  line = StringTools::newCleanString(line);
  tok_line = StringTools::tokenize(line, ",");
  int charge = atoi(tok_line[0].c_str());
  int multi = atoi(tok_line[1].c_str());

  int k=0;
  vector<string> tmp2;

  tmp3.natoms = natoms;
  tmp3.charge = charge;
  tmp3.multi = multi;
 
  size_t natom = unsigned(natoms);
  tmp3.coords.reserve(3*natom);
  tmp3.element.reserve(natom);

  for(;;)
  {
 //   line.clear();
    getline(infile,line);
    line = StringTools::newCleanString(line);
    if (line.length()!=0 && line.find("$") == string::npos)
    {
      tok_line = StringTools::tokenize(line," ");
      k++;
      tmp3.element.push_back(tok_line[0]);

      for (int j=1; j< 4; j++)
      {
        tmp3.coords.push_back(atof(tok_line[j].c_str()));
      }

//      tmp.clear();
//      tok_line.clear();
    }
    else break;    
 } 
// while(line.find("XYZ") == string::npos || line.find ("xyz") == string :: npos);
 /*end of for loop*/

  if (k != natoms)
  {
    cout << "The number of atoms in " << filename << " doesn't match" << endl;
    exit(-1);
  }

  string block;

  while (!infile.eof())
  {
//    line.clear();
    getline(infile, line);

//    tok_line.clear();
//    tok_line2.clear();
      tmp2.clear();

   line = StringTools::newCleanString(line);
 
   if (line.empty()) continue;
   else if (line[0] == '$')
   {
     line = StringTools::lowerCase(line);
     if (line.find("end") == string::npos)
        block = line.substr(1);
//     else block.clear();
   }

   else
   {
     tok_line = StringTools::tokenize(line, " ");
     line = StringTools::lowerCase(line);
   
     if (block.find("label")!=string::npos)
     {
       nlabels++;
       tmp2.push_back(tok_line[0]);
       tok_line2 = StringTools::tokenize(tok_line[1],",");
       for (size_t i=0; i < tok_line2.size(); i++)
       {
       if (tok_line2[i].find("-")!=string::npos)
       {
           tok_line3 = StringTools::tokenize(tok_line2[i],"-");
          for (int j = atoi(tok_line3[0].c_str()); j < atoi(tok_line3[1].c_str())+1; j++)
          {
             tmp2.push_back(StringTools::int2str(j,4," "));     /*for c++11*/
          }
//          tok_line3.clear();
       }  
       else  
         tmp2.push_back(tok_line2[i]);
       }      /*end of for loop*/

       tmp3.label.push_back(tmp2);
        
   
       tmp3.nlabels = nlabels;
     }
     #if 0
     else if (block.find("connect")!=string::npos)
     {
       tmp2.push_back(tok_line[0]);
       tok_line2 = StringTools::tokenize(tok_line[1],":");
       tmp2.push_back(tok_line2[0]);   /*order may vary, scaff:R or R:scaff*/
       tmp2.push_back(tok_line2[1]);
       tmp3.nodes.push_back(tmp2);
     }
     #endif
   
//      xyzf_temp.push_back(tmp3);
//      tmp3.reset();   
     }
   }  /*end of while loop*/

/**************************************************/
//gen joint list: |V1|V2|=|3 4 5 6|1 2 7 8|  pos = label value-1, pos in element array
  vector<string>::iterator it;
  int pos;
  vector<int> Vn;
//  string var;
// 
for  (size_t j =0; j < r_sym.size(); j++)
 {  
  for (it = r_sym[j].begin()+1; it != r_sym[j].end(); it++)
  {
    for (int i = 0; i< tmp3.nlabels; i++)
    {
      if (tmp3.label[i][0] == *it)
      {
        pos = atoi(tmp3.label[i][1].c_str())-1;
        cout << pos ;
        Vn.push_back(pos); 
        break;
      }
    }
   }
  tmp3.joint.push_back(Vn);
  vector<int>().swap(Vn);
  cout << endl;
  }     
/**************************************************/
   xyzf_temp.push_back(tmp3);
   tmp3.reset();
   n_temp++;  //number of templates

   cout << "n_temp= " << n_temp << endl;   
   infile.close();

 } /*end of for loop*/
   
   return;

}  /*end of function*/  

void readinp::get_basisfile(string filedir)
{
  DIR *dpdf;
  struct dirent *epdf;
  vector<string> rlist; 
  string str;
  n_rfiles = 0;
  string filename;
  string newfilename;
  string cmd;

  dpdf = opendir(filedir.c_str());
  if (dpdf != NULL)
  {
   while ((epdf = readdir(dpdf)) != NULL)
    {
      str = string(epdf->d_name);
      if ((str != ".") && (str != ".."))
      {
      cout << "R filename:" << str << endl;
      rlist.push_back(str);
      n_rfiles++;
      }
    }
   closedir(dpdf);
   }
  
  else 
 {
  cout << "cannot open R folder" << endl;
  exit(-1);
 }


  int status;
  string path = filedir + "data/";
  status = mkdir(path.c_str(),S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);

  for (int i=0; i < n_rfiles; i++)
  {
     filename = filedir + rlist[i];
     newfilename = StringTools::int2str(i+1, 1, "0");
     newfilename = StringTools::newCleanString(newfilename);
     newfilename = newfilename + ".basis";
     newfilename = path  +  newfilename;
     cmd = "cp " + filename + " " + newfilename;
     status = system(cmd.c_str());     
  } 

  cout << "Finish generating a list of rbasis file" << endl;

  return;  
}

void readinp::read_mopac(string inpfile)
{
  cout << "Reading MOPAC settings" << endl;
  ifstream infile;
  string filename = "./inp/";
  filename = filename + inpfile;
  infile.open(filename.c_str());

  if(!infile)
  {
  cout << "Cannot read " << inpfile << endl;
  exit(1);
  }
  
  string block;
  string line;
  vector<string> tok_line;

  while(getline(infile,line))
  {    
//    getline(infile,line);
    line = StringTools::newCleanString(line);

    cout << line << endl;
  
    if (line.length() == 0) continue;
    else
    {
    tok_line = StringTools::tokenize(line," ");
    line = StringTools::lowerCase(line);

    if (line.find("method") != string::npos)
     {
       mopac.method = tok_line[1];
       cout << "method for mopac is" << mopac.method << endl;
     }
    else if (line.find("walltime") != string::npos)
     { 
      mopac.walltime = tok_line[1];
      cout << "walltime for mopac is " << mopac.walltime << endl;
     }
//    else if (line.find("use_aux") != string::npos)
//      {
//      string tmp = tok_line[1];
//      if (tmp.find("on") != string::npos)
//         mopac.use_aux = true;
//      else  mopac.use_aux = false;
//      }   
    else if (line.find("gnorm") != string::npos)
      {
        mopac.gnorm = atof(tok_line[1].c_str());
        cout << "gnorm for mopac is " << mopac.gnorm << endl;
      }
    else if (line.find("unrestricted") != string::npos)
      {
      string tmp = tok_line[1];
      if(tmp.find("on") != string::npos)  
      mopac.uhf = true;  
      else mopac.uhf = false;
      }
//more keywords could be added later!
    }

}

    infile.close();

}

