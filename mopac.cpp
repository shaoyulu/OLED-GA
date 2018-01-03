#include "mopac.h"
#include <stdio.h>
#include <sys/stat.h>

using namespace std;

int Mopac::opt_fragment(string inpfile, readinp* pInp,
                vector<string> & anames, vector<double> & coords)
 {
     ofstream outfile;
     outfile.open(inpfile.c_str());
     outfile.setf(ios::fixed);
     outfile << setprecision(6);

///Geometric Optimization Parameters//
//       outfile << " " << pInp->mopac.method << "XYZ UHF T=" << pInp->mopac.walltime << " SINGLET BONDS AUX ";

//// without GO/////////
   outfile << " " << pInp->mopac.method << " NOSYM T=" << pInp->mopac.walltime << " GNORM=    " << pInp->mopac.gnorm << " LOCALIZE ";




     outfile << endl;
     outfile << "    MOPAC run   " << endl;
     outfile << "    ready to go?    " << endl;

     for (int i=0; i < anames.size(); i++)
     outfile << " " << anames[i] << "  " << coords[3*i+0] << " 1 " << coords[3*i+1] << " 1 "    << coords[3*i+2] << " 1 " << endl;

     outfile.close();

     string cmd = "~/MOPAC2016.exe " + inpfile;
     system(cmd.c_str());

     if(check_success(inpfile))
     {
     xyz_read(inpfile, anames, coords);
     return 1;
     }
     else
     return 0;
 }

void Mopac::opt_header(ofstream & inpfile, readinp* pInp, Node<GANode>* pNode)
{

/////GOPT///////
//  inpfile << " " << pInp->mopac.method << " XYZ T=" << pInp->mopac.walltime << " GNORM=    " << pInp->mopac.gnorm << " SINGLET BONDS AUX  CHARGE=" << pNode->data.charge;

// inpfile << " " << pInp->mopac.method << " GEO-OK T=" << pInp->mopac.walltime << " GNORM=    " << pInp->mopac.gnorm << " SINGLET BONDS AUX  CHARGE=" << pNode->data.charge;

////Original//
inpfile << " " << pInp->mopac.method << " NOSYM T=" << pInp->mopac.walltime << " GNORM="  << pInp->mopac.gnorm << " CHARGE=" << pNode->data.charge;



  if(pInp->mopac.uhf)
  inpfile << " UHF";

//  if(pInp->mopac.aux)
//  inpfile << " AUX";

  if(pNode->data.multi == 3)
  inpfile << " TRIPLET";

  inpfile << endl;
  inpfile << "    MOPAC run  " << endl;
  inpfile << "    ready to go?  " << endl;

  return;
}

int Mopac::gen_inp(string mopac_folder, readinp* pInp, Node<GANode>* pNode)
{
//cout << "pre-opt new cat with mopac (core freeze)" << endl;
  string filename = "m"+ pNode->data.comment;

  cout << "generating inpfile: " << filename << endl;
//cout << "folder: " << mopac_folder << endl;

  ofstream inpfile;
  string inpfile_string = mopac_folder + filename;
  inpfile.open(inpfile_string.c_str());
  inpfile.setf(ios::fixed);
  inpfile.setf(ios::left);
  inpfile << setprecision(6);
  Mopac::opt_header(inpfile, pInp, pNode);

//cout << "pNode natoms = " << pNode->data.natoms << endl;

  for (size_t i=0; i< pNode->data.natoms; i++)
  inpfile << " " << pNode->data.anames[i] << " " << pNode->data.coords[3*i+0] << " 1 " << pNode->data.coords[3*i+1] << " 1 " << pNode->data.coords[3*i+2] << " 1 " << endl;

  inpfile.close();
 
  cout << "successfully generate a mopac input" << endl;

   return 1;
  
}
#if 0
int Mopac::opt(string mopac_folder, readinp* pInp, xyznode* pNode, vector<int> frezlist)
{
  cout << "pre-opt new cat with mopac (core freeze)" << endl;
//  string fname = StringTools::int2str(count, 4,"0");
//  string filedir = "./Work/mopac/"; 
  string filename = "m"+ pNode->comment;

  cout << "generating inpfile: " << filename << endl;

  ofstream inpfile;
  string inpfile_string = mopac_folder + filename;
//  ofstream xyzfile;
//  string xyzfile_string = "testmopac.xyz";

  inpfile.open(inpfile_string.c_str());
  inpfile.setf(ios::fixed);
  inpfile.setf(ios::left);
  inpfile << setprecision(6);

//  xyzfile.open(xyzfile_string.c_str());
//  xyzfile.setf(ios::fixed);   
//  xyzfile.setf(ios::left);
//  xyzfile << setprecision(6);
//  xyzfile << " " << natoms >> endl << endl;

  Mopac::opt_header(inpfile, pInp, pNode);

/*
  int ntemp = pNode->ntemp-1;
  int temp_natoms = pInp->xyzf_temp[ntemp].natoms;
  vector<vector<int> > joint = pInp->xyzf_temp[ntemp].joint;
  
  vector<int> frezlist(pNode->natoms);

  for (size_t i=0; i < pNode->natoms; i++) frezlist[i]=0;
  for (int i=0; i < temp_natoms; i++)  frezlist[i] = 1;
  for (size_t i =0; i< joint.size(); i++)
  {
     for(size_t j=0; j< joint[i].size(); j++)
        frezlist[joint[i][j]] = 0;
  }
*/
    
  for (size_t i=0; i< pNode->natoms; i++)
  if (frezlist[i])
  inpfile << " " << pNode->element[i] << " " << pNode->coords[3*i+0] << " 0 " << pNode->coords[3*i+1] << " 0 " << pNode->coords[3*i+2] << " 0 " << endl;
  else
  inpfile << " " << pNode->element[i] << " " << pNode->coords[3*i+0] << " 1 " << pNode->coords[3*i+1] << " 1 " << pNode->coords[3*i+2] << " 1 " << endl;

  inpfile.close();
  
  string cmd = "/tmp/MOPAC2012.exe "+ inpfile_string;
  system(cmd.c_str());  
  
//  int success = 1;
  int success =  Mopac::check_success(inpfile_string); 
 
  if(success)  
  {
  cout << "mopac is finished successfully!" << endl;
  Mopac::xyz_read(mopac_folder,pNode);
  return 1;
  }
  else
  return 0;
}
#endif

#if 1
int Mopac::save_all_into_list(string filedir, LinkedList<Node<GANode> >* pList)
{
    Node<GANode>* pNode;
    pNode = pList->gethead();
    int count =0;
  
    while (count < pList->listlength)
    {
       cout << endl;
       cout << "*******The " << count+1 << "-th individual in this generation*******" << endl;
       if(pNode->data.flag == 0)
       {
//       cout << "start reading mopac output..." << endl;
         Mopac::readout(filedir, pNode);
       }

       pNode = pNode->NextNode();
       count++;
    }
 
     return 1;
}
#endif
#if 0
int Mopac::hard_code_save_alltolist(string filedir, LinkedList<Node<GANode> >* pList)
{
    Node<GANode>* pNode;
    pNode = pList->gethead();
    int count = 0;
    vector<int> boundry_list;
    int off_set = 1;
    int old_label;

    boundry_list.push_back(9999);
    
    for (int m = 12; m>=0; m--)
    {
        for (int n = 12; n >=0; n--)
        {
            boundry_list.push_back(4+n*13+m*13*13);
        }
    }
    while (count < pList->listlength)
    {
        cout << endl;
        cout << "*******The " << count+1 << "-th individual in this generation*******" << endl;
        if(pNode->data.flag == 0)
        {
            
          if (count == boundry_list.back()) 
          {
            off_set++;
            old_label = count + off_set;
            boundry_list.pop_back();
          }
          else if (count < boundry_list.back())
          {
             old_label = count + off_set;
          }    
            
          cout << "offset: " << off_set << endl; 
          cout << " go to " << old_label << ".out :" << pNode->data.gene << endl;
          Mopac::readout(filedir, pNode, old_label);
        }

        pNode = pNode->NextNode();
        count++;
    }
}
#endif
void Mopac::readout(string mopac_folder,Node<GANode>* pNode)
{
//  string filedir = "Work/mopac/";

//  string filelabel = StringTools::int2str(old_label, 4, "0");  
//hard_coded one:
//  string oname = mopac_folder + "m" + filelabel + ".out";
//cout << "mopac output : "<<  mopac_folder << endl;
  string oname = mopac_folder + "m" + pNode->data.comment + ".out";
  ifstream output(oname.c_str(),ios::in);

//cout<< oname << endl;

  if(output.fail())
  {
    cout << "Cannot read inpfile" << endl;
    pNode->data.HOMO = 0.00000;
    pNode->data.LUMO = 0.00000; 
    return;
//    exit(-1);
  }

  string line;
  vector<string> tok_line;
  int count = 0;
  int i = 0;
  while(getline(output,line))
  {

  if (line.find("HOMO LUMO ENERGIES (EV)") != string::npos)
  {
    cout << pNode->data.gene  << endl;

  tok_line = StringTools::tokenize(line, " ");
//  cout << tok_line[5] << endl;
//  cout << tok_line[6] << endl;

  pNode->data.HOMO = atof(tok_line[5].c_str());
  cout << "HOMO(ev) : "<<pNode->data.HOMO << endl;
  pNode->data.LUMO = atof(tok_line[6].c_str());
  cout << "LUMO(eV) : "<<pNode->data.LUMO << endl;
  count = 1;
  }
 
  if (line.find("CARTESIAN COORDINATES")!= string::npos && count == 1)
   {
//     cout << "find optimized coordinates" << endl;
     getline(output, line); 
//   cout <<"natoms = " << pNode->data.natoms << endl;
     for(int i=0; i < pNode->data.natoms; i++)
     {
     getline(output, line);
//     cout << line << endl;
     tok_line = StringTools::tokenize(line, " ");
     pNode->data.coords[3*i+0] = atof(tok_line[2].c_str());
     pNode->data.coords[3*i+1] = atof(tok_line[3].c_str());
     pNode->data.coords[3*i+2] = atof(tok_line[4].c_str());
     }
     count++;
   }
   
  if (count == 2) break;

  }
//cout << "Done reading output of mopac." << endl;

  output.close();
  return;
}

 void Mopac::xyz_read(string filename, vector<string> & anames, vector<double> & coords)
 {
    ifstream output(filename.c_str(),ios::in);

    if(output.fail())
    {
      cout << "Cannot read inpfile" << endl;
      exit(-1);
    }

    string line;
    vector<string> tok_line;
    int count = 0;
    int i = 0;
    while(!output.eof())
    {
    getline(output, line);

    if (count == 2)
    {
       for(size_t i=0; i < anames.size(); i++)
       {
       getline(output, line);
       tok_line = StringTools::tokenize(line, " ");
       coords[3*i+0] = atof(tok_line[2].c_str());
       coords[3*i+1] = atof(tok_line[3].c_str());
       coords[3*i+2] = atof(tok_line[4].c_str());
       }
       break;
    }

    if (line.find("CARTESIAN COORDINATES")!= string::npos && count == 0)
    {
      count++;
    }

    else if (line.find("CARTESIAN COORDINATES")!= string::npos && count == 1)
    {
     count++;
    }
    }

    output.close();
    return;
  }

/*
void Mopac::xyz_read(string mopac_folder,xyznode* pNode)
{
//  string filedir = "Work/mopac/";
  string oname = mopac_folder + "m" + pNode->comment + ".out";
  ifstream output(oname.c_str(),ios::in);
  
  if(output.fail())
  {
    cout << "Cannot read inpfile" << endl;
    exit(-1);
  }

  string line;
  vector<string> tok_line;
  int count = 0;
  int i = 0;
  while(!output.eof())
  {
  getline(output, line);

  if (count == 2)
  {  
     for(size_t i=0; i < pNode->natoms; i++)
     {
     getline(output, line);
     tok_line = StringTools::tokenize(line, " ");
     pNode->coords[3*i+0] = atof(tok_line[2].c_str());
     pNode->coords[3*i+1] = atof(tok_line[3].c_str());
     pNode->coords[3*i+2] = atof(tok_line[4].c_str());
     }
     break;
  }

  if (line.find("CARTESIAN COORDINATES")!= string::npos && count == 0)
  {
    count++;
  }
  
  else if (line.find("CARTESIAN COORDINATES")!= string::npos && count == 1)
  {
   count++;
  }
  }

  output.close();
  return;
}
*/     
int Mopac::check_success(string filename) 
{
//  string filedir = "/Work/mopac";
//cout << "test here" << endl;
   string oname = filename + ".out";
  ifstream output(oname.c_str(),ios::in);
  string line;
  if(output.fail())
  {
  cout << "Cannot open " << oname << endl;
  }

  string tmp;  
  while(getline(output, tmp))
  {
   line = tmp;
  }

  if (line.find("MOPAC DONE") != string::npos)
  {
  output.close();
  return 1;
  }
  else
  {
  output.close();
  return 0;
  }
 
} 

  void Mopac::mopac_genslurm(string filedir, vector<int> & list)
{

   ofstream slurmfile;
    string slurmfile_string = filedir + "go_slurm";

    slurmfile.open(slurmfile_string.c_str());
    cout << "Generating SBATCH  script for mopac jobs ..." << endl;

    slurmfile << "#!/bin/bash" << endl;
    slurmfile << "#SBATCH --array=";

    vector<int>::iterator it;
    it = list.begin();
    slurmfile << *it << "-";
    it = list.end()-1;
    slurmfile << *it << endl<< endl;

 /*Hard coded!!!*/
    slurmfile << "#SBATCH -N 1  --ntasks-per-node=1 --time=12:00:00" << endl;
    slurmfile << "#SBATCH -A zimmerman  " << endl;
 /*Hard coded!!*/
    slurmfile << "#SBATCH -p guest --job-name=go_slurm" << endl;
    slurmfile << "cp ~/MOPAC2016.exe /shaoyulu/" << endl;
    slurmfile << "export MOPAC_LICENSE=/home/paulzim/mopac" << endl;
    slurmfile << "export OMP_NUM_THREADS=1" << endl;
    slurmfile << ". /etc/profile.d/slurm.sh" << endl;
    slurmfile << endl;
    slurmfile << "ID=`printf \"%0*d\\n\" 4 ${SLURM_ARRAY_TASK_ID}`" << endl << endl;
    slurmfile << "cd " << filedir << endl << endl;
    slurmfile << "~/MOPAC2016.exe m$ID " << endl << endl;
    slurmfile << "wait" << endl << endl;
    slurmfile << "echo \"done with m$ID\" > mdone${SLURM_ARRAY_TASK_ID}" << endl << endl;
    slurmfile.close();

}


 void Mopac::mopac_genqsh(string filedir, vector<int> & list)
 {
    ofstream qshfile;
    string qshfile_string = filedir + "go_mopac.qsh";

    qshfile.open(qshfile_string.c_str());
    cout << "Generating PBS script for mopac jobs ..." << endl;

    qshfile << "#PBS -t ";

    vector<int>::iterator it;
    it = list.begin();
    qshfile << *it << "-";
    it = list.end()-1;
    qshfile << *it << endl<< endl;

 /*Hard coded!!!*/
    qshfile << "#PBS -l nodes=1:ppn=2 -l walltime=02:00:00" << endl;

    qshfile << "#PBS -e " << filedir << " -o " << filedir << " " << endl;
 /*Hard coded!!*/

//    qshfile << "#PBS -q guest    " << endl;
    qshfile << "#PBS -q zimmerman" << endl;
//    qshfile << "PBS -q parallel" << endl;
    
   qshfile << "cp ~/MOPAC2016.exe ~/" << endl;
   qshfile << "export MOPAC_LICENSE=/home/paulzim/mopac" << endl;
   qshfile << "export OMP_NUM_THREADS=2" << endl;
 
    qshfile << endl;
    qshfile << "ID=`printf \"%0*d\\n\" 4 ${PBS_ARRAYID}`" << endl << endl;

    qshfile << "cd " << filedir << endl << endl;
    qshfile << "~/MOPAC2016.exe m$ID " << endl << endl;

    qshfile << "wait" << endl << endl;

//    qshfile << "rm $QCSCRATCH" << endl << endl;

    qshfile << "echo \"done with m$ID\" > mdone${PBS_ARRAYID}" << endl << endl;

    qshfile.close();
 }

bool Mopac::run_opt_jobs(string filedir, vector<int> & list, int jobtype)
{
  
   string qshfile = "go_mopac.qsh";
  string slurmfile ="go_slurm";

  string cmdslurm = "sbatch " + filedir + slurmfile;
  string cmd = "qsub "+ filedir + qshfile;


  if(jobtype==1){
  system(cmdslurm.c_str());
}
else if(jobtype==2){ system(cmd.c_str());}



  size_t size = list.size();
  size_t done = 0;
  vector<int> mopacdone(size);
  vector<int>::iterator it;

  for (it=mopacdone.begin(); it!=mopacdone.end(); ++it)
  {
    *it = 0;
  }

  int max_wait = MAX_TIME_WAIT;
  int tc = 0;
  do
  {
    tc++;
    if (tc > max_wait)
    {
    cout << "Done waiting!" <<endl;
    break;
    }

  sleep(30);

    for (int i=0; i<size; i++)
    {
     if(!mopacdone[i])
     {
      string nstr = StringTools::int2str(list[i],1,"0");
      string dftfile_string = filedir + "mdone" + nstr;
      struct stat sts;
      if (stat(dftfile_string.c_str(), &sts) != -1)
      {
        mopacdone[i] = 1;
        cout  << "mopac" << nstr << " is done  " ;
        done++;
      }
     }
    }
  } while(done != size);
  cout << "Done with MOPAC calculation for this generation!" << endl;

  sleep(5);

  return true;

}  
