// Author: Simon Taylor 
//
//
// bfield_finemesh.cc
//
// Program to create a fine-mesh map of the magnetic field (including 
// gradients) in the form of an evio file using the coarse-mesh map computed
// by TOSCA/Ansys/.. as input.

#include "DANA/DApplication.h"
using namespace std;

#include <vector>
#include <HDGEOMETRY/DMagneticFieldMap.h>
#include <evioFileChannel.hxx>
#include <evioUtil.hxx>
using namespace evio;

void Usage(void);

//-----------
// main
//-----------
int main(int narg, char *argv[])
{
  if(narg<=1)Usage();
    
  double rmax=100.;
  double rmin=0.;
  double zmax=600.;
  double zmin=0.;
  double dr=0.1;
  double dz=0.1;
  string evioFileName="finemesh.evio";

  // parse command line arguments
  for(int i=1; i<narg; i++){
    string arg(argv[i]);
    string next(i<narg-1 ? argv[i+1]:"");
    float argf = atof(next.c_str());
    bool used_next = false; 
    if(arg=="-zmin"){used_next=true; zmin = argf;}
    if(arg=="-zmax"){used_next=true; zmax = argf;}
    if(arg=="-rmin"){used_next=true; rmin = argf;}
    if(arg=="-rmax"){used_next=true; rmax = argf;}
    if(arg=="-dr"){used_next=true;dr=argf;}
    if(arg=="-dz"){used_next=true;dz=argf;}  
    if(arg=="-file"){evioFileName=next;};
    if(used_next){
      // skip to next argument
      i++;
      
      //If i is now > than narg, then no value was passed to an "-XXX" argument
      if(i>narg){
	_DBG_<<"No argument given for \""<<arg<<"\" argument!"<<endl;
	Usage();
	exit(-1);
      }
    }    
  }
  cout << "Generation fine-mesh B-field map with the following parameters:" 
       <<endl;
  cout << "   rmin = " << rmin <<endl;
  cout << "   rmax = " << rmax << endl;
  cout << "   dr   = " << dr <<endl;
  cout << "   zmin = " << zmin << endl;
  cout << "   zmax = " << zmax << endl;
  cout << "   dz   = " << dz <<endl;

  // Instantiate an event loop object
  DApplication app(narg, argv);
  
  // Initialize
  app.Init();

  // Read the Ansys/TOSCA/Poisson map from CalibDB
  DMagneticFieldMap *bfield = app.GetBfield();

  // Determine the number of points in x and z from the input parameters
  unsigned int Nr=(unsigned int)floor((rmax-rmin)/dr+0.5);
  unsigned int Nz=(unsigned int)floor((zmax-zmin)/dz+0.5);

  // Create the finer-mesh map 
  vector<float>Br_;
  vector<float>Bz_;
  vector<float>dBrdr_;
  vector<float>dBrdz_;  
  vector<float>dBzdr_;
  vector<float>dBzdz_;
  for (unsigned int i=0;i<Nr;i++){
    double r=rmin+dr*double(i);
    for (unsigned int j=0;j<Nz;j++){
      double z=zmin+dz*double(j);
      double Bx,By,Bz;
      double dBxdx,dBxdy,dBxdz,dBydx,dBydy,dBydz,dBzdx,dBzdy,dBzdz;
      
      bfield->GetFieldAndGradient(r,0.,z,Bx,By,Bz,dBxdx,dBxdy,dBxdz,
				  dBydx,dBydy,dBydz,dBzdx,dBzdy,dBzdz);
      
      Br_.push_back(Bx);
      Bz_.push_back(Bz);
      dBrdr_.push_back(dBxdx); 
      dBrdz_.push_back(dBxdz);
      dBzdr_.push_back(dBzdx); 
      dBzdz_.push_back(dBzdz);
    }
  }

  // Open the evio file channel
  unsigned long bufsize=Nr*Nz*6*sizeof(float)+6;
  evioFileChannel chan(evioFileName,"w",bufsize);
  chan.open();
  
  // create an event tree, root node has (tag=1,num=0)
  evioDOMTree tree(1,0);

  float minmaxdelta[6]={rmin,rmax,dr,zmin,zmax,dz};
  tree.addBank(2,0,minmaxdelta,6);

  // Add the banks
  tree.addBank(3,0,Br_);
  tree.addBank(3,1,Bz_);
  tree.addBank(3,2,dBrdr_);
  tree.addBank(3,3,dBrdz_);
  tree.addBank(3,4,dBzdr_);
  tree.addBank(3,5,dBzdz_);

  chan.write(tree);
  chan.close();

  return 0;
}

//-----------
// Usage
//-----------
void Usage(void)
{
	cout<<endl;
	cout<<"Usage:"<<endl;
	cout<<"   bfield_finemesh -rmin [rmin] -rmax [rmax] -dr [dr] -zmin [zmin] -zmax [zmax] -dz [dz] [-file filename]"<<endl;
	cout<<endl;
	
	exit(0);
}

