// Author: David Lawrence  Dec. 20, 2013
//
//
// hd_ana.cc
//

#include <fstream>

#include <DANA/DApplication.h>
#include <HDGEOMETRY/DGeometry.h>
using namespace std;

void Usage(JApplication &app);


//-----------
// main
//-----------
int main(int narg, char *argv[])
{
	// Instantiate an DApplication object
	DApplication app(narg, argv);

	//if(narg<=1)Usage(app);

	app.Init();
	DGeometry *dgeom = app.GetDGeometry(1);
	if(!dgeom){
		jerr << "Couldn't get DGeometry pointer!!" << endl;
		return -1;
	}

	// Get CDC wires
	vector<vector<DCDCWire *> > cdcwires;
	dgeom->GetCDCWires(cdcwires);

	// Get FDC wires
	vector<vector<DFDCWire *> > fdcwires;
	dgeom->GetFDCWires(fdcwires);
	
	cout << endl;
	cout << "Writing GlueX wire definitions to gluex_wires.txt ... " << endl;

	ofstream ofs("gluex_wires.txt");
	ofs << "#" << endl;
	ofs << "# CDC wires" <<endl;
	ofs << "# id  length pos_x pos_y pos_z dir_x dir_y dir_z" << endl;
	ofs << "#" << endl;

	ofs << "# ----- CDC ------" << endl;
	int Ncdc = 0;
	for(unsigned int i=0; i<cdcwires.size(); i++){
		vector<DCDCWire *> &wires = cdcwires[i];
		for(unsigned int j=0; j<wires.size(); j++){

			DCDCWire *w = wires[j];

			int id = 100000 + (i+1)*1000 + j;
			ofs << id << " " << w->L << " ";

			ofs << w->origin.X() << " ";
			ofs << w->origin.Y() << " ";
			ofs << w->origin.Z() << " ";

			ofs << w->udir.X() << " ";
			ofs << w->udir.Y() << " ";
			ofs << w->udir.Z() << " ";
			ofs << endl;
			Ncdc++;
		}
	}

	ofs << "#" << endl;
	ofs << "# ----- FDC ------" << endl;
	int Nfdc = 0;
	for(unsigned int i=0; i<fdcwires.size(); i++){
		vector<DFDCWire *> &wires = fdcwires[i];
		for(unsigned int j=0; j<wires.size(); j++){

			DFDCWire *w = wires[j];

			int id = 200000 + (i+1)*1000 + j;
			ofs << id << " " << w->L << " ";

			ofs << w->origin.X() << " ";
			ofs << w->origin.Y() << " ";
			ofs << w->origin.Z() << " ";

			ofs << w->udir.X() << " ";
			ofs << w->udir.Y() << " ";
			ofs << w->udir.Z() << " ";
			ofs << endl;
			Nfdc++;
		}
	}

	ofs.close();

	cout << "Done." << endl;
	cout << "Num. CDC wires: " << Ncdc << endl;
	cout << "Num. FDC wires: " << Nfdc << endl;


	return 0;
}

//-----------
// Usage
//-----------
void Usage(JApplication &app)
{
	cout<<endl;
	cout<<"Usage:"<<endl;
	cout<<"    hd_ana [options] source1 source2 source3 ..."<<endl;
	cout<<endl;
	app.Usage();
	cout<<endl;
	
	exit(0);
}

