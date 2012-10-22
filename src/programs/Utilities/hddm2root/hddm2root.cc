

#include <iostream>
#include <iomanip>
#include <string>
#include <fstream>
#include <stack>
using namespace std;

#include <sys/stat.h>
#include <time.h>
#include <expat.h>

#include "DClassDef.h"

string XML_FILENAME = "/Users/davidl/HallD/builds/sim-recon_latest/src/libraries/HDDM/event.xml";
string HDDM_CLASS = "x";

map<string, DClassDef> CLASSES;
unsigned int MAX_DEPTH=0;

void startElement(void *data, const char *el, const char **attr);
void endElement(void *data, const char *el);
void CreateCopyRoutines(void);
string NameToHDDMname(string name);
string HDDMPlural(string name_hddm);

//--------------
// main
//--------------
int main(int narg, char *argv[])
{
	// Parse XML file to get class definitions
	stack<string> ancestory;

	XML_Parser parser = XML_ParserCreate(NULL);
	XML_SetUserData(parser, &ancestory);
	XML_SetElementHandler(parser, startElement, endElement);

	ifstream ifs(XML_FILENAME.c_str());

	// get length of file:
	ifs.seekg (0, ios::end);
	int length = ifs.tellg();
	ifs.seekg (0, ios::beg);

	char *buff = new char[length];
	ifs.read (buff,length);
	ifs.close();
	
	int done = 0;
	XML_Parse(parser, buff, length, done);

	delete[] buff;


	// Create a header file for each class. This is needed to make
	// the ROOT dictionaries (cint doesn't seem to like it when they
	// are all in the same file).
	
	// Put all headers in sub-directory just to keep things a little cleaner
	mkdir("hddm_root_generated", S_IRWXU | S_IRWXG | S_IRWXO);
	
	// Print class definitions
	// Here we must print them in reverse-depth order so the classes
	// that show up as members of others are defined first.
	map<string, DClassDef>::iterator iter;
	for(iter=CLASSES.begin(); iter!=CLASSES.end(); iter++){
		string name = iter->first;
		string type = iter->second.name;
		map<string,string> &members = iter->second.members;

		//if(iter->second.depth != depth)continue; // printing in reverse-depth order
		if(name == "HDDM_t")continue;  // Skip outer HDDM tags
		
		// Open output file
		string fname = string("hddm_root_generated/") + name + ".h";
		cout<<"writing "<<fname<<endl;
		ofstream ofs(fname.c_str());

		// Write header to file
		ofs << "#include<vector>" << endl;
		ofs << "using namespace std;" << endl;
		ofs << endl;
		ofs << "#include<RTypes.h>" << endl;
		ofs << "#include<TObject.h>" << endl;
		ofs << endl;
		
		ofs << "#ifndef _" << name << "_"<<endl;
		ofs << "#define _" << name << "_"<<endl;
		
		ofs << "typedef int Particle_t;" << endl;
		ofs << endl;
		
		// Add relevant includes for other classes this depends on
		set<string> &include_types = iter->second.include_types;
		set<string>::iterator iter_inc;
		for(iter_inc=include_types.begin(); iter_inc!=include_types.end(); iter_inc++){
			ofs << "#include \"" << *iter_inc << "_t.h\""<<endl;
		}
		ofs << endl;
		
		// Write class definition to file
		ofs<<"class " << name << ":public TObject{" << endl;
		ofs<<"	public:"<<endl;
		ofs<<endl;
		
		map<string,string>::iterator iter2;
		for(iter2=members.begin(); iter2!=members.end(); iter2++){
			ofs << "		" << iter2->second << " " << iter2->first << ";" << endl;
		}
		ofs<<endl;
		ofs<<"		ClassDef("<<name<<",1)"<<endl;
		ofs<<"};"<<endl;
		ofs<<endl;

		ofs << "#endif" << endl;
		
		// Close file
		ofs.close();
	}
	
	// Create script that will generate ROOT dictionaries for each header file
	ofstream ofs("hddm_root_generated/mk_dictionaries.csh");
	ofs << "#!/bin/tcsh -f "<<endl;
	ofs << endl;
	for(iter=CLASSES.begin(); iter!=CLASSES.end(); iter++){
		string name = iter->first;
		if(name == "HDDM_t")continue;  // Skip outer HDDM tags

		ofs << "echo making ROOT dictionary for " << name << " ..." << endl;
		string cmd = "rootcint -f " + name + "_Dict.cc -c " + name + ".h";
		ofs << cmd << endl;
	}
	ofs << endl;

	for(iter=CLASSES.begin(); iter!=CLASSES.end(); iter++){
		string name = iter->first;
		if(name == "HDDM_t")continue;  // Skip outer HDDM tags

		ofs << "echo compiling ROOT dictionary for " << name << " ..." << endl;
		string cmd = "c++ `root-config --cflags` -c " + name + "_Dict.cc";
		ofs << cmd << endl;
	}
	ofs << endl;
	
	ofs << "ar -r libhddm_root.a *_t_Dict.o" << endl;
	ofs << endl;
	
	ofs << "echo compiling hddm_root_CopyRoutines.cc ..." << endl;
	string cmd = "c++ `root-config --cflags` -c -I${HALLD_HOME}/include -I. hddm_root_CopyRoutines.cc";
	ofs << cmd << endl;
	ofs << endl;

	ofs << "ar -r libhddm_root.a hddm_root_CopyRoutines.o" << endl;
	ofs << endl;

	ofs << "echo compiling hddm_root_main.cc ..." << endl;
	cmd = "c++ `root-config --cflags --libs` -I${HALLD_HOME}/include -I. hddm_root_main.cc ./libhddm_root.a -lbz2 -lz -L${HALLD_HOME}/lib/${BMS_OSNAME} -lHDDM";
	ofs << cmd << endl;
	ofs << endl;

	ofs.close();
	
	// Make the dictionary generator script executable
	chmod("hddm_root_generated/mk_dictionaries.csh", S_IRWXU | S_IRWXG | S_IRWXO);
	
	// Create Copy Routines
	CreateCopyRoutines();
	
	return 0;
}

//--------------
// startElement
//--------------
void startElement(void *data, const char *el, const char **attr)
{
	// Setup reference to our ancestory stack
	stack<string> &ancestory = *(stack<string>*)data;
	
	// Convert element tag to string
	string name(el);
	bool is_HDDM_tag = (name == "HDDM");
	
	// Loop over attributes, temporarily storing them into map
	map<string, string> members;
	bool is_unbounded = false;
	for (int i = 0; attr[i]; i += 2) {
		string name(attr[i]);
		string type(attr[i+1]);
		
		// Some attributes are HDDM-specific
		if(name == "class") HDDM_CLASS = type;
		if(name == "minOccurs")continue;
		if(name == "version")continue;
		if(name == "maxOccurs"){
			if(type == "unbounded"){
				is_unbounded = true;
				continue;
			}
			if(atoi(type.c_str())>1){
				is_unbounded = true;
				continue;
			}
			continue;
		}
		
		if(type=="boolean") type = "bool"; // (sigh...)
		
		// Add member to list
		members[name] = type;
	}

	
	// If there is a parent, add this element as a member
	if(!ancestory.empty()){
		if(is_unbounded){
			CLASSES[ancestory.top()].members[name+"s"] = string("vector<") + name + "_t>";
		}else{
			CLASSES[ancestory.top()].members[name] = name + "_t";
		}
		CLASSES[ancestory.top()].include_types.insert(name);
	}

	// Get reference to class definition, creating it if necessary
	DClassDef &class_def = CLASSES[name + "_t"];
	class_def.members.insert(members.begin(), members.end());
	class_def.depth = ancestory.size();
	if(class_def.depth > MAX_DEPTH) MAX_DEPTH = class_def.depth;

	ancestory.push(name + "_t");
}

//--------------
// endElement
//--------------
void endElement(void *data, const char *el)
{
	stack<string> &ancestory = *(stack<string>*)data;
	
	ancestory.pop();
}

//--------------
// CreateCopyRoutines
//--------------
void CreateCopyRoutines(void)
{
	// Open header file to write into
	cout << "Writing hddm_root_CopyRoutines.h" << endl;
	ofstream ofs("hddm_root_generated/hddm_root_CopyRoutines.h");

	// Add file header
	time_t t = time(NULL);
	ofs << "//" << endl;
	ofs << "// Auto-generated HDDM to ROOT copy routines DO NOT EDIT" << endl;
	ofs << "//" << endl;
	ofs << "// "<<ctime(&t);
	ofs << "// Generated from: "<<XML_FILENAME<<endl;
	ofs << "//" << endl;
	ofs << endl;
	
	ofs << "#define Particle_t Particle_tt"<<endl;
	map<string, DClassDef>::iterator iter;
	for(int depth=0; depth<MAX_DEPTH; depth++){
		for(iter=CLASSES.begin(); iter!=CLASSES.end(); iter++){
			if(iter->second.depth != depth) continue;
			if(iter->first == "HDDM_t")continue;
			ofs << "#include \"" << iter->first << ".h\"" << endl;
		}
	}
	ofs << "#undef Particle_t" << endl;
	ofs << endl;
	ofs << endl;
	
	ofs << "#include <hddm_"<<HDDM_CLASS<<".hpp>"<<endl;
	ofs << "using namespace hddm_"<<HDDM_CLASS<<";"<<endl;
	ofs << endl;

	// Copy Routine declarations
	for(iter=CLASSES.begin(); iter!=CLASSES.end(); iter++){
		string name = iter->first;
		string varname = name.substr(0, name.length()-2);
		string name_hddm = NameToHDDMname(varname);

		if(varname == "HDDM")continue;
		
		ofs << "void Copy"<<name_hddm<<"("<<name<<" &"<<varname<<", class "<<name_hddm<<" &"<<varname<<"_hddm, bool hddm_invalid=false);"<<endl;
	}
	ofs << endl;
	ofs << endl;
	ofs.close();  // header file
	
	// Open implementation file to write into
	cout << "Writing hddm_root_CopyRoutines.cc" << endl;
	ofs.open("hddm_root_generated/hddm_root_CopyRoutines.cc");
	
	// Add file header
	ofs << "//" << endl;
	ofs << "// Auto-generated HDDM to ROOT copy routines DO NOT EDIT" << endl;
	ofs << "//" << endl;
	ofs << "// "<<ctime(&t);
	ofs << "// Generated from: "<<XML_FILENAME<<endl;
	ofs << "//" << endl;
	ofs << endl;
	ofs << "#include \"hddm_root_CopyRoutines.h\"" << endl;
	ofs << endl;
	ofs << endl;

	// Copy Routines
	for(iter=CLASSES.begin(); iter!=CLASSES.end(); iter++){
		string name = iter->first;
		string varname = name.substr(0, name.length()-2);
		string name_hddm = NameToHDDMname(varname);
		string name_hddm_plural = HDDMPlural(name_hddm);
		map<string, string> &members = iter->second.members;
		
		if(varname == "HDDM")continue;
		
		ofs << "void Copy"<<name_hddm<<"("<<name<<" &"<<varname<<", class "<<name_hddm<<" &"<<varname<<"_hddm, bool hddm_invalid)"<<endl;
		ofs << "{"<<endl;
		ofs << "	if(hddm_invalid)return; // FIXME!!!" << endl;
		ofs << endl;
		
		map<string, string>::iterator iter_members;
		for(iter_members=members.begin(); iter_members!=members.end(); iter_members++){
			string myname = iter_members->first;
			string mytype = iter_members->second;
			string myvarname = myname;
			string myname_hddm = NameToHDDMname(myvarname);
			
			if(   mytype == "int"
			   || mytype == "float"
			   || mytype == "bool"
			   || mytype == "Particle_t"
			   || mytype == "string"){
			   
			   // Atomic type uses equals
				ofs << "	" << varname<<"."<<myvarname<<" = "<<varname<<"_hddm.get"<<myname_hddm<<"();" << endl;
			}else if(mytype.find("vector<") == 0){
				// STL vector requires iteration
				mytype = myname.substr(0, myname.length()-1); // chop off "s" that was added
				myname_hddm = NameToHDDMname(mytype);
				string myname_hddm_plural = HDDMPlural(mytype);
				string iter_name = "iter_" + mytype;
				ofs << endl;
				//ofs << "--- skipping "<<mytype<<":"<<myname<<":"<<myname_hddm_plural<<endl;
				ofs << "	"<<myname_hddm<<"List &"<<myname_hddm_plural<<" = "<<varname<<"_hddm.get"<<HDDMPlural(myname_hddm)<<"();"  <<endl;
				ofs << "	"<<myname_hddm<<"List::iterator "<<iter_name<<";"<<endl;
				ofs << "	for("<<iter_name<<"="<<myname_hddm_plural<<".begin(); "<<iter_name<<"!="<<myname_hddm_plural<<".end(); "<<iter_name<<"++){"<<endl;
				ofs << "		"<<mytype<<"_t a;"<<endl;
				ofs << "		Copy"<<myname_hddm<<"(a, *"<<iter_name<<");"<<endl;
				ofs << "		"<<varname<<"."<<myname<<".push_back(a);"<<endl;
				ofs << "	}"<<endl;
			}else{
				// Non-vector calls copy routine
				string myname_hddm_plural = HDDMPlural(myname_hddm);
				if(myname_hddm_plural == "Properties") myname_hddm_plural = "PropertiesList";
				ofs << "	Copy"<<myname_hddm<<"("<<varname<<"."<<myvarname<<", "<<varname<<"_hddm.get"<<myname_hddm<<"(), "<<varname<<"_hddm.get"<< myname_hddm_plural<< "().empty());" << endl;
			}
		}
		
		ofs << "}" << endl;
		ofs << endl;
	}
	ofs << endl;
	ofs << endl;
	
	ofs.close();
}

//--------------
// NameToHDDMname
//--------------
string NameToHDDMname(string name)
{
	// Convert our class name used for ROOT
	// into the HDDM class name
	string name_hddm = name;
	name_hddm[0] = toupper(name[0]);
	
	return name_hddm;
}

//--------------
// HDDMPlural
//--------------
string HDDMPlural(string name_hddm)
{
	// Re-create HDDM-style plurals

	if(name_hddm.size()<3) return name_hddm + "s";

	if(name_hddm.substr(name_hddm.length()-2, 2) == "ex"){
		name_hddm.erase(name_hddm.length()-2);
		name_hddm += "ices";
	}else if(name_hddm.substr(name_hddm.length()-4, 4) == "ties"){
		// do nothing
	}else if(name_hddm.substr(name_hddm.length()-2, 2) == "ty"){
		name_hddm.erase(name_hddm.length()-2);
		name_hddm += "ties";
	}else if(name_hddm.substr(name_hddm.length()-3, 3) == "tum"){
		name_hddm.erase(name_hddm.length()-3);
		name_hddm += "ta";
	}else if(name_hddm.substr(name_hddm.length()-1, 1) == "s"){
		name_hddm += "es";
	}else{
		name_hddm += "s";
	}
	
	return name_hddm;
}


