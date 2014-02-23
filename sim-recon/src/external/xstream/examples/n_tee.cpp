#include <xstream/tee.h>

#include <iostream>
#include <fstream>


using namespace std;;
using namespace xstream;;


int main(int argc, char* argv[]){

	try{
		if(1<argc){
			tee::ostreambuf Tee;
			for(int i=1; i<argc;++i){
				ofstream* file = new ofstream(argv[i]);
				Tee.add(file->rdbuf());
			}

			ostream tee(&Tee);
			
			//raise exceptions
			tee.exceptions(ios::badbit);

			tee<<cin.rdbuf();
			tee.flush();
			
		}else{
			cerr<<argv[0]<<" file1 file2 ... filen"<<endl;
		}
	}
	catch(exception& e){
		cerr<<"Error: "<<e.what()<<endl;
	}

	return 0;
}
