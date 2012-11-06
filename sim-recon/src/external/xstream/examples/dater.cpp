#include <xstream/dater.h>

#include <iostream>
#include <fstream>


using namespace std;;
using namespace xstream;;


int main(int argc, char* argv[]){
	const size_t len = 4*1024;
	char buf[len];

	ostream* out;

	try{
		if(2==argc){
			out = new ofstream(argv[1]);
		}else{
			out=&cout;
		}

		dater dsb(out->rdbuf());
		//dater dsb(out->rdbuf(),"{%F %T %z}-> ");
		ostream dout(&dsb);

		//raise exceptions
		dout.exceptions(ios::badbit);

		while(cin.good()){
			cin.read(buf,len);
			dout.write(buf,cin.gcount());
		}

		dout.flush();

	}
	catch(exception& e){
		cerr<<"Error: "<<e.what()<<endl;
	}


	//pointers not being freed

	return 0;
}
