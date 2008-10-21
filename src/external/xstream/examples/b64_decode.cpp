/*
 * Base64 decode data
 *
 * without parameters reads encoded data from standard input and writes decoded data to standard output
 * with 1 parameter reads encoded data from the file given and writes decoded data to standard output
 * with 2 parameters reads encoded  data from the file in the 1st argument and writes decoded data to the file in the 2nd argument
 *
 */

#include <xstream/base64.h>

#include <iostream>
#include <fstream>

using namespace std;;
using namespace xstream;;


int main(int argc, char* argv[]){

	const size_t len = 4*1024;
	char buf[len];

	istream* readfrom;
	ostream* writeto;

	try{
		if(1<argc){
			readfrom = new ifstream(argv[1]);
			if(2<argc){
				writeto = new ofstream(argv[2]);
			}else{
				writeto = &cout;
			}
		}else{
			readfrom = &cin;
			writeto  = &cout;
		}

		base64::istreambuf b64sb(readfrom->rdbuf());
		istream b64(&b64sb);
		
		//raise exceptions
		b64.exceptions(ios::badbit);

		while(b64.good()){
			b64.read(buf,len);
			writeto->write(buf,b64.gcount());
		}

		b64.sync();
		writeto->flush();

	}

	catch(exception& e){
		cerr<<"Error: "<<e.what()<<endl;
	}


	//pointers not being freed

	return 0;
}
