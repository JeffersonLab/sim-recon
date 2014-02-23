/*
 * Base64 encode data
 *
 * without parameters reads data from standard input and writes encoded data to standard output
 * with 1 parameter reads data from the file given and writes encoded data to standard output
 * with 2 parameters reads data from the file in the 1st argument and writes encoded data to the file in the 2nd argument
 *
 */

#include <xstream/base64.h>

#include <iostream>
#include <fstream>

using namespace std;
using namespace xstream;


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

		base64::ostreambuf b64sb(writeto->rdbuf());
		ostream b64(&b64sb);
		
		//raise exceptions
		b64.exceptions(ios::badbit);


		while(readfrom->good()){
			readfrom->read(buf,len);
			b64.write(buf,readfrom->gcount());
		}

		b64.flush();

	}
	catch(exception& e){
		cerr<<"Error: "<<e.what()<<endl;
	}


	//pointers not being freed

	return 0;
}
