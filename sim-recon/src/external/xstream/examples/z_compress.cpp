/*
 * Compresses data in Zlib format
 *
 * without parameters reads data from standard input and writes compressed data to standard output
 * with 1 parameter reads data from the file given and writes compressed data to standard output
 * with 2 parameters reads data from the file in the 1st argument and writes compressed data to the file in the 2nd argument
 *
 */

#include <xstream/z.h>

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
		/* Binary flags mandatory on win32 systems, otherwise files will be corrupt */
		if(1<argc){
			readfrom = new ifstream(argv[1], ios::binary);
			if(2<argc){
				writeto = new ofstream(argv[2], ios::binary);
			}else{
				writeto = &cout;
			}
		}else{
			readfrom = &cin;
			writeto  = &cout;
		}

		z::ostreambuf z_o(writeto->rdbuf());
		ostream zout(&z_o);

		//raise exceptions
		zout.exceptions(ios::badbit);

		while(readfrom->good()){
			readfrom->read(buf,len);
			zout.write(buf,readfrom->gcount());
		}

		zout.flush();
		//clog<<"Adler32 = "<< (z_o.checksum() )<<endl;

	}
	catch(exception& e){
		cerr<<"Error: "<<e.what()<<endl;
	}


	//pointers not being freed

	return 0;
}
