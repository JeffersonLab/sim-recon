/*
 * Decompresses Bz2 data
 *
 * without parameters reads compressed data from standard input and writes decompressed data to standard output
 * with 1 parameter reads compressed data from the file given and writes decompressed data to standard output
 * with 2 parameters reads compressed data from the file in the 1st argument and writes decompressed data to the file in the 2nd argument
 *
 */

#include <xstream/bz.h>

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
		/* Binary flags mandtory on win32 systems, otherwise files will be corrupt */
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

		bz::istreambuf bz_i(readfrom->rdbuf());
		istream bzin(&bz_i);

		//raise exceptions
		bzin.exceptions(ios::badbit);

		while(bzin.good()){
			bzin.read(buf,len);
			writeto->write(buf,bzin.gcount());
		}

	}
	catch(exception& e){
		cerr<<"Error: "<<e.what()<<endl;
	}

	writeto->flush();

	//pointers not being freed

	return 0;
}
