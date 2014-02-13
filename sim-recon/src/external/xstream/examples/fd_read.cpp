#include <xstream/fd.h>

#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <stdio.h>
#include <errno.h>

#include <iostream>
#include <string>

using namespace std;
using namespace xstream;


int main(int argc, char* argv[]){

	try{
		int fd;
		if(2==argc){
			fd = open(argv[1], O_RDONLY);
			if(-1==fd){
				perror("open");
				return errno;
			}else{

				string line;

				fd::streambuf fs(fd,true);
				istream is(&fs);
				is.exceptions(ios::badbit);

				while(is.good()){
					getline(is,line);
					cout<<line<<endl;
				}
			}
		}else{
			cerr<<argv[0]<<" file"<<endl;
		}
	}
	catch(exception& e){
		cerr<<"Error: "<<e.what()<<endl;
	}

	return 0;
}
