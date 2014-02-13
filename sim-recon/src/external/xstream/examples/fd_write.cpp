#include <xstream/fd.h>

#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <stdio.h>
#include <errno.h>

#include <iostream>

using namespace std;
using namespace xstream;


int main(int argc, char* argv[]){

	try{
		int fd;
		if(2==argc){
			fd = open(argv[1], O_WRONLY|O_CREAT);
			if(-1==fd){
				perror("open");
				return errno;
			}else{
				fd::streambuf fs(fd,true);
				ostream os(&fs);
				os.exceptions(ios::badbit);

				os<<cin.rdbuf();
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
