#include <xstream/xdr.h>
#include "xdr_test.h"

#include <iostream>
#include <stdexcept>

using namespace std;;
using namespace xstream;;

int main(int argc, char* argv[]){
	try{
		xdr::istream xdr_i(cin);

		struct test o;
		xdr_i>>o.i>>o.ui>>o.f>>o.d>>o.li>>o.uli>>o.s>>o.vi;;
		cout<<o.i<<","<<o.ui<<","<<o.f<<","<<o.d<<","<<o.li<<","<<o.uli<<",["<<o.s<<"]"<<endl;
		size_t len = o.vi.size();
		cout<<"vector len = "<<len<<endl;
		while(len--!=0)
			cout<<o.vi[len]<<endl;
	}
	catch(exception& e){
		cerr<<"Error: "<<e.what()<<endl;
	}

	return 0;
}
