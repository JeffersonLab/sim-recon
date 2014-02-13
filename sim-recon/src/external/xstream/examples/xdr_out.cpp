#include <xstream/xdr.h>
#include "xdr_test.h"

#include <iostream>
#include <stdexcept>

using namespace std;;
using namespace xstream;;

int main(int argc, char* argv[]){
	try{
		xdr::ostream xdr_o(cout);

		struct test o=get_obj();
		xdr_o<<o.i<<o.ui<<o.f<<o.d<<o.li<<o.uli<<o.s<<o.vi;;
	}
	catch(exception& e){
		cerr<<"Error: "<<e.what()<<endl;
	}

	return 0;
}
