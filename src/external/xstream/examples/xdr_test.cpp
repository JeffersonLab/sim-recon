#include "xdr_test.h"

bool test::operator=(const test &t2){
	if(!(
		i==t2.i && ui==t2.ui && f==t2.f && d==t2.d && li==t2.li && uli==t2.uli && s == t2.s
	)){
		return false;
	}else{
		size_t len=vi.size();
		if(len!=t2.vi.size()){
			return false;
		}else{
			for(size_t i=0; i<len; ++i)
				if(vi[i]!=t2.vi[i])
					return false;
		}
	}
	return true;
}

struct test get_obj(void){
	test t;
	t.i=-512;
	t.ui=512;
	t.f=-3.1415;
	t.d=-3.1415;
	t.li=-123456789;
	t.uli=123456789;

	t.s="This is a string 1";

	t.vi.push_back(-1);
	t.vi.push_back(-2);
	t.vi.push_back(-3);
	t.vi.push_back(-4);
	t.vi.push_back(-5);

	return t;
}

