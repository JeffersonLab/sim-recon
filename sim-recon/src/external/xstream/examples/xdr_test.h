#include <vector>
#include <string>

using std::vector;
using std::string;

struct test{
	int i;
	unsigned int ui;
	float f;
	double d;
	long long int li;
	unsigned long long int uli;
	string s;

	vector<int> vi;

	bool operator=(const test &t2);
};

struct test get_obj();
