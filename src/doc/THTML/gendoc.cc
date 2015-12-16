
#include <stdlib.h>

#include <THtml.h>

int main(int argc, char *argv[])
{
	// Create THtml object and generate documentation
	THtml html;
	//html.SetInputDir("JANA:/usr/local/root/PRO/include");
	html.MakeAll();
	
	return 0;
}
