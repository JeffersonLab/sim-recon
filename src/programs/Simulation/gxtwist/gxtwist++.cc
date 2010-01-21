


extern "C" int gxint_(void);
extern "C" int getwebfile(const char *url);

int main(int narg, char *argv[])
{
	// Make sure magnetic field file exists in local directory
	const char *url = "http://zeus.phys.uconn.edu/halld/tagger/simulation/taggerBfield-quad-map.gz";
	getwebfile(url);

	return gxint_();
}

