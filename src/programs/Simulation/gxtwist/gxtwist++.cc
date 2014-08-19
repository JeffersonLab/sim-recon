extern "C" int gxint_(void);
extern "C" int getwebfile(const char *url);

int main(int narg, char *argv[])
{
   // Make sure magnetic field file exists in local directory
   getwebfile("http://zeus.phys.uconn.edu/halld/tagger/simulation/TOSCA_tagger_dipole-15000G.map.gz");
   getwebfile("http://zeus.phys.uconn.edu/halld/tagger/simulation/TOSCA_tagger_quadrupole-nominal.map.gz");
   return gxint_();
}

