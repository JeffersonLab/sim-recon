	If you're reading this, you must be attempting to develop some tracking code for the FDC.
Since I was an undergraduate when I began this project, I'm writing this for a C++-literate undergraduate. I also assume you've looked through the paper I wrote for this project that 
gives a good introduction to the chambers and my strategy for reconstruction, found at
<http://www.jlab.org/~craigb/HallD/FDCreport.pdf>. That being said, since it was my strategy, it's totally up for reconsideration.
	To start, HDGeant creates FDC hits (along with hits in all of the other detectors) and
stores them in an HDDM file. Currently, the data model HDGeant uses to place hits in the FDC 
doesn't use wire and strip elements that detect collisions--it calculates a wire or strip number
based on the X-Y position of the track when it reaches the Z position of a layer. Coincidentally, the data model does not save the wire/strip number it calculates to the data
file when it saves the "truth" point. This must be fixed--truth points must include wire/strip
number and layer information to plot the difference between truth points and reconstructed points. Also, there are bugs in the data model that cause analysis programs to crash after
a random number of events (always the same event for a given data file), so don't attempt
to debug those events by questioning your DANA code. It is most likely something askew in
HDGeant.  
	Next, in a DANA program like hd_ana, the DFactory_DFDCHit class will extract the FDC data 
created by HDGeant from the data file supplied on the command line. DFactory_DFDCHit is relatively complete, but any changes to the data model or the DFDCHit class will require equivalent changes in DFactory_DFDCHit. 
	From there, DFactory_FDCCathodeCluster takes the hit data, sorts out the cathode hits, and
then associates them into clusters. The algorithm right now merely looks for consecutive strip
numbers to place them into a cluster. It's feasible and even likely that clusters will overlap, and that peaks will need to be distinguished with a sophisticated peak-finding algorithm. If one
can ascertain the position of the peak, that can be used in tracking to improve the resolution.
	After clusters are made, the DFactory_DFDCPseudo class attempts to intersect the clusters
and reconstruct points by using anode hits to thin out the cathode candidates. Right now, this
doesn't work. It is true that required an anode hit makes us "picky" about our points, but we should be finding at least some sizable fraction of points placed by the simulator, and we're finding none. So, this should be the first order of business for anyone adopting this project--
get this part working and have a look at the quality of the reconstructed points.
	After that works, drift time will certainly be incorporated for the anode wires. HDGeant must use some distance-to-time relationship to calculate its simulated drift times, so implementation of the reverse of that relationship to calculate drift distance will definitely improve resolution. 
	So in summary, I humbly suggest these improvements (in order):
	1) Implement the changes to HDGeant that would add wire/strip and layer numbers to the
	truth points.
	2) Get the intersection algorithm working at a respectable rate, and get an idea of the
	resolution using no drift and no peak-finding.
	3) Implement drift time relationships. 
	4) Implement a peak-finding algorithm for the cathode strip clusters.

Have fun!

- Craig
