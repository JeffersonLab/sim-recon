#if !defined(DECAYCHANNELGENERATOR)
#define DECAYCHANNELGENERATOR

#include <vector>

using namespace std;

class DecayChannelGenerator {
    
public:
    
    DecayChannelGenerator();
    
    void addChannel( unsigned int channelNum, double bf );
    unsigned int operator()();
    
    const vector< unsigned int >& availableChannels() const { return m_index; }
    
    double getProb( unsigned int channelNum );
    
private:
    
    double m_bfTotal;
    
    vector< double > m_upperBound;
    vector< double > m_prob;
    vector< unsigned int > m_index;
    
    bool m_probRenormalized;
};

#endif
