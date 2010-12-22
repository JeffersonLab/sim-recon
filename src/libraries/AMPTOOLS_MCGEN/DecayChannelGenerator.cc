
#include <iostream>
#include <math.h>
#include <stdlib.h>
#include <cassert>

#include "AMPTOOLS_MCGEN/DecayChannelGenerator.h"

DecayChannelGenerator::DecayChannelGenerator() :
m_bfTotal( 0 ),
m_upperBound( 0 ),
m_index( 0 ), 
m_probRenormalized( false )
{}

using namespace std;

void
DecayChannelGenerator::addChannel( unsigned int channel, double bf ) 
{
    
    m_bfTotal += bf;
    
    m_index.push_back( channel );
    m_prob.push_back( bf );
    m_upperBound.push_back( m_bfTotal );

    if( m_probRenormalized ){
        
        cout << "ERROR:  channels cannot be added after BF normalization!"
        << endl;
        
        assert( false );
    }
}

unsigned int
DecayChannelGenerator::operator()(){
 
    if( fabs( m_bfTotal - 1 ) > 0.001 ){
        
        cout << "WARNING:  sum of branching fractions: " << m_bfTotal 
        << "\n\t is not within .1% of 1 -- renormalizing" << endl;
    
        for( vector< double >::iterator val = m_upperBound.begin();
            val != m_upperBound.end();
            ++val ){
            
            (*val) /= m_bfTotal;
        }

        for( vector< double >::iterator val = m_prob.begin();
            val != m_prob.end();
            ++val ){
            
            (*val) /= m_bfTotal;
        }
        
        m_bfTotal = 1;
        m_probRenormalized = true;
    }
    
    double rand = drand48();
    for( unsigned int i = 0; i < m_upperBound.size(); ++i ){
        
        if( rand < m_upperBound[i] ){
            
            return m_index[i];
        }
    }
    
    return m_index[m_index.size()-1];
}

double
DecayChannelGenerator::getProb( unsigned int chan ){
 
    if( fabs( m_bfTotal - 1 ) > 0.001 ){
        
        cout << "WARNING:  sum of branching fractions: " << m_bfTotal 
        << "\n\t is not within .1% of 1 -- renormalizing" << endl;
        
        for( vector< double >::iterator val = m_upperBound.begin();
            val != m_upperBound.end();
            ++val ){
            
            (*val) /= m_bfTotal;
        }
        
        for( vector< double >::iterator val = m_prob.begin();
            val != m_prob.end();
            ++val ){
            
            (*val) /= m_bfTotal;
        }
        
        m_bfTotal = 1;
        m_probRenormalized = true;
    }
    
    return m_prob.at( chan );
}

