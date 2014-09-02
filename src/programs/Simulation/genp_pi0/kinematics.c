
/****************************************************************/
/* kinematics.c                                                 */
/*                                                              */
/* Some basic kinematics routines and tools                     */
/*                                                              */
/* D. Lawrence                                                  */
/* 3/18/99                                                      */
/****************************************************************/


#include <kinematics.h>
#include <stdio.h>



vect4 vect4_add(vect4 v1,vect4 v2)
{
   v1.E+=v2.E;
   v1.x+=v2.x;
   v1.y+=v2.y;
   v1.z+=v2.z;
   
   return v1;
}

vect4 vect4_sub(vect4 v1,vect4 v2)
{
   v1.E-=v2.E;
   v1.x-=v2.x;
   v1.y-=v2.y;
   v1.z-=v2.z;
   
   return v1;
}

double vect4_mul(vect4 v1,vect4 v2)
{
   v1.E*=v2.E;
   v1.x*=v2.x;
   v1.y*=v2.y;
   v1.z*=v2.z;
   
   return v1.E - v1.x - v1.y - v1.z;
}

double vect4_sq(vect4 v)
{
   return vect4_mul(v,v);
}

vect4 vect4_boost(vect4 p,double beta)
{
   vect4 u;
   double gamma=1.0/sqrt(1.0-(beta*beta));
   
   u.E=(p.E*gamma) + (p.z*gamma*beta);
   u.x=p.x;
   u.y=p.y;
   u.z=(p.E*gamma*beta) + (p.z*gamma);

   return u;
}

double vect4_mag2(vect4 v)
{
	return (v.x*v.x) + (v.y*v.y) + (v.z*v.z);
}

double vect4_mag(vect4 v)
{
	return sqrt(vect4_mag2(v));
}

double vect4_theta(vect4 v)
{
	return (180.0/M_PI)*acos(v.z/vect4_mag(v));
}

double vect4_phi(vect4 v)
{
	double phi;

	phi = atan2(v.y, v.x);
	if(phi<0.0)phi = 2.0*M_PI + phi; 
	return (180.0/M_PI)*phi;
}


/*****************************************************/
/* This enum is taken from cern.h since I don't want */
/* to have this file #including anything other than  */
/* kinematics.h                                      */
/*****************************************************/
/* GEANT Particle types */
enum geant_particles{
   ptype_none,
   ptype_gamma,
   ptype_positron,
   ptype_electron,
   ptype_neutrino,
   ptype_muon_plus,
   ptype_muon_minus,
   ptype_pion_zero,
   ptype_pion_plus,
   ptype_pion_minus,
   ptype_kaon_zero_long,
   ptype_kaon_plus,
   ptype_kaon_minus,
   ptype_neutron,
   ptype_proton,
   ptype_antiproton,
   ptype_kaon_zero_short,
   ptype_eta   
};

char* part_type_str(int type)
{
   static char str[256];
   switch(type){
   
      case ptype_gamma:
         return "photon";
         break;
      case ptype_positron:
         return "positron";
         break;
      case ptype_electron:
         return "electron";
         break;
      case ptype_neutrino:
         return "neutrino";
         break;
      case ptype_muon_plus:
         return "mu+";
         break;
      case ptype_muon_minus:
         return "mu-";
         break;
      case ptype_pion_zero:
         return "pi0";
         break;
      case ptype_pion_plus:
         return "pi+";
         break;
      case ptype_pion_minus:
         return "pi-";
         break;
      case ptype_neutron:
         return "neutron";
         break;
      case ptype_proton:
         return "proton";
         break;   
      default:
         sprintf(str,"%d",type);
         return str;
   }
}

int chargeof(int type)
{
	if(type==ptype_none           )return  0;
	if(type==ptype_gamma          )return  0;
	if(type==ptype_positron       )return  1;
	if(type==ptype_electron       )return -1;
	if(type==ptype_neutrino       )return  0;
	if(type==ptype_muon_plus      )return  1;
	if(type==ptype_muon_minus     )return -1;
	if(type==ptype_pion_zero      )return  0;
	if(type==ptype_pion_plus      )return  1;
	if(type==ptype_pion_minus     )return -1;
	if(type==ptype_kaon_zero_long )return  0;
	if(type==ptype_kaon_plus      )return  1;
	if(type==ptype_kaon_minus     )return -1;
	if(type==ptype_neutron        )return  0;
	if(type==ptype_proton         )return  1;
	if(type==ptype_antiproton     )return -1;
	if(type==ptype_kaon_zero_short)return  0;
	if(type==ptype_eta            )return  0;      
}




