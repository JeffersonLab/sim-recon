#include "GPDs.hh"

GPDs::GPDs(const char *file_name, int nn_q2, int nn_t, int nn_eta, double q2_, double t_, double eta_)
{
  fname = (char*)file_name;
  n_q2 = nn_q2;
  n_t = nn_t;
  n_eta = nn_eta;
  Q2 = q2_;
  tM = -t_; // t -s real t, tM is -t
  eta = eta_;
  
  arr_q2  = new double[n_q2];
  arr_t   = new double[n_t];
  arr_eta = new double[n_eta];
  
  ImH_2 = new double[n_q2*n_t*n_eta];
  ReH_2 = new double[n_q2*n_t*n_eta];
  ImE_2 = new double[n_q2*n_t*n_eta];
  ReE_2 = new double[n_q2*n_t*n_eta];
  ImHtild_2 = new double[n_q2*n_t*n_eta];
  ReHtild_2 = new double[n_q2*n_t*n_eta];
  Dterm_2 = new double[n_q2*n_t*n_eta];
  
  ReadFile();
  DefineValues();
}

bool GPDs::ReadFile()
{
  ifstream inp(fname);
  
  if( inp.is_open() )
    {
      for( int i = 0; i < n_q2; i++ )
	{
	  for( int j = 0; j < n_t; j++ )
	    {
	      for( int k = 0; k < n_eta; k++ )
		{
		  inp>>arr_q2[i];
		  inp>>arr_t[j];
		  arr_t[j] = -arr_t[j];
		  inp>>arr_eta[k];
		  inp>>ImH_2[i*n_t*n_eta + j*n_eta + k];
		  inp>>ReH_2[i*n_t*n_eta + j*n_eta + k];
		  inp>>ImE_2[i*n_t*n_eta + j*n_eta + k];
		  inp>>ReE_2[i*n_t*n_eta + j*n_eta + k];
		  inp>>ImHtild_2[i*n_t*n_eta + j*n_eta + k];
		  inp>>ReHtild_2[i*n_t*n_eta + j*n_eta + k];
		  inp>>Dterm_2[i*n_t*n_eta + j*n_eta + k];
		  
		  //cout<<i<<" "<<j<<" "<<k<<" Q2="<<arr_q2[i]<<" t="<<arr_t[j]<<" eta="<<arr_eta[k]<<endl;
		}
	    }
	}
      fstat = 1;
    }
  else
    {
      //cout<<"Th input File "<<fname<<" was not found"<<endl;
      //GPDs::~GPDs();
      fstat = 0;
    }
  
  inp.close();
  
  for( int i = 0; i < n_q2; i++ )
    {
      v_q2.push_back(arr_q2[i]);
    }
  
  for( int i = 0; i < n_t; i++ )
    {
      v_t.push_back(arr_t[i]);
    }
  
  for( int i = 0; i < n_eta; i++ )
    {
      v_eta.push_back(arr_eta[i]);
    }
  return fstat;
}


void GPDs::DefineValues()
{
  int ind_q2 = int(lower_bound(v_q2.begin(), v_q2.end(), Q2) - v_q2.begin()) - 1;
  int ind_tM = int(lower_bound(v_t.begin(), v_t.end(), tM) - v_t.begin()) - 1;
  int ind_eta = int(lower_bound(v_eta.begin(), v_eta.end(), eta) - v_eta.begin()) - 1;
  
  //cout<<"ind_q2="<<ind_q2<<"  ind_tM="<<ind_tM<<"  ind_eta="<<ind_eta<<endl;
  //cout<<"Q2="<<Q2<<"  tM="<<tM<<"  eta="<<eta<<endl;
  
  //assert( ind_q2 >= 0 && ind_tM >=0 && ind_eta >= 0 && ind_q2 < n_q2 - 1 && ind_tM < n_t - 1 && ind_eta < n_eta - 1);  
  assert( ind_q2 >= 0 );
  assert( ind_q2 < n_q2 - 1);
  assert( ind_tM >=0 );
  assert( ind_tM < n_t - 1 );
  assert( ind_eta >= 0 );
  assert( ind_eta < n_eta - 1);  
  
  q21 =  v_q2[ind_q2];
  q22 = v_q2[ind_q2 + 1];
  t1 = v_t[ind_tM];
  t2 = v_t[ind_tM + 1];
  eta1 = v_eta[ind_eta];
  eta2 = v_eta[ind_eta + 1];
  
  double f111, f112, f121, f122, f211, f212, f221, f222; // f(q2,t,eta)

  
  //cout<<"Q2_1="<<q21<<"\t q2_2-"<<q22<<endl;
  
  //==================ImH=========================
  f111 = ImH_2[ind_q2*n_t*n_eta + ind_tM*n_eta + ind_eta];
  f112 = ImH_2[ind_q2*n_t*n_eta + ind_tM*n_eta + ind_eta + 1];
  f121 = ImH_2[ind_q2*n_t*n_eta + (ind_tM + 1)*n_eta + ind_eta];
  f122 = ImH_2[ind_q2*n_t*n_eta + (ind_tM + 1)*n_eta + ind_eta + 1];
  f211 = ImH_2[(ind_q2 + 1)*n_t*n_eta + ind_tM*n_eta + ind_eta];
  f212 = ImH_2[(ind_q2 + 1)*n_t*n_eta + ind_tM*n_eta + ind_eta + 1];
  f221 = ImH_2[(ind_q2 + 1)*n_t*n_eta + (ind_tM + 1)*n_eta + ind_eta];
  f222 = ImH_2[(ind_q2 + 1)*n_t*n_eta + (ind_tM + 1)*n_eta + ind_eta + 1];
  
  ImH = (1/((q22-q21)*(t2-t1)*(eta2-eta1)))*
    (f111*(q22 - Q2)*(t2 - tM)*(eta2 - eta) + f112*(q22 - Q2)*(t2 - tM)*(eta - eta1)
    + f121*(q22 - Q2)*(tM - t1)*(eta2 - eta) + f122*(q22 - Q2)*(tM - t1)*(eta - eta1)
    + f211*(Q2 - q21)*(t2 - tM)*(eta2 - eta) + f212*(Q2 - q21)*(t2 - tM)*(eta - eta1)
    + f221*(Q2 - q21)*(tM - t1)*(eta2 - eta) + f222*(Q2 - q21)*(tM - t1)*(eta - eta1));
  

  //==================ReH=========================
  f111 = ReH_2[ind_q2*n_t*n_eta + ind_tM*n_eta + ind_eta];
  f112 = ReH_2[ind_q2*n_t*n_eta + ind_tM*n_eta + ind_eta + 1];
  f121 = ReH_2[ind_q2*n_t*n_eta + (ind_tM + 1)*n_eta + ind_eta];
  f122 = ReH_2[ind_q2*n_t*n_eta + (ind_tM + 1)*n_eta + ind_eta + 1];
  f211 = ReH_2[(ind_q2 + 1)*n_t*n_eta + ind_tM*n_eta + ind_eta];
  f212 = ReH_2[(ind_q2 + 1)*n_t*n_eta + ind_tM*n_eta + ind_eta + 1];
  f221 = ReH_2[(ind_q2 + 1)*n_t*n_eta + (ind_tM + 1)*n_eta + ind_eta];
  f222 = ReH_2[(ind_q2 + 1)*n_t*n_eta + (ind_tM + 1)*n_eta + ind_eta + 1];

  ReH = (1/((q22-q21)*(t2-t1)*(eta2-eta1)))*
    (f111*(q22 - Q2)*(t2 - tM)*(eta2 - eta) + f112*(q22 - Q2)*(t2 - tM)*(eta - eta1)
    + f121*(q22 - Q2)*(tM - t1)*(eta2 - eta) + f122*(q22 - Q2)*(tM - t1)*(eta - eta1)
    + f211*(Q2 - q21)*(t2 - tM)*(eta2 - eta) + f212*(Q2 - q21)*(t2 - tM)*(eta - eta1)
    + f221*(Q2 - q21)*(tM - t1)*(eta2 - eta) + f222*(Q2 - q21)*(tM - t1)*(eta - eta1));
  
  //==================ImE=========================
  f111 = ImE_2[ind_q2*n_t*n_eta + ind_tM*n_eta + ind_eta];
  f112 = ImE_2[ind_q2*n_t*n_eta + ind_tM*n_eta + ind_eta + 1];
  f121 = ImE_2[ind_q2*n_t*n_eta + (ind_tM + 1)*n_eta + ind_eta];
  f122 = ImE_2[ind_q2*n_t*n_eta + (ind_tM + 1)*n_eta + ind_eta + 1];
  f211 = ImE_2[(ind_q2 + 1)*n_t*n_eta + ind_tM*n_eta + ind_eta];
  f212 = ImE_2[(ind_q2 + 1)*n_t*n_eta + ind_tM*n_eta + ind_eta + 1];
  f221 = ImE_2[(ind_q2 + 1)*n_t*n_eta + (ind_tM + 1)*n_eta + ind_eta];
  f222 = ImE_2[(ind_q2 + 1)*n_t*n_eta + (ind_tM + 1)*n_eta + ind_eta + 1];

  ImE = (1/((q22-q21)*(t2-t1)*(eta2-eta1)))*
    (f111*(q22 - Q2)*(t2 - tM)*(eta2 - eta) + f112*(q22 - Q2)*(t2 - tM)*(eta - eta1)
    + f121*(q22 - Q2)*(tM - t1)*(eta2 - eta) + f122*(q22 - Q2)*(tM - t1)*(eta - eta1)
    + f211*(Q2 - q21)*(t2 - tM)*(eta2 - eta) + f212*(Q2 - q21)*(t2 - tM)*(eta - eta1)
    + f221*(Q2 - q21)*(tM - t1)*(eta2 - eta) + f222*(Q2 - q21)*(tM - t1)*(eta - eta1));

  
  //==================ReE=========================
  f111 = ReE_2[ind_q2*n_t*n_eta + ind_tM*n_eta + ind_eta];
  f112 = ReE_2[ind_q2*n_t*n_eta + ind_tM*n_eta + ind_eta + 1];
  f121 = ReE_2[ind_q2*n_t*n_eta + (ind_tM + 1)*n_eta + ind_eta];
  f122 = ReE_2[ind_q2*n_t*n_eta + (ind_tM + 1)*n_eta + ind_eta + 1];
  f211 = ReE_2[(ind_q2 + 1)*n_t*n_eta + ind_tM*n_eta + ind_eta];
  f212 = ReE_2[(ind_q2 + 1)*n_t*n_eta + ind_tM*n_eta + ind_eta + 1];
  f221 = ReE_2[(ind_q2 + 1)*n_t*n_eta + (ind_tM + 1)*n_eta + ind_eta];
  f222 = ReE_2[(ind_q2 + 1)*n_t*n_eta + (ind_tM + 1)*n_eta + ind_eta + 1];

  ReE = (1/((q22-q21)*(t2-t1)*(eta2-eta1)))*
    (f111*(q22 - Q2)*(t2 - tM)*(eta2 - eta) + f112*(q22 - Q2)*(t2 - tM)*(eta - eta1)
    + f121*(q22 - Q2)*(tM - t1)*(eta2 - eta) + f122*(q22 - Q2)*(tM - t1)*(eta - eta1)
    + f211*(Q2 - q21)*(t2 - tM)*(eta2 - eta) + f212*(Q2 - q21)*(t2 - tM)*(eta - eta1)
    + f221*(Q2 - q21)*(tM - t1)*(eta2 - eta) + f222*(Q2 - q21)*(tM - t1)*(eta - eta1));
  
  
  //==================ImHtild=========================
  f111 = ImHtild_2[ind_q2*n_t*n_eta + ind_tM*n_eta + ind_eta];
  f112 = ImHtild_2[ind_q2*n_t*n_eta + ind_tM*n_eta + ind_eta + 1];
  f121 = ImHtild_2[ind_q2*n_t*n_eta + (ind_tM + 1)*n_eta + ind_eta];
  f122 = ImHtild_2[ind_q2*n_t*n_eta + (ind_tM + 1)*n_eta + ind_eta + 1];
  f211 = ImHtild_2[(ind_q2 + 1)*n_t*n_eta + ind_tM*n_eta + ind_eta];
  f212 = ImHtild_2[(ind_q2 + 1)*n_t*n_eta + ind_tM*n_eta + ind_eta + 1];
  f221 = ImHtild_2[(ind_q2 + 1)*n_t*n_eta + (ind_tM + 1)*n_eta + ind_eta];
  f222 = ImHtild_2[(ind_q2 + 1)*n_t*n_eta + (ind_tM + 1)*n_eta + ind_eta + 1];

  ImHtild = (1/((q22-q21)*(t2-t1)*(eta2-eta1)))*
    (f111*(q22 - Q2)*(t2 - tM)*(eta2 - eta) + f112*(q22 - Q2)*(t2 - tM)*(eta - eta1)
    + f121*(q22 - Q2)*(tM - t1)*(eta2 - eta) + f122*(q22 - Q2)*(tM - t1)*(eta - eta1)
    + f211*(Q2 - q21)*(t2 - tM)*(eta2 - eta) + f212*(Q2 - q21)*(t2 - tM)*(eta - eta1)
    + f221*(Q2 - q21)*(tM - t1)*(eta2 - eta) + f222*(Q2 - q21)*(tM - t1)*(eta - eta1));
  
  //==================ReHtild=========================
  f111 = ReHtild_2[ind_q2*n_t*n_eta + ind_tM*n_eta + ind_eta];
  f112 = ReHtild_2[ind_q2*n_t*n_eta + ind_tM*n_eta + ind_eta + 1];
  f121 = ReHtild_2[ind_q2*n_t*n_eta + (ind_tM + 1)*n_eta + ind_eta];
  f122 = ReHtild_2[ind_q2*n_t*n_eta + (ind_tM + 1)*n_eta + ind_eta + 1];
  f211 = ReHtild_2[(ind_q2 + 1)*n_t*n_eta + ind_tM*n_eta + ind_eta];
  f212 = ReHtild_2[(ind_q2 + 1)*n_t*n_eta + ind_tM*n_eta + ind_eta + 1];
  f221 = ReHtild_2[(ind_q2 + 1)*n_t*n_eta + (ind_tM + 1)*n_eta + ind_eta];
  f222 = ReHtild_2[(ind_q2 + 1)*n_t*n_eta + (ind_tM + 1)*n_eta + ind_eta + 1];

  ReHtild = (1/((q22-q21)*(t2-t1)*(eta2-eta1)))*
    (f111*(q22 - Q2)*(t2 - tM)*(eta2 - eta) + f112*(q22 - Q2)*(t2 - tM)*(eta - eta1)
    + f121*(q22 - Q2)*(tM - t1)*(eta2 - eta) + f122*(q22 - Q2)*(tM - t1)*(eta - eta1)
    + f211*(Q2 - q21)*(t2 - tM)*(eta2 - eta) + f212*(Q2 - q21)*(t2 - tM)*(eta - eta1)
    + f221*(Q2 - q21)*(tM - t1)*(eta2 - eta) + f222*(Q2 - q21)*(tM - t1)*(eta - eta1));
  
  //==================Dterm=========================
  f111 = Dterm_2[ind_q2*n_t*n_eta + ind_tM*n_eta + ind_eta];
  f112 = Dterm_2[ind_q2*n_t*n_eta + ind_tM*n_eta + ind_eta + 1];
  f121 = Dterm_2[ind_q2*n_t*n_eta + (ind_tM + 1)*n_eta + ind_eta];
  f122 = Dterm_2[ind_q2*n_t*n_eta + (ind_tM + 1)*n_eta + ind_eta + 1];
  f211 = Dterm_2[(ind_q2 + 1)*n_t*n_eta + ind_tM*n_eta + ind_eta];
  f212 = Dterm_2[(ind_q2 + 1)*n_t*n_eta + ind_tM*n_eta + ind_eta + 1];
  f221 = Dterm_2[(ind_q2 + 1)*n_t*n_eta + (ind_tM + 1)*n_eta + ind_eta];
  f222 = Dterm_2[(ind_q2 + 1)*n_t*n_eta + (ind_tM + 1)*n_eta + ind_eta + 1];

  Dterm = (1/((q22-q21)*(t2-t1)*(eta2-eta1)))*
    (f111*(q22 - Q2)*(t2 - tM)*(eta2 - eta) + f112*(q22 - Q2)*(t2 - tM)*(eta - eta1)
    + f121*(q22 - Q2)*(tM - t1)*(eta2 - eta) + f122*(q22 - Q2)*(tM - t1)*(eta - eta1)
    + f211*(Q2 - q21)*(t2 - tM)*(eta2 - eta) + f212*(Q2 - q21)*(t2 - tM)*(eta - eta1)
    + f221*(Q2 - q21)*(tM - t1)*(eta2 - eta) + f222*(Q2 - q21)*(tM - t1)*(eta - eta1));

}

void GPDs::Set_q2_t_eta(double q2_, double t_, double eta_)
{
  tM = -t_;
  Q2 = q2_;
  eta = eta_;
  
  DefineValues();
}

double GPDs::GetImH() const
{
  return ImH;
}

double GPDs::GetReH() const
{
  return ReH;
}

double GPDs::GetImE() const
{
  return ImE;
}

double GPDs::GetReE() const
{
  return ReE;
}

double GPDs::GetImHtild() const
{
  return ImHtild;
}

double GPDs::GetReHtild() const
{
  return ReHtild;
}

double GPDs::GetDterm() const
{
  return Dterm;
}

GPDs::~GPDs()
{
  delete arr_q2;
  delete arr_t;
  delete arr_eta;
  
  delete ImH_2;
  delete ReH_2;
  delete ImE_2;
  delete ReE_2;
  delete ImHtild_2;
  delete ReHtild_2;
  delete Dterm_2;
}
