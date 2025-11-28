//#include <stdio.h>
//#include <nr3.h>
//#include "Vector.C"



#define P2 6.28318530717959 

//using namespace std;

class KinematicModel3D
{
private:
  size_t nk; // numero di modi
  VecDoub k, A, SS;
  MatDoub omega, phi;
  Doub epsilon, AspectRatio, eta;

  
public:
  KinematicModel3D(size_t n): k(n), A(n), SS(n), omega(n,3), phi(n,2)
  {
    nk = n;
  }

  KinematicModel3D(Doub L_start, Doub L_end, Doub fact, Doub eps, Doub decay)
  {
    if(!global.Kin3D)
      {
	nk = 0;
	return;
      }
   
    /* version isotropica */
    size_t i,j, n;
    Doub L;

    epsilon = eps;
    eta = decay;

    n=0; L = L_end;
    while(L > L_start) 
      {
	n++;
	L /= fact;
      }

    nk = n;
    k.resize(n);
    A.resize(n);
    SS.resize(n);
    omega.resize(n,3);
    phi.resize(n,2);

    L *= fact;

    for(i=0; i < nk; i++)
      {
	k[i] = P2 / L;
	A[i] = pow(epsilon / k[i], 0.3333333);
	SS[i] = 0.2 * L;

	omega[i][0] = 2. * k[i] * A[i];
	omega[i][1] = P2 / 6. * omega[i][0];
	omega[i][2] = P2 / 4. / sqrt(2.5) * omega[i][0];

	phi[i][0] = ranf_arr_next() * P2;
	phi[i][1] = ranf_arr_next() * P2;

	//  fprintf(stderr,"%lf %lf %lf\n",L, A[i], P2 * L / A[i]);
	
	L *= fact;
      }
    SS[nk-1] = 0; // per fissare l'ultima scala in modo che faccia da barriera
    AspectRatio = sqrt(2.0);     
  }

  Vec3<double> u(Vec3<double> X, double t, size_t i)
  {
    double arg1, arg2, arg3, expo, uk3;
    
    expo = exp(X.Z() / eta); 

    uk3 = 1. / k[i] / AspectRatio / eta;
    
    arg1 = k[i] * (X.X() - SS[i] * sin(omega[i][0] * t));
    arg2 = k[i] * (X.Y() - SS[i] * sin(omega[i][1] * t + phi[i][0]));
    // arg3 = AspectRatio * k[i] * (X.Z() - SS[i] * cos(omega[i][2] * t + phi[i][1]));

    arg3 = AspectRatio * k[i] * (X.Z());

    Vec3<double> V(\
		   A[i] * expo * ( sin(arg1) * cos(arg3) + uk3 * sin(arg1) * sin(arg3) ), \
		  -A[i] * expo * ( cos(arg3) * sin(arg2) + uk3 * sin(arg2) * sin(arg3) ), \
		   A[i] * expo / AspectRatio * ( -cos(arg1) * sin(arg3) + cos(arg2) * sin(arg3) ) \
		   );


    return V;
  }

  Vec3<double> u(Vec3<double> X, double t)
  {
    size_t i;
    Vec3<double> V;

    for(i=0; i<nk; i++) V += u(X, t, i);

    return V;
  }
  
  friend std::ostream& operator<< (std::ostream&, const KinematicModel3D&);
};

class KinematicModel2D
{
private:
  size_t nk; // numero di modi
  VecDoub k, A, SS, omega, phi;
  Doub enstrophy;

  
public:
  KinematicModel2D(size_t n): k(n), A(n), SS(n), omega(n), phi(n)
  {
    nk = n;
  }

  KinematicModel2D(Doub L_start, Doub L_end, Doub fact, Doub ens)
  {
    if(!global.Kin2D)
      {
	nk = 0;
	return;
      }


    size_t i,j, n;
    Doub L;

    enstrophy = ens;

    n=0; L = L_end;
    while(L > L_start) 
      {
	n++;
	L /= fact;
      }

    nk = n;
    k.resize(n);
    A.resize(n);
    SS.resize(n);
    omega.resize(n);
    phi.resize(n);

    L *= fact;

    for(i=0; i < nk; i++)
      {
	k[i] = P2 / L;
	//A[i] = enstrophy; // da completare!
	A[i] = pow(enstrophy*L/2., 0.333333); // da completare!
	SS[i] = 0.2 * L;

	omega[i] =  k[i] * A[i];
	phi[i] = ranf_arr_next() * P2;

	L *= fact;
      }

  }

  Vec3<double> u(Vec3<double> X, double t, size_t i)
  {
    double arg1, arg2;
    

    arg1 = k[i] * (X.X() - SS[i] * sin(omega[i] * t));
    arg2 = k[i] * (X.Y() - SS[i] * sin(omega[i] * t + phi[i]));

    Vec3<double> V(\
		    A[i] * sin(arg1) * cos(arg2), \
		    -A[i] * cos(arg1) * sin(arg2),\
		    0.\
		   );
    return V;
  }

  double psi(Vec3<double> X, double t, size_t i)
  {
    double arg1, arg2;
    

    arg1 = k[i] * (X.X() - SS[i] * sin(omega[i] * t));
    arg2 = k[i] * (X.Y() - SS[i] * sin(omega[i] * t + phi[i]));

    return sin(arg1) * sin(arg2) / k[i];
  }


  Vec3<double> u(Vec3<double> X, double t)
  {
    size_t i;
    Vec3<double> V(0.,0.,0.);

    for(i=0; i<nk; i++) V += u(X, t, i);

    return V;
  }

  double psi(Vec3<double> X, double t)
  {
    size_t i;
    double p = 0.;

    for(i=0; i<nk; i++) p += psi(X, t, i);

    return p;
  }


  friend std::ostream& operator<< (std::ostream&, const KinematicModel2D&);
  
};


inline std::ostream& operator << (std::ostream& ostream, const KinematicModel3D& C)
{
  size_t i; // numero di modi

  if(!global.Kin3D) return ostream;  

  ostream << "\nNMODI " << C.nk << endl;
  
  ostream << "k | L | \t A_k | \t omega_k[0]  omega_k[1]   omega_k[2] | \t phi_k[0] phi_k[1] \n" << endl;
  ostream << "---------------------------\n" << endl;

  for(i=0; i<C.nk; i++)
    {
      ostream << C.k[i] << "\t" << P2 / C.k[i] << "\t" <<C.A[i] << "\t" << " ";
      ostream << C.omega[i][0] << " " << C.omega[i][1] << " " << C.omega[i][2] << " ";
      ostream << C.phi[i][0] << " " << C.phi[i][1] << endl;
    }

  ostream << endl;

  ostream << "Aspect Ratio " << (double)C.AspectRatio << endl;
  ostream << "Eta " << (double)C.eta << endl;
  ostream << "Epsilon " << (double)C.epsilon << endl << endl;
 

  return ostream;
}

inline std::ostream& operator << (std::ostream& ostream, const KinematicModel2D& C)
{
  size_t i; // numero di modi

  if(!global.Kin2D) return ostream;  

  ostream << "\nNMODI " << C.nk << endl;
  
  ostream << "k | L | \t A_k | \t omega_k | \t phi_k \n" << endl;
  ostream << "---------------------------\n" << endl;

  for(i=0; i<C.nk; i++)
    {
      ostream << C.k[i] << "\t" << P2 / C.k[i]  << "\t" << C.A[i] << "\t" << " ";
      ostream << C.omega[i] << " ";
      ostream << C.phi[i] << endl;
    }

  ostream << endl;

  ostream << "Enstrophy " << (double)C.enstrophy << endl << endl;
 

  return ostream;
}
