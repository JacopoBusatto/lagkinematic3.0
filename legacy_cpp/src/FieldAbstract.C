/**********************************************************************************************/
/*  Questa e' la classe astratta, ogni implementazione deve definire                         */
/*  void Init(....)                                                                          */
/*  void Load(....)                                                                          */	
/*  double Value(X,Y,Z)                                                                      */
/*  double Value(X,Y,Z, Neighbour)                                                           */
/*                                                                                           */
/*  La classe derivata Mask ha bisogno                                                       */
/*  bool Inside(Vec3<double>& p)                                                             */
/*  double orografia(Vec3<double>& p)                                                        */
/*  double toplevel(Vec3<double>& p)                                                         */
/**********************************************************************************************/

//#include <iostream>
//#include <stdio.h>
//#include <string.h>


#include "netcdf.h"

#include <nr3.h>
#include <ludcmp.h>
#include <interp_1d.h>
#include <interp_linear.h>
#include <interp_2d.h>

enum KindOfFile_t {MFS_MyOcean = 0, Mercator = 1, HyCom = 2, AVISO = 3, WRF = 4, \
		   SST_L4 = 5, SST_L3S = 6, GERRIS_SW = 7, OCEAN_COLOR = 8, OceanColor_MyOcean2D = 9};

#include "Vector.C"
#include "Global.C"
#include "CerchioMassimo.C"
#include "Neighbours.C"
#include "interp_cerchio_massimo.h"


using namespace std;


class ScalarField3D {
protected:

  double time; // in secondi da una data di riferimento?
  bool extension; // dice se c'e' l'estensione al bottom del campo di velocita'm tipo u10, v10


public:
  size_t n_punti, nz;  size_t ntime; /* l'implementazione di Init si occupa di assegnarlo*/
  VecDoub z; /* z serve sia da z che da livelli sigma */  
 
 ScalarField3D(){}
       
  virtual void Init(const char netcdf_file[], const char nome_lon[], \
		    const char nome_lat[], const char nome_depth[], const char nomeExtension[]) = 0;

  inline void Init(const char netcdf_file[], const char nome_lon[],	\
		    const char nome_lat[], const char nome_depth[])
  {
    cout << "No extension\n";
    Init(netcdf_file, nome_lon, nome_lat, nome_depth,"NULL");
  }
  

  virtual void Load(const char netcdf_file[], const char nome_var[], size_t T, const char nomeExtension[]) = 0;

  inline void Load(const char netcdf_file[], const char nome_var[], size_t T)
  {
    Load(netcdf_file, nome_var, T, "NULL");
  }

  virtual void CloseDivergence(ScalarField3D& U, ScalarField3D& V)=0;

  virtual double Value(double, double, double) = 0;

  virtual double Value(double, double, double, Neighbours&)=0;

  inline double Value(Vec3<double> X) { return Value(X.X(), X.Y(), X.Z()); }
  
  inline double Value(Vec3<double> X, Neighbours &vicini) 
  { return Value( X.X(), X.Y(), X.Z(), vicini);}


  virtual size_t NTimes(const char netcdf_file[])=0;

  virtual MatDoub& Points4Neighbour(bool& isUsable) = 0;

  size_t Nz(){ return nz;}

  double Time(){return time;}
  
  void Time(double t){time = t;}

  

  Vec3<double> Gradient(Vec3<double> X)
  {
    Vec3<double> g, dx(1,0,0), dy(0,1,0);
    Doub gx, gy, h = 0.02, N; /*h circa 1 km */
    
    try
      {
	gx = 0.5 / h * (Value(X + dx * h) - Value(X + dx * (-h)));
	gy = 0.5 / h * (Value(X + dy * h) - Value(X + dy * (-h)));
      }    
    catch( int ercode )
      {
	cout << "Calcolo del gradiente fuori dominio" << endl;
	global.er << "Calcolo del gradiente fuori dominio" << endl;
	global.er << X.X() << " " << X.Y() << " " << X.Z() << endl;
	throw ercode;
      }
    
    g.setX(gx); g.setY(gy);
    g.setZ(0.);

    return g;
  }

  size_t Ntot(){ return n_punti; }

  Vec3<double> Gradient(Vec3<double> X, Neighbours& neigh)
  {
    Vec3<double> g, dx(1,0,0), dy(0,1,0);
    Doub gx, gy, h = 0.1, N;
    
    try
      {
	gx = 0.5 / h * (Value(X + dx * h, neigh) - Value(X + dx * (-h), neigh ));
	gy = 0.5 / h * (Value(X + dy * h, neigh) - Value(X + dy * (-h), neigh ));
      }    
    catch( int ercode )
      {
	global.er << "Calcolo del gradiente fuori dominio" << endl;
	global.er << X.X() << " " << X.Y() << " " << X.Z() << endl;
	throw ercode;
      }
    
    g.setX(gx); g.setY(gy);
    g.setZ(0.);

    return g;
  }

};

/*%%%%%%%%%%%%%%%%%%%%%% Classe Mask %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

class ScalarField3D_Mask : public ScalarField3D{


public:

  ScalarField3D_Mask() : ScalarField3D(){}

  virtual void Mask(const char netcdf_file[], const char nome_var[])=0;
  /* carica la maschera invede della variabile   */

  virtual bool Inside(Vec3<double>& p)=0; /* Inside e' solo laterale al confine del dominio*/
  
  virtual double bottom(Vec3<double>& p)=0; // batimetria nel caso Ocean

  virtual double bottom(Vec3<double>& p, Neighbours& vicini)=0; // batimetria nel caso Ocean
  
  virtual double toplevel(Vec3<double>& p)=0; 

  virtual double toplevel(Vec3<double>& p, Neighbours& vicini)=0; 

  virtual void LoadOrografia(const char netcdf_file[], const char nome_var[])=0;

  void SetExtension()
  {
    extension=true;
  }

  bool InDomain(Vec3<double>& p)
  {
    double top = toplevel(p), e=1.e-5;
    //  double bot = bottom(p);

    //    if(p.X() <= -180. || p.X() >= 180.) return true;

    if( Value(p.X(), p.Y(), top - e) < global.tt ) return false;        

    //    cout << "controllino " << top << " " << bot << endl;

    // if( top <= bot) return false;

    //    if(Value(p) < global.tt) return false;

    return true;
  }
  

  bool AbsorptionTOP(Vec3<double>& p) 
  {
    if(global.buoyant) return false;

    double top = toplevel(p);
    if(p.Z() >= top) return true;
    else
      return false;
  }

  bool AbsorptionBOTTOM(Vec3<double>& p)
  {
    /* sono lontano dal fondo? */ 
    // if( Value(p) > global.tt ) return false;    

    if(global.buoyant) return false;

    double bot = bottom(p);

    if(p.Z() <= bot) 
   	return true; // avvisa se si e' toccato il fondo
    else
      return false;
  }

  bool ReboundTOP(Vec3<double>& p)
  {
    if(global.buoyant) return false;

    const double e = 1.e-5;
    double old;

    old = p.Z();

    double top = toplevel(p);
    if(p.Z() >= top)
      {
	p.setZ(2. * top - p.Z() - e);
	//	global.er << "top " << p.Z() << " " << old << " " << top << endl;
	//	global.er.flush();
	return true;
      }

    cout.flush();
    return false;
  }

  bool ReboundBOTTOM(Vec3<double>& p)
  {
    /* sono lontano dal fondo? */ 
    
    if(global.buoyant) return false;

    try{
      if( Value(p) > global.tt ) return false;    
    }
    catch (int ercode)
      {
	global.er << "Rimbalzo Bottom problematico " << p.Z() << endl;
	return false;
      }

    double bot = bottom(p);
    const double e = 1.e-5;
    double old;

    old = p.Z();

    if(p.Z() <= bot) 
      {
	p.setZ(2. * bot - p.Z() + e);
	//	global.er << "bottom " << p.Z() << " " << old << " " << bot << endl;
	//	global.er.flush();

	return true; // avvisa se si e' toccato il fondo
      }
    /*             in tutti gli altri casi ritorna false */


    cout.flush();

    return false;
  }

  /* Lateral si intende sul terreno o sulla costa, alla fine del dominio 
   vengono sempre killed */


  bool ReboundLateral(Vec3<double>&p)
  {
    double V;
    const double security = 0.05; /* di quanto voglio stare sopra global.tt  */

    try{
    V = Value(p);
    }
    catch (int ercode)
      {
	cout << "l'errore sgrunt e' qui \n\n"; 
	return false;
	//throw int (-99);
      }

    if( V > global.kt3d ) return false;    

    Vec3<double> H, dx;
    double hh, fac;

    /*    try{
      H = Gradient(p);
    }
    catch(int ercode)
      {
	cout << "l'errore porco giuda e' qui nel gradiente\n";
	}*/


    /*  p e' ortogonale alla costa/boundary e diretto verso il mare/aria*/

    long count = 0;

    try{
      while( (V = Value(p)) <= global.kt3d  && Inside(p) && count < 20) 
	{
	  H = Gradient(p);
	  hh = H.magnitude();

	  if (hh < 1.e-6) return false;
	  H.normalise();  
	  
	  //	  fac = -0.3 * (V - global.kt3d ) / hh;
	  fac = - (V - 1.) / hh;

	  //	  global.absorptionFile << "iter " << p.X() << " " << p.Y() << " " << Value(p) << " " << H.X() << " " << H.Y() << endl;
	  p += fac * H;
	  count ++;
	}

      //      global.absorptionFile << "iterend " << p.X() << " " << p.Y() << " " << Value(p) << " " << H.X() << " " << H.Y() << endl;
    }
    catch(int ercode)
      {
	global.er << "era in questo loop dentro rebound lateral" << endl;
	global.er.flush();
      }

    return true;
  }

  //  size_t NTimes(){ return ScalarField3D::ntime;}

  bool AbsorptionLateral(Vec3<double>&p)
  {
    //    cout << Value(p) << " " << global.tt << endl;
    
    if( Value(p) > global.tt ) return false;    
    else return true;
  }

};


/* include i dettagli delle implementazioni */

#include "MFS_MyOcean.C" 
#include "MFS_MyOcean2D.C" //\RC

#include "Mercator.C" 

#include "WRF.C" 

/* coming soon 

#include "HyCom.C" 

#include "Altimetry.C" 

*/

#include "GerrisLuigi.C" 

class ScalarField3D* Switcher(KindOfFile_t type)
{
  //  cout << "Il tipo di file e' " << type << endl;
  // cout << "corrispondente a ";


  switch( type ){

  case MFS_MyOcean:
    cout << " MFS_MyOcean" << endl;
    global.sigma = false;
    Lista.ocean = true;
    return new class ScalarField3D_MFS;
    break;

    //\RC
  case OceanColor_MyOcean2D:
    cout << " OceanColor_MyOcean 2D" << endl;
    global.sigma = false;
    Lista.ocean = true;
    return new class ScalarField2D_OcColor;
    break;

  case Mercator:
    global.sigma = false;
    global.ocean = true;
    return new ScalarField3D_Mercator;
    break;

  case WRF:
    //  cout << " WRF" << endl;
    global.sigma = true;
    Lista.ocean = false;
    return new class ScalarField3D_WRF;
    break;
  
  case GERRIS_SW:
    global.sigma = false;
    Lista.ocean = true;
    return new class ScalarField3D_Gerris;
    break;
  }
}


class ScalarField3D_Mask* SwitcherMask(KindOfFile_t type)
{
  switch( type ){

  case MFS_MyOcean:
    return new class ScalarField3D_Mask_MFS;
    break;

    //\RC
  case OceanColor_MyOcean2D:
    return new class ScalarField2D_Mask_OcColor;
    break;

  case Mercator:
    global.sigma = false;
    global.ocean = true;
    return new class ScalarField3D_Mask_Mercator;
    break;

  case WRF:
    global.sigma = true;
    Lista.ocean = false;
    return new class ScalarField3D_Mask_WRF;
    break;

  case GERRIS_SW:
    global.sigma = false;
    Lista.ocean = true;
    return new class ScalarField3D_Mask_Gerris;
    break;


  }
}

