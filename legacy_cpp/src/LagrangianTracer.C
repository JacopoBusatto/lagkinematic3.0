//#include <iostream>
//#include <math.h>

#define Lo0 12.
#define La0 36.


//#include "CerchioMassimo.C"
//#include "Neighbours.C"
//#include "interp_cerchio_massimo.h" //fmg

//#define ScalarField3D ScalarGerris
//#define VectorField3D VectorGerris

//#define R_earth 6378137 // equatoriale

class LagrangianTracer 
{
private:
  Vec3<double> X, X0; 
  // x, y, z in metri dal punto di origine (Long0, Lat0, 0)
  double LONG, LAT, temp, S;
  double age0, age, start; // start e' il tempo a cui iniziare ad evolvere (plume)
  bool flag2D, kinematic2D, kinematic3D, alive;
  short kind_of_tracer;  
  double orografia, top;
  static const int n_vicini=16;

public:
  // ------------ Constructors ------------

  //Per tartarughe e acciughe  
  double Psi, Chl;
  
  Neighbours vicinix, viciniy, viciniz;
  
  LagrangianTracer():vicinix(n_vicini), viciniy(n_vicini), viciniz(n_vicini)
  {
    X0.set(0., 0., 0.);
    X = X0;
    LONG = Lo0; LAT = La0;

    flag2D = 0;
    age = 0.; start = 0.;
    alive = 1;
    Psi = 1.; // probabilità di sopravvivenza
    Chl = 0.; // clorofilla incontrata durante il cammino

    // todo da controllare, luigi dice che nn fungerebbe vicini = Neighbours(1, n_vicini); //fmg
  }

  LagrangianTracer(const char KindOfTracer[]):vicinix(n_vicini), viciniy(n_vicini), viciniz(n_vicini)
  {
    X0.set(0., 0., 0.);
    X = X0;
    LONG = Lo0; LAT = La0;

    if(!strcmp(KindOfTracer,"passive")) kind_of_tracer = flag2D = 0;
    if(!strcmp(KindOfTracer,"buoyant")) kind_of_tracer = flag2D = 1;
    if(!strcmp(KindOfTracer,"anchovy")) 
      {
	kind_of_tracer = 2;
	flag2D = 1; // all'inizio
      }

    age = 0.; start = 0.;
    alive = 1;
  }

  void Set(const char KindOfTracer[])
  {
    if(!strcmp(KindOfTracer,"passive")) kind_of_tracer = flag2D = 0;
    if(!strcmp(KindOfTracer,"buoyant")) kind_of_tracer = flag2D = 1;
    if(!strcmp(KindOfTracer,"anchovy")) 
      {
	kind_of_tracer = 2;
	flag2D = 1; // all'inizio
      }
  }

  inline void SetBuoyant(){ flag2D = 1;}
  
  ~LagrangianTracer(){}

  void SetAlive() { alive = 1;}
  void SetDead() { alive = 0;}
  bool Alive() { return alive;}
  inline bool IsAnchovy(){ 
    if(kind_of_tracer == 2) return true;
    else return false;
  }


  //---- metodi gestione coordinate -------------

  Vec3<double> Xv(){ return X;}

  void AggiornaLongLat(Vec3<double> dx)
  {
    double dLat, dLon;
    
    //    dLat = ( dx.Y() / R_earth ) / P2 * 360.;
    dLat = dx.Y() / 111039;
    LAT += dLat;

    dLon = dx.X() / 111039 / cos(LAT * P2 / 360.);
    LONG += dLon;

    /* LP 30 nov 2019 */
    // if(LONG > 179.75) LONG -= 359.75;                    
    // if(LONG <= -180) LONG += 359.75;
    if(LONG >   179.875) LONG -= 359.75;                    
    if(LONG <= -179.875) LONG += 359.75;
    //    if(LONG > 180.) LONG -= 360.;
    // if(LONG < -180.) LONG += 360.;

    /* .........  */
  }

  Vec3<double> LongLatZ()
  {
    Vec3<double> out;
    out.set(LONG, LAT, X.Z());
    return out;
  }

  Vec3<double> LongLatSigma(class ScalarField3D* gpm, class ScalarField3D* gpa)
  {
    if(!Lista.sigma) return LongLatZ();

    Vec3<double> out;
    out.set(LONG, LAT, X.Z());

    out.setZ(Zeta(gpm, gpa));

    return out;
  }



  void SetLongLatZ(Vec3<double>& xl)
  {
    LONG = xl.X();
    LAT = xl.Y();
    
    const size_t nsteps = 100;
    double dlong, dlat, dx, dy, xt = 0., yt = 0.;

    yt = (LAT - La0) * 111039;
    xt = (LONG - Lo0) * 111039 *  cos( ( LAT + La0) * P2 / 360.);

    X.setX(xt);
    X.setY(yt);
    X.setZ(xl.Z());
    X0 = X;
    
  }
  
  double Long(){ return LONG; }

  double Lat(){ return LAT; }

  double T(){ return temp;}
  void T(double x){ temp = x;}

  //  void Set(const Vec3<double> X1){ X = X1;}

  void IncrSet(Vec3<double> dx)
  { 
    X += dx;
    AggiornaLongLat(dx);
  }

  double Age(){return age/86400.;}

  void Age(double& xx){ this->age = xx; this->age0 = xx;} // in secondi

  bool IsStarted()
  {
    if(global.Time() < start) return false;
    else return true;
  }

  void Start(double st){ start = st;}
  double Start(){ return start;}

  double Zeta(class ScalarField3D* gpm, class ScalarField3D* gpa)
  {
    /* data la z restituisce il valore di zeta o sigma compreso fra 0 e 1   */

    if(!Lista.sigma) return X.Z();

    /* gmp e' il campo medio, gpa la fluttuazione */
    double lon, lat;
    const double g = 9.81;

    lon = Long(); lat = Lat();

    VecDoub geopot(gpm->Nz()-1);

    // Nz()-1 perche' gli ultimi due valori sono uguali 

    for(size_t i=0; i<gpm->Nz()-1; i++)
      {
	try{
	  geopot[i] = ( gpm->Value(lon, lat, gpm->z[i], viciniz) + \
			gpa->Value(lon, lat, gpm->z[i], viciniz) ) / g;
	}
	catch (int ercode){
	  global.er << "Codice errore in Value: " << ercode << endl;
	  global.er << "Lanciato da LagragianTracer::Zeta\n";
	}	
      }

    if(X.Z() < geopot[0]) return gpm->z[0];
    if(X.Z() > geopot[gpm->Nz() - 2]) return gpm->z[gpm->Nz() - 2];

    try{
      Poly_interp GeopotVsZ(geopot, gpm->z, 5);

      return GeopotVsZ.interp(X.Z());
    }
    catch(int ercode){
      global.er << "Errore nell'interpolazione Zeta lanciato da LagrangianTracer::Zeta\n";
      for(size_t j=0; j<gpm->Nz() -1; j++) global.er << geopot[j] << " " << gpm->z[j] << endl;
      
      throw int(7777);
    }
  }


  // Evoluzione con campo risolto
  //! fa muovere i traccianti facendo muovere da un campo risolto ( tipo Gerris )
  // fieldNee e fielOld -> fa l'iterpolazione lineare fra nuovo e vecchio


void EvolveCampoRisolto(class VectorField3D& field_old, class VectorField3D& field_new, \
			class ScalarField3D_Mask* mask, class ScalarField3D* gpm, class ScalarField3D* gpa)
  {
    Vec3<double> U1, U2, U, xlonglat, M, dX;
    double f1, f2, dotP, dm, modM;

    if(!IsStarted()) return; // non e' ancora giunto il suo momento

    if(!alive) return; // se non e' viva non viene evoluta

    /* if( Age() < 0. ) 
      {
	xlonglat = LongLatZ();
	SetDead();
	global.absorptionFile << "Birth " << xlonglat.X() << " " << xlonglat.Y() << " "\
			      << xlonglat.Z() << " " << Age() << " " << global.Time() / 86400. \
			      << " " << age0/86400. << endl;
	return;
	}   LP 17 dic 2019, questo serviva per le acciughe */ 


    /**********************  BOUNDARY CONDITION **************************/

    ApplyBoundaryCondition(mask, gpm, gpa);

    xlonglat = LongLatZ();

    
    try{

      if(Lista.sigma)
	xlonglat.setZ( Zeta(gpm, gpa) );

      U1 = field_old.Value(xlonglat, vicinix, viciniy, viciniz);
      U2 = field_new.Value(xlonglat, vicinix, viciniy, viciniz);
    }
    catch(int ercode)
      {
	global.er << "Lanciata da EvolveCampoRisolto " << ercode << endl;
	global.er << "sul vettore "; 
	xlonglat.display(global.er);
	global.er << endl << endl;
	U1.zero(); U2.zero();
      }

    f2 = (global.Time() - field_old.Time() ) / Lista.DT;
    f1 = 1. - f2;

    U = (f1 * U1 + f2 * U2 ) ; 


    //    if(global.divergence)
    // if(true)
    if (false) /* spengo la chiusura della divergenza  */
      {
	double dwdz, wz1, wz2, hx, hy, hz = 50., ddx = 1./24;
	double XX, YY; 

	wz1 = wz2 = 0.;
	
	XX = xlonglat.X();
	YY = xlonglat.Y();

	hy = 110000. * ddx; // 1/24 di grado
	hx = 110000. * ddx * cos(YY*6.28318530717959 / 180.); // 1/24 di grado


	for(double zz = mask->bottom(xlonglat); zz < xlonglat.Z(); zz += hz )
	    {
	      wz1 -= 0.5 * hz * ( field_old.X->Value(XX+ddx,YY,zz) - field_old.X->Value(XX-ddx,YY,zz) ) / hx;
	      wz1 -= 0.5 * hz * ( field_old.Y->Value(XX,YY+ddx,zz) - field_old.Y->Value(XX,YY-ddx,zz) ) / hy;
	      
	      wz2 -= 0.5 * hz * ( field_new.X->Value(XX+ddx,YY,zz) - field_new.X->Value(XX-ddx,YY,zz) ) / hx;
	      wz2 -= 0.5 * hz * ( field_new.Y->Value(XX,YY+ddx,zz) - field_new.Y->Value(XX,YY-ddx,zz) ) / hy;

	      U.setZ(f1 * wz1 + f2 * wz2);

	      //	      cout << zz << " " << wz1 << " " << endl;
	      //              cout.flush();
	    }
      }


    

    if(global.buoyant) U.setZ(0.); // spengo la velocita' verticale risolta
    
    U.setZ(U.Z() - global.Vdep); // metto la velocita' di deposizione

    if(global.backtraj) 
      {
	U = -1. * U; 
	age -= global.dt();
      }
    else age += global.dt();

    IncrSet(U * global.dt()); // questo aggiorna anche LONG e LAT
    if(global.anchovy){
      // da palatella et al. 2014 timeeggs=1.5 days, time_up = 6am, time_down = 6pm, zup = -3m, zdown = -100m , gam=0.1 h^-1 (12/2/25 gam=1 h^-1)
      DialVerticalMigration(1.5, 6, 18, -3, -50, 1./3600, mask);
    }
  }


  inline void DialVerticalMigration(double timeeggs, double time_up, double time_down, double zup, double zdown, double gam, 	class ScalarField3D_Mask* mask)
  {
    long int hour, secs;
    double z0, z;
    Vec3<double> Uz(0., 0., 0.); 

    // 20250220 commentato JB
    // if (global.Time() <= 2 * global.dt()) {
    //   global.log <<"STO FACENDO DAVVERO LE ACCIUGHE\n";
    // }
    if(global.Time() < start) return; // non e' ancora giunto il suo momento
    

    // JACOPO BUSATTO 20-12-2024 flag2D > flag3D
    if(age < timeeggs){ flag2D = 1; return;}
    
    flag2D = 0; // non piu' galleggiante!
    
    secs = ((long int)global.Time()) % 86400;
    hour = secs / 3600;
    
    if(hour > time_up && hour < time_down) z0 = zup;
    else z0 = zdown;

    z = X.Z();

    Uz.setZ( -gam * (X.Z() - z0)  );
    
    IncrSet(Uz * global.dt());       

    //    ApplyBoundaryCondition(mask);
  }

  inline void DialVerticalMigration(class Anchovy& A, 	class ScalarField3D_Mask* mask)
  {
    DialVerticalMigration(A.timeeggs, A.time_up, A.time_down, A.zup, A.zdown, A.gam, mask);
  }


  void ApplyBoundaryCondition(ScalarField3D_Mask *mask, class ScalarField3D* gpm, class ScalarField3D* gpa)
  {
    /* questa routine applica tutte le boundary condition al Lagrangian tracer*/
    
    bool esito;
    Vec3<double> tmp, p, pt;

    p = LongLatZ(); pt = LongLatSigma(gpm, gpa);

    if(!mask->Inside(pt))
      {
	SetDead();
	global.absorptionFile << " e' stato inside  "<< p.X() << " " << p.Y() << " " << p.Z() << " " <<Age() << endl;
	return;
      }

    if(!Alive()) return;
    
    if(!mask->InDomain(pt))
      {
	SetDead();
	global.absorptionFile << " e' stato in domain" << p.X() << " " << p.Y() << " " << p.Z() << " " << Age() << endl;
	return;
      }



    /* TOP E BOTTOM SONO DATI IN METRI QUANDO I LIVELLI SONO SIGMA*/

    if(global.absorptionTOP) 
      {  
	esito = mask->AbsorptionTOP(p);
	
	if(esito) 
	  {
	    SetDead();
	    global.absorptionFile << p.X() << " " << p.Y() << " " << p.Z() << " " <<Age() << endl;
	    return; 
	  }
      }
    else
      if(mask->ReboundTOP(p)) X.setZ(p.Z());    


    p = LongLatZ(); pt = LongLatSigma(gpm, gpa);

    /***************************************/

    if(global.absorptionBOTTOM) 
      {  
	esito = mask->AbsorptionBOTTOM(p);
	
	if(esito) 
	  {
	    SetDead();
	    global.absorptionFile << p.X() << " " << p.Y() << " " << p.Z() << " " <<Age() << endl;
	    return;
	  }
      }
    else 
      if(mask->ReboundBOTTOM(p)) X.setZ(p.Z());

    p = LongLatZ(); pt = LongLatSigma(gpm, gpa);

    /***************************************/

    if(global.absorptionLATERAL) 
      {  
	esito = mask->AbsorptionLateral(pt);
	
	if(esito) 
	  {
	    SetDead();
	    global.absorptionFile << p.X() << " " << p.Y() << " " << p.Z() << " " <<Age() << endl;
	    return;
	  }
      }
    else
      if(mask->ReboundLateral(pt))
	{ 
	  if(Lista.sigma)
	    {
	      tmp = pt;
	      tmp.setZ(p.Z());
	      SetLongLatZ(tmp); /* controllare !!!*/
	    }
	  else SetLongLatZ(pt);
	}


  }

  friend class CoupleLagrangianTracer;
};
