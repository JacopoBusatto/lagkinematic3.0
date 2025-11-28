
//#include "Field3D.C"
//#include "KinematicModel.C"
//#include "LagrangianTracer.C"


class CoupleLagrangianTracer 
{
private:

public:
  // ------------ Constructors ------------
  double mask1, mask2; // servono solo per l'output

  LagrangianTracer part1, part2;  

  CoupleLagrangianTracer(): part1(), part2() {}

  CoupleLagrangianTracer(const char KindOfTracer[]): part1(KindOfTracer), part2(KindOfTracer) {}
  
  void Set(const char KindOfTracer[])
  {
    part1.Set(KindOfTracer);
    part2.Set(KindOfTracer);
  }

  inline void SetBuoyant()
  {
    part1.SetBuoyant();
    part2.SetBuoyant();
  }

  Vec3<double> CM()
  {
    Vec3<double> x;

    x = part1.Xv() + part2.Xv();
    return 0.5 * x;
  }

  Vec3<double> ModifiedCM()
  {
    Vec3<double> cm, xcm0;
    double zz;

    xcm0 = part1.X0 + part2.X0;
    //    zz = xcm0.Z(); // lungo z uso coordinate assolute

    cm =  ( part1.Xv() + part2.Xv() - xcm0 ); 
    // cm.setZ( 0.5 * ( part1.orografia + part2.orografia ) );
    cm.setZ(0.);

    return cm / 2.;
  }
  

  void EvolveAnalytical(class KinematicModel3D& subgrid, class ScalarField3D_Mask* mask,  \
			class ScalarField3D* gpm, class ScalarField3D* gpa)	
  {
    // subgrid esprime il campo in coordinate relatite x-y-z!
    Vec3<double> k1, k2, k3, k4, cm, x1, x2, x1ll, x2ll, dummy, xlonglat;
    double t, h, m1, m2, oro, top;

    //    cout << "cinematico 3D" << endl;

    if(!part1.Alive() || !part2.Alive()) return;

    if(!part1.IsStarted() || !part2.IsStarted()) return;


    part1.ApplyBoundaryCondition(mask, gpm, gpa);
    part2.ApplyBoundaryCondition(mask, gpm, gpa);


    t = global.Time(); // tempo corrente
    h = global.dt();

    x1 = part1.Xv();
    x2 = part2.Xv();

    x1ll = part1.LongLatZ();
    x2ll = part2.LongLatZ();

    
    cm = ModifiedCM();

    // qui inizio Runge-Kutta al quarto ordine

    k1 = subgrid.u(x1 - cm, t ) * h; if(part1.flag2D) k1.setZ(0.);
    k2 = subgrid.u(x1 - cm + 0.5 * k1, t + h/2.) * h; if(part1.flag2D) k2.setZ(0.);
    k3 = subgrid.u(x1 - cm + 0.5 * k2, t + h/2.) * h; if(part1.flag2D) k3.setZ(0.);
    k4 = subgrid.u(x1 - cm + k3,t + h) * h; if(part1.flag2D) k4.setZ(0.);
    
    part1.IncrSet(k1/6. + k2/3. + k3/3. + k4/6.);

    k1 = subgrid.u(x2 - cm, t ) * h; if(part2.flag2D) k1.setZ(0.);
    k2 = subgrid.u(x2 - cm + 0.5 * k1, t + h/2.) * h; if(part2.flag2D) k2.setZ(0.);
    k3 = subgrid.u(x2 - cm + 0.5 * k2, t + h/2.) * h; if(part2.flag2D) k3.setZ(0.);
    k4 = subgrid.u(x2 - cm + k3,t + h) * h; if(part2.flag2D) k4.setZ(0.);
    
    part2.IncrSet(k1/6. + k2/3. + k3/3. + k4/6.);

   part1.ApplyBoundaryCondition(mask, gpm, gpa);
   part2.ApplyBoundaryCondition(mask, gpm, gpa);


  }


  //  void EvolveAnalytical(class KinematicModel2D& subgrid,		\
  //			class KinematicModel2D& subgrid2, class ScalarField3D& mask)
  //

  void EvolveAnalytical(class KinematicModel2D& subgrid, class ScalarField3D_Mask* mask,\
			class ScalarField3D* gpm, class ScalarField3D* gpa)	
  {
    // subgrid esprime il campo in coordinate relatite x-y-z!
    Vec3<double> k1, k2, k3, k4, cm, x1, x2, x1ll, x2ll, grad, g;
    double t, h, m1, m2;

    if(!part1.Alive() || !part2.Alive()) return;

    if(!part1.IsStarted() || !part2.IsStarted()) return;

    //  part1.ApplyBoundaryCondition(mask, gpm, gpa);
    //  part2.ApplyBoundaryCondition(mask, gpm, gpa);


    t = global.Time(); // tempo corrente
    h = global.dt();

    x1 = part1.Xv();
    x2 = part2.Xv();

    x1ll = part1.LongLatZ();
    x2ll = part2.LongLatZ();

    //    if(mask.Value(x1ll) < global.kt2d || mask.Value(x2ll) < global.kt2d) return;

    cm = ModifiedCM();

    // qui inizio Runge-Kutta al quarto ordine
    
    try{
      //      global.er << "CoupleLagrangianTracer:EvolveAnalytical\n";
      m1 = mask->Value(x1ll, part1.vicinix);     
      m2 = mask->Value(x2ll, part2.vicinix);
    }
    catch(int ercode)
      {
	cout << "Lanciata da EvolveCinematico2D: mask " << ercode << endl;
	global.er << "Lanciata da EvolveCinematico2D: mask " << ercode << endl;
	global.er << "sul vettore ";
	x1ll.display(global.er);	global.er << " o ";
	x2ll.display(global.er);
	global.er << endl << endl;
	m1 = m2 = 0.;
      }

   if(m1 < global.kt2d || m2 < global.kt2d)   
    {
	try{
	  grad = mask->Gradient(x1ll);
	}
	catch(int ercode)
	  {
	    global.er << "Lanciata da EvolveCinematico2D: mask.Gradient " << ercode << endl;
	    global.er << "sul vettore "; x1ll.display(global.er);
	    global.er << endl << endl;
	    grad.zero();
	  }	

	g.setX(grad.Y() / R_earth * 360. / P2 );
	g.setY(grad.X() / R_earth * 360. / P2 / cos( x1ll.Y() * P2 / 360. ) );
	
	k1 = ( subgrid.u(x1 - cm, t ) * m1 + g * subgrid.psi(x1 - cm, t) ) * h;
	k2 = ( subgrid.u(x1 - cm + 0.5 * k1, t + h/2.) * m1 + g * subgrid.psi(x1 - cm + 0.5 * k1, t + h/2)) * h;
	k3 = ( subgrid.u(x1 - cm + 0.5 * k2, t + h/2.) * m1 + g * subgrid.psi(x1 - cm + 0.5 * k2, t + h/2)) * h;
	k4 = ( subgrid.u(x1 - cm + k3,t + h) * m1 + g * subgrid.psi(x1 - cm + 0.5 * k3, t + h)) * h;
	
	part1.IncrSet(k1/6. + k2/3. + k3/3. + k4/6.);
	
	try{
	  grad = mask->Gradient(x2ll);
	}
	catch(int ercode)
	  {
	    global.er << "Lanciata da EvolveCinematico2D: mask.Gradient " << ercode << endl;
	    global.er << "sul vettore ";
	    x2ll.display(global.er);
	    global.er << endl << endl;
	    grad.zero();
	  }

	g.setX(grad.Y() / R_earth * 360. / P2 );
	g.setY(grad.X() / R_earth * 360. / P2 / cos( x1ll.Y() * P2 / 360. ) );
		
	k1 = ( subgrid.u(x2 - cm, t ) * m2 + g * subgrid.psi(x2 - cm, t) ) * h;
	k2 = ( subgrid.u(x2 - cm + 0.5 * k1, t + h/2.) * m2 + g * subgrid.psi(x2 - cm + 0.5 * k1, t + h/2)) * h;
	k3 = ( subgrid.u(x2 - cm + 0.5 * k2, t + h/2.) * m2 + g * subgrid.psi(x2 - cm + 0.5 * k2, t + h/2)) * h;
	k4 = ( subgrid.u(x2 - cm + k3,t + h) * m2 + g * subgrid.psi(x2 - cm + 0.5 * k3, t + h)) * h;

	part2.IncrSet(k1/6. + k2/3. + k3/3. + k4/6.);
      } 

    // qui inizio Runge-Kutta al quarto ordine

    else{
      
      k1 = subgrid.u(x1 - cm, t ) * h;
      k2 = subgrid.u(x1 - cm + 0.5 * k1, t + h/2.) * h;
      k3 = subgrid.u(x1 - cm + 0.5 * k2, t + h/2.) * h;
      k4 = subgrid.u(x1 - cm + k3,t + h) * h;
      
      part1.IncrSet(k1/6. + k2/3. + k3/3. + k4/6.);
      
      k1 = subgrid.u(x2 - cm, t ) * h;
      k2 = subgrid.u(x2 - cm + 0.5 * k1, t + h/2.) * h;
      k3 = subgrid.u(x2 - cm + 0.5 * k2, t + h/2.) * h;
      k4 = subgrid.u(x2 - cm + k3,t + h) * h;    
      
      part2.IncrSet(k1/6. + k2/3. + k3/3. + k4/6.);

    }

   part1.ApplyBoundaryCondition(mask, gpm, gpa);
   part2.ApplyBoundaryCondition(mask, gpm, gpa);
      

  }


  //  void EvolveCampoRisolto(class VectorField3D& field_old, class VectorField3D& field_new, class ScalarField3D& mask)
void EvolveCampoRisolto(class VectorField3D& field_old, class VectorField3D& field_new, class ScalarField3D_Mask* mask, class ScalarField3D* gpm, class ScalarField3D* gpa)
  {
    part1.EvolveCampoRisolto(field_old, field_new, mask, gpm, gpa);
    part2.EvolveCampoRisolto(field_old, field_new, mask, gpm, gpa);
  }

  inline double Age(){ return part1.Age();}
  inline void Age(double xx){ part1.Age(xx); part2.Age(xx);}
  inline bool IsAnchovy(){ return part1.IsAnchovy();}
  
  inline bool ToBePrinted()
  {
    if(!part1.IsStarted() || !part2.IsStarted()) return false;
    if(!part1.Alive() || !part2.Alive()) return false; 
    return true;
  }

};



std::ostream& operator << (std::ostream& ostrm, CoupleLagrangianTracer C)
{

  Vec3<double> l1, l2;

  l1 = C.part1.LongLatZ();
  l2 = C.part2.LongLatZ();

      //      ostrm << global.Time() << " ";
  l1.display(ostrm);
  l2.display(ostrm);
 
  //  ostrm <<" " << C.part1.Age() << " " << C.part1.Start() << endl;
  ostrm <<" " << C.part1.Age() << endl;

  return ostrm;
}

/*
Sostituisce la funzione già presente in CoupleLagrangianTracer.C
Se non ci sono le colonne con Age e Start assegna valori di default a
queste due variabili (righe 21-30)
*/

std::istream& operator >> (std::istream& istream, CoupleLagrangianTracer& C)
{
  double x1,y1,z1,x2,y2,z2,age,start; // \RC
  
  std::string line, pippo; 
  std::istringstream ss;

  // Ci sono sempre almeno 6 valori da leggere (coordinate delle coppie)
  istream >> x1; istream >> y1; istream >> z1;
  istream >> x2; istream >> y2; istream >> z2;

  /*  if(!ss.eof()) // se ci sono anche le colonne 7 e 8
    {
      ss >> start;
      ss >> age;
      age *= 86400.;
    }
  else 
    {
      //      start = 0.;
      // age = 0.; i defaults li stabilisce il costruttore
      }*/


  Vec3<double> l1(x1,y1,z1), l2(x2,y2,z2);
 
  C.part1.SetLongLatZ(l1);
  C.part2.SetLongLatZ(l2);
  C.part1.Start(start);
  C.part2.Start(start);
  C.part1.Age(age);
  C.part2.Age(age);

  return istream;
}


/*std::istream& operator >> (std::istream& istream, CoupleLagrangianTracer& C)
{
  double x,y,z,x1,y1,z1;

  istream >> x;   istream >> y;   istream >> z;
  istream >> x1;  istream >> y1;  istream >> z1;

  //  cout << x << " " << " " << y << " " << z << endl;

  Vec3<double> l1(x,y,z), l2(x1,y1,z1);
 
  C.part1.SetLongLatZ(l1);
  C.part2.SetLongLatZ(l2);

  return istream;
  }*/
