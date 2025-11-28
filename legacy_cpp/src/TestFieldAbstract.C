#include <iostream>
#include <stdio.h>
#include <string.h>
#include <cstdlib>
#include <fstream>

#include "FieldAbstract.C"
#include "VectorField.C"

int main(int argc, char *argv[])
{
  class ScalarField3D* field1 = Switcher(MFS_MyOcean);
  class ScalarField3D* field2 = Switcher(WRF);

  Lista.ocean = true;

  cout << "Carico il file tipo MFS_MyOcean\n\n"; cout.flush();
 
  field1->Init(argv[1],"lon","lat","depth","NULL");
  field1->Load(argv[1],"vozocrtx",0,"NULL");

  double x,y,z;

  for(long i=0; i<20; i++)
    {
      x = 10. + 1. * rand() / RAND_MAX * 2.;
      y = 36 + 1. * rand() / RAND_MAX * 2.;
      z = -20 + 1. * rand() / RAND_MAX;

      cout << x << " " << y << " " << z << " " << field1->Value(x,y,z) << endl;
    }

  //\RC
  class ScalarField3D* field2D = Switcher(OceanColor_MyOcean2D);

  cout << "Carico il file tipo OceanColor_MyOcean 2D\n\n"; cout.flush();
 
  field2D->Init(argv[2],"lon","lat","depth","NULL");
  cout << "# Ho inizializzato per CHL\n";
  field2D->Load(argv[2],"CHL",0,"NULL");

  for(long i=0; i<20; i++)
    {
      x = 10. + 1. * rand() / RAND_MAX * 2.;
      y = 36 + 1. * rand() / RAND_MAX * 2.;

      cout << x << " " << y << " " << field2D->Value(x,y,z) << endl;
    }
  // Fine \RC

  /*
  Lista.ocean = false;
  global.sigma = true;

  cout << "Carico il file tipo WRF\n\n"; cout.flush();
 
  global.er.open("pippo.er");

  field2->Init(argv[3],"XLONG_U","XLAT_U","ZNU","U10");
  field2->Load(argv[3],"U",0,"U10");

  for(long i=0; i<20; i++)
    {
      x = 10. + 1. * rand() / RAND_MAX * 2.;
      y = 36 + 1. * rand() / RAND_MAX * 2.;
      z = 0. + .99 * rand() / RAND_MAX;

      cout << x << " " << y << " " << z << " " << field2->Value(x,y,z) << endl;
    }

  global.er.close();

  cout << "\nTEST MASK \n\n";

  class ScalarField3D_Mask* pippo;

  pippo = new ScalarField3D_Mask_WRF;

  Lista.orografia = false;

  pippo -> Init(argv[3],"XLONG","XLAT","ZNW","NULL");
  pippo -> Mask(argv[3],"W");

  for(long i=0; i<20; i++)
    {
      x = 10. + 1. * rand() / RAND_MAX * 2.;
      y = 36 + 1. * rand() / RAND_MAX * 2.;
      z = 0. + .99 * rand() / RAND_MAX;

      cout << x << " " << y << " " << z << " " << pippo->Value(x,y,z) << endl;
    }

  cout << "\nTEST OROGRAPHY \n\n";

  pippo->LoadOrografia(argv[3],"HGT");

  Vec3<double> p;

  */

  /*

  ofstream oro("oro");

  for(x=7; x<18; x+=0.2)
    {
      for(y=37; y<46; y+=0.2)
	{
	  p.setX(x);
	  p.setY(y);
	  p.setZ(0.2);
	  
	  oro << x << " " << y << " " << pippo->bottom(p) << " estra"<< endl;
	}
      oro << "\n";
    }
  
    oro.close(); */

  VectorField3D U(MFS_MyOcean);

}
