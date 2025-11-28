	
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
