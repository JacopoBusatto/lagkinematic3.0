namespace palatella_math_functions
{

  double CerchioMassimo(double long1, double lat1, double long2, double lat2)
  {
    const double g2rad = 0.0174532925199433; // converte gradi in radianti
    const double unosuR = 0.00015696123057604771; // l'inverso del raggio della terra
    const double soglia = 1.5;

    double phi1, phi2, theta1, theta2, dist;
    
    //    if(fabs(phi1 - phi2) > soglia || fabs(theta2 - theta1) > soglia) return 0.;

    phi1 = g2rad * long1;
    phi2 = g2rad * long2;
    
    theta1 = g2rad * lat1;
    theta2 = g2rad * lat2;
	
    dist = cos(phi1 - phi2) * cos(theta1) * cos(theta2) + sin(theta1) * sin(theta2);
    
    if(fabs(dist) >=  1.) return 1.e6; // per motivi numerici quando coincidono a volte restituisce > 1
    return  unosuR / acos(dist); // restituisce in km^-1 
    
  }

  double CerchioMassimo(double long1, double lat1, double z1, double long2, double lat2, double z2) // z in metri
  {
    double d2, d3, dz = 1000. * (z2-z1);
    
    d2 = CerchioMassimo(long1, lat1, long2, lat2);
    
    //cout << long1-long2 << " " << lat1 - lat2 << " " <<   d2 << " " << dz << endl;

    return sqrt(dz*dz + d2*d2);
  }
  

  double CerchioMassimo(Vec3<double>& X1, Vec3<double>& X2)
  {
    return CerchioMassimo(X1.X(), X1.Y(), X1.Z(), X2.X(), X2.Y(), X2.Z());
  }

  double CerchioMassimo(const Doub *p1,const  Doub *p2)
  {
    //   return CerchioMassimo(p1[0], p1[1], p1[2], p2[0], p2[1], p2[2]);
   
      return CerchioMassimo(p1[0], p1[1], p2[0], p2[1]);
  }




  
}
