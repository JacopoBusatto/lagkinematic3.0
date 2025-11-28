void DensityT(CoupleLagrangianTracer **P, size_t density_time, size_t Nbinx, \
	      size_t Nbiny ,size_t NTOT, size_t NTIME, const char* nomedir)
{
  /* density_time e' un intero = giorni * 86400 / global.snap*/
  ofstream density; 
  char nomefile[500];
  Vec3<double> x1, x2;
  size_t snap_time, I, J;
  MatDoub PDF(Nbinx, Nbiny);
  double dx, dy, lon_max = 37., lon_min=-6., lat_min=30., lat_max = 46.;
  double T1, T2, T3;

  for(size_t k=0; k<NTOT; k++)
    {
      for(size_t itime=0; itime < NTIME - 2; itime++)
	{
	  T1 = (*(P[itime]+k)).part1.T(); 
	  T2 = (*(P[itime+1]+k)).part1.T(); 
	  T3 = (*(P[itime+2]+k)).part1.T();

	  if(T1 < 15. && T2 < 15. && T3 < 15.)
	    {
	      if((itime + k)%2) (*(P[itime]+k)).part1.SetDead(); 
	    }

	  T1 = (*(P[itime]+k)).part2.T(); 
	  T2 = (*(P[itime+1]+k)).part2.T(); 
	  T3 = (*(P[itime+2]+k)).part2.T();

	  if(T1 < 15. && T2 < 15. && T3 < 15.)
	    {
	      if((itime + k)%2) (*(P[itime]+k)).part2.SetDead(); 
	    }
	}
    }


  dx = (lon_max - lon_min) / Nbinx;
  dy = (lat_max - lat_min) / Nbiny;

  sprintf(nomefile,"%s/densityT_%ld.dat",nomedir,density_time);
  density.open(nomefile);

  for(size_t k=0; k<NTOT; k++)
    {  
      if( (*(P[density_time]+k)).part1.Alive() )
	   {
	     x1 = (*(P[density_time]+k)).part1.LongLatZ(); 
	     
	     I = (x1.X() - lon_min) / dx;
	     J = (x1.Y() - lat_min) / dy;
	     
	     if(I<0 || I>=Nbinx || J<0 || J>= Nbiny)
	       cout << I << " " << J << " " << x1.X() << " " << x1.Y() << " " << k << endl;
	     PDF[I][J]++;
	  
	   }
      if( (*(P[density_time]+k)).part2.Alive() )
	{
	  x2 = (*(P[density_time]+k)).part2.LongLatZ(); 
	  
	  I = (x2.X() - lon_min) / dx;
	  J = (x2.Y() - lat_min) / dy;

	  if(I<0 || I>=Nbinx || J<0 || J>= Nbiny)
	     cout << I << " " << J << " " << x2.X() << " " << x2.Y() << " " << k << endl;


	  PDF[I][J]++;
	}
    }

  for(size_t i=0; i<Nbinx; i++)
    {
      for(size_t j=0; j<Nbiny; j++)
	density << lon_min + (i+0.5)*dx << " " << lat_min + (j+0.5)*dy << " " << PDF[i][j]/2./NTOT << endl;
      density << endl;
    }

}
