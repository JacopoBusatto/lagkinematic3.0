
/*

Post Processer 
versione ottimizzata per l'utilizzo della memoria

*/

#include <iostream>
#include <stdio.h>
#include <string.h>
#include <cstdlib>
#include <fstream>

#include "FieldAbstract.C"
#include "VectorField.C"

#define R_earth 6378137 // equatoriale

#include "randoknuth.c"

#include "KinematicModel.C"
#include "LagrangianTracer.C"
#include "CoupleLagrangianTracer2.C"


// Funzioni usate da PostProcesser

size_t Nalive(CoupleLagrangianTracer **P, size_t i_time, size_t NTOT)
{
  size_t count;

  count = 0;

  for(size_t i=0; i<NTOT; i++) 
    {
      if( (*(P[i_time]+i)).part1.Alive() ) count++;
    }
  return count;
}


Vec3<double> Xavg(CoupleLagrangianTracer **P, size_t i_time, size_t NTOT, size_t c)
{
  Vec3<double> xavg; // azzerato di default
  double count;

  count = 0.;

  if(c==1)
    {
      for(size_t i=0; i<NTOT; i++)
	{
	  if( (*(P[i_time]+i)).part1.Alive() ) 
	    {
	      count++;
	      xavg += (*(P[i_time]+i)).part1.LongLatZ(); 
	    }
	}
    }

  if(c==2)
    {
      for(size_t i=0; i<NTOT; i++)
	{
	  if( (*(P[i_time]+i)).part2.Alive() ) 
	    {
	      count++;
	      xavg += (*(P[i_time]+i)).part2.LongLatZ(); 
	    }
	}
    }


  return (1./count) * xavg;
}


Vec3<double> Xavg(CoupleLagrangianTracer **P, size_t i_time, size_t NTOT)
{
  Vec3<double> x;

  x = Xavg(P, i_time, NTOT,1) + Xavg(P, i_time, NTOT,2);
  return 0.5 * x;
}

MatDoub Cov(CoupleLagrangianTracer **P, size_t i_time, size_t NTOT)
{ 
  Vec3<double> xm = Xavg(P, i_time, NTOT), x1, x2;
  MatDoub C(3,3);
  double e, count1, count2;

  count1 = count2 = 0.;

  for(size_t k=0; k<NTOT; k++)
    {  
      if( (*(P[i_time]+k)).part1.Alive() )
	{
	  x1 = (*(P[i_time]+k)).part1.LongLatZ(); 
	  C[0][0] += (x1.X() - xm.X()) * (x1.X() - xm.X());
	  C[1][1] += (x1.Y() - xm.Y()) * (x1.Y() - xm.Y());
	  C[0][1] += (x1.Y() - xm.Y()) * (x1.X() - xm.X());
	  C[2][2] += (x1.Z() - xm.Z()) * (x1.Z() - xm.Z());
	  count1++;
	}

      if( (*(P[i_time]+k)).part2.Alive() )
	{
	  x2 = (*(P[i_time]+k)).part2.LongLatZ(); 
	  C[0][0] += (x2.X() - xm.X()) * (x2.X() - xm.X());
	  C[1][1] += (x2.Y() - xm.Y()) * (x2.Y() - xm.Y());
	  C[0][1] += (x2.Y() - xm.Y()) * (x2.X() - xm.X());
	  C[2][2] += (x1.Z() - xm.Z()) * (x1.Z() - xm.Z());
	  count1++;
	}
    }
  
  C[0][0] /= count1; C[1][1] /= count1;  C[2][2] /= count1;
  C[0][1] /= count1;

    C[1][0] = C[0][1];

  return C;
}


void Movie(CoupleLagrangianTracer **P, long int *indice_t,double *starttime, double *endtime, size_t NTOT, double timestep, double TMAX  ,const char* nomedir)
{
  size_t nsnap, TE,j;
  ofstream gnuplot, snap; 
  char nomefile[500];
  Vec3<double> x1, x2;
 

  sprintf(nomefile,"%s/grafico.gnu",nomedir);
  gnuplot.open(nomefile);

  gnuplot << "set xra [-20:50]\n set yra[-50:-25]\n";
  gnuplot << "set terminal gif small animate optimize" << endl;
  gnuplot << "set output 'grafico.gif'" << endl;
  gnuplot << "set palette rgb 33,13,10;" << endl;
  gnuplot << "set cbrange [0:500]" << endl;
  gnuplot << "unset key" << endl;
  gnuplot << "set size ratio -1" << endl;

  long int snapt;

  for(double time = 0; time<TMAX; time += timestep)
    {

      snapt = time / timestep;
      sprintf(nomefile,"%s/snap%ld.dat",nomedir, snapt);
      snap.open(nomefile);

      gnuplot << "plot  \"" << nomefile << "\" u 1:2:7 title \"" << TE  << "\" w d  lc palette,\"\" u 4:5:7 notitle w d  lc palette,'/archivio/palatella/world_10m.txt' u ($1<0?360+$1:$1):2  w d" << endl;

      for(size_t i=0; i<NTOT; i++)
	{
	  if(time >= starttime[i] && time < endtime[i])
	    {
	      j = (time-starttime[i]) / timestep;
	      if(j<0 || j> indice_t[i])
		{
		  cout << "errore particella " << i << endl;
		  //throw int(137);
		}

	      if( (*(P[i]+j)).part1.Alive() &&  (*(P[i]+j)).part1.Alive())
		{
		  x1 = (*(P[i]+j)).part1.LongLatZ(); 
		  x2 = (*(P[i]+j)).part2.LongLatZ(); 
		  x1.display(snap); x2.display(snap);	      
		  snap << (*(P[i]+j)).part1.Age();
		  snap << endl;
		} 
	    }
	}
      snap.close();
    }

  gnuplot.close();
}


void DensityAge(CoupleLagrangianTracer **P, size_t age, size_t Nbinx, \
		size_t Nbiny ,size_t NTOT, size_t MaxTime, const char* nomedir, bool coast)
{
  /* density_time e' un intero = giorni * 86400 / global.snap*/
  ofstream density;
  char nomefile[500];
  Vec3<double> x1, x2;
  double xxc, yyc;
  size_t snap_time, I, J;
  MatDoub PDF(Nbinx, Nbiny);
  double dx, dy, lon_max = 37., lon_min=-6., lat_min=30., lat_max = 46., eta;
  double *xcoast, *ycoast;

  dx = (lon_max - lon_min) / Nbinx;
  dy = (lat_max - lat_min) / Nbiny;


  long jj, nc = 0;

  if(coast)
    {
      ifstream coste;
      coste.open("/archivio/palatella/medcoastBig.dat");

      while(!coste.eof())
	{
	  coste >> xxc >> yyc;
	  nc++;
	} 

      coste.close();

      xcoast = new double[nc];
      ycoast = new double[nc];

      coste.open("/archivio/palatella/medcoastBig.dat");
      
      for(jj=0; jj<nc; jj++)
	{
	  coste >> xxc >> yyc;
	  xcoast[jj] = xxc;
	  ycoast[jj] = yyc;
	}
      
      coste.close();

    }
  


  sprintf(nomefile,"%s/densityAge_%ld.dat",nomedir,age);
  density.open(nomefile);

  for(size_t i=0; i<Nbinx; i++)
    for(size_t j=0; j<Nbiny; j++)
      PDF[i][j]=0.;

  for(size_t k=0; k<NTOT; k++)
    {
      for(size_t run_time =0 ; run_time< MaxTime; run_time++)
	{
	  eta = -int( (*(P[run_time]+k)).part1.Age() );

	  //	  if(k<3) cout << eta <<" " << age << " " << (*(P[run_time]+k)).part1.Age() << endl;

	  if( eta == age)
	    {
	      if( (*(P[run_time]+k)).part1.Alive())
		{
		  x1 = (*(P[run_time]+k)).part1.LongLatZ();
		  
		  I = (x1.X() - lon_min) / dx;
		  J = (x1.Y() - lat_min) / dy;
		  
		  if(I<0 || I>=Nbinx || J<0 || J>= Nbiny)
		    cout << I << " " << J << " " << x1.X() << " " << x1.Y() << " " << k << endl;
		  else
		    {
		      if(!coast)
			PDF[I][J]++;
		      if(coast)
			{
			  jj=0;
			  while(jj<nc && ( fabs( x1.X() - xcoast[jj] ) > 1*dx || fabs( x1.Y() - ycoast[jj] ) > 1*dy)) jj++;
			  if(jj != nc)
			    PDF[I][J]++;
			}
		    }
		}

	      if( (*(P[run_time]+k)).part2.Alive())
                {
                  x1 = (*(P[run_time]+k)).part2.LongLatZ();

                  I = (x1.X() - lon_min) / dx;
                  J = (x1.Y() - lat_min) / dy;

                  if(I<0 || I>=Nbinx || J<0 || J>= Nbiny)
                    cout << I << " " << J << " " << x1.X() << " " << x1.Y() << " " << k << endl;
		  else
		    {
		      if(!coast)
                        PDF[I][J]++;
                      if(coast)
                        {
                          jj=0;
                          while(jj<nc && ( fabs( x1.X() - xcoast[jj] ) > 1*dx || fabs( x1.Y() - ycoast[jj] ) > 1*dy)) jj++;
                          if(jj != nc)
                            PDF[I][J]++;
                        }                      
                    }
                }
	    }
	}
    }
     
  for(size_t i=0; i<Nbinx; i++)
	{
	  for(size_t j=0; j<Nbiny; j++)
	    density << lon_min + (i+0.5)*dx << " " << lat_min + (j+0.5)*dy << " " << PDF[i][j]/2./NTOT << endl;
	  density << endl;
	}
}


void Density(CoupleLagrangianTracer **P, size_t density_time, size_t Nbinx, \
	     size_t Nbiny ,size_t NTOT, const char* nomedir)
{
  /* density_time e' un intero = giorni * 86400 / global.snap*/
  ofstream density; 
  char nomefile[500];
  Vec3<double> x1, x2;
  size_t snap_time, I, J;
  MatDoub PDF(Nbinx, Nbiny);
  double dx, dy, lon_max = 37., lon_min=-6., lat_min=30., lat_max = 46.;

  dx = (lon_max - lon_min) / Nbinx;
  dy = (lat_max - lat_min) / Nbiny;

  sprintf(nomefile,"%s/density_%ld.dat",nomedir,density_time);
  density.open(nomefile);

  for(size_t i=0; i<Nbinx; i++)
    for(size_t j=0; j<Nbiny; j++)
      PDF[i][j]=0.;
      
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


void Movie3D(CoupleLagrangianTracer **P, double *time,size_t NTIME, size_t snap_time, size_t NTOT, const char* nomedir)
{
  size_t nsnap, TE;
  ofstream gnuplot, snap; 
  char nomefile[500];
  Vec3<double> x1, x2;

  sprintf(nomefile,"%s/grafico3D.gnu",nomedir);
  gnuplot.open(nomefile);

  gnuplot << "set xra [-1:4]\n set yra[36:38]\n";
  gnuplot << "set terminal gif small animate optimize" << endl;
  gnuplot << "set pm3d map" << endl;
  gnuplot << "set view 50,270,1.2" << endl;
  gnuplot << "set output 'grafico3D.gif'" << endl;
  gnuplot << "unset border ; unset xtics; unset ytics; unset ztics" << endl; 
  //  gnuplot << "rgb(x)=int(-x/4)" << endl;

  nsnap = NTIME / snap_time;

  for(size_t i=0; i<nsnap; i++)
    {
      sprintf(nomefile,"%s/snap%ld.dat",nomedir,i);
      snap.open(nomefile);
      
      TE = time[i * snap_time]/86400.;

      gnuplot << "splot \"../batimetria.dat\" u 1:2:($3>0?0:$3) notitle w l lt -1,\"" << nomefile << "\" u 1:2:3:3 title \"" << TE  << "\" w p pt 7 ps 0.7 lc palette,\"\" u 4:5:6:6 notitle w p pt 7 ps 0.7 lc palette" << endl;
 
      for(size_t k=0; k<NTOT; k++)
	{  
	  if( (*(P[i * snap_time]+k)).part1.Alive() &&  (*(P[i * snap_time]+k)).part2.Alive() )
	    {
	      x1 = (*(P[i * snap_time]+k)).part1.LongLatZ(); 
	      x2 = (*(P[i * snap_time]+k)).part2.LongLatZ(); 

	      x1.display(snap); x2.display(snap); snap << endl;
	    }
	}
      snap.close();
    }

  gnuplot.close();
}


bool IsNearStation(Vec3<double> X, Vec3<double> X0, double distanza)
{
  Vec3<double> d = X - X0;

  d.setZ(0.);
  
  if(d.magnitude() < distanza) return true;
  else return false;
}

VecDoub PDF(CoupleLagrangianTracer **P, double *time,size_t NTIME, size_t NTOT, Vec3<double> S, size_t bin)
{
  VecDoub pdf(NTIME);
  Vec3<double> X1, X2;
  size_t count = 0, i;

  for(size_t k=0; k<NTIME; k++) pdf[k] = 0.;

  for(size_t k=0; k<NTOT; k++)
    {
      i=0;
      while( !IsNearStation( (*(P[i]+k)).part1.LongLatZ(), S, 0.5) &&  (*(P[i]+k)).part1.Alive() && i < NTIME - 1 ) i++;
      if( (*(P[i]+k)).part1.Alive() && i < NTIME -1 ) pdf[i/bin]++; 

      i=0;
      while( !IsNearStation( (*(P[i]+k)).part2.LongLatZ(), S, 0.5) &&  (*(P[i]+k)).part2.Alive() && i < NTIME - 1) i++;
      if( (*(P[i]+k)).part2.Alive() && i < NTIME - 1 ) pdf[i/bin]++;
    }

  /*  for(size_t k=0; k<NTIME; k++) count += pdf[k];

  if(count)
  for(size_t k=0; k<NTIME; k++) pdf[k] = pdf[k] / count;*/

  return pdf;

}


VecDoub ATL(CoupleLagrangianTracer **P, double *time,size_t NTIME, size_t snap_time, size_t NTOT, double lat_ATL, double lon_ATL)
{
  //VecDoub InATL(NTIME), flag(2*NTOT);
  VecDoub InATL(NTIME), flag1(NTOT), flag2(NTOT);
  Vec3<double> x1, x2;
  size_t count = 0, i;

  for(size_t i=0; i<NTIME; i++) InATL[i] = 0;
  //for(size_t k=0; k<2 * NTOT; k++) flag[k] = 0;
  for(size_t k=0; k< NTOT; k++) {
    flag1[k] = 0;
    flag2[k] = 0;
  }


  for(size_t i=0; i<NTIME; i++) 
    {
      for(size_t k=0; k<NTOT; k++)
	{
	  if( (*(P[i * snap_time]+k)).part1.Alive())
	    {
	      x1 = (*(P[i * snap_time]+k)).part1.LongLatZ();
	      if(x1.Y() > lat_ATL && x1.X() < lon_ATL && flag1[k] == 0) 
		{
		  flag1[k] = 1;
		  InATL[i]++;
		}
	    }

	  if( (*(P[i * snap_time]+k)).part2.Alive())
            {
              x2 = (*(P[i * snap_time]+k)).part2.LongLatZ();
              if(x2.Y() > lat_ATL && x2.X() < lon_ATL && flag2[k] == 0) 
		{
		  flag2[k] = 1;
		  InATL[i]++;
		}
            }

	}
    }

  return InATL;
}


VecDoub CrossATL(CoupleLagrangianTracer **P, double *time,size_t NTIME, size_t snap_time, size_t NTOT, double lat_ATL, double lon_ATL)
{
  VecDoub InATL(NTIME);
  Vec3<double> x1, x1old, x2, x2old;
  size_t count = 0, i;

  for(size_t i=0; i<NTIME; i++) InATL[i] = 0;

  for(size_t i=1; i<NTIME; i++)
    {
      for(size_t k=0; k<NTOT; k++)
        {
	  //          if( (*(P[i * snap_time]+k)).part1.Alive() &&  (*(P[i * snap_time]+k)).part2.Alive() )
          if( (*(P[i * snap_time]+k)).part1.Alive())
            {
              x1 = (*(P[i * snap_time]+k)).part1.LongLatZ();
              x1old = (*(P[(i - 1) * snap_time]+k)).part1.LongLatZ();
	      if( (x1.Y() > lat_ATL && x1old.Y() <= lat_ATL) && x1.X() < lon_ATL ) InATL[i]++;
            }

	  if( (*(P[i * snap_time]+k)).part2.Alive())
            {
              x2 = (*(P[i * snap_time]+k)).part2.LongLatZ();
              x2old = (*(P[(i - 1) * snap_time]+k)).part2.LongLatZ();
	      if( (x2.Y() > lat_ATL && x2old.Y() <= lat_ATL) && x2.X() < lon_ATL ) InATL[i]++;
            }

        }
    }

  return InATL;
}


//********************************************************************************
// 
//********************************************************************************

int main(int argc, char *argv[])
{

  if(!argv[1])
    {
      cout << "\n\n# PostProcesser 1.0\n\n " << endl;
      printf("\nUsage %s <nome directory simulazione> <numero particelle> <numero processori> <file parametri> [movie snap] <nome out file>\n\n", argv[0]);
      return 0;
    }

  size_t part, Nline, n, Npart, np, NTOT, k, NTIME;
  double time_max, time_old, age, time;
  char nomefile[500];
  CoupleLagrangianTracer dummy;
  string line;
  size_t* linexfile;
  long nt;
  char *in_directory = NULL;
  char *out_file     = NULL;
  VecDoub InATL(NTIME), cross(NTIME);
  in_directory = argv[1];
  out_file     = argv[6];
  size_t i_part, count, nnt;

  sscanf(argv[2],"%ld", &NTOT);
  sscanf(argv[3],"%ld", &np);

  Npart = NTOT / np; // particelle per processore
  
  for(size_t tt = 0; tt < NTIME; tt++) {
    InATL[tt]=0;
    cross[tt]=0;
  }
 
  linexfile = new size_t[np];
  for(size_t i = 0; i<np; i++) linexfile[i]=0;

  ifstream *filestream = new ifstream[np];
  ifstream logstream(argv[4]);

  logstream >> global;
  logstream.close();

  global.log.open("paperino.log");

  Nline = 0;
  time_max = 0.;

  long startfile, movie_snap;

  if(argv[5])
    sscanf(argv[5],"%ld",&movie_snap);
  else movie_snap = 1;

  startfile = 0;

  double *starttime = new double[NTOT]; // l'array con i tempi di start
  double *endtime = new double[NTOT]; // l'array con i tempi di start
  size_t *alive = new size_t[NTOT];

  for(size_t i=0;i<NTOT;i++)
    {
      starttime[i] = endtime[i] = -1.; // valore di default
      alive[i] = 0;
    }

  double timeold;

  timeold = -1;


  for(size_t i=0; i<np; i++) 
    {
      if (i + startfile < 100)
	sprintf(nomefile,"%s/trajectories_LongLatZ_%.2ld.dat",argv[1],i + startfile);
      else
	sprintf(nomefile,"%s/trajectories_LongLatZ_%.3ld.dat",argv[1],i + startfile);
      cout << "# apro il file " << nomefile << endl;
      filestream[i].open(nomefile);

      while(filestream[i] >> time >> part >> dummy >> age)
	{ 
	  //	  cout << Nline << endl;
	  i_part = part + i * Npart;	
	  if(starttime[i_part]<0) starttime[i_part] = time;
	  if(starttime[i_part]>0) endtime[i_part] = time; 

	  //  cout << i_part << " " << starttime[i_part] << " " << endtime[i_part] << endl;

	  Nline++;
	  linexfile[i]++;
	  if(time > time_max) time_max = time;
	}
      filestream[i].close();
    }
   
  cout << "# nei file ci sono Nline = " << Nline << endl;

  NTIME = time_max/global.snap; // massimo valore di NTIME

  NTIME /= movie_snap;
  NTIME+=5; // per sicurezza

  cout << "# NTOT = " << NTOT << endl;
  cout << "# NTIME = " << NTIME << endl;


  long int NNEW;
  long int*  indice_t = new long int[NTOT];

  NNEW = 0;

  for(size_t i=0;i<NTOT;i++) 
    {
      indice_t[i] = ( endtime[i] - starttime[i] ) / global.snap / movie_snap;    
      if(indice_t[i]>0)
	{
	  NNEW += indice_t[i];
	  cout << i << " " << starttime[i] << " " << endtime[i] << " " << indice_t[i] << endl;
	}
    }

  long int STIMA;

  STIMA = NNEW * sizeof( CoupleLagrangianTracer  ); 

  cout << "NNEW = " << NNEW << " Stima RAM necessaria = " << STIMA / pow(2,30) << " Gb\n"; /* boh, tutta da controllre */

  CoupleLagrangianTracer** P = new CoupleLagrangianTracer*[NTOT];

  if(!P)
    {
      cout << "Memoria insufficiente " << endl;
      cout.flush();
      return 1;
    }

  double *Time = new double[NTIME];

  cout << "Allocando...\n\n";
  cout.flush();
 

  for(size_t i=0; i<NTOT; i++)
    {
      cout << i << " " << indice_t[i] << endl;
      if(indice_t[i]>0)
	{
	  try{
	    P[i] = new CoupleLagrangianTracer[indice_t[i]];
	  }
	  catch (std::bad_alloc& ba)
	    {
	      cout << "Memoria insufficiente alla particella : " << i << " / " << NTOT << " "  << ba.what() << '\n';
	      cout.flush();
	      throw int(123);
	    }
	}
      
      // Time[i] = -1; // valore non inizializzato
    }


  /*  for(size_t j=0; j<NTOT; j++)
    for(size_t i=0; i<indice_t[j]; i++)
      { 
	(*(P[j] + i)).part1.SetDead();
	(*(P[j] + i)).part2.SetDead();
	}*/
    
 
  cout << "Allocate le traiettorie con successo\n\n";
  cout.flush();

  nt = k = 0; time = time_old = 0.;
  long int ttime;
 

  for(size_t i=0; i<np; i++) 
    {
      count = 0;
      if (i + startfile < 100)
	sprintf(nomefile,"%s/trajectories_LongLatZ_%.2ld.dat",argv[1],i + startfile);
      else
	sprintf(nomefile,"%s/trajectories_LongLatZ_%.3ld.dat",argv[1],i + startfile);
      filestream[i].open(nomefile);
      cout << "# Riapro  " << nomefile << endl;

      nnt = 0;

      while(filestream[i] >> time >> part >> dummy >> age) 
      {	
	i_part = part + i * Npart;	
	ttime = ( time - starttime[i] ) / global.snap / movie_snap;    

	if(ttime > 0 && ttime < indice_t[i_part])
	  {
	    *(P[i_part] + ttime) = dummy;
	    (*(P[i_part] + ttime)).Age(age * 86400.);
	    (*(P[i_part] + ttime)).part1.SetAlive();
	    (*(P[i_part] + ttime)).part2.SetAlive();
	  }
      }
      filestream[i].close();
    }

  cout << "# lettura terminata con successo\n " << endl;

  //cout << "# Faccio il movie....\n " << endl;

  //  Movie(P, indice_t, NTOT, 1, argv[1]);

  Movie(P,indice_t, starttime, endtime, NTOT, global.snap*movie_snap, time_max, argv[1]);

  /*  cout << "# Controllo chi è in Atlantico...\n " << endl;
  InATL =      ATL(P, Time, NTIME,1,  NTOT, -30, 20);
  cout << "# Controllo chi ha attraversato la separazione per l'Atlantico...\n " << endl;
  cross = CrossATL(P, Time, NTIME, 1, NTOT, -30, 20);

  in_directory = strcat(in_directory,"/");
  ofstream atl( strcat(in_directory, out_file));

  for(size_t tt = 0; tt<NTIME;tt++)
    atl << tt << " " << InATL[tt] << " " << cross[tt] << endl;


    atl.close();*/

  return 0;
} 

