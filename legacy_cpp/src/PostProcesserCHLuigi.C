
/*

Post Processer 

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


void Movie(CoupleLagrangianTracer **P, double *time,size_t NTIME, size_t snap_time, size_t NTOT, const char* nomedir)
{
  size_t nsnap, TE;
  ofstream gnuplot, snap; 
  char nomefile[500];
  Vec3<double> x1, x2;

  sprintf(nomefile,"%s/grafico.gnu",nomedir);
  gnuplot.open(nomefile);

  gnuplot << "set xra [-5:32]\n set yra[33:46]\n";
  gnuplot << "set terminal gif small animate optimize" << endl;
  gnuplot << "set output 'grafico.gif'" << endl;
  gnuplot << "set palette rgb 33,13,10;" << endl;
  gnuplot << "set cbrange [-500:0]" << endl;
  gnuplot << "unset key" << endl;


  nsnap = NTIME / snap_time;

  for(size_t i=0; i<nsnap; i++)
    {
      //      sprintf(nomefile,"%s/snap%ld.dat",nomedir,i);
      sprintf(nomefile,"snap%ld.dat",i);
      snap.open(nomefile);
      
      //      cout << i << " " << endl;
      
      TE = time[i * snap_time]/86400.;

      gnuplot << "plot  \"" << nomefile << "\" u 1:2:3 title \"" << TE  << "\"w d  lc palette,\"\" u 4:5:6 notitle w d  lc palette,'/cantiere/palatella/medcoastBig.dat' w d" << endl;
 
      for(size_t k=0; k<NTOT; k++)
	{  
	  if( (*(P[i * snap_time]+k)).part1.Alive() &&  (*(P[i * snap_time]+k)).part2.Alive() )
	    {
	      x1 = (*(P[i * snap_time]+k)).part1.LongLatZ(); 
	      x2 = (*(P[i * snap_time]+k)).part2.LongLatZ(); 


	      //	      cout << x1.X() << " " << x1.Y() << endl;

	      x1.display(snap); x2.display(snap); snap << endl;
	    }
	}
      snap.close();
    }

  gnuplot.close();
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


//********************************************************************************
// 
//********************************************************************************

int main(int argc, char *argv[])
{

  if(!argv[1])
    {
      cout << "\n\n# PostProcesser 1.0\n\n " << endl;
      printf("\nUsage %s <nome directory simulazione> <numero particelle> <numero processori> <file parametri> [elencoNetCDF_CHL(opz)] [startFileTraj (opz)]\n\n", argv[0]);
      return 0;
    }

  size_t part, Nline, n, Npart, np, NTOT, k, NTIME;
  double time_max, time_old, age, time;
  char nomefile[500];
  CoupleLagrangianTracer dummy;
  string line;
  size_t* linexfile;
  long nt;

  sscanf(argv[2],"%ld", &NTOT);
  sscanf(argv[3],"%ld", &np);

  Npart = NTOT / np; // particelle per processore


  linexfile = new size_t[np];
  for(size_t i = 0; i<np; i++) linexfile[i]=0;

  ifstream *filestream = new ifstream[np];
  ifstream logstream(argv[4]);

  logstream >> global;
  logstream.close();

  global.log.open("paperino.log");

  Nline = 0;
  time_max = 0.;

  long startfile;

  if(argv[6])
    sscanf(argv[6],"%ld",&startfile);
  else startfile = 0;

  startfile = 0;

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
	  Nline++;
	  linexfile[i]++;
	  if(time > time_max) time_max = time;
	}
      filestream[i].close();
    }
   
  cout << "# nei file ci sono Nline = " << Nline << endl;

  NTIME = time_max/global.snap + 1; // massimo valore di NTIME

  cout << "# NTOT = " << NTOT << endl;
  cout << "# NTIME = " << NTIME << endl;

  CoupleLagrangianTracer** P = new CoupleLagrangianTracer*[NTIME];



  if(!P)
    {
      cout << "Memoria insufficiente " << endl;
      cout.flush();
      return 1;
    }

  double *Time = new double[NTIME];

  cout << "Allocando...\n\n";
  cout.flush();
 

  for(size_t i=0; i<NTIME; i++)
    {
      //      cout << i << " " << NTIME << endl;
      try{
	P[i] = new CoupleLagrangianTracer[NTOT];
      }
      catch (std::bad_alloc& ba)
	{
	  cout << "bad_alloc caught: " << ba.what() << '\n';
	  cout.flush();
	  throw int(123);
	}
      
      Time[i] = -1; // valore non inizializzato
    }

  for(size_t i=0; i<NTIME; i++)
    {
      for(size_t j=0; j<NTOT; j++)
	{ 
	  (*(P[i] + j)).part1.SetDead();
	  (*(P[i] + j)).part2.SetDead();
	}
    }
 
  cout << "Allocate le traiettorie con successo\n\n";
  cout.flush();

  nt = k = 0; time = time_old = 0.;
  size_t i_part, count;

  for(size_t i=0; i<np; i++) 
    {
      count = 0;
      if (i + startfile < 100)
	sprintf(nomefile,"%s/trajectories_LongLatZ_%.2ld.dat",argv[1],i + startfile);
      else
	sprintf(nomefile,"%s/trajectories_LongLatZ_%.3ld.dat",argv[1],i + startfile);
      filestream[i].open(nomefile);
      cout << "# Riapro  " << nomefile << endl;

      while(filestream[i] >> time >> part >> dummy >> age && count < linexfile[i] - 1) 
      {	
	nt = (long)(time / global.snap);
	i_part = part + i * Npart;
	Time[nt] = time;
	count ++;

	//	cout << nt << " " << time << " " << age << endl;
	
	if(i_part >= NTOT)
	  {
	    cout << "errore nell'indice di particella " << i_part << " " << part << endl;
	    return 0;
	  }

	if(nt >= NTIME)
	  {
	    cout << "errore nell'indice di tempo " << nt << " " << NTIME << endl;
	    return 0;
	  }

	*(P[nt]+i_part) = dummy;
	(*(P[nt]+i_part)).Age(age * 86400.);
	(*(P[nt]+i_part)).part1.SetAlive();
	(*(P[nt]+i_part)).part2.SetAlive();
      }
      filestream[i].close();
    }

  cout << "# lettura terminata\n " << endl;

  cout << "# Faccio il movie....\n " << endl;

  Movie(P, Time, NTIME, 6, NTOT, argv[1]);

  return 0;
} 

