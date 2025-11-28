
/*
Pensata per includere mortalita' delle tartarughe causata da basse temperature. 
Include DensityT.C con tre diverse funzioni che calcolano la mortalita'.  
*/

#include <iostream>
#include <stdio.h>
#include <string.h>
#include <cstdlib>
#include <fstream>

#include "FieldAbstract.C"
#include "VectorField.C"

//Global global;
#define R_earth 6378137 // equatoriale, in metri

//Altri include necessari
#include "randoknuth.c"
#include "KinematicModel.C"
#include "LagrangianTracer.C"
#include "CoupleLagrangianTracer.C"

//#include "DensityT.C" // Tutte le funzioni che considerano la Temperatura
//#include "DensityCHL.C" // Tutte le funzioni che considerano la Clorofilla


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


//***********************************************************************
//***********************************************************************

int main(int argc, char *argv[])
{
  if(!argv[1])
    {
      cout << "\n\n# PostProcesser Anchovies Chlorofilla - TESTING\n\n " << endl;
      printf("\nUsage %s <nome directory simulazione> <numero particelle> <numero processori> <file parametri> Ora Legge Fin Qui \n\n<elencoNetCDF_CHL>\n\n", argv[0]);
      return 0;
    }

  size_t part, Nline, n, Npart, np, NTOT, k, NTIME, nt;
  double time_max, time_old, age, time;
  char nomefile[500];
  CoupleLagrangianTracer dummy;
  string line;
  size_t* linexfile;

  sscanf(argv[2],"%ld", &NTOT);
  sscanf(argv[3],"%ld", &np);

  Npart = NTOT / np; // particelle per processore
  size_t loadedFiles; 
  size_t maxLoadedFiles = 16; // Per non esaurire RAM - deve essere divisore di 128
  if (np < 16) 
    maxLoadedFiles = np;

  cout << "\n\n# Carico " << maxLoadedFiles << " file con le traiettorie;\n\n";

  linexfile = new size_t[np];
  for(size_t i = 0; i<np; i++) linexfile[i]=0;

  ifstream *filestream = new ifstream[np];
  ifstream logstream(argv[4]);

  logstream >> global;
  logstream.close();

  //global.log.open("donaldDuck.log");

  Nline = 0;
  time_max = 0.;

  // Qui legge tutti i file ma non memorizza - OK
  for(size_t i=0; i<np; i++) 
    {
      if (i<100)
	sprintf(nomefile,"%s/trajectories_LongLatZ_%.2ld.dat",argv[1],i); // Primo problema - risolto
      else
	sprintf(nomefile,"%s/trajectories_LongLatZ_%.3ld.dat",argv[1],i);
      cout << "# apro il file " << nomefile << endl;
      filestream[i].open(nomefile);

      while(filestream[i] >> time >> part >> dummy >> age)
	{ 
	  Nline++;
	  linexfile[i]++;
	  if(time > time_max) time_max = time;
	}
      filestream[i].close();
    }
   
  cout << "# nei file ci sono Nline = " << Nline << endl;

  // Inizializza P
  NTIME = time_max/global.snap + 1; // massimo valore di NTIME

  cout << "# NTOT = " << NTOT << endl;
  cout << "# NTIME = " << NTIME << endl;

  CoupleLagrangianTracer** P = new CoupleLagrangianTracer*[NTIME];

  if(!P) cout << "Memoria insufficiente " << endl;

  double *Time = new double[NTIME];
 
  ////Attento, non NTOT ma al massimo Npart*maxLoadedFiles
  size_t actualNTOT; 
  if (NTOT > Npart*maxLoadedFiles)
    actualNTOT = Npart*maxLoadedFiles;
  else
    actualNTOT = NTOT;

  for(size_t i=0; i<NTIME; i++)
    {
      P[i] = new CoupleLagrangianTracer[actualNTOT];
      if(!P[i]) {
	  cout << "Memoria insufficiente " << i << endl;
	  return 0;
	}
    }

  for(size_t i=0; i<NTIME; i++)
    {
      for(size_t j=0; j<actualNTOT; j++)
	{ 
	  (*(P[i] + j)).part1.SetDead();
	  (*(P[i] + j)).part2.SetDead();
	}
    }
  //Tutto pronto

  nt = k = 0; time = time_old = 0.;
  size_t i_part, count;

  loadedFiles = 0;

  ofstream outFile;

  //cout << "\nClorofilliamo tutti insieme!\n" << endl << endl;

  // Dichiaro e inizializzo i campi 2D che conterrano le info sulla chl  
  class  ScalarField3D* CHL1 = Switcher(OceanColor_MyOcean2D);
  class  ScalarField3D* CHL2 = Switcher(OceanColor_MyOcean2D);
  ListOfFile ListaCHL;
  ListaCHL.Init(argv[5]);
  ListaCHL.LoadFieldParameters();

  CHL1->Init(ListaCHL.Nome(0),"lon","lat","depth","NULL");
  CHL2->Init(ListaCHL.Nome(0),"lon","lat","depth","NULL");

  size_t day, dayPrev;
  double kappa = 0.9, intChl_11, intChl_12, intChl_21, intChl_22, age1=0, age2=0.;
  double chl1, chl2, CHLinterp, f1;
  double x,y; //coordinate

  // Processo tutti i files a blocchi di 'maxloadedFiles'
  while (loadedFiles < np)
    {

      for(size_t j=0; j<NTIME; j++)
	{
	  for(size_t k=0; k<actualNTOT; k++)
	    { 
	      (*(P[j] + k)).part1.SetDead();
	      (*(P[j] + k)).part2.SetDead();
	    }
	}
      
      cout << "Inizializzato i tracers\n";
      for(size_t i=loadedFiles; i<maxLoadedFiles; i++) // sostituisce for(i=0;i<np;i++)
	{
	  count = 0;
	  if (i<100)
	    sprintf(nomefile,"%s/trajectories_LongLatZ_%.2ld.dat",argv[1],i); 
	  else
	    sprintf(nomefile,"%s/trajectories_LongLatZ_%.3ld.dat",argv[1],i);
	  
	  filestream[i].open(nomefile);
	  
	  cout << "# Riapro  " << nomefile << endl;
	  
	  while(filestream[i] >> time >> part >> dummy >> age && count < linexfile[i] - 1) 
	    {	
	      nt = (size_t)(time / global.snap);
	      i_part = part + i * Npart;
	      Time[nt] = time;
	      count ++;
	      
	      if(i_part >= actualNTOT)
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

	  // Ho caricato un file con Npart coppie
	  cout << "\n# Azioni...Clorofilla...\n\n";
	  //File Output
	  sprintf(nomefile,"/home/palatella/Corrado/TestingCHL/fileOutput_%.3ld.dat",i);	  
	  outFile.open(nomefile);
	  outFile << "# Sono il file " << i << endl;      	  
	  /*
	  day = floor(Time[NTIME-1]/86400);
	  for (nt = NTIME -1; nt > 0; nt++)
	    {
	      dayPrev = day -1;	      
	      //cout << nt/1./NTIME << "\r"; cout.flush(); // Lo tengo?
	      
	      if(day < ListaCHL.nfile)
		{
		  if(day != dayPrev)
		    {
		      CHL1->Load(ListaCHL.Nome(day),"CHL",0,"NULL");
		      CHL2->Load(ListaCHL.Nome(dayPrev),"CHL",0,"NULL");
		      day = dayPrev;
		      dayPrev = floor(Time[nt]/86400);
		    }
		  
		  for(size_t k = 0; k < actualNTOT; k++) //Per ogni coppia
		    {  
		      f1 = (Time[nt] - day*86400.) / 86400.;
		      intChl_11 = intChl_12 = intChl_21 = intChl_22 = 0.;
		      
		      if( (*(P[nt]+k)).part1.Alive())
			{
			  x=(*(P[nt]+k)).part1.Long();
			  y=(*(P[nt]+k)).part1.Lat();
			  chl1 = CHL1->Value(x,y,-1.); //
			  chl2 = CHL2->Value(x,y,-1.);
			  CHLinterp = chl1 * (1. - f1) + chl2 * f1;
			  intChl_11 += CHLinterp / kappa;
			  intChl_12 += CHLinterp / (CHLinterp + kappa);
			  age1=(*(P[nt]+k)).part1.Age();
			}
		      else 
			outFile << x << "  " << y << "  " << age1 << "  " << intChl_11 << "  " << intChl_12 << endl;
		      
		      if( (*(P[nt]+k)).part2.Alive())
			{
			  x=(*(P[nt]+k)).part2.Long();
			  y=(*(P[nt]+k)).part2.Lat();
			  chl1 = CHL1->Value(x,y,-1.); //
			  chl2 = CHL2->Value(x,y,-1.);
			  CHLinterp = chl1 * (1. - f1) + chl2 * f1;
			  intChl_21 += CHLinterp / kappa;
			  intChl_22 += CHLinterp / (CHLinterp + kappa);
			  age2=(*(P[nt]+k)).part2.Age();
			}
		      else 
			outFile << x << "  " << y << "  " << age2 << "  " << intChl_21 << "  " << intChl_22 << endl;

		    }
		}
	    }

	  // Fine clorofilla
	  */
	  outFile.close();
	  loadedFiles++;
	}
      
      cout << "# caricamento da " << loadedFiles-maxLoadedFiles << " a " << loadedFiles-maxLoadedFiles+1 << " terminato\n " << endl;
      
    } // end while block

  return 0;

} 

