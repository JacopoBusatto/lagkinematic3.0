#include <iostream>
#include <ctime>
#include <string>

#include <sys/types.h>
#include <sys/stat.h>


#include <mpi.h>

#define R_earth 6378137 // equatoriale

#include "randoknuth.c"

#include "FieldAbstract.C"
#include "VectorField.C"

#include "KinematicModel.C"
#include "LagrangianTracer.C"
#include "CoupleLagrangianTracer.C"
#include <fstream>



#define AMIROOT (!me)

//using namespace std;



int main(int argc, char *argv[])
{

  if(!argv[1])
       {
      cout << "\n\n# L A G K I N E M A T I C 1.0\n\n# LAGrangian with KINematic deterministic subgrid\n\n# by L. Palatella, G. Lacorata, R. Corrado, F. Grasso\n# 31 Jan 2015\n\n" << endl;
      printf("\nUsage %s file_elenco_netcdf file parametri condizioniiniziali filelog\n\nLa prima riga di file elenco deve contenere il nome del file MFS/MyOcean/POM con la descrizione dei file\n\n", argv[0]);
      return 0;
    }

  Vec3<double> Ui, X, tmp;
  Doub x,y,z, M;
  long int SEED = 123435, npa;
  size_t j;
  char nomelog[500], nomeout[500], nomep[500], nomeer[500],nomef[500], nomea[500];
  int me, np, i;
  
  MPI::Init(argc,argv);
  np = MPI::COMM_WORLD.Get_size();
  me = MPI::COMM_WORLD.Get_rank();

  ifstream parametri, condizioni_iniziali, anchovy_file, ele;
  char *stringa1, *stringa2, dumm[500], dname[500];
  size_t ndname, i0, i1;

  if(AMIROOT)
    {
      ele.open(argv[1]);
      ele >> dumm; ele >> dumm;
      ele.close();

      stringa1 = CopiaDaSlashaPunto(dumm);
      stringa2 = CopiaDaSlashaPunto(argv[3]);

      //ndname = sprintf(dname,"/store3/JacopoBusatto/LAGKINEMATIC_OUTS/lagkinematic_%s",argv[4]);
      ndname = sprintf(dname,"%s",argv[4]);


      // ndname = sprintf(dname,"/store3/JacopoBusatto/LAGKINEMATIC_OUTS/lagkinematic_%s_%s_%s",argv[4],stringa1,stringa2);
      // ndname = sprintf(dname,"/cantiere/palatella/lagkinematic_%s_%s_%s",argv[4],stringa1,stringa2);
      //ndname = sprintf(dname,"/Users/luigipalatella/Work/MareDiRoss/OUT/lagkinematic_%s_%s_%s",argv[4],stringa1,stringa2);
      fprintf(stdout,"Output nella directory : %s\n",dname);
    
      mkdir(dname,S_IWUSR | S_IXUSR | S_IRUSR);

      sprintf(nomep,"%s/parametri.dat",dname);
      sprintf(nomef,"%s/tipofile.dat",dname);

      global.param.open(nomep);
      global.tfile.open(nomef);
    }

  MPI::COMM_WORLD.Bcast(dname,500,MPI_CHAR,0);

  sprintf(nomelog,"%s/log%.2d.dat",dname, me);
  sprintf(nomeout,"%s/trajectories_LongLatZ_%.2d.dat",dname, me);
  sprintf(nomeer,"%s/errori_%.2d.dat",dname, me);
  sprintf(nomea,"%s/absorption_%.2d.dat",dname, me);

  global.er.open(nomeer);

  fprintf(stdout,"%s %s\n", nomelog, nomeout);

  ofstream OUT(nomeout);
  global.log.open(nomelog); // apro il fil di log
  global.absorptionFile.open(nomea);


  global.log << "\n---------------------------------------------------------" << endl;
  global.log << "\n L A G K I N E M A T I C 1.0\n\n LAGrangian with KINematic deterministic subgrid\n\n by L. Palatella, G. Lacorata, R. Corrado, f. Grasso\n 31 Jan 2015\n\n" << endl;
  
  global.log << "\n command line: ";
  for(j=0; j<argc; j++ ) global.log << argv[j] << " ";
  global.log << "\n";
  
  global.log << "\n---------------------------------------------------------" << endl << endl;

  Vec3<double> a1(1,35,-100.), a2(5,40,-1.);

  for(i=0;i<np;i++)
    {
      if(i==me) 
  {
    condizioni_iniziali.open(argv[3]);
    
    if(!condizioni_iniziali)
      {
        global.GestoreErrori("File condizioni iniziali non trovato");
        cout << "File condizioni iniziali non trovato" << endl;

        MPI::Finalize();
        return 0;
      }
    condizioni_iniziali >> npa;
    condizioni_iniziali.close();

    parametri.open(argv[2]);

    if(!parametri)
      {
        global.GestoreErrori("File parametri non trovato");
        MPI::Finalize();
        return 0;
      }

    fprintf(stdout,"Proc. %d legge i files NetCDF per impostare i Fields\n\n",me); 
    fflush(stdout);
    parametri >> global; // leggo i parametri dinamici della simulazione   

    Lista.Init(argv[1]); // leggo la lista dei file

    cout << "Carico i parametri dei file " << endl;
    cout.flush();

    Lista.LoadFieldParameters(); 

    parametri.close();

    cout.flush();
  }
      MPI_Barrier(MPI_COMM_WORLD);
    }
  
  global.log << "Inizializzo i campi vettoriali U e Unew" << endl << endl;
  cout << "Inizializzo i campi vettoriali U e Unew" << endl << endl;
  cout.flush();
  global.log.flush();

  /************************ U E Unew ********************************/

  //  ScalarField3D Mask, GPM, GPA;
  VectorField3D U(Lista.KindOfFile), Unew(Lista.KindOfFile); /* serve solo a dare i puntatori del tipo giusto*/

  cout << "Init() dei  campi vettoriali " << endl;
  global.log << "Init() dei  campi vettoriali " << endl;

  U.Init(Lista); Unew.Init(Lista);

  /************************ MASK ********************************/

  cout << "Init & load mask and orography " << endl;
  global.log << "Init & load mask and orography " << endl;

  ScalarField3D_Mask* Mask =  SwitcherMask(Lista.KindOfFile);

  if(Lista.IsThereW)
    {
      Mask->Init(Lista.Nome(0), Lista.nomelonW, Lista.nomelatW, Lista.nomezW);
      Mask->Mask(Lista.Nome(0), Lista.nomeW);
    }

  else 
    {
      Mask->Init(Lista.Nome(0),Lista.nomelonV, Lista.nomelatU, Lista.nomezU);
      Mask->Mask(Lista.Nome(0), Lista.nomeU);

      cout << "qui e' " << Mask->NTimes(Lista.Nome(0)) << endl;
    }  

  if(global.sigma || strcmp(Lista.nomeorografia,"NULL"))
    Mask->LoadOrografia(Lista.Nome(0), Lista.nomeorografia);


  /************************ PB e PBH ********************************/

  ScalarField3D* GPM =  Switcher(Lista.KindOfFile);
  ScalarField3D* GPA =  Switcher(Lista.KindOfFile);


  if(Lista.sigma)
    {
      cout << "Init PB e PBH " << endl;
            
      GPM->Init(Lista.Nome(0),Lista.nomelonW, Lista.nomelatW, Lista.nomezW, "NULL");
      GPA->Init(Lista.Nome(0),Lista.nomelonW, Lista.nomelatW, Lista.nomezW, "NULL");
    }

  long int nfields =  Mask->NTimes(Lista.Nome(0));
  
  cout << "numero tempi per file: " << nfields << endl;

  if(nfields != Lista.nfieldsxfile) Lista.nfieldsxfile = nfields;

  /*************** Fine inizializzazioni *******************/

  
  i0 = global.ifile_start;
 
  if(nfields==1)
    i1 = i0 + 1;
  else i1 = i0;

  U.X->Load(Lista.Nome(i0),Lista.nomeU, 0,Lista.nomeextensionU);
  U.Y->Load(Lista.Nome(i0),Lista.nomeV, 0,Lista.nomeextensionV);
  


  if(Lista.IsThereW && strcmp(Lista.nomeW,"vo"))
    U.Z->Load(Lista.Nome(i0),Lista.nomeW, 0, Lista.nomeextensionW);


  /* se uso di nuovo vo vuol dire che usa close divergence*/
  if(!strcmp(Lista.nomeW,"vo")) 
    {
      cout << "chiudo con la divergenza\n";
      U.Z->CloseDivergence(*(U.X),*(U.Y));
    }

  
  if(nfields == 1)
    {
  
      if(Lista.sigma)
  {
    GPM->Load(Lista.Nome(i1), "PHB",0);
    GPA->Load(Lista.Nome(i1), "PH",0);
  }

      Unew.X->Load(Lista.Nome(i1),Lista.nomeU,0,Lista.nomeextensionU);
      Unew.Y->Load(Lista.Nome(i1),Lista.nomeV,0, Lista.nomeextensionV);

      if(Lista.IsThereW)
  U.Z->Load(Lista.Nome(i1),Lista.nomeW, 0 , Lista.nomeextensionW);
    }
  else
    {
      if(Lista.sigma)
  {
    GPM->Load(Lista.Nome(i0), "PHB",1);
    GPA->Load(Lista.Nome(i0), "PH",1);
  }

      Unew.X->Load(Lista.Nome(i0),Lista.nomeU, 1, Lista.nomeextensionU);
      Unew.Y->Load(Lista.Nome(i0),Lista.nomeV, 1, Lista.nomeextensionV);

      if(Lista.IsThereW)
  U.Z->Load(Lista.Nome(i0),Lista.nomeW,1, Lista.nomeextensionW);
    }


  ranf_start(SEED);

  U.Time(0.);
  Unew.Time(Lista.DT);

  global.log <<"# " << U.Time() << " " << Unew.Time() << endl;

  KinematicModel3D sgs3D(global.Lstart3D, global.Lend3D, global.Fact3D, global.Eps3D, global.Eta);
  KinematicModel2D sgs2D(global.Lstart2D, global.Lend2D, global.Fact2D, global.Enstrophy);


  if(npa%np) 
    {
      cout << "numero di particelle non-multiplo del numero di processori\n\n";
      global.er << "numero di particelle non-multiplo del numero di processori\n\nSto uscendo\n\n";
      MPI::Finalize();
      return 0;
    }

  global.nparticles = npa / np; 

  fprintf(stdout,"\nnumero coppie = %ld\n", npa);

  global.log << "\n\nnumero coppie totali " << npa << " divise su " << np << " processori" << endl << endl;

  global.log << global << endl;

  global.log << "Parametri NetCDF letti dal file " << Lista.tipo_files << endl;
  global.log << Lista;

  if(AMIROOT)
    {
      global.param << global;
      global.tfile << Lista;
    }

  global.log << sgs3D;
  global.log << sgs2D;

  CoupleLagrangianTracer *particles, dummy;
  particles = new CoupleLagrangianTracer[global.nparticles];
  double dtstart;

  global.log << "\nLeggo condizioni iniziali dal file " << argv[3] << endl << endl; 
  cout << me <<": Leggo condizioni iniziali dal file " << argv[3] << endl; 
  
  if(global.plume)
    {
      dtstart = global.time_end_plume / global.nparticles;    
      global.log << "Plume parameters: le particelle partono ogni " << dtstart << " sec"<< endl << endl;
    }
  
  for(i=0;i<np;i++)
    {
      if(i==me) 
  {    
    condizioni_iniziali.open(argv[3]);
    condizioni_iniziali >> npa;

    for(j= 0; j< global.nparticles * me; j++) condizioni_iniziali >> dummy;
    for(j= me * global.nparticles ; j< global.nparticles * (me + 1); j++)
      { 
        condizioni_iniziali >> particles[j%global.nparticles];     
        //global.log << " carico particella " << j << " " << particles[j%global.nparticles] << endl;
      }
    condizioni_iniziali.close();
    
    if(global.buoyant)
      for(j= 0; j< global.nparticles; j++) particles[j].SetBuoyant();

    if(global.plume)
      {
        for(j=0; j< global.nparticles; j++) 
    {
      particles[j].part1.Start(dtstart * j);
      particles[j].part2.Start(dtstart * j);
    }
      }

  }
      //      MPI_Barrier(MPI_COMM_WORLD);     LP + JB 17 feb 2022 
    }

  if(global.anchovy) 
    {
      global.log << "\nParametri dal file ANCHOVY\n\n";
      global.log << anchovy;
    }

  /***********  INIZIALIZZA I VICINI *************/

  global.log << "Inizializzo i vicini....\n\n";

  bool isUsable;

  if(Lista.irregular)
    {
      
      U.X->Points4Neighbour(isUsable);
      
      if(isUsable)
  {
    for(j=0; j<global.nparticles; j++)
      {
        tmp = particles[j].part1.LongLatZ();
        
        particles[j].part1.vicinix.Init( U.X->Points4Neighbour(isUsable), tmp);
        particles[j].part1.viciniy.Init( U.Y->Points4Neighbour(isUsable), tmp);
        particles[j].part1.viciniz.Init( U.Z->Points4Neighbour(isUsable), tmp);
        
        tmp = particles[j].part2.LongLatZ();
        
        particles[j].part2.vicinix.Init( U.X->Points4Neighbour(isUsable), tmp);
        particles[j].part2.viciniy.Init( U.Y->Points4Neighbour(isUsable), tmp);
        particles[j].part2.viciniz.Init( U.Z->Points4Neighbour(isUsable), tmp);
      }
  }
    }


  /* Si controlla se le particelle partono dentro il dominio altrimenti si ammazzano */

  for(i=0; i<global.nparticles; i++)
    {
      particles[i].part1.ApplyBoundaryCondition(Mask, GPM, GPA);      
      particles[i].part2.ApplyBoundaryCondition(Mask, GPM, GPA);      
    }

  /* INIZIO LOOP SUL TEMPO */

  MPI_Barrier(MPI_COMM_WORLD);      

  global.log << "\nINIZIO SIMULAZIONE.... \n\n " << endl; 

  size_t nkf, nin, nin2;
  double inner_time, time1, time2;

  /* Si parte! */

  nkf = (Lista.nfieldsxfile == 1 ? global.ifile_start + 2 : global.ifile_start);
  nin = (Lista.nfieldsxfile == 1 ? 0 : 2);

  for(; nkf < Lista.nfile; nkf++)
    {
      for(; nin < Lista.nfieldsxfile; nin++)  
  {
    for(inner_time=0.; inner_time < Lista.DT; inner_time += global.dt())
      {        
        for(i=0; i<global.nparticles; i++)
    {
      /* UNCOMMENTED JB 20240921 */
      /* COMMENTED JACOPO BUSATTO 20240514
      ADDED BY JACOPO BUSATTO 20231213 */
      if(global.Snapshot()) 
        {
          if(particles[i].ToBePrinted())
      OUT << global.Time() << " " << i <<" " << particles[i];
        }
      /* ADDED BY JACOPO BUSATTO 20231213 */
      /* */
      particles[i].EvolveCampoRisolto(U, Unew, Mask, GPM, GPA);
      /* Evolve campo risolto applica le condizioni al contorno */

      if(global.Kin3D) particles[i].EvolveAnalytical(sgs3D, Mask, GPM, GPA);
      if(global.Kin2D) particles[i].EvolveAnalytical(sgs2D, Mask, GPM, GPA);


      /* anche gli evolutori cinematici dopo applicano le condizioni al contorno*/

      /* Ogni evolutore applica al suo interno le corrette boundary conditions*/
      
      /* COMMENTED JB 20240921 
      /* UNCOMMENTED BY JACOPO BUSATTO 20231213 
      /* COMMENTED BY JACOPO BUSATTO 20231213   
      if(global.Snapshot()) 
        {
          if(particles[i].ToBePrinted())
      OUT << global.Time() << " " << i <<" " << particles[i];
        }
      /* COMMENTED BY JACOPO BUSATTO 20231213 
      /* UNCOMMENTED BY JACOPO BUSATTO 20231213 
      */
    }
        global.IncrTime(global.dt());    
        } 
      
    U ^ Unew; // li scambio senza copiarli!!
    
    time2 = Lista.DT * ( nin + (nkf - global.ifile_start) * Lista.nfieldsxfile );
    time1 = time2 - Lista.DT;

    U.Time(time1);
    Unew.Time(time2);

    try
      {
        Unew.X->Load(Lista.Nome(nkf), Lista.nomeU, nin, Lista.nomeextensionU);
        Unew.Y->Load(Lista.Nome(nkf), Lista.nomeV, nin, Lista.nomeextensionV);
        if(Lista.IsThereW)
    Unew.Z->Load(Lista.Nome(nkf), Lista.nomeW, nin, Lista.nomeextensionW);
       
        if(Lista.sigma)
    {
      GPM->Load(Lista.Nome(nkf), "PHB",nin);
      GPA->Load(Lista.Nome(nkf), "PH",nin);
    }
      }
    catch(int ercode)
      {
        global.er << "Lanciata nal main loop sul time " << ercode << endl;
        global.er << "file = " << Lista.Nome(nkf) << " nin = " << nin << endl; 
        cout << "Lanciata nal main loop sul time " << ercode << endl;
        MPI::Finalize();
        return 0;
      }
    global.log <<"# " << U.Time() << " " << Unew.Time() << endl;
  }
      nin = 0;
    }

  global.er.close();
  global.log.close();
  global.absorptionFile.close();

  OUT.close();

  if(AMIROOT)
    {

      global.param.close();
      global.tfile.close();
    }

  MPI::Finalize();
}
