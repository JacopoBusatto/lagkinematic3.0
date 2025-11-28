//#include <iostream>
//#include <math.h>


using namespace std;

class Global
{
private:
  double time, time_step;
public:
  double kt3d;
  double kt2d, tt, snap;
  size_t nparticles, ifile_start, numero_errori;
  ofstream log, param, tfile, er, absorptionFile;

  bool absorptionTOP, absorptionBOTTOM, absorptionLATERAL;
  bool Kin3D, Kin2D, plume, backtraj, buoyant, anchovy, sigma, ocean;
  double Lstart3D, Lend3D, Fact3D, Eps3D, Eta;
  double Lstart2D, Lend2D, Fact2D, Enstrophy, time_end_plume;
  double Vdep;
  
  Global(){
    /* Set default value che rimangono in caso non vengano specificati nel file parametri*/
    numero_errori = 0;
    time = time_end_plume = 0.;
    time_step = 120.;
    kt3d = 0.1; kt2d = 0.1; tt = 0.02; snap = 3600.;
    nparticles = 0;
    ifile_start = 0; 
    Kin3D = Kin2D = plume = backtraj = buoyant = anchovy = false; 
    absorptionTOP = absorptionBOTTOM = absorptionLATERAL = false;
    // valori di default

    Lstart3D =  Lend3D =  Fact3D =  Eps3D =  Eta = 1.;
    Lstart2D =  Lend2D =  Fact2D =  Enstrophy = 1.;
  }

  size_t GestoreErrori(const char errore[])
  {
    //cout << "Errore: " << errore << " " << numero_errori << endl;
    // cout.flush();

    er << "Errore: " << errore << " " << numero_errori << endl;
    er.flush();

    numero_errori++;
    return numero_errori;
  }
  
  double Time(){ return time;}
  void Time(double t){time = t;}
  void IncrTime(double t){time += t;}
  double dt(){ return time_step;}
  void dt(double ddt){time_step = ddt;}
  //  bool Plume(){ return plume;}

  bool Snapshot(double sn)
  {
    size_t ttt, ts;
    /* JACOPO BUSATTO 20240418*/
    if (time == 0) return true;
    /* JACOPO BUSATTO 20240418*/
    ttt = (size_t)time; ts = (size_t)sn;
    if(!(ttt % ts) ) return true;
    else return false;
  }

  bool Snapshot()
  {
    return Snapshot(snap);
  }

};

Global global;

class ListOfFile
{
public:
  long int nfile, nfieldsxfile, ndimlon;
  char **nome_files;
  char tipo_files[20], nome[200];
  char nomelatU[20], nomelonU[20], nomezU[20];
  char nomelatV[20], nomelonV[20], nomezV[20];
  char nomelatW[20], nomelonW[20], nomezW[20];
  char nomeU[20], nomeV[20], nomeW[20], nomeOrografia[20];
  char nomeextensionU[20];
  char nomeextensionV[20];
  char nomeextensionW[20];
  char nomeorografia[20];
  KindOfFile_t KindOfFile; 
  double DT;
  bool ocean, IsThereW, sigma, irregular, orografia, extension, extensionW;

  /*
    ocean 0 => atmosfera
    sigma 1 => livelli sigma, necessaria orografia
    irregular => griglia irregolare, anche WRF e' irregolare in long e lat
   */


  ListOfFile()
  {
    nfile = 1;
    nfieldsxfile = 1; // default values
  }

  ListOfFile(const char elenco[])
  {
    Init(elenco);
  }
 
  inline void Init(const char elenco[])
  {
    FILE *in;
    size_t dum;

    /* default values*/

    ocean = true;

    DT = 86400.; IsThereW = false;

    in = fopen(elenco,"r");

    dum = fscanf(in,"%s",tipo_files);

    size_t i;
    i=0;

    while(fscanf(in,"%s\n",nome)!= EOF){i++;}

    nfile = i;
    fclose(in);
    fprintf(stderr,"Ci sono %ld files\n",nfile);

    nome_files = new char *[nfile];
   
    in = fopen(elenco,"r");
    dum = fscanf(in,"%s\n",tipo_files);


    for(i=0; i<nfile; i++)
      {
	nome_files[i] = new char[300];
      }


    for(i=0; i<nfile; i++)
      {
	if(!global.backtraj) dum = fscanf(in,"%s\n",nome_files[i]); 
	if(global.backtraj) dum = fscanf(in,"%s\n",nome_files[nfile - 1 - i]); // al rovescio!
      }
    fclose(in);
  }

  char *Nome(size_t k)
  {
    if(k<0 || k >= nfile) 
      {
	fprintf(stderr,"Indice file out of range\n");
	global.GestoreErrori("Indice file out of range");
	return nome_files[0];
      }
    return nome_files[k];
  }

  void LoadFieldParameters()
  {
    LoadFieldParameters(tipo_files);
  }


  void LoadFieldParameters(const char *nomefiletipo)
  {
    /*
      tipo_files e' il nome del file parametri 
      non e' l'intero che indica esattamente il tipo
     */

    ifstream p(nomefiletipo);
    char tipo[20];
    double vd;
    long int vi;
    int kind;
    bool bo;

    cout << "Setting File parameters....\n\n";

    while(!p.eof())
    {
      p >> tipo; cout << tipo << endl;
      if(!strcmp(tipo,"KINDOFFILE")) { p >> kind; 
	KindOfFile = (KindOfFile_t)kind;
      }
      if(!strcmp(tipo,"TIMESTEP_O")) { p >> vd; DT = vd;}
      if(!strcmp(tipo,"NUM_FIELDS_FILE")) { p >> vi; nfieldsxfile = vi;}
      if(!strcmp(tipo,"N_DIM_LONLAT")) { p >> vi;  ndimlon = vi; }
      if(!strcmp(tipo,"LONGITUDE_U")) { p >> nomelonU; }
      if(!strcmp(tipo,"LATITUDE_U")) { p >> nomelatU; }
      if(!strcmp(tipo,"DEPTH_U")) { p >> nomezU; }
      if(!strcmp(tipo,"LONGITUDE_V")) { p >> nomelonV; }
      if(!strcmp(tipo,"LATITUDE_V")) { p >> nomelatV; }
      if(!strcmp(tipo,"DEPTH_V")) { p >> nomezV; }
      if(!strcmp(tipo,"LONGITUDE_W")) { p >> nomelonW; }
      if(!strcmp(tipo,"LATITUDE_W")) { p >> nomelatW; }
      if(!strcmp(tipo,"DEPTH_W")) { p >> nomezW; }

      if(!strcmp(tipo,"EXTENSION_U")) { p >> nomeextensionU; }
      if(!strcmp(tipo,"EXTENSION_V")) { p >> nomeextensionV; }
      if(!strcmp(tipo,"EXTENSION_W")) { p >> nomeextensionW; }

      if(!strcmp(tipo,"SIGMA_LEVEL")) { p >> bo; sigma = bo;}
      if(!strcmp(tipo,"OROGRAPHY")) { p >> nomeorografia;}
      if(!strcmp(tipo,"OCEAN")) { p >> bo; ocean = bo;}
      if(!strcmp(tipo,"IRREGULAR")) { p >> bo; irregular = bo;}

      if(!strcmp(tipo,"U")) { p >> nomeU; }
      if(!strcmp(tipo,"V")) { p >> nomeV; }
      if(!strcmp(tipo,"W")) { p >> nomeW; }
    }


    if(!strcmp(nomeW,"NULL")) IsThereW = false;
    else IsThereW = true;

    if(!strcmp(nomeorografia,"NULL")) orografia = false;
    else orografia = true;

    if(!strcmp(nomeextensionU,"NULL")) extension = false;
    else extension = true;
    
    if(!strcmp(nomeextensionW,"NULL")) extensionW = false;
    else extensionW = true;
    
    
    p.close();
  }
};

ListOfFile Lista; 

std::istream& operator >> (std::istream& istream, Global& C)
{
  char tipo[20];
  double vd;
  long int vi;
  bool b;

  while(!istream.eof())
    {
      istream >> tipo;

      cout << " Setting " << tipo << endl; 
      //      if(!strcmp(tipo,"NPARTICLES")) { istream >> vi; C.nparticles = vi;}
      if(!strcmp(tipo,"TIMESTEP")) { istream >> vd; C.dt(vd);}
      if(!strcmp(tipo,"SOGLIA_CINEMATICO3D")) { istream >> vd; C.kt3d = vd;}
      if(!strcmp(tipo,"SOGLIA_CINEMATICO2D")) { istream >> vd; C.kt2d = vd;}
      if(!strcmp(tipo,"SOGLIA_ALIVE")) { istream >> vd; C.tt = vd;}
      if(!strcmp(tipo,"IFILE_START")) { istream >> vi; C.ifile_start = vi;}

      if(!strcmp(tipo,"LSTART3D")) { istream >> vd; C.Lstart3D = vd;}
      if(!strcmp(tipo,"LEND3D")) { istream >> vd; C.Lend3D = vd;}
      if(!strcmp(tipo,"FACT3D")) { istream >> vd; C.Fact3D = vd;}
      if(!strcmp(tipo,"Eps3D")) { istream >> vd; C.Eps3D = vd;}
      if(!strcmp(tipo,"VERTICAL_DECAY")) { istream >> vd; C.Eta = vd;}

      if(!strcmp(tipo,"BUOYANT")) { istream >> b; C.buoyant = b;}
      if(!strcmp(tipo,"ANCHOVY")) { istream >> b; C.anchovy = b;}
      if(!strcmp(tipo,"LSTART2D")) { istream >> vd; C.Lstart2D = vd;}
      if(!strcmp(tipo,"LEND2D")) { istream >> vd; C.Lend2D = vd;}
      if(!strcmp(tipo,"FACT2D")) { istream >> vd; C.Fact2D = vd;}
      if(!strcmp(tipo,"ENSTROPHY")) { istream >> vd; C.Enstrophy = vd;}
      
      if(!strcmp(tipo,"FLAG3D")) { istream >> b; C.Kin3D = b;}
      if(!strcmp(tipo,"FLAG2D")) { istream >> b; C.Kin2D = b;}
      if(!strcmp(tipo,"SNAPSHOT")) { istream >> vd; C.snap = vd;}
      if(!strcmp(tipo,"PLUME")) { istream >> b; C.plume = b;}
      if(!strcmp(tipo,"TIME_END_PLUME")) { istream >> vd; C.time_end_plume = vd;}
      if(!strcmp(tipo,"BACKTRAJECTORY")) { istream >> b; C.backtraj = b;}

      if(!strcmp(tipo,"ABSORPTIONTOP")) { istream >> b; C.absorptionTOP = b;}
      if(!strcmp(tipo,"ABSORPTIONBOTTOM")) { istream >> b; C.absorptionBOTTOM = b;}
      if(!strcmp(tipo,"ABSORPTIONLATERAL")) { istream >> b; C.absorptionLATERAL = b;}


    }

  return istream;
}

std::ostream& operator << (std::ostream& ostream, Global& C)
{
  ostream << "\nNPARTICLES " << C.nparticles << endl;
  ostream << "TIMESTEP " << C.dt() << endl;
  ostream << "SOGLIA_CINEMATICO3D " << C.kt3d << endl;
  ostream << "SOGLIA_CINEMATICO2D " << C.kt2d << endl;
  ostream << "SOGLIA_ALIVE " << C.tt << endl;
  ostream << "SNAPSHOT " << C.snap << endl << endl;
  ostream << "IFILE_START " << C.ifile_start << endl;
  
  ostream << "BUOYANT " << C.buoyant << endl;
  ostream << "ANCHOVY " << C.anchovy << endl;
  ostream << "PLUME " << C.plume << endl;
  if(C.plume) ostream << "TIME_END_PLUME " << C.time_end_plume << endl;

  ostream << "BACKTRAJECTORY " << C.backtraj << endl;

  ostream << "ABSORPTIONTOP " << C.absorptionTOP << endl;
  ostream << "ABSORPTIONBOTTOM " << C.absorptionBOTTOM << endl;
  ostream << "ABSORPTIONLATERAL " << C.absorptionLATERAL << endl;


  ostream << "FLAG3D " << C.Kin3D << endl;

  if(C.Kin3D)
    {
      ostream << "LSTART3D " << C.Lstart3D << endl;
      ostream << "LEND3D " << C.Lend3D << endl;
      ostream << "FACT3D " << C.Fact3D << endl;
      ostream << "Eps3D " << C.Eps3D << endl;
      ostream << "VERTICAL_DECAY " << C.Eta << endl << endl;
    }

  ostream << "FLAG2D " << C.Kin2D << endl;
  
  if(C.Kin2D)
    {
      ostream << "LSTART2D " << C.Lstart2D << endl;
      ostream << "LEND2D " << C.Lend2D << endl;
      ostream << "FACT2D " << C.Fact2D << endl;
      ostream << "ENSTROPHY " << C.Enstrophy << endl;
    }

  return ostream;
}

std::ostream& operator << (std::ostream& ostream, ListOfFile& C)
{

  ostream << "\nKINDOFFILE " << C.KindOfFile << endl;
  ostream << "\nTIMESTEP_O " << C.DT << endl;
  ostream << "NUM_FIELDS_FILE " << C.nfieldsxfile << endl;
  ostream << "LONGITUDE_U " << C.nomelonU << endl;
  ostream << "LATITUDE_U " << C.nomelatU << endl;
  ostream << "DEPTH_U " << C.nomezU << endl;
  ostream << "LONGITUDE_V " << C.nomelonV << endl;
  ostream << "LATITUDE_V " << C.nomelatV << endl;
  ostream << "DEPTH_V " << C.nomezV << endl;
  ostream << "LONGITUDE_W " << C.nomelonW << endl;
  ostream << "LATITUDE_W " << C.nomelatW << endl;
  ostream << "DEPTH_W " << C.nomezW << endl;
  ostream << "U " << C.nomeU << endl;
  ostream << "V " << C.nomeV << endl;
  ostream << "W " << C.nomeW << endl;
  ostream << "EXTENSION_U " << C.nomeextensionU << endl;
  ostream << "EXTENSION_V " << C.nomeextensionV << endl;
  ostream << "EXTENSION_W " << C.nomeextensionW << endl;
  ostream << "SIGMA_LEVEL " << C.sigma << endl;
  ostream << "OROGRAPHY " << C.nomeorografia << endl;
  ostream << "OCEAN " << C.ocean << endl;
  ostream << "IRREGULAR " << C.irregular << endl;
  ostream << "nfile  " << C.nfile << endl;

  return ostream;
}


class Anchovy
{
public:
  double timeeggs, time_up, time_down, zup, zdown, gam;
};

Anchovy anchovy; // globale come global, sono parametri fissi!

std::istream& operator >> (std::istream& istream, Anchovy& A)
{
  char tipo[20];
  double vd;
 
  while(!istream.eof())
    {
      istream >> tipo;

      cout << " Setting " << tipo << endl; 
  
      if(!strcmp(tipo,"TIMEEGGS")) { istream >> vd; A.timeeggs = vd;}
      if(!strcmp(tipo,"TIMEUP")) { istream >> vd; A.time_up = vd;}
      if(!strcmp(tipo,"TIMEDOWN")) { istream >> vd; A.time_down = vd;}
      if(!strcmp(tipo,"Z_UP")) { istream >> vd; A.zup = vd;}
      if(!strcmp(tipo,"Z_DOWN")) { istream >> vd; A.zdown = vd;}
      if(!strcmp(tipo,"GAMMA")) { istream >> vd; A.gam = vd;}
    }
  return istream;
}


std::ostream& operator << (std::ostream& ostream, Anchovy& A)
{
  ostream << "TIMEEGGS " << A.timeeggs << endl;
  ostream << "TIMEUP " << A.time_up << endl;
  ostream << "TIMEDOWN " << A.time_down << endl;
  ostream << "Z_UP " << A.zup << endl;
  ostream << "Z_DOWN " << A.zdown << endl;
  ostream << "GAMMA " << A.gam << endl << endl;

  return ostream;
}



char *CopiaDaSlashaPunto(char *src)
{
  size_t i,j,k,l,len;

  len = strlen(src);

  i = len;
  while(i>0 && src[i]!='/') i--;
  if(i!=0) i++;

  j=len;
  while(j>0 && src[j]!='.') j--; 

  printf("%ld %ld\n",i,j);

  if(j<i) j = len;

  k = j - i;

  char *out = (char *)malloc(sizeof(char)*(k+1));

  for(l=i;l<j;l++) out[l-i] = src[l];
  out[k]='\0';

  return out;
}


