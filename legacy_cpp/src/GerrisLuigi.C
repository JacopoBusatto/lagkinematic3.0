/** 
 * @file Gerris.C 
 * It manages the output of the simulator Gerris(R).
 * @author Fabio Grasso, Luigi Palatella, mia riedizione del 26 mar 2015 LP
 * @version 1.0
 */

/**
 * @class ScalarField3D_Gerris 
 * It loads a text file, column based, as e.g.: \n 
 * 1:x 2:y 3:z 4:P 5:Pmac 6:U .........    
 *  
 */
class ScalarField3D_Gerris : public ScalarField3D {

protected:
  static const size_t n_dim = 2;  /*!< Number of point coordinates  */
  MatDoub PtiGriglia;
  double x_min, y_min, z_min, x_max, y_max, z_max, time;
  VecDoub data;  
  vector<string> Header;
  
public:

  ScalarField3D_Gerris() : ScalarField3D(){}
  
   /** Memory allocation and loading of coordinates
   *  
   * @param netcdf_file output file by Gerris(R) or a bathymetry
   * @param nome_lon header label of longitude column
   * @param nome_lat header label of latitude column
   * @param nome_depth **not used, you pass NULL
   * @param nomeExtension **not used, you pass NULL
   * @return void 
   */

  void CloseDivergence(ScalarField3D& U, ScalarField3D& V){};

  inline void Init(const char netcdf_file[], const char nome_lon[],	\
		    const char nome_lat[], const char nome_depth[], const char nomeExtension[]) {
    string s,ss;
    string s1, s2, s3;
    int i=0, ii, ind_lon, ind_lat, primo_indice, secondo_indice;
    const char* cc_nf;
    ifstream f;
    std::stringstream stream(std::stringstream::in | std::stringstream::out);
   
    f.open(netcdf_file);

    cout << "Sto inizializzando campi Gerris\n";
  
    if(!f){
      cout << "Il file non esiste!\n" << endl;  //todo eccezione
      return ;
    }

    if (f.good()) getline(f,s);       
    if (s.length()>10){  // -- Acquisizione dell'header
      stream << s;
      stream >> ss;
      while (stream.good()){
	if(ss.compare("#")!=0)
	  Header.push_back(ss);
	stream >> ss;
      }
    }
    else{
      return; //todo eccezione! il file è vuoto
    }

    stream.str("");  // -- Azzero la parte rimanente dello stream
    stream.clear();  // -- visto che non è più 'good()', cancello il flag di errore
    
    n_punti=0; // -- n_punti e' della classe astratta ScalarField3D
    nz = 1; // -- nz e' della classe astratta ScalarField3D
    z.resize(1); // -- z e' della classe astratta ScalarField3D
    z[0] = 0.; // -- per convenzione profondità zero


    // -- Scorro tutto il file per determinare il numero di righe che lo compongono
    while(f.good()){  
      std:: getline(f,s);
      if(s.length()>10){
	n_punti++;
      }
    }
    f.close();

    ind_lon = -1; ind_lat = -1;
    // -- Scorro header e prendo gli indici di lon e lat in modo da renderli parametri effettivi
    for(i=0;i<Header.size();i++){
      if(Header[i].compare(nome_lon)==0) 
	ind_lon=i;
      if(Header[i].compare(nome_lat)==0)
	ind_lat=i;
    }
    
    if(ind_lon == -1 || ind_lat == -1){
      //todo throw eccezione
      cout << " Errore nei parametri ! " << endl;
      return;
    }
    
    if(ind_lon > ind_lat){
      primo_indice = ind_lat;
      secondo_indice = ind_lon;
    }
    else{
      primo_indice = ind_lon;
      secondo_indice = ind_lat;
    }

    // -- Adesso scorro di nuovo il file dall'inizio per caricare la matrice Ptigriglia 
    data.resize(n_punti);
    PtiGriglia.resize(n_punti,n_dim);    
    f.open(netcdf_file);
    if (f.good()) getline(f,s);   // -- Butto la prima linea, quella di intestazione
    
    i=0; ii=0;
    while(f.good() && i<n_punti){   
      std:: getline(f,s);
      if(s.length()>10)	{
	stream << s;
	stream >> ss;
	for(ii=0;ii<primo_indice;ii++){
	  stream >> ss;
	}
	if(ind_lon == ii){
	  PtiGriglia[i][0]=atof(ss.c_str());
	} else  if(ind_lat == ii){
	  PtiGriglia[i][1]=atof(ss.c_str());
	}
	for(;ii<secondo_indice;ii++){
	  stream >> ss;
	}
	if( ind_lon == ii){
	  PtiGriglia[i][0]=atof(ss.c_str());
	} else  if( ind_lat == ii){
	  PtiGriglia[i][1]=atof(ss.c_str());
	}
	stream.str(""); // -- Azzero la parte rimanente dello stream
	  
	//cout << "pto#: " << i << " crd: "<<PtiGriglia[i][0]<< " "<<PtiGriglia[i][1]<<" "<<PtiGriglia[i][2]<<endl;
	//cout << "U: " << U[i] << " V: " << V[i]   << " \n" << endl;
	
	// -- Aggiorno min e max
	if (i==0) {
	  x_min=PtiGriglia[i][0]; 
	  y_min=PtiGriglia[i][1]; 
	  if (n_dim==3)	  z_min=PtiGriglia[i][2]; 
	  x_max=PtiGriglia[i][0]; 
	  y_max=PtiGriglia[i][1]; 
	  if (n_dim==3) z_max=PtiGriglia[i][2]; 
	}
	else {
	  if(PtiGriglia[i][0]<x_min) x_min=PtiGriglia[i][0]; 
	  if(PtiGriglia[i][1]<y_min) y_min=PtiGriglia[i][1]; 
	  if (n_dim==3) if(PtiGriglia[i][2]<z_min) z_min=PtiGriglia[i][2]; 
	  if(PtiGriglia[i][0]>x_max) x_max=PtiGriglia[i][0]; 
	  if(PtiGriglia[i][1]>y_max) y_max=PtiGriglia[i][1]; 
	  if (n_dim==3) if(PtiGriglia[i][2]>z_max) z_max=PtiGriglia[i][2]; 
	}
	i++;
      }
    }
    f.close();
    
    // -- dbg
    /*for(i=0;i<Header.size();i++){
      cout << "header[" << i << "] = " << Header[i] << endl;
      }*/
    
    cout << "n_punti: " <<n_punti <<"\n" << endl;
    if (n_dim==3) cout << "x_min: " << x_min << " y_min: " << y_min << " z_min: " << z_min << endl;
    else cout << "x_min: " << x_min << " y_min: " << y_min << endl;
    if (n_dim==3)    cout << "x_max: " << x_max << " y_max: " << y_max << " z_max: " << z_max << endl;
    else cout << "x_max: " << x_max << " y_max: " << y_max << endl; 
    
  }
  
  /** Load file data fields in the class. 
   *  
   * @param netcdf_file Gerris(R) file.
   * @param nome_var Data column header name. 
   * @param T **not used
   * @param nomeExtension **not used, is NULL
   * @return void
   */
  void Load(const char netcdf_file[], const char nome_var[], size_t T, const char nomeExtension[])
  {
    const char* cc_nf;
    ifstream f;
    std::stringstream stream(std::stringstream::in | std::stringstream::out);
    int i, ii ;
    string s, ss;
    
    f.open(netcdf_file);
    if(!f){
      cout << "Il file non esiste!\n" << endl;  //todo eccezione
      return ;
    }

    global.log << "Loading field: " << netcdf_file <<" : " << nome_var << " : " << T << endl;

    //data.resize(n_punti);

    if (f.good()) getline(f,s);       // -- Prima riga è l'header, nn lo uso

    for(i=0;i<n_punti && f.good();i++) //todo intertire i cicli
      {  
	//todo gestione eccezione su f
	std::getline(f,s);
	stream << s; 
	for(ii=0; ii<Header.size() && stream.good(); ii++)
	  {
	    stream >> ss;

	    //cout << Header[ii] << " " << campo << endl;

	    if(Header[ii].compare(nome_var)==0)
	      {
		//cout << Header[ii] << " bis " << nome_var << endl;
		data[i]=atof(ss.c_str());
		ii=Header.size();
		stream.str(""); // -- Azzero la parte rimanente dello stream
	      }
	  }
      }
    f.close();

    /* dbg --
    cout << " dbg 0 nome_var: " << campo << endl; 
    for(i=0;i<100;i++){
      cout << " data[" << i << "]= " << data[i] << " " ;
      }
    */

  }

  double getXmax(){ return x_max;};  
  double getXmin(){ return x_min;};

  /** It gives a interpolation of data-value on the input point
   * @return data interpolated
   */
  double Value(double x, double y, double z){
    //  double Value(Vec3<double>& ptoIn){
    
    double v1;
    VecDoub p(2);

    //p[0] = ptoIn.X();
    p[0] = x;
    //p[1] = ptoIn.Y(); 
    p[1] = y;
       
    Shep_interp xx(PtiGriglia, data, 4); //todo il 4 deve diventare un parametro ! (riguarda 
    v1=xx.interp(p);
    return v1;
  }

  double Value(double x, double y, double z, Neighbours &vic)
  {  
    double v1;
    VecDoub p(2);

    //p[0] = ptoIn.X();
    p[0] = x;
    //p[1] = ptoIn.Y();
    p[1] = y;

    Shep_interp xx(PtiGriglia, data, 4); //todo il 4 deve diventare un parametro ! 
    v1=xx.interp(p, vic);
       
    return v1;
  }

  size_t NTimes(const char netcdf_file[]) { return 1;} // i file Gerris DEVONO essere monotemporali

  MatDoub& Points4Neighbour(bool& isUsable) 
  {
    isUsable = true;
    return PtiGriglia;
  }
  
};

/* @class ScalarFiels3D_Mask_Gerris
 */

class ScalarField3D_Mask_Gerris : public ScalarField3D_Mask, public ScalarField3D_Gerris{

public:

  //size_t nx_orography, ny_orography;
  size_t n_p;
  bool orography;
  VecDoub Orografia;
  static const double THRESHOLD = -0.1; //todo su Global.C ?
  static const int N_VICINI = 16;
  
  ScalarField3D_Mask_Gerris() : ScalarField3D_Mask(), ScalarField3D_Gerris(){} 

  void Init(const char netcdf_file[], const char nome_lon[], const char nome_lat[], \
       const char nome_depth[], const char nomeExtension[])  
  {
    ScalarField3D_Gerris::Init(netcdf_file, nome_lon, nome_lat, nome_depth, nomeExtension);
    
    n_p=ScalarField3D_Gerris::n_punti;
  }

  void CloseDivergence(ScalarField3D& U, ScalarField3D& V){};

  void Load(const char netcdf_file[], const char nome_var[], size_t T, const char nomeExtension[])
  {
    cout << "Mask non permetta il Load standard, utilizzare ScalarField3D_Mask_[type]::Mask(....)\n\n";
    throw int(90);
  }


  double Value(double X, double Y, double Z) 
  {
    return ScalarField3D_Gerris::Value(X, Y, Z);
  }

  double Value(double X, double Y, double Z, Neighbours &vicini)
  {
    return ScalarField3D_Gerris::Value(X, Y, Z, vicini);
  }

  MatDoub& Points4Neighbour(bool& isUsable) 
  {
    return ScalarField3D_Gerris::Points4Neighbour(isUsable);
  }

  size_t NTimes(const char netcdf_file[])
  {
    return ScalarField3D_Gerris::NTimes(netcdf_file);
  }


  
  /*  bool Inside(Vec3<double>& p)
  {
    // -- tanto l'interpolazione si puo' fare sempre LP e FMG 26 mar 2015 
    return true;
    }*/


  void LoadOrografia(const char netcdf_file[], const char nome_var[]="3:z")
  {
    long int i,ii,k, ToT;
    int retval, status;
    ifstream f;
    std::stringstream stream(std::stringstream::in | std::stringstream::out);
    string s, ss;
    
    f.open(netcdf_file);
    if(!f){
      cout << "Il file non esiste!\n" << endl;  //todo eccezione
      return ;
    }
    
    global.log << "# Loading bathymetry from: " << netcdf_file << " :  " << nome_var << endl;  
    cout << "# Loading bathymetry from: " << netcdf_file << " :  " << nome_var << endl;  

    Orografia.resize(n_p);
    if (f.good()) getline(f,s);       // -- Prima riga è l'header, nn lo uso
    
    for(i=0;i<n_p && f.good();i++) //todo intertire i cicli
      {  
	//todo gestione eccezione su f
	std::getline(f,s);
	stream << s; 
	for(ii=0; ii<Header.size() && stream.good(); ii++)
	  {
	    stream >> ss;
	    if(Header[ii].compare(nome_var)==0)
	      {
		Orografia[i]=atof(ss.c_str());
		ii=Header.size();
		stream.str(""); 
	      }
	  }
      }
    f.close();
  }


  void Mask(const char netcdf_file[], const char nome_var[])
  {
    long int i,ii,k;
    ifstream f;
    std::stringstream stream(std::stringstream::in | std::stringstream::out);
    string s, ss;


    global.log << "# Loading Mask from: " << netcdf_file << " :  " << nome_var << endl;  
    cout << "# Loading Mask from: " << netcdf_file << " :  " << nome_var << endl; 
    
    f.open(netcdf_file);
    if(!f){
      cout << "Il file non esiste!\n" << endl;  //todo eccezione
      return ;
    }
    if (f.good()) getline(f,s);       // -- Prima riga è l'header, nn lo uso

    cout << " *** n_p: " << n_p << endl;

    for(i=0;i<n_p && f.good();i++) //todo intertire i cicli
      {  
	//todo gestione eccezione su f
	std::getline(f,s);
	stream << s; 
	for(ii=0; ii<Header.size() && stream.good(); ii++)
	  {
	    stream >> ss;
	    if(Header[ii].compare("3:z")==0)
	      {
		if(atof(ss.c_str()) > THRESHOLD){
		  data[i] = 0 ;
		  cout << " mask 1 , i: " << i << endl; //kiki 
		}
		else{
		  data[i] = 1;
		  cout << " mask 0 , i: " << i << endl;
		}
		ii=Header.size();
		stream.str(""); 
	      }
	  }
      }
    f.close();
  }
  
  
  double toplevel(Vec3<double>& p, Neighbours &vicini)
  { 
    return 0.;
  } 

  double toplevel(Vec3<double>& p)
  { 
    /*Neighbours dummy;
      return toplevel(p,dummy);*/
    return 0.;
  } 

  double bottom(Vec3<double>& p)
  {
    Neighbours vic(N_VICINI);
    vic.Init(ScalarField3D_Gerris::PtiGriglia, p); 
    return bottom(p, vic);
  }


  double bottom(Vec3<double>& p, Neighbours& vicini)
  { 
    double v1;
    VecDoub pp(2);
    
    pp[0] = p.X();
    pp[1] = p.Y();
    
    Shep_interp xx(ScalarField3D_Gerris::PtiGriglia, Orografia, 4); //todo il 4 deve diventare un parametro !  
    v1=xx.interp(pp, vicini);
    if(vicini.checkNews())
      vicini.Init(ScalarField3D_Gerris::PtiGriglia, p);
       
    return v1;
  }

  bool Inside(Vec3<double> &p)
  {
    if(p.X() < x_min || p.X() > x_max) return false;
    if(p.Y() < y_min || p.Y() > y_max) return false;
    
    return true;
  }

}; 

