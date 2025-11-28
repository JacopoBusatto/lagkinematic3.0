/** 
 * @file Gerris.C 
 * It manages the output of the simulator Gerris(R).
 * @author Fabio Grasso, Luigi Palatella
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
  inline void Init(const char netcdf_file[], const char nome_lon[],	\
		    const char nome_lat[], const char nome_depth[], const char nomeExtension[]) {

    string s,ss;
    string s1, s2, s3;
    int i=0, ii, ind_lon, ind_lat, primo_indice, secondo_indice;
    const char* cc_nf;
    ifstream f;
    std::stringstream stream(std::stringstream::in | std::stringstream::out);
    f.open(netcdf_file);
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
    
    n_punti=0;

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

    //cout << "ind_lon: " << ind_lon << " ind_lat: " << ind_lat << " primo_indice: " << primo_indice << " secondo_indice: " << secondo_indice << endl;

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
  void Load(const char netcdf_file[], const char nome_var[], size_t T, const char nomeExtension[]){
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

  double Value(double x, double y, double z, Neighbours &vic){
    //double Value(Vec3<double>& ptoIn, Neighbours &vic){
    
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
  
  void operator ^ (ScalarField3D_Gerris &B)
  {
    /* operatore di swap senza lavoro, scambia i puntatori*/
    //MatDoub *a,*temp;
    VecDoub a,temp;

    temp = data;
    data = B.data;
    B.data = temp;
  } 

  void ReboundZ(Vec3<double>& p)
  {
    const double e = 1.e-5;
    
    /*if(p.Z() >= z[0]) p.setZ(2. * z[0] - p.Z() - e);
      if(p.Z() <= z[nz - 1]) p.setZ(2. * z[nz -1] - p.Z() + e);*/
    if(p.Z() >= z_max) p.setZ(2. * z_max - p.Z() - e);
    if(p.Z() <= z_min) p.setZ(2. * z_min - p.Z() + e);

  }
  
  bool Inside(Vec3<double>& p)
  {
    /*if(p.X() < Lon[0] || p.X() > Lon[nx-1]) return false;
      if(p.Y() < Lat[0] || p.Y() > Lat[ny-1]) return false;*/
    /*if(p.X() < PtiGriglia[0][0] || p.X() > PtiGriglia[n_punti-1][0]) return false;
    if(p.Y() < PtiGriglia[0][1] || p.Y() > PtiGriglia[n_punti-1][1]) return false;
    */
    //fmg todo return false se la distanza minima di p da tutti i punti è maggiore di una soglia (20.0)
    //altrimnti true
    return true;
  }

  double Time(){return time;}

  void Time(double t){time = t;}

  Vec3<double> Gradient(Vec3<double> X)
  {
    Vec3<double> g, dx(1,0,0), dy(0,1,0);
    Doub gx, gy, h = 0.1, N;
    Vec3<double> xxx1, xxx2,  yyy1, yyy2;  //fmg
    
    /*
    //    if(X.X() - h < Lon[0] || X.X() + h > Lon[nx-1])   // todo se è dentro
    if(X.X() - h < PtiGriglia[0][0] || X.X() + h > PtiGriglia[n_punti-1][0])
      {
	global.GestoreErrori("Gradiente troppo vicino al bordo in longitudine\n");
	throw int(17);
      }
      
    //if(X.Y() - h < Lat[0] || X.Y() + h > Lat[ny-1])
    if(X.Y() - h < PtiGriglia[0][1] || X.Y() + h > PtiGriglia[n_punti-1][1])
      {
	global.GestoreErrori("Gradiente troppo vicino al bordo in latitudine\n");
	throw int(18);
      }
    */

    /*gx = 0.5 / h * (Value(X + dx * h) - Value(X + dx * (-h)));
      gy = 0.5 / h * (Value(X + dy * h) - Value(X + dy * (-h))); */ 
    xxx1 = X + dx * h; xxx2 = X + dx * (-h) ;
    yyy1 = X + dy * h; yyy2 = X + dy * (-h) ;
    //gx = 0.5 / h * (Value(xxx1) - Value(xxx2));
    //gy = 0.5 / h * (Value(yyy1) - Value(yyy2));
    gx = 0.5 / h * (Value(xxx1.X(),xxx1.Y(),xxx1.Z()) - Value(xxx2.X(),xxx2.Y(),xxx2.Z()));
    gy = 0.5 / h * (Value(yyy1.X(),yyy1.Y(),yyy1.Z()) - Value(yyy2.X(),yyy2.Y(),yyy2.Z()));
    
    //    cout << gx << " " << gy << endl;

    //    N = sqrt(gx * gx + gy * gy);

    g.setX(gx); g.setY(gy);
    g.setZ(0.);

    return g;
  }

  /*size_t Nz(){ return nz;}
  size_t Ny(){ return ny;}
  size_t Nx(){ return nx;}
  */

  //todo non servono!?
  size_t Nz(){ return n_punti;}
  
  size_t Ny(){ return n_punti;}
  
  size_t Nx(){ return n_punti;}

  size_t Nxslab()
  {
    //return nx * ny;
    return n_punti;
  }
  
  //size_t Ntot(){ return nx * ny * nz;}
  size_t Ntot(){ return n_punti; }
  
  void CloseDivergence(ScalarField3D_Gerris& U, ScalarField3D_Gerris& V)
  {
    size_t k, i, j;
    double dwdz, wz, hx, hy, hz;
    
    //    for(i=1; i<Nx(); i++) 
    for(i=1; i<n_punti; i++) 
      { 
	// for(j=1; j<Ny(); j++) 
	for(j=1; j<n_punti; j++) 
	  {
	    wz = 0;
	    /*hx = 1./(U.Lon[i] - U.Lon[i-1]) / 111000. / cos (6.28318530717959 * Lat[j] / 360.) ;
	      hy = 1./(V.Lat[j] - V.Lat[j-1]) / 111000.; */
	    hx = 1./(U.PtiGriglia[i][0] - U.PtiGriglia[i-1][0]) / 111000. / cos (6.28318530717959 * PtiGriglia[j][1] / 360.) ;
	    hy = 1./(V.PtiGriglia[j][1] - V.PtiGriglia[j-1][1]) / 111000.;
	    
	    //todo controllare con Luigi !!! per ora azzardo
	    //(*(data+nz-1))[i][j] = 0.;
	    data[i] = 0.;
	    //todo CHIEDERE A LUIGI !	    
	    //for(k=nz-1; k>0; k--)
	    for(k=n_punti-1; k>0; k--)
	      {
		//hz = z[k] - z[k - 1];
		hz = PtiGriglia[k][2] - PtiGriglia[k - 1][2];
		//wz -= 0.5 * hz * hx * ( U.Value(i,j,k) - U.Value(i-1,j,k) );
		wz -= 0.5 * hz * hx *  U.data[k];
		//wz -= 0.5 * hz * hy * ( V.Value(i,j,k) - V.Value(i,j-1,k) );
		wz -= 0.5 * hz * hy *  V.data[k] - V.data[k] ;
		//(*(data+k-1))[i][j] = wz;
		data[k] = wz;
		
		//	        cout << hx << " " << hy << " " << hz << " " << wz << endl;
		//		cout.flush();
	      }
	  }
      }
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
  static const double THRESHOLD = -0.1; //todo su Global.C
  static const int N_VICINI = 16;
  
  ScalarField3D_Mask_Gerris() : ScalarField3D_Mask(), ScalarField3D_Gerris(){} 

  void Init(const char netcdf_file[], const char nome_lon[], const char nome_lat[], \
       const char nome_depth[], const char nomeExtension[])  
  {
    ScalarField3D_Gerris::Init(netcdf_file, nome_lon, nome_lat, nome_depth, nomeExtension);
    
    n_p=ScalarField3D_Gerris::n_punti;
  }

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


  void Mask(const char netcdf_file[], const char nome_var[]="6:U")
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
	    if(Header[ii].compare(nome_var)==0)
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

