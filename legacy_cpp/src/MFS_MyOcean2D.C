/* Interfaccia campi clorofilla OceanColor */

// Derivata da ScalarField3D
class ScalarField2D_OcColor : public ScalarField3D
{

private:

public:
  VecDoub Lon, Lat; 
  MatDoub *data;

  size_t nx, ny, mx, my;
  static const size_t nz = 1; // Da qui in poi dovrebbe essere valido

  ScalarField2D_OcColor() : ScalarField3D(){} 

  void Init(const char netcdf_file[], const char nome_lon[], \
	    const char nome_lat[], const char nome_depth[], const char nomeExtension[]);

  void Load(const char netcdf_file[], const char nome_var[], size_t T, const char nomeExtension[]);

  void CloseDivergence(ScalarField3D& U, ScalarField3D& V){};

  double Value(double, double, double); //Qui è solo dichiarata

  double Value(double x, double y, double z, Neighbours& neigh) 
  { return Value(x,y,z); }

  /* quando la griglia e' regolare il load con Neighbours e' banalmente overloaded */

  MatDoub& Points4Neighbour(bool& isUsable) 
  {
    isUsable = false;
    return data[0];
  }
  
  size_t Ny() { return nx; }
  size_t Nx() { return ny; }
  //Modificata \RC
  size_t Nz() { return nz;}

  //Serve? non credo
  double Sigma(size_t k)
  { 
    if(Lista.sigma) 
      {
	cout << "il campo sigma non puo' essere utilizzato\n";
	throw int(2343);
      }
    return 0.;
  }
  
  //Modificare?
  double Value(size_t i, size_t j, size_t k)
  {
    //k = 0;
    if(i>=0 && i<nx && j>=0 && j<ny) 
      return (*(data))[i][j];
    else 
      {
	global.GestoreErrori("Valore oltre il range 0 - nx-1 in Field3D.value(i,j,k): Fields3D.C");	
	throw int(30);
      }
  }

  size_t NTimes(const char netcdf_file[])
  {
    int status, ncid, latid, recid;
    size_t length, recs;
    char recname[NC_MAX_NAME+1];
        
    status = nc_open(netcdf_file, NC_SHARE, &ncid);

    status = nc_inq_dimid(ncid, "time", &recid); 

    status = nc_inq_dimlen(ncid, recid, &length);

    cout << "numero of timesssss " << length << endl;

    nc_close(ncid);

    return length;
  }
    
};

void ScalarField2D_OcColor::Init(const char netcdf_file[], const char nome_lon[], const char nome_lat[],\
			      const char nome_depth[], const char nomeExtension[])
  {
    long int i,j,k;
    int ncid, var_id, long_id, lat_id, z_id, Time_id;
    int retval, status;
    
    mx = my = 5;

    cout << "Initializing file of OcColor_MyOcean kind - Only 2D\n\n";
    cout.flush();

    if ((retval = nc_open(netcdf_file, NC_SHARE, &ncid))) 
      {
	global.er << "File NetCDF non trovato o di formato non corretto\n\n";
	throw int(1);
      }

    status = nc_inq_varid(ncid, nome_lon, &long_id);
    if(status) { global.er << "Variabile " << nome_lon << " non trovata\n"; throw int(2);}

    status = nc_inq_varid(ncid, nome_lat, &lat_id);
    if(status) { global.er << "Variabile " << nome_lat << " non trovata\n"; throw int(2);}

    //    status = nc_inq_varid(ncid, nome_depth, &z_id);
    //    if(status) { global.er << "Variabile " << nome_depth << " non trovata\n"; throw int(2);}

    //    Variabili spaziali lon e lat settate, non c'è depth
    
    double *lats, *lons, v;
    nc_type type;   
    int ndims, natts, ndims2, ndims3;
    int  dimids[NC_MAX_VAR_DIMS]; 
    size_t *count, *count2, lenghtp, *index;

    count = (size_t*)malloc(sizeof(size_t)*ndims);
    count2 = (size_t*)malloc(sizeof(size_t)*2);

    status = nc_inq_var(ncid, long_id, 0, &type, &ndims, dimids, &natts);

    cout << "ndims = " << ndims << endl;

    if(ndims == 1) 
      {
	nc_inq_dimlen(ncid, dimids[0], &nx);

	Lon.resize(nx);

	for(size_t i=0; i<nx; i++)
	  {
	    status = nc_get_var1_double(ncid, long_id, &i, &v);
	    Lon[i] = v;
	  }

	status = nc_inq_var(ncid, lat_id, 0, &type, &ndims2, dimids, &natts);
	nc_inq_dimlen(ncid, dimids[0], &ny);

	Lat.resize(ny);

	for(size_t i=0; i<ny; i++)
	  {
	    status = nc_get_var1_double(ncid, lat_id, &i, &v);
	    Lat[i]=v;
	  }
      }
    
    if(ndims == 2) 
      {
	index = new size_t[2];

	nc_inq_dimlen(ncid, dimids[1], &nx);
	nc_inq_dimlen(ncid, dimids[0], &ny);

	Lon.resize(nx); Lat.resize(ny);

	for(size_t l=0; l<nx; l++)
	  {
	    index[0]=0; index[1]=l;
	    status = nc_get_var1_double(ncid, long_id, index, &v);
	    Lon[l] = v;
	  }

	for(size_t j=0; j<ny; j++)
	  {
	    index[0]=j; index[1]=0;
	    status = nc_get_var1_double(ncid, lat_id, index, &v);
	    Lat[j] = v;
	  }
      }

    /*
    status = nc_inq_var(ncid, z_id, 0, &type, &ndims3, dimids, &natts);
    nc_inq_dimlen(ncid, dimids[0], &nz);

    z.resize(nz);

    //Modifico qui, tolto un ciclo inutile \RC

    status = nc_get_var1_double(ncid, z_id, &i, &v);
    
    if(Lista.ocean)
      z[0] = -fabs(v);
    else z[0] = v;

    z[0] = 0.;
    */    

    cout << "# griglia: " << nx << " x " << ny << " x " << nz << endl;
    global.log << "# dimensioni griglia: nx x ny x nz: " << nx << " " <<  ny << " " << nz << endl;

    data = new MatDoub[nz];

    for(i=0; i<nz; i++) (data + i)->resize(nx,ny);  //Potrei togliere il ciclo \RC

    n_punti = nx * ny * nz;

    cout << "# numero totale punti: " << n_punti << endl;
    cout << "# Completato caricamento caratteristiche griglia\n";

    free(count);
    free(count2);

    nc_close(ncid);
  }


void ScalarField2D_OcColor::Load(const char netcdf_file[], const char nome_var[], size_t T, const char nomeExtension[])
  {
    long int i,j,k, ToT;
    size_t nt;
    double v;
    int ncid, var_id, long_id, lat_id, z_id, Time_id, x_id;
    int retval, status;
    nc_type type;   
    int ndims, natts, ndims2, ndims3;
 
    int  dimids[NC_MAX_VAR_DIMS]; 
    size_t *count, *count2, *start;
    ptrdiff_t *stride;
   
    //nz = 1; //\RC
 
    count =  (size_t *)malloc(sizeof(size_t)*3);
    start = (size_t *)malloc(sizeof(size_t)*3);
    stride = (ptrdiff_t *)malloc(sizeof(ptrdiff_t)*3);

    status =  nc_open(netcdf_file, NC_SHARE, &ncid);
    if(status) { global.er << "File NetCDF non trovato \n"; throw int(1);}

    status = nc_inq_varid(ncid, nome_var, &x_id);
    if(status) { global.er << "Variabile " << nome_var << " non trovata\n"; throw int(3);}

    status = nc_inq_var(ncid, x_id, 0, &type, &ndims, dimids, &natts);
    if(status) { global.er << "Variabile " << nome_var << " non trovata\n"; throw int(3);}

    status = nc_inq_dimlen(ncid, dimids[0], &nt);
    if(status) { global.er << "Dimensione variabile " << nome_var << " non trovata\n"; throw int(3);}

    ntime = nt;

    if(global.backtraj) ToT = ScalarField2D_OcColor::ntime - T - 1;
    else ToT = T;

    //    cout << "# checking " << nt << " " << ntime << " " << ScalarField2D_OcColor::ntime << endl;

    if(ScalarField2D_OcColor::ntime <= ToT)
      {
	cout << ScalarField2D_OcColor::ntime << " " << ToT << endl;
	global.GestoreErrori("Indice tempo T out of range: Fields3D.C");
	throw int(20);
      }

    global.log << "# Loading: " << netcdf_file << " :  " << nome_var << " : " << ToT << endl;  
    cout << "# Loading: " << netcdf_file << " :  " << nome_var << " : " << ToT << endl;  

    start[0]=ToT;
    start[1]=start[2]=0;

    count[0]=1; count[1]= ny; count[2]= nx;
    stride[0] = stride[1] = stride[2] = 1;     

    double *f = new double[nz* ny * nx];

    status = nc_get_vars_double(ncid, x_id, start, count, stride, f);
    
    if(status)
      {
	global.log << "Probabili problemi nel caricamento dati \n" << endl;
	cout << "Probabili problemi nel caricamento dati \n" << endl;
	throw int(21);
      }

    //\RC
    k = 0;
    for(i=0; i<nx; i++) 
      {
	for(j=0; j<ny; j++) 
	  {
	    v = *(f + i + nx*j + nx*ny*k);
	    if(v < 0.)  
	      {
		if (j>0 && (*(data+k))[i][j-1]>0.)
		  (*(data+k))[i][j] =(*(data+k))[i][j-1];
		else (*(data+k))[i][j]=0.;
	      }
	    else (*(data+k))[i][j] = v;
	  }
      }
    
    delete[] f;

    free(count);
    free(start);
    free(stride);

    nc_close(ncid);
  }

//Modificata \RC
double ScalarField2D_OcColor :: Value(double X, double Y, double Z) 
{
  size_t k;
  double v1;

  k=0;
  
  if(X < Lon[0] || X > Lon[nx-1] || Y < Lat[0] || Y > Lat[ny-1]) 
    {
      //      global.GestoreErrori("Valore oltre il range 0 - nx-1 in Field3D.value(x,y,z): Fields3D.C");
      cout << "Valore ooooltre il range 0 - nx-1 in Field2D.value(x,y): Fields2D.C\n" ;
      cout << X << " " << Lon[0] << endl;
      cout << Y << " " << Lat[0] << endl;
      return 0.; /* todo DA CONTROLLARE */
      //      throw int(33);
    }

  Poly2D_interp suz1(Lon, Lat, *(data + k), mx, my);

  v1 = suz1.interp(X,Y);
  
  return v1;
}


/* %%%%%%%%%%%% QUI INIZIA LA MASK %%%%%%%%%%%%%%%%%%%%%%%%%%*/

class ScalarField2D_Mask_OcColor : public ScalarField3D_Mask, public ScalarField2D_OcColor
{

public:
  static const size_t nz = 1;

  ScalarField2D_Mask_OcColor() : ScalarField3D_Mask(), ScalarField2D_OcColor(){} 

  void Init(const char netcdf_file[], const char nome_lon[], const char nome_lat[], \
       const char nome_depth[], const char nomeExtension[])  
  {
    ScalarField2D_OcColor::Init(netcdf_file, nome_lon, nome_lat, nome_depth, "NULL");
  }

  void CloseDivergence(ScalarField3D& U, ScalarField3D& V){};

  void Load(const char netcdf_file[], const char nome_var[], size_t T, const char nomeExtension[])
  {
    cout << "Mask non permette il Load standard, utilizzare ScalarField3D_Mask_[type]::Mask(....)\n\n";
    throw int(90);
  }

  double Value(double X, double Y, double Z) 
  {
    return ScalarField2D_OcColor::Value(X, Y, Z);
  }

  double Value(double X, double Y, double Z, Neighbours &vicini)
  {
    return ScalarField2D_OcColor::Value(X, Y, Z);
  }

  size_t NTimes(const char netcdf_file[])
  {
    return ScalarField2D_OcColor::NTimes(netcdf_file);
  }

  //Non mi serve?
  MatDoub& Points4Neighbour(bool& isUsable) 
  {
    return ScalarField2D_OcColor::Points4Neighbour(isUsable);
  }

  void LoadOrografia(const char netcdf_file[], const char nome_var[]) {}; /* MFS non ha l'orografia a parte */

  void Mask(const char netcdf_file[], const char nome_var[])
  {
    long int i,j,k, ToT;
    double v;
    int ncid, var_id, long_id, lat_id, z_id, Time_id, x_id;
    int retval, status;
    nc_type type;   
    int ndims, natts, ndims2, ndims3;
    size_t nt;
    int  dimids[NC_MAX_VAR_DIMS]; 
    size_t *count, *count2, *start;
    ptrdiff_t *stride;
    
    //nz = 1; //\RC

    count =  (size_t *)malloc(sizeof(size_t)*3);
    start = (size_t *)malloc(sizeof(size_t)*3);
    stride = (ptrdiff_t *)malloc(sizeof(ptrdiff_t)*3);

    status =  nc_open(netcdf_file, NC_SHARE, &ncid);
    if(status) { global.er << "File NetCDF non trovato \n"; throw int(1);}

    status = nc_inq_varid(ncid, nome_var, &x_id);
    if(status) { global.er << "Variabile " << nome_var << " non trovata\n"; throw int(3);}

    status = nc_inq_var(ncid, x_id, 0, &type, &ndims, dimids, &natts);
    if(status) { global.er << "Variabile " << nome_var << " non trovata\n"; throw int(3);}

    status = nc_inq_dimlen(ncid, dimids[0], &nt);
    if(status) { global.er << "Dimensione variabile " << nome_var << " non trovata\n"; throw int(3);}

    ScalarField2D_OcColor::ntime = nt;

    cout << "checking in mask" << nt << " " << ScalarField2D_OcColor::ntime << endl;

    global.log << "# Loading mask from: " << netcdf_file << " based on  " << nome_var << endl;  

    start[0]=0;
    start[1]=start[2]=0;

    count[0]=1; count[1]=ny; count[2]= nx;
    stride[0] = stride[1] = stride[2] = 1;

    double *f = new double[nz * ny * nx];

    status = nc_get_vars_double(ncid, x_id, start, count, stride, f);
    
    if(status)
      {
	global.log << "Probabili problemi nel caricamento dati \n" << endl;
	cout << "Probabili problemi nel caricamento dati \n" << endl;
	throw int(21);
      }

    k=0;
    for(i=0; i<nx; i++) 
      {
	for(j=0; j<ny; j++) 
	  {
	    v = *(f + i + nx*j + nx*ny*k);
	    if(v > 1.e6)  (*(data+k))[i][j] = 0.;
	    else (*(data+k))[i][j] = 1.;
	  }
      }
  
    delete[] f;

    free(count);
    free(start);
    free(stride);

    nc_close(ncid);
  }

  double toplevel(Vec3<double>& p)
  {     
    //    return ScalarField2D_OcColor::z[0]; 
    return 0.;
  } /* nel caso oceano e' sempre il primo livello*/

  double toplevel(Vec3<double>& p, Neighbours& vicini)
  {
    return toplevel(p);
  }

  double bottom(Vec3<double> &p)
  {
    return -10000.;
  }

  double bottom(Vec3<double> &p, Neighbours& vicini)
  {
    return -10000.;
  }



  //Tolto bottom
  bool Inside(Vec3<double> &p)
  {
    if(p.X() < Lon[0] || p.X() > Lon[nx-1]) return false;
    if(p.Y() < Lat[0] || p.Y() > Lat[ny-1]) return false;
    
    return true;
  }

}; 
