// Interfaccia Mercator \RC

class ScalarField3D_Mercator : public ScalarField3D {

private: 

public: 
  VecDoub Lon, Lat;
  MatDoub *data;

  size_t nx, ny, mx, my; 
  double scale_factor;

  ScalarField3D_Mercator() : ScalarField3D(){}

  void Init(const char netcdf_file[], const char nome_lon[], \
	    const char nome_lat[], const char nome_depth[], const char nomeExtension[]);

  void Load(const char netcdf_file[], const char nome_var[], size_t T, const char nomeExtension[]);
  
  double Value(double, double, double); 

  void CloseDivergence(ScalarField3D& U, ScalarField3D& V);
  
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
  
  double Sigma(size_t k){ 
    if(Lista.sigma) 
      {
	cout << "il campo sigma non puo' essere utilizzato\n";
	throw int(2343);
      }
    return 0.;
  }
  
  double Value(size_t i, size_t j, size_t k)
  {
    if(i>=0 && i<nx && j>=0 && j<ny && k >= 0 && k < nz) 
      return (*(data+k))[i][j];

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

    //    status = nc_inq_unlimdim(ncid, &recid); 

    //    status = nc_inq_dimid(ncid, "time_counter", &recid); 
    status = nc_inq_dimid(ncid, "time", &recid); 

    status = nc_inq_dimlen(ncid, recid, &length);

    cout << "numero of timesssss " << length << endl;

    return length;
  }
  
};

////////////////////

void ScalarField3D_Mercator::Init(const char netcdf_file[], const char nome_lon[], const char nome_lat[], \
				  const char nome_depth[], const char nomeExtension[])
  {
    /*nome_var serve per calcolare la mask*/
    long int i,j,k;
    int ncid, var_id, long_id, lat_id, z_id, Time_id;
    int retval, status;
    
    mx = my = 5;

    cout << "Initializing file of Mercator kind - Test\n\n";
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

    status = nc_inq_varid(ncid, nome_depth, &z_id);
    if(status) { global.er << "Variabile " << nome_depth << " non trovata\n"; throw int(2);}

    //mask = false;

    double *lats, *lons, *zz, v;
    nc_type type;   
    int ndims, natts, ndims2, ndims3;
    int  dimids[NC_MAX_VAR_DIMS]; 
    size_t *count, *count2, lenghtp, *index;

    count = (size_t*)malloc(sizeof(size_t)*ndims);
    count2 = (size_t*)malloc(sizeof(size_t)*2);

    status = nc_inq_var(ncid, long_id, 0, &type, &ndims, dimids, &natts);

    cout << "mercator ndims = " << ndims << endl; 
    //global.log << "ndims = " << ndims << endl; 

    if(ndims == 1) 
      {
	nc_inq_dimlen(ncid, dimids[0], &nx);

	Lon.resize(nx);

	for(size_t i=0; i<nx; i++)
	  {
	    status = nc_get_var1_double(ncid, long_id, &i, &v);
	    Lon[i] = v; 
	    //if(v>=0) Lon[i] = v;
	    //if(v<0) Lon[i] = v + 360.;
	    
	    /* LP per gestire la il taglio con simulazioni vicine a 180 gradi di longitudine*/ 
	    
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
	    Lon[i] = v; 
	    //if(v>=0) Lon[i] = v;
            //if(v<0) Lon[i] = v + 360.;
	    
	    /* LP per gestire la il taglio con simulazioni vicine a 180 gradi di longitudine*/ 
	  }

	for(size_t j=0; j<ny; j++)
	  {
	    index[0]=j; index[1]=0;
	    status = nc_get_var1_double(ncid, lat_id, index, &v);
	    Lat[j] = v;
	  }
      }
 
    status = nc_inq_var(ncid, z_id, 0, &type, &ndims3, dimids, &natts);
    nc_inq_dimlen(ncid, dimids[0], &nz);

    z.resize(nz); // tanto nz = 1 const

    for(size_t i=0; i<nz; i++)
      {
	status = nc_get_var1_double(ncid, z_id, &i, &v);
	if(Lista.ocean)
	  z[i] = -fabs(v);
	else z[i] = v;
      }
    
    cout << "# griglia: " << nx << " x " << ny << " x " << nz << endl;
    global.log << "# dimensioni griglia: nx x ny x nz: " << nx << " " <<  ny << " " << nz << endl;

    data = new MatDoub[nz];

    for(i=0; i<nz; i++) (data + i)->resize(nx,ny);

    n_punti = nx *ny * nz;

    cout << "# numero totale punti: " << n_punti << endl;
    cout << "# Completato caricamento caratteristiche griglia\n";

    free(count);
    free(count2);

    nc_close(ncid);
  }



void ScalarField3D_Mercator::Load(const char netcdf_file[], const char nome_var[], size_t T, const char nomeExtension[])
  {
    long int i,j,k, ToT;
    size_t nt;
    double v, sc;
    int ncid, var_id, long_id, lat_id, z_id, Time_id, x_id;
    int retval, status;
    nc_type type;   
    int ndims, natts, ndims2, ndims3;

    int  dimids[NC_MAX_VAR_DIMS]; 
    size_t *count, *count2, *start;
    ptrdiff_t *stride;
    


    count =  (size_t *)malloc(sizeof(size_t)*4);
    start = (size_t *)malloc(sizeof(size_t)*4);
    stride = (ptrdiff_t *)malloc(sizeof(ptrdiff_t)*4);

    status = nc_open(netcdf_file, NC_SHARE, &ncid);
    if(status) { global.er << "File NetCDF non trovato \n"; throw int(1);}

    status = nc_inq_varid(ncid, nome_var, &x_id);
    if(status) { global.er << "Variabile " << nome_var << " non trovata\n"; throw int(3);}

    status = nc_get_att_double(ncid, x_id,"scale_factor",&sc);
    cout << "scale_factor = " << sc << endl;

    status = nc_inq_var(ncid, x_id, 0, &type, &ndims, dimids, &natts);
    if(status) { global.er << "Variabile " << nome_var << " non trovata\n"; throw int(3);} 

    status = nc_inq_dimlen(ncid, dimids[0], &nt);
    if(status) { global.er << "Dimensione variabile " << nome_var << " non trovata\n"; throw int(3);} 
    ntime = nt;

    if(global.backtraj) ToT = ScalarField3D_Mercator::ntime - T - 1;
    else ToT = T;

    if(ScalarField3D_Mercator::ntime <= ToT)
      {
	cout << ScalarField3D_Mercator::ntime << " " << ToT << endl;
	global.GestoreErrori("Indice tempo T out of range: Fields3D.C");
	throw int(20);
      }

    global.log << "# Loading: " << netcdf_file << " :  " << nome_var << " : " << ToT << endl;  

    start[0]=ToT;
    start[1]=start[2]=start[3]=0;

    count[0]=1; count[1]= nz; count[2]= ny; count[3]= nx;
    stride[0] = stride[1] = stride[2] = stride[3] = 1;     

    short *f = new short[nx * ny * nz];

    status = nc_get_vars_short(ncid, x_id, start, count, stride, f);
    
    if(status)
      {
	global.log << "Probabili problemi nel caricamento dati \n" << endl;
	cout << "Probabili problemi nel caricamento dati \n" << endl;
	throw int(22);
      }

    for(k=0; k<nz; k++)
      { 
	for(i=0; i<nx; i++) 
	  {
	    for(j=0; j<ny; j++) 
	      {
		v = sc * (*(f + i + nx*j + nx*ny*k));
		if(fabs(v) > 1.e5)  (*(data+k))[i][j] = 0.;
		else (*(data+k))[i][j] = v;
	      }
	  }
      }

    delete[] f;
    free(count);
    free(start);
    free(stride);

    nc_close(ncid);
  }



void ScalarField3D_Mercator::CloseDivergence(ScalarField3D& U, ScalarField3D& V)
{
  size_t k, i, j;
  double dwdz, wz, hx, hy, hz, ddx = 1./24;
  double XX, YY;

  for(i=1; i<nx-1; i++) 
    { 
      cout << 100.*i/nx << endl;

      for(j=1; j<ny-1; j++) 
	{
	  wz = 0;
	  
	  (*(data+nz-1))[i][j] = 0.;
	  
	  XX = Lon[i];
	  YY = Lat[j];

	  hy = 110000. * ddx; // 1/24 di grado
	  hx = 110000. * ddx * cos(YY*6.28318530717959 / 180.); // 1/24 di grado


	  for(k=nz-1; k>0; k--)
	    {
	      hz = z[k] - z[k - 1];
	      wz -= 0.5 * hz * ( U.Value(XX+ddx,YY,z[k]) - U.Value(XX-ddx,YY,z[k]) ) / hx;
	      wz -= 0.5 * hz * ( V.Value(XX,YY+ddx,z[k]) - V.Value(XX,YY-ddx,z[k]) ) / hy;
	      (*(data+k-1))[i][j] = wz;
	      
	      //cout << z[k] << " " << wz << " " << endl;
	      //              cout.flush();
	    }
	}
    }
}




double ScalarField3D_Mercator::Value(double X, double Y, double Z) 
{
  size_t k, k1;
  double v1, v2;

  k=0;
  
  if(Z >= z[0] || Z < z[nz - 1] || X < Lon[0] || X > Lon[nx-1] || Y < Lat[0] || Y > Lat[ny-1]) 
    {
      //      cout << "Valore ooooltre il range 0 - nx-1 in Field3D.value(x,y,z): Fields3D.C\n" ;
      // cout << Z << " " << z[0] << endl;

      //      cout << X << " " << Lon[0] << endl;
      //  cout << Y << " " << Lat[0] << endl;

      return 0.; /* todo DA CONTROLLARE */
      //      throw int(33);
    }

  while(z[k] >= Z && k < nz-1) k++;

  Poly2D_interp suz1(Lon, Lat, *(data + k), mx, my), suz2(Lon, Lat, *(data + k - 1), mx, my);

  v1 = suz1.interp(X,Y);
  v2 = suz2.interp(X,Y);
  
  return v1 + (Z - z[k])/(z[k-1] - z[k])*(v2 - v1);
}


// %%%%%%%%%%%% QUI INIZIA LA MASK %%%%%%%%%%%%%%%%%%%%%%%%%%

class ScalarField3D_Mask_Mercator : public ScalarField3D_Mask, public ScalarField3D_Mercator
{

public:

 ScalarField3D_Mask_Mercator() : ScalarField3D_Mask(), ScalarField3D_Mercator(){} 

  void Init(const char netcdf_file[], const char nome_lon[], const char nome_lat[], \
	    const char nome_depth[], const char nomeExtension[])  
  {
    ScalarField3D_Mercator::Init(netcdf_file, nome_lon, nome_lat, nome_depth, "NULL");
  }
  
  void Load(const char netcdf_file[], const char nome_var[], size_t T, const char nomeExtension[])
  {
    cout << "Mask non permetta il Load standard, utilizzare ScalarField3D_Mask_[type]::Mask(....)\n\n";
    throw int(90);
  }

  void CloseDivergence(ScalarField3D& U, ScalarField3D& V){};

  double Value(double X, double Y, double Z) 
  {
    return ScalarField3D_Mercator::Value(X, Y, Z);
  }

  double Value(double X, double Y, double Z, Neighbours &vicini)
  {
    return ScalarField3D_Mercator::Value(X, Y, Z);
  }

  size_t NTimes(const char netcdf_file[])
  {
    return ScalarField3D_Mercator::NTimes(netcdf_file);
  }

  MatDoub& Points4Neighbour(bool& isUsable) 
  {
    return ScalarField3D_Mercator::Points4Neighbour(isUsable);
  }

  // Mercator non ha l'orografia a parte ????? \RC 
  void LoadOrografia(const char netcdf_file[], const char nome_var[]) {}; 

  void Mask(const char netcdf_file[], const char nome_var[])
  {
    long int i,j,k, ToT;
    double v, sc;
    int ncid, var_id, long_id, lat_id, z_id, Time_id, x_id;
    int retval, status;
    nc_type type;   
    int ndims, natts, ndims2, ndims3;
    size_t nt;
    int  dimids[NC_MAX_VAR_DIMS]; 
    size_t *count, *count2, *start;
    ptrdiff_t *stride;
    
    count =  (size_t *)malloc(sizeof(size_t)*4);
    start = (size_t *)malloc(sizeof(size_t)*4);
    stride = (ptrdiff_t *)malloc(sizeof(ptrdiff_t)*4);

    status =  nc_open(netcdf_file, NC_SHARE, &ncid);
    if(status) { global.er << "File NetCDF non trovato \n"; throw int(1);}

    status = nc_inq_varid(ncid, nome_var, &x_id);
    if(status) { global.er << "Variabile " << nome_var << " non trovata\n"; throw int(3);}

    status = nc_get_att_double(ncid, x_id,"scale_factor",&sc); // scale factor

    status = nc_inq_var(ncid, x_id, 0, &type, &ndims, dimids, &natts);
    if(status) { global.er << "Variabile " << nome_var << " non trovata\n"; throw int(3);}

    status = nc_inq_dimlen(ncid, dimids[0], &nt);
    if(status) { global.er << "Dimensione variabile " << nome_var << " non trovata\n"; throw int(3);}

    ScalarField3D_Mercator::ntime = nt;

    //    cout << "numero tempi = " <<     ScalarField3D_Mask_Mercator::ntime << endl;

    cout << "checking in mask" << nt << " " << ScalarField3D_Mercator::ntime << endl;

    global.log << "# Loading mask from: " << netcdf_file << " based on  " << nome_var << endl;  

    start[0]=0;
    start[1]=start[2]=start[3]=0;

    count[0]=1; count[1]= ScalarField3D_Mercator::nz; count[2]= ny; count[3]= nx;
    stride[0] = stride[1] = stride[2] = stride[3] = 1;     

    short *f = new short[ ScalarField3D_Mercator::nz * ny * nx];

    status = nc_get_vars_short(ncid, x_id, start, count, stride, f);
    
    if(status)
      {
	global.log << "Probabili problemi nel caricamento dati \n" << endl;
	cout << "Probabili problemi nel caricamento dati \n" << endl;
	throw int(21);
      }

    for(k=0; k<ScalarField3D_Mercator::nz; k++)
      { 
	for(i=0; i<nx; i++) 
	  {
	    for(j=0; j<ny; j++) 
	      {
		v = sc * (*(f + i + nx*j + nx*ny*k));
		if(fabs(v) > 1.e5)  (*(data+k))[i][j] = 0.;
		else (*(data+k))[i][j] = 1.;
	      }
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
    return ScalarField3D_Mercator::z[0]; 
    
  } /* nel caso oceano e' sempre il primo livello*/

  double toplevel(Vec3<double>& p, Neighbours& vicini)
  {
    return toplevel(p);
  }

  double bottom(Vec3<double>& p)
  { 
    double z1, z2, zm, v1, v2, vm, e=1.e-5; 
    
    size_t k =  ScalarField3D_Mercator::nz -1;

    try{
      while( Value(p.X(), p.Y(), ScalarField3D_Mercator::z[k] + e) < global.tt && k > 0 ) k--;  
    }
    catch (int ercode)
      {
	cout << "e' proprio qui in bottom\n\n";
	cout << "z = " << ScalarField3D_Mercator::z[k] << endl;
      }

    return ScalarField3D_Mercator::z[k];
  }

  double bottom(Vec3<double>& p, Neighbours& vicini)
  {
    return bottom(p);
  }

  bool Inside(Vec3<double> &p)
  {

    if(p.Y() < Lat[0] || p.Y() > Lat[ny-1]) return false;

    //    return true; /* per ora */

    //if(fabs(Lon[nx-1]-Lon[0]) - 360. < 0.5 ) return true;
    /*LP 3 dic 2019, se su longitudine si chiude non può uscire */
    if(p.X() < Lon[0] || p.X() > Lon[nx-1]) return false;    
    return true;
  }


};
