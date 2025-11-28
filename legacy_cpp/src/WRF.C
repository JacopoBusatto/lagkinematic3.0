/* Interfaccia WRF */

class ScalarField3D_WRF : public ScalarField3D {

  /* Questa e' signa e atmosferica */

private:
  size_t nx, ny;
  double missing_value;
  

public:
  static const size_t degree = 4;
  double x_min, y_min, z_min, x_max, y_max, z_max;
  VecDoub *data;
  MatDoub *PtiGriglia;

  
  ScalarField3D_WRF() : ScalarField3D(){} 

  void Init(const char netcdf_file[], const char nome_lon[], const char nome_lat[], \
	    const char nome_depth[], const char nomeExtension[]);

  void Load(const char netcdf_file[], const char nome_var[], size_t T, const char nomeExtension[]);

  void CloseDivergence(ScalarField3D& U, ScalarField3D& V){};

  double Value(double, double, double); 

  double Value(double x, double y, double z, Neighbours& neigh);

  /* quando la griglia e' regolare il load con Neighbours e' banalmente overloaded */

  MatDoub& Points4Neighbour(bool& isUsable) 
  {
    isUsable = true;
    return PtiGriglia[0];
  }
  
  
  size_t Ny() { return ny; }
 
  size_t Nx() { return nx; }

  //  size_t Nz() { return nz; }

  double Value(size_t i, size_t j, size_t k)
  {
    size_t index;


    if(i>=0 && i<nx && j>=0 && j<ny && k >= 0 && k < nz) 
      {
	index = i + j*ny;      
	return (*(data+k))[index];
      }

    global.GestoreErrori("Valore oltre il range 0 - nx-1 in Field3D.value(i,j,k): Fields3D.C");	
    throw int(30);
  }


  size_t NTimes(const char netcdf_file[])
  {
    int status, ncid, latid, recid;
    size_t length, recs;
    char recname[NC_MAX_NAME+1];
    
    status = nc_open(netcdf_file, NC_SHARE, &ncid);

    //    status = nc_inq_unlimdim(ncid, &recid); 

    status = nc_inq_dimid(ncid, "Time", &recid); 

    status = nc_inq_dimlen(ncid, recid, &length);

    cout << "numero of timesssss " << length << endl;

    return length;
  }

};



void ScalarField3D_WRF::Init(const char netcdf_file[], const char nome_lon[], const char nome_lat[],\
			     const char nome_depth[], const char nomeExtension[])
 {
    long int i,j,k, iter;
    int ncid, var_id, long_id, lat_id, z_id, Time_id;
    int retval, status;

    cout << "Initializing file of  WRF kind\n\n";
    
    if(!strcmp(nomeExtension,"NULL")) extension = false;
    else extension = true;


    if ((retval = nc_open(netcdf_file, NC_SHARE, &ncid))) 
      {
        global.er << "File NetCDF non trovato o di formato non corretto\n\n";
        cout << "File NetCDF non trovato o di formato non corretto\n\n";
	throw int(2);
      }


    status = nc_inq_varid(ncid, nome_lon, &long_id);
    if(status) { global.er << "Variabile " << nome_lon << " non trovata\n"; global.er.close(); throw int(2);}

    status = nc_inq_varid(ncid, nome_lat, &lat_id);
    if(status) { global.er << "Variabile " << nome_lat << " non trovata\n"; throw int(2);}

    status = nc_inq_varid(ncid, nome_depth, &z_id);
    if(status) { global.er << "Variabile " << nome_depth << " non trovata\n"; throw int(2);}



    double *lats, *lons, *zz, v;
    nc_type type;   
    int ndims, natts, ndims2, ndims3;
    size_t *count, *count2;
    size_t lenghtp, *index, *index2 , ntime;
    int  dimids[NC_MAX_VAR_DIMS]; 
    

    status = nc_inq_var(ncid, long_id, 0, &type, &ndims, dimids, &natts);

    index = new size_t[3];
    index2 = new size_t[2];

    nc_inq_dimlen(ncid, dimids[2], &nx);
    nc_inq_dimlen(ncid, dimids[1], &ny);

    nc_inq_dimlen(ncid, dimids[0], &ntime);

    status = nc_inq_var(ncid, z_id, 0, &type, &ndims3, dimids, &natts);
    nc_inq_dimlen(ncid, dimids[1], &nz);

    cout << "# dimensioni griglia: nx x ny x nz: " << nx << " " <<  ny << " " << nz << endl;
    global.log << "# dimensioni griglia: nx x ny x nz: " << nx << " " <<  ny << " " << nz << endl;

    if(extension)
      {
	nz++;
	global.log << "# Attivata l'estensione a eta = 1 nz = " << nz << endl;
      }
    
    z.resize(nz);

    n_punti = nx * ny * nz;
     
    PtiGriglia = new MatDoub[nz];
    data = new VecDoub[nz];

    for(long k=0; k<nz; k++)
      {
	PtiGriglia[k].resize(nx * ny,2);
	data[k].resize(nx * ny);
      }


    for(long k=0; k<nx * ny; k++)
      {
	/* Qui si riempe solo PtiGriglia[0], ovvero il livello suolo, 
	   con i valori di longitudine e latitudine */
	i = k%nx;
	j = k/nx;
	index[0]=0; index[1]=j; index[2] = i;
	status = nc_get_var1_double(ncid, long_id, index, &v);
	if(status) throw int(21273);
	(*PtiGriglia)[k][0] = v;

	status = nc_get_var1_double(ncid, lat_id, index, &v);
	if(status) throw int(17276);
	(*PtiGriglia)[k][1] = v;
      }

    if(!extension)
      for(size_t k=0; k<nz; k++)
	{
	  index2[0]=0; index2[1]=k;
	  status = nc_get_var1_double(ncid, z_id, index2, &v);
	  
	  z[k] = v;
	}
    else
      {
	for(size_t k=0; k<nz-1; k++)
	  {
	    index2[0]=0; index2[1]=k;
	    status = nc_get_var1_double(ncid, z_id, index2, &v);
	    
	    z[k+1] = v;
	  }
	z[0] = 1.;
      }


    for(size_t k=0; k<nz; k++)
      {
	index2[0]=0; index2[1]=k;
	status = nc_get_var1_double(ncid, z_id, index2, &v);

	z[k] = v;
	//	cout << "z = " << v << endl;

	for(i=0; i< nx * ny; i++)
	    {
	      if(k!=0)
		{
		  (*(PtiGriglia+k))[i][0] = (*PtiGriglia)[i][0];
		  (*(PtiGriglia+k))[i][1] = (*PtiGriglia)[i][1];
		}
	    }
       }

    // -- Aggiorno min e max
    i=0;

    x_min=(*PtiGriglia)[i][0]; 
    y_min=(*PtiGriglia)[i][1]; 
    
    x_max=(*PtiGriglia)[i][0]; 
    y_max=(*PtiGriglia)[i][1]; 
      
    for(i=1; i<nx * ny; i++)
      {	
	if((*PtiGriglia)[i][0]<x_min) x_min=(*PtiGriglia)[i][0]; 
	if((*PtiGriglia)[i][1]<y_min) y_min=(*PtiGriglia)[i][1]; 

	
	if((*PtiGriglia)[i][0]>x_max) x_max=(*PtiGriglia)[i][0]; 
	if((*PtiGriglia)[i][1]>y_max) y_max=(*PtiGriglia)[i][1]; 
      }


    z_min = 0.;
    z_max = 20000.;

    cout << "n_punti: " <<n_punti <<"\n" << endl;
    cout << "x_min: " << x_min << " y_min: " << y_min << endl;
    cout << "x_max: " << x_max << " y_max: " << y_max << endl;
    cout << "sigma_max: " << z_max << " sigma_min: " << z_min << endl;

    free(index);
    free(index2);

    nc_close(ncid);


 }
 

void ScalarField3D_WRF::Load(const char netcdf_file[], const char nome_var[], size_t T, const char nomeExtension[])
  {
    long int i,j,k, ToT;
    double v, missingvalue;
    int ncid, var_id, long_id, lat_id, z_id, Time_id, x_id;
    int retval, status;
    nc_type type;   
    int ndims, natts, ndims2, ndims3;
    size_t ntimes;
    int  dimids[NC_MAX_VAR_DIMS]; 
    size_t *count, *count2, *start;
    ptrdiff_t *stride;
    
    if(!strcmp(nomeExtension,"NULL")) extension = false;
    else extension = true;

    count =  (size_t *)malloc(sizeof(size_t)*4);
    start = (size_t *)malloc(sizeof(size_t)*4);
    stride = (ptrdiff_t *)malloc(sizeof(ptrdiff_t)*4);

    status =  nc_open(netcdf_file, NC_SHARE, &ncid);
    if(status) { global.er << "File NetCDF non trovato \n"; throw int(1);}

    status = nc_inq_varid(ncid, nome_var, &x_id);
    if(status) { global.er << "Variabile " << nome_var << " non trovata\n"; throw int(3);}

    status = nc_inq_var(ncid, x_id, 0, &type, &ndims, dimids, &natts);
    if(status) { global.er << "Variabile " << nome_var << " non trovata\n"; throw int(3);}

    status = nc_inq_dimlen(ncid, dimids[0], &ntimes);
    if(status) { global.er << "Dimensione variabile " << nome_var << " non trovata\n"; throw int(3);}

    ScalarField3D_WRF::ntime = ntimes;

    status = nc_get_att_double(ncid, x_id,"missing_value",&missingvalue);

    if(global.backtraj) ToT = ntimes - T - 1;
    else ToT = T;

    if(ntimes <= ToT)
      {
	global.er  << "Indice tempo T out of range: Fields3D.C\n";
	throw int(20);
      }

    global.log << "# Loading: " << netcdf_file << " :  " << nome_var << " : " << ToT << endl;  
    cout << "# Loading: " << netcdf_file << " :  " << nome_var << " : " << nomeExtension << " : "  << ToT << endl;  

    start[0]=ToT;
    start[1]=start[2]=start[3]=0;

    count[0]=1; count[1]= nz; count[2]= ny; count[3]= nx;
    stride[0] = stride[1] = stride[2] = stride[3] = 1;     

    if(extension) count[1]= nz - 1; // uno meno di quanto allocato! 

    double *f; 

    size_t n_punti_bulk;

    if(extension)
      n_punti_bulk = (nz-1) * ny * nx;
    else 
      n_punti_bulk = nz * ny * nx;
      
    f = new double[n_punti_bulk];
    
    status = nc_get_vars_double(ncid, x_id, start, count, stride, f);
    
    if(status)
      {
	global.er << "Probabili problemi nel caricamento dati \n" << endl;
	cout << "Probabili problemi nel caricamento dati \n" << endl;
	throw int(21);
      }

    size_t level, nxlevel;

    nxlevel = nx * ny;
 
    for(k=0; k<n_punti_bulk; k++)
      { 
	if(f[k] > 1.e6)  f[k]=0.;

	level = k / nxlevel;
	if(extension) level++; // sposto in su di un livello per dar posto all'estensione
	i = k % nxlevel;
	(*(data+level))[i] = f[k];
      }


    delete f;

    free(count);
    free(start);
    free(stride);

    if(extension)
      {
	count =  (size_t *)malloc(sizeof(size_t)*3);
	start = (size_t *)malloc(sizeof(size_t)*3);
	stride = (ptrdiff_t *)malloc(sizeof(ptrdiff_t)*3);
	
	start[0]=ToT;
	start[1]=start[2]=0;
	
	count[0]=1; 
	
	count[1]= ny; count[2]= nx;
	stride[0] = stride[1] = stride[2] = 1;     
	
	status = nc_inq_varid(ncid, nomeExtension, &x_id);
	if(status) throw int(280714);
	
	f = new double[ny * nx];

	status = nc_get_vars_double(ncid, x_id, start, count, stride, f);
	
	for(k=0; k<ny * nx; k++)
	  { 
	    if(f[k] > 1.e6)  f[k]=0.;     
	    (*data)[k] = f[k];
	  }

	delete f;
	free(count);  
	free(start); free(stride);
      }

    nc_close(ncid);
  }


double ScalarField3D_WRF :: Value(double X, double Y, double Z) 
{


  double vup, vdown;
  VecDoub p(2);
  
  p[0] = X;
  p[1] = Y;

  
  if(Z > 1. || Z < z[nz-1])
      {
	cout << "Valore oltre il range z in Field3D.value(x,y,z): Fields3D.C" << endl;
	cout << Z << " " << z[0]  << endl;
	global.er <<  "Valore oltre il range z in Field3D.value(x,y,z)\n";
        throw int(13);
      }
    
    size_t k;
    k=0;

    while(z[k] >= Z && k < nz - 1) k++;
    
    /* i campi li estendo a eta = 1 ponendo il valore del campo = 0 (U,V,W) */

    Shep_interp xx_up(PtiGriglia[k], data[k], degree);

    if(k==0) /* quindi Z > sigma[0]*/
      {
	vup = xx_up.interp(p);
	return vup * (1. - Z) / (1 -  z[0]);
      }

    Shep_interp xx_down(PtiGriglia[k-1], data[k-1], degree);
    

    if(z[k] == Z) return xx_up.interp(p); 
    if(z[k-1] == Z) return xx_down.interp(p);

    vup = xx_up.interp(p);
    vdown = xx_down.interp(p);
    
    return vup + (Z - z[k])/(z[k-1] - z[k])*(vdown - vup);
}


double  ScalarField3D_WRF::Value(double X, double Y, double Z, Neighbours &vicini)
{    
  double vup, vdown;
  VecDoub p(2);
  Vec3<double> pt;
  
  pt.setX(X);
  pt.setY(Y);
  pt.setZ(0.); // serve per vicini.Init(....)

  p[0] = X;
  p[1] = Y;
  
  /* i campi li estendo a eta = 1 ponendo il valore del campo = 0 (U,V,W) */
  if(Z > 1. || Z < z[nz-1])
    {
      cout << "Valore oltre il range z in Field3D.value(x,y,z, vicini): Fields3D.C" << endl;
      cout << Z << " " << z[0]  << endl;
      global.GestoreErrori("Valore oltre il range z in Field3D.value(x,y,z): Fields3D.C");
      throw int(13);
    }
  
  size_t k;
  k=0;
  
  while(z[k] >= Z && k < nz - 1) k++;
  
  /* i campi li estendo a eta = 1 ponendo il valore del campo = 0 (U,V,W) */
  
  Shep_interp xx_up(PtiGriglia[k], data[k], degree);
  
  long vici;      

  if(k==0) /* quindi Z > sigma[0]*/
    {
      vup = xx_up.interp(p, vicini);
      
      if(vicini.checkNews()) 
	vicini.Init(PtiGriglia[k], pt);
      
      return vup * (1. - Z) / (1 -  z[0]);
    }

  Shep_interp xx_down(PtiGriglia[k-1], data[k-1], degree);
  
  if(z[k] == Z) 
    {
      vup = xx_up.interp(p, vicini); 
      if(vicini.checkNews()) 
	vicini.Init(PtiGriglia[k], pt);
      return vup;
    }
  
  if(z[k-1] == Z)
    {
      vup =  xx_down.interp(p, vicini);
      if(vicini.checkNews()) 
	vicini.Init(PtiGriglia[k-1], pt);
      return vup;
    }
  
  vup = xx_up.interp(p, vicini);
  vdown = xx_down.interp(p, vicini);
  
  if(vicini.checkNews()) 
    vicini.Init(PtiGriglia[k], pt);
  
    return vup + (Z - z[k])/(z[k-1] - z[k])*(vdown - vup);
}



/* %%%%%%%%%%%% QUI INIZIA LA MASK %%%%%%%%%%%%%%%%%%%%%%%%%%*/

class ScalarField3D_Mask_WRF : public ScalarField3D_Mask, public ScalarField3D_WRF{

public:

  size_t nx_orography, ny_orography;
  bool orography;
  VecDoub Orografia;

  ScalarField3D_Mask_WRF() : ScalarField3D_Mask(), ScalarField3D_WRF(){} 

  void Init(const char netcdf_file[], const char nome_lon[], const char nome_lat[], \
       const char nome_depth[], const char nomeExtension[])  
  {
    ScalarField3D_WRF::Init(netcdf_file, nome_lon, nome_lat, nome_depth, nomeExtension);
  }

  void CloseDivergence(ScalarField3D& U, ScalarField3D& V){};

  void Load(const char netcdf_file[], const char nome_var[], size_t T, const char nomeExtension[])
  {
    cout << "Mask non permetta il Load standard, utilizzare ScalarField3D_Mask_[type]::Mask(....)\n\n";
    throw int(90);
  }


  double Value(double X, double Y, double Z) 
  {
    return ScalarField3D_WRF::Value(X, Y, Z);
  }

  double Value(double X, double Y, double Z, Neighbours &vicini)
  {
    return ScalarField3D_WRF::Value(X, Y, Z, vicini);
  }

  MatDoub& Points4Neighbour(bool& isUsable) 
  {
    return ScalarField3D_WRF::Points4Neighbour(isUsable);
  }

  size_t NTimes(const char netcdf_file[])
  {
    return ScalarField3D_WRF::NTimes(netcdf_file);
  }

  void LoadOrografia(const char netcdf_file[], const char nome_var[])
  {
    long int i,j,k, ToT;
    double v, missingvalue;
    int ncid, var_id, long_id, lat_id, z_id, Time_id, x_id;
    int retval, status;
    nc_type type;   
    int ndims, natts, ndims2, ndims3;
    size_t ntimes;
    int  dimids[NC_MAX_VAR_DIMS]; 
    size_t *count, *count2, *start;
    ptrdiff_t *stride;
    
    
    count =  (size_t *)malloc(sizeof(size_t)*3);
    start = (size_t *)malloc(sizeof(size_t)*3);
    stride = (ptrdiff_t *)malloc(sizeof(ptrdiff_t)*3);

    status =  nc_open(netcdf_file, NC_SHARE, &ncid);
    if(status) { cout << "File NetCDF non trovato \n"; throw int(1);}

    status = nc_inq_varid(ncid, nome_var, &x_id);
    if(status) { cout << "Variabile " << nome_var << " non trovata\n"; throw int(3);}

    status = nc_inq_var(ncid, x_id, 0, &type, &ndims, dimids, &natts);
    if(status) { cout << "Variabile " << nome_var << " non trovata\n"; throw int(3);}

    status = nc_inq_dimlen(ncid, dimids[1], &ny_orography);
    status = nc_inq_dimlen(ncid, dimids[2], &nx_orography);

    cout << "Orografia " << nx_orography << " x " << ny_orography << endl; 

    global.log << "# Loading orography from: " << netcdf_file << " :  " << nome_var << endl;  
    cout << "# Loading orography from: " << netcdf_file << " :  " << nome_var << endl;  

    start[0]=0; start[1]=0; start[2]=0;

    count[0]= 1; count[1]=ny_orography; count[2]= nx_orography;
    stride[0] = stride[1] = stride[2] = 1;

    double *f = new double[nx_orography * ny_orography];

    status = nc_get_vars_double(ncid, x_id, start, count, stride, f);
    
    if(status)
      {
	global.er << "Probabili problemi nel caricamento dati \n" << endl;
	cout << "Probabili problemi nel caricamento dati \n" << endl;
	throw int(21);
      }

    Orografia.resize(nx_orography * ny_orography);

    for(i=0; i<nx_orography; i++) 
      for(j=0; j<ny_orography; j++) 
	{
	  Orografia[i + nx_orography * j] = *(f + i + nx_orography *j);    
 	}

    delete f;

    free(count);
    free(start);
    free(stride);

    nc_close(ncid);
    
  }


void Mask(const char netcdf_file[], const char nome_var[])
  {
    long int i,j,k;
    double v;
    int ncid, x_id;
    int retval, status;
    nc_type type;   
    int ndims, natts, ndims2, ndims3;
    size_t ntimes;
    int  dimids[NC_MAX_VAR_DIMS]; 
    size_t *count,  *start;
    ptrdiff_t *stride;
    
    count =  (size_t *)malloc(sizeof(size_t)*4);
    start = (size_t *)malloc(sizeof(size_t)*4);
    stride = (ptrdiff_t *)malloc(sizeof(ptrdiff_t)*4);

    status =  nc_open(netcdf_file, NC_SHARE, &ncid);
    if(status) { cout << "File NetCDF non trovato \n"; throw int(1);}

    status = nc_inq_varid(ncid, nome_var, &x_id);
    if(status) { cout << "Variabile " << nome_var << " non trovata\n"; throw int(3);}

    status = nc_inq_var(ncid, x_id, 0, &type, &ndims, dimids, &natts);
    if(status) { cout << "Variabile " << nome_var << " non trovata\n"; throw int(3);}

    /*    status = nc_inq_dimlen(ncid, dimids[3], &nx);
    cout << ntimes << " " << Nx() << endl;

    status = nc_inq_dimlen(ncid, dimids[2], &ny);
    cout << ntimes << " " << Ny() << endl;

    status = nc_inq_dimlen(ncid, dimids[1], &nz);
    cout << nz << " " << endl; */

    status = nc_inq_dimlen(ncid, dimids[0], &ntimes);
    cout << ntimes << endl; 

    ScalarField3D_WRF::ntime = ntimes;

    if(status) { cout << "Dimensione variabile " << nome_var << " non trovata\n"; throw int(3);}


    global.log << "# Loading Mask from: " << netcdf_file << " :  " << nome_var << endl;  
    cout << "# Loading Mask from: " << netcdf_file << " :  " << nome_var << endl;  

   
    start[0]=0;
    start[1]=start[2]=start[3]=0;

    count[0]=1; count[1]= ScalarField3D_WRF::Nz(); count[2]= Ny(); count[3]= Nx();
    //count[0]=1; count[1]= Nz(); count[2]= Ny(); count[3]= 10;
    stride[0] = stride[1] = stride[2] = stride[3] = 1;     

    cout << count[1] << endl;
    cout << count[2] << endl;
    cout << count[3] << endl;

    size_t n_punti_bulk;

    n_punti_bulk = ScalarField3D_WRF::Nz() * Ny() * Nx();
      
    cout << n_punti_bulk << endl;

    float *ff;

    ff = new float[n_punti_bulk];

    status = 0;

    status = nc_get_vars_float(ncid, x_id, start, count, stride, ff);
    
    if(status)
      {
	global.er << "Probabili problemi nel caricamento dati \n" << endl;
	cout << "Probabili problemi nel caricamento dati in Mask \n" << endl;
	throw int(21);
      }

    size_t level, nxlevel;

    nxlevel = Nx() * Ny();
 
    for(k=0; k<n_punti_bulk; k++)
      { 
	if(ff[k] > 1.e6)  ff[k]=0.;

	level = k / nxlevel;
	i = k % nxlevel;
	(*(data+level))[i] = 1.;
      }


    delete ff;

    free(count);
    free(start);
    free(stride);

    nc_close(ncid);
  }


  double toplevel(Vec3<double>& p, Neighbours &vicini)
  { 
    /*    if(global.ocean && !global.sigma) return ScalarField3D_WRF::Sigma(0); 
    if(!global.sigma) return ScalarField3D_WRF::sigma[ScalarField3D_WRF::Nz() -1];
    if(global.ocean && global.sigma) return 0.;
    if(!global.ocean && global.sigma) return 12000.; */ // TO CHECK

    return 20000.;
  } 

  double toplevel(Vec3<double>& p)
  { 
    Neighbours dummy;
    return toplevel(p,dummy);
  } 

  double bottom(Vec3<double>& p)
  {
    Neighbours dummy(16);

    dummy.Dummy(Nx() * Ny() );

    return bottom(p,dummy);
  }


  double bottom(Vec3<double>& p, Neighbours& vicini)
  { 
    double z1, z2, zm, v1, v2, vm; 

    VecDoub pp(2);
    Shep_interp xx(PtiGriglia[0], Orografia, 2);

    pp[0] = p.X();  pp[1] = p.Y();
    
    return xx.interp(pp, vicini);	
  }


  bool Inside(Vec3<double> &p)
  {
    if(p.X() < x_min || p.X() > x_max) return false;
    if(p.Y() < y_min || p.Y() > y_max) return false;
    
    return true;
  }

}; 
