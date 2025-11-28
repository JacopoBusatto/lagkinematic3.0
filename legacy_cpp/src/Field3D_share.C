//#include <netcdf.h>

/* #include <nr3.h>
#include <interp_1d.h>
#include <interp_linear.h>
#include <interp_2d.h>
#include "Vector.C" */

using namespace std;

//static const int NC_ERR2 = 2;

/****** CODICI ERRORI NETCDF FATTI IN CASA *******/

/*
  0  file non trovato
  1  dimensione non trovata
  2  variabile non trovata
  3  lettura dati
  4  posizionamento indice

  7 dimensione longitudine > 1
 */


/***********************************************/

void errore(int code){
  char codice_errore[100];

  sprintf(codice_errore,"Errore gestione file NetCDF: code %d",code);
  global.GestoreErrori(codice_errore);
}

class ScalarField3D
{
private:
  VecDoub Lon, Lat, z; // le coordinate dei punti griglia
  size_t nx, ny, nz;
  size_t mx, my;
  double time;
  bool filled_long_lat, mask;

public:
  MatDoub *data;
  // ------------ Constructors ------------
  
  ScalarField3D(){mx = my = 15;}

  ScalarField3D(size_t n, size_t m, size_t k) : Lon(n), Lat(m), z(k)
  {
    size_t i;
    data = new MatDoub[k];

    for(i=0; i<k; i++) (data + i)->resize(n,m);

    nx = n; ny = m; nz = k;

    filled_long_lat = false;
    mask = false;
    mx = my = 15;
  }
  
  // Three parameter constructor, non esiste il copy contructor, da
  // scrivere in modo che la cosa venga evitata!

  ScalarField3D(size_t Nx, size_t Ny, size_t Nz, double *x, double *y, double*zz, double* f) : Lon(Nx), Lat(Ny), z(Nz)
  {
    Int i,j,k;
    Doub v;
    
    mask = false;
    nx = Nx; ny = Ny; nz = Nz; 
    mx = my = 15;

    data = new MatDoub[nz];
    for(i=0; i<nz; i++) (data + i)->resize(nx,ny);
    
    for(i=0; i<nx; i++) Lon[i] = x[i];
    for(j=0; j<ny; j++) Lat[j] = y[i];
    for(k=0; k<nz; k++) z[j] = zz[i];

    for(k=0; k<nz; k++)
      {
	for(j=0; j<ny; j++) 
	  {
	    for(i=0; i<nx; i++) 
	      {
		v = *(f + i + nx*j + nx*ny*k);
		(*(data+k))[i][j] = v;
	      }
	  }
      }
    filled_long_lat = true;
  }

  ScalarField3D(const char netcdf_file[], const char nome_lon[], const char nome_lat[], const char nome_depth[], int dimlon)
  {
    Init(netcdf_file,nome_lon, nome_lat, nome_depth, dimlon);
  }

  inline void Init(const char netcdf_file[], const char nome_lon[], const char nome_lat[], const char nome_depth[], int dimlon)
  {
    /*nome_var serve per calcolare la mask*/
    long int i,j,k;
    int ncid, var_id, long_id, lat_id, z_id, Time_id;
    int retval, status;
    
    mx = my = 15;
    mask = false;

    if ((retval = nc_open(netcdf_file, NC_SHARE, &ncid))) errore(0);

    status = nc_inq_varid(ncid, nome_lon, &long_id);
    if(status) errore(2);
    status = nc_inq_varid(ncid, nome_lat, &lat_id);
    if(status) errore(2);
    status = nc_inq_varid(ncid, nome_depth, &z_id);
    if(status) errore(2);

    double *lats, *lons, *zz, v;
    nc_type type;   
    int ndims, natts, ndims2, ndims3;
    int  dimids[NC_MAX_VAR_DIMS]; 
    size_t *count, *count2, lenghtp, *index;

    count = (size_t*)malloc(sizeof(size_t)*ndims);
    count2 = (size_t*)malloc(sizeof(size_t)*2);

    status = nc_inq_var(ncid, long_id, 0, &type, &ndims, dimids, &natts);

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
 
    status = nc_inq_var(ncid, z_id, 0, &type, &ndims3, dimids, &natts);
    nc_inq_dimlen(ncid, dimids[0], &nz);

    z.resize(nz);

    for(size_t i=0; i<nz; i++)
      {
	status = nc_get_var1_double(ncid, z_id, &i, &v);
	z[i] = -fabs(v);
      }

    fprintf(stdout,"# griglia: %ld %ld %ld\n",nx, ny, nz);
    global.log << "# dimensioni griglia: nx x ny x nz: " << nx << " " <<  ny << " " << nz << endl;

    data = new MatDoub[nz];

    for(i=0; i<nz; i++) (data + i)->resize(nx,ny);

    filled_long_lat = true;

    fprintf(stdout,"# Completato caricamento caratteristiche griglia\n");
    fflush(stdout);

    //    for(j=0; j<nx; j++) cout <<  "lon " << Lon[j] << "\n ";
    //   for(j=0; j<ny; j++) cout <<  Lat[j] << "\n ";
    nc_close(ncid);
  }

  inline void Init(ListOfFile& Lista, char component)
  {
    /* di default prende la griglia della U*/
    if(component=='U')
      Init(Lista.Nome(0), Lista.nomelonU,Lista.nomelatU,Lista.nomezU,Lista.ndimlon);
    if(component=='V')
      Init(Lista.Nome(0), Lista.nomelonV,Lista.nomelatV,Lista.nomezV,Lista.ndimlon);
    if(component=='W')
      Init(Lista.Nome(0), Lista.nomelonW,Lista.nomelatW,Lista.nomezW,Lista.ndimlon);
  }

  void Load(double *f)
  {
    size_t k, i, j;
    double v;

    for(k=0; k<nz; k++)
      { 
	for(i=0; i<nx; i++) 
	  {
	    for(j=0; j<ny; j++) 
	      {
		v = *(f + i + nx*j + nx*ny*k);
		if(v > 1.e6)  (*(data+k))[i][j] = 0.;
		else (*(data+k))[i][j] = v;
	      }
	  }
      }
  }



  void Load(const char netcdf_file[], const char nome_var[], size_t T)
  {
    double *f = new double[nz* ny * nx];
    Load(netcdf_file, nome_var, T, f);
    free(f);
  }

  void Load(const char netcdf_file[], const char nome_var[], size_t T, double *f)
  {
    if(!filled_long_lat)
      {
	global.GestoreErrori("File non pronto con lat e lon");
	return;
      }

    long int i,j,k, ToT;
    double v;
    int ncid, var_id, long_id, lat_id, z_id, Time_id, x_id;
    int retval, status;
    nc_type type;   
    int ndims, natts, ndims2, ndims3;
    size_t ntimes;
    int  dimids[NC_MAX_VAR_DIMS]; 
    size_t *count, *count2, *start;
    ptrdiff_t *stride;
    
    count =  (size_t *)malloc(sizeof(size_t)*4);
    start = (size_t *)malloc(sizeof(size_t)*4);
    stride = (ptrdiff_t *)malloc(sizeof(ptrdiff_t)*4);

    if ((retval = nc_open(netcdf_file, NC_SHARE, &ncid))) errore(0);

    status = nc_inq_varid(ncid, nome_var, &x_id);
    if(status) errore(2);

    status = nc_inq_var(ncid, x_id, 0, &type, &ndims, dimids, &natts);
    if(status) errore(2);

    status = nc_inq_dimlen(ncid, dimids[0], &ntimes);
    if(status) errore(2);

    if(global.backtraj) ToT = ntimes - T - 1;
    else ToT = T;

    if(ntimes <= ToT)
      {
	global.GestoreErrori("Indice tempo T out of range: Fields3D.C");
	return;
      }

    global.log << "# Loading: " << netcdf_file << " :  " << nome_var << " : " << ToT << endl;  

    start[0]=ToT;
    start[1]=start[2]=start[3]=0;

    count[0]=1; count[1]= nz; count[2]= ny; count[3]= nx;
    stride[0] = stride[1] = stride[2] = stride[3] = 1;     

    status = nc_get_vars_double(ncid, x_id, start, count, stride, f);
    
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
		v = *(f + i + nx*j + nx*ny*k);
		if(v > 1.e6)  (*(data+k))[i][j] = 0.;
		else (*(data+k))[i][j] = v;
	      }
	  }
      }

    nc_close(ncid);
  }

  void Mask(const char netcdf_file[], const char nome_var[])
  {
    mx = my = 3;

    mask = true;

    if(!filled_long_lat)
      {
	global.GestoreErrori("File non pronto con lat e lon: Mask");
	return;
      }

    double v;
    long int i,j,k, ToT;
    int ncid, var_id, long_id, lat_id, z_id, Time_id, x_id;
    int retval, status;
    nc_type type;   
    int ndims, natts, ndims2, ndims3, ntimes;
    int  dimids[NC_MAX_VAR_DIMS]; 
    size_t *count, *count2, *start;
    ptrdiff_t *stride;
    
    count =  (size_t *)malloc(sizeof(size_t)*4);
    start = (size_t *)malloc(sizeof(size_t)*4);
    stride = (ptrdiff_t *)malloc(sizeof(ptrdiff_t)*4);

    if ((retval = nc_open(netcdf_file, NC_SHARE, &ncid))) errore(0);

    status = nc_inq_varid(ncid, nome_var, &x_id);
    if(status) errore(2);

    start[0]=0;
    start[1]=start[2]=start[3]=0;

    count[0]=1; count[1]= nz; count[2]= ny; count[3]= nx;
    stride[0] = stride[1] = stride[2] = stride[3] = 1;     

    double *f = new double[nz* ny * nx];

    status = nc_get_vars_double(ncid, x_id, start, count, stride, f);

    for(k=0; k<nz; k++)
      { 
	for(i=0; i<nx; i++) 
	  {
	    for(j=0; j<ny; j++) 
	      {
		v = *(f + i + nx*j + nx*ny*k);
		if(v > 1.e6)  (*(data+k))[i][j] = 0.;
		else (*(data+k))[i][j] = 1.;
	      }
	  }
      }
    nc_close(ncid);
  }

  double profondita(Vec3<double>& p)
  {
    /* solo se e' un campo mask!!! pensare ad un flag mask*/
    double zp;

    if(!mask) 
      {
	global.GestoreErrori("Campo non mask, impossibile calcolare la profondita");
	return -15000.; // in caso mettessi il rimbalzo
      }
    
    zp = z[0];
    
    while(Value(p.X(), p.Y(), zp) > 0.1 && p.Z() > z[nz -1] ) zp -= 0.5;
    return zp;
  }

  void DummyMask()
  {
    if(!mask)
      {
	global.GestoreErrori("Campo non mask, impossibile eseguire DummyMask");
	return ;
      }

    global.param << endl;

    size_t i,j;

    for(j=0; j<ny; j++) 
      {
	for(i=0; i<nx; i++) 
	  {
	    global.param << i << " " << j << " " << Value(i,j,0) << endl;
	  }
	global.param << endl;
      }
    
  }



  void ReboundZ(Vec3<double>& p)
  {
    const double e = 1.e-5;

    if(p.Z() >= z[0]) p.setZ(2. * z[0] - p.Z() - e);
    if(p.Z() <= z[nz - 1]) p.setZ(2. * z[nz -1] - p.Z() + e);
  }

  bool Inside(Vec3<double>& p)
  {
    if(p.X() < Lon[0] || p.X() > Lon[nx-1]) return false;
    if(p.Y() < Lat[0] || p.Y() > Lat[ny-1]) return false;
    
    return true;
  }


  double Time(){ return time;}
  void Time(double t){  time = t;}

  void operator ^ (ScalarField3D &B)
  {
    /* operatore di swap senza lavoro, scambia i puntatori*/
    MatDoub *a,*temp;

    temp = data;
    data = B.data;
    B.data = temp;
  } 

  Doub Value(size_t i, size_t j, size_t k)
  {
    if(i>=0 && i<nx && j>=0 && j<ny && k >= 0 && k < nz) return (*(data+k))[i][j];
    else 
      {
	//	cout << i << " " <<  j <<" " << k << " " << nx << " " << ny << " " << nz << endl; 
	global.GestoreErrori("Valore oltre il range 0 - nx-1 in Field3D.value(i,j,k): Fields3D.C");	
	return 0.;
      }
  }

  Doub Value(double X, double Y, double Z) 
  {
    size_t k, k1;
    Doub v1, v2;

    k=0;
    
    if(Z >= z[0] || Z < z[nz - 1] || X < Lon[0] || X > Lon[nx-1] || Y < Lat[0] || Y > Lat[ny-1]) 
      {
	global.GestoreErrori("Valore oltre il range 0 - nx-1 in Field3D.value(x,y,z): Fields3D.C");
	throw int(13);
	//	return 0.;
      }

    while(z[k] >= Z && k < nz) k++;

    Poly2D_interp suz1(Lon, Lat, *(data + k), mx, my), suz2(Lon, Lat, *(data + k - 1), mx, my);

    v1 = suz1.interp(X,Y);
    v2 = suz2.interp(X,Y);

    return v1 + (Z - z[k])/(z[k-1] - z[k])*(v2 - v1);
  }


  Doub Value(Vec3<double> X)
  {
    return Value(X.X(), X.Y(), X.Z());
  }

  Vec3<double> Gradient(Vec3<double> X)
  {
    Vec3<double> g, dx(1,0,0), dy(0,1,0);
    Doub gx, gy, h = 0.1, N;
    
    if(X.X() - h < Lon[0] || X.X() + h > Lon[nx-1])
      {
	global.GestoreErrori("Gradiente troppo vicino al bordo in longitudine\n");
	throw int(17);
      }
      
    if(X.Y() - h < Lat[0] || X.Y() + h > Lat[ny-1])
      {
	global.GestoreErrori("Gradiente troppo vicino al bordo in latitudine\n");
	throw int(18);
      }
    
    gx = 0.5 / h * (Value(X + dx * h) - Value(X + dx * (-h)));
    gy = 0.5 / h * (Value(X + dy * h) - Value(X + dy * (-h)));
    
    //    cout << gx << " " << gy << endl;

    //    N = sqrt(gx * gx + gy * gy);

    g.setX(gx); g.setY(gy);
    g.setZ(0.);

    return g;
  }

  size_t Nz(){ return nz;}
  size_t Ny(){ return ny;}
  size_t Nx(){ return nx;}

  size_t Nxslab()
  {
    return nx * ny;
  }

  size_t Ntot(){ return nx * ny * nz;}

  void CloseDivergence(ScalarField3D& U, ScalarField3D& V)
  {
    size_t k, i, j;
    double dwdz, wz, hx, hy, hz;
    
    for(i=1; i<Nx(); i++) 
      { 
	for(j=1; j<Ny(); j++) 
	  {
	    wz = 0;
	    hx = 1./(U.Lon[i] - U.Lon[i-1]) / 111000. / cos (6.28318530717959 * Lat[j] / 360.) ;
	    hy = 1./(V.Lat[j] - V.Lat[j-1]) / 111000.;

	    (*(data+nz-1))[i][j] = 0.;

	    for(k=nz-1; k>0; k--)
	      {
		hz = z[k] - z[k - 1];
		wz -= 0.5 * hz * hx * ( U.Value(i,j,k) - U.Value(i-1,j,k) );
		wz -= 0.5 * hz * hy * ( V.Value(i,j,k) - V.Value(i,j-1,k) );
		(*(data+k-1))[i][j] = wz;

		//	        cout << hx << " " << hy << " " << hz << " " << wz << endl;
		//		cout.flush();
	      }
	  }
      }
  }


  void PrintOut(Vec3<double> X1, Vec3<double> X2, const char nomefile[])
  {
    ofstream out;
    double x,y,z;

    out.open(nomefile);

    for(z = X1.Z(); z < X2.Z(); z += 1)
      {
	for(y = X1.Y(); y < X2.Y(); y += 0.1)
	  {
	    for(x = X1.X(); x < X2.X(); x += 0.1)
	      out << x << " " << y << " " << z  << " " << this->Value(x,y,z) << endl; 
	    out << endl; 
	  }
	out << endl;
      }
    out.close();
  }

  friend class VectorField3D;

};


/* ================================================================= */

class VectorField3D
{
  
public:
  ScalarField3D X,Y,Z;
  bool zflag;

  VectorField3D(){}

  VectorField3D(ScalarField3D& x, ScalarField3D& y, ScalarField3D& z)
  {
    X = x; Y = y; Z = z;
    zflag = true;
  }

  VectorField3D(ScalarField3D& x, ScalarField3D& y)
  {
    X = x; Y = y; 
    zflag = false;
  }

  void Init(ListOfFile& Lista)
  {
    X.Init(Lista,'U');
    Y.Init(Lista,'V');
    if(Lista.IsThereW)
      Z.Init(Lista,'W');
    else
      {
	Z.Init(Lista.Nome(0), Lista.nomelonV,Lista.nomelatU,Lista.nomezU,Lista.ndimlon);
      }
    zflag=true;
  }

  VectorField3D(ListOfFile& Lista){  Init(Lista); }

  void operator ^ (VectorField3D &B)
  {
    /* operatore di swap senza lavoro, scambia i puntatori*/
    X ^ B.X;
    Y ^ B.Y;
    if(zflag) Z ^ B.Z;
  } 

  Vec3<double> Value(size_t i, size_t j, size_t k)
  {
    Vec3<double> out;
    
    if(zflag) out.set(X.Value(i,j,k), Y.Value(i,j,k), Z.Value(i,j,k) );
    else out.set(X.Value(i,j,k), Y.Value(i,j,k), 0. );
    
    return out;
  }

  Vec3<double> Value(Vec3<double>& P)
  {
    Vec3<double> out;
    
    if(zflag) out.set(X.Value(P), Y.Value(P), Z.Value(P) );
    else out.set(X.Value(P), Y.Value(P), 0. );

    return out;
  }

  double Time(){ return X.Time(); }

  void Time(double t){ X.Time(t);  Y.Time(t); Z.Time(t); }

};



long int NFieldsXFile(const char netcdf_file[], const char nome_var[])
{
  size_t out;
  int ncid, var_id, x_id;
  int retval, status;
  nc_type type;   
  int ndims, natts, ndims2, ndims3;
  int  dimids[NC_MAX_VAR_DIMS]; 
  
  if ((retval = nc_open(netcdf_file, NC_SHARE, &ncid))) errore(0);
  
  status = nc_inq_varid(ncid, nome_var, &x_id);
  if(status) errore(2);
  
  status = nc_inq_var(ncid, x_id, 0, &type, &ndims, dimids, &natts);

  nc_inq_dimlen(ncid, dimids[0], &out);
  
  cout << "numero campi per file "<< out << endl;

  nc_close(ncid);

  return out;
}



	
