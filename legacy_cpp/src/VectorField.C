


/*****************
void Wdivergence(ScalarField3D& U, ScalarField3D& V)
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
	      
	      cout << z[k] << " " << wz << " " << endl;
	      //              cout.flush();
	    }
	}
    }
}********/









class VectorField3D
{
  
public:
  ScalarField3D *X,*Y,*Z;
  bool zflag;

  VectorField3D(){}

  VectorField3D(ScalarField3D* x, ScalarField3D* y, ScalarField3D* z)
  {
    X = x; Y = y; Z = z;
    zflag = true;
  }

  VectorField3D(ScalarField3D* x, ScalarField3D* y)
  {
    X = x; Y = y; 
    zflag = false;
  }

  VectorField3D(KindOfFile_t type)
  {
    /* da lanciare come prima cosa in modo che poi partano gli Init e Load giusti */
    
    X = Switcher(type);
    Y = Switcher(type);

    if(Lista.IsThereW)
      {
	Z = Switcher(type);
	zflag=true;
      }
    else zflag=false;
  }

  void Init(ListOfFile& Lista)
  {
    X->Init(Lista.Nome(0), Lista.nomelonU, Lista.nomelatU, Lista.nomezU, Lista.nomeextensionU);
    Y->Init(Lista.Nome(0), Lista.nomelonV, Lista.nomelatV, Lista.nomezV, Lista.nomeextensionV);

    if(zflag)
      Z->Init(Lista.Nome(0), Lista.nomelonW, Lista.nomelatW, Lista.nomezW, Lista.nomeextensionW);   
  }

  VectorField3D(ListOfFile& Lista){  Init(Lista); }


  

  void operator ^ (VectorField3D &B)
  {
    ScalarField3D *temp;
    /* operatore di swap senza lavoro, scambia i puntatori*/
    temp = X; X = B.X; B.X = temp;
    temp = Y; Y = B.Y; B.Y = temp;

    if(zflag) 
      {
	temp = Z; Z = B.Z; B.Z = temp;
      }
  } 

  Vec3<double> Value(size_t i, size_t j, size_t k)
  {
    Vec3<double> out;
    
    if(zflag) out.set(X->Value(i,j,k), Y->Value(i,j,k), Z->Value(i,j,k) );
    else out.set(X->Value(i,j,k), Y->Value(i,j,k), 0. );
    
    return out;
    return out;
  }

  Vec3<double> Value(Vec3<double>& P)
  {
    Vec3<double> out;

    try{
    
    if(zflag) out.set(X->Value(P), Y->Value(P), Z->Value(P) );
    else out.set(X->Value(P), Y->Value(P), 0. );
    }
    catch( int ercode ){
      cout << "Errore lanciato da VectorField.Value(x,y,z) " << ercode << endl;
      throw int(ercode);
    }

    return out;
  }

  Vec3<double> Value(Vec3<double>& P, Neighbours& vicinix, Neighbours& viciniy, Neighbours& viciniz)
  {
    Vec3<double> out;
    
    if(zflag) out.set(X->Value(P, vicinix), Y->Value(P, viciniy), Z->Value(P, viciniz) );
    else out.set(X->Value(P, vicinix), Y->Value(P, viciniy), 0. );

    return out;
  }

  double Time(){ return X->Time(); }

  void Time(double t)
  { 
    X->Time(t);  Y->Time(t); 
    if(zflag)
      Z->Time(t); 
  }

};
