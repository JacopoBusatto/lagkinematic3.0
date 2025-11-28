
struct RBF_fn {
	virtual Doub rbf(Doub r) = 0;
};

struct RBF_interp {
	Int dim, n;
	const MatDoub &pts;
	const VecDoub &vals;
	VecDoub w;
	RBF_fn &fn;
	Bool norm;

	RBF_interp(MatDoub_I &ptss, VecDoub_I &valss, RBF_fn &func, Bool nrbf=false)
	: dim(ptss.ncols()), n(ptss.nrows()) , pts(ptss), vals(valss),
	w(n), fn(func), norm(nrbf) {
		Int i,j;
		Doub sum;
		MatDoub rbf(n,n);
		VecDoub rhs(n);
		for (i=0;i<n;i++) {
		  sum = 0.;
   		  for (j=0;j<n;j++) {
		    sum += (rbf[i][j] = fn.rbf(palatella_math_functions::CerchioMassimo(&pts[i][0],&pts[j][0])));
		  }
		  
		  if (norm) rhs[i] = sum*vals[i];
		  else rhs[i] = vals[i];
		}
		LUdcmp lu(rbf);
		lu.solve(rhs,w);
	}

	Doub interp(VecDoub_I &pt) {
		Doub fval, sum=0., sumw=0.;
		if (pt.size() != dim) throw("RBF_interp bad pt size");
		for (Int i=0;i<n;i++) {
			fval = fn.rbf(palatella_math_functions::CerchioMassimo(&pt[0],&pts[i][0]));
			sumw += w[i]*fval;
			sum += fval;
		}
		return norm ? sumw/sum : sumw;
	}

	Doub rad(const Doub *p1, const Doub *p2) {
		Doub sum = 0.;
		for (Int i=0;i<dim;i++) sum += SQR(p1[i]-p2[i]);
		return sqrt(sum);
	}
};
struct RBF_multiquadric : RBF_fn {
	Doub r02;
	RBF_multiquadric(Doub scale=1.) : r02(SQR(scale)) {}
	Doub rbf(Doub r) { return sqrt(SQR(r)+r02); }
};

struct RBF_thinplate : RBF_fn {
	Doub r0;
	RBF_thinplate(Doub scale=1.) : r0(scale) {}
	Doub rbf(Doub r) { return r <= 0. ? 0. : SQR(r)*log(r/r0); }
};

struct RBF_gauss : RBF_fn {
	Doub r0;
	RBF_gauss(Doub scale=1.) : r0(scale) {}
	Doub rbf(Doub r) { return exp(-0.5*SQR(r/r0)); }
};

struct RBF_inversemultiquadric : RBF_fn {
	Doub r02;
	RBF_inversemultiquadric(Doub scale=1.) : r02(SQR(scale)) {}
	Doub rbf(Doub r) { return 1./sqrt(SQR(r)+r02); }
};
struct Shep_interp {
  Int dim, n;
  const MatDoub &pts;
  const VecDoub &vals;
  Doub pneg;

Shep_interp(MatDoub_I &ptss, VecDoub_I &valss, Doub p)
: dim(ptss.ncols()), n(ptss.nrows()) , pts(ptss),
    vals(valss), pneg(-p) {}
  
  Doub interp(VecDoub_I &pt) {
    Doub r, w, sum=0., sumw=0.;
    if (pt.size() != dim) 
      {
	cout << pt.size() << " " << dim << endl;
	throw("RBF_interp bad pt size");
      }
    for (Int i=0;i<n;i++) {
      if ((r=palatella_math_functions::CerchioMassimo(&pt[0],&pts[i][0])) > 1.e4) return vals[i]; 
      if(r > 0.)
	{
	  sum += (w = pow(r,-pneg));
	  sumw += w*vals[i];
	}
      //	if(w > 1.e-4)
      //	  cout <<"gasp " << pt[0] << " " << pt[1]		\
      //    << " " << pt[2] << " " << w << " " << vals[i] << endl; 
    }
    return sumw/sum;
  }
  
  /*! Interpolate only on a sub-set of grid points, pointed by iv[]
   */
  Doub interp(VecDoub_I &pt, Neighbours &vic){
    Doub r, w, sum=0., sumw=0.;
    Int ii;
    Doub first_dist=-1;
    Int doCheck=0;
    if (pt.size() != dim) 
      {
	cout << pt.size() << " " << dim << endl;
	global.er << " RBF_interp error  "<< pt.size() << " " << dim << endl;
	throw("RBF_interp bad pt size");
      }

    for (Int i=0;i<vic.ind.size();i++) 
      {
	ii=vic.ind[i];

	if(ii >= pts.nrows()  || ii < 0)
	  {
	    global.er << "indice vicini fuori range!\n\n";
	    throw int(1357);
	  }

	if ((r=palatella_math_functions::CerchioMassimo(&pt[0],&pts[ii][0])) > 1.e4) return vals[ii];   

	/* r è l'inverso della distanza ! , quindi equivale a dire r prossimo a 0, 
	   quindi se il punto è sostanzialmente coincidente con uno della griglia, 
	   resituisci quel valore, senza interpolare */

	if (i==0) first_dist=r;  // il primo punto , quello più vicino ( e quindi con r maggiore  ! )
	
	else{
	  if(r>first_dist) doCheck=1;    // Il vettore dei vicini va ridefinito
	}
	if(r > 0.)
	  {
	    sum += (w = pow(r,-pneg));
	    sumw += w*vals[ii];
	  }

	//	cout << i << " " << ii << " " << r << " " <<  
      }
    vic.doCheck=doCheck;
    return sumw/sum;
  }
  
  
  Doub rad(const Doub *p1, const Doub *p2) {
    Doub sum = 0.;
    for (Int i=0;i<dim;i++) sum += SQR(p1[i]-p2[i]);
    return sqrt(sum);
  }
};


