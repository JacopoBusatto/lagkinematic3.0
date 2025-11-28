#include "sort.h"

using namespace std;

/** A index array of grid points 
 *
 */
class Neighbours
{
private:
  int doCheck; 

public:
  VecInt ind;  // indexes of them  // todo diventa vettore 
  
  Neighbours(){}
  
  Neighbours(int n_vicini){
    ind.resize(n_vicini);
  }

  void Dummy(int n_all)
  {
    ind.resize(n_all);

    for (int i=0;i<ind.size();i++)
      ind[i]=i;
  }
  
  int Nvicini(){ return ind.size();}
  

  /** Calcola i vicini unicamente della riga richiesta
   *
   */
  inline void Init(const MatDoub_I &gridPts, const Vec3<double> &pts) 
  { 
    doCheck=0; // Reset the flag 0 aggiornato, 1 da aggiornare
    int nn = gridPts.nrows();
   
    double *vD = new double[nn];
    double *ptss = new double[2];

    ptss[0]=pts.X(); 
    ptss[1]=pts.Y(); 


    for (int i =0; i<nn;i++)
      {
	vD[i]=palatella_math_functions::CerchioMassimo(pts.X(), pts.Y(), gridPts[i][0], gridPts[i][1]);
	//      cout << "i_v= " << i_v << " vD[" << i << "]= " << vD[i] << endl;
      }
    
    /** Sorting !!
     */
    Indexx indexx;
    indexx.n=nn;
    indexx.index(vD, nn);
    

    /** Copio in ind[] component 0-x, 1-y, 2-z     */
    int yyy = indexx.indx.size();
    
    for (int i=0;i<ind.size();i++)
      ind[i]=indexx.indx[yyy-1-i];
    
    free(vD);
    free(ptss);
      
  }
  
  
  /*  inline void Init(const ScalarField3D& field, const Vec3<double> &pts) 
  {
    Init(field.PtiGriglia[0], pts);
    }*/

  /** It checks 
   * Controlla che il primo valore dell'array sia il più vicino, fra questi, al punto
   *@param iv : indice di ind[][] su cui effettuare la ricerca
   */
  
  bool checkNews(){
    return doCheck;
  }
  
  //dbg
  /*
  int nC(){
    return ind.ncols();
  }
  */

  friend class Shep_interp;
  
};


