#include "crpropa/magneticField/FITSMagneticField.h"

namespace crpropa {
  FITSMagneticField::FITSMagneticField(const std::string& fitsFileName){

    // Where reading these parameters?
    long _fNx = 1;
    long _fNy = 1;
    long _fNz = 1;
    double _fBmax;
    double _fpB[10000];
    //
    
    int lStatus = 0 ;
    int lAnyNull = 0 ;
    int lNaxis = 0 ;
    long *lNaxes = new long[3] ; 
    
    double *lpBx = new double[_fNx*_fNy*_fNz] ; 
    fitsfile *lpBxFits ;
	
    const char *lExtBx = (fitsFileName).c_str() ;
    
    fits_open_file(&lpBxFits, lExtBx, READONLY, &lStatus) ;
    int lDummy = 0 ;
    fits_movabs_hdu(lpBxFits,1,&lDummy,&lStatus);
    fits_get_img_dim(lpBxFits, &lNaxis, &lStatus) ;
    //////if (lNaxis != 3) throw TFitsErr("NAXIS incorrect for Bx.") ;
    fits_get_img_size(lpBxFits, 3, lNaxes, &lStatus) ;
    ///if (lNaxes[0]!=_fNx || lNaxes[1]!=_fNy || lNaxes[2]!=_fNz) 
    ///  throw TFitsErr("size of Bx array not ok with xml specifications.") ;
    fits_read_img(lpBxFits, TDOUBLE, 1, _fNx*_fNy*_fNz,
		  NULL, lpBx, &lAnyNull, &lStatus) ;
    fits_close_file(lpBxFits,&lStatus) ;
    ////if (lStatus) throw TFitsErr(lStatus) ;
    ////if (lAnyNull) throw TFitsErr("Undefined elements in Bx extension.") ;
    
    double *lpBy = new double[_fNx*_fNy*_fNz] ; 
    fitsfile *lpByFits ;
    const char *lExtBy = (fitsFileName).c_str() ;
    fits_open_file(&lpByFits, lExtBy, READONLY, &lStatus) ;
    fits_movabs_hdu(lpByFits,2,&lDummy,&lStatus);
    fits_get_img_dim(lpByFits, &lNaxis, &lStatus) ;
    //if (lNaxis != 3) throw TFitsErr("NAXIS incorrect for By.") ;
    fits_get_img_size(lpByFits, 3, lNaxes, &lStatus) ;
    //if (lNaxes[0]!=_fNx || lNaxes[1]!=_fNy || lNaxes[2]!=_fNz) 
    //  throw TFitsErr("size of By array not ok with xml specifications.") ;
    fits_read_img(lpByFits, TDOUBLE, 1, _fNx*_fNy*_fNz,
		  NULL, lpBy, &lAnyNull, &lStatus) ;
    fits_close_file(lpByFits,&lStatus) ;
    //if (lStatus) throw TFitsErr(lStatus) ;
    //if (lAnyNull) throw TFitsErr("Undefined elements in Bx extension.") ;

    double *lpBz = new double[_fNx*_fNy*_fNz] ; 
    fitsfile *lpBzFits ;
    const char *lExtBz = (fitsFileName).c_str() ;
    fits_open_file(&lpBzFits, lExtBz, READONLY, &lStatus) ;
    fits_movabs_hdu(lpBzFits,3,&lDummy,&lStatus);
    fits_get_img_dim(lpBzFits, &lNaxis, &lStatus) ;
    //if (lNaxis != 3) throw TFitsErr("NAXIS incorrect for Bz.") ;
    fits_get_img_size(lpBzFits, 3, lNaxes, &lStatus) ;
    //if (lNaxes[0]!=_fNx || lNaxes[1]!=_fNy || lNaxes[2]!=_fNz) 
    //  throw TFitsErr("size of Bx array not ok with xml specifications.") ;
    fits_read_img(lpBzFits, TDOUBLE, 1, _fNx*_fNy*_fNz,
		  NULL, lpBz, &lAnyNull, &lStatus) ;
    fits_close_file(lpBzFits,&lStatus) ;
    //if (lStatus) throw TFitsErr(lStatus) ;
    //if (lAnyNull) throw TFitsErr("Undefined elements in Bz extension.") ;

    _fBmax = 0. ;

    long lInd, lIndFits ;
    for (int i=0; i<_fNx; i++) {
      for (int j=0; j<_fNy; j++) {
	for (int k=0; k<_fNz; k++) {
	  lIndFits = i+j*_fNx+k*_fNx*_fNy ;
	  lInd = i*_fNy*_fNz+j*_fNz+k ;
	  //(_fpB[lInd]).set(lpBx[lIndFits], lpBy[lIndFits], lpBz[lIndFits]);
	  //(_fpB[lInd]) *= gauss ; // unit in the fits file
	  //_fBmax = std::max(_fBmax, (_fpB[lInd]).mag() ) ;
	}
      }
    }
    delete[] lpBx;
    delete[] lpBy;
    delete[] lpBz ;
    delete[] lNaxes ;
    
    return;
  }

  Vector3d FITSMagneticField::getField(const Vector3d &position) const {
    return Vector3d(0,0,0);
  }
}

    // AMRMagneticField(const string& fitsFilename)
    //   { 
    // 	// 
    // 	long _fNx = 1;
    // 	long _fNy = 1;
    // 	long _fNz = 1;
    // 	//
	
	

	

	
	
    // 	_fBmax = 0. ;
    // 	long lInd, lIndFits ;
    // 	for (int i=0; i<_fNx; i++) {
    // 	  for (int j=0; j<_fNy; j++) {
    // 	    for (int k=0; k<_fNz; k++) {
    // 	      lIndFits = i+j*_fNx+k*_fNx*_fNy ;
    // 	      lInd = i*_fNy*_fNz+j*_fNz+k ;
    // 	      (_fpB[lInd]).set(lpBx[lIndFits], lpBy[lIndFits], lpBz[lIndFits]);
    // 	      (_fpB[lInd]) *= gauss ; // unit in the fits file
    // 	      //if( _fBmax < (_fpB[lInd]).mag() ){ 
    // 	      //std::cout<<"lInd"<<lInd<<"\t_fpB[lInd]).mag()="<<(_fpB[lInd]).mag()<<std::endl;
    // 	      //}
    // 	      _fBmax = max(_fBmax, (_fpB[lInd]).mag() ) ;
    // 	    }
    // 	  }
    // 	}
    // 	delete[] lpBx;
    // 	delete[] lpBy;
    // 	delete[] lpBz ;
    // 	delete[] lNaxes ;
	
    //   }
    
    // Vector3d getField(const Vector3d &position) const {

    //   /*double x = position.x/cfLength;
    //     double y = position.y/cfLength;
    //     double z = position.z/cfLength;

    //     std::vector<double> b;
    //     #ifdef _OPENMP
    //         #pragma omp critical 
    //         {    
    // 		        b = field->getField(x, y, z);
    //         }
    //     #else 
    //         b = field->getField(x, y, z);
    //     #endif

    //     for(int i=0; i<3; i++)
    //         b[i]*=cfMagneticField;
    //     //std::cout << x*cfLength/Mpc << "  " << y*cfLength/Mpc << "  " << z*cfLength/Mpc << "  ||   " << b[0] << " " << b[1] << "  " << b[0] << std::endl;

    // 	return Vector3d(b[0], b[1], b[2]) * tesla;   */
    // }
