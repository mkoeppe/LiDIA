//==============================================================================================
//
//	This file is part of LiDIA --- a library for computational number theory
//
//	Copyright (c) 1994--2001 the LiDIA Group.  All rights reserved.
//
//	See http://www.informatik.tu-darmstadt.de/TI/LiDIA/
//
//----------------------------------------------------------------------------------------------
//
//	$Id$
//
//	Author	: Volker Mueller (VM)
//	Changes	: See CVS log
//
//==============================================================================================


#ifndef LIDIA_POINT_OPERATIONS_CC_GUARD_
#define LIDIA_POINT_OPERATIONS_CC_GUARD_


#ifndef LIDIA_POINT_OPERATIONS_H_GUARD_
# include	"LiDIA/elliptic_curves/point_operations.h"
#endif
#ifndef LIDIA_BASE_POINT_H_GUARD_
# include	"LiDIA/elliptic_curves/base_point.h"
#endif
#ifndef LIDIA_ELLIPTIC_CURVE_FLAGS_H_GUARD_
# include	"LiDIA/elliptic_curve_flags.h"
#endif



#ifdef LIDIA_NAMESPACE
# ifndef IN_NAMESPACE_LIDIA
namespace LiDIA {
# endif
#endif



//
// constructor / destructor
//
template< class T >
point_operations< T >::point_operations ()
{
	debug_handler("point_operations< T >", "point_operations()");

	_add = NULL;
	_negate = NULL;
	_mult_by_2 = NULL;
}



template< class T >
point_operations< T >::~point_operations ()
{
	debug_handler("point_operations< T >", "~point_operations()");
}



//
// assigment
//
template< class T >
void
point_operations< T >::set (elliptic_curve_flags::curve_parametrization cp,
			     elliptic_curve_flags::curve_model cm,
			     const bigint& characteristic)
{
	// Check for the model: affine or projective
	//
	if (cm == elliptic_curve_flags::AFFINE) {
		// Verify the parametrization: SHORT_W, LONG_W
		//
		switch (cp) {
		case elliptic_curve_flags::SHORT_W:
			_add = add_swnf_affine;
			_negate = negate_swnf_affine;
			_mult_by_2 = mult_by_2_swnf_affine;
			break;

		case  elliptic_curve_flags::LONG_W:
			_add = add_lwnf_affine;
			_negate = negate_lwnf_affine;
			_mult_by_2 = mult_by_2_lwnf_affine;
			break;
		default: lidia_error_handler("point_operations< T >::set",
					     "Wrong type.");
		}
	}
	else {
		switch (cp) {
		case elliptic_curve_flags::SHORT_W:
			_add = add_swnf_projective;
			_negate = negate_swnf_projective;
			_mult_by_2 = mult_by_2_swnf_projective;
			break;

		case  elliptic_curve_flags::LONG_W:
			_add = add_lwnf_projective;
			_negate = negate_lwnf_projective;
			_mult_by_2 = mult_by_2_lwnf_projective;
			break;
		default: lidia_error_handler("point_operations< T >::set",
					     "Wrong type.");
		}
	}
}



template< class T >
point_operations< T > &
point_operations< T >::operator = (const point_operations< T > & p)
{
	if (this != &p) {
		_add = p._add;
		_negate = p._negate;
		_mult_by_2 = p._mult_by_2;
	}
	return *this;
}



//
// Implementation of basic arithmetic on points
// used in base_point.cc.
//

//***************************************************************
//  Addition of different points (this is NOT checked !!)
//**************************************************************


template < class T >
void add_swnf_affine (base_point < T > &erg,
		      const base_point < T > &wp1,
		      const base_point < T > &wp2)
{
  debug_handler ("point_operations.c", "add_swnf_affine(R,P,Q)");
  T lambda, h1 (wp1.y), h2 (wp1.x);
  subtract (h2, wp2.x, h2); if (h2.is_zero ())
    {
      erg.assign_zero (wp1.ec); return;
    }
  
  h2.invert (); 
  subtract (h1, wp2.y, h1); 
  multiply (lambda, h1, h2);	// gradient

  square (h1, lambda); 
  subtract (h1, h1, wp1.x); 
  subtract (h1, h1, wp2.x);	// x-coordinate

  subtract (h2, wp1.x, h1); 
  multiply (h2, h2, lambda); 
  subtract (h2, h2, wp1.y);	// y-coordinate

#ifdef DEBUG
     erg.assign (h1, h2, wp1.ec);
#else
     erg.assign_no_test (h1, h2, wp1.ec);
#endif
     }

//------------------------------------------------------------------------

template < class T >
void add_lwnf_affine (base_point < T > &erg, const base_point < T > &wp1,
		      const base_point < T > &wp2)
{
  debug_handler ("point_operations.c", "add_lwnf_affine(R,P,Q)");

  T lambda, x3, y3;
  
  subtract(x3, wp2.x, wp1.x);

  if (x3.is_zero())
    {
      erg.assign_zero(wp1.ec);
      return;
    }
  
  x3.invert();
  subtract(lambda, wp2.y, wp1.y);
  multiply(lambda, lambda, x3);  // gradient
  
  add(x3, lambda, wp1.get_curve().get_a1());
  multiply(x3, lambda, x3);
  subtract(x3, x3, wp1.get_curve().get_a2());
  subtract(x3, x3, wp1.x);
  subtract(x3, x3, wp2.x);       // x-coordinate

  subtract(y3, wp1.x, x3);
  multiply(y3, y3, lambda);
  subtract(y3, y3, wp1.y);
  multiply(lambda, x3, wp1.get_curve().get_a1());
  subtract(y3, y3, lambda);
  subtract(y3, y3, wp1.get_curve().get_a3());  // y-coordinate

#ifdef DEBUG
  erg.assign(x3, y3, wp1.ec);
#else
  erg.assign_no_test(x3, y3, wp1.ec);
#endif

}

//----------------------------------------------------------------------

template < class T >
void add_swnf_projective (base_point < T > &erg,
			  const base_point < T > &wp1,
			  const base_point < T > &wp2)
{
  debug_handler ("point_operations.c", "add_swnf_projective(R,P,Q)");

    if(wp1.is_zero())
    {
      erg.assign(wp2);   
      return;
    }
  
  if (wp2.is_zero())
    {
      erg.assign(wp1);     
      return;
    }

  T h1 (wp1.x), h2(wp1.y), h3(wp1.z), h4(wp2.x), h5(wp2.y), h6, h7;

  if (!wp2.z.is_one()) 
    {
      square(h7, wp2.z);              
      multiply(h1, wp1.x, h7);        // U_0 = X_0 * Z_1^2
      multiply(h7, wp2.z, h7);        
      multiply(h2, wp1.y, h7);        // S_0 = y_0 * Z_1^3
    }

  if (!wp1.z.is_one())
    {
      square(h7, wp1.z);      
      multiply(h4, wp2.x, h7);   // U_1 = X_1 * Z_0^2
      multiply(h7, wp1.z, h7); 
      multiply(h5, wp2.y, h7);   // S_1 = Y_1 * Z_0^3
    }

  subtract(h4, h1, h4);         // W = U_0 - U_1    
  subtract(h5, h2, h5);         // R = S_0 - S_1

  if (h4.is_zero()) 
    {                   
      if (h5.is_zero())
        mult_by_2_swnf_projective(erg, wp1); 
      else 
        erg.assign_zero(wp1.ec); 
      return;
    }

  h1.multiply_by_2();
  subtract(h1, h1, h4);       // T = U_0 + U_1
  h2.multiply_by_2();             
  subtract(h2, h2, h5);       // M = S_0 + S_1

  if (!wp2.z.is_one())
    multiply(h3, h3, wp2.z);     
  multiply(h3, h3, h4);      // Z_2 = Z_0 * Z_1 * W

  square(h7, h4);    
  multiply(h4, h4, h7);   
  multiply(h7, h1, h7);     
  square(h1, h5);           
  subtract(h1, h1, h7);     // X_2 = R^2 - T * W^2

  add(h6, h1, h1);          
  subtract(h7, h7, h6);     // V = T * W^2 - 2 * X_2

  multiply(h5, h5, h7);     
  multiply(h4, h2, h4);     
  subtract(h2, h5, h4);    // 2*Y_2 = V*R - M * W^3 
  h2.divide_by_2();        

#ifdef DEBUG
  erg.assign(h1, h2, h3, wp1.ec);  
#else
  erg.assign_no_test(h1, h2, h3, wp1.ec);
#endif
}


//------------------------------------------------------------------

template < class T >
void add_lwnf_projective (base_point < T > &erg,
			  const base_point < T > &wp1,
			  const base_point < T > &wp2)
{
  debug_handler ("point_operations.c", "add_lwnf_projective(R,P,Q)");

  if(wp1.is_zero())
    {
      erg.assign(wp2);   
      return;
    }
  
  if (wp2.is_zero())
    {
      erg.assign(wp1);                           
      return;
    }

  T h1, h2, h3, h4, h5, h6, h7, h8, h9, z31;
  
  if (!wp2.z.is_one()) 
    {
      square(h7, wp2.z);              
      multiply(h1, wp1.x, h7);        // U_0 = X_0 * Z_1^2
      multiply(h7, wp2.z, h7);        
      multiply(h2, wp1.y, h7);        // S_0 = y_0 * Z_1^3
    }
  else
    {
      h1.assign(wp1.x);
      h2.assign(wp1.y);
    }

  if (!wp1.z.is_one())
    {
      square(h7, wp1.z);      
      multiply(h4, wp2.x, h7);   // U_1 = X_1 * Z_0^2
      multiply(h7, wp1.z, h7); 
      multiply(h5, wp2.y, h7);   // S_1 = Y_1 * Z_0^3
    }
  else
    {
      h4.assign(wp2.x);
      h5.assign(wp2.y);
    }

  subtract(h3, h1, h4);         // W = U_0 - U_1    
  subtract(h6, h2, h5);         // R = S_0 - S_1

  if (h3.is_zero()) 
    {                   
      if (h6.is_zero())
        mult_by_2_lwnf_projective(erg, wp1); 
      else 
        erg.assign_zero(wp1.ec); 
      return;
    }

  add(h1, h1, h4);        // T = U_0 + U_1
  add(h2, h2, h5);        // M = S_0 + S_1

  multiply(z31, h3, wp2.z);     // Z_2 / Z_0
  multiply(h4, z31, wp1.z);      // Z_2 = Z_0 * Z_1 * W

  square(h8, h3);    
  multiply(h7, h1, h8);     // T*W^2
  multiply(h8, h8, h3);     // W^3
  square(h5, h6);           
  subtract(h5, h5, h7);     // R^2 - T * W^2
  multiply(h1, wp1.get_curve().get_a1(), h6);
  multiply(h3, wp1.get_curve().get_a2(), h4);
  subtract(h1, h1, h3);
  multiply(h1, h4, h1);
  add(h5, h5, h1);

           // X_2 = R^2 - T * W^2 + Z_2 * (a1 * R - a2 * Z_2)

  multiply(h1, wp1.get_curve().get_a1(), h5);
  square(h3, h4);
  multiply(h2, wp1.get_curve().get_a3(), h3);
  add(h1, h1, h2);
  multiply(h1, h1, h4);

  square(h2, z31);
  multiply(h3, h2, wp1.x);
  subtract(h3, h3, h5);
  multiply(h3, h3, h6);
  subtract(h3, h3, h1);
  
  multiply(h2, h2, z31);
  multiply(h2, h2, wp1.y);
  subtract(h3, h3, h2);
  
  // Y_2 = R * (X_0 * (Z_2/Z_0)^2 - X_2) - Y_0*(Z_2/Z_0)^3 - Z_2 * 
  //       (a1 *X_2 + a3 * Z_2^2)

#ifdef DEBUG
  erg.assign(h5, h3, h4, wp1.ec);  
#else
  erg.assign_no_test(h5, h3, h4, wp1.ec);
#endif
}



//***************************************************************
// Doubling routines 
//**************************************************************/


template < class T >
void mult_by_2_swnf_affine (base_point < T > &erg,
			    const base_point < T > &wp)
{
  debug_handler ("point_operations.c", "mult_by_2_swnf_affine");

   T lambda(wp.y), hx(wp.x), hy(wp.y);

  hy.multiply_by_2();
  
  if (hy.is_zero())
    {
      erg.assign_zero(wp.get_curve());
      return;
    }

  hy.invert();
  square(hx, hx);
  multiply(hx, hx, bigint(3));
  add(hx, hx, wp.get_curve().get_a4());
  multiply(lambda, hx, hy);            // gradient
  
  square(hx, lambda);
  hy.assign(wp.x);
  hy.multiply_by_2();
  subtract(hy, hx, hy);                // x-coordinate

  subtract(hx, wp.x, hy);
  multiply(hx, lambda, hx);
  subtract(hx, hx, wp.y);              // y-coordinate

#ifdef DEBUG  
  erg.assign(hy, hx, wp.get_curve());
#else
  erg.assign_no_test(hy, hx, wp.get_curve());
#endif
}


//-------------------------------------------------------------------

template < class T >
void mult_by_2_lwnf_affine (base_point < T > &erg,
			    const base_point < T > &wp)
{
  debug_handler ("point_operations.c", "mult_by_2_lwnf_affine");

   T lambda(wp.get_curve().get_a3()), x3(wp.y), y3(wp.x);
  
  x3.multiply_by_2();
  add(lambda, lambda, x3);
  multiply(x3, wp.get_curve().get_a1(), y3);
  add(lambda, lambda, x3);

  if (lambda.is_zero())
    {
      erg.assign_zero(wp.get_curve());
      return;
    }
  
  lambda.invert();
  x3.assign(wp.x);
  x3.multiply_by_2();
  add(x3, x3, wp.x);
  add(y3, wp.get_curve().get_a2(), wp.get_curve().get_a2());
  add(x3, x3, y3);
  multiply(x3, x3, wp.x);
  add(x3, x3, wp.get_curve().get_a4());
  multiply(y3, wp.get_curve().get_a1(), wp.y);
  subtract(x3, x3, y3);
  multiply(lambda, x3, lambda);               // gradient
  
  add(x3, lambda, wp.get_curve().get_a1());
  multiply(x3, x3, lambda);
  subtract(x3, x3, wp.get_curve().get_a2());
  add(y3, wp.x, wp.x);
  subtract(x3, x3, y3);           // x-coordinate

  subtract(y3, wp.x, x3);
  multiply(y3, y3, lambda);
  subtract(y3, y3, wp.y);
  subtract(y3, y3, wp.get_curve().get_a3());
  multiply(lambda, wp.get_curve().get_a1(), x3);
  subtract(y3, y3, lambda);       // y-coordinate

#ifdef DEBUG
  erg.assign(x3, y3, wp.get_curve());
#else
  erg.assign_no_test(x3, y3, wp.get_curve());
#endif
}

//------------------------------------------------------------------

template < class T >
void mult_by_2_swnf_projective (base_point < T > &erg,
				const base_point < T > &wp)
{
  debug_handler ("point_operations.c", "mult_by_2_swnf_projective");

  
  if (wp.is_zero())
    {           
      erg.assign_zero(wp.get_curve());
      return;
    }

  T  h1(wp.x), h2(wp.y), h3(wp.z), h4, h5, h;
  
  square(h5, h3); 
  square(h5, h5);                  
  multiply(h5, wp.get_curve().get_a4(), h5); 
  square(h4, h1);                
  add(h, h4, h4);
  add(h4, h, h4);                   
  add(h4, h4, h5);             // M = 3 x^2 + a4 * z^4
  
  multiply(h3, h2, h3);                            
  h3.multiply_by_2();          // Z_2 = h3 = 2 * y * z

  square(h2, h2);    
  multiply(h5, h1, h2);                        
  h5.multiply_by_2();                        
  h5.multiply_by_2();         // S = 4 * x * y^2
                
  square(h1, h4);                            
  add(h, h5, h5);  
  subtract(h1, h1, h);        // X_2 = M^2 - 2 * S

  square(h2, h2);       
  h2.multiply_by_2();   
  h2.multiply_by_2();   
  h2.multiply_by_2();         // T = 8 * y^4

  subtract(h5, h5, h1);
  multiply(h5, h4, h5);  
  subtract(h2, h5, h2);       // Y_2 = M * (S-X_2) - T

#ifdef DEBUG
    erg.assign(h1, h2, h3, wp.get_curve());                    
#else
    erg.assign_no_test(h1, h2, h3, wp.get_curve());
#endif
} 


//--------------------------------------------------------------------

template < class T >
void mult_by_2_lwnf_projective (base_point < T > &erg,
				const base_point < T > &wp)
{
  debug_handler ("point_operations.c", "mult_by_2_lwnf_projective");

    if (wp.is_zero())
    {
      erg.assign_zero(wp.get_curve());
      return;
    }

  T m, t, x2, y2, z2, h1, h2, h3;

  square(h1, wp.z);
  multiply(t, wp.get_curve().get_a1(), wp.x);
  multiply(h2, wp.get_curve().get_a3(), h1);
  add(t, t, h2);
  multiply(t, t, wp.z);
  add(t, t, wp.y);
  add(t, t, wp.y);    // T = 2*Y + (a1 * X + a3 * Z^2) * Z

  if (t.is_zero())
    {
      erg.assign_zero(wp.get_curve());
      return;
    }

  multiply(z2, t, wp.z);   // Z_2

  multiply(h2, wp.get_curve().get_a2(), wp.x);
  add(h2, h2, h2);
  multiply(m, wp.get_curve().get_a4(), h1);
  add(m, m, h2);
  multiply(m, m, wp.z);
  multiply(h2, wp.get_curve().get_a1(), wp.y);
  subtract(m, m, h2);
  multiply(m, m, wp.z);
  square(h2, wp.x);
  add(m, m, h2);
  add(m, m, h2);
  add(m, m, h2); // M = 3*X^2 + Z*( Z*(2*a2*X + a4*Z^2) - a1*Y)

  square(h1, z2);
  multiply(h2, wp.get_curve().get_a1(), z2);
  add(x2, m, h2);
  multiply(x2, m, x2);
  multiply(h2, wp.get_curve().get_a2(), h1);
  subtract(x2, x2, h2);
  square(h2, t);
  multiply(h3, h2, wp.x);
  subtract(x2, x2, h3);
  subtract(x2, x2, h3);    // X_2 =  M * (M + a1*Z_2) - a2*Z_2^2 - 2*X*T^2

  subtract(y2, h3, x2);
  multiply(y2, y2, m);
  multiply(m, h2, t);
  multiply(m, m, wp.y);
  subtract(y2, y2, m);
  multiply(m, wp.get_curve().get_a1(), x2);
  multiply(h1, h1, wp.get_curve().get_a3());
  add(m, m, h1);
  multiply(m, m, z2);
  subtract(y2, y2, m);  // Y_2 = (T^2*X - X_2) * M - Y*T^3 - Z_2 * (a1*X_2
                        //       + a3*Z_2^2)
  
#ifdef DEBUG
    erg.assign(x2, y2, z2, wp.get_curve());                    
#else
    erg.assign_no_test(x2, y2, z2, wp.get_curve());
#endif
  
}



/***************************************************************
  Negation
  **************************************************************/

template < class T >
void negate_swnf_affine (base_point < T > &erg,
			 const base_point < T > &P)
{
  debug_handler ("point_operations.c", "negate_swnf_affine");

 erg = P;

  if (!P.is_zero())
    negate(erg.y, erg.y);
}


//---------------------------------------------------------------------

template < class T >
void negate_lwnf_affine (base_point < T > &erg,
			 const base_point < T > &P)
{
  debug_handler ("point_operations.c", "negate_lwnf_affine");

 T h;
 erg = P;
 
 if (P.is_zero())
   return;
 
 multiply(h, P.get_curve().get_a1(), P.x);
 add(erg.y, h, P.y);
 add(erg.y, erg.y, P.get_curve().get_a3());
 negate(erg.y, erg.y);
}


//--------------------------------------------------------------------

template < class T >
void negate_swnf_projective (base_point < T > &erg,
			     const base_point < T > &P)
{
  debug_handler ("point_operations.c", "negate_swnf_projective");

  erg = P;
  if (!P.is_zero())
    negate(erg.y, erg.y);
}


//--------------------------------------------------------------------

template < class T >
void negate_lwnf_projective (base_point < T > &erg,
			     const base_point < T > &P)
{
  debug_handler ("point_operations.c", "negate_lwnf_projective");

  T h, h2;

  erg = P;
  if (P.is_zero())
    return;

  multiply(h, P.get_curve().get_a1(), P.x);
  square(h2, P.z);
  multiply(h2, h2, P.get_curve().get_a3());
  add(h2, h2, h);
  multiply(h2, h2, P.z);
  add(h2, P.y, h2);
  negate(erg.y, h2);     // -(X:Y:Z) = (X: -Y - a1*X*Z - a3*Z^3 : Z) 
}






#ifdef LIDIA_NAMESPACE
# ifndef IN_NAMESPACE_LIDIA
}	// end of namespace LiDIA
# endif
#endif



#endif	// LIDIA_POINT_OPERATIONS_CC_GUARD_
