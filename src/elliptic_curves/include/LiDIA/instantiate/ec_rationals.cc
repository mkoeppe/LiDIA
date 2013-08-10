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
//	Author	: Volker Mueller (VM), Markus Maurer  (MM)
//	Changes	: See CVS log
//
//==============================================================================================


//
//  usage:
//
//	include all headers that are necessary to describe the type TYPE
//
//	define the type TYPE
//
//	include this file
//


#define HAS_NO_GET_FIELD

#include	"LiDIA/elliptic_curves/base_elliptic_curve_rep.h"
#include	"LiDIA/elliptic_curves/elliptic_curve_rep.h"
#include	"LiDIA/elliptic_curves/base_elliptic_curve.h"
#include	"LiDIA/elliptic_curves/base_point.h"
#include	"LiDIA/elliptic_curves/point_operations.h"

#include	"LiDIA/elliptic_curve.h"
#include	"LiDIA/point.h"

#include	"LiDIA/elliptic_curves/base_elliptic_curve_rep.cc"
#include	"LiDIA/elliptic_curves/elliptic_curve_rep.cc"
#include	"LiDIA/elliptic_curves/base_elliptic_curve.cc"
#include	"LiDIA/elliptic_curves/base_point.cc"
#include	"LiDIA/elliptic_curves/point_operations.cc"

#include	"LiDIA/elliptic_curve.cc"
#include	"LiDIA/point.cc"



#ifdef LIDIA_NAMESPACE
namespace LiDIA {
#endif



//****************** instantiation of classes *****************************

//
// base_elliptic_curve_rep
//

template class base_elliptic_curve_rep< TYPE >;


//
// base_elliptic_curve
//

template class base_elliptic_curve< TYPE >;


#ifndef BASE

//
// elliptic_curve_rep
//

template class elliptic_curve_rep< TYPE >;


//
// elliptic_curve
//

template class elliptic_curve< TYPE >;


#endif



//
// base_point
//

template class base_point< TYPE >;


template void negate(base_point< TYPE > & r,
		     const base_point< TYPE > & x);

template void add(base_point< TYPE > & R,
		  const base_point< TYPE > & P1,
		  const base_point< TYPE > & Q);

template void subtract(base_point< TYPE > & r,
		       const base_point< TYPE > & x,
		       const base_point< TYPE > & y);

template void multiply_by_2(base_point< TYPE > & R,
			    const base_point< TYPE > & P);

template void multiply (base_point< TYPE > & ,
			const bigint &,
			const base_point< TYPE > &);

template void swap(base_point< TYPE > & P,
		   base_point< TYPE > & Q);



#ifndef BASE

//
// point
//

template class point< TYPE >;

#endif


//
// point_operations
//

template class point_operations< TYPE >;


//
// add
//

template void add_swnf_affine(base_point< TYPE > &,
			      const base_point< TYPE > &,
			      const base_point< TYPE > &);

template void add_lwnf_affine(base_point< TYPE > &,
			      const base_point< TYPE > &,
			      const base_point< TYPE > &);

template void add_swnf_projective(base_point< TYPE > &,
				  const base_point< TYPE > &,
				  const base_point< TYPE > &);

template void add_lwnf_projective(base_point< TYPE > &,
				  const base_point< TYPE > &,
				  const base_point< TYPE > &);


//
// negate
//

template void negate_swnf_affine(base_point< TYPE > &,
				 const base_point< TYPE > &);


template void negate_lwnf_affine(base_point< TYPE > &,
				 const base_point< TYPE > &);

template void negate_swnf_projective(base_point< TYPE > &,
				     const base_point< TYPE > &);


template void negate_lwnf_projective(base_point< TYPE > &,
				     const base_point< TYPE > &);



//
// mult_by_2
//

template void mult_by_2_swnf_affine(base_point< TYPE > &,
				    const base_point< TYPE > &);

template void mult_by_2_lwnf_affine(base_point< TYPE > &,
				    const base_point< TYPE > &);

template void mult_by_2_swnf_projective(base_point< TYPE > &,
					const base_point< TYPE > &);

template void mult_by_2_lwnf_projective(base_point< TYPE > &,
					const base_point< TYPE > &);





#ifdef LIDIA_NAMESPACE
}	// end of namespace LiDIA
#endif
