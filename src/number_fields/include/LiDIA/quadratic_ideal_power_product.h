// -*- C++ -*-
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
//	Author	: Markus Maurer (MM)
//	Changes	: See CVS log
//
//==============================================================================================


#ifndef LIDIA_QUADRATIC_IDEAL_POWER_PRODUCT_H_GUARD_
#define LIDIA_QUADRATIC_IDEAL_POWER_PRODUCT_H_GUARD_


#ifndef LIDIA_BASE_POWER_PRODUCT_H_GUARD_
#include	"LiDIA/base_power_product.h"
#endif
#ifndef LIDIA_QUADRATIC_IDEAL_H_GUARD_
#include	"LiDIA/quadratic_ideal.h"
#endif
#ifndef LIDIA_BIGINT_H_GUARD_
#include	"LiDIA/bigint.h"
#endif
#ifndef LIDIA_BASE_VECTOR_H_GUARD_
#include	"LiDIA/base_vector.h"
#endif



#ifdef LIDIA_NAMESPACE
namespace LiDIA {
# define IN_NAMESPACE_LIDIA
#endif



class quadratic_number_standard;



class quadratic_ideal_power_product : public base_power_product< quadratic_ideal, bigint >
{
public:
	//
	// constructor / destructor
	//

	quadratic_ideal_power_product ();
	~quadratic_ideal_power_product ();

	//
	// assignment
	//
	quadratic_ideal_power_product & operator = (
		const quadratic_ideal_power_product &);

	void assign (const base_vector< quadratic_number_standard > & b,
		     const base_vector< bigint > & e);
	//
	// factor refinement
	//
	int factor_refinement();

	bool is_gcd_free() const;
};



#ifdef LIDIA_NAMESPACE
}	// end of namespace LiDIA
# undef IN_NAMESPACE_LIDIA
#endif



#endif	// LIDIA_QUADRATIC_IDEAL_POWER_PRODUCT_H_GUARD_
