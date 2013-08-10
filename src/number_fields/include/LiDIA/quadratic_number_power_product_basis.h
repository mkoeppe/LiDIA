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


#ifndef LIDIA_QUADRATIC_NUMBER_POWER_PRODUCT_BASIS_H_GUARD_
#define LIDIA_QUADRATIC_NUMBER_POWER_PRODUCT_BASIS_H_GUARD_


#ifndef LIDIA_BIGINT_H_GUARD_
# include	"LiDIA/bigint.h"
#endif
#ifndef LIDIA_XBIGFLOAT_H_GUARD_
# include	"LiDIA/xbigfloat.h"
#endif
#ifndef LIDIA_BASE_VECTOR_H_GUARD_
# include	"LiDIA/base_vector.h"
#endif
#ifndef LIDIA_QUADRATIC_NUMBER_STANDARD_H_GUARD_
# include	"LiDIA/quadratic_number_standard.h"
#endif
#ifndef LIDIA_MATRIX_H_GUARD_
# include	"LiDIA/matrix.h"
#endif



#ifdef LIDIA_NAMESPACE
namespace LiDIA {
# define IN_NAMESPACE_LIDIA
#endif



class quadratic_order;



class quadratic_number_power_product_basis
{
private:

	static xbigfloat dummy;

	base_vector< quadratic_number_standard > basis;
	// elements of the basis

	mutable base_vector< xbigfloat > ln_basis;
	// approx. to their logarithms

	mutable base_vector< long > prec_ln_basis;
	// their precisions

	mutable base_vector< char > init;
	// ln_basis[i] valid ?

public:

	//
	//  constructor
	//
	quadratic_number_power_product_basis();

	//
	//  destructor
	//
	~quadratic_number_power_product_basis();


	//
	//  assignment
	//
	quadratic_number_power_product_basis & operator = (const quadratic_number_power_product_basis & b);

	void assign (const quadratic_number_power_product_basis & b);

	void set_basis (const base_vector< quadratic_number_standard > & b);
	void set_basis (const quadratic_number_standard & b);

	void concat (const quadratic_number_power_product_basis & b,
		     const quadratic_number_standard & q);

	void concat (const quadratic_number_power_product_basis & b,
		     const quadratic_number_power_product_basis & c);

	//
	//  access
	//

	const base_vector< quadratic_number_standard > & get_basis () const;
	lidia_size_t get_size() const;

	const quadratic_number_standard & operator[](lidia_size_t i) const;
	quadratic_number_standard & operator[](lidia_size_t i);

	const quadratic_order & get_order () const;


	//
	//  comparison. Only compares the addresses !
	//
	friend bool operator == (const quadratic_number_power_product_basis &,
				 const quadratic_number_power_product_basis &);

	//
	//  Assigns to l an absolute k-approximation of the Ln of the
	//  i-th base element. If the i is out of range the error handler
	//  will be called. If ln_prec_basis[i] is sufficient accurate,
	//  the computation is only cutting this element to precision k.
	//  Otherwise the Ln is recomputed and stored in ln_prec_basis[i].
	//

	const xbigfloat & get_absolute_Ln_approximation (long, lidia_size_t) const;

	//
	//
	//

	void absolute_Ln_approximations (const matrix< bigint > & M, long k);

	bool check_Ln_correctness (lidia_size_t trials = 5) const;


	//
	//  conjugate
	//

	void conjugate();


	//
	//   output operator
	//
	friend std::istream & operator >> (std::istream & in, quadratic_number_power_product_basis & b);

	//
	//   output operator
	//
	friend std::ostream & operator << (std::ostream & out, const quadratic_number_power_product_basis & b);
};



#ifdef LIDIA_NAMESPACE
}	// end of namespace LiDIA
# undef IN_NAMESPACE_LIDIA
#endif



#endif	// LIDIA_QUADRATIC_NUMBER_POWER_PRODUCT_BASIS_H_GUARD_
