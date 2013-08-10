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
//	Author	: Markus Maurer (MM), Patrick Theobald (PT)
//	Changes	: See CVS log
//
//==============================================================================================


#ifndef LIDIA_BASE_POWER_PRODUCT_CC_GUARD_
#define LIDIA_BASE_POWER_PRODUCT_CC_GUARD_


#ifndef LIDIA_BASE_POWER_PRODUCT_H_GUARD_
# include	"LiDIA/base_power_product.h"
#endif



#ifdef LIDIA_NAMESPACE
# ifndef IN_NAMESPACE_LIDIA
namespace LiDIA {
# endif
#endif



//************************************************************
// Sometimes, we have to assign '1'. In those cases, there must
// be an reference element in the base_power_product.
// Otherwise, if the
// base_power_product is empty, the lidia_error_handler is called.
// This is necessary, because it is possible that the element '1'
// does not only depend on the type T, but also on internal
// variables of type T elements (i.e., the modulus like in
// multi_bigmod, see LiDIA manual).
//
// Examples:
// 
// Correct:
//
// T x;
// x = base(0);
// x.assign_one(); // now x contains the modulus of base(0)
//
// Incorrect:
//
// T x;
// x.assign:one(); // What is the modulus of x ??
//*************************************************************

//
// input / output
//

template< class T, class exp_type >
void
base_power_product< T, exp_type >::read (std::istream &is)
{
	//
	// I/O - format :
	//
	// [ [base_1, exp_1], [base_2, exp_2], ..., [base_i, exp_i] ]
	//
	// [] yields an empty base_power_product
	// (must be possible, because this is the output
	//  in case of an empty base_power_product, s. write)
	//

	exp_type exp;
	T base;
	char c;

	this->reset();

	// read outermost opening '['
	is >> c;
	if (c != '[')
		lidia_error_handler("base_power_product< T, exp_type >",
				    "read(std::istream&)::[ expected");

	is >> c;
	if (c == ']')
		is.putback(c);
	else {
		is.putback(c);

		do {
			is >> c;
			if (c != '[')
				lidia_error_handler("base_power_product< T, exp_type >",
						    "read(std::istream&)::[ expected");

			is >> base; // base of the component
			is >> c;
			switch (c) {
			case (',') : is >> exp; // exp of the component
				this->append(base, exp);
				is >> c;
				if (c != ']')
					lidia_error_handler("base_power_product< T, exp_type >",
							    "read(std::istream&)::']' expected");
				break;
			default :
				lidia_error_handler("base_power_product< T, exp_type >",
						    "read(std::istream&)::',' or ']' expected");
			}

			is >> c;
		} while (c == ',');

	}
	// end if (c == ']')

	// read outermost closing ']'
	if (c != ']')
		lidia_error_handler("base_power_product< T, exp_type >",
				    "read(std::istream&)::']' expected");
}



template< class T, class exp_type >
void
base_power_product< T, exp_type >::write (std::ostream &os) const
{
	//
	// I/O - format :
	//
	// [ [base_1, exp_1], [base_2, exp_2], ..., [base_i, exp_i] ]
	//
	// [] in case of an empty base_power_product
	//

	lidia_size_t i;
	lidia_size_t nc;

	nc = get_no_of_components();

	os << "[";

	for (i = 0; i < nc; i++) {
		os << " [";
		os << this->get_base(i);
		os << ", ";
		os << this->get_exponent(i);
		os << "] ";
	}

	os << "]";
}



//
// high-level functions
//

template< class T, class exp_type >
void
base_power_product< T, exp_type >::append (const T &a, exp_type exp)
{
	// append T a and exp_type exp
	// to the base_power_product

	if (exp != 0) {
		lidia_size_t nc = get_no_of_components();

		component[nc].left() = a;
		component[nc].right() = exp;

		clear_attributes();
	}
}



#ifdef LIDIA_NAMESPACE
# ifndef IN_NAMESPACE_LIDIA
}	// end of namespace LiDIA
# endif
#endif



#endif	// LIDIA_BASE_POWER_PRODUCT_CC_GUARD_
