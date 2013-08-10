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
//	Author	:
//	Changes	: See CVS log
//
//==============================================================================================


#ifndef LIDIA_COEFF_SPARSE_POWER_SERIES_CC_GUARD_
#define LIDIA_COEFF_SPARSE_POWER_SERIES_CC_GUARD_


#ifndef LIDIA_COEFF_SPARSE_POWER_SERIES_H_GUARD_
# include	"LiDIA/finite_fields/coeff_sparse_power_series.h"
#endif



#ifdef LIDIA_NAMESPACE
# ifndef IN_NAMESPACE_LIDIA
namespace LiDIA {
# endif
#endif



//
//  *****  I/O  *****
//

template< class T >
void
spc< T >::read (std::istream & in)
{
	debug_handler ("sparse_powser::spc" , "read");

	char ch;

	in >> ch;
	if (ch != '[')
		lidia_error_handler ("sparse_powser::spc" , "read::character ']' expected");

	in >> this->coeff;
	in >> this->exp;

	in >> ch;
	if (ch != ']')
		lidia_error_handler ("sparse_powser::spc" , "read::character ']' expected");
}



template< class T >
void
spc< T >::write (std::ostream & out) const
{
	debug_handler ("sparse_powser::spc" , "print");
	out << "[ " << this->coeff << " " << this->exp << " ]";
}



template< class T >
int
spc< T >::coeff_cmp_zero (const spc< T > & b) const
{
	debug_handler ("sparse_powser::spc" , "coeff_cmp_zero()");

	T zero_T(0);
	int rc;

	// This relation says :  1.) 0 is greater than any other element
	//                       2.) any non-zero elements will be identified
	//
	// Sorting a vector of spcs according to this function will
	// move the zeros into the rear part of the vector


	if (this->coeff == zero_T && b.coeff == zero_T)
		rc = 0;
	else if (this->coeff == zero_T)
		rc = 1;
	else if (b.coeff == zero_T)
		rc = -1;
	else
		rc = 0;

	return rc;
}



#ifdef LIDIA_NAMESPACE
# ifndef IN_NAMESPACE_LIDIA
}	// end of namespace LiDIA
# endif
#endif



#endif	// LIDIA_COEFF_SPARSE_POWER_SERIES_CC_GUARD_
