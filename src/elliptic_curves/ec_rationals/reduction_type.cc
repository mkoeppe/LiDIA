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
//	Author	: Nigel Smart
//		  Adaption of John Cremona's code
//                (some of which itself is based on code of
//                Oisin McGuiness and Richard Pinch)
//	Changes	: See CVS log
//
//==============================================================================================


#ifdef HAVE_CONFIG_H
# include	"config.h"
#endif
#include	"LiDIA/kodaira_code.h"
#include	"LiDIA/reduction_type.h"



#ifdef LIDIA_NAMESPACE
namespace LiDIA {
#endif



reduction_type
::reduction_type(int opd, int opN, int opj, int kc, int cp)

	: ord_p_discr(opd), ord_p_N(opN), ord_p_j_denom(opj), Kcode(kc),
	  c_p(cp)
{
	debug_handler("reduction_type",
		      "reduction_type(int, int, int, int, int)");
}



reduction_type & reduction_type
::operator = (const reduction_type& r)
{
	debug_handler("reduction_type",
		      "operator = (const reduction_type&)");

	if (this != &r) {
		ord_p_discr = r.ord_p_discr; ord_p_N = r.ord_p_N;
		ord_p_j_denom = r.ord_p_j_denom; Kcode = r.Kcode; c_p = r.c_p;
	}

	return *this;
}



//
// Access Functions
//

int reduction_type
::get_ord_p_discriminant() const
{
	debug_handler("reduction_type", "get_ord_p_discriminant()");
	return ord_p_discr;
}



int reduction_type
::get_ord_p_N() const
{
	debug_handler("reduction_type", "get_ord_p_N()");
	return ord_p_N;
}



int reduction_type
::get_ord_p_j_denom() const
{
	debug_handler("reduction_type", "ord_p_j_denom");
	return ord_p_j_denom;
}



int reduction_type
::get_c_p() const
{
	debug_handler("reduction_type", "get_c_p()");
	return c_p;
}



Kodaira_code reduction_type
::get_code()	const
{
	debug_handler("reduction_type", "get_code()");
	return Kcode;
}



std::ostream& operator << (std::ostream& os, const reduction_type& R)
{
	debug_handler("reduction_type",
		      "operator << (std::ostream&, const reduction_type&)");

	os << R.ord_p_discr << "\t";
	os << R.ord_p_N << "\t";
	os << R.ord_p_j_denom << "\t\t";
	os << R.Kcode << "\t";
	os << R.c_p;
	return os;
}



#ifdef LIDIA_NAMESPACE
}	// end of namespace LiDIA
#endif
