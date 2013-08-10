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
//	Author	: Nigel Smart
//		  Adaption of John Cremona's code
//                (some of which itself is based on code of
//                Oisin McGuiness and Richard Pinch)
//	Changes	: See CVS log
//
//==============================================================================================


#ifndef LIDIA_REDUCTION_TYPE_H_GUARD_
#define LIDIA_REDUCTION_TYPE_H_GUARD_


#ifndef LIDIA_KODAIRA_CODE_H_GUARD_
# include	"LiDIA/kodaira_code.h"
#endif



#ifdef LIDIA_NAMESPACE
namespace LiDIA {
# define IN_NAMESPACE_LIDIA
#endif



// Class for holding the reduction type at p

class reduction_type
{
	int ord_p_discr;
	int ord_p_N;
	int ord_p_j_denom;
	Kodaira_code Kcode; // NB the constructor makes this from an int
	int c_p;

public:

	reduction_type(int opd = 0, int opN = 0, int opj = 0, int kc = -1, int cp = 0);

	reduction_type & operator = (const reduction_type& r);

	// Access Functions
	int get_ord_p_discriminant() const;
	int get_ord_p_N() const;
	int get_ord_p_j_denom() const;
	int get_c_p()	const;
	Kodaira_code get_code() const;

	friend std::ostream& operator << (std::ostream& os, const reduction_type& R);
};



#ifdef LIDIA_NAMESPACE
}	// end of namespace LiDIA
# undef IN_NAMESPACE_LIDIA
#endif



#endif	// LIDIA_REDUCTION_TYPE_H_GUARD_
