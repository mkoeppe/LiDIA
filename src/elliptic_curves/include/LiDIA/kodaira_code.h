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


#ifndef LIDIA_KODAIRA_CODE_H_GUARD_
#define LIDIA_KODAIRA_CODE_H_GUARD_


#ifndef LIDIA_LIDIA_H_GUARD_
# include	"LiDIA/LiDIA.h"
#endif



#ifdef LIDIA_NAMESPACE
namespace LiDIA {
# define IN_NAMESPACE_LIDIA
#endif



// class Kodaira_code just holds an int which "codes" a
// Kodaira code as follows:
// (this coding originally from R.G.E.Pinch)
//
// Im                 ->10*m
// I*m                ->10*m+1
// I, II, III, IV     ->1, 2, 3, 4
// I*, II*. III*, IV* ->5, 6, 7, 8
//

class Kodaira_code
{
	int code;

public:

	Kodaira_code(int k = -1);
	Kodaira_code(const Kodaira_code &);
	~Kodaira_code();

	Kodaira_code & operator = (int k);
	Kodaira_code & operator = (const Kodaira_code& c);

	// Access Functions
	int get_code() const;

	// output
	friend std::ostream& operator << (std::ostream& os, const Kodaira_code& c);
};



#ifdef LIDIA_NAMESPACE
}	// end of namespace LiDIA
# undef IN_NAMESPACE_LIDIA
#endif



#endif	// LIDIA_KODAIRA_CODE_H_GUARD_
