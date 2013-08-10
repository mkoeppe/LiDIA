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
//	Author	: Nigel Smart (NS)
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



#ifdef LIDIA_NAMESPACE
namespace LiDIA {
#endif



Kodaira_code
::Kodaira_code(int k) : code(k)
{
	debug_handler("Kodaira_code", "Kodaira_code(int k = -1)");
}



Kodaira_code
::Kodaira_code(const Kodaira_code & k)
{
	debug_handler("Kodaira_code", "Kodaira_code(const Kodaira_code &)");
	code = k.code;
}



Kodaira_code
::~Kodaira_code()
{
	debug_handler("Kodaira_code", "~Kodaira_code()");
}



Kodaira_code & Kodaira_code
::operator = (int k)
{
	debug_handler("Kodaira_code", "operator = (int)");
	code = k;
	return *this;
}



Kodaira_code & Kodaira_code
::operator = (const Kodaira_code& c)
{
	debug_handler("Kodaira_code", "operator = (const Kodaira_code&)");
	code = c.code;
	return *this;
}



//
// Access Functions
//
int Kodaira_code
::get_code() const
{
	debug_handler("Kodaira_code", "get_code() const");
	return code;
}



//
// output
//

std::ostream& operator << (std::ostream& os, const Kodaira_code& c)
{
	int code = c.code;
	switch (code%10) {
	case 0:
		os << "I" << ((code)/10); break;
	case 1:
		os << "I*" << ((code - 1)/10); break;
	case 2:
		os << "II   "; break;
	case 3:
		os << "III  "; break;
	case 4:
		os << "IV   "; break;
	case 5:
		os << "IV*  "; break;
	case 6:
		os << "III* "; break;
	case 7:
		os << "II*  "; break;
	default:
		os << "???? "; break;
	};
	return os;
}



#ifdef LIDIA_NAMESPACE
}	// end of namespace LiDIA
#endif
