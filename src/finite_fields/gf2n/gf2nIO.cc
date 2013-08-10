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


#ifdef HAVE_CONFIG_H
# include	"config.h"
#endif
#include	"LiDIA/gf2n.h"
#include	"LiDIA/finite_fields/gf2nIO.h"



#ifdef LIDIA_NAMESPACE
namespace LiDIA {
#endif



//
// STATIC VARIABLES
// ----------------------------------------------------------------------

gf2nIO::base      gf2nIO::iobase;
char*             gf2nIO::ioprefix;



//
// CONSTRUCTORS
// ----------------------------------------------------------------------

gf2nIO::gf2nIO()
{
	gf2nIO::iobase = Dec;
	delete[] gf2nIO::ioprefix;
	gf2nIO::ioprefix = new char[5];
	strcpy(gf2nIO::ioprefix, "Dec");
}



gf2nIO::gf2nIO(gf2nIO::base base, char* prefix)
{
	gf2nIO::iobase = base;
	delete[] gf2nIO::ioprefix;
	gf2nIO::ioprefix = new char [strlen(prefix) + 1];
	strcpy(gf2nIO::ioprefix, prefix);
}



//
// MEMBERFUNCTIONS
// ----------------------------------------------------------------------

void gf2nIO::setbase(gf2nIO::base base)
{
	if ((base == Dec) || (base == Hex)) {
		delete[] gf2nIO::ioprefix;
		gf2nIO::ioprefix = new char[5];

		if (base == Dec)
			strcpy(gf2nIO::ioprefix, "Dec");

		if (base == Hex)
			strcpy(gf2nIO::ioprefix, "Hex");

		gf2nIO::iobase = base;
	}
	else {
		lidia_error_handler("gf2nIO", "setbase::nonexistent base");
		return;
	}
}



// ---------------------------------------------------------------------

gf2nIO::base gf2nIO::showbase()
{
	return gf2nIO::iobase;
}



// ---------------------------------------------------------------------

void gf2nIO::setprefix(char* prefix)
{
	delete[] gf2nIO::ioprefix;
	gf2nIO::ioprefix = new char [strlen(prefix) + 2];
	strcpy(gf2nIO::ioprefix, prefix);
}



// ---------------------------------------------------------------------

void gf2nIO::setprefix(gf2nIO::base base)
{
	switch(base) {
	case Dec :
		delete[] gf2nIO::ioprefix;
		gf2nIO::ioprefix = new char[5];
		strcpy(gf2nIO::ioprefix, "Dec");
		break;

	case Hex :
		delete[] gf2nIO::ioprefix;
		gf2nIO::ioprefix = new char[5];
		strcpy(gf2nIO::ioprefix, "Hex");
		break;

	default: lidia_error_handler("gf2nIO", "setprefix::nonexistent base");
		break;
	};
}



// ---------------------------------------------------------------------

void gf2nIO::noprefix()
{
	delete[]  gf2nIO::ioprefix;
	gf2nIO::ioprefix = NULL;
}



// ---------------------------------------------------------------------

char* gf2nIO::showprefix()
{
	return gf2nIO::ioprefix;
}



#ifdef LIDIA_NAMESPACE
}	// end of namespace LiDIA
#endif
