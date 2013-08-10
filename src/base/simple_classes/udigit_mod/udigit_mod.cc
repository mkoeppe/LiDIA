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
//	Author	: Thorsten Rottschaefer
//	Changes	: See CVS log
//
//==============================================================================================


#ifdef HAVE_CONFIG_H
# include	"config.h"
#endif
#include	"LiDIA/udigit_mod.h"



#ifdef LIDIA_NAMESPACE
namespace LiDIA {
#endif



udigit udigit_mod::prime = 0;
udigit udigit_mod::h = 0;



#ifdef LIDIA_NAMESPACE
}	// end of namespace LiDIA
#endif
