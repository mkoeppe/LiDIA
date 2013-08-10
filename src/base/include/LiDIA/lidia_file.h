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
//	Author	:
//	Changes	: See CVS log
//
//==============================================================================================


#ifndef LIDIA_LIDIA_FILE_H_GUARD_
#define LIDIA_LIDIA_FILE_H_GUARD_


#ifndef LIDIA_LIDIA_H_GUARD_
# include	"LiDIA/LiDIA.h"
#endif



#ifdef LIDIA_NAMESPACE
namespace LiDIA {
# define IN_NAMESPACE_LIDIA
#endif



class lidia_file
{
private:

	// static residue_class_list < char* > *file_names;
	// static long file_number;

public:

	// static bool get_filename (char* & fn, char* prefix);
	// static void delete_file  (char* fn);
	static bool file_exists  (const char* fn);
};



#ifdef LIDIA_NAMESPACE
}	// end of namespace LiDIA
# undef IN_NAMESPACE_LIDIA
#endif



#endif	// LIDIA_LIDIA_FILE_H_GUARD_
