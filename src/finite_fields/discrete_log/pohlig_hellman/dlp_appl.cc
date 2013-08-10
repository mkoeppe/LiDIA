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
//	Author	: Damian Weber (DW)
//	Changes	: See CVS log
//
//==============================================================================================


#include	"LiDIA/dlp.h"
#include	"LiDIA/timer.h"

#include        <cstdlib>



#ifdef LIDIA_NAMESPACE
using namespace LiDIA;
#endif



#define PROGNAME "dlp_appl"



void help()
{
	std::cerr << "Solve the Discrete Logarithm Problem a^x = b mod p, p prime" << "\n";
	std::cerr << "\n" << "usage: dlp_appl a b p" << "\n" << "\n";
	exit(0);
}



int  flag_set(int argc, char **argv, char *s)
{
	int j;

	for (j = 1; j <= argc-1; j++)
		if (!strcmp(argv[j], s))
			return j;
	return 0;
}



int main_LiDIA(int argc, char** argv)
{
	bigint a, b, p;
	rational_factorization f;

	if (flag_set(argc, argv, "-h") || argc != 4)
		help();

	if (!string_to_bigint(argv[1], a)) help();
	if (!string_to_bigint(argv[2], b)) help();
	if (!string_to_bigint(argv[3], p)) help();

	dl(a, b, p, 1);

	return 0;
}


int main(int argc, char** argv) {

#if defined(LIDIA_EXCEPTIONS)
    try {
#endif

	main_LiDIA(argc, argv);
	
#if defined(LIDIA_EXCEPTIONS)
    }
    catch(basic_error const& ex) {
	ex.traditional_error_handler();
	return 1;
    }
    catch(std::exception const& ex) {
	std::cerr << "unexpected exception: " << ex.what() << "\n";
	return 1;
    }
#endif
    
}
