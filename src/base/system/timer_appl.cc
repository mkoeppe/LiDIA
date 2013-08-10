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
//	Author	: Thomas Papanikolaou (TP)
//	Changes	: See CVS log
//
//==============================================================================================


#include	"LiDIA/timer.h"



#ifdef LIDIA_NAMESPACE
using namespace LiDIA;
#endif



int main_LiDIA(int argc, char** argv)
{
	timer x;
	long i;

	x.set_print_mode(HMS_MODE); // choose h m s ms

	x.start_timer();
#if defined(__hppa)
	for (i = 0; i < 0x7ffffff; i++);
#else
	for (i = 0; i < 0x7fffff; i++);
#endif
	x.stop_timer();
	std::cout << std::endl;
	x.print();
	std::cout << std::endl;

	x.set_print_mode(TIME_MODE); // choose time mode

	x.cont_timer();
#if defined(__hppa)
	for (i = 0; i < 0x7ffffff; i++);
#else
	for (i = 0; i < 0x7fffff; i++);
#endif
	x.stop_timer();
	std::cout << std::endl;
	x.print();
	std::cout << std::endl;

	return (0);
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
