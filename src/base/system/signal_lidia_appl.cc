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
//	Author	: Markus Maurer (MM)
//	Changes	: See CVS log
//
//==============================================================================================


#include	"LiDIA/lidia_signal.h"
#include	<unistd.h>



#ifdef LIDIA_NAMESPACE
using namespace LiDIA;
#endif



LIDIA_SIGNAL_FUNCTION (myhandler)
{
	std::cout << " Hello, this is the myhandler function.\n";
	std::cout << " I have caught your LIDIA_SIGINT.\n";
	std::cout << std::endl;
}



void
foo ()
{
	// install myhandler for LIDIA_SIGINT
	lidia_signal s(LIDIA_SIGINT, myhandler);

	std::cout << "\n";
	std::cout << " If you type Ctrl-C within the next 15 seconds, \n";
	std::cout << " the function myhandler is called.\n";
	std::cout << " sleep(15) ...\n" << std::endl;

#if defined (_MSC_VER)
	while (true)
		;
#else
	sleep(15);
#endif


	// The automatic call of the destructor reinstalls
	// the previous handler for LIDIA_SIGINT.
}



int main_LiDIA (int, char**)
{
	std::cout << "\n";
	std::cout << "Going into function foo ...\n";

	foo();

	std::cout << "Back from foo.\n";
	std::cout << "The default handler for LIDIA_SIGINT is installed again.\n";
	std::cout << "sleep(15) ...\n" << std::endl;

#if defined (_MSC_VER)
	while (true)
		;
#else
	sleep(15);
#endif

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
