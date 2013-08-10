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
//	Changes	: See CVS log
//
//==============================================================================================


#include	"LiDIA/quartic.h"



#ifdef LIDIA_NAMESPACE
using namespace LiDIA;
#endif



int main_LiDIA(int argc, char** argv)
{
	long lidia_precision = 40;
	bigfloat::set_precision(lidia_precision);

	// First a test of quartics...
	quartic q(-1, 0, 54, -256, 487);
	std::cout << q << std::endl;
	std::cout << "Testing " << std::endl;
	std::cout << "2 : " << q.qp_soluble(2) << std::endl;
	std::cout << "3 : " << q.qp_soluble(3) << std::endl;
	std::cout << "11 : " << q.qp_soluble(11) << std::endl;
	std::cout << "941 : " << q.qp_soluble(941) << std::endl;

	q.assign(1, 2, -9, -42, -43);
	std::cout << q << std::endl;
	std::cout << "Testing " << std::endl;
	std::cout << "2 : " << q.qp_soluble(2) << std::endl;
	std::cout << "3 : " << q.qp_soluble(3) << std::endl;
	std::cout << "11 : " << q.qp_soluble(11) << std::endl;
	std::cout << "941 : " << q.qp_soluble(941) << std::endl;

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
