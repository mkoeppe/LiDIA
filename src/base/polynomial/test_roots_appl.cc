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
#include	"LiDIA/polynomial.h"
#include	"LiDIA/bigcomplex.h"



#ifdef LIDIA_NAMESPACE
using namespace LiDIA;
#endif



void
shell(bigcomplex * v, int n)
{
	int gap, i, j, c;
	bigcomplex temp;

	for (gap = n / 2; gap > 0; gap /= 2)
		for (i = gap; i < n; i++) {
			c = (abs(v[i - gap]) > abs(v[i]));
			for (j = i - gap; j >= 0 && c; j -= gap) {
				c = (abs(v[j]) > abs(v[j + gap]));
				temp = v[j];
				v[j] = v[j + gap];
				v[j + gap] = temp;
			}
		}
}



int main_LiDIA(int argc, char** argv)
{

	polynomial< bigcomplex > Q;
	bigcomplex x, *rt;
	timer T;
	int i = 0;
	char ch = 'y';

	while (ch == 'y') {
		std::cout << "\f\n\n";
		std::cout << " A root finding algorithm over C\n";
		std::cout << " -------------------------------\n\n";
		std::cout << " Solve for x in C          \n\n";
		std::cout << "             n               \n";
		std::cout << "           _____             \n";
		std::cout << "           \\                \n";
		std::cout << "            \\        i      \n";
		std::cout << " f(x) = )   a  x = 0 \n";
		std::cout << "            /     i          \n";
		std::cout << "           /____             \n";
		std::cout << "           i = 0             \n";
		std::cout << "                             \n";
		std::cout << " Enter the working precision : = " << std::flush;
		std::cin >> i;
		bigfloat::set_precision(i);
		std::cout << "\n";
		std::cout << " Enter a polynomial\n\n";
		std::cout << " f(x) : = " << std::flush;
		std::cin >> Q;
		rt = new bigcomplex[Q.degree()];
		std::cout << "\n The roots of f(x) : = " << Q;
		std::cout << "\n (in a random order) are :" << std::endl;
		T.start_timer();
		roots(Q, rt);
		T.stop_timer();
		std::cout << T;
		shell(rt, Q.degree());
		for (i = 0; i < Q.degree(); i++) {
			std::cout << "\n %" << i + 1 << ": " << rt[i] << std::endl;
			std::cout << " eval f(%" << i + 1 << ") : " << Q(rt[i]) << std::endl;
		}

		delete[] rt;
		std::cout << "\n\n One more time (y/n) ? : " << std::endl;
		std::cin >> ch;
	}

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
