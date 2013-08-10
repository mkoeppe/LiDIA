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



#ifdef LIDIA_NAMESPACE
using namespace LiDIA;
#endif



void help()
{

	std::cout << "\n" << "usage: fpc <prime>\n\n";
	std::cout << "Once started, you can enter the following "
		  << "expressions:\n";
	std::cout << "  q\n";
	std::cout << "  h\n";
	std::cout << "  a+b\n";
	std::cout << "  a-b\n";
	std::cout << "  a*b\n";
	std::cout << "  a/b\n";
	std::cout << "  a^b\n";
	std::cout << "where a, b are either integers or the character '&'.\n";
	std::cout << "'&' is a shorthand for the result of the "
		  << "previous line.\n";
	std::cout << "'q' terminates the program, 'h' prints this help.\n";
	std::cout << std::endl;
}



int main_LiDIA(int argc, char** argv)
{
        typedef std::string::size_type StringIndex;

	std::string inBuffer;
	std::string const prompt("fpc> ");
	std::string const operators("+-*/^");

	bigint p;
	if (argc > 1) {
		string_to_bigint(argv[1], p);
	}
	else {
		help();
		return (0);
	}

	if (!is_prime(p, 6)) {
		std::cerr << "fpc: " << p << " is not a prime number " << "\n";
		return(1);
	}

	std::cout << "Computing in F" << p << " ..." << "\n" << "\n";

	int count = 0;
	bigint y = 0;

	while (true) {
		std::cout << prompt;
		std::getline(std::cin, inBuffer);

		if (inBuffer.empty() || inBuffer[0] == 'h') {
			help();
			continue;
		}

		if (inBuffer[0] == 'q')
			break;

		StringIndex opPos = inBuffer.find_first_of(operators);
// 		op = strpbrk(string, "+-*/^");

//		if (op == 0) {
		if(opPos == std::string::npos) {
			std::cerr << "no operand found (op = '+-*/^')" << "\n";
			continue;
		}
		std::string lhs = inBuffer.substr(0, opPos);
		char opc = inBuffer[opPos];
		std::string rhs = inBuffer.substr(opPos+1);

//		opc = *op;

//		*op = 0;

		bool check1 = !lhs.empty();
		bigint a;
		if (check1 && lhs[0] == '&') {
		    a = y;
		}
		else {
		    check1 = check1 && string_to_bigint(lhs.c_str(), a);
		}
		
		bool check2 = !rhs.empty();
		bigint b;
		if (check2 && rhs[0] == '&') {
		    b = y;
		}
		else {
		    check2 = check2 && string_to_bigint(rhs.c_str(), b);
		}

		if (!check1) {
			std::cerr << "operand 1 not numeric" << "\n";
			continue;
		}

		if (!check2) {
			std::cerr << "operand 2 not numeric" << "\n";
			continue;
		}

		if (a.is_lt_zero()) {
			a %= p;
			a += p;
		}

		if (b.is_lt_zero()) {
			if (opc == '^') {
				b %= (p-bigint(1));
				b += (p-bigint(1));
			}
			else {
				b %= p;
				b += p;
			}
		}

		switch (opc) {
		case '^':
		    power_mod(y, a, b, p);
		    break;
		case '+':
		    y = (a+b) % p;
		    break;
		case '-':
		    y = (a-b) % p;
		    break;
		case '*':
		    y = (a*b) % p;
		    break;
		case '/':
		    bigint inv_b;
		    power_mod(inv_b, b, p-bigint(2), p);
		    y = (a*inv_b) % p;
		}

		std::cout << "out(" << ++count << ") = " << y << "\n";
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
