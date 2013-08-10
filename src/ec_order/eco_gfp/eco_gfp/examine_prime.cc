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
//	Author	: Frank Lehmann, Markus Maurer, Peter Noss, Volker Mueller
//	Changes	: See CVS log
//
//==============================================================================================



#ifdef HAVE_CONFIG_H
# include	"config.h"
#endif
#include	"LiDIA/meq_prime.h"



#ifdef LIDIA_NAMESPACE
namespace LiDIA {
#endif



// modus == 0 : normal
//          1 : XONHK1
//          2 : X2NHK1
//          3 : XONHK2
//          4 : theta - Reihen



bool
meq_prime::examine_prime ()
{
	debug_handler ("meq_prime", "examine_prime ()");

	//
	// For type 1: choose integer s such that
	//	         v = s * (l-1)/12 is an integer.
	//
	//		 The number of coefficients of the corresponding
	//	 	 bivariate polynomial is (l+2) * (v+1).
	//
	// For type 2: choose v for each l.
	//
	//		 The number of coefficients of the corresponding
	//	 	 bivariate polynomial is (l+2) * (v+1).
	//
	// They are stored as
	//
	//  [ a_00 a_01 ... a_0v ]
	//  [ a_10 a_11     a_1v ]
	//  :
	//  [ a_(l+1)0 a_(l+1)1... a_(l+1)v ]
	//

	s = 0;
	v = 0;
	mode = 0;


	// type is set to 2 by default.
	// If l is of type 1, this must be set with a case statement;
	// then s and v are computed at the end of the function.
	//
	type = 2;

	bool found = true;

	switch (l) {
	case 3   :
	case 5   :
	case 7   :
	case 11  :
	case 13  :
	case 17  :
	case 19  :
	case 23  : type = 1;
		break;

		// Valence == 1

	case 29  : type = 1;
		// mode = 1; hecke_op1 = 5; v = 1;
		break;

	case 31  : type = 1;
		// mode = 1; hecke_op1 = 4; v = 1;
		break;

	case 37  : type = 1;
		// mode = 3; hecke_op1 =; v   = 7;
		break;

	case 41  : type = 1;
		// mode = 1; hecke_op1 = 5; v = 1;
		break;

	case 43  : type = 1;
		// mode = 3; hecke_op1 =; v   = 4;
		break;

	case 47  : type = 1;
		// mode = 1; hecke_op1 = 2; v   = 1;
		break;

	case 53  : mode = 1;
		hecke_op1 = 9;
		v   = 2;
		break;

	case 59  : type = 1;
		// mode = 1; hecke_op1 = 3; v   = 1;
		break;

	case 61  : type = 1;
		// mode = 1; hecke_op1 = 13; v   = 2;
		break;

	case 67  : type = 1;
		// mode = 3; hecke_op1 =; v   = 2;
		break;

	case 71 :  mode = 1;
		hecke_op1 = 2;
		v = 1;
		break;

	case 73  : type = 1;
		// mode = 1; hecke_op1 = 37; v = 3;
		break;

	case 79  : type = 1;
		// mode = 1; hecke_op1 = 4; v = 2;
		break;

	case 83  : type = 1;
		// mode = 1; hecke_op1 = 3; v = 2;
		break;

	case 89  : type = 1;
		// mode = 1; hecke_op1 = 9; v = 2;
		break;

	case 97  : type = 1;
		// mode = 1; hecke_op1 = 49; v = 4;
		break;

	case 101 : type = 1;
		// mode = 2; hecke_op1 = 5; hecke_op2 = 9; v = 2;
		break;

	case 103 : mode = 1;
		hecke_op1 = 4;
		v = 3;
		break;

	case 107 : type = 1;
		// mode = 1; hecke_op1 = 3; v = 3;
		break;

	case 109 : type = 1;
		// mode = 1; hecke_op1 = 25; v = 3;
		break;

	case 113 : type = 1;
		// mode = 1; hecke_op1 = 9; v = 4;
		break;

	case 127 : type = 1;
		// mode = 1; hecke_op1 = 4; v = 4;
		break;

	case 131 : type = 1;
		// mode = 2; hecke_op1 = 5; hecke_op2 = 7; v = 2;
		break;

	case 137 : mode = 1;
		hecke_op1 = 9;
		v = 5;
		break;

	case 139 : type = 1;
		// mode = 2; hecke_op1 = 7; hecke_op2 = 13; v = 4;
		break;

	case 149 : type = 1;
		// mode = 2; hecke_op1 = 5; hecke_op2 = 9; v = 4;
		break;

	case 151 : type = 1;
		// mode = 1; hecke_op1 = 4; v = 3;
		break;

	case 157 : type = 1;
		// mode = 1; hecke_op1 = 13; v = 6;
		break;

	case 163 : type = 1;
		// mode = 3; hecke_op1 =; v = 2;
		break;

	case 167 : mode = 1;
		hecke_op1 = 2;
		v = 3;
		break;

	case 173 : type = 1;
		// mode = 1; hecke_op1 = 9; v = 4;
		break;

	case 179 : mode = 1;
		hecke_op1 = 3;
		v = 5;
		break;

	case 181 : type = 1;
		// mode = 1; hecke_op1 = 37; v = 6;
		break;

	case 191 : mode = 1;
		hecke_op1 = 6;
		v = 3;
		break;

	case 193 : type = 1;
		// mode = 1; hecke_op1 = 97; v = 8;
		break;

	case 197 : type = 1;
		// mode = 1; hecke_op1 = 9; v = 6;
		break;

	case 199 : type = 1;
		// mode = 2; hecke_op1 = 13; hecke_op2 = 16; v = 5;
		break;

	case 211 : type = 1;
		// mode = 2; hecke_op1 = 19; hecke_op2 = 25; v = 7;
		break;

	case 223 : type = 1;
		// mode = 1; hecke_op1 = 4; v = 7;
		break;

	case 227 : mode = 1;
		hecke_op1 = 3;
		v = 6;
		break;

	case 229 : type = 1;
		// mode = 1; hecke_op1 = 49; v = 7;
		break;

	case 233 : type = 1;
		// mode = 1; hecke_op1 = 9; v = 7;
		break;

	case 239 : mode = 1;
		hecke_op1 = 2;
		v = 5;
		break;

	case 241 : type = 1;
		// mode = 1; hecke_op1 = 25; v = 6;
		break;

	case 251 : mode = 1;
		hecke_op1 = 49;
		v = 6;
		break;

	case 257 : mode = 1;
		hecke_op1 = 57;
		v = 6;
		break;

	case 263 : mode = 1;
		hecke_op1 = 2;
		v = 5;
		break;

	case 269 : mode = 1;
		hecke_op1 = 81;
		v = 6;
		break;

	case 271 : type = 1;
		// mode = 2; hecke_op1 = 7; hecke_op2 = 22; v = 6;
		break;

	case 277 : type = 1;
		// mode = 1; hecke_op1 = 13; v = 9;
		break;

	case 281 : mode = 1;
		hecke_op1 = 29;
		v = 7;
		break;

	case 283 : type = 1;
		// mode = 2; hecke_op1 = 7; hecke_op2 = 13; v = 8;
		break;

	case 293 : mode = 1;
		hecke_op1 = 33;
		v = 7;
		break;

	case 307 : mode = 2;
		hecke_op1 = 7;
		hecke_op2 = 19;
		v = 10;
		break;

	case 311 : mode = 1;
		hecke_op1 = 243;
		v = 5;
		break;

        case 313 : mode = 1;
		hecke_op1 = 13;
		v = 12;
		break;

        case 317 : mode = 1;
		hecke_op1 = 9;
		v = 11;
		break;

        case 331 : mode = 1;
		hecke_op1 = 67;
		v = 11;
		break;

        case 337 : mode = 1;
		hecke_op1 = 13;
		v = 13;
		break;

        case 347 : mode = 1;
		hecke_op1 = 3;
		v = 9;
		break;

        case 349 : mode = 1;
		hecke_op1 = 49;
		v = 10;
		break;

        case 353 : mode = 1;
		hecke_op1 = 81;
		v = 8;
		break;

        case 359 : mode = 1;
		hecke_op1 = 153;
		v = 6;
		break;

        case 367 : mode = 1;
		hecke_op1 = 4;
		v = 11;
		break;

	case 373 : type = 1;
		break;

        case 379 : mode = 1;
		hecke_op1 = 79;
		v = 12;
		break;

        case 383 : mode = 1;
		hecke_op1 = 2;
		v = 8;
		break;

        case 389 : mode = 1;
		hecke_op1 = 9;
		v = 12;
		break;

        case 397 : type = 1;
		break;

        case 401 : mode =  1;
		hecke_op1 = 41;
		v = 10;
		break;

        case 409 : mode =  1;
		hecke_op1 = 49;
		v = 10;
		break;

        case 419 : mode =  1;
		hecke_op1 =  3;
		v = 11;
		break;

        case 421 : mode =  1;
		hecke_op1 = 97;
		v = 14;
		break;

        case 431 : mode =  1;
		hecke_op1 =  2;
		v =  9;
		break;

        case 433 : mode =  1;
		hecke_op1 = 13;
		v = 15;
		break;

        case 439 : mode =  1;
		hecke_op1 = 112;
		v = 11;
		break;

        case 443 : mode =  1;
		hecke_op1 =  3;
		v = 12;
		break;

        case 449 : mode =  1;
		hecke_op1 = 53;
		v = 11;
		break;

        case 457 : mode =  1;
		hecke_op1 = 73;
		v = 13;
		break;

        case 461 : mode =  1;
		hecke_op1 = 49;
		v = 13;
		break;

        case 463 : mode =  1;
		hecke_op1 =  4;
		v = 14;
		break;

        case 467 : mode =  1;
		hecke_op1 = 27;
		v = 11;
		break;

        case 479 : mode =  1;
		hecke_op1 =  2;
		v = 10;
		break;

        case 487 : mode =  1;
		hecke_op1 =  4;
		v = 15;
		break;

#if 0
	case 491 : mode =   1;
		hecke_op1 = 243;
		v =   8;
		break;
#endif

        case 491 : mode =   2;
		hecke_op1 =   3;
		hecke_op2 =  33;
		v =   8;
		break;

        case 499 : mode =   1;
		hecke_op1 = 103;
		v =  16;
		break;

        case 503 : mode =   1;
		hecke_op1 =   2;
		v =  10;
		break;

        case 509 : mode =   1;
		hecke_op1 =  37;
		v =  15;
		break;

        case 521 : mode =   1;
		hecke_op1 =  53;
		v =  13;
		break;

        case 523 : mode =   1;
		hecke_op1 =  73;
		v =  15;
		break;

        case 541 : mode =    2;
		hecke_op1 =   25;
		hecke_op2 =  109;
		v =   17;
		break;

        case 547 : mode =   1;
		hecke_op1 =  73;
		v =  18;
		break;

        case 557 : mode =   1;
		hecke_op1 =  69;
		v =  16;
		break;

        case 563 : mode =   1;
		hecke_op1 =  27;
		v =  13;
		break;

        case 569 : mode =    1;
		hecke_op1 =  101;
		v =   14;
		break;

        case 571 : mode =    1;
		hecke_op1 =  127;
		v =   19;
		break;

        case 577 : mode =   1;
		hecke_op1 =  49;
		v =  20;
		break;

        case 587 : mode =   1;
		hecke_op1 =   3;
		v =  16;
		break;

        case 593 : mode =    1;
		hecke_op1 =  153;
		v =   14;
		break;

        case 599 : mode =   1;
		hecke_op1 =   2;
		v =  12;
		break;

        case 601 : type = 1;
		// mode = 1; hecke_op1 =  61; v =  15;
		break;

        case 607 : mode =   1;
		hecke_op1 =  52;
		v =  18;
		break;

        case 613 : mode =   1;
		hecke_op1 =  61;
		v =  19;
		break;

        case 617 : mode =   1;
		hecke_op1 =   9;
		v =  20;
		break;

        case 619 : mode =    2;
		hecke_op1 =   61;
		hecke_op2 =  127;
		v =   17;
		break;

        case 631 : mode =   1;
		hecke_op1 =  16;
		v =  15;
		break;

        case 641 : mode =    1;
		hecke_op1 =  137;
		v =   16;
		break;

        case 643 : mode =   1;
		hecke_op1 =  49;
		v =  22;
		break;

        case 647 : mode =   1;
		hecke_op1 =   2;
		v =  13;
		break;

        case 653 : mode =   1;
		hecke_op1 =   9;
		v =  24;
		break;

        case 659 : mode =   1;
		hecke_op1 =   3;
		v =  18;
		break;

        case 661 : mode =   1;
		hecke_op1 = 157;
		v =  21;
		break;

        case 673 : mode =   1;
		hecke_op1 =  13;
		v =  25;
		// not available
		break;

        case 677 : mode =   1;
		hecke_op1 =  81;
		v =  16;
		break;

        case 683 : mode =   1;
		hecke_op1 =   3;
		v =  19;
		break;

        case 691 : mode =   1;
		hecke_op1 =  79;
		v =  21;
		break;

        case 701 : mode =   2;
		hecke_op1 = 101;
		hecke_op2 = 153;
		v =  21;
		break;

        case 709 : mode =   1;
		hecke_op1 = 181;
		v =  23;
		break;

        case 719 : mode =   1;
		hecke_op1 =   2;
		v =  15;
		break;

        case 727 : mode =   1;
		hecke_op1 =   4;
		v =  21;
		break;

        case 733 : mode =   1;
		hecke_op1 =  37;
		v =  25;
		break;

        case 739 : mode =   1;
		hecke_op1 = 151;
		v =  24;
		break;

        case 743 : mode =   1;
		hecke_op1 =   2;
		v =  15;
		break;

        case 751 : mode =   1;
		hecke_op1 = 208;
		v =  18;
		break;

        case 757 : mode =   1;
		hecke_op1 =  73;
		v =  26;
		break;

        case 761 : mode =   1;
		hecke_op1 = 489;
		v =  23;
		break;

        case 769 : mode =   1;
		hecke_op1 = 121;
		v =  19;
		break;

        case 773 : mode =   1;
		hecke_op1 = 513;
		v =  18;
		break;

        case 787 : mode =   1;
		hecke_op1 =  67;
		v =  26;
		break;

        case 797 : mode =   1;
		hecke_op1 = 621;
		v =  19;
		break;

        case 809 : mode =   1;
		hecke_op1 =  89;
		v =  20;
		break;

        case 811 : mode =   1;
		hecke_op1 = 163;
		v =  27;
		break;

        case 821 : mode =   1;
		hecke_op1 = 405;
		v =  24;
		break;

        case 823 : mode =   1;
		hecke_op1 =   4;
		v =  24;
		break;

        case 827 : mode =   1;
		hecke_op1 =   3;
		v =  23;
		break;

        case 829 : mode =   1;
		hecke_op1 = 169;
		v =  27;
		break;

        case 839 : mode =   1;
		hecke_op1 = 243;
		v =  14;
		break;

        case 853 : mode =   1;
		hecke_op1 =  49;
		v =  29;
		// not available
		break;

        case 857 : mode =   1;
		hecke_op1 = 117;
		v =  20;
		break;

        case 859 : mode =   1;
		hecke_op1 = 103;
		v =  27;
		// not available
		break;

        case 863 : mode =   1;
		hecke_op1 =   2;
		v =  18;
		break;

        case 877 : mode =   1;
		hecke_op1 =  73;
		v =  28;
		break;

        case 881 : mode =   1;
		hecke_op1 =  89;
		v =  22;
		break;

        case 883 : mode =   1;
		hecke_op1 =  67;
		v =  32;
		// not available
		break;

        case 887 : mode =   1;
		hecke_op1 =   2;
		v =  18;
		break;

        case 911 : mode =   1;
		hecke_op1 = 243;
		v =  15;
		break;

        case 929 : mode =   1;
		hecke_op1 =  101;
		v =   23;
		break;

        case 971 : mode =   1;
		hecke_op1 =  515;
		v =   23;
		// not available
		break;

        case 983 : mode =   1;
		hecke_op1 =   2;
		v =  20;
		break;

        case 1013 : mode =   1;
		hecke_op1 = 837;
		v =  24;
		break;

        case 1019 : mode =   1;
		hecke_op1 =  81;
		v =  17;
		break;

	default: found = false;
		break;
	} // switch


	//
	// Determine s and v for type 1 primes.
	//
	if (type == 1) {
		s = 1;
		v = l-1;

		while (v % 12 != 0) {
			v += l-1;
			s += 1;
		}

		v /= 12;
		mode = 0;
	}

	//
	// For type 2 primes double v, such that v+1
	// is the number of cols in the coefficient
	// array of the "modular" polynomial.
	//
	if (type == 2)
		v *= 2;

#ifdef AIX_DEBUG
	std::cout << "meq_prime::examine_prime:: l = " << l << std::endl;
	std::cout << "meq_prime::examine_prime:: type = " << type << std::endl;
	std::cout << "meq_prime::examine_prime:: mode = " << mode << std::endl;
	std::cout << "meq_prime::examine_prime:: v = " << v << std::endl;
	std::cout << "meq_prime::examine_prime:: s = " << s << std::endl;
#endif

	return found;
}



#ifdef LIDIA_NAMESPACE
}	// end of namespace LiDIA
#endif
