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
//	Author	: Patrick Theobald (PT)
//	Changes	: See CVS log
//
//==============================================================================================


#include	"LiDIA/matrix/crt_and_prime_handling.h"
#include	"LiDIA/prime_list.h"
#include	"LiDIA/modular_operations.inl"
#include        "LiDIA/isstream.h"
#include	"LiDIA/path.h"

#include        <vector>
#include        <string>
#include	<cstdlib>
#include        <cstring>   // declares strcpy()


#ifdef LIDIA_NAMESPACE
namespace LiDIA {
#endif



//
// debug defines / error defines
//

    extern const char *PRT;
    extern const char *matrix_error_msg[];

#define DMESSAGE "crt_and_prime_handling"
#define DVALUE   300
#define EMESSAGE matrix_error_msg

#ifndef LIDIA_PRIMELIST_SIZE
# define LIDIA_PRIMELIST_SIZE 300
#endif
#ifndef LIDIA_PRIMELIST_MAXSIZE
# define LIDIA_PRIMELIST_MAXSIZE 86416
#endif



//
// CRT
//

    bigint
    chinrest(const bigint *values, const bigint *prim)
    {
	bigint RES;
	debug_handler_l(DMESSAGE, "in function "
			"chinrest(bigint &, const bigint *, const bigint *)", DVALUE);

	register lidia_size_t i;
	bigint M, mod;
	bigint TMP, TMP0, TMP1, TMP2;

	long len;
	prim[0].longify(len);

	bigint e;

	// Step 1
	bigint L = 1;
	for (i = 0; i < len; i++)
	    LiDIA::multiply(L, L, prim[i + 1]);

	bigint X;
	for (i = 0; i < len; i++) {
	    mod.assign(prim[i + 1]);
	    div_rem(M, TMP, L, mod);
	    best_remainder(TMP, M, mod);
	    xgcd_right(TMP1, mod, TMP);
	    LiDIA::multiply(e, TMP1, M);
	    LiDIA::multiply(TMP, e, values[i]);
	    LiDIA::add(X, X, TMP);
	}

	best_remainder(X, X, L);
	LiDIA::subtract(TMP0, L, bigint(1));
	shift_left(TMP1, X, 1);
	shift_left(TMP2, L, 1);

	if (TMP1 >= -TMP0 && TMP1 <= TMP0)
	    RES.assign(X);
	else if (TMP1 - TMP2 <= TMP0 || TMP1 - TMP2 >= -TMP0)
	    LiDIA::subtract(RES, X, L);
	else
	    lidia_error_handler("void chinrest(bigint & RES, const bigint * values, const bigint * prim)",
				DMESSAGE, EMESSAGE[8]);
	return RES;
    }



    void
    chinrest(bigint & RES, const bigint * values, const bigint * prim)
    {
	//
	// DESCRIPTION: chinrest(Res, v, prim);
	// =  > Res = solution of the Chinese remaindering theorem
	//                 with parameters v and prim
	// VERSION: 2.0
	//

	debug_handler_l(DMESSAGE, "in function "
			"chinrest(bigint &, const bigint *, const bigint *)", DVALUE);

	register lidia_size_t i;
	bigint M, mod;
	bigint TMP, TMP0, TMP1, TMP2;

	long len;
	prim[0].longify(len);

	bigint e;

	// Step 1
	bigint L = 1;
	for (i = 0; i < len; i++)
	    LiDIA::multiply(L, L, prim[i + 1]);

	bigint X;
	for (i = 0; i < len; i++) {
	    mod.assign(prim[i + 1]);
	    div_rem(M, TMP, L, mod);
	    best_remainder(TMP, M, mod);
	    xgcd_right(TMP1, mod, TMP);
	    LiDIA::multiply(e, TMP1, M);
	    LiDIA::multiply(TMP, e, values[i]);
	    LiDIA::add(X, X, TMP);
	}

	best_remainder(X, X, L);
	LiDIA::subtract(TMP0, L, bigint(1));
	shift_left(TMP1, X, 1);
	shift_left(TMP2, L, 1);

	if (TMP1 >= -TMP0 && TMP1 <= TMP0)
	    RES.assign(X);
	else if (TMP1 - TMP2 <= TMP0 || TMP1 - TMP2 >= -TMP0)
	    LiDIA::subtract(RES, X, L);
	else
	    lidia_error_handler("void chinrest(bigint & RES, const bigint * values, const bigint * prim)",
				DMESSAGE, EMESSAGE[8]);
    }



//
// prime handling
//

    static void
    load_primes(std::string const& filename, lidia_size_t number_of_primes,
		std::vector<bigint>& list) {
	list.resize(0);
	list.reserve(number_of_primes + 1);
        list.push_back(1);         // list[0] is a dummy

	prime_list pl;
	pl.load_from_file(filename.c_str(), number_of_primes);

	PRIME_LIST_NUMBER p = pl.get_last_prime();
	for (lidia_size_t i = 1; i <= number_of_primes; ++i) {
	    list.push_back(p);
	    p = pl.get_prev_prime();
	}
    }



    bigint *
    get_primes(const bigint & C, const bigint & m, const bool SW_COPY) {
	//
	// DESCRIPTION: get_primes(C, m) = v;
	// =  > v = List of primes
	// =  > !(m | v[i]) for all i = 1, ..., v[0]
	// =  > v[1]*...*v[v[0]] > C
	// VERSION: 2.0
	//

	debug_handler_l(DMESSAGE, "in function "
			"get_primes(bigint, bigint)", DVALUE);

	static std::vector<bigint> static_list;
	static lidia_size_t bitlength = 0;

	// GET ENV Variable
	static lidia_size_t anzahl = 0;
	static std::string lidia_primes_name = "";

	if (anzahl == 0) {
	    char* lidia_primes_size = getenv("LIDIA_PRIMES_SIZE");
	    if (lidia_primes_size == NULL) {
		anzahl = LIDIA_PRIMELIST_SIZE;
	    }
	    else {
		isstream is(lidia_primes_size);
		is >> anzahl;
		if(anzahl <= 0 || is.fail()) {
		    lidia_error_handler("get_primes()",
					"Error while evaluating "
					"environment variable "
					"LIDIA_PRIMES_SIZE");
		}
	    }
	}

	if (lidia_primes_name.length() == 0) {
	    char* primes_name_env = getenv("LIDIA_PRIMES_NAME");
	    if (primes_name_env == NULL || primes_name_env[0] == 0) {
		lidia_primes_name = LIDIA_PRIMELIST_NAME;
	    }
	    else {
		lidia_primes_name = primes_name_env;
	    }
	}

	if (static_list.size() == 0) {
	    debug_handler_c(DMESSAGE, "in function "
			    "get_primes(bigint, bigint)", DVALUE,
			    std::cout << " Auffuellen der Primzahlliste " << std::endl;);
	    load_primes(lidia_primes_name, anzahl, static_list);
	    for (lidia_size_t i = 1; i <= anzahl; i++) {
		bitlength += static_list[i].bit_length();
	    }
	}

	lidia_size_t PRObits = 0;
	lidia_size_t count = 0;
	lidia_size_t Hbits = C.bit_length();
	bool SW = true;
	do {
	    if (bitlength >= Hbits) {
		debug_handler_c(DMESSAGE, "in function "
				"get_primes(bigint, bigint)", DVALUE,
				std::cout << " Test der Primzahlliste " << std::endl;);

		for (lidia_size_t i = 1; i <= anzahl && PRObits <= Hbits; i++) {
		    if (m % static_list[i] != 0) {
			PRObits += static_list[i].bit_length();
			count++;
		    }
		}
	    }
	    else {
		PRObits = -1;
	    }

	    if (PRObits < Hbits) {
		debug_handler_c(DMESSAGE, "in function "
				"get_primes(bigint, bigint)",
				DVALUE,
				std::cout << " Erweitern der Primzahlliste "
				          << std::endl;);
		// fix anzahl
		if (anzahl == LIDIA_PRIMELIST_MAXSIZE) {
		    lidia_error_handler(DMESSAGE,
					"get_primes :: Primlist too small !!");
		}
		anzahl *= 2;
		if (anzahl > LIDIA_PRIMELIST_MAXSIZE) {
		    anzahl = LIDIA_PRIMELIST_MAXSIZE;
		}

		std::cout << DMESSAGE << ":: get_primes :: "
			  << "Size of primelist increased !\n"
			  << "new size: " << anzahl << std::endl;

		load_primes(lidia_primes_name, anzahl, static_list);
		bitlength = 0;
		for (lidia_size_t i = 1; i <= anzahl; i++) {
		    bitlength += static_list[i].bit_length();
		}
	    }
	    else {
		SW = false;
	    }
	} while (SW);

	static_list[0].assign(count);
	debug_handler_c(DMESSAGE, "in function "
			"get_primes(bigint, bigint)", DVALUE,
			std::cout << "Meldung aus get_primes:\n";
			std::cout << "=======================\n";
			std::cout << "Bitlaenge = " << bitlength << '\n';
			std::cout << "Hbits = " << Hbits << '\n';
			std::cout << "Anzahl der Primzahlen = "
			<< static_list[0] << std::endl;);

	if (SW_COPY) {
	    bigint *PRIM = new bigint[count + 1];
	    for (lidia_size_t i = 1, j = 1; i <= count ; j++)
		if (m % static_list[j] != 0) {
		    PRIM[i] = static_list[j];
		    i++;
		}
	    PRIM[0] = count;
	    return PRIM;
	}

	return &static_list[0];
    }



#undef DMESSAGE
#undef DVALUE
#undef EMESSAGE



#ifdef LIDIA_NAMESPACE
}	// end of namespace LiDIA
#endif
