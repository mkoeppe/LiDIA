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
//	Author	: Markus Maurer, Volker Mueller
//	Changes	: See CVS log
//
//==============================================================================================


#ifndef LIDIA_MEQ_PRIME_H_GUARD_
#define LIDIA_MEQ_PRIME_H_GUARD_



#ifndef LIDIA_FP_POLYNOMIAL_H_GUARD_
# include	"LiDIA/Fp_polynomial.h"
#endif
#ifndef LIDIA_UDIGIT_H_GUARD_
# include	"LiDIA/udigit.h"
#endif
#ifndef LIDIA_SORT_VECTOR_H_GUARD_
# include	"LiDIA/sort_vector.h"
#endif
#include        <string>
#include        <vector>



#ifdef LIDIA_NAMESPACE
namespace LiDIA {
# define IN_NAMESPACE_LIDIA
#endif



typedef udigit  meq_prime_t;
typedef int     meq_size_t;
#define MEQ_PREFIX "meq"         // prefix for name of data files



class bigmod;



class meq_prime
{
public:

	static const char ASCII;
	static const char BINARY;
	static const char GZIP;
	static const char BIN_GZIP;
	static const char BZIP;
	static const char BIN_BZIP;

	static const udigit MAX_MEQ_PRIME; // the maximal possible prime

	meq_prime_t l; // prime number
	meq_prime_t v; // valence of modular equation
	meq_prime_t s; // s == s(l, v)

	char        temp_file[1024]; // file name for .gz files
	int         type; // type of modular function
	int         mode; // mode of computation
	meq_prime_t hecke_op1, hecke_op2; // Hecke-Operators for computation

	char       file_format; // format of individual file
	FILE      *fp_me; // file-pointer to individual modular equation
	meq_size_t me_row; // index of row, fp_me is currently pointing to

	base_vector< bigint >  row; // vector to store coeff. of modular equations


private:

	//  *****  private member functions  *****

	bool examine_prime (); // determine necess. information about a prime.
	// returns true if information are available, false otherwise.

	//  *****  public  member functions  *****

public :

	meq_prime ();
	meq_prime (const meq_prime & p);
	~meq_prime ();

	meq_prime & operator = (const meq_prime & p);
	meq_prime & operator = (meq_prime_t p);

	bool set_prime (meq_prime_t p);

	meq_prime_t get_prime  () const
	{
		debug_handler ("meq_prime", "get_prime()");
		return l;
	}

	meq_prime_t valence    ();
	int         trans_type ();


private:
        static std::string const bin_suffix;
        static std::string const gz_suffix;
        static std::string const bin_gz_suffix;
        static std::string const bz2_suffix;
        static std::string const bin_bz2_suffix;

	static int fexist  (char const* fname);
        void get_filename_vector(std::vector<std::string>& fn_vec);
	bool meq_file_is_available();

	int ropen_meq_file ();
	int close_meq_file ();

	int read_row_ascii (base_vector< bigint > & vec);
	int read_row_ascii (base_vector< bigint > & vec, const bigint & p);

	int read_row_bin   (base_vector< bigint > & vec);
	int read_row_bin   (base_vector< bigint > & vec, const bigint & p);

public:

	int reset ();
	int read_row ();
	int read_row (const bigint & p);

	void build_poly_in_X (Fp_polynomial & p, const bigmod & y);
	void build_poly_in_Y (Fp_polynomial & p, const bigmod & x);
};



#ifdef LIDIA_NAMESPACE
}	// end of namespace LiDIA
# undef IN_NAMESPACE_LIDIA
#endif



#endif	// LIDIA_MEQ_PRIME_H_GUARD_
