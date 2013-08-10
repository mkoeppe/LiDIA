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
//	$Id: crt.h,v 2.6 2004/06/15 10:19:16 lidiaadm Exp $
//
//	Author	: Frank J. Lehmann (FL), Thomas Pfahler (TPF)
//	Changes	: See CVS log
//
//==============================================================================================


#ifndef LIDIA_CRT_H_GUARD_
#define LIDIA_CRT_H_GUARD_



#ifndef LIDIA_BIGINT_H_GUARD_
# include	"LiDIA/bigint.h"
#endif
#ifndef LIDIA_BIGMOD_H_GUARD_
# include	"LiDIA/bigmod.h"
#endif
#ifndef LIDIA_BASE_VECTOR_H_GUARD_
# include	"LiDIA/base_vector.h"
#endif
#ifndef LIDIA_BASE_MATRIX_H_GUARD_
# include	"LiDIA/base_matrix.h"
#endif
#ifndef LIDIA_UDIGIT_H_GUARD_
# include	"LiDIA/udigit.h"
#endif
#ifndef LIDIA_CRT_TABLE_H_GUARD_
# include	"LiDIA/crt_table.h"
#endif



#ifdef LIDIA_NAMESPACE
namespace LiDIA {
# define IN_NAMESPACE_LIDIA
#endif



//***************************************************************
//
//					class crt
//
//***************************************************************


class crt
{
private :

	static const char INIT;
        static const char SINGLE;
        static const char VECTOR;
        static const char POINTER;
        static const char MATRIX;

        static const char MODULAR;
        static const char INTEGER;


	crt_table *m;

	lidia_size_t used; // number of used primes when combining
	unsigned char *ix; // ix[i] == 1 iff primes[i] has been used already

	char mode;

	bigint H;
	double d;

	lidia_size_t size;
	bigint* Hv;
	double* dv;

	base_vector< bigint > Hvec;
	base_vector< double > dvec;


	base_matrix< bigint > Hmat;
	base_matrix< double > dmat;



	void init_single ();
	void init_pointer(lidia_size_t sz);
	void init_vector (lidia_size_t sz);
	void init_matrix (lidia_size_t row, lidia_size_t col);

	void clear_pointer();
	void clear_vector ();
	void clear_matrix ();

	bool test_mode (char m);


public:

	crt ();
	crt (crt_table &M);
	~crt ();


	void init  (crt_table & M);
	void reset ();
	void clear ();

	sdigit get_prime (lidia_size_t index) const;


	void reduce (sdigit &small_value, const bigint &big, lidia_size_t index) const;
	void reduce (sdigit* small_vec, const bigint* big_vec, lidia_size_t bsize, lidia_size_t index) const;
	void reduce (base_vector< sdigit > & small_vec, const base_vector< bigint > & big_vec, lidia_size_t index) const;

	void reduce (base_matrix< sdigit > & small_mat, const base_matrix< bigint > & big_mat, lidia_size_t index) const;


	void reduce (sdigit &small_value, const bigmod &big, lidia_size_t index) const;
	void reduce (sdigit* small_vec, const bigmod* big_vec, lidia_size_t bsize, lidia_size_t index) const;
	void reduce (base_vector< sdigit > & small_vec, const base_vector< bigmod > & big_vec, lidia_size_t index) const;

	void reduce (base_matrix< sdigit > & small_mat, const base_matrix< bigmod > & big_mat, lidia_size_t index) const;


	void combine (const sdigit &  small_vec, lidia_size_t index);
	void combine (const sdigit* small_vec, lidia_size_t l, lidia_size_t index);
	void combine (const udigit* small_vec, lidia_size_t l, lidia_size_t index);
	void combine (const base_vector< sdigit > & small_vec, lidia_size_t index);

	void combine (const base_matrix< sdigit > & small_vec, lidia_size_t index);


	void get_result (bigint & big);
	void get_result (bigint*& big_vec, lidia_size_t & l);
	void get_result (base_vector< bigint > & big_vec);

	void get_result (base_matrix< bigint > & big_mat);


	void get_result (bigmod & big);
	void get_result (bigmod*& big_vec, lidia_size_t & l);
	void get_result (base_vector< bigmod > & big_vec);

	void get_result (base_matrix< bigmod > & big_mat);


	lidia_size_t number_of_primes() const;
	lidia_size_t how_much_to_use (const bigint &B) const;

	void info () const;

};
// class crt



#ifdef LIDIA_NAMESPACE
}	// end of namespace LiDIA
# undef IN_NAMESPACE_LIDIA
#endif



#endif	// LIDIA_CRT_H_GUARD_
