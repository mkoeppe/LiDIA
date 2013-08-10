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
//	Author	: Werner Backes (WB), Thorsten Lauer (TL)
//	Changes	: See CVS log
//
//==============================================================================================


#ifndef LIDIA_BIGINT_LATTICE_H_GUARD_
#define LIDIA_BIGINT_LATTICE_H_GUARD_


#ifndef LIDIA_BIGINT_H_GUARD_
# include	"LiDIA/bigint.h"
#endif
#ifndef LIDIA_BIGRATIONAL_H_GUARD_
# include	"LiDIA/bigrational.h"
#endif
#ifndef LIDIA_BIGFLOAT_H_GUARD_
# include	"LiDIA/bigfloat.h"
#endif
#ifndef LIDIA_BASE_VECTOR_H_GUARD_
# include	"LiDIA/base_vector.h"
#endif
#ifndef LIDIA_BIGINT_MATRIX_H_GUARD_
# include	"LiDIA/bigint_matrix.h"
#endif
#ifndef LIDIA_LATTICE_DEFS_H_GUARD_
# include	"LiDIA/lattices/lattice_defs.h"
#endif
#ifndef LIDIA_LATTICE_KERNEL_H_GUARD_
# include	"LiDIA/lattices/lattice_kernel.h"
#endif
#ifndef LIDIA_LATTICE_MODULES_H_GUARD_
# include	"LiDIA/lattices/lattice_modules.h"
#endif



#ifdef LIDIA_NAMESPACE
namespace LiDIA {
# define IN_NAMESPACE_LIDIA
#endif



//
// class bigint_lattice
// ~~~~~~~~~~~~~~~~~~~~
//
// A class derived form the class bigint_matrix, because we use
// a matrix for the representation of the lattice. This class
// allows you to work with lattices of bigint entries.
//
// You try to reduce the bigfloat problem to a kind of bigint problem,
// then use the bigint algorithm, and go back to bigfloats
// thats the reason for the friend declaration
//
class bigint_lattice:public bigint_matrix
{
private :
	friend class bigfloat_lattice;

public :

	//
	// Constructors / Destructor
	//
	bigint_lattice():bigint_matrix()
	{ bitfield.lattice_mode = DEFAULT_LATTICE_MODE; }
//    bigint_lattice(lidia_size_t n):bigint_matrix(n)
//      { bitfield.lattice_mode = DEFAULT_LATTICE_MODE; }
	bigint_lattice(lidia_size_t n, lidia_size_t m):bigint_matrix(n, m)
	{ bitfield.lattice_mode = DEFAULT_LATTICE_MODE; }
	bigint_lattice(lidia_size_t n, lidia_size_t m, const bigint** PPbin):
		bigint_matrix(n, m, PPbin)
	{ bitfield.lattice_mode = DEFAULT_LATTICE_MODE; }
	bigint_lattice(const bigint_matrix& b):bigint_matrix(b)
	{ bitfield.lattice_mode = DEFAULT_LATTICE_MODE; }

	~bigint_lattice() {};

	//
	// Special experimantal functions
	//
	//  randomize = permutation
	//  sort by length
	//
	void randomize_vectors();
	void sort_big_vectors();
	void sort_small_vectors();

	//
	// int bin_cmp_function(bigint *a, bigint *b, long)
	// returns
	// > 0     if    a > b
	// = 0     if    a = b
	// < 0     if    a < b
	//
	void sort_vectors(bin_cmp_func);

	//
	// Set special lattice characteristics
	//
	void set_red_orientation_columns();
	void set_red_orientation_rows();
	void set_gram_flag();
	void set_basis_flag();

	void delete_gram_flag();
	void delete_basis_flag();

	//
	// get the special lattice characteristics
	//
	bool get_red_orientation();
	bool get_basis_flag();
	bool get_gram_flag();

	//
	// check for special lattice characteristics
	//
	bool check_basis();
	bool check_gram();

	//
	// Test for lll - reduction
	//
	bool lll_check(sdigit, sdigit);
	bool lll_check(double y)
	{
		if (fabs(y) > 2.0) // allow correct checkings even if y*1000000
			y = 2.0; // forces overflow
		return(lll_check(static_cast<sdigit>(y*1000000), 1000000));
	}

	void lll_check_search(sdigit&, sdigit&);
	double lll_check_search()
	{
		sdigit nom, denom;
		lll_check_search(nom, denom);
		return(static_cast<double>(nom)/static_cast<double>(denom));
	}

	//
	// Interface to Algorithms
	//
	// Buchmann - Kessler   (only gensys/basis, no gram matrices)
	//
	// member functions
	void buchmann_kessler(math_matrix< bigint > &, double, lattice_info&);
	void buchmann_kessler(math_matrix< bigint > & T, double y)
	{
		lattice_info tli;
		buchmann_kessler(T, y, tli);
	}

	//
	// friend functions
	//
	friend bigint_lattice buchmann_kessler(const bigint_lattice&,
					       math_matrix< bigint > &,
					       lattice_info&);
	friend bigint_lattice buchmann_kessler(const bigint_lattice&,
					       math_matrix< bigint > &,
					       double);

	//
	// Schnorr - Euchner
	//
	// member functions
	//
	void lll_schnorr_euchner_orig(double, lattice_info&, sdigit = 1);
	void lll_schnorr_euchner_orig(math_matrix< bigint > &, double,
				      lattice_info&, sdigit = 1);

	void lll_schnorr_euchner_orig(double y, sdigit factor = 1)
	{
		lattice_info tli;
		lll_schnorr_euchner_orig(y, tli, factor);
	}
	void lll_schnorr_euchner_orig(math_matrix< bigint > & T, double y,
				      sdigit factor = 1)
	{
		lattice_info tli;
		lll_schnorr_euchner_orig(T, y, tli, factor);
	}
	//
	// user defined Scalar Product
	//
	void lll_schnorr_euchner_orig(double, lattice_info&, user_SP, sdigit = 1);
	void lll_schnorr_euchner_orig(math_matrix< bigint > &, double,
				      lattice_info&, user_SP, sdigit = 1);

	void lll_schnorr_euchner_orig(double y, user_SP, sdigit factor = 1)
	{
		lattice_info tli;
		lll_schnorr_euchner_orig(y, tli, factor);
	}
	void lll_schnorr_euchner_orig(math_matrix< bigint > & T, double y,
				      user_SP, sdigit factor = 1)
	{
		lattice_info tli;
		lll_schnorr_euchner_orig(T, y, tli, factor);
	}
	//
	// friend functions
	//
	friend bigint_lattice lll_schnorr_euchner_orig(const bigint_lattice&, double,
						       lattice_info&, sdigit);
	friend bigint_lattice lll_schnorr_euchner_orig(const bigint_lattice&,
						       math_matrix< bigint > &, double,
						       lattice_info&, sdigit);

	friend bigint_lattice lll_schnorr_euchner_orig(const bigint_lattice&,
						       double, sdigit);
	friend bigint_lattice lll_schnorr_euchner_orig(const bigint_lattice&,
						       math_matrix< bigint > &,
						       double, sdigit);
	//
	// user defined Scalar Product
	//
	friend bigint_lattice lll_schnorr_euchner_orig(const bigint_lattice&, double,
						       lattice_info&, user_SP, sdigit);
	friend bigint_lattice lll_schnorr_euchner_orig(const bigint_lattice&,
						       math_matrix< bigint > &, double,
						       lattice_info&, user_SP, sdigit);

	friend bigint_lattice lll_schnorr_euchner_orig(const bigint_lattice&, double,
						       user_SP, sdigit);
	friend bigint_lattice lll_schnorr_euchner_orig(const bigint_lattice&,
						       math_matrix< bigint > &, double,
						       user_SP, sdigit);

	//
	// Variation of Schnorr - Euchner
	//
	// member functions
	//
	void lll_schnorr_euchner_fact(double, lattice_info&, sdigit = 1);
	void lll_schnorr_euchner_fact(math_matrix< bigint > &, double,
				      lattice_info&, sdigit = 1);

	void lll_schnorr_euchner_fact(double y, sdigit factor = 1)
	{
		lattice_info tli;
		lll_schnorr_euchner_fact(y, tli, factor);
	}
	void lll_schnorr_euchner_fact(math_matrix< bigint > & T, double y, sdigit factor = 1)
	{
		lattice_info tli;
		lll_schnorr_euchner_fact(T, y, tli, factor);
	}
	//
	// user defined Scalar Product
	//
	void lll_schnorr_euchner_fact(double, lattice_info&, user_SP, sdigit = 1);
	void lll_schnorr_euchner_fact(math_matrix< bigint > &, double,
				      lattice_info&, user_SP, sdigit = 1);

	void lll_schnorr_euchner_fact(double y, user_SP, sdigit factor = 1)
	{
		lattice_info tli;
		lll_schnorr_euchner_fact(y, tli, factor);
	}
	void lll_schnorr_euchner_fact(math_matrix< bigint > & T, double y,
				      user_SP, sdigit factor = 1)
	{
		lattice_info tli;
		lll_schnorr_euchner_fact(T, y, tli, factor);
	}
	//
	// friend functions
	//
	friend bigint_lattice lll_schnorr_euchner_fact(const bigint_lattice&, double,
						       lattice_info&, sdigit);
	friend bigint_lattice lll_schnorr_euchner_fact(const bigint_lattice&,
						       math_matrix< bigint > &, double,
						       lattice_info&, sdigit);

	friend bigint_lattice lll_schnorr_euchner_fact(const bigint_lattice&, double,
						       sdigit);
	friend bigint_lattice lll_schnorr_euchner_fact(const bigint_lattice&,
						       math_matrix< bigint > &, double,
						       sdigit);
	//
	// user defined Scalar Product
	//
	friend bigint_lattice lll_schnorr_euchner_fact(const bigint_lattice&, double,
						       lattice_info&, user_SP, sdigit);
	friend bigint_lattice lll_schnorr_euchner_fact(const bigint_lattice&,
						       math_matrix< bigint > &, double,
						       lattice_info&, user_SP, sdigit);

	friend bigint_lattice lll_schnorr_euchner_fact(const bigint_lattice&, double,
						       user_SP, sdigit);
	friend bigint_lattice lll_schnorr_euchner_fact(const bigint_lattice&,
						       math_matrix< bigint > &, double,
						       user_SP, sdigit);

	//
	// Schnorr - Euchner (interface to fastest algorithm)
	//
	// member functions
	//
	void lll(double y, lattice_info& li, sdigit factor = 1)
	{
		lll_schnorr_euchner_orig(y, li, factor);
	}
	void lll(math_matrix< bigint > & T, double y,
		 lattice_info& li, sdigit factor = 1)
	{
		lll_schnorr_euchner_orig(T, y, li, factor);
	}

	void lll(double y, sdigit factor = 1)
	{
		lll_schnorr_euchner_orig(y, factor);
	}
	void lll(math_matrix< bigint > & T, double y, sdigit factor = 1)
	{
		lll_schnorr_euchner_orig(T, y, factor);
	}
	//
	// user defined Scalar Product
	//
	void lll(double y, lattice_info& li, user_SP SP, sdigit factor = 1)
	{
		lll_schnorr_euchner_orig(y, li, SP, factor);
	}
	void lll(math_matrix< bigint > & T, double y,
		 lattice_info& li, user_SP SP, sdigit factor = 1)
	{
		lll_schnorr_euchner_orig(T, y, li, SP, factor);
	}

	void lll(double y, user_SP SP, sdigit factor = 1)
	{
		lll_schnorr_euchner_orig(y, SP, factor);
	}
	void lll(math_matrix< bigint > & T, double y, user_SP SP, sdigit factor = 1)
	{
		lll_schnorr_euchner_orig(T, y, SP, factor);
	}
	//
	// friend functions
	//
	friend bigint_lattice lll(const bigint_lattice& A, double y,
				  lattice_info& li, sdigit factor);
	friend bigint_lattice lll(const bigint_lattice& A,
				  math_matrix< bigint > & T, double y,
				  lattice_info& li, sdigit factor);

	friend bigint_lattice lll(const bigint_lattice& A,
				  double y, sdigit factor);
	friend bigint_lattice lll(const bigint_lattice& A,
				  math_matrix< bigint > & T, double y,
				  sdigit factor);
	//
	// user defined Scalar Product
	//
	friend bigint_lattice lll(const bigint_lattice& A, double y,
				  lattice_info& li, user_SP SP,
				  sdigit factor);
	friend bigint_lattice lll(const bigint_lattice& A,
				  math_matrix< bigint > & T, double y,
				  lattice_info& li, user_SP SP,
				  sdigit factor);
	friend bigint_lattice lll(const bigint_lattice& A,
				  double y, user_SP SP, sdigit factor);
	friend bigint_lattice lll(const bigint_lattice& A,
				  math_matrix< bigint > & T, double y,
				  user_SP SP, sdigit factor);

	//
	// Benne de Weger
	//
	// member functions
	//
	void lll_benne_de_weger(sdigit, sdigit, lattice_info&);
	void lll_benne_de_weger(math_matrix< bigint > &, sdigit, sdigit,
				lattice_info&);

	void lll_benne_de_weger(sdigit y_nom, sdigit y_denom)
	{
		lattice_info tli;
		lll_benne_de_weger(y_nom, y_denom, tli);
	}
	void lll_benne_de_weger(math_matrix< bigint > &T, sdigit y_nom,
				sdigit y_denom)
	{
		lattice_info tli;
		lll_benne_de_weger(T, y_nom, y_denom, tli);
	}
	//
	// friend functions
	//
	friend bigint_lattice lll_benne_de_weger(const bigint_lattice&,
						 sdigit, sdigit, lattice_info&);
	friend bigint_lattice lll_benne_de_weger(const bigint_lattice&,
						 math_matrix< bigint > &,
						 sdigit, sdigit, lattice_info&);

	friend bigint_lattice lll_benne_de_weger(const bigint_lattice&,
						 sdigit, sdigit);
	friend bigint_lattice lll_benne_de_weger(const bigint_lattice&,
						 math_matrix< bigint > &,
						 sdigit, sdigit);

	//
	// Modified lll
	//
	// member functions
	//
	void mlll(double, bigint*&, lattice_info&);
	void mlll(double, base_vector< bigint > &, lattice_info&);

	void mlll(double y, bigint*& v)
	{
		lattice_info tli;
		mlll(y, v, tli);
	}
	void mlll(double y, base_vector< bigint > & v)
	{
		lattice_info tli;
		mlll(y, v, tli);
	}
	//
	// friend functions
	//
	friend bigint_lattice mlll(const bigint_lattice&, double,
				   bigint*&, lattice_info&);
	friend bigint_lattice mlll(const bigint_lattice&, double,
				   base_vector< bigint > &, lattice_info&);

	friend bigint_lattice mlll(const bigint_lattice&, double, bigint*&);
	friend bigint_lattice mlll(const bigint_lattice&, double,
				   base_vector< bigint > &);

	//
	// Shortest vector
	//
	lidia_size_t shortest_vector(math_matrix< bigint > &, lidia_size_t = 20);
	friend lidia_size_t shortest_vector(const bigint_lattice&,
					    math_matrix< bigint > &, lidia_size_t);

	//
	// Short vectors
	//
	void short_vector_coeff(math_matrix< bigint > &,
				const bigrational &, lidia_size_t = 20);
	friend void short_vector_coeff(const bigint_lattice&,
				       math_matrix< bigint > &,
				       const bigrational &, lidia_size_t);
	//
	// Close vector
	//
	void close_vector(const base_vector< bigint > &, base_vector< bigint > &,
			  sdigit = 1);
	friend base_vector< bigint > close_vector(const bigint_lattice& A,
						  const base_vector< bigint > & v,
						  sdigit x_factor);

protected :
	//
	// Create needed structure
	// Restore back in matrix structure (maybe different)
	//
	void Tr_dense_create(dense_alg< bigint > &);
	void Tr_dense_extract(const dense_alg< bigint > &);

	//
	// some flag checkings (...)
	//
	bool chk_basis();
	bool chk_gram();
	bool chk_trans();
	bool chk_reduce_columns();

	//
	// Dimension checking
	//
	bool chk_mlll_dim();
	bool chk_lll_dim();

	//
	// Other checkings
	//
	void chk_corr_param(double& y);
	void chk_corr_param(sdigit& nom, sdigit& denom);


	//
	// Helpful functions
	//
	sdigit TrD_search_factor(const dense_alg< bigint > &);
	bool TrD_gram_schmidt_orth_bdw(const dense_alg< bigint > &, bigint**,
				       bigint**, bigint*);
	bool TrD_lll_check(const dense_alg< bigint > &, sdigit, sdigit);
	void TrD_lll_check_search(const dense_alg< bigint > &, sdigit&, sdigit&);
	bool TrD_check_basis(const dense_alg< bigint > &);
	bool TrD_check_gram(const dense_alg< bigint > &);

	//
	//
	// Algorithms, real Implementaion
	//
	void TrD_randomize_vectors(dense_alg< bigint > &);
	void TrD_sort_vectors(dense_alg< bigint > &, sdigit);
	void TrD_sort_vectors(dense_alg< bigint > &, bin_cmp_func);

	//
	// Benne de Weger
	//
	void TrD_lll(dense_alg< bigint > &, lattice_info&);
	void TrD_lll_trans(dense_alg< bigint > &, lattice_info&);
	void TrD_lll_rand(dense_alg< bigint > &, lattice_info&);
	void TrD_lll_gram(dense_alg< bigint > &, lattice_info&);

	//
	// Schnorr - Euchner
	// the implementation is done in the class schnorr_euchner_bigint
	//

	//
	// Modified lll
	//
	bigint* TrD_mlll_bfl(dense_alg< bigint > &, lattice_info&);

	//
	// Buchmann - Kessler
	//
	void TrD_buchmann_kessler(dense_alg< bigint > &, lattice_info&);

	//
	// Useful functions needed by Buchmann - Kessler and mlll
	// for precision problems and needed constants
	//
	sdigit compute_read_precision(const dense_alg< bigint > &);
	sdigit compute_precision(const dense_alg< bigint > &);

	void alpha_compute(const dense_alg< bigint > &, bigfloat&);
	void gamma_compute(bigfloat&, sdigit);
	void zwei_pot_q_compute(const dense_alg< bigint > &, bigfloat&, sdigit&,
				bigfloat&);
};



// friend functions

bigint_lattice buchmann_kessler(const bigint_lattice&,
				math_matrix< bigint > &,
				lattice_info&);
bigint_lattice buchmann_kessler(const bigint_lattice&,
				math_matrix< bigint > &,
				double);

//
// Schnorr - Euchner
//

bigint_lattice lll_schnorr_euchner_orig(const bigint_lattice&, double,
					lattice_info&, sdigit = 1);
bigint_lattice lll_schnorr_euchner_orig(const bigint_lattice&,
					math_matrix< bigint > &, double,
					lattice_info&, sdigit = 1);

bigint_lattice lll_schnorr_euchner_orig(const bigint_lattice&,
					double, sdigit = 1);
bigint_lattice lll_schnorr_euchner_orig(const bigint_lattice&,
					math_matrix< bigint > &,
					double, sdigit = 1);
//
// user defined Scalar Product
//
bigint_lattice lll_schnorr_euchner_orig(const bigint_lattice&, double,
					lattice_info&, user_SP, sdigit = 1);
bigint_lattice lll_schnorr_euchner_orig(const bigint_lattice&,
					math_matrix< bigint > &, double,
					lattice_info&, user_SP, sdigit = 1);

bigint_lattice lll_schnorr_euchner_orig(const bigint_lattice&, double,
					user_SP, sdigit = 1);
bigint_lattice lll_schnorr_euchner_orig(const bigint_lattice&,
					math_matrix< bigint > &, double,
					user_SP, sdigit = 1);

//
// Variation of Schnorr - Euchner
//

bigint_lattice lll_schnorr_euchner_fact(const bigint_lattice&, double,
					lattice_info&, sdigit = 1);
bigint_lattice lll_schnorr_euchner_fact(const bigint_lattice&,
					math_matrix< bigint > &, double,
					lattice_info&, sdigit = 1);

bigint_lattice lll_schnorr_euchner_fact(const bigint_lattice&, double,
					sdigit = 1);
bigint_lattice lll_schnorr_euchner_fact(const bigint_lattice&,
					math_matrix< bigint > &, double,
					sdigit = 1);
//
// user defined Scalar Product
//
bigint_lattice lll_schnorr_euchner_fact(const bigint_lattice&, double,
					lattice_info&, user_SP, sdigit = 1);
bigint_lattice lll_schnorr_euchner_fact(const bigint_lattice&,
					math_matrix< bigint > &, double,
					lattice_info&, user_SP, sdigit = 1);

bigint_lattice lll_schnorr_euchner_fact(const bigint_lattice&, double,
					user_SP, sdigit = 1);
bigint_lattice lll_schnorr_euchner_fact(const bigint_lattice&,
					math_matrix< bigint > &, double,
					user_SP, sdigit = 1);

//
// Schnorr - Euchner (interface to fastest algorithm)
//

inline
bigint_lattice lll(const bigint_lattice& A, double y,
		   lattice_info& li, sdigit factor = 1) {
	return (LiDIA::lll_schnorr_euchner_orig(A, y, li, factor));
}

inline
bigint_lattice lll(const bigint_lattice& A,
		   math_matrix< bigint > & T, double y,
		   lattice_info& li, sdigit factor = 1) {
	return (LiDIA::lll_schnorr_euchner_orig(A, T, y, li, factor));
}

inline
bigint_lattice lll(const bigint_lattice& A,
		   double y, sdigit factor = 1) {
	return (LiDIA::lll_schnorr_euchner_orig(A, y, factor));
}

inline
bigint_lattice lll(const bigint_lattice& A,
		   math_matrix< bigint > & T, double y,
		   sdigit factor = 1) {
	return (LiDIA::lll_schnorr_euchner_orig(A, T, y, factor));
}

//
// user defined Scalar Product
//

inline
bigint_lattice lll(const bigint_lattice& A, double y,
		   lattice_info& li, user_SP SP,
		   sdigit factor = 1) {
	return (LiDIA::lll_schnorr_euchner_orig(A, y, li, SP, factor));
}

inline
bigint_lattice lll(const bigint_lattice& A,
		   math_matrix< bigint > & T, double y,
		   lattice_info& li, user_SP SP,
		   sdigit factor = 1) {
	return (LiDIA::lll_schnorr_euchner_orig(A, T, y, li, SP, factor));
}

inline
bigint_lattice lll(const bigint_lattice& A,
		   double y, user_SP SP, sdigit factor = 1) {
	return (LiDIA::lll_schnorr_euchner_orig(A, y, SP, factor));
}

inline
bigint_lattice lll(const bigint_lattice& A,
		   math_matrix< bigint > & T, double y,
		   user_SP SP, sdigit factor = 1) {
	return (LiDIA::lll_schnorr_euchner_orig(A, T, y, SP, factor));
}

//
// Benne de Weger
//

bigint_lattice lll_benne_de_weger(const bigint_lattice&,
				  sdigit, sdigit, lattice_info&);
bigint_lattice lll_benne_de_weger(const bigint_lattice&,
				  math_matrix< bigint > &,
				  sdigit, sdigit, lattice_info&);

bigint_lattice lll_benne_de_weger(const bigint_lattice&,
				  sdigit, sdigit);
bigint_lattice lll_benne_de_weger(const bigint_lattice&,
				  math_matrix< bigint > &,
				  sdigit, sdigit);

//
// Modified lll
//

bigint_lattice mlll(const bigint_lattice&, double,
		    bigint*&, lattice_info&);
bigint_lattice mlll(const bigint_lattice&, double,
		    base_vector< bigint > &, lattice_info&);

bigint_lattice mlll(const bigint_lattice&, double, bigint*&);
bigint_lattice mlll(const bigint_lattice&, double,
		    base_vector< bigint > &);

//
// Shortest vector
//

lidia_size_t shortest_vector(const bigint_lattice&,
			     math_matrix< bigint > &, lidia_size_t = 20);

//
// Short vectors
//

void short_vector_coeff(const bigint_lattice&,
			math_matrix< bigint > &,
			const bigrational &, lidia_size_t = 20);
//
// Close vector
//

inline
base_vector< bigint > close_vector(const bigint_lattice& A,
				   const base_vector< bigint > & v,
				   sdigit x_factor = 1) {
	bigint_lattice tA(A);
	base_vector< bigint > cv;
	tA.close_vector(v, cv, x_factor);
	return (cv);
}



#ifdef LIDIA_NAMESPACE
}	// end of namespace LiDIA
# undef IN_NAMESPACE_LIDIA
#endif



#endif	// LIDIA_BIGINT_LATTICE_H_GUARD_
