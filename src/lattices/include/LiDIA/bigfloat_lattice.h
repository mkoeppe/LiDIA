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


#ifndef LIDIA_BIGFLOAT_LATTICE_H_GUARD_
#define LIDIA_BIGFLOAT_LATTICE_H_GUARD_


#ifndef LIDIA_BIGFLOAT_H_GUARD_
# include	"LiDIA/bigfloat.h"
#endif
#ifndef LIDIA_BIGFLOAT_MATRIX_H_GUARD_
# include	"LiDIA/bigfloat_matrix.h"
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
#ifndef LIDIA_BIGINT_LATTICE_H_GUARD_
# include	"LiDIA/bigint_lattice.h"
#endif



#ifdef LIDIA_NAMESPACE
namespace LiDIA {
# define IN_NAMESPACE_LIDIA
#endif



//
// class bigfloat_lattice
// ~~~~~~~~~~~~~~~~~~~~~~
//
// A class derived form the class bigfloat_matrix, because we use
// a matrix for the representation of the lattice. This class
// allows you to work with lattices of bigfloat entries.
//
class bigfloat_lattice : public field_matrix< bigfloat >
{
private :
	static bigint_lattice CallBiL;

public :

	//
	// Constructors / Destructor
	//
	bigfloat_lattice():field_matrix< bigfloat > ()
	{ bitfield.lattice_mode = DEFAULT_LATTICE_MODE; };
//    bigfloat_lattice(lidia_size_t n):field_matrix < bigfloat > (n)
//      { bitfield.lattice_mode = DEFAULT_LATTICE_MODE; };
	bigfloat_lattice(lidia_size_t n, lidia_size_t m):field_matrix< bigfloat > (n, m)
	{ bitfield.lattice_mode = DEFAULT_LATTICE_MODE; };
	bigfloat_lattice(lidia_size_t n, lidia_size_t m, const bigfloat** PPbin):
		field_matrix< bigfloat > (n, m, PPbin)
	{ bitfield.lattice_mode = DEFAULT_LATTICE_MODE; };
	bigfloat_lattice(const field_matrix< bigfloat > & b):field_matrix< bigfloat > (b)
	{ bitfield.lattice_mode = DEFAULT_LATTICE_MODE; };

	~bigfloat_lattice() {};

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
	// int bin_cmp_function(bigfloat *a, bigfloat *b, long)
	// returns
	// > 0     if    a > b
	// = 0     if    a = b
	// < 0     if    a < b
	//
	void sort_vectors(bfl_cmp_func);

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
	void buchmann_kessler(ring_matrix< bigint > &, double, lattice_info&);
	void buchmann_kessler(field_matrix< bigfloat > &, double, lattice_info&);
	void buchmann_kessler(ring_matrix< bigint > & T, double y)
	{
		lattice_info tli;
		buchmann_kessler(T, y, tli);
	}
	void buchmann_kessler(field_matrix< bigfloat > & T, double y)
	{
		lattice_info tli;
		buchmann_kessler(T, y, tli);
	}

	//
	// friend functions
	//
	friend bigfloat_lattice buchmann_kessler(const bigfloat_lattice&,
						 ring_matrix< bigint > &,
						 lattice_info&);
	friend bigfloat_lattice buchmann_kessler(const bigfloat_lattice&,
						 field_matrix< bigfloat > &,
						 lattice_info&);
	friend bigfloat_lattice buchmann_kessler(const bigfloat_lattice&,
						 ring_matrix< bigint > &,
						 double);
	friend bigfloat_lattice buchmann_kessler(const bigfloat_lattice&,
						 field_matrix< bigfloat > &,
						 double);

	//
	// Schnorr - Euchner
	//
	// member functions
	//
	void lll_schnorr_euchner_orig(double, lattice_info&, sdigit = 1);
	void lll_schnorr_euchner_orig(ring_matrix< bigint > &, double,
				      lattice_info&, sdigit = 1);
	void lll_schnorr_euchner_orig(field_matrix< bigfloat > &, double,
				      lattice_info&, sdigit = 1);

	void lll_schnorr_euchner_orig(double y, sdigit factor = 1)
	{
		lattice_info tli;
		lll_schnorr_euchner_orig(y, tli, factor);
	}
	void lll_schnorr_euchner_orig(ring_matrix< bigint > & T, double y,
				      sdigit factor = 1)
	{
		lattice_info tli;
		lll_schnorr_euchner_orig(T, y, tli, factor);
	}
	void lll_schnorr_euchner_orig(field_matrix< bigfloat > & T, double y,
				      sdigit factor = 1)
	{
		lattice_info tli;
		lll_schnorr_euchner_orig(T, y, tli, factor);
	}
	//
	// user defined Scalar Product
	//
	void lll_schnorr_euchner_orig(double, lattice_info&, user_SP, sdigit = 1);
	void lll_schnorr_euchner_orig(ring_matrix< bigint > &, double,
				      lattice_info&, user_SP, sdigit = 1);
	void lll_schnorr_euchner_orig(field_matrix< bigfloat > &, double,
				      lattice_info&, user_SP, sdigit = 1);

	void lll_schnorr_euchner_orig(double y, user_SP, sdigit factor = 1)
	{
		lattice_info tli;
		lll_schnorr_euchner_orig(y, tli, factor);
	}
	void lll_schnorr_euchner_orig(ring_matrix< bigint > & T, double y,
				      user_SP, sdigit factor = 1)
	{
		lattice_info tli;
		lll_schnorr_euchner_orig(T, y, tli, factor);
	}
	void lll_schnorr_euchner_orig(field_matrix< bigfloat > & T, double y,
				      user_SP, sdigit factor = 1)
	{
		lattice_info tli;
		lll_schnorr_euchner_orig(T, y, tli, factor);
	}
	//
	// friend functions
	//
	friend bigfloat_lattice lll_schnorr_euchner_orig(const bigfloat_lattice&, double,
							 lattice_info&, sdigit);
	friend bigfloat_lattice lll_schnorr_euchner_orig(const bigfloat_lattice&,
							 ring_matrix< bigint > &, double,
							 lattice_info&, sdigit);
	friend bigfloat_lattice lll_schnorr_euchner_orig(const bigfloat_lattice&,
							 field_matrix< bigfloat > &, double,
							 lattice_info&, sdigit);

	friend bigfloat_lattice lll_schnorr_euchner_orig(const bigfloat_lattice&,
							 double, sdigit);
	friend bigfloat_lattice lll_schnorr_euchner_orig(const bigfloat_lattice&,
							 ring_matrix< bigint > &, double,
							 sdigit);
	friend bigfloat_lattice lll_schnorr_euchner_orig(const bigfloat_lattice&,
							 field_matrix< bigfloat > &, double,
							 sdigit);
	//
	// user defined Scalar Product
	//
	friend bigfloat_lattice lll_schnorr_euchner_orig(const bigfloat_lattice&, double,
							 lattice_info&, user_SP,
							 sdigit);
	friend bigfloat_lattice lll_schnorr_euchner_orig(const bigfloat_lattice&,
							 ring_matrix< bigint > &, double,
							 lattice_info&, user_SP,
							 sdigit);
	friend bigfloat_lattice lll_schnorr_euchner_orig(const bigfloat_lattice&,
							 field_matrix< bigfloat > &, double,
							 lattice_info&, user_SP,
							 sdigit);

	friend bigfloat_lattice lll_schnorr_euchner_orig(const bigfloat_lattice&, double,
							 user_SP, sdigit);
	friend bigfloat_lattice lll_schnorr_euchner_orig(const bigfloat_lattice&,
							 ring_matrix< bigint > &, double,
							 user_SP, sdigit);
	friend bigfloat_lattice lll_schnorr_euchner_orig(const bigfloat_lattice&,
							 field_matrix< bigfloat > &, double,
							 user_SP, sdigit);

	//
	// Variation of Schnorr - Euchner
	//
	// member functions
	//
	void lll_schnorr_euchner_fact(double, lattice_info&, sdigit = 1);
	void lll_schnorr_euchner_fact(ring_matrix< bigint > &, double,
				      lattice_info&, sdigit = 1);
	void lll_schnorr_euchner_fact(field_matrix< bigfloat > &, double,
				      lattice_info&, sdigit = 1);

	void lll_schnorr_euchner_fact(double y, sdigit factor = 1)
	{
		lattice_info tli;
		lll_schnorr_euchner_fact(y, tli, factor);
	}
	void lll_schnorr_euchner_fact(ring_matrix< bigint > & T, double y,
				      sdigit factor = 1)
	{
		lattice_info tli;
		lll_schnorr_euchner_fact(T, y, tli, factor);
	}
	void lll_schnorr_euchner_fact(field_matrix< bigfloat > & T, double y,
				      sdigit factor = 1)
	{
		lattice_info tli;
		lll_schnorr_euchner_fact(T, y, tli, factor);
	}
	//
	// user defined Scalar Product
	//
	void lll_schnorr_euchner_fact(double, lattice_info&, user_SP, sdigit = 1);
	void lll_schnorr_euchner_fact(ring_matrix< bigint > &, double,
				      lattice_info&, user_SP, sdigit = 1);
	void lll_schnorr_euchner_fact(field_matrix< bigfloat > &, double,
				      lattice_info&, user_SP, sdigit = 1);

	void lll_schnorr_euchner_fact(double y, user_SP, sdigit factor = 1)
	{
		lattice_info tli;
		lll_schnorr_euchner_fact(y, tli, factor);
	}
	void lll_schnorr_euchner_fact(ring_matrix< bigint > & T, double y,
				      user_SP, sdigit factor = 1)
	{
		lattice_info tli;
		lll_schnorr_euchner_fact(T, y, tli, factor);
	}
	void lll_schnorr_euchner_fact(field_matrix< bigfloat > & T, double y,
				      user_SP, sdigit factor = 1)
	{
		lattice_info tli;
		lll_schnorr_euchner_fact(T, y, tli, factor);
	}
	//
	// friend functions
	//
	friend bigfloat_lattice lll_schnorr_euchner_fact(const bigfloat_lattice&, double,
							 lattice_info&, sdigit);
	friend bigfloat_lattice lll_schnorr_euchner_fact(const bigfloat_lattice&,
							 ring_matrix< bigint > &, double,
							 lattice_info&, sdigit);
	friend bigfloat_lattice lll_schnorr_euchner_fact(const bigfloat_lattice&,
							 field_matrix< bigfloat > &, double,
							 lattice_info&, sdigit);

	friend bigfloat_lattice lll_schnorr_euchner_fact(const bigfloat_lattice&,
							 double, sdigit);
	friend bigfloat_lattice lll_schnorr_euchner_fact(const bigfloat_lattice&,
							 ring_matrix< bigint > &, double,
							 sdigit);
	friend bigfloat_lattice lll_schnorr_euchner_fact(const bigfloat_lattice&,
							 field_matrix< bigfloat > &, double,
							 sdigit);
	//
	// user defined Scalar Product
	//
	friend bigfloat_lattice lll_schnorr_euchner_fact(const bigfloat_lattice&, double,
							 lattice_info&, user_SP,
							 sdigit);
	friend bigfloat_lattice lll_schnorr_euchner_fact(const bigfloat_lattice&,
							 ring_matrix< bigint > &, double,
							 lattice_info&, user_SP,
							 sdigit);
	friend bigfloat_lattice lll_schnorr_euchner_fact(const bigfloat_lattice&,
							 field_matrix< bigfloat > &, double,
							 lattice_info&, user_SP,
							 sdigit);

	friend bigfloat_lattice lll_schnorr_euchner_fact(const bigfloat_lattice&, double,
							 user_SP, sdigit);
	friend bigfloat_lattice lll_schnorr_euchner_fact(const bigfloat_lattice&,
							 ring_matrix< bigint > &, double,
							 user_SP, sdigit);
	friend bigfloat_lattice lll_schnorr_euchner_fact(const bigfloat_lattice&,
							 field_matrix< bigfloat > &, double,
							 user_SP, sdigit);

	//
	// Schnorr - Euchner (interface to fastest(stable) algorithm)
	//
	// member functions
	//
	void lll(double y, lattice_info& li, sdigit factor = 1)
	{ lll_schnorr_euchner_fact(y, li, factor); }
	void lll(ring_matrix< bigint > & T, double y,
		 lattice_info& li, sdigit factor = 1)
	{ lll_schnorr_euchner_fact(T, y, li, factor); }
	void lll(field_matrix< bigfloat > & T, double y,
		 lattice_info& li, sdigit factor = 1)
	{ lll_schnorr_euchner_fact(T, y, li, factor); }

	void lll(double y, sdigit factor = 1)
	{ lll_schnorr_euchner_fact(y, factor); }
	void lll(ring_matrix< bigint > & T, double y, sdigit factor = 1)
	{ lll_schnorr_euchner_fact(T, y, factor); }
	void lll(field_matrix< bigfloat > & T, double y, sdigit factor = 1)
	{ lll_schnorr_euchner_fact(T, y, factor); }
	//
	// user defined Scalar Product
	//
	void lll(double y, lattice_info& li, user_SP SP, sdigit factor = 1)
	{ lll_schnorr_euchner_fact(y, li, SP, factor); }
	void lll(ring_matrix< bigint > & T, double y,
		 lattice_info& li, user_SP SP, sdigit factor = 1)
	{ lll_schnorr_euchner_fact(T, y, li, SP, factor); }
	void lll(field_matrix< bigfloat > & T, double y,
		 lattice_info& li, user_SP SP, sdigit factor = 1)
	{ lll_schnorr_euchner_fact(T, y, li, SP, factor); }

	void lll(double y, user_SP SP, sdigit factor = 1)
	{ lll_schnorr_euchner_fact(y, SP, factor); }
	void lll(ring_matrix< bigint > & T, double y, user_SP SP, sdigit factor = 1)
	{ lll_schnorr_euchner_fact(T, y, SP, factor); }
	void lll(field_matrix< bigfloat > & T, double y, user_SP SP, sdigit factor = 1)
	{ lll_schnorr_euchner_fact(T, y, SP, factor); }
	//
	// friend functions
	//
	friend bigfloat_lattice lll(const bigfloat_lattice& A, double y,
				    lattice_info& li, sdigit factor);
	friend bigfloat_lattice lll(const bigfloat_lattice& A,
				    ring_matrix< bigint > & T, double y,
				    lattice_info& li, sdigit factor);
	friend bigfloat_lattice lll(const bigfloat_lattice& A,
				    field_matrix< bigfloat > & T, double y,
				    lattice_info& li, sdigit factor);

	friend bigfloat_lattice lll(const bigfloat_lattice& A,
				    double y, sdigit factor);
	friend bigfloat_lattice lll(const bigfloat_lattice& A,
				    ring_matrix< bigint > & T, double y,
				    sdigit factor);
	friend bigfloat_lattice lll(const bigfloat_lattice& A,
				    field_matrix< bigfloat > & T, double y,
				    sdigit factor);
	//
	// user defined Scalar Product
	//
	friend bigfloat_lattice lll(const bigfloat_lattice& A, double y,
				    lattice_info& li, user_SP SP,
				    sdigit factor);
	friend bigfloat_lattice lll(const bigfloat_lattice& A,
				    ring_matrix< bigint > & T, double y,
				    lattice_info& li, user_SP SP,
				    sdigit factor);
	friend bigfloat_lattice lll(const bigfloat_lattice& A,
				    field_matrix< bigfloat > & T, double y,
				    lattice_info& li, user_SP SP,
				    sdigit factor);

	friend bigfloat_lattice lll(const bigfloat_lattice& A,
				    double y, user_SP SP, sdigit factor);
	friend bigfloat_lattice lll(const bigfloat_lattice& A,
				    ring_matrix< bigint > & T, double y,
				    user_SP SP, sdigit factor);
	friend bigfloat_lattice lll(const bigfloat_lattice& A,
				    field_matrix< bigfloat > & T, double y,
				    user_SP SP, sdigit factor);

	//
	// Modified lll
	//
	// member functions
	//
	void mlll(double, bigint*&, lattice_info&);
	void mlll(double, bigfloat*&, lattice_info&);
	void mlll(double, base_vector< bigint > &, lattice_info&);
	void mlll(double, base_vector< bigfloat > &, lattice_info&);

	void mlll(double y, bigint*& v)
	{
		lattice_info tli;
		mlll(y, v, tli);
	}
	void mlll(double y, bigfloat*& v)
	{
		lattice_info tli;
		mlll(y, v, tli);
	}
	void mlll(double y, base_vector< bigint > & v)
	{
		lattice_info tli;
		mlll(y, v, tli);
	}
	void mlll(double y, base_vector< bigfloat > & v)
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
				   bigfloat*&, lattice_info&);
	friend bigint_lattice mlll(const bigint_lattice&, double,
				   base_vector< bigint > &, lattice_info&);
	friend bigint_lattice mlll(const bigint_lattice&, double,
				   base_vector< bigfloat > &, lattice_info&);

	friend bigint_lattice mlll(const bigint_lattice&, double, bigint*&);
	friend bigint_lattice mlll(const bigint_lattice&, double, bigfloat*&);
	friend bigint_lattice mlll(const bigint_lattice&, double,
				   base_vector< bigint > &);
	friend bigint_lattice mlll(const bigint_lattice&, double,
				   base_vector< bigfloat > &);

	//
	// Close vector
	//
	void close_vector(const base_vector< bigfloat > &, base_vector< bigfloat > &,
			  sdigit = 1);
	friend base_vector< bigfloat > close_vector(const bigfloat_lattice& A,
						    const base_vector< bigfloat > & v,
						    sdigit x_factor);


protected :
	//
	// Create needed structure
	// Restore back in matrix structure (maybe different)
	//
	void Tr_dense_create(dense_alg< bigfloat > &);
	void Tr_dense_extract(const dense_alg< bigfloat > &);

	//
	// Generate bigint_lattice / factor
	//
	void Tr_dense_create_f(dense_alg< bigint > &);
	void Tr_dense_extract_f(const dense_alg< bigint > &);

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
	// Helpfull functions
	//
	sdigit TrD_search_factor(const dense_alg< bigfloat > &);
	void conv_bfl_bin(field_matrix< bigfloat > &, ring_matrix< bigint > &);
	void conv_bin_bfl(ring_matrix< bigint > &, field_matrix< bigfloat > &);

	//
	// Algorithms, real Implementaion
	//
	void TrD_randomize_vectors(dense_alg< bigfloat > &);
	void TrD_sort_vectors(dense_alg< bigfloat > &, sdigit);
	void TrD_sort_vectors(dense_alg< bigfloat > &, bfl_cmp_func);

	//
	// Schnorr - Euchner
	// the implementation is done in the class schnorr_euchner_bigfloat
	//

	//
	// Modified lll
	//
	bigint* TrD_mlll_bfl(dense_alg< bigfloat > &, lattice_info&);

	//
	// Buchmann - Kessler
	//
	void TrD_buchmann_kessler(dense_alg< bigfloat > &, lattice_info&);

	//
	// Useful functions needed by Buchmann - Kessler and mlll
	// for precision problems and needed constants
	//
	sdigit compute_read_precision(const dense_alg< bigfloat > &);
	sdigit compute_precision(const dense_alg< bigfloat > &);

	void alpha_compute(const dense_alg< bigfloat > &, bigfloat&);
	void gamma_compute(bigfloat&, sdigit);
	void zwei_pot_q_compute(const dense_alg< bigfloat > &, bigfloat&, sdigit&,
				bigfloat&);
};


// friend functions

bigfloat_lattice buchmann_kessler(const bigfloat_lattice&,
				  ring_matrix< bigint > &,
				  lattice_info&);
bigfloat_lattice buchmann_kessler(const bigfloat_lattice&,
				  field_matrix< bigfloat > &,
				  lattice_info&);
bigfloat_lattice buchmann_kessler(const bigfloat_lattice&,
				  ring_matrix< bigint > &,
				  double);
bigfloat_lattice buchmann_kessler(const bigfloat_lattice&,
				  field_matrix< bigfloat > &,
				  double);

//
// Schnorr - Euchner
//

bigfloat_lattice lll_schnorr_euchner_orig(const bigfloat_lattice&, double,
					  lattice_info&, sdigit = 1);
bigfloat_lattice lll_schnorr_euchner_orig(const bigfloat_lattice&,
					  ring_matrix< bigint > &, double,
					  lattice_info&, sdigit = 1);
bigfloat_lattice lll_schnorr_euchner_orig(const bigfloat_lattice&,
					  field_matrix< bigfloat > &, double,
					  lattice_info&, sdigit = 1);

bigfloat_lattice lll_schnorr_euchner_orig(const bigfloat_lattice&,
					  double, sdigit = 1);
bigfloat_lattice lll_schnorr_euchner_orig(const bigfloat_lattice&,
					  ring_matrix< bigint > &, double,
					  sdigit = 1);
bigfloat_lattice lll_schnorr_euchner_orig(const bigfloat_lattice&,
					  field_matrix< bigfloat > &, double,
					  sdigit = 1);
//
// user defined Scalar Product
//
bigfloat_lattice lll_schnorr_euchner_orig(const bigfloat_lattice&, double,
					  lattice_info&, user_SP,
					  sdigit = 1);
bigfloat_lattice lll_schnorr_euchner_orig(const bigfloat_lattice&,
					  ring_matrix< bigint > &, double,
					  lattice_info&, user_SP,
					  sdigit = 1);
bigfloat_lattice lll_schnorr_euchner_orig(const bigfloat_lattice&,
					  field_matrix< bigfloat > &, double,
					  lattice_info&, user_SP,
					  sdigit = 1);

bigfloat_lattice lll_schnorr_euchner_orig(const bigfloat_lattice&, double,
					  user_SP, sdigit = 1);
bigfloat_lattice lll_schnorr_euchner_orig(const bigfloat_lattice&,
					  ring_matrix< bigint > &, double,
					  user_SP, sdigit = 1);
bigfloat_lattice lll_schnorr_euchner_orig(const bigfloat_lattice&,
					  field_matrix< bigfloat > &, double,
					  user_SP, sdigit = 1);

//
// Variation of Schnorr - Euchner
//

bigfloat_lattice lll_schnorr_euchner_fact(const bigfloat_lattice&, double,
					  lattice_info&, sdigit = 1);
bigfloat_lattice lll_schnorr_euchner_fact(const bigfloat_lattice&,
					  ring_matrix< bigint > &, double,
					  lattice_info&, sdigit = 1);
bigfloat_lattice lll_schnorr_euchner_fact(const bigfloat_lattice&,
					  field_matrix< bigfloat > &, double,
					  lattice_info&, sdigit = 1);
    
bigfloat_lattice lll_schnorr_euchner_fact(const bigfloat_lattice&,
					  double, sdigit = 1);
bigfloat_lattice lll_schnorr_euchner_fact(const bigfloat_lattice&,
					  ring_matrix< bigint > &, double,
					  sdigit = 1);
bigfloat_lattice lll_schnorr_euchner_fact(const bigfloat_lattice&,
					  field_matrix< bigfloat > &, double,
					  sdigit = 1);
//
// user defined Scalar Product
//
bigfloat_lattice lll_schnorr_euchner_fact(const bigfloat_lattice&, double,
					  lattice_info&, user_SP,
					  sdigit = 1);
bigfloat_lattice lll_schnorr_euchner_fact(const bigfloat_lattice&,
					  ring_matrix< bigint > &, double,
					  lattice_info&, user_SP,
					  sdigit = 1);
bigfloat_lattice lll_schnorr_euchner_fact(const bigfloat_lattice&,
					  field_matrix< bigfloat > &, double,
					  lattice_info&, user_SP,
					  sdigit = 1);

bigfloat_lattice lll_schnorr_euchner_fact(const bigfloat_lattice&, double,
					  user_SP, sdigit = 1);
bigfloat_lattice lll_schnorr_euchner_fact(const bigfloat_lattice&,
					  ring_matrix< bigint > &, double,
					  user_SP, sdigit = 1);
bigfloat_lattice lll_schnorr_euchner_fact(const bigfloat_lattice&,
					  field_matrix< bigfloat > &, double,
					  user_SP, sdigit = 1);

//
// Schnorr - Euchner (interface to fastest(stable) algorithm)
//

inline
bigfloat_lattice lll(const bigfloat_lattice& A, double y,
		     lattice_info& li, sdigit factor = 1) {
    return (LiDIA::lll_schnorr_euchner_fact(A, y, li, factor));
}

inline
bigfloat_lattice lll(const bigfloat_lattice& A,
		     ring_matrix< bigint > & T, double y,
		     lattice_info& li, sdigit factor = 1) {
    return (LiDIA::lll_schnorr_euchner_fact(A, T, y, li, factor));
}

inline
bigfloat_lattice lll(const bigfloat_lattice& A,
		     field_matrix< bigfloat > & T, double y,
		     lattice_info& li, sdigit factor = 1) {
    return (LiDIA::lll_schnorr_euchner_fact(A, T, y, li, factor));
}

inline
bigfloat_lattice lll(const bigfloat_lattice& A,
		     double y, sdigit factor = 1) {
    return (LiDIA::lll_schnorr_euchner_fact(A, y, factor));
}

inline
bigfloat_lattice lll(const bigfloat_lattice& A,
		     ring_matrix< bigint > & T, double y,
		     sdigit factor = 1) {
    return (LiDIA::lll_schnorr_euchner_fact(A, T, y, factor));
}

inline
bigfloat_lattice lll(const bigfloat_lattice& A,
		     field_matrix< bigfloat > & T, double y,
		     sdigit factor = 1) {
    return (LiDIA::lll_schnorr_euchner_fact(A, T, y, factor));
}

//
// user defined Scalar Product
//
    
inline
bigfloat_lattice lll(const bigfloat_lattice& A, double y,
		     lattice_info& li, user_SP SP,
		     sdigit factor = 1) {
    return (LiDIA::lll_schnorr_euchner_fact(A, y, li, SP, factor));
}

inline
bigfloat_lattice lll(const bigfloat_lattice& A,
		     ring_matrix< bigint > & T, double y,
		     lattice_info& li, user_SP SP,
		     sdigit factor = 1) {
    return (LiDIA::lll_schnorr_euchner_fact(A, T, y, li, SP, factor));
}

inline
bigfloat_lattice lll(const bigfloat_lattice& A,
		     field_matrix< bigfloat > & T, double y,
		     lattice_info& li, user_SP SP,
		     sdigit factor = 1) {
    return (LiDIA::lll_schnorr_euchner_fact(A, T, y, li, SP, factor));
}

inline
bigfloat_lattice lll(const bigfloat_lattice& A,
		     double y, user_SP SP, sdigit factor = 1) {
    return (LiDIA::lll_schnorr_euchner_fact(A, y, SP, factor));
}

inline
bigfloat_lattice lll(const bigfloat_lattice& A,
		     ring_matrix< bigint > & T, double y,
		     user_SP SP, sdigit factor = 1) {
    return (LiDIA::lll_schnorr_euchner_fact(A, T, y, SP, factor));
}

inline
bigfloat_lattice lll(const bigfloat_lattice& A,
		     field_matrix< bigfloat > & T, double y,
		     user_SP SP, sdigit factor = 1) {
    return (LiDIA::lll_schnorr_euchner_fact(A, T, y, SP, factor));
}

//
// Modified lll
//

bigint_lattice mlll(const bigint_lattice&, double,
		    bigint*&, lattice_info&);
bigint_lattice mlll(const bigint_lattice&, double,
		    bigfloat*&, lattice_info&);
bigint_lattice mlll(const bigint_lattice&, double,
		    base_vector< bigint > &, lattice_info&);
bigint_lattice mlll(const bigint_lattice&, double,
		    base_vector< bigfloat > &, lattice_info&);

bigint_lattice mlll(const bigint_lattice&, double, bigint*&);
bigint_lattice mlll(const bigint_lattice&, double, bigfloat*&);
bigint_lattice mlll(const bigint_lattice&, double,
		    base_vector< bigint > &);
bigint_lattice mlll(const bigint_lattice&, double,
		    base_vector< bigfloat > &);

//
// Close vector
//

inline
base_vector< bigfloat > close_vector(const bigfloat_lattice& A,
					 const base_vector< bigfloat > & v,
					 sdigit x_factor = 1)
{
	bigfloat_lattice tA(A);
	base_vector< bigfloat > cv;
	tA.close_vector(v, cv, x_factor);
	return (cv);
}



#ifdef LIDIA_NAMESPACE
}	// end of namespace LiDIA
# undef IN_NAMESPACE_LIDIA
#endif



#endif	// LIDIA_BIGFLOAT_LATTICE_H_GUARD_
