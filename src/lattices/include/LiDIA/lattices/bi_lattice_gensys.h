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
//	Changes	: see CVS log
//
//==============================================================================================


#ifndef LIDIA_BI_LATTICE_GENSYS_H_GUARD_
#define LIDIA_BI_LATTICE_GENSYS_H_GUARD_


#ifndef LIDIA_BIGFLOAT_MATRIX_H_GUARD_
# include	"LiDIA/bigfloat_matrix.h"
#endif
#ifndef LIDIA_BASE_VECTOR_H_GUARD_
# include	"LiDIA/base_vector.h"
#endif
#ifndef LIDIA_BIGINT_H_GUARD_
# include	"LiDIA/bigint.h"
#endif
#ifndef LIDIA_BIGFLOAT_H_GUARD_
# include	"LiDIA/bigfloat.h"
#endif
#ifndef LIDIA_TIMER_H_GUARD_
# include	"LiDIA/timer.h"
#endif
#ifndef LIDIA_RANDOM_GENERATOR_H_GUARD_
# include	"LiDIA/random_generator.h"
#endif
#include	<cstdlib>



#ifdef LIDIA_NAMESPACE
namespace LiDIA {
# define IN_NAMESPACE_LIDIA
#endif



#ifndef rint
#define rint(x) (((x)-std::floor(x)) < 0.5 ? std::floor(x) : std::ceil(x))
#endif



typedef int (*bin_cmp)(const bigint*, const bigint*, long);

//
// you try to reduce  the bigfloat problem to a kind of bigint problem,
// then use the bigint algorithm, and go back to bigfloats
// thats the reason for the friend declaration
//

class bigfloat_lattice_gensys;
class bigfloat_lattice_basis;



class bigint_lattice_gensys : public ring_matrix< bigint >
{
protected :
	friend          class bigfloat_lattice_gensys;
	friend          class bigfloat_lattice_basis;

//
// Global configuration
//
	sdigit          y_nom;
	sdigit          y_denom;
	double          y_par;

	bool            trans_flag;

//
// algorithm - information
//
	sdigit          reduction_steps;
	sdigit          swapping_steps;
	sdigit          correction_steps;

//
// static variables for temporary use to increase speed
//
	static bigint tempmz0;
	static bigint tempmz1;
	static bigint tempmz2;
	static bigint tempmz3;
	static bigint tempmz4;
	static bigint ergmz;

	static double   vectdblz;
	static bigint   vectbinz;
	static bigfloat vectbflz;

	static lidia_size_t vectsize;

public :

//
// Constructor / Destructors
//
	bigint_lattice_gensys();
//    bigint_lattice_gensys(lidia_size_t);
	bigint_lattice_gensys(lidia_size_t, lidia_size_t);
	bigint_lattice_gensys(lidia_size_t, lidia_size_t, const bigint**);
	bigint_lattice_gensys(const ring_matrix< bigint > &);
	bigint_lattice_gensys(const bigint_lattice_gensys&);
	~bigint_lattice_gensys();

//
//  Assignments
//
	void assign(const bigint_lattice_gensys&);
	bigint_lattice_gensys& operator = (const bigint_lattice_gensys&);
	friend void assign(bigint_lattice_gensys&, const bigint_lattice_gensys&);

//
// Set or get reduction info, paramters
//
	void set_reduction_parameter(double);
	void set_reduction_parameter(sdigit, sdigit);
	double get_reduction_parameter();
	void get_reduction_parameter(sdigit&, sdigit&);

//
// Set representation direction, reduce
// column by column or row by row
//
	void set_representation_columns();
	void set_representation_rows();
	bool get_representation()
	{
		return (trans_flag);
	}

//
// Get detailed information about lll - algorithm
//
	sdigit get_no_of_reduction_steps()
	{
		return (reduction_steps);
	}
	sdigit get_no_of_swaps()
	{
		return (swapping_steps);
	}
	sdigit get_no_of_corrections()
	{
		return (correction_steps);
	}

//
// Algorithms - Interface
//

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
// int bfl_cmp_function(bigint *a, bigint *b, long)
// returns
// > 0     if    a > b
// = 0     if    a = b
// < 0     if    a < b
//
	void sort_vectors(bin_cmp);

//
// Buchmann - Kessler Algorthm for
// linaer generating systems
//
	void lin_gen_system(ring_matrix< bigint > &, lidia_size_t&);

	friend bigint_lattice_gensys lin_gen_system(const bigint_lattice_gensys& L,
						    lidia_size_t& rank);

//
// Variation of the Schnorr - Euchner - lllfp using doubles for
// paramter = 1 and x*doubleprecision for parameter = x
//
	void lll_gensys(lidia_size_t&, sdigit = 1);
	void lll_gensys(const ring_matrix< bigint > &, lidia_size_t&, sdigit = 1);
	friend bigint_lattice_gensys lll_gensys(const bigint_lattice_gensys& L,
						lidia_size_t& rank,
						sdigit factor);


	void lll_trans_gensys(ring_matrix< bigint > &, lidia_size_t&, sdigit = 1);
	friend bigint_lattice_gensys lll_trans_gensys(const bigint_lattice_gensys& L,
						      ring_matrix< bigint > & T,
						      lidia_size_t& rank,
						      sdigit factor);


//
// Modified lll for genertating systems
// of the form n x n+1
//
// Returns vector of dimension n
// where you can find the relations
//
	void mlll(base_vector< bigint > &);
	void mlll(const ring_matrix< bigint > &, base_vector< bigint > &);
	friend bigint_lattice_gensys mlll(const bigint_lattice_gensys& L,
					  base_vector< bigint > & v);

	void mlll(bigint *&);
	void mlll(const ring_matrix< bigint > &, bigint *&);
	friend bigint_lattice_gensys mlll(const bigint_lattice_gensys& L,
					  bigint*& v);

//
// Compute a basis from lattice A and Transformation lattice T
//
	void compute_basis(const ring_matrix< bigint > & A, const ring_matrix< bigint > & T)
	{
		debug_handler("bigint_lattice_gensys", "compute_basis(A, T)");
		if (trans_flag)
			LiDIA::multiply(*this, A, T);
		else
			LiDIA::multiply(*this, T, A);
	}

protected :

//
// Functions needed by every constructor
//
	void init_parameters();
	void do_some_other_stuff() {};
	void assign_the_rest(const bigint_lattice_gensys&);

	sdigit compute_read_precision();
	sdigit compute_precision();

	void alpha_compute(bigfloat&);
	void gamma_compute(bigfloat&, sdigit);
	void zwei_pot_q_compute(bigfloat&, sdigit&, bigfloat&);

//
// very funny, but helpfull, speedup
//
	void Tr_trans_swap(bigint_lattice_gensys&);
	void bigintify(sdigit&, const bigfloat_matrix &);

//
// Dimension checking
//
	bool Tr_check_mlll()
	{
		if (trans_flag == true) {
			if (rows+1 != columns)
				return(false);
			else
				return(true);
		}
		else {
			if (columns+1 != rows)
				return(false);
			else
				return(true);
		}
	}


//
// Algorithms, real Implementaion
//
	void Tr_randomize_vectors();
	void Tr_sort_vectors(sdigit);
	void Tr_sort_vectors(bin_cmp);

	void Tr_lll();
	void Tr_lll_dbl();
	void Tr_lll_bfl(sdigit);

	void Tr_lll_deep_insert_dbl(sdigit); // not working yet
	void Tr_lll_deep_insert_bfl(sdigit, sdigit); // not yet implemented (s.o.)

	void Tr_lll_trans(ring_matrix< bigint > &);
	void Tr_lll_trans_dbl(ring_matrix< bigint > &);
	void Tr_lll_trans_bfl(ring_matrix< bigint > &, sdigit);

	void Tr_lll_rand();

	void Tr_mlll_dbl(bigint*); // not yet implemented
	void Tr_mlll_bfl(bigint*);

	void Tr_lin_gen_system(ring_matrix< bigint > &, lidia_size_t&);

	void Tr_lll_dbl_gensys(lidia_size_t&);
	void Tr_lll_bfl_gensys(lidia_size_t&, sdigit = 2);

	void Tr_lll_trans_dbl_gensys(ring_matrix< bigint > &, lidia_size_t&);
	void Tr_lll_trans_bfl_gensys(ring_matrix< bigint > &, lidia_size_t&, sdigit = 2);


//
// "Vector" - Operations
//

//
// For double vectors
//
	friend void bin_assign_dbl(double*, double*);
	friend void bin_assign_zero_dbl(double*);
	friend void bin_add_dbl(double*, double*, double*);
	friend void bin_subtract_dbl(double*, double*, double*);

	friend void bin_scalmul_dbl(double*, const double&, double*);
	friend void bin_scalprod_dbl(double&, double*, double*);
	friend void bin_scalquad_dbl(double&, double*);

	friend void bin_l1_norm_dbl(double&, double*);
	friend void bin_l2_norm_dbl(double&, double*);

//
// for bigfloat vectors
//
	friend void bin_assign_bfl(bigfloat*, bigfloat*);
	friend void bin_assign_zero_bfl(bigfloat*);
	friend void bin_add_bfl(bigfloat*, bigfloat*, bigfloat*);
	friend void bin_subtract_bfl(bigfloat*, bigfloat*, bigfloat*);

	friend void bin_scalmul_bfl(bigfloat*, const bigfloat&, bigfloat*);
	friend void bin_scalprod_bfl(bigfloat&, bigfloat*, bigfloat*);
	friend void bin_scalquad_bfl(bigfloat&, bigfloat*);

	friend void bin_l1_norm_bfl(bigfloat&, bigfloat*);
	friend void bin_l2_norm_bfl(bigfloat&, bigfloat*);


//
// for bigint vectors
//
	friend void bin_assign_bin(bigint*, bigint*);
	friend void bin_assign_zero_bin(bigint*);
	friend void bin_add_bin(bigint*, bigint*, bigint*);
	friend void bin_subtract_bin(bigint*, bigint*, bigint*);

	friend void bin_scalmul_bin(bigint*, const bigint&, bigint*);
	friend void bin_scalprod_bin(bigint&, bigint*, bigint*);
	friend void bin_scalquad_bin(bigint&, bigint*);

	friend void bin_l1_norm_bin(bigint&, bigint*);
	friend void bin_l2_norm_bin(bigint&, bigint*);


	friend void bin_swap_dbl(double*& a, double*& b);
	friend void bin_swap_bfl(bigfloat*& a, bigfloat*& b);
	friend void bin_swap_bin(bigint*& a, bigint*& b);


};


// friend functions

inline
bigint_lattice_gensys lin_gen_system(const bigint_lattice_gensys& L,
				     lidia_size_t& rank)
{
	bigint_lattice_gensys TL(L);
	bigint_lattice_gensys T;
	TL.lin_gen_system(T, rank);
	return (T);
}

//
// Variation of the Schnorr - Euchner - lllfp using doubles for
// paramter = 1 and x*doubleprecision for parameter = x
//

inline
bigint_lattice_gensys lll_gensys(const bigint_lattice_gensys& L,
				 lidia_size_t& rank,
				 sdigit factor = 1)
{
	bigint_lattice_gensys TL(L);
	TL.lll_gensys(rank, factor);
	return (TL);
}

inline
bigint_lattice_gensys lll_trans_gensys(const bigint_lattice_gensys& L,
				       ring_matrix< bigint > & T,
				       lidia_size_t& rank,
				       sdigit factor = 1)
{
	bigint_lattice_gensys LT(L);
	LT.lll_trans_gensys(T, rank, factor);
	return (LT);
}


//
// Modified lll for genertating systems
// of the form n x n+1
//
// Returns vector of dimension n
// where you can find the relations
//

inline
bigint_lattice_gensys mlll(const bigint_lattice_gensys& L,
			       base_vector< bigint > & v)
{
	bigint_lattice_gensys TL(L);
	TL.mlll(v);
	return(TL);
}

inline
bigint_lattice_gensys mlll(const bigint_lattice_gensys& L,
			   bigint*& v)
{
	bigint_lattice_gensys TL(L);
	TL.mlll(v);
	return(TL);
}


//
// For double vectors
//
void bin_assign_dbl(double*, double*);
void bin_assign_zero_dbl(double*);
void bin_add_dbl(double*, double*, double*);
void bin_subtract_dbl(double*, double*, double*);

void bin_scalmul_dbl(double*, const double&, double*);
void bin_scalprod_dbl(double&, double*, double*);
void bin_scalquad_dbl(double&, double*);

void bin_l1_norm_dbl(double&, double*);
void bin_l2_norm_dbl(double&, double*);

//
// for bigfloat vectors
//
void bin_assign_bfl(bigfloat*, bigfloat*);
void bin_assign_zero_bfl(bigfloat*);
void bin_add_bfl(bigfloat*, bigfloat*, bigfloat*);
void bin_subtract_bfl(bigfloat*, bigfloat*, bigfloat*);

void bin_scalmul_bfl(bigfloat*, const bigfloat&, bigfloat*);
void bin_scalprod_bfl(bigfloat&, bigfloat*, bigfloat*);
void bin_scalquad_bfl(bigfloat&, bigfloat*);

void bin_l1_norm_bfl(bigfloat&, bigfloat*);
void bin_l2_norm_bfl(bigfloat&, bigfloat*);


//
// for bigint vectors
//
void bin_assign_bin(bigint*, bigint*);
void bin_assign_zero_bin(bigint*);
void bin_add_bin(bigint*, bigint*, bigint*);
void bin_subtract_bin(bigint*, bigint*, bigint*);

void bin_scalmul_bin(bigint*, const bigint&, bigint*);
void bin_scalprod_bin(bigint&, bigint*, bigint*);
void bin_scalquad_bin(bigint&, bigint*);

void bin_l1_norm_bin(bigint&, bigint*);
void bin_l2_norm_bin(bigint&, bigint*);


inline
void bin_swap_dbl(double*& a, double*& b)
{
	debug_handler("bigint_lattice_gensys", "void bin_swap_dbl(a, b)");
	register double* temp = a;
	a = b;
	b = temp;
}

inline
void bin_swap_bfl(bigfloat*& a, bigfloat*& b)
{
	debug_handler("bigint_lattice_gensys", "bin_swap_bfl(a, b)");
	register bigfloat* temp = a;
	a = b;
	b = temp;
}

inline
void bin_swap_bin(bigint*& a, bigint*& b)
{
	debug_handler("bigint_lattice_gensys", "bin_swap_bin(a, b)");
	register bigint* temp = a;
	a = b;
	b = temp;
}



//
// for double vectors
//

inline void
bin_assign_dbl(double* a, double* b)
{
	register lidia_size_t	i;

	debug_handler("bigint_lattice_gensys", "void bin_assign_dbl(a, b)");
	for (i = bigint_lattice_gensys::vectsize-1; i >= 0; i--)
		a[i] = b[i];
}



inline void
bin_assign_zero_dbl(double* a)
{
	register lidia_size_t	i;

	debug_handler("bigint_lattice_gensys", "void bin_assign_zero_dbl(a)");
	for (i = bigint_lattice_gensys::vectsize-1; i >= 0; i--)
		a[i] = 0;
}



inline void
bin_add_dbl(double* c, double* a, double* b)
{
	register lidia_size_t	i;

	debug_handler("bigint_lattice_gensys", "void bin_add_dbl(c, a, b)");
	for (i = bigint_lattice_gensys::vectsize-1; i >= 0; i--)
		c[i] = a[i]+b[i];
}



inline void
bin_subtract_dbl(double* c, double* a, double* b)
{
	register lidia_size_t	i;

	debug_handler("bigint_lattice_gensys", "void bin_subtract_dbl(c, a, b)");
	for (i = bigint_lattice_gensys::vectsize-1; i >= 0; i--)
		c[i] = a[i]-b[i];
}



inline void
bin_scalmul_dbl(double* b, const double& d, double* a)
{
	register lidia_size_t	i;

	debug_handler("bigint_lattice_gensys", "void bin_scalmul_dbl(b, d, a)");
	for (i = bigint_lattice_gensys::vectsize-1; i >= 0; i--)
		b[i] = a[i]*d;
}



inline void
bin_scalprod_dbl(double& res, double* a, double* b)
{
	register lidia_size_t	i;

	debug_handler("bigint_lattice_gensys", "void bin_scalprod_dbl(res, a, b)");
	res = 0;
	for (i = bigint_lattice_gensys::vectsize-1; i >= 0; i--)
		res += a[i]*b[i];
}



inline void
bin_scalquad_dbl(double& res, double* a)
{
	register lidia_size_t	i;

	debug_handler("bigint_lattice_gensys", "void bin_scalquad_dbl(res, a, b)");
	res = 0;
	for (i = bigint_lattice_gensys::vectsize-1; i >= 0; i--)
		res += a[i]*a[i];
}



inline void
bin_l2_norm_dbl(double& norm, double* v)
{
	register lidia_size_t	i;

	debug_handler("bigint_lattice_gensys", "void bin_l2_norm_dbl(norm, v)");
	norm = 0; // Initialisation
	for (i = bigint_lattice_gensys::vectsize-1; i >= 0; i--)
		norm += v[i]*v[i];
	norm = std::sqrt(norm);
}



inline void
bin_l1_norm_dbl(double& norm, double* v)
{
	register lidia_size_t	i;

	debug_handler("bigint_lattice_gensys", "void bin_l1_norm_dbl(norm, v)");
	norm = 0; // Initialisation
	for (i = bigint_lattice_gensys::vectsize-1; i >= 0; i--)
		norm += std::fabs(v[i]);
}

//
// for bigfloat vectors
//

inline void
bin_assign_bfl(bigfloat* a, bigfloat* b)
{
	register lidia_size_t	i;

	debug_handler("bigint_lattice_gensys", "bin_assign_bfl(a, b)");
	for (i = bigint_lattice_gensys::vectsize-1; i >= 0; i--)
		a[i].assign(b[i]);
}



inline void
bin_assign_zero_bfl(bigfloat* a)
{
	register lidia_size_t	i;

	debug_handler("bigint_lattice_gensys", "bin_assign_zero_bfl(a)");
	for (i = bigint_lattice_gensys::vectsize-1; i >= 0; i--)
		a[i].assign_zero();
}



inline void
bin_add_bfl(bigfloat* c, bigfloat* a, bigfloat* b)
{
	register lidia_size_t	i;

	debug_handler("bigint_lattice_gensys", "bin_add_bfl(c, a, b)");
	for (i = bigint_lattice_gensys::vectsize-1; i >= 0; i--)
		LiDIA::add(c[i], a[i], b[i]);
}



inline void
bin_subtract_bfl(bigfloat* c, bigfloat* a, bigfloat* b)
{
	register lidia_size_t	i;

	debug_handler("bigint_lattice_gensys", "bin_subtract_bfl(c, a, b)");
	for (i = bigint_lattice_gensys::vectsize-1; i >= 0; i--)
		LiDIA::subtract(c[i], a[i], b[i]);
}



inline void
bin_scalmul_bfl(bigfloat* b, const bigfloat& d, bigfloat* a)
{
	register lidia_size_t	i;

	debug_handler("bigint_lattice_gensys", "bin_scalmul_bfl(b, d, a)");
	for (i = bigint_lattice_gensys::vectsize-1; i >= 0; i--)
		LiDIA::multiply(b[i], a[i], d);
}



inline void
bin_scalprod_bfl(bigfloat& res, bigfloat* a, bigfloat* b)
{
	register lidia_size_t	i;

	debug_handler("bigint_lattice_gensys", "bin_scalprod_bfl(res, a, b)");
	res.assign_zero();
	for (i = bigint_lattice_gensys::vectsize-1; i >= 0; i--) {
		LiDIA::multiply(bigint_lattice_gensys::vectbflz, a[i], b[i]);
		LiDIA::add(res, res, bigint_lattice_gensys::vectbflz);
	}
}



inline void
bin_scalquad_bfl(bigfloat& res, bigfloat* a)
{
	register lidia_size_t	i;

	debug_handler("bigint_lattice_gensys", "bin_scalquad_bfl(res, a)");
	res.assign_zero();
	for (i = bigint_lattice_gensys::vectsize-1; i >= 0; i--) {
		LiDIA::square(bigint_lattice_gensys::vectbflz, a[i]);
		LiDIA::add(res, res, bigint_lattice_gensys::vectbflz);
	}
}



inline void
bin_l1_norm_bfl(bigfloat& norm, bigfloat* v)
{
	register lidia_size_t	i;

	debug_handler("bigint_lattice_gensys", "bin_l1_norm_bfl(norm, v)");
	norm.assign_zero(); // Initialisation
	for (i = bigint_lattice_gensys::vectsize-1; i >= 0; i--) {
		bigint_lattice_gensys::vectbflz.assign(abs(v[i]));
		LiDIA::add(norm, norm, bigint_lattice_gensys::vectbflz);
	}
}



inline void
bin_l2_norm_bfl(bigfloat& norm, bigfloat* v)
{
	register lidia_size_t	i;

	debug_handler("bigint_lattice_gensys", "bin_l2_norm_bfl(norm, v)");
	norm.assign_zero(); // Initialisation
	for (i = bigint_lattice_gensys::vectsize-1; i >= 0; i--) {
		LiDIA::square(bigint_lattice_gensys::vectbflz, v[i]);
		LiDIA::add(norm, norm, bigint_lattice_gensys::vectbflz);
	}
	sqrt(norm, norm);
}



//
// for bigint vectors
//

inline void
bin_assign_bin(bigint* a, bigint* b)
{
	register lidia_size_t	i;

	debug_handler("bigint_lattice_gensys", "bin_assign_bin(a, b)");
	for (i = bigint_lattice_gensys::vectsize-1; i >= 0; i--)
		a[i].assign(b[i]);
}



inline void
bin_assign_zero_bin(bigint* a)
{
	register lidia_size_t	i;

	debug_handler("bigint_lattice_gensys", "bin_assign_zero_bin(a)");
	for (i = bigint_lattice_gensys::vectsize-1; i >= 0; i--)
		a[i].assign_zero();
}



inline void
bin_add_bin(bigint* c, bigint* a, bigint* b)
{
	register lidia_size_t	i;

	debug_handler("bigint_lattice_gensys", "bin_add_bin(c, a, b)");
	for (i = bigint_lattice_gensys::vectsize-1; i >= 0; i--)
		LiDIA::add(c[i], a[i], b[i]);
}



inline void
bin_subtract_bin(bigint* c, bigint* a, bigint* b)
{
	register lidia_size_t	i;

	debug_handler("bigint_lattice_gensys", "bin_subtract_bin(c, a, b)");
	for (i = bigint_lattice_gensys::vectsize-1; i >= 0; i--)
		LiDIA::subtract(c[i], a[i], b[i]);
}



inline void
bin_scalmul_bin(bigint* b, const bigint& d, bigint* a)
{
	register lidia_size_t	i;

	debug_handler("bigint_lattice_gensys", "bin_scalmul_bin(b, d, a)");
	for (i = bigint_lattice_gensys::vectsize-1; i >= 0; i--)
		LiDIA::multiply(b[i], a[i], d);
}



inline void
bin_scalprod_bin(bigint& res, bigint* a, bigint* b)
{
	register lidia_size_t	i;

	debug_handler("bigint_lattice_gensys", "bin_scalprod_bin(res, a, b)");
	res.assign_zero();
	for (i = bigint_lattice_gensys::vectsize-1; i >= 0; i--) {
		LiDIA::multiply(bigint_lattice_gensys::vectbinz, a[i], b[i]);
		LiDIA::add(res, res, bigint_lattice_gensys::vectbinz);
	}
}



inline void
bin_scalquad_bin(bigint& res, bigint* a)
{
	register lidia_size_t	i;

	debug_handler("bigint_lattice_gensys", "bin_scalquad_bin(res, a)");
	res.assign_zero();
	for (i = bigint_lattice_gensys::vectsize-1; i >= 0; i--) {
		LiDIA::square(bigint_lattice_gensys::vectbinz, a[i]);
		LiDIA::add(res, res, bigint_lattice_gensys::vectbinz);
	}
}



inline void
bin_l2_norm_bin(bigint& norm, bigint* v)
{
	register lidia_size_t	i;

	debug_handler("bigint_lattice_gensys", "bin_l2_norm_bin(norm, v)");
	norm.assign_zero(); // Initialisation
	for (i = bigint_lattice_gensys::vectsize-1; i >= 0; i--) {
		LiDIA::square(bigint_lattice_gensys::vectbinz, v[i]);
		LiDIA::add(norm, norm, bigint_lattice_gensys::vectbinz);
	}
	//  sqrt(norm, norm);
}



inline void
bin_l1_norm_bin(bigint& norm, bigint* v)
{
	register lidia_size_t	i;

	debug_handler("bigint_lattice_gensys", "bin_l1_norm_bin(norm, v)");
	norm.assign_zero(); // Initialisation
	for (i = bigint_lattice_gensys::vectsize-1; i >= 0; i--) {
		bigint_lattice_gensys::vectbinz.assign(abs(v[i]));
		LiDIA::add(norm, norm, bigint_lattice_gensys::vectbinz);
	}
}



#ifdef LIDIA_NAMESPACE
}	// end of namespace LiDIA
# undef IN_NAMESPACE_LIDIA
#endif



#endif	// LIDIA_BI_LATTICE_GENSYS_H_GUARD_
