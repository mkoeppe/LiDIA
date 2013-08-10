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


#ifndef LIDIA_BF_LATTICE_GENSYS_H_GUARD_
#define LIDIA_BF_LATTICE_GENSYS_H_GUARD_


#ifndef LIDIA_MATH_MATRIX_H_GUARD_
# include	"LiDIA/math_matrix.h"
#endif
#ifndef LIDIA_BI_LATTICE_GENSYS_H_GUARD_
# include	"LiDIA/lattices/bi_lattice_gensys.h"
#endif
#ifndef LIDIA_BIGINT_H_GUARD_
# include	"LiDIA/bigint.h"
#endif
#ifndef LIDIA_BIGFLOAT_H_GUARD_
# include	"LiDIA/bigfloat.h"
#endif
#ifndef LIDIA_BASE_VECTOR_H_GUARD_
# include	"LiDIA/base_vector.h"
#endif
#ifndef LIDIA_TIMER_H_GUARD_
# include	"LiDIA/timer.h"
#endif



#ifdef LIDIA_NAMESPACE
namespace LiDIA {
# define IN_NAMESPACE_LIDIA
#endif



#ifndef rint
# define rint(x) (((x)-std::floor(x)) < 0.5 ? std::floor(x) : std::ceil(x))
#endif

typedef int (*bfl_cmp)(const bigfloat*, const bigfloat*, long);



class bigfloat_lattice_gensys:public math_matrix< bigfloat >
{
protected :
//
// configuration of lattice
//
	sdigit         y_nom;
	sdigit         y_denom;
	double         y_par;

	bool           trans_flag;

//
// algorithm - information
//
	sdigit         reduction_steps;
	sdigit         swapping_steps;
	sdigit         correction_steps;

//
// static variables for temporary use to increase speed
//
	static bigfloat tempmz0;
	static bigfloat tempmz1;
	static bigfloat tempmz2;
	static bigfloat tempmz3;
	static bigfloat tempmz4;
	static bigfloat ergmz;

	static double   vectdblz;
	static bigint   vectbinz;
	static bigfloat vectbflz;

	static lidia_size_t vectsize;

public :

//
// Constructor / Destructors
//
	bigfloat_lattice_gensys();
//    bigfloat_lattice_gensys(lidia_size_t);
	bigfloat_lattice_gensys(lidia_size_t, lidia_size_t);
	bigfloat_lattice_gensys(lidia_size_t, lidia_size_t, const bigfloat**);
	bigfloat_lattice_gensys(const math_matrix< bigfloat > &);
	bigfloat_lattice_gensys(const bigfloat_lattice_gensys&);
	~bigfloat_lattice_gensys();

//
// Assignments
//
	void assign(const bigfloat_lattice_gensys&);
	bigfloat_lattice_gensys& operator = (const bigfloat_lattice_gensys&);
	friend void assign(bigfloat_lattice_gensys&, const bigfloat_lattice_gensys&);

//
// Set or get reduction info, paramters
//
	void set_reduction_parameter(double);
	void set_reduction_parameter(sdigit, sdigit);
	double get_reduction_parameter();
	void get_reduction_parameter(sdigit&, sdigit&);

//
// Set representation direction
// reduce column by column or row by row
//
	void set_representation_columns();
	void set_representation_rows();
	bool get_representation()
	{
		return (trans_flag);
	}

//
// get detailed information about lll - algorithm
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
	void sort_vectors(bfl_cmp);

//
// Buchmann - Kessler Algorthm for
// linaer generating systems
//
	void lin_gen_system(math_matrix< bigint > &, lidia_size_t&);
	void lin_gen_system(math_matrix< bigfloat > &, lidia_size_t&);
	friend bigint_lattice_gensys lin_gen_system(const bigfloat_lattice_gensys& L,
						    lidia_size_t& rank);

//
// Variation of the Schnorr - Euchner - lllfp using doubles for
// parameter = 1 and x*doubleprec for parameter = x
//
	void lll_gensys(lidia_size_t&, sdigit = 1);
	void lll_gensys(const math_matrix< bigfloat > &, lidia_size_t&, sdigit = 1);
	friend bigfloat_lattice_gensys lll_gensys(const bigfloat_lattice_gensys& L,
						  lidia_size_t& rank,
						  sdigit factor );

	void lll_trans_gensys(math_matrix< bigint > &, lidia_size_t&, sdigit = 1);
	void lll_trans_gensys(math_matrix< bigfloat > &, lidia_size_t&, sdigit = 1);
	friend bigfloat_lattice_gensys lll_trans_gensys(const bigfloat_lattice_gensys& L,
							math_matrix< bigint > & T,
							lidia_size_t& rank,
							sdigit factor);

	friend bigfloat_lattice_gensys lll_trans_gensys(const bigfloat_lattice_gensys& L,
							math_matrix< bigfloat > & T,
							lidia_size_t& rank,
							sdigit factor);

	void lll_var_gensys(lidia_size_t&, sdigit = 1);
	void lll_var_gensys(const math_matrix< bigfloat > &, lidia_size_t&, sdigit = 1);
	friend bigfloat_lattice_gensys lll_var_gensys(const bigfloat_lattice_gensys& L,
						      lidia_size_t& rank,
						      sdigit factor );


	void lll_trans_var_gensys(math_matrix< bigint > &, lidia_size_t&, sdigit = 1);
	void lll_trans_var_gensys(math_matrix< bigfloat > &, lidia_size_t&, sdigit = 1);
	friend bigfloat_lattice_gensys lll_trans_var_gensys(const bigfloat_lattice_gensys& L,
							    math_matrix< bigint > & T,
							    lidia_size_t& rank,
							    sdigit factor);

	friend bigfloat_lattice_gensys lll_trans_var_gensys(const bigfloat_lattice_gensys& L,
							    math_matrix< bigfloat > & T,
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
	void mlll(const math_matrix< bigfloat > &, base_vector< bigint > &);
	friend bigfloat_lattice_gensys mlll(const bigfloat_lattice_gensys& L,
					    base_vector< bigint > & v);

	void mlll(bigint *&);
	void mlll(const math_matrix< bigfloat > &, bigint *&);
	friend bigfloat_lattice_gensys mlll(const bigfloat_lattice_gensys& L,
					    bigint*& v);

	void mlll(base_vector< bigfloat > &);
	void mlll(const math_matrix< bigfloat > &, base_vector< bigfloat > &);
	friend bigfloat_lattice_gensys mlll(const bigfloat_lattice_gensys& L,
					    base_vector< bigfloat > & v);

	void mlll(bigfloat *&);
	void mlll(const math_matrix< bigfloat > &, bigfloat *&);
	friend bigfloat_lattice_gensys mlll(const bigfloat_lattice_gensys& L,
					    bigfloat*& v);

//
// Compute a basis from lattice A and Transformation lattice T
//
	void compute_basis(const math_matrix< bigfloat > & A, const math_matrix< bigint > & T)
	{
		debug_handler("bigfloat_lattice_gensys", "compute_basis(A, T)");
		math_matrix< bigfloat > Bi(columns, columns);
		bigfloatify(Bi, T);
		compute_basis(A, Bi);
	}

	void compute_basis(const math_matrix< bigfloat > & A, const math_matrix< bigfloat > & T)
	{
		debug_handler("bigfloat_lattice_gensys", "compute_basis(A, T)");
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
	void assign_the_rest(const bigfloat_lattice_gensys&);

	void bigfloatify(math_matrix< bigfloat > &, const math_matrix< bigint > &);

	sdigit compute_read_precision();
	sdigit compute_precision();

	void alpha_compute(bigfloat&);
	void gamma_compute(bigfloat&, sdigit);
	void zwei_pot_q_compute(bigfloat&, sdigit&, bigfloat&);

//
// Very funny, helpful, speedup
//
	void Tr_trans_swap(bigfloat_lattice_gensys&);

//
// Dimension Checking
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
	void Tr_sort_vectors(bfl_cmp);

	void Tr_lll();
	void Tr_lll_dbl();
	void Tr_lll_bfl(sdigit);

	void Tr_lll_var_dbl(math_matrix< bigint > &, sdigit);
	void Tr_lll_var_bfl(math_matrix< bigint > &, sdigit, sdigit);

	void Tr_lll_deep_insert_dbl(sdigit);
	void Tr_lll_deep_insert_bfl(sdigit, sdigit);

	void Tr_lll_trans(math_matrix< bigint > &);
	void Tr_lll_trans_dbl(math_matrix< bigint > &);
	void Tr_lll_trans_bfl(math_matrix< bigint > &, sdigit);

	void Tr_lll_trans_var_dbl(math_matrix< bigint > &, math_matrix< bigint > &,
				  sdigit);
	void Tr_lll_trans_var_bfl(math_matrix< bigint > &, math_matrix< bigint > &,
				  sdigit, sdigit);

	void Tr_mlll_bfl(bigint*);

	void Tr_lin_gen_system(math_matrix< bigint > &, lidia_size_t&);

	void Tr_lll_dbl_gensys(lidia_size_t&);
	void Tr_lll_bfl_gensys(lidia_size_t&, sdigit);

	void Tr_lll_var_dbl_gensys(math_matrix< bigint > &,
				   sdigit, lidia_size_t&);
	void Tr_lll_var_bfl_gensys(math_matrix< bigint > &,
				   sdigit, lidia_size_t&, sdigit);

	void Tr_lll_trans_dbl_gensys(math_matrix< bigint > &, lidia_size_t&);
	void Tr_lll_trans_bfl_gensys(math_matrix< bigint > &, lidia_size_t&, sdigit);

	void Tr_lll_trans_var_dbl_gensys(math_matrix< bigint > &,
					 math_matrix< bigint > &,
					 sdigit, lidia_size_t&);
	void Tr_lll_trans_var_bfl_gensys(math_matrix< bigint > &,
					 math_matrix< bigint > &,
					 sdigit, lidia_size_t&, sdigit);

//
// "Vector" - Operations
//


//
// For double vectors
//
	friend void bfl_assign_dbl(double*, double*);
	friend void bfl_assign_zero_dbl(double*);
	friend void bfl_add_dbl(double*, double*, double*);
	friend void bfl_subtract_dbl(double*, double*, double*);

	friend void bfl_scalmul_dbl(double*, const double&, double*);
	friend void bfl_scalprod_dbl(double&, double*, double*);
	friend void bfl_scalquad_dbl(double&, double*);

	friend void bfl_l1_norm_dbl(double&, double*);
	friend void bfl_l2_norm_dbl(double&, double*);

//
// for bigfloat vectors
//
	friend void bfl_assign_bfl(bigfloat*, bigfloat*);
	friend void bfl_assign_zero_bfl(bigfloat*);
	friend void bfl_add_bfl(bigfloat*, bigfloat*, bigfloat*);
	friend void bfl_subtract_bfl(bigfloat*, bigfloat*, bigfloat*);

	friend void bfl_scalmul_bfl(bigfloat*, const bigfloat&, bigfloat*);
	friend void bfl_scalprod_bfl(bigfloat&, bigfloat*, bigfloat*);
	friend void bfl_scalquad_bfl(bigfloat&, bigfloat*);

	friend void bfl_l1_norm_bfl(bigfloat&, bigfloat*);
	friend void bfl_l2_norm_bfl(bigfloat&, bigfloat*);


//
// for bigint vectors
//
	friend void bfl_assign_bin(bigint*, bigint*);
	friend void bfl_assign_zero_bin(bigint*);
	friend void bfl_add_bin(bigint*, bigint*, bigint*);
	friend void bfl_subtract_bin(bigint*, bigint*, bigint*);

	friend void bfl_scalmul_bin(bigint*, const bigint&, bigint*);
	friend void bfl_scalprod_bin(bigint&, bigint*, bigint*);
	friend void bfl_scalquad_bin(bigint&, bigint*);

	friend void bfl_l1_norm_bin(bigint&, bigint*);
	friend void bfl_l2_norm_bin(bigint&, bigint*);



	friend void bfl_swap_dbl(double*& a, double*& b);

	friend void bfl_swap_bfl(bigfloat*& a, bigfloat*& b);

	friend void bfl_swap_bin(bigint*& a, bigint*& b);

};


// friend functions
 
void assign(bigfloat_lattice_gensys&, const bigfloat_lattice_gensys&);

//
// Buchmann - Kessler Algorthm for
// linaer generating systems
//

inline
bigint_lattice_gensys lin_gen_system(const bigfloat_lattice_gensys& L,
					 lidia_size_t& rank)
{
	bigfloat_lattice_gensys TL(L);
	bigint_lattice_gensys T;
	TL.lin_gen_system(T, rank);
	return (T);
}

//
// Variation of the Schnorr - Euchner - lllfp using doubles for
// parameter = 1 and x*doubleprec for parameter = x
//

inline
bigfloat_lattice_gensys lll_gensys(const bigfloat_lattice_gensys& L,
				       lidia_size_t& rank,
				       sdigit factor = 1)
{
	bigfloat_lattice_gensys TL(L);
	TL.lll_gensys(rank, factor);
	return (TL);
}

inline
bigfloat_lattice_gensys lll_trans_gensys(const bigfloat_lattice_gensys& L,
					     math_matrix< bigint > & T,
					     lidia_size_t& rank,
					     sdigit factor = 1)
{
	bigfloat_lattice_gensys LT(L);
	LT.lll_trans_gensys(T, rank, factor);
	return (LT);
}

inline
bigfloat_lattice_gensys lll_trans_gensys(const bigfloat_lattice_gensys& L,
					     math_matrix< bigfloat > & T,
					     lidia_size_t& rank,
					     sdigit factor = 1)
{
	bigfloat_lattice_gensys LT(L);
	LT.lll_trans_gensys(T, rank, factor);
	return (LT);
}

inline
bigfloat_lattice_gensys lll_var_gensys(const bigfloat_lattice_gensys& L,
					   lidia_size_t& rank,
					   sdigit factor = 1)
{
	bigfloat_lattice_gensys TL(L);
	TL.lll_var_gensys(rank, factor);
	return (TL);
}

inline
bigfloat_lattice_gensys lll_trans_var_gensys(const bigfloat_lattice_gensys& L,
						 math_matrix< bigint > & T,
						 lidia_size_t& rank,
						 sdigit factor = 1)
{
	bigfloat_lattice_gensys LT(L);
	LT.lll_trans_var_gensys(T, rank, factor);
	return (LT);
}

inline
bigfloat_lattice_gensys lll_trans_var_gensys(const bigfloat_lattice_gensys& L,
						 math_matrix< bigfloat > & T,
						 lidia_size_t& rank,
						 sdigit factor = 1)
{
	bigfloat_lattice_gensys LT(L);
	LT.lll_trans_var_gensys(T, rank, factor);
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
bigfloat_lattice_gensys mlll(const bigfloat_lattice_gensys& L,
				 base_vector< bigint > & v)
{
	bigfloat_lattice_gensys TL(L);
	TL.mlll(L, v);
	return (TL);
}

inline
bigfloat_lattice_gensys mlll(const bigfloat_lattice_gensys& L,
				 bigint*& v)
{
	bigfloat_lattice_gensys TL(L);
	TL.mlll(L, v);
	return (TL);
}

inline
bigfloat_lattice_gensys mlll(const bigfloat_lattice_gensys& L,
				 base_vector< bigfloat > & v)
{
	bigfloat_lattice_gensys TL(L);
	TL.mlll(L, v);
	return (TL);
}

inline
bigfloat_lattice_gensys mlll(const bigfloat_lattice_gensys& L,
				 bigfloat*& v)
{
	bigfloat_lattice_gensys TL(L);
	TL.mlll(L, v);
	return (TL);
}


//
// For double vectors
//
void bfl_assign_dbl(double*, double*);
void bfl_assign_zero_dbl(double*);
void bfl_add_dbl(double*, double*, double*);
void bfl_subtract_dbl(double*, double*, double*);

void bfl_scalmul_dbl(double*, const double&, double*);
void bfl_scalprod_dbl(double&, double*, double*);
void bfl_scalquad_dbl(double&, double*);

void bfl_l1_norm_dbl(double&, double*);
void bfl_l2_norm_dbl(double&, double*);

//
// for bigfloat vectors
//
void bfl_assign_bfl(bigfloat*, bigfloat*);
void bfl_assign_zero_bfl(bigfloat*);
void bfl_add_bfl(bigfloat*, bigfloat*, bigfloat*);
void bfl_subtract_bfl(bigfloat*, bigfloat*, bigfloat*);

void bfl_scalmul_bfl(bigfloat*, const bigfloat&, bigfloat*);
void bfl_scalprod_bfl(bigfloat&, bigfloat*, bigfloat*);
void bfl_scalquad_bfl(bigfloat&, bigfloat*);

void bfl_l1_norm_bfl(bigfloat&, bigfloat*);
void bfl_l2_norm_bfl(bigfloat&, bigfloat*);


//
// for bigint vectors
//
void bfl_assign_bin(bigint*, bigint*);
void bfl_assign_zero_bin(bigint*);
void bfl_add_bin(bigint*, bigint*, bigint*);
void bfl_subtract_bin(bigint*, bigint*, bigint*);

void bfl_scalmul_bin(bigint*, const bigint&, bigint*);
void bfl_scalprod_bin(bigint&, bigint*, bigint*);
void bfl_scalquad_bin(bigint&, bigint*);

void bfl_l1_norm_bin(bigint&, bigint*);
void bfl_l2_norm_bin(bigint&, bigint*);



inline
void bfl_swap_dbl(double*& a, double*& b)
{
	debug_handler("bigint_lattice_gensys", "void bfl_swap_dbl(a, b)");
	register double* temp = a;
	a = b;
	b = temp;
}

inline
void bfl_swap_bfl(bigfloat*& a, bigfloat*& b)
{
	debug_handler("bigint_lattice_gensys", "bfl_swap_bfl(a, b)");
	register bigfloat* temp = a;
	a = b;
	b = temp;
}

inline
void bfl_swap_bin(bigint*& a, bigint*& b)
{
	debug_handler("bigfloat_lattice_gensys", "bfl_swap_bin(a, b)");
	register bigint* temp = a;
	a = b;
	b = temp;
}


//
// for double vectors
//

inline void
bfl_assign_dbl(double* a, double* b)
{
	register lidia_size_t i;

	debug_handler("bigfloat_lattice_gensys", "void bfl_assign_dbl(a, b)");
	for (i = bigfloat_lattice_gensys::vectsize-1; i >= 0; i--)
		a[i] = b[i];
}



inline void
bfl_assign_zero_dbl(double* a)
{
	register lidia_size_t i;

	debug_handler("bigfloat_lattice_gensys", "void bfl_assign_zero_dbl(a)");
	for (i = bigfloat_lattice_gensys::vectsize-1; i >= 0; i--)
		a[i] = 0;
}



inline void
bfl_add_dbl(double* c, double* a, double* b)
{
	register lidia_size_t i;

	debug_handler("bigfloat_lattice_gensys", "void bfl_add_dbl(c, a, b)");
	for (i = bigfloat_lattice_gensys::vectsize-1; i >= 0; i--)
		c[i] = a[i]+b[i];
}



inline void
bfl_subtract_dbl(double* c, double* a, double* b)
{
	register lidia_size_t i;

	debug_handler("bigfloat_lattice_gensys", "void bfl_subtract_dbl(c, a, b)");
	for (i = bigfloat_lattice_gensys::vectsize-1; i >= 0; i--)
		c[i] = a[i]-b[i];
}



inline void
bfl_scalmul_dbl(double* b, const double& d, double* a)
{
	register lidia_size_t i;

	debug_handler("bigfloat_lattice_gensys", "void bfl_scalmul_dbl(b, d, a)");
	for (i = bigfloat_lattice_gensys::vectsize-1; i >= 0; i--)
		b[i] = a[i]*d;
}



inline void
bfl_scalprod_dbl(double& res, double* a, double* b)
{
	register lidia_size_t i;

	debug_handler("bigfloat_lattice_gensys", "void bfl_scalprod_dbl(res, a, b)");
	res = 0;
	for (i = bigfloat_lattice_gensys::vectsize-1; i >= 0; i--)
		res += a[i]*b[i];
}



inline void
bfl_scalquad_dbl(double& res, double* a)
{
	register lidia_size_t i;

	debug_handler("bigfloat_lattice_gensys", "void bfl_scalquad_dbl(res, a, b)");
	res = 0;
	for (i = bigfloat_lattice_gensys::vectsize-1; i >= 0; i--)
		res += a[i]*a[i];
}



inline void
bfl_l2_norm_dbl(double& norm, double* v)
{
	register lidia_size_t i;

	debug_handler("bigfloat_lattice_gensys", "void bfl_l2_norm_dbl(norm, v)");
	norm = 0; // Initialisation
	for (i = bigfloat_lattice_gensys::vectsize-1; i >= 0; i--)
		norm += v[i]*v[i];
	norm = std::sqrt(norm);
}



inline void
bfl_l1_norm_dbl(double& norm, double* v)
{
	register lidia_size_t i;

	debug_handler("bigfloat_lattice_gensys", "void bfl_l1_norm_dbl(norm, v)");
	norm = 0; // Initialisation
	for (i = bigfloat_lattice_gensys::vectsize-1; i >= 0; i--)
		norm += std::fabs(v[i]);
}



//
// for bigfloat vectors
//

inline void
bfl_assign_bfl(bigfloat* a, bigfloat* b)
{
	register lidia_size_t i;

	debug_handler("bigfloat_lattice_gensys", "bfl_assign_bfl(a, b)");
	for (i = bigfloat_lattice_gensys::vectsize-1; i >= 0; i--)
		a[i].assign(b[i]);
}



inline void
bfl_assign_zero_bfl(bigfloat* a)
{
	register lidia_size_t i;

	debug_handler("bigfloat_lattice_gensys", "bfl_assign_zero_bfl(a)");
	for (i = bigfloat_lattice_gensys::vectsize-1; i >= 0; i--)
		a[i].assign_zero();
}



inline void
bfl_add_bfl(bigfloat* c, bigfloat* a, bigfloat* b)
{
	register lidia_size_t i;

	debug_handler("bigfloat_lattice_gensys", "bfl_add_bfl(c, a, b)");
	for (i = bigfloat_lattice_gensys::vectsize-1; i >= 0; i--)
		LiDIA::add(c[i], a[i], b[i]);
}



inline void
bfl_subtract_bfl(bigfloat* c, bigfloat* a, bigfloat* b)
{
	register lidia_size_t i;

	debug_handler("bigfloat_lattice_gensys", "bfl_subtract_bfl(c, a, b)");
	for (i = bigfloat_lattice_gensys::vectsize-1; i >= 0; i--)
		LiDIA::subtract(c[i], a[i], b[i]);
}



inline void
bfl_scalmul_bfl(bigfloat* b, const bigfloat& d, bigfloat* a)
{
	register lidia_size_t i;

	debug_handler("bigfloat_lattice_gensys", "bfl_scalmul_bfl(b, d, a)");
	for (i = bigfloat_lattice_gensys::vectsize-1; i >= 0; i--)
		LiDIA::multiply(b[i], a[i], d);
}



inline void
bfl_scalprod_bfl(bigfloat& res, bigfloat* a, bigfloat* b)
{
	register lidia_size_t i;

	debug_handler("bigfloat_lattice_gensys", "bfl_scalprod_bfl(res, a, b)");
	res.assign_zero();
	for (i = bigfloat_lattice_gensys::vectsize-1; i >= 0; i--) {
		LiDIA::multiply(bigfloat_lattice_gensys::vectbflz, a[i], b[i]);
		LiDIA::add(res, res, bigfloat_lattice_gensys::vectbflz);
	}
}



inline void
bfl_scalquad_bfl(bigfloat& res, bigfloat* a)
{
	register lidia_size_t i;

	debug_handler("bigfloat_lattice_gensys", "bfl_scalquad_bfl(res, a)");
	res.assign_zero();
	for (i = bigfloat_lattice_gensys::vectsize-1; i >= 0; i--) {
		LiDIA::square(bigfloat_lattice_gensys::vectbflz, a[i]);
		LiDIA::add(res, res, bigfloat_lattice_gensys::vectbflz);
	}
}



inline void
bfl_l1_norm_bfl(bigfloat& norm, bigfloat* v)
{
	register lidia_size_t i;

	debug_handler("bigfloat_lattice_gensys", "bfl_l1_norm_bfl(norm, v)");
	norm.assign_zero(); // Initialisation
	for (i = bigfloat_lattice_gensys::vectsize-1; i >= 0; i--) {
		bigfloat_lattice_gensys::vectbflz.assign(abs(v[i]));
		LiDIA::add(norm, norm, bigfloat_lattice_gensys::vectbflz);
	}
}



inline void
bfl_l2_norm_bfl(bigfloat& norm, bigfloat* v)
{
	register lidia_size_t i;

	debug_handler("bigfloat_lattice_gensys", "bfl_l2_norm_bfl(norm, v)");
	norm.assign_zero(); // Initialisation
	for (i = bigfloat_lattice_gensys::vectsize-1; i >= 0; i--) {
		LiDIA::square(bigfloat_lattice_gensys::vectbflz, v[i]);
		LiDIA::add(norm, norm, bigfloat_lattice_gensys::vectbflz);
	}
	LiDIA::sqrt(norm, norm);
}



//
// for bigint vectors
//

inline void
bfl_assign_bin(bigint* a, bigint* b)
{
	register lidia_size_t i;

	debug_handler("bigfloat_lattice_gensys", "bfl_assign_bin(a, b)");
	for (i = bigfloat_lattice_gensys::vectsize-1; i >= 0; i--)
		a[i].assign(b[i]);
}



inline void
bfl_assign_zero_bin(bigint* a)
{
	register lidia_size_t i;

	debug_handler("bigfloat_lattice_gensys", "bfl_assign_zero_bin(a)");
	for (i = bigfloat_lattice_gensys::vectsize-1; i >= 0; i--)
		a[i].assign_zero();
}



inline void
bfl_add_bin(bigint* c, bigint* a, bigint* b)
{
	register lidia_size_t i;

	debug_handler("bigfloat_lattice_gensys", "bfl_add_bin(c, a, b)");
	for (i = bigfloat_lattice_gensys::vectsize-1; i >= 0; i--)
		LiDIA::add(c[i], a[i], b[i]);
}



inline void
bfl_subtract_bin(bigint* c, bigint* a, bigint* b)
{
	register lidia_size_t i;

	debug_handler("bigfloat_lattice_gensys", "bfl_subtract_bin(c, a, b)");
	for (i = bigfloat_lattice_gensys::vectsize-1; i >= 0; i--)
		LiDIA::subtract(c[i], a[i], b[i]);
}



inline void
bfl_scalmul_bin(bigint* b, const bigint& d, bigint* a)
{
	register lidia_size_t i;

	debug_handler("bigfloat_lattice_gensys", "bfl_scalmul_bin(b, d, a)");
	for (i = bigfloat_lattice_gensys::vectsize-1; i >= 0; i--)
		LiDIA::multiply(b[i], a[i], d);
}



inline void
bfl_scalprod_bin(bigint& res, bigint* a, bigint* b)
{
	register lidia_size_t i;

	debug_handler("bigint_lattice_gensys", "bfl_scalprod_bin(res, a, b)");
	res.assign_zero();
	for (i = bigfloat_lattice_gensys::vectsize-1; i >= 0; i--) {
		LiDIA::multiply(bigfloat_lattice_gensys::vectbinz, a[i], b[i]);
		LiDIA::add(res, res, bigfloat_lattice_gensys::vectbinz);
	}
}



inline void
bfl_scalquad_bin(bigint& res, bigint* a)
{
	register lidia_size_t i;

	debug_handler("bigfloat_lattice_gensys", "bfl_scalquad_bin(res, a)");
	res.assign_zero();
	for (i = bigfloat_lattice_gensys::vectsize-1; i >= 0; i--) {
		LiDIA::square(bigfloat_lattice_gensys::vectbinz, a[i]);
		LiDIA::add(res, res, bigfloat_lattice_gensys::vectbinz);
	}
}



inline void
bfl_l2_norm_bin(bigint& norm, bigint* v)
{
	register lidia_size_t i;

	debug_handler("bigfloat_lattice_gensys", "bfl_l2_norm_bin(norm, v)");
	norm.assign_zero(); // Initialisation
	for (i = bigfloat_lattice_gensys::vectsize-1; i >= 0; i--) {
		LiDIA::square(bigfloat_lattice_gensys::vectbinz, v[i]);
		LiDIA::add(norm, norm, bigfloat_lattice_gensys::vectbinz);
	}
	//  LiDIA::sqrt(norm, norm);
}



inline void
bfl_l1_norm_bin(bigint& norm, bigint* v)
{
	register lidia_size_t i;

	debug_handler("bigfloat_lattice_gensys", "bfl_l1_norm_bin(norm, v)");
	norm.assign_zero(); // Initialisation
	for (i = bigfloat_lattice_gensys::vectsize-1; i >= 0; i--) {
		bigfloat_lattice_gensys::vectbinz.assign(abs(v[i]));
		LiDIA::add(norm, norm, bigfloat_lattice_gensys::vectbinz);
	}
}



#ifdef LIDIA_NAMESPACE
}	// end of namespace LiDIA
# undef IN_NAMESPACE_LIDIA
#endif



#endif	// LIDIA_BF_LATTICE_GENSYS_H_GUARD_
