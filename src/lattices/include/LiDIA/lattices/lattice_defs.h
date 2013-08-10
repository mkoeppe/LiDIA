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


#ifndef LIDIA_LATTICE_DEFS_H_GUARD_
#define LIDIA_LATTICE_DEFS_H_GUARD_


#ifndef LIDIA_XDOUBLE_H_GUARD_
# include	"LiDIA/xdouble.h"
#endif
#ifndef LIDIA_BIGINT_H_GUARD_
# include	"LiDIA/bigint.h"
#endif
#ifndef LIDIA_BIGFLOAT_H_GUARD_
# include	"LiDIA/bigfloat.h"
#endif
#ifndef LIDIA_P_VECTOR_H_GUARD_
# include	"LiDIA/lattices/p_vector.h"
#endif
#ifndef LIDIA_BIGINT_MATRIX_H_GUARD_
# include	"LiDIA/bigint_matrix.h"
#endif
#ifndef LIDIA_RANDOM_GENERATOR_H_GUARD_
# include	"LiDIA/random_generator.h"
#endif



#ifdef LIDIA_NAMESPACE
namespace LiDIA {
# define IN_NAMESPACE_LIDIA
#endif



class bigint_lattice;
class bigfloat_lattice;



//
// typedef for sortfunction
//

typedef sdigit (*bin_cmp_func)(const bigint*, const bigint*, lidia_size_t);
typedef sdigit (*bfl_cmp_func)(const bigfloat*, const bigfloat*, lidia_size_t);



//
// typedef for scalarproduct
// the bigfloat - version is not allowed to change precision
// (if it's changed then restore after computation)
//
typedef void (*scal_dbl)(double&, double*, double*, lidia_size_t);
typedef void (*scal_xdbl)(xdouble&, xdouble*, xdouble*, lidia_size_t);
typedef void (*scal_bin)(bigint&, bigint*, bigint*, lidia_size_t);
typedef void (*scal_bfl)(bigfloat&, bigfloat*, bigfloat*, lidia_size_t);

typedef struct {
	scal_dbl dbl;
	scal_xdbl xdbl;
	scal_bin bin;
	scal_bfl bfl;
} user_SP;

//
// Information about the algorithms
//
typedef union {
	struct {
		lidia_size_t rank;
		double y;
		sdigit y_nom;
		sdigit y_denom;
		sdigit reduction_steps;
		sdigit correction_steps;
		sdigit swaps;
	} lll;
} lattice_info;

//
// structure definition for internal class use
// to simplify parameters for algorithms
//
typedef struct {
	double y;
	sdigit y_nom;
	sdigit y_denom;
	bool transpose;
	bool alg_trans;
	lidia_size_t rows;
	lidia_size_t columns;
	lidia_size_t rank;
	lidia_size_t real_columns;
} base_alg;


typedef struct {
	sdigit bit_prec;
	sdigit cut_bit_prec;
	sdigit bit_factor;
	sdigit approx_prec;
	sdigit exact_prec;
	sdigit old_prec;
	sdigit read_prec;
} info_dgts;

template< class T >
struct dense_struc
{
	T **value;
	T *delvalue;
	math_matrix< T > *TMatrix;
};

//
// Parameters for dense lattice algorithms
//
template< class T >
struct dense_alg
{
	dense_struc< T > s;
	info_dgts d;
	base_alg b;
};

//
// Parameters for sparse lattice algorithms
//
struct sparse_alg
{
	base_alg b;
};


//
// Vector operations for the needed types using
// the new class p_vector. These are put together to
// have one template argument instead of three
//

template< class E, class A >
struct vector_op
{
	p_vector< A > approx;
	p_vector< E > exact;
};

//
// See above, using function pointer defined scalarproduct
//
template< class E, class A >
struct vector_op_SP
{
	p_vector_SP< A > approx;
	p_vector_SP< E > exact;
};

const sdigit DOUBLE_MANTISSA_BITS = 52;
//const sdigit MANTISSA_CUT = bigint::bits_per_digit();
//MM: Replaced this by a define for SunPro
#ifndef MANTISSA_CUT
#define MANTISSA_CUT bigint::bits_per_digit()
#endif


#define ALG_DEF_BASIS(E, Vect, Var)                                      \
struct LIDIA_CONCAT4(E, Vect, basis, Var) {                               \
  lll_kernel_op< E, double, Vect < E, double >, \
                basis_modules< E, double, Var > > dbl; \
  lll_kernel_op< E, xdouble, Vect < E, xdouble >, \
                basis_modules< E, xdouble, Var > > xdbl; \
  lll_kernel_fu< E, bigfloat, Vect < E, bigfloat >, \
                basis_modules< E, bigfloat, Var > > bfl; \
}

#define ALG_DEF_GENSYS(E, Vect, Var)                                     \
struct LIDIA_CONCAT4(E, Vect, igensys, Var) {                             \
  lll_kernel_op< E, double, Vect < E, double >, \
                gensys_modules< E, double, Var > > dbl; \
  lll_kernel_op< E, xdouble, Vect < E, xdouble >, \
                gensys_modules< E, xdouble, Var > > xdbl; \
  lll_kernel_fu< E, bigfloat, Vect < E, bigfloat >, \
                gensys_modules< E, bigfloat, Var > > bfl; \
}

#define ALG_POINTER(alg, SP, ex)                                                              \
{                                                                                           \
  alg.dbl.vector_operations().approx.set_pointer(SP.dbl); \
  alg.dbl.vector_operations().exact.set_pointer(SP.ex); \
  alg.xdbl.vector_operations().approx.set_pointer(SP.xdbl); \
  alg.xdbl.vector_operations().exact.set_pointer(SP.ex); \
  alg.bfl.vector_operations().approx.set_pointer(SP.bfl); \
  alg.bfl.vector_operations().exact.set_pointer(SP.ex); \
}

#define ALG_CALL(cl, alg, da, li, factor)     \
{                                             \
  if (factor < 2)                             \
    cl.dbl.alg(da, li);                       \
  else                                        \
    {                                         \
      if (factor == 2)                        \
        cl.xdbl.alg(da, li);                  \
      else                                    \
        cl.bfl.alg(da, li);                   \
    }                                         \
}

//
// Const for structure_info
//
const sdigit GRAM_MATRIX = 0x80000000;

//
// Const for lattice_info
//
const sdigit REDUCE_ROWS = 0x80000000;

// __dummy__ is never referenced and it is an illegal identifier anyway. (see C++ standard 17.4.3.1.2)
// static const bigfloat __dummy__ = 0.0;



#ifdef LIDIA_NAMESPACE
}	// end of namespace LiDIA
# undef IN_NAMESPACE_LIDIA
#endif



#endif	// LIDIA_LATTICE_DEFS_H_GUARD_
