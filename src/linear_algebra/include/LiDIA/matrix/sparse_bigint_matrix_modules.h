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


#ifndef LIDIA_SPARSE_BIGINT_MATRIX_MODULES_H_GUARD_
#define LIDIA_SPARSE_BIGINT_MATRIX_MODULES_H_GUARD_



#ifndef LIDIA_LIDIA_H_GUARD_
# include	"LiDIA/LiDIA.h"
#endif



#ifdef LIDIA_NAMESPACE
namespace LiDIA {
# define IN_NAMESPACE_LIDIA
#endif



class bigint;
class bigfloat;

template <>
class matrix< bigint >;



#define DELTA 2

#define column_oriented_sparse_matrix_modules COSMM



template< class T >
class column_oriented_sparse_matrix_modules
{
public:

	//
	// constructor
	//

	column_oriented_sparse_matrix_modules() {}

	//
	// destructor
	//

	~column_oriented_sparse_matrix_modules() {}


	//
	// member
	//

	const T & member(const MR< T > &A, lidia_size_t x, lidia_size_t y) const;

	//
	// sto
	//

	void sto(MR< T > &A, lidia_size_t y, lidia_size_t x, const T &e) const;

	//
	// swap rows
	//

	void swap_columns(MR< T > &A, lidia_size_t pos1, lidia_size_t pos2) const;

	//
	// max
	//

	void max_abs(const MR< T > &RES, T &MAX) const;

	//
	// subtract
	//

	void subtract_multiple_of_column(MR< T > &A,
					 lidia_size_t index1,
					 const T &q,
					 lidia_size_t index2,
					 lidia_size_t startr) const;

	void subtract_multiple_of_column_mod(MR< T > &A,
					     lidia_size_t index1,
					     const T &q,
					     lidia_size_t index2,
					     lidia_size_t startr,
					     const T &mod) const;

	void normalize_column_mod(MR< T > &A,
				  lidia_size_t index1,
				  const T &q,
				  lidia_size_t index2,
				  lidia_size_t startr,
				  const T &mod) const;

	//
	// negate
	//

	void negate_column(MR< T > &A,
			   lidia_size_t pos,
			   lidia_size_t len) const;

	void negate_column_mod(MR< T > &A,
			       lidia_size_t pos,
			       lidia_size_t len,
			       const T &mod) const;

	//
	// hadamard bound
	//

	void hadamard(const MR< T > &A, bigint &H) const;

	//
	// Init and update
	//

	void update_max_array(const MR< T > &A,
			      lidia_size_t i,
			      T *max_array) const;

	T *init_max_array(const MR< T > &A) const;

	//
	// Pivot search
	//

	lidia_size_t normal_pivot(const MR< T > &A,
				  lidia_size_t startr,
				  lidia_size_t startc,
				  lidia_size_t &index) const;

	lidia_size_t min_abs_of_row(const MR< T > &A,
				    lidia_size_t startr,
				    lidia_size_t startc,
				    lidia_size_t &index) const;

	lidia_size_t pivot_sorting_gcd(const MR< T > &A,
				       lidia_size_t startr,
				       lidia_size_t startc,
				       lidia_size_t &index) const;

	lidia_size_t minimal_norm(const MR< T > &A,
				  lidia_size_t startr,
				  lidia_size_t startc,
				  lidia_size_t &index,
				  lidia_size_t norm) const;


	lidia_size_t min_abs_of_row_plus_minimal_norm(const MR< T > &A,
						      lidia_size_t startr,
						      lidia_size_t startc,
						      lidia_size_t &index,
						      lidia_size_t norm) const;


	lidia_size_t minimal_norm_plus_min_abs_of_row(const MR< T > &A,
						      lidia_size_t startr,
						      lidia_size_t startc,
						      lidia_size_t &index,
						      lidia_size_t norm) const;

	lidia_size_t minimal_norm_plus_sorting_gcd(const MR< T > &A,
						   lidia_size_t startr,
						   lidia_size_t startc,
						   lidia_size_t &index,
						   lidia_size_t norm) const;

	lidia_size_t minimal_norm_plus_min_no_of_elements(const MR< T > &A,
							  lidia_size_t startr,
							  lidia_size_t startc,
							  lidia_size_t &index,
							  lidia_size_t norm) const;


	lidia_size_t sorting_gcd_plus_minimal_norm(const MR< T > &A,
						   lidia_size_t startr,
						   lidia_size_t startc,
						   lidia_size_t &index,
						   lidia_size_t norm) const;

	lidia_size_t min_no_of_elements_plus_minimal_norm(const MR< T > &A,
							  lidia_size_t startr,
							  lidia_size_t startc,
							  lidia_size_t &index,
							  lidia_size_t norm) const;

	//
	// Norm
	//

	bigint column_norm(const MR< T > &, lidia_size_t, lidia_size_t) const;

	void kennwerte(MR< T > &, T &, lidia_size_t &, T &) const;
	void max(MR< T > &, T &) const;

	//
	// Elimination
	//

	bool normalize_one_element_of_row(MR< T > &A,
					  lidia_size_t startr,
					  lidia_size_t startc,
					  lidia_size_t index,
					  lidia_size_t len) const;

	bool normalize_one_element_of_row(MR< T > &A,
					  matrix< bigint > &TR,
					  lidia_size_t startr,
					  lidia_size_t startc,
					  lidia_size_t index,
					  lidia_size_t len) const;

	bool normalize_one_element_of_row(MR< T > &A,
					  math_vector< bigfloat > &,
					  lidia_size_t startr,
					  lidia_size_t startc,
					  lidia_size_t index,
					  lidia_size_t len) const;

	bool mgcd_linear(MR< T > &A,
			 lidia_size_t startr,
			 lidia_size_t startc,
			 lidia_size_t len) const;
	bool mgcd_linear(MR< T > &A,
			 matrix< bigint > &TR,
			 lidia_size_t startr,
			 lidia_size_t startc,
			 lidia_size_t len) const;

	bool mgcd_bradley(MR< T > &A,
			  lidia_size_t startr,
			  lidia_size_t startc,
			  lidia_size_t len) const;
	bool mgcd_bradley(MR< T > &A,
			  matrix< bigint > &TR,
			  lidia_size_t startr,
			  lidia_size_t startc,
			  lidia_size_t len) const;

	bool mgcd_ilio(MR< T > &A,
		       lidia_size_t startr,
		       lidia_size_t startc,
		       lidia_size_t len) const;
	bool mgcd_ilio(MR< T > &A,
		       matrix< bigint > &TR,
		       lidia_size_t startr,
		       lidia_size_t startc,
		       lidia_size_t len) const;

	bool mgcd_opt(MR< T > &A,
		      lidia_size_t startr,
		      lidia_size_t startc,
		      lidia_size_t len) const;
	bool mgcd_opt(MR< T > &A,
		      matrix< bigint > &TR,
		      lidia_size_t startr,
		      lidia_size_t startc,
		      lidia_size_t len) const;


	bool xgcd_elimination(MR< T > &A,
			      lidia_size_t startr,
			      lidia_size_t startc,
			      lidia_size_t index,
			      lidia_size_t len) const;
	bool xgcd_elimination(MR< T > &A,
			      matrix< bigint > &,
			      lidia_size_t startr,
			      lidia_size_t startc,
			      lidia_size_t index,
			      lidia_size_t len) const;

	bool xgcd_elimination_mod(MR< T > &A,
				  lidia_size_t startr,
				  lidia_size_t startc,
				  lidia_size_t index,
				  lidia_size_t len,
				  const T &) const;

	bool mgcd_storjohann(MR< T > &A,
			     lidia_size_t startr,
			     lidia_size_t startc,
			     lidia_size_t len) const;
	bool mgcd_storjohann(MR< T > &A,
			     matrix< bigint > &TR,
			     lidia_size_t startr,
			     lidia_size_t startc,
			     lidia_size_t len) const;

	//
	// normalize
	//

	bool normalize_row(MR< T > &A,
			   lidia_size_t startr,
			   lidia_size_t startc,
			   lidia_size_t index,
			   lidia_size_t len) const;

	bool normalize_row_mod(MR< T > &A,
			       lidia_size_t startr,
			       lidia_size_t startc,
			       lidia_size_t index,
			       lidia_size_t len,
			       const T &) const;

	bool normalize_row(MR< T > &A,
			   matrix< bigint > &TR,
			   lidia_size_t startr,
			   lidia_size_t startc,
			   lidia_size_t index,
			   lidia_size_t len) const;

	bool normalize_row(MR< T > &A,
			   math_vector< bigfloat > &vec,
			   lidia_size_t startr,
			   lidia_size_t startc,
			   lidia_size_t index,
			   lidia_size_t len) const;

	bool normalize_row(MR< T > &A,
			   lidia_size_t startr,
			   lidia_size_t startc,
			   lidia_size_t index,
			   lidia_size_t len,
			   T *max_array,
			   const T &BOUND) const;

	bool normalize_row(MR< T > &A,
			   matrix< bigint > &TR,
			   lidia_size_t startr,
			   lidia_size_t startc,
			   lidia_size_t index,
			   lidia_size_t len,
			   T *max_array,
			   const T &BOUND) const;

	bool normalize_row(MR< T > &A,
			   math_vector< bigfloat > &vec,
			   lidia_size_t startr,
			   lidia_size_t startc,
			   lidia_size_t index,
			   lidia_size_t len,
			   T *max_array,
			   const T &BOUND) const;

};



#undef column_oriented_sparse_matrix_modules



#ifdef LIDIA_NAMESPACE
}	// end of namespace LiDIA
# undef IN_NAMESPACE_LIDIA
#endif



#endif	// LIDIA_SPARSE_BIGINT_MATRIX_MODULES_H_GUARD_
