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


#ifndef LIDIA_DENSE_BIGINT_MATRIX_MODULES_H_GUARD_
#define LIDIA_DENSE_BIGINT_MATRIX_MODULES_H_GUARD_



#ifndef LIDIA_DENSE_RING_MATRIX_KERNEL_H_GUARD_
# include	"LiDIA/matrix/dense_ring_matrix_kernel.h"
#endif



#ifdef LIDIA_NAMESPACE
namespace LiDIA {
# define IN_NAMESPACE_LIDIA
#endif



#define row_oriented_dense_matrix_modules RODMM



template< class T >
class row_oriented_dense_matrix_modules : public DRMK< T >
{
public:

	//
	// constructor
	//

	row_oriented_dense_matrix_modules() {}

	//
	// destructor
	//

	~row_oriented_dense_matrix_modules() {}

	//
	// max
	//

	void max_abs(const MR< T > &RES, T &MAX) const;

	//
	// special functions
	//

	void subtract_multiple_of_column(MR< T > &A,
					 lidia_size_t l1,
					 const T &q,
					 lidia_size_t l2,
					 lidia_size_t l) const;

	void subtract_multiple_of_column_mod(MR< T > &A,
					     lidia_size_t l1,
					     const T &q,
					     lidia_size_t l2,
					     lidia_size_t l,
					     const T &mod) const;

	void normalize_column_mod(MR< T > &A,
				  lidia_size_t l1,
				  const T &q,
				  lidia_size_t l2,
				  lidia_size_t l,
				  const T &mod) const;

	//
	// negate
	//

	void negate_column_mod(MR< T > &A,
			       lidia_size_t index,
			       lidia_size_t l,
			       const T&mod) const;

	void negate_column(MR< T > &A,
			   lidia_size_t index,
			   lidia_size_t l) const;

	//
	// init and update
	//

	T * init_max_array(const MR< T > &A) const;

	void update_max_array(const MR< T > &A,
			      lidia_size_t i,
			      T *max_array) const;

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

	lidia_size_t minimal_norm(MR< T > &A,
				  lidia_size_t startr,
				  lidia_size_t startc,
				  lidia_size_t &index,
				  lidia_size_t norm) const;


	lidia_size_t min_abs_of_row_plus_minimal_norm(MR< T > &A,
						      lidia_size_t startr,
						      lidia_size_t startc,
						      lidia_size_t &index,
						      lidia_size_t norm) const;


	lidia_size_t minimal_norm_plus_min_abs_of_row(MR< T > &A,
						      lidia_size_t startr,
						      lidia_size_t startc,
						      lidia_size_t &index,
						      lidia_size_t norm) const;

	lidia_size_t minimal_norm_plus_sorting_gcd(MR< T > &A,
						   lidia_size_t startr,
						   lidia_size_t startc,
						   lidia_size_t &index,
						   lidia_size_t norm) const;

	lidia_size_t minimal_norm_plus_min_no_of_elements(MR< T > &A,
							  lidia_size_t startr,
							  lidia_size_t startc,
							  lidia_size_t &index,
							  lidia_size_t norm) const;

	lidia_size_t sorting_gcd_plus_minimal_norm(MR< T > &A,
						   lidia_size_t startr,
						   lidia_size_t startc,
						   lidia_size_t &index,
						   lidia_size_t norm) const;

	lidia_size_t min_no_of_elements_plus_minimal_norm(MR< T > &A,
							  lidia_size_t startr,
							  lidia_size_t startc,
							  lidia_size_t &index,
							  lidia_size_t norm) const;

	//
	// Norm
	//

	bigint column_norm(MR< T > &,
			   lidia_size_t,
			   lidia_size_t) const;

	void kennwerte(MR< T > &,
		       T &,
		       lidia_size_t &,
		       T &) const;
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
			      matrix< bigint > &TR,
			      lidia_size_t startr,
			      lidia_size_t startc,
			      lidia_size_t index,
			      lidia_size_t len) const;

	bool xgcd_elimination_mod(MR< T > &A,
				  lidia_size_t startr,
				  lidia_size_t startc,
				  lidia_size_t index,
				  lidia_size_t len,
				  const T &DET) const;


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
			       const T &DET) const;

	bool normalize_row(MR< T > &A,
			   matrix< bigint > &TR,
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


	//
	// mgcd
	//


};



#undef row_oriented_dense_matrix_modules
#define column_oriented_dense_matrix_modules CODMM


template< class T >
class column_oriented_dense_matrix_modules : public DRMK< T >
{
public:

	//
	// constructor
	//

	column_oriented_dense_matrix_modules() {}

	//
	// destructor
	//

	~column_oriented_dense_matrix_modules() {}

	//
	// special functions
	//

	const T & member(const MR< T > &A,
			 lidia_size_t x,
			 lidia_size_t y) const;
	void swap_columns(const MR< T > &A,
			  lidia_size_t l1,
			  lidia_size_t l2) const;

	void subtract_multiple_of_column(MR< T > &A,
					 lidia_size_t l1,
					 const T &q,
					 lidia_size_t l2,
					 lidia_size_t l) const;

	void subtract_multiple_of_column_mod(MR< T > &A,
					     lidia_size_t l1,
					     const T &q,
					     lidia_size_t l2,
					     lidia_size_t l,
					     const T &mod) const;

	void normalize_column_mod(MR< T > &A,
				  lidia_size_t l1,
				  const T &q,
				  lidia_size_t l2,
				  lidia_size_t l,
				  const T &mod) const;

	//
	// negate
	//

	void negate_column(MR< T > &A,
			   lidia_size_t index,
			   lidia_size_t l) const;

	void negate_column_mod(MR< T > &A,
			       lidia_size_t index,
			       lidia_size_t l,
			       const T &mod) const;

	//
	// Init and update
	//

	T * init_max_array(const MR< T > &A) const;

	void update_max_array(const MR< T > &A,
			      lidia_size_t i,
			      T *max_array) const;

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

	lidia_size_t sorting_gcd_plus_minimal_norm(const MR< T > &A,
						   lidia_size_t startr,
						   lidia_size_t startc,
						   lidia_size_t &index,
						   lidia_size_t norm) const;

	//
	// Norm
	//

	bigint column_norm(const MR< T > &,
			   lidia_size_t,
			   lidia_size_t) const;

	void kennwerte(MR< T > &,
		       T &,
		       lidia_size_t &,
		       T &) const;
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

	bool xgcd_elimination_mod(MR< T > &A,
				  lidia_size_t startr,
				  lidia_size_t startc,
				  lidia_size_t index,
				  lidia_size_t len,
				  const T &mod) const;


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
			       const T &mod) const;

	bool normalize_row(MR< T > &A,
			   matrix< bigint > &TR,
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

};



#undef column_oriented_dense_matrix_modules



#ifdef LIDIA_NAMESPACE
}	// end of namespace LiDIA
# undef IN_NAMESPACE_LIDIA
#endif



#endif	// LIDIA_DENSE_BIGINT_MATRIX_MODULES_H_GUARD_
