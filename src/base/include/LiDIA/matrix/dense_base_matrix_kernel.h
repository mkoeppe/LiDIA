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


#ifndef LIDIA_DENSE_BASE_MATRIX_KERNEL_H_GUARD_
#define LIDIA_DENSE_BASE_MATRIX_KERNEL_H_GUARD_



#ifndef LIDIA_MATH_VECTOR_H_GUARD_
# include	"LiDIA/math_vector.h"
#endif
#ifndef LIDIA_MATRIX_REPRESENTATION_H_GUARD_
# include	"LiDIA/matrix/matrix_representation.h"
#endif



#ifdef LIDIA_NAMESPACE
namespace LiDIA {
# define IN_NAMESPACE_LIDIA
#endif



#define dense_base_matrix_kernel DBMK



template< class T >
class dense_base_matrix_kernel
{

public:

	//
	// constructor
	//

	dense_base_matrix_kernel () {}

	//
	// destructor
	//

	~dense_base_matrix_kernel () {}

	//
	// constructor kernel
	//

	void constructor(MR< T > &, lidia_size_t, lidia_size_t) const;
	void constructor(MR< T > &, lidia_size_t, lidia_size_t, const T **) const;
	void constructor(MR< T > &, const MR< T > &) const;

	//
	// destructor kernel
	//

	void destructor(MR< T > &) const;

	//
	// element access kernel
	//

	void sto(MR< T > &, lidia_size_t, lidia_size_t, const T &) const;

	T & member(const MR< T > &, lidia_size_t, lidia_size_t) const;

	//
	// column access kernel
	//

	void sto_column(MR< T > &, const T *, lidia_size_t,
			lidia_size_t, lidia_size_t) const;

	void get_column(const MR< T > &, T *, lidia_size_t) const;

	//
	// row access kernel
	//

	void sto_row(MR< T > &, const T *, lidia_size_t,
		     lidia_size_t, lidia_size_t) const;

	void get_row(const MR< T > &, T *, lidia_size_t) const;

	//
	// insert and remove kernel
	//

	void insert_columns(MR< T > &, lidia_size_t *, const T **) const;
	void insert_column_at(MR< T > &, lidia_size_t, const T *, lidia_size_t) const;
	void remove_columns(MR< T > &, lidia_size_t *) const;

	void insert_rows(MR< T > &, lidia_size_t *, const T **) const;
	void insert_row_at(MR< T > &, lidia_size_t, const T *, lidia_size_t) const;
	void remove_rows(MR< T > &, lidia_size_t *) const;

	//
	// exchange functions / swap functions
	//

	void swap_columns(MR< T > &, lidia_size_t, lidia_size_t) const;

	void swap_rows(MR< T > &, lidia_size_t, lidia_size_t) const;

	//
	// assignment kernel
	//

	void assign(MR< T > &, const MR< T > &) const;

	//
	// diagonal function
	//

	void diag(MR< T > &, const T &, const T &) const;

	//
	// transpose function
	//

	void trans(MR< T > &, const MR< T > &) const;

	//
	// stream handling
	//

	void write_to_beauty(const MR< T > &, std::ostream &) const;

	void write_to_stream(const MR< T > &, std::ostream &) const;
	void read_from_stream(MR< T > &, std::istream &) const;

	void write_to_mathematica(const MR< T > &, std::ostream &) const;
	void read_from_mathematica(MR< T > &, std::istream &) const;

	void write_to_maple(const MR< T > &, std::ostream &) const;
	void read_from_maple(MR< T > &, std::istream &) const;

	void write_to_gp(const MR< T > &, std::ostream &) const;
	void read_from_gp(MR< T > &, std::istream &) const;

	void write_to_kash(const MR< T > &, std::ostream &) const;
	void read_from_kash(MR< T > &, std::istream &) const;

	void write_to_latex(const MR< T > &, std::ostream &) const;

	void write_to_magma(const MR< T > &, std::ostream &) const;

	//
	// structure functions kernel
	//

	void set_no_of_rows(MR< T > &, lidia_size_t) const;
	void set_no_of_columns(MR< T > &, lidia_size_t) const;

	void kill(MR< T > &) const;

	//
	// boolean functions
	//

	bool is_column_zero(const MR< T > &, lidia_size_t) const;

	bool is_row_zero(const MR< T > &, lidia_size_t) const;

	bool is_matrix_zero(const MR< T > &) const;

	//
	// change orientation
	//

	void change_orientation(MR< T > &, unsigned long) const;
};



#ifdef LIDIA_NAMESPACE
}	// end of namespace LiDIA
# undef IN_NAMESPACE_LIDIA
#endif



#ifdef LIDIA_INCLUDE_CC
# include	"LiDIA/matrix/dense_base_matrix_kernel.cc"
#endif



#undef dense_base_matrix_kernel



#endif	// LIDIA_DENSE_BASE_MATRIX_KERNEL_H_GUARD_
