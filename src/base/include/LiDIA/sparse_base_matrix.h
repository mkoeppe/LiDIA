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


#ifndef LIDIA_SPARSE_BASE_MATRIX_H_GUARD_
#define LIDIA_SPARSE_BASE_MATRIX_H_GUARD_



#ifndef LIDIA_BASE_VECTOR_H_GUARD_
# include	"LiDIA/base_vector.h"
#endif
#ifndef LIDIA_MATRIX_REPRESENTATION_H_GUARD_
# include	"LiDIA/matrix/matrix_representation.h"
#endif
#ifndef LIDIA_SPARSE_BASE_MATRIX_KERNEL_H_GUARD_
# include	"LiDIA/matrix/sparse_base_matrix_kernel.h"
#endif
#ifndef LIDIA_BASE_MATRIX_ALGORITHMS_H_GUARD_
# include	"LiDIA/matrix/base_matrix_algorithms.h"
#endif



#ifdef LIDIA_NAMESPACE
namespace LiDIA {
# define IN_NAMESPACE_LIDIA
#endif



template< class T >
class sparse_base_matrix : public MR< T >
{
	//
	// friend classes
	//

	friend class base_matrix< T >;

	//
	// modul definitions
	//

public:

	const SBMK< T > S_base_modul;

	const BMA< T, SBMK < T >, SBMK< T > > SS_base_modul;

	//
	// constructors
	//

public:

	sparse_base_matrix();
	sparse_base_matrix(lidia_size_t, lidia_size_t);
	sparse_base_matrix(const base_vector< T > &);
	sparse_base_matrix(const sparse_base_matrix< T > &);
	sparse_base_matrix(lidia_size_t, lidia_size_t, const T **);

	//
	// destructor
	//

public:

	~sparse_base_matrix();

	//
	// Input / Output
	//

	std::ostream & write(std::ostream &) const;
	std::istream & read(std::istream &);


	/////////////////////////////
	// BEGIN: access functions //
	/////////////////////////////

	//
	// set zero_element / get zero_element
	//

public:

	void set_zero_element(const T &a)
	{
		this->Zero = a;
	}

	const T & get_zero_element() const
	{
		return this->Zero;
	}

	//
	// set print_mode / get print_mode
	//

public:

	void set_print_mode(unsigned long art)
	{
		this->bitfield.set_print_mode(art);
	}

	unsigned long get_print_mode()
	{
		return this->bitfield.get_print_mode();
	}

	//
	// set storage_mode / get storage_mode
	//

public:

	void set_storage_mode(unsigned long art)
	{
		this->bitfield.set_storage_mode(art);
	}

	unsigned long get_storage_mode()
	{
		return this->bitfield.get_storage_mode();
	}

	// set orientation / get orientation
public:

	void set_orientation(unsigned long art)
	{
		S_base_modul.change_orientation(*this, art);
		this->bitfield.set_orientation(art);
	}

	unsigned long get_orientation()
	{
		return this->bitfield.get_orientation();
	}

	//
	// set structure_mode / get structure_mode
	//

public:

	void set_structure_mode(unsigned long art)
	{
		this->bitfield.set_structure_mode(art);
	}

	unsigned long get_structure_mode()
	{
		return this->bitfield.get_structure_mode();
	}

	//
	// set info_mode / get info_mode
	//

public:

	void set_info_mode(unsigned long art)
	{
		this->bitfield.set_info_mode(art);
	}

	unsigned long get_info_mode()
	{
		return this->bitfield.get_info_mode();
	}

	//
	// set lattice_mode / get lattice_mode
	//

public:

	void set_lattice_mode(unsigned long art)
	{
		this->bitfield.set_lattice_mode(art);
	}

	unsigned long get_lattice_mode()
	{
		return this->bitfield.get_lattice_mode();
	}

	//
	// element access
	//

public:

	void sto(lidia_size_t, lidia_size_t, const T &);

	const T & operator() (lidia_size_t x, lidia_size_t y) const
	{
		return member(x, y);
	}

	const T & member(lidia_size_t, lidia_size_t) const;

	//
	// column access
	//

public:

	base_vector< T > operator() (lidia_size_t i) const
	{
		base_vector< T > RES;
		get_column_vector(RES, i);
		return RES;
	}

	void sto_column(const T *, lidia_size_t, lidia_size_t, lidia_size_t from = 0);
	void sto_column_vector(const base_vector< T > &, lidia_size_t, lidia_size_t,
			       lidia_size_t from = 0);

	T * get_column(lidia_size_t i) const
	{
		register T *RES = new T[this->rows];
		memory_handler(RES, DMESSAGE, "get_column(lidia_size_t) :: "
			       "Error in memory allocation (RES)");
		get_column(RES, i);
		return RES;
	}

	base_vector< T > get_column_vector(lidia_size_t i) const
	{
		base_vector< T > RES;
		get_column_vector(RES, i);
		return RES;
	}

	void get_column(T *, lidia_size_t) const;
	void get_column_vector(base_vector< T > &, lidia_size_t) const;

	//
	// row access
	//

public:

	base_vector< T > operator[] (lidia_size_t i) const
	{
		base_vector< T > RES;
		get_row_vector(RES, i);
		return RES;
	}

	void sto_row(const T *, lidia_size_t, lidia_size_t, lidia_size_t from = 0);
	void sto_row_vector(const base_vector< T > &, lidia_size_t, lidia_size_t,
			    lidia_size_t from = 0);

	T *get_row(lidia_size_t i) const
	{
		register T *RES = new T[this->columns];
		memory_handler(RES, DMESSAGE, "get_column(lidia_size_t) :: "
			       "Error in memory allocation (RES)");
		get_row(RES, i);
		return RES;
	}

	base_vector< T > get_row_vector(lidia_size_t i) const
	{
		base_vector< T > RES;
		get_row_vector(RES, i);
		return RES;
	}

	void get_row(T *, lidia_size_t) const;
	void get_row_vector(base_vector< T > &, lidia_size_t) const;

	//
	// array access
	//

	//protected:
	//Remove public again, if the lattice code does not
	//need this function anymore.
public:

	T ** get_data_address() const
	{
		return this->value;
	}

public:

	T ** get_data() const;
	void set_data(const T **, lidia_size_t, lidia_size_t);

	///////////////////////////
	// END: access functions //
	///////////////////////////

	//
	// insert_at
	//

public:

	void insert_at(lidia_size_t, lidia_size_t,
		       const sparse_base_matrix< T > &, lidia_size_t,
		       lidia_size_t, lidia_size_t, lidia_size_t);

	//
	// insert_columns, insert_rows, remove_columns, remove_rows
	//

public:

	void insert_columns(lidia_size_t *, const T **);
	void remove_columns(lidia_size_t *);

	void insert_rows(lidia_size_t *, const T **);
	void remove_rows(lidia_size_t *);

	//
	// split functions
	//

public:

	void split_t(sparse_base_matrix< T > &, sparse_base_matrix< T > &,
		     sparse_base_matrix< T > &, sparse_base_matrix< T > &) const;

	void split_h(sparse_base_matrix< T > &, sparse_base_matrix< T > &) const;
	void split_h(T *, sparse_base_matrix< T > &) const;
	void split_h(sparse_base_matrix< T > &, T *) const;
	void split_h(base_vector< T > &, sparse_base_matrix< T > &) const;
	void split_h(sparse_base_matrix< T > &, base_vector< T > &) const;

	void split_v(sparse_base_matrix< T > &, sparse_base_matrix< T > &) const;
	void split_v(T *, sparse_base_matrix< T > &) const;
	void split_v(sparse_base_matrix< T > &, T *) const;
	void split_v(base_vector< T > &, sparse_base_matrix< T > &) const;
	void split_v(sparse_base_matrix< T > &, base_vector< T > &) const;

	//
	// compose functions
	//

public:

	void compose_t(const sparse_base_matrix< T > &, const sparse_base_matrix< T > &,
		       const sparse_base_matrix< T > &, const sparse_base_matrix< T > &);

	void compose_h(const sparse_base_matrix< T > &, const sparse_base_matrix< T > &);
	void compose_h(const T *, const sparse_base_matrix< T > &);
	void compose_h(const sparse_base_matrix< T > &, const T *);
	void compose_h(const base_vector< T > &, const sparse_base_matrix< T > &);
	void compose_h(const sparse_base_matrix< T > &, const base_vector< T > &);

	void compose_v(const sparse_base_matrix< T > &, const sparse_base_matrix< T > &);
	void compose_v(const T *, const sparse_base_matrix< T > &);
	void compose_v(const sparse_base_matrix< T > &, const T *);
	void compose_v(const base_vector< T > &, const sparse_base_matrix< T > &);
	void compose_v(const sparse_base_matrix< T > &, const base_vector< T > &);

	//
	// exchange functions / swap functions
	//

public:

	void swap(sparse_base_matrix< T > &);

	void swap_columns(lidia_size_t, lidia_size_t);
	void swap_rows(lidia_size_t, lidia_size_t);


	//
	// structur functions
	//

public:

	lidia_size_t get_no_of_rows() const
	{
		return this->rows;
	}

	lidia_size_t get_no_of_columns() const
	{
		return this->columns;
	}

	void set_no_of_rows(lidia_size_t);
	void set_no_of_columns(lidia_size_t);

	void resize(lidia_size_t, lidia_size_t);

	void reset()
	{
		this->kill();
	}

	void kill();

	//
	// assignment
	//

public:

	sparse_base_matrix< T > & operator = (const sparse_base_matrix< T > &B)
	{
		this->assign(B);
		return *this;
	}

	void assign(const sparse_base_matrix< T > &);

	//
	// diagonal function
	//

public:

	void diag(const T &, const T &);

	//
	// transpose function
	//

public:

	void trans(const sparse_base_matrix< T > &);

	sparse_base_matrix< T > trans() const;

	//
	// stream handling
	//

public:

	void write_to_beauty(std::ostream &) const;

	void write_to_lidia(std::ostream &out) const
	{
		write_to_stream(out);
	}

	void read_from_lidia(std::istream &in)
	{
		read_from_stream(in);
	}

	void write_to_stream(std::ostream &) const;
	void read_from_stream(std::istream &);

	void write_to_mathematica(std::ostream &) const;
	void read_from_mathematica(std::istream &);

	void write_to_maple(std::ostream &) const;
	void read_from_maple(std::istream &);

	void write_to_gp(std::ostream &) const;
	void read_from_gp(std::istream &);

	void write_to_kash(std::ostream &) const;
	void read_from_kash(std::istream &);

	void write_to_latex(std::ostream &) const;

	void write_to_magma(std::ostream &) const;

	//
	// boolean functions
	//

public:

	bool is_column_zero(lidia_size_t) const;

	bool is_row_zero(lidia_size_t) const;

	bool is_matrix_zero() const;

	//
	// change orientation
	//

	void change_orientation(unsigned long);

	//
	// status report
	//

	void status_report();

};

//
// Input / Output
//

template< class T >
inline std::ostream &
operator << (std::ostream &out, const sparse_base_matrix< T > &M)
{
	return M.write(out);
}



template< class T >
inline std::istream &
operator >> (std::istream &in, sparse_base_matrix< T > &M)
{
	return M.read(in);
}



//
// exchange functions / swap functions
//

template< class T >
inline void
swap(sparse_base_matrix< T > &A, sparse_base_matrix< T > &B)
{
	A.swap(B);
}



//
// assignment
//

template< class T >
inline void
assign(sparse_base_matrix< T > &RES, const sparse_base_matrix< T > &M)
{
	RES.assign(M);
}



//
// diagonal function
//

template< class T >
inline void
diag(sparse_base_matrix< T > &A, const T &a, const T &b)
{
	A.diag(a, b);
}



//
// transpose function
//

template< class T >
inline void
trans(sparse_base_matrix< T > &AT, const sparse_base_matrix< T > &A)
{
	AT.trans(A);
}



template< class T >
inline sparse_base_matrix< T >
trans(const sparse_base_matrix< T > &A)
{
	return A.trans();
}



#ifdef LIDIA_NAMESPACE
}	// end of namespace LiDIA
# undef IN_NAMESPACE_LIDIA
#endif



#ifdef LIDIA_INCLUDE_CC
# include	"LiDIA/sparse_base_matrix.cc"
#endif



#endif	// LIDIA_SPARSE_BASE_MATRIX_H_GUARD_
