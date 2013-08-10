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


#ifndef LIDIA_BASE_MATRIX_H_GUARD_
#define LIDIA_BASE_MATRIX_H_GUARD_



#ifndef LIDIA_LIDIA_H_GUARD_
# include	"LiDIA/LiDIA.h"
#endif
#ifndef LIDIA_BASE_VECTOR_H_GUARD_
# include	"LiDIA/base_vector.h"
#endif
#ifndef LIDIA_DENSE_BASE_MATRIX_H_GUARD_
# include	"LiDIA/dense_base_matrix.h"
#endif
#ifndef LIDIA_SPARSE_BASE_MATRIX_H_GUARD_
# include	"LiDIA/sparse_base_matrix.h"
#endif
#ifndef LIDIA_MATRIX_REPRESENTATION_H_GUARD_
# include	"LiDIA/matrix/matrix_representation.h"
#endif
#ifndef LIDIA_DENSE_BASE_MATRIX_KERNEL_H_GUARD_
# include	"LiDIA/matrix/dense_base_matrix_kernel.h"
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



#define DBMKex DBMK< T >
#define SBMKex SBMK< T >



template< class T >
class base_matrix : public MR< T >
{
	//
	// modul definitions
	//

protected:

	const DBMKex D_base_modul;
	const SBMKex S_base_modul;

	const BMA< T, DBMKex, DBMKex > DD_base_modul;
	const BMA< T, DBMKex, SBMKex > DS_base_modul;
	const BMA< T, SBMKex, DBMKex > SD_base_modul;
	const BMA< T, SBMKex, SBMKex > SS_base_modul;

	//
	// constructors
	//

public:

	base_matrix();
	base_matrix(const matrix_flags &);
	base_matrix(lidia_size_t, lidia_size_t);
	base_matrix(lidia_size_t, lidia_size_t, const matrix_flags &);
	base_matrix(const base_vector< T > &);
	base_matrix(const base_vector< T > &, const matrix_flags &);
	base_matrix(const base_matrix< T > &);
	base_matrix(const base_matrix< T > &, const matrix_flags &);
	base_matrix(const dense_base_matrix< T > &);
	base_matrix(const sparse_base_matrix< T > &);
	base_matrix(lidia_size_t, lidia_size_t, const T **);
	base_matrix(lidia_size_t, lidia_size_t, const T **, const matrix_flags &);

	//
	// destructor
	//

public:

	~base_matrix();

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
		set_representation(art & matrix_flags::representation);
		set_orientation(art & matrix_flags::orientation);
	}

	unsigned long get_storage_mode()
	{
		return this->bitfield.get_storage_mode();
	}

	//
	// set representation / get representation
	//

public:

	void set_representation(unsigned long art)
	{
		this->change_representation(art);
		this->bitfield.set_representation(art);
	}

	unsigned long get_representation()
	{
		return this->bitfield.get_representation();
	}

	// set orientation / get orientation
public:

	void set_orientation(unsigned long art)
	{
		if (this->bitfield.get_representation() == matrix_flags::dense_representation)
			this->D_base_modul.change_orientation(*this, art);
		else
			this->S_base_modul.change_orientation(*this, art);

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
		memory_handler(RES, DMESSAGE, "get_row(lidia_size_t) :: "
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
		       const base_matrix< T > &, lidia_size_t,
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

	void split_t(base_matrix< T > &, base_matrix< T > &,
		     base_matrix< T > &, base_matrix< T > &) const;

	void split_h(base_matrix< T > &, base_matrix< T > &) const;
	void split_h(T *, base_matrix< T > &) const;
	void split_h(base_matrix< T > &, T *) const;
	void split_h(base_vector< T > &, base_matrix< T > &) const;
	void split_h(base_matrix< T > &, base_vector< T > &) const;

	void split_v(base_matrix< T > &, base_matrix< T > &) const;
	void split_v(T *, base_matrix< T > &) const;
	void split_v(base_matrix< T > &, T *) const;
	void split_v(base_vector< T > &, base_matrix< T > &) const;
	void split_v(base_matrix< T > &, base_vector< T > &) const;

	//
	// compose functions
	//

public:

	void compose_t(const base_matrix< T > &, const base_matrix< T > &,
		       const base_matrix< T > &, const base_matrix< T > &);

	void compose_h(const base_matrix< T > &, const base_matrix< T > &);
	void compose_h(const T *, const base_matrix< T > &);
	void compose_h(const base_matrix< T > &, const T *);
	void compose_h(const base_vector< T > &, const base_matrix< T > &);
	void compose_h(const base_matrix< T > &, const base_vector< T > &);

	void compose_v(const base_matrix< T > &, const base_matrix< T > &);
	void compose_v(const T *, const base_matrix< T > &);
	void compose_v(const base_matrix< T > &, const T *);
	void compose_v(const base_vector< T > &, const base_matrix< T > &);
	void compose_v(const base_matrix< T > &, const base_vector< T > &);

	//
	// exchange functions / swap functions
	//

public:

	void swap(base_matrix< T > &);

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

	base_matrix< T > & operator = (const base_matrix< T > &B)
	{
		this->assign(B);
		return *this;
	}

	void assign(const base_matrix< T > &);

	//
	// diagonal function
	//

public:

	void diag(const T &, const T &);

	//
	// transpose function
	//

public:

	void trans(const base_matrix< T > &);
	base_matrix< T > trans() const;

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
	// internal functions
	//

public:

	void status_report();

private:

	void change_representation(unsigned long);
};



//
// Input / Output
//

template< class T >
inline std::ostream &
operator << (std::ostream & out, const base_matrix< T > & M)
{
	return M.write(out);
}



template< class T >
inline std::istream &
operator >> (std::istream & in, base_matrix< T > & M)
{
	return M.read(in);
}



//
// exchange functions / swap functions
//

template< class T >
inline void
swap(base_matrix< T > &A, base_matrix< T > &B)
{
	A.swap(B);
}



//
// assignment
//

template< class T >
inline void
assign(base_matrix< T > &RES, const base_matrix< T > &M)
{
	RES.assign(M);
}



//
// diagonal function
//

template< class T >
inline void
diag(base_matrix< T > &A, const T &a, const T &b)
{
	A.diag(a, b);
}



//
// transpose function
//

template< class T >
inline void
trans(base_matrix< T > &AT, const base_matrix< T > &A)
{
	AT.trans(A);
}



template< class T >
inline base_matrix< T >
trans(const base_matrix< T > &A)
{
	return A.trans();
}



//
// specialization for type bigmod
//

template <>
class base_matrix< bigmod >
{
	//
	// data structure
	//

private:

	MR< bigint > bigint_value;
	bigint bigint_modulus;

	MR< long > long_value;
	long long_modulus;

	bool BigintLong;

	//
	// modul definitions
	//

protected:

	const DBMK< bigint > D_base_modul_bigint;
	const SBMK< bigint > S_base_modul_bigint;

	const BMA< bigint, DBMK < bigint >, DBMK< bigint > > DD_base_modul_bigint;
	const BMA< bigint, DBMK < bigint >, SBMK< bigint > > DS_base_modul_bigint;
	const BMA< bigint, SBMK < bigint >, DBMK< bigint > > SD_base_modul_bigint;
	const BMA< bigint, SBMK < bigint >, SBMK< bigint > > SS_base_modul_bigint;

	const DBMK< long > D_base_modul_long;
	const SBMK< long > S_base_modul_long;

	const BMA< long, DBMK < long >, DBMK< long > > DD_base_modul_long;
	const BMA< long, DBMK < long >, SBMK< long > > DS_base_modul_long;
	const BMA< long, SBMK < long >, DBMK< long > > SD_base_modul_long;
	const BMA< long, SBMK < long >, SBMK< long > > SS_base_modul_long;

	//
	// constructors
	//

public:

	base_matrix< bigmod > ();
	base_matrix< bigmod > (const bigint &);
	base_matrix< bigmod > (const bigint &, const matrix_flags &);
	base_matrix< bigmod > (lidia_size_t, lidia_size_t, const bigint &);
	base_matrix< bigmod > (lidia_size_t, lidia_size_t, const bigint &,
			       const matrix_flags &);
	base_matrix< bigmod > (const base_matrix< bigmod > &);
	base_matrix< bigmod > (const base_matrix< bigmod > &,
			       const matrix_flags &);

	base_matrix< bigmod > (const base_matrix< bigint > &, const bigint &);
	base_matrix< bigmod > (const base_matrix< long > &, const long &);

	//
	// destructor
	//

public:

	~base_matrix();

	//
	// Input / Output
	//

public:

	std::ostream & write(std::ostream &) const;
	std::istream & read(std::istream &);


	/////////////////////////////
	// BEGIN: access functions //
	/////////////////////////////

	//
	// set zero_element / get zero_element
	//

public:

	void set_zero_element(const bigmod &a)
	{
		if (BigintLong)
			bigint_value.Zero = a.mantissa();
		else {
			bigint TMP = a.mantissa();
			TMP.longify(long_value.Zero);
		}
	}

	bigmod get_zero_element() const
	{
		bigmod TMP;
		if (BigintLong) {
			TMP.set_modulus(bigint_modulus);
			TMP = bigint_value.Zero;
		}
		else {
			TMP.set_modulus(long_modulus);
			TMP = long_value.Zero;
		}
		return TMP;
	}

	//
	// set print_mode / get print_mode
	//

public:

	void set_print_mode(unsigned long art)
	{
		if (BigintLong)
			bigint_value.bitfield.set_print_mode(art);
		else
			long_value.bitfield.set_print_mode(art);
	}

	unsigned long get_print_mode()
	{
		if (BigintLong)
			return bigint_value.bitfield.get_print_mode();
		else
			return long_value.bitfield.get_print_mode();
	}

	// set representation / get representation
public:

	void set_representation(unsigned long art)
	{
		change_representation(art);
		if (BigintLong)
			bigint_value.bitfield.set_representation(art);
		else
			long_value.bitfield.set_representation(art);
	}

	unsigned long get_representation()
	{
		if (BigintLong)
			return bigint_value.bitfield.get_representation();
		else
			return long_value.bitfield.get_representation();
	}

	// set orientation / get orientation
public:

	void set_orientation(unsigned long art)
	{
		if (BigintLong) {
			if (bigint_value.bitfield.get_representation() ==
			    matrix_flags::dense_representation)
				D_base_modul_bigint.change_orientation(bigint_value, art);
			else
				S_base_modul_bigint.change_orientation(bigint_value, art);

			bigint_value.bitfield.set_orientation(art);
		}
		else {
			if (long_value.bitfield.get_representation() ==
			    matrix_flags::dense_representation)
				D_base_modul_long.change_orientation(long_value, art);
			else
				S_base_modul_long.change_orientation(long_value, art);

			long_value.bitfield.set_orientation(art);
		}
	}

	unsigned long get_orientation()
	{
		if (BigintLong)
			return bigint_value.bitfield.get_orientation();
		else
			return long_value.bitfield.get_orientation();
	}

	//
	// set storage_mode / get storage_mode
	//

public:

	void set_storage_mode(unsigned long art)
	{
		set_representation(art & matrix_flags::representation);
		set_orientation(art & matrix_flags::orientation);
	}

	unsigned long get_storage_mode()
	{
		if (BigintLong)
			return bigint_value.bitfield.get_storage_mode();
		else
			return long_value.bitfield.get_storage_mode();
	}

	//
	// set structure_mode / get structure_mode
	//

public:

	void set_structure_mode(unsigned long art)
	{
		if (BigintLong)
			bigint_value.bitfield.set_structure_mode(art);
		else
			long_value.bitfield.set_structure_mode(art);
	}

	unsigned long get_structure_mode()
	{
		if (BigintLong)
			return bigint_value.bitfield.get_structure_mode();
		else
			return long_value.bitfield.get_structure_mode();
	}

	//
	// set info_mode / get info_mode
	//

public:

	void set_info_mode(unsigned long art)
	{
		if (BigintLong)
			bigint_value.bitfield.set_info_mode(art);
		else
			long_value.bitfield.set_info_mode(art);
	}

	unsigned long get_info_mode()
	{
		if (BigintLong)
			return bigint_value.bitfield.get_info_mode();
		else
			return long_value.bitfield.get_info_mode();
	}

	//
	// set lattice_mode / get lattice_mode
	//

public:

	void set_lattice_mode(unsigned long art)
	{
		if (BigintLong)
			bigint_value.bitfield.set_lattice_mode(art);
		else
			long_value.bitfield.set_lattice_mode(art);
	}

	unsigned long get_lattice_mode()
	{
		if (BigintLong)
			return bigint_value.bitfield.get_lattice_mode();
		else
			return long_value.bitfield.get_lattice_mode();
	}

	//
	// element access
	//

public:

	void sto(lidia_size_t, lidia_size_t, const bigmod &);

	const bigmod operator() (lidia_size_t x, lidia_size_t y) const
	{
		return member(x, y);
	}

	const bigmod member(lidia_size_t, lidia_size_t) const;

	//
	// column access
	//

public:

	base_vector< bigmod > operator() (lidia_size_t i) const
	{
		base_vector< bigmod > RES;
		get_column_vector(RES, i);
		return RES;
	}

	void sto_column(const bigmod *, lidia_size_t, lidia_size_t,
			lidia_size_t from = 0);
	void sto_column_vector(const base_vector< bigmod > &, lidia_size_t,
			       lidia_size_t, lidia_size_t from = 0);

	bigmod * get_column(lidia_size_t i) const
	{
		register bigmod *RES;
		if (BigintLong) {
			RES = new bigmod[bigint_value.rows];
			memory_handler(RES, DMESSAGE, "get_column(lidia_size_t) :: "
				       "Error in memory allocation (RES)");
		}
		else {
			RES = new bigmod[long_value.rows];
			memory_handler(RES, DMESSAGE, "get_column(lidia_size_t) :: "
				       "Error in memory allocation (RES)");
		}
		get_column(RES, i);
		return RES;
	}

	base_vector< bigmod > get_column_vector(lidia_size_t i) const
	{
		base_vector< bigmod > RES;
		get_column_vector(RES, i);
		return RES;
	}

	void get_column(bigmod *, lidia_size_t) const;
	void get_column_vector(base_vector< bigmod > &, lidia_size_t) const;

	//
	// row access
	//

public:

	base_vector< bigmod > operator[] (lidia_size_t i) const
	{
		base_vector< bigmod > RES;
		get_row_vector(RES, i);
		return RES;
	}

	void sto_row(const bigmod *, lidia_size_t, lidia_size_t, lidia_size_t from = 0);
	void sto_row_vector(const base_vector< bigmod > &, lidia_size_t, lidia_size_t,
			    lidia_size_t from = 0);

	bigmod *get_row(lidia_size_t i) const
	{
		register bigmod *RES;
		if (BigintLong) {
			RES = new bigmod[bigint_value.columns];
			memory_handler(RES, DMESSAGE, "get_row(lidia_size_t) :: "
				       "Error in memory allocation (RES)");
		}
		else {
			RES = new bigmod[long_value.columns];
			memory_handler(RES, DMESSAGE, "get_row(lidia_size_t) :: "
				       "Error in memory allocation (RES)");
		}

		get_row(RES, i);
		return RES;
	}

	base_vector< bigmod > get_row_vector(lidia_size_t i) const
	{
		base_vector< bigmod > RES;
		get_row_vector(RES, i);
		return RES;
	}

	void get_row(bigmod *, lidia_size_t) const;
	void get_row_vector(base_vector< bigmod > &, lidia_size_t) const;

	//
	// array access
	//

public:

	bigmod ** get_data() const;
	void set_data(const bigint **, lidia_size_t, lidia_size_t);

	//
	// modulus access
	//

	bigint get_modulus() const
	{
		return (BigintLong ? bigint_modulus : bigint(long_modulus));
	}

	void set_modulus(const bigint &);

	///////////////////////////
	// END: access functions //
	///////////////////////////

	//
	// insert_at
	//

public:

	void insert_at(lidia_size_t, lidia_size_t,
		       const base_matrix< bigmod > &, lidia_size_t,
		       lidia_size_t, lidia_size_t, lidia_size_t);

	//
	// insert_columns, insert_rows, remove_columns, remove_rows
	//

public:

	void insert_columns(lidia_size_t *, const bigint **);
	void remove_columns(lidia_size_t *);

	void insert_rows(lidia_size_t *, const bigint **);
	void remove_rows(lidia_size_t *);

	//
	// split functions
	//

public:

	void split_t(base_matrix< bigmod > &, base_matrix< bigmod > &,
		     base_matrix< bigmod > &, base_matrix< bigmod > &) const;

	void split_h(base_matrix< bigmod > &, base_matrix< bigmod > &) const;
	void split_h(bigmod *, base_matrix< bigmod > &) const;
	void split_h(base_matrix< bigmod > &, bigmod *) const;
	void split_h(base_vector< bigmod > &, base_matrix< bigmod > &) const;
	void split_h(base_matrix< bigmod > &, base_vector< bigmod > &) const;

	void split_v(base_matrix< bigmod > &, base_matrix< bigmod > &) const;
	void split_v(bigmod *, base_matrix< bigmod > &) const;
	void split_v(base_matrix< bigmod > &, bigmod *) const;
	void split_v(base_vector< bigmod > &, base_matrix< bigmod > &) const;
	void split_v(base_matrix< bigmod > &, base_vector< bigmod > &) const;

	//
	// compose functions
	//

public:

	void compose_t(const base_matrix< bigmod > &,
		       const base_matrix< bigmod > &,
		       const base_matrix< bigmod > &,
		       const base_matrix< bigmod > &);

	void compose_h(const base_matrix< bigmod > &,
		       const base_matrix< bigmod > &);
	void compose_h(const bigmod *, const base_matrix< bigmod > &);
	void compose_h(const base_matrix< bigmod > &, const bigmod *);
	void compose_h(const base_vector< bigmod > &,
		       const base_matrix< bigmod > &);
	void compose_h(const base_matrix< bigmod > &,
		       const base_vector< bigmod > &);

	void compose_v(const base_matrix< bigmod > &,
		       const base_matrix< bigmod > &);
	void compose_v(const bigmod *, const base_matrix< bigmod > &);
	void compose_v(const base_matrix< bigmod > &, const bigmod *);
	void compose_v(const base_vector< bigmod > &,
		       const base_matrix< bigmod > &);
	void compose_v(const base_matrix< bigmod > &,
		       const base_vector< bigmod > &);

	//
	// exchange functions / swap functions
	//

public:

	void swap(base_matrix< bigmod > &);

public:

	void swap_columns(lidia_size_t, lidia_size_t);
	void swap_rows(lidia_size_t, lidia_size_t);

	//
	// structur functions
	//

public:

	lidia_size_t get_no_of_rows() const
	{
		if (BigintLong)
			return bigint_value.rows;
		else
			return long_value.rows;
	}

	lidia_size_t get_no_of_columns() const
	{
		if (BigintLong)
			return bigint_value.columns;
		else
			return long_value.columns;
	}

	void set_no_of_rows(lidia_size_t);
	void set_no_of_columns(lidia_size_t);

	void resize(lidia_size_t, lidia_size_t);

	void reset()
	{
		kill();
	}

	void kill();

	//
	// assignment
	//

public:

	base_matrix< bigmod > & operator = (const base_matrix< bigmod > &B)
	{
		this->assign(B);
		return *this;
	}

	void assign(const base_matrix< bigmod > &);

	//
	// diagonal function
	//

public:

	void diag(const bigmod &, const bigmod &);

	//
	// transpose function
	//

public:

	void trans(const base_matrix< bigmod > &);

	base_matrix< bigmod > trans() const;

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

	//
	// boolean functions
	//

public:

	bool is_column_zero(lidia_size_t) const;

	bool is_row_zero(lidia_size_t) const;

	bool is_matrix_zero() const;

	//
	// internal functions
	//

private:

	void status_report();
	void change_representation(unsigned long);
};



#ifdef LIDIA_NAMESPACE
}	// end of namespace LiDIA
# undef IN_NAMESPACE_LIDIA
#endif



#undef DBMKex
#undef SBMKex



#ifdef LIDIA_INCLUDE_CC
# include	"LiDIA/base_matrix.cc"
#endif



#endif	// LIDIA_BASE_MATRIX_H_GUARD_
