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


#ifndef LIDIA_SPARSE_BASE_MATRIX_CC_GUARD_
#define LIDIA_SPARSE_BASE_MATRIX_CC_GUARD_



#include "LiDIA/sparse_base_matrix.h"
#include "LiDIA/precondition_error.h"



#ifdef LIDIA_NAMESPACE
# ifndef IN_NAMESPACE_LIDIA
namespace LiDIA {
# endif
#endif



//
// debug defines / error defines
//

extern const char *matrix_error_msg[];
extern const char *PRT;

#define DVALUE   LDBL_MATRIX             // Debug value
#define DMESSAGE "sparse_base_matrix"    // Debug message
#define EMESSAGE matrix_error_msg        // Error message array

//
// debug level
//
//   0 : constructors / destructor
//   1 : input / output
//   2 : access functions
//   3 : split / compose
//   4 : swap functions
//   5 : structur functions
//   6 : assignments
//   7 : diagonal and transpose functions
//   8 : stream handling
//   9 : internal functions
//  10 : boolean functions
//

//
// constructors
//

template< class T >
sparse_base_matrix< T >::sparse_base_matrix ()
{
	debug_handler_l(DMESSAGE, "sparse_base_matrix()", DVALUE);

	this->bitfield.set_representation(matrix_flags::sparse_representation);
	this->SS_base_modul.constructor(*this);
}



template< class T >
sparse_base_matrix< T >::sparse_base_matrix (lidia_size_t r, lidia_size_t c)
{
	debug_handler_l(DMESSAGE, "sparse_base_matrix(lidia_size_t, lidia_size_t)",
			DVALUE);

	if (r < 0 || c < 0)
		precondition_error_handler(r, "r", "r >= 0",
				    c, "c", "c >= 0",
				    "sparse_base_matrix< T >::"
				    "sparse_base_matrix(lidia_size_t r, lidia_size_t c)",
				    DMESSAGE, EMESSAGE[0]);

	this->bitfield.set_representation(matrix_flags::sparse_representation);
	this->S_base_modul.constructor(*this, r, c);
}



template< class T >
sparse_base_matrix< T >::sparse_base_matrix (const base_vector< T > &v)
{
	debug_handler_l(DMESSAGE, "sparse_base_matrix(const base_vector< T > &)", DVALUE);

	this->bitfield.set_representation(matrix_flags::sparse_representation);
	this->SS_base_modul.constructor(*this, v);
}



template< class T >
sparse_base_matrix< T >::sparse_base_matrix (const sparse_base_matrix< T > &M)
{
	debug_handler_l(DMESSAGE, "sparse_base_matrix(const MR< T > &)",
			DVALUE);

	this->bitfield.set_representation(matrix_flags::sparse_representation);
	this->S_base_modul.constructor(*this, M);
}



template< class T >
sparse_base_matrix< T >::sparse_base_matrix (lidia_size_t r, lidia_size_t c, const T **A)
{
	debug_handler_l(DMESSAGE,
			"sparse_base_matrix(lidia_size_t, lidia_size_t, const T **)",
			DVALUE);

	if (r < 0 || c < 0 || A == NULL)
		precondition_error_handler(r, "r", "r >= 0",
				    c, "c", "c >= 0",
				    PRT, "A", "A != NULL",
				    "sparse_base_matrix< T >::"
				    "sparse_base_matrix(lidia_size_t r, lidia_size_t c, "
				    "const T **A)",
				    DMESSAGE, EMESSAGE[1]);

	this->bitfield.set_representation(matrix_flags::sparse_representation);
	this->S_base_modul.constructor(*this, r, c, A);
}



//
// destructor
//

template< class T >
sparse_base_matrix< T >::~sparse_base_matrix()
{
	debug_handler_l(DMESSAGE, "~sparse_base_matrix()", DVALUE);

	this->S_base_modul.destructor(*this);
}



//
// Input / Output
//

template< class T >
std::ostream &
sparse_base_matrix< T >::write (std::ostream &out) const
{
	//
	//    Task: Output to stream
	// Version: 2.0
	//

	debug_handler_l(DMESSAGE, "print(std::ostream &)", DVALUE + 1);

	switch(this->bitfield.get_print_mode()) {
	case matrix_flags::beauty_mode:
		this->S_base_modul.write_to_beauty(*this, out);
		break;
	case matrix_flags::lidia_mode:
		this->S_base_modul.write_to_stream(*this, out);
		break;
	case matrix_flags::gp_mode:
		this->S_base_modul.write_to_gp(*this, out);
		break;
	case matrix_flags::maple_mode:
		this->S_base_modul.write_to_maple(*this, out);
		break;
	case matrix_flags::mathematica_mode:
		this->S_base_modul.write_to_mathematica(*this, out);
		break;
	case matrix_flags::kash_mode:
		this->S_base_modul.write_to_kash(*this, out);
		break;
	case matrix_flags::latex_mode:
		this->S_base_modul.write_to_latex(*this, out);
		break;
	case matrix_flags::magma_mode:
		this->S_base_modul.write_to_magma(*this, out);
		break;
	}
	return out;
}



template< class T >
std::istream &
sparse_base_matrix< T >::read (std::istream &in)
{
	//
	//       Task: Input from stream
	// Conditions: first character not in {[, {, 1, 2, 3, 4, 5, 6, 7, 8, 9, M, a}
	//    Version: 2.0
	//

	debug_handler_l(DMESSAGE, "read(std::istream &)", DVALUE + 1);

	char c;
	in >> c;
	in.putback(c);

	switch(c) {
	case '1':
	case '2':
	case '3':
	case '4':
	case '5':
	case '6':
	case '7':
	case '8':
	case '9':
		this->S_base_modul.read_from_stream(*this, in);
		break;
	case '[':
		this->S_base_modul.read_from_gp(*this, in);
		break;
	case 'a':
		this->S_base_modul.read_from_maple(*this, in);
		break;
	case '{':
		this->S_base_modul.read_from_mathematica(*this, in);
		break;
	case 'M':
		this->S_base_modul.read_from_kash(*this, in);
		break;
	default:
		precondition_error_handler(c, "c", "Character c have to be in "
				    "{[, {, a, M, 1, 2, 3, 4, 5, 6, 7, 8, 9}",
				    "std::istream & sparse_base_matrix< T >::"
				    "read(std::istream &in)",
				    DMESSAGE, EMESSAGE[2]);
	}
	return in;
}



/////////////////////////////
// BEGIN: access functions //
/////////////////////////////

//
// element access
//

template< class T >
void
sparse_base_matrix< T >::sto (lidia_size_t x, lidia_size_t y, const T &e)
{
	//
	//       Task: A.sto(x, y, e) stores e at position (x, y) in matrix A.
	// Conditions: 0 <= x < A.rows and
	//             0 <= y < A.columns
	//    Version: 2.0
	//

	debug_handler_l(DMESSAGE, "sto(lidia_size_t, lidia_size_t, const T &)", DVALUE + 2);

	if (x< 0 || x >= this->rows || y< 0 || y >= this->columns)
		precondition_error_handler(x, "x", "0 <= x < rows",
				    this->rows, "rows", "",
				    y, "y", "0 <= y < columns",
				    this->columns, "", "",
				    "void sparse_base_matrix< T >::"
				    "sto(lidia_size_t x, lidia_size_t y, const T &e)",
				    DMESSAGE, EMESSAGE[3]);

	this->S_base_modul.sto(*this, x, y, e);
}



template< class T >
const T &
sparse_base_matrix< T >::member (lidia_size_t x, lidia_size_t y) const
{
	//
	//       Task: A.member(x, y) = value at position (x, y) in matrix A
	// Conditions: 0 <= x < A.rows and
	//             0 <= y < A.columns
	//    Version: 2.0
	//

	debug_handler_l(DMESSAGE, "member(lidia_size_t, lidia_size_t)", DVALUE + 2);

	if (x< 0 || y < 0 || x >= this->rows || y >= this->columns)
		precondition_error_handler(x, "x", "0 <= x < rows",
				    this->rows, "rows", "",
				    y, "y", "0 <= y < columns",
				    this->columns, "columns", "",
				    "T sparse_base_matrix< T >::"
				    "member(lidia_size_t x, lidia_size_t y) const",
				    DMESSAGE, EMESSAGE[3]);

	return this->S_base_modul.member(*this, x, y);
}



//
// column access
//

template< class T >
void
sparse_base_matrix< T >::sto_column (const T *v, lidia_size_t l, lidia_size_t j, lidia_size_t from)
{
	//
	//       Task: A.sto_column(v, l, j, p) stores v[0], ..., v[l-1] in column j
	//             of matrix A starting at position p.
	// Conditions: 0 <= j < A.columns and
	//             0 <= l <= A.rows and
	//             0 <= from < A.rows and
	//             from + l <= A.rows and
	//             v != NULL
	//    Version: 2.0
	//

	debug_handler_l(DMESSAGE, "sto_column(const T *, lidia_size_t, lidia_size_t, "
			"lidia_size_t)", DVALUE + 2);

	if (j >= this->columns || j< 0 || l > this->rows || l < 0 ||
	    from< 0 || from >= this->rows || from + l > this->rows || v == NULL)
		precondition_error_handler(j, "j", "0 <= j < columns",
				    this->columns, "columns", "",
				    l, "l", "0 <= l <= rows and l+from <= rows",
				    from, "from", "0 <= from < rows and l + from <= rows",
				    this->rows, "rows", "",
				    PRT, "v", " v != NULL",
				    "void sparse_base_matrix< T >::"
				    "sto_column(const T *v, lidia_size_t l, "
				    "lidia_size_t j, lidia_size_t from)",
				    DMESSAGE, EMESSAGE[1]);

	this->S_base_modul.sto_column(*this, v, l, j, from);
}



template< class T >
void
sparse_base_matrix< T >::sto_column_vector (const base_vector< T > &v, lidia_size_t l,
					    lidia_size_t j, lidia_size_t from)
{
	//
	//       Task: A.sto_column_vector(v, l, j, p) stores v[0], ..., v[l-1] in column j
	//             of matrix A starting at position p.
	// Conditions: 0 <= j < A.columns and
	//             0 <= l <= min(v.size, A.rows) and
	//             0 <= from < A.rows and
	//             from + l <= A.rows
	//    Version: 2.0
	//

	debug_handler_l(DMESSAGE,
			"sto_column_vector(const base_vector< T > &, "
			"lidia_size_t, lidia_size_t, lidia_size_t)",
			DVALUE + 2);

	lidia_size_t min = (this->rows < v.get_size() ? this->rows : v.get_size());

	if (j >= this->columns || j< 0 || l < 0 || l > min
	    || from< 0 || from >= this->rows || l + from > this->rows)
		precondition_error_handler(j , "j", "0 <= j < columns",
				    this->columns, "columns", "",
				    l, "l", "0 < l <= min(v.size, rows) and from + l <= rows",
				    from, "from", "0 <= from < rows and from + l <= rows",
				    this->rows, "rows", "",
				    v.get_size(), "v.size", "0 < l <= min(v.size, rows)",
				    "void sparse_base_matrix< T >::"
				    "sto_column_vector(const base_vector< T > &v, "
				    "lidia_size_t l, lidia_size_t j, lidia_size_t from)",
				    DMESSAGE, EMESSAGE[3]);

	this->S_base_modul.sto_column(*this, v.get_data_address(), l, j, from);
}



template< class T >
void
sparse_base_matrix< T >::get_column (T *RES, lidia_size_t i) const
{
	//
	//       Task: A.get_column(RES, i);
	// =  > RES[0] = A.value[0][i], ..., RES[A.rows-1] = A.value[A.rows-1][i]
	// Conditions: 0 <= i < A.columns and
	//             RES != NULL
	//    Version: 2.0
	//

	debug_handler_l(DMESSAGE, "get_column(T *, lidia_size_t)", DVALUE + 2);

	if (i >= this->columns || i < 0 || RES == NULL)
		precondition_error_handler(i, "i", "0 <= i < columns",
				    this->columns, "columns", "",
				    PRT, "RES", "RES != NULL",
				    "void sparse_base_matrix< T >::"
				    "get_column(T *RES, lidia_size_t i) const",
				    DMESSAGE, EMESSAGE[1]);

	this->S_base_modul.get_column(*this, RES, i);
}



template< class T >
void
sparse_base_matrix< T >::get_column_vector (base_vector< T > &RES, lidia_size_t i) const
{
	//
	//       Task: A.get_column_vector(RES, i);
	// =  > RES[0] = A.value[0][i], ..., RES[A.rows-1] = A.value[A.rows-1][i]
	// Conditions: 0 <= i < A.columns
	//    Version: 2.0
	//

	debug_handler_l(DMESSAGE, "get_column_vector(base_vector< T > &, lidia_size_t)",
			DVALUE + 2);

	if (i >= this->columns || i < 0)
		precondition_error_handler(i, "i", "0 <= i < columns",
				    this->columns, "columns", "",
				    "void sparse_base_matrix< T >::"
				    "get_column_vector(base_vector< T > &RES, lidia_size_t i) const",
				    DMESSAGE, EMESSAGE[3]);

	if (RES.get_capacity() < this->rows)
		RES.set_capacity(this->rows);
	if (RES.size() != this->rows)
		RES.set_size(this->rows);

	this->S_base_modul.get_column(*this, RES.get_data_address(), i);
}



//
// row access
//

template< class T >
void
sparse_base_matrix< T >::sto_row (const T *v, lidia_size_t l, lidia_size_t i, lidia_size_t from)
{
	//
	//       Task: A.sto_row(v, l, j, p) stores v[0], ..., v[l-1] in row i
	//             of matrix A starting at position p.
	// Conditions: 0 <= i < A.rows and
	//             0 <= l <= A.columns and
	//             0 <= from < A.columns and
	//             v != NULL and
	//             from + l <= A.rows
	//    Version: 2.0
	//

	debug_handler_l(DMESSAGE, "sto_row(const T *, lidia_size_t, "
			"lidia_size_t, lidia_size_t)", DVALUE + 2);

	if (i< 0 || i >= this->rows || 0 > from || from >= this->columns ||
	    l >= this->columns || l <= 0 || from + l > this->columns || v == NULL)
		precondition_error_handler(i, "i", "0 <= i < rows",
				    this->rows, "rows", "",
				    l, "l", "0 <= l <= columns and l + from <= columns" ,
				    from, "from", "0 <= from < columns and l + from <= columns",
				    this->columns, "columns", "",
				    PRT, "v", "v != NULL",
				    "void sparse_base_matrix< T >::"
				    "sto_row(const T *v, lidia_size_t l, "
				    "lidia_size_t i, lidia_size_t)",
				    DMESSAGE, EMESSAGE[1]);

	this->S_base_modul.sto_row(*this, v, l, i, from);
}



template< class T >
void
sparse_base_matrix< T >::sto_row_vector (const base_vector< T > &v, lidia_size_t l,
					 lidia_size_t i, lidia_size_t from)
{
	//
	//       Task: A.sto_row_vector(v, l, j, p) stores v[0], ..., v[l-1] in row i
	//             of matrix A starting at position p.
	// Conditions: 0 <= i < A.rows and
	//             0 <= l <= min(v.size, A.columns) and
	//             0 <= from < A.rows and
	//             from + l <= A.columns
	//    Version: 2.0
	//

	debug_handler_l(DMESSAGE, "sto_row_vector(const base_vector< T > &, "
			"lidia_size_t, lidia_size_t, lidia_size_t)", DVALUE + 2);

	lidia_size_t min = (this->columns < v.get_size() ? this->columns : v.get_size());

	if (i >= this->rows || i< 0 || l < 0 || l > min
	    || from< 0 || from >= this->columns || l + from > this->columns)
		precondition_error_handler(i, "i", "0 <= i < rows",
				    this->rows, "rows", "",
				    l, "l", "0 < l <= min(columns, v.size) and from + l <= columns",
				    from, "from", "0 <= from < columns and from + l <= columns",
				    this->columns, "columns", "",
				    v.size(), "v.size", "0 < l <= min(columns, v.size)",
				    "void sparse_base_matrix< T >::"
				    "sto_row_vector(const base_vector< T > &v, lidia_size_t l, "
				    "lidia_size_t i, lidia_size_t from)",
				    DMESSAGE, EMESSAGE[3]);

	this->S_base_modul.sto_row(*this, v.get_data_address(), l, i, from);
}



template< class T >
void
sparse_base_matrix< T >::get_row (T *RES, lidia_size_t i) const
{
	//
	//       Task: A.get_row(RES, i);
	// =  > RES[0] = A.value[i][0], ...,
	//                RES[A.columns-1] = A.value[i][A.columns-1]
	// Conditions: 0 <= i < A.rows and
	//             RES != NULL
	//    Version: 2.0
	//

	debug_handler_l(DMESSAGE, "get_row(T *, lidia_size_t)", DVALUE + 2);

	if (i >= this->rows || i < 0 || RES == NULL)
		precondition_error_handler(i, "i", "0 <= i < rows",
				    this->rows, "rows", "",
				    PRT, "RES", "RES != NULL",
				    "void sparse_base_matrix< T >::"
				    "get_row(T *RES, lidia_size_t i) const",
				    DMESSAGE, EMESSAGE[3]);

	this->S_base_modul.get_row(*this, RES, i);
}



template< class T >
void
sparse_base_matrix< T >::get_row_vector (base_vector< T > &RES, lidia_size_t i) const
{
	//
	//       Task: A.get_row_vector(RES, i);
	// =  > RES[0] = A.value[i][0], ...,
	//                RES[A.columns-1] = A.value[i][A.columns-1]
	// Conditions: 0 <= i < A.rows
	//    Version: 2.0
	//

	debug_handler_l(DMESSAGE, "get_row_vector(base_vector< T > &, lidia_size_t)",
			DVALUE + 2);

	if (i >= this->rows || i < 0)
		precondition_error_handler(i, "i", "0 <= i < rows",
				    this->rows, "rows", "",
				    "void sparse_base_matrix< T >::"
				    "row_vector(base_vector< T > &RES, lidia_size_t i) const",
				    DMESSAGE, EMESSAGE[3]);

	if (RES.get_capacity() < this->columns)
		RES.set_capacity(this->columns);
	if (RES.size() != this->columns)
		RES.set_size(this->columns);

	this->S_base_modul.get_row(*this, RES.get_data_address(), i);
}



//
// array access
//

template< class T >
T **
sparse_base_matrix< T >::get_data () const
{
	//
	//    Task: prt = A.get_data()
	// =  > prt = pointer to a copy of array A.value
	// Version: 2.0
	//

	debug_handler_l(DMESSAGE, "get_value()", DVALUE + 2);

	register T **copy = NULL;
	if (this->rows != 0) {
		copy = new T *[this->rows];
		memory_handler(copy, DM_BM, "get_value :: "
			       "Error in memory allocation (copy)");
	}

	register T *tmp, *tmp1;
	for (register lidia_size_t i = 0; i < this->rows; i++) {
		tmp1 = this->value[i];
		tmp = new T[this->columns];
		memory_handler(tmp, DM_BM, "get_value :: "
			       "Error in memory allocation (tmp)");

		register lidia_size_t l = 0;
		for (register lidia_size_t j = 0; l < this->value_counter[i] && j < this->columns; j++) {
			if (j == this->index[i][l]) {
				tmp[j] = tmp1[l];
				l++;
			}
			else
				tmp[j] = this->Zero;
		}
		copy[i] = tmp;
	}
	return copy;
}



template< class T >
void
sparse_base_matrix< T >::set_data (const T **v, lidia_size_t r, lidia_size_t c)
{
	//
	//       Task: A.set_data(v, r, c)
	// =  > A.value[i][j] = v[i][j]
	//                0 <= i < r and 0 <= j < c
	// Conditions: r >= 0 and
	//             c >= 0 and
	//             v != NULL
	//    Version: 2.0
	//

	debug_handler_l(DMESSAGE,
			"set_value(const T **, lidia_size_t, lidia_size_t)",
			DVALUE + 2);

	if (r < 0 || c < 0 || v == NULL)
		precondition_error_handler(r, "r", "r >= 0",
				    c, "c", "c >= 0",
				    PRT, "v", "v != NULL",
				    "void sparse_base_matrix< T >::"
				    "set_data(const T **v, lidia_size_t r, lidia_size_t c)",
				    DMESSAGE, EMESSAGE[1]);

	if (r != this->rows)
		this->S_base_modul.set_no_of_rows(*this, r);
	if (c != this->columns)
		this->S_base_modul.set_no_of_columns(*this, c);

	for (register lidia_size_t len1 = this->rows - 1; len1 >= 0; len1--)
		for (register lidia_size_t len = this->columns - 1; len >= 0; len--)
			this->S_base_modul.sto(*this, len1, len, v[len1][len]);
}



///////////////////////////
// END: access functions //
///////////////////////////

//
// insert_at
//

template< class T >
void
sparse_base_matrix< T >::insert_at (lidia_size_t r, lidia_size_t c, const sparse_base_matrix< T > &A,
				    lidia_size_t startr, lidia_size_t startc,
				    lidia_size_t l1, lidia_size_t l2)
{
	//
	//       Task: A.insert_at(r, c, B, startr, startc, l1, l2)
	//             stores the submatrix B[startr, startc]...B[startr+l1][startc+l2]
	//             in matrix A starting at position (r, c)
	// Conditions: 0 <= r < A.rows and
	//             0 <= c < A.columns and
	//             0 <= startr < B.rows and
	//             0 <= startc < B.columns and
	//             0 <= l1 <= B.rows and
	//             0 <= l2 <= B.columns and
	//             startr + l1 <= B.rows and
	//             startc + l2 <= B.columns and
	//             r + l1 <= A.rows and
	//             c + l2 <= A.columns
	//             v != NULL
	// Version: 2.0
	//

	if (r< 0 || r >= this->rows || c< 0 || c >= this->columns
	    || startr< 0 || startr >= A.rows || startc< 0 || startc >= A.columns
	    || l1< 0 || l1 > A.rows || l2< 0 || l2 > A.columns
	    || startr + l1 > A.rows || startc + l2 > A.columns
	    || r + l1 > this->rows || c + l2 > this->columns)
		precondition_error_handler(r, "r", "0 <= r < this->rows and r + l1 <= rows",
				    this->rows, "rows", "",
				    c, "c", "0 <= c < columns and c + l2 <= columns",
				    this->columns, "columns", "",
				    startr, "startr", "0 <= startr < A.rows and startr + l1 <= A.rows",
				    A.rows, "A.rows", "",
				    startc, "startc", "0 <= startc < A.columns and startc + l2 <= A.columns",
				    A.columns, "A.columns", "",
				    l1, "l1", "r + l1 <= rows and startr + l1 <= A.rows",
				    l2, "l2", "c + l2 <= columns and startc + l2 <= A.columns",
				    "void sparse_base_matrix< T >::"
				    "insert_at(lidia_size_t r, lidia_size_t c, const sparse_base_matrix< T > &A, "
				    "lidia_size_t startr, lidia_size_t startc, lidia_size_t l1, lidia_size_t l2)",
				    DMESSAGE, EMESSAGE[1]);

	this->SS_base_modul.insert_at(*this, r, c, A, startr, startc, l1, l2);
}



//
// insert_columns, insert_rows, remove_columns, remove_rows
//

template< class T >
void
sparse_base_matrix< T >::insert_columns (lidia_size_t *ind, const T **news)
{
	//
	//       Task: A.insert_columns(ind, news);
	// =  > insert the columns of array in matrix A
	//    Version: 2.0
	//

	debug_handler_l(DMESSAGE, "insert_columns(lidia_size_t *, const T **)", DVALUE + 2);

	this->S_base_modul.insert_columns(*this, ind, news);
}



template< class T >
void
sparse_base_matrix< T >::remove_columns (lidia_size_t *rem)
{
	//
	//       Task: A.remove_column(rem);
	// =  > remove columns
	//    Version: 2.0
	//

	debug_handler_l(DMESSAGE, "remove_column(lidia_size_t *)", DVALUE + 2);

	this->S_base_modul.remove_columns(*this, rem);
}



template< class T >
void
sparse_base_matrix< T >::insert_rows (lidia_size_t *ind, const T **news)
{
	//
	//       Task: A.insert_rows(ind, news);
	// =  > insert the rows of array in matrix A
	//    Version: 2.0
	//

	debug_handler_l(DMESSAGE, "insert_rows(lidia_size_t *, const T **)", DVALUE + 2);

	this->S_base_modul.insert_rows(*this, ind, news);
}



template< class T >
void
sparse_base_matrix< T >::remove_rows (lidia_size_t *rem)
{
	//
	//       Task: A.remove_rows(rem);
	// =  > remove rows
	//    Version: 2.0
	//

	debug_handler_l(DMESSAGE, "remove_rows(lidia_size_t *)", DVALUE + 2);

	this->S_base_modul.remove_rows(*this, rem);
}



//
// split functions
//

template< class T >
void
sparse_base_matrix< T >::split_t (sparse_base_matrix< T > &A, sparse_base_matrix< T > &B,
				  sparse_base_matrix< T > &C, sparse_base_matrix< T > &D) const
{
	//
	//       Task: E.split(A, B, C, D);
	// =  > splits matrix E into the submatrices A, B, C, D
	//                    (A B)
	// =  > E = (C D)
	// Conditions: A.columns, B.columns, C.columns, D.columns <= columns and
	//             A.rows, B.rows, C.rows, D.rows <= rows
	//    Version: 2.0
	//

	debug_handler_l(DMESSAGE,
			"split(sparse_base_matrix< T > &, sparse_base_matrix< T > &, "
			"sparse_base_matrix< T > &, sparse_base_matrix< T > &)",
			DVALUE + 3);

	if (A.rows > this->rows || A.columns > this->columns ||
	    B.rows > this->rows || B.columns > this->columns ||
	    C.rows > this->rows || C.columns > this->columns ||
	    D.rows > this->rows || D.columns > this->columns)
		precondition_error_handler(A.rows, "A.rows", "A.rows <= rows",
				    A.columns, "A.columns", "A.columns <= columns",
				    B.rows, "B.rows", "B.rows <= rows",
				    B.columns, "B.columns", "B.columns <= columns",
				    C.rows, "C.rows", "C.rows <= rows",
				    C.columns, "C.columns", "C.columns <= columns",
				    D.rows, "D.rows", "D.rows <= rows",
				    D.columns, "D.columns", "D.columns <= columns",
				    this->rows, "rows", "",
				    this->columns, "columns", "",
				    "void sparse_base_matrix< T >::"
				    "split(sparse_base_matrix< T > &A, sparse_base_matrix< T > &B, "
				    "sparse_base_matrix< T > &C, sparse_base_matrix< T > &D) "
				    "const",
				    DMESSAGE, EMESSAGE[4]);

	this->SS_base_modul.insert_at(A, 0, 0, *this, 0, 0, A.rows, A.columns);
	this->SS_base_modul.insert_at(B, 0, 0, *this, 0, this->columns - B.columns, B.rows, B.columns);
	this->SS_base_modul.insert_at(C, 0, 0, *this, this->rows - C.rows, 0, C.rows, C.columns);
	this->SS_base_modul.insert_at(D, 0, 0, *this, this->rows - D.rows, this->columns - D.columns, D.rows, D.columns);
}



template< class T >
void
sparse_base_matrix< T >::split_h (sparse_base_matrix< T > &A, sparse_base_matrix< T > &B) const
{
	//
	//       Task: C.split_h(A, B);
	// =  > splits matrix C into the submatrices A, B
	// =  > C = (A B)
	// Conditions: A.columns, B.columns <= columns and
	//             A.rows, B.rows <= rows
	//    Version: 2.0
	//

	debug_handler_l(DMESSAGE,
			"split_h(sparse_base_matrix< T > &, sparse_base_matrix< T > &)",
			DVALUE + 3);

	if (A.rows > this->rows || A.columns > this->columns ||
	    B.rows > this->rows || B.columns > this->columns)
		precondition_error_handler(A.rows, "A.rows", "A.rows <= rows",
				    A.columns, "A.columns", "A.columns <= coluns",
				    B.rows, "B.rows", "B.rows <= rows",
				    B.columns, "B.columns", "B.columns <= coluns",
				    this->rows, "rows", "",
				    this->columns, "columns", "",
				    "void sparse_base_matrix< T >::"
				    "split_h(sparse_base_matrix< T > &A, "
				    "sparse_base_matrix< T > &B) const",
				    DMESSAGE, EMESSAGE[4]);

	this->SS_base_modul.insert_at(A, 0, 0, *this, 0, 0, A.rows, A.columns);
	this->SS_base_modul.insert_at(B, 0, 0, *this, 0, this->columns - B.columns, B.rows, B.columns);
}



template< class T >
void
sparse_base_matrix< T >::split_h (T *v, sparse_base_matrix< T > &A) const
{
	//
	//       Task: C.split_h(v, A);
	// =  > splits matrix C into the array v and the submatrix A
	// =  > C = (v A)
	// Conditions: A.columns <= columns and
	//             A.rows <= rows and
	//             v != NULL
	//    Version: 2.0
	//

	debug_handler_l(DMESSAGE,
			"split_h(T *, sparse_base_matrix< T > &)",
			DVALUE + 3);

	if (A.rows > this->rows || A.columns > this->columns || v == NULL)
		precondition_error_handler(A.rows, "A.rows", "A.rows <= rows",
				    A.columns, "A.columns", "A.columns <= columns",
				    this->rows, "rows", "",
				    this->columns, "columns", "",
				    PRT, "v", "v != NULL",
				    "void sparse_base_matrix< T >::"
				    "split_h(T *v, sparse_base_matrix< T > &A) const",
				    DMESSAGE, EMESSAGE[4]);

	this->SS_base_modul.insert_at(A, 0, 0, *this, 0, this->columns - A.columns, A.rows, A.columns);
	this->S_base_modul.get_column(*this, v, 0);
}



template< class T >
void
sparse_base_matrix< T >::split_h (base_vector< T > &v, sparse_base_matrix< T > &A) const
{
	//
	//       Task: C.split_h(v, A);
	// =  > splits matrix C into the vector v and the submatrix A
	// =  > C = (v A)
	// Conditions: v.size <= rows and
	//             A.columns <= columns and
	//             A.rows <= rows
	//    Version: 2.0
	//

	debug_handler_l(DMESSAGE,
			"split_h(base_vector< T > &, sparse_base_matrix< T > &)",
			DVALUE + 3);

	if (A.rows > this->rows || A.columns > this->columns || v.size() > this->rows)
		precondition_error_handler(A.rows, "A.rows", "A.rows <= rows",
				    A.columns, "A.columns", "A.columns <= columns",
				    this->rows, "rows", "",
				    this->columns, "columns", "",
				    v.size(), "v.size", "v.size <= rows",
				    "void sparse_base_matrix< T >::"
				    "split_h(base_vector< T > &v, "
				    "sparse_base_matrix< T > &A) const",
				    DMESSAGE, EMESSAGE[4]);

	this->SS_base_modul.insert_at(A, 0, 0, *this, 0, this->columns - A.columns, A.rows, A.columns);
	this->S_base_modul.get_column(*this, v.get_data_address(), 0);
}



template< class T >
void
sparse_base_matrix< T >::split_h (sparse_base_matrix< T > &A, T *v) const
{
	//
	//       Task: C.split_h(A, v);
	// =  > splits matrix C into the submatrix A and the array v
	// =  > C = (A v)
	// Conditions: A.columns <= this->columns and
	//             A.rows <= rows and
	//             v != NULL
	//    Version: 2.0
	//

	debug_handler_l(DMESSAGE,
			"split_h(sparse_base_matrix< T > &, T *)",
			DVALUE + 3);

	if (A.rows > this->rows || A.columns > this->columns || v == NULL)
		precondition_error_handler(A.rows, "A.rows", "A.rows <= rows",
				    A.columns, "A.columns", "A.columns <= columns",
				    this->rows, "rows", "",
				    this->columns, "columns", "",
				    PRT, "v", "v != NULL",
				    "void sparse_base_matrix< T >::"
				    "split_h(sparse_base_matrix< T > &A, T *v) const",
				    DMESSAGE, EMESSAGE[4]);

	this->SS_base_modul.insert_at(A, 0, 0, *this, 0, 0, A.rows, A.columns);
	this->S_base_modul.get_column(*this, v, this->columns - 1);
}



template< class T >
void
sparse_base_matrix< T >::split_h (sparse_base_matrix< T > &A, base_vector< T > &v) const
{
	//
	//       Task: C.split_h(A, v);
	// =  > splits matrix C into the submatrix A and the vector v
	// =  > C = (A v)
	// Conditions: A.columns <= columns and
	//             A.rows <= rows and
	//             v.size <= rows
	//    Version: 2.0
	//

	debug_handler_l(DMESSAGE,
			"split_h(sparse_base_matrix< T > &, base_vector< T > &)",
			DVALUE + 3);

	if (A.rows > this->rows || A.columns > this->columns || v.size() > this->rows)
		precondition_error_handler(A.rows, "A.rows", "A.rows <= rows",
				    A.columns, "A.columns", "A.columns <= columns",
				    this->rows, "rows", "",
				    this->columns, "columns", "",
				    v.size(), "v.size", "v.size <= rows",
				    "void sparse_base_matrix< T >::"
				    "split_h(sparse_base_matrix< T > &A, "
				    "base_vector< T > &v) const",
				    DMESSAGE, EMESSAGE[4]);

	this->SS_base_modul.insert_at(A, 0, 0, *this, 0, 0, A.rows, A.columns);
	this->S_base_modul.get_column(*this, v.get_data_address(), this->columns - 1);
}



template< class T >
void
sparse_base_matrix< T >::split_v (sparse_base_matrix< T > &A, sparse_base_matrix< T > &B) const
{
	//
	//       Task: C.split_v(A, B);
	// =  > splits matrix C into the submatrices A, B
	//                    (A)
	// =  > C = (B)
	// Conditions: A.columns, B.columns <= this->columns and
	//             A.rows, B.rows <= this->rows
	//    Version: 2.0
	//

	debug_handler_l(DMESSAGE,
			"split_v(sparse_base_matrix< T > &, sparse_base_matrix< T > &)",
			DVALUE + 3);

	if (A.rows > this->rows || A.columns > this->columns ||
	    B.rows > this->rows || B.columns > this->columns)
		precondition_error_handler(A.rows, "A.rows", "A.rows <= rows",
				    A.columns, "A.columns", "A.columns <= columns",
				    B.rows, "B.rows", "B.rows <= rows",
				    B.columns, "B.columns", "B.columns <= columns",
				    this->rows, "rows", "",
				    this->columns, "columns", "",
				    "void sparse_base_matrix< T >::"
				    "split_v(sparse_base_matrix< T > &A, "
				    "sparse_base_matrix< T > &B) const",
				    DMESSAGE, EMESSAGE[4]);

	this->SS_base_modul.insert_at(A, 0, 0, *this, 0, 0, A.rows, A.columns);
	this->SS_base_modul.insert_at(B, 0, 0, *this, this->rows - B.rows, 0, B.rows, B.columns);
}



template< class T >
void
sparse_base_matrix< T >::split_v (T *v, sparse_base_matrix< T > &A) const
{
	//
	//       Task: C.split_v(v, A);
	// =  > splits matrix C into the array v and the submatrix A
	//                    (v)
	// =  > C = (A)
	// Conditions: A.columns <= columns and
	//             A.rows <= rows and
	//             v != NULL
	//    Version: 2.0
	//

	debug_handler_l(DMESSAGE,
			"split_v(T *, sparse_base_matrix< T > &)",
			DVALUE + 3);

	if (A.rows > this->rows || A.columns > this->columns || v == NULL)
		precondition_error_handler(A.rows, "A.rows", "A.rows <= rows",
				    A.columns, "A.columns", "A.columns <= columns",
				    this->rows, "rows", "",
				    this->columns, "columns", "",
				    PRT, "v", "v != NULL",
				    "void sparse_base_matrix< T >::"
				    "split_v(T *v, sparse_base_matrix< T > &A) const",
				    DMESSAGE, EMESSAGE[4]);

	this->SS_base_modul.insert_at(A, 0, 0, *this, this->rows - A.rows, 0, A.rows, A.columns);
	this->S_base_modul.get_row(*this, v, 0);
}



template< class T >
void
sparse_base_matrix< T >::split_v (base_vector< T > &v, sparse_base_matrix< T > &A) const
{
	//
	//       Task: C.split_v(v, A);
	// =  > splits matrix C into the vector v and the submatrix A
	//                    (v)
	// =  > C = (A)
	// Conditions: A.columns <= columns and
	//             A.rows <= rows and
	//             v.size <= columns
	//    Version: 2.0
	//

	debug_handler_l(DMESSAGE,
			"split_v(base_vector< T > &, sparse_base_matrix< T > &)",
			DVALUE + 3);

	if (A.rows > this->rows || A.columns > this->columns || v.size() > this->columns)
		precondition_error_handler(A.rows, "A.rows", "A.rows <= rows",
				    A.columns, "A.columns", "A.columns <= columns",
				    this->rows, "rows", "",
				    this->columns, "columns", "",
				    v.size(), "v.size", "v.size <= columns",
				    "void sparse_base_matrix< T >::"
				    "split_v(base_vector< T > &v, "
				    "sparse_base_matrix< T > &A) const",
				    DMESSAGE, EMESSAGE[4]);

	this->SS_base_modul.insert_at(A, 0, 0, *this, this->rows - A.rows, 0, A.rows, A.columns);
	this->S_base_modul.get_row(*this, v.get_data_address(), 0);
}



template< class T >
void
sparse_base_matrix< T >::split_v (sparse_base_matrix< T > &A, T *v) const
{
	//
	//       Task: C.split_v(A, B);
	// =  > splits matrix C into the submatrix A and the vector v
	//                    (A)
	// =  > C = (v)
	// Conditions: A.columns <= columns and
	//             A.rows <= rows and
	//             v != NULL
	//    Version: 2.0
	//

	debug_handler_l(DMESSAGE,
			"split_v(sparse_base_matrix< T > &, T *)",
			DVALUE + 3);

	if (A.rows > this->rows || A.columns > this->columns || v == NULL)
		precondition_error_handler(A.rows, "A.rows", "A.rows <= rows",
				    A.columns, "A.columns", "A.columsn <= columns",
				    this->rows, "rows", "",
				    this->columns, "columns", "",
				    PRT, "v", "v != NULL",
				    "void sparse_base_matrix< T >::"
				    "split_v(sparse_base_matrix< T > &A, T *v) const",
				    DMESSAGE, EMESSAGE[4]);

	this->SS_base_modul.insert_at(A, 0, 0, *this, 0, 0, A.rows, A.columns);
	this->S_base_modul.get_row(*this, v, this->rows - 1);
}



template< class T >
void
sparse_base_matrix< T >::split_v (sparse_base_matrix< T > &A, base_vector< T > &v) const
{
	//
	//       Task: C.split_v(A, B);
	// =  > splits matrix C into the submatrix A and the vector v
	//                    (A)
	// =  > C = (v)
	// Conditions: A.columns <= columns and
	//             A.rows <= rows and
	//             v.size <= columns
	//    Version: 2.0
	//

	debug_handler_l(DMESSAGE,
			"split_v(sparse_base_matrix< T > &, base_vector< T > &)",
			DVALUE + 3);

	if (A.rows > this->rows || A.columns > this->columns || v.size() > this->columns)
		precondition_error_handler(A.rows, "A.rows", "A.rows <= rows",
				    A.columns, "A.columns", "A.columns <= columns",
				    this->rows, "rows", "",
				    this->columns, "columns", "",
				    v.size(), "v.size", "v.size <= columns",
				    "void sparse_base_matrix< T >::"
				    "split_v(sparse_base_matrix< T > &A, "
				    "base_vector< T > &v) const",
				    DMESSAGE, EMESSAGE[4]);

	this->SS_base_modul.insert_at(A, 0, 0, *this, 0, 0, A.rows, A.columns);
	this->S_base_modul.get_row(*this, v.get_data_address(), this->rows - 1);
}



//
// compose functions
//

template< class T >
void
sparse_base_matrix< T >::compose_t (const sparse_base_matrix< T > &A, const sparse_base_matrix< T > &B,
				    const sparse_base_matrix< T > &C, const sparse_base_matrix< T > &D)
{
	//
	//       Task: E.compose(A, B, C, D)
	// =  > compose the submatrices A, B, C, D to E
	//                    (A B)
	// =  > E = (C D)
	// Conditions: A.columns, B.columns, C.columns, D.columns <= columns and
	//             A.rows, B.rows, C.rows, D.rows <= rows
	//    Version: 2.0
	//

	debug_handler_l(DMESSAGE,
			"compose(const sparse_base_matrix< T > &, "
			"const sparse_base_matrix< T > &, "
			"const sparse_base_matrix< T > &, "
			"const sparse_base_matrix< T > &)",
			DVALUE + 3);

	if (A.rows > this->rows || A.columns > this->columns ||
	    B.rows > this->rows || B.columns > this->columns ||
	    C.rows > this->rows || C.columns > this->columns ||
	    D.rows > this->rows || D.columns > this->columns)
		precondition_error_handler(A.rows, "A.rows", "A.rows <= rows",
				    A.columns, "A.columns", "A.columns <= columns",
				    B.rows, "B.rows", "B.rows <= rows",
				    B.columns, "B.columns", "B.columns <= columns",
				    C.rows, "C.rows", "C.rows <= rows",
				    C.columns, "C.columns", "C.columns <= columns",
				    D.rows, "D.rows", "D.rows <= rows",
				    D.columns, "D.columns", "D.columns <= columns",
				    this->rows, "rows", "",
				    this->columns, "columns", "",
				    "void sparse_base_matrix< T >::"
				    "compose(const sparse_base_matrix< T > &A, "
				    "const sparse_base_matrix< T > &B, "
				    "const sparse_base_matrix< T > &C, "
				    "const sparse_base_matrix< T > &D)",
				    DMESSAGE, EMESSAGE[4]);

	this->SS_base_modul.insert_at(*this, 0, 0, A, 0, 0, A.rows, A.columns);
	this->SS_base_modul.insert_at(*this, 0, this->columns - B.columns, B, 0, 0, B.rows, B.columns);
	this->SS_base_modul.insert_at(*this, this->rows - C.rows, 0, C, 0, 0, C.rows, C.columns);
	this->SS_base_modul.insert_at(*this, this->rows - D.rows, this->columns - D.columns, D, 0, 0, D.rows, D.columns);
}



template< class T >
void
sparse_base_matrix< T >::compose_h (const sparse_base_matrix< T > &A, const sparse_base_matrix< T > &B)
{
	//
	//       Task: C.compose_h(A, B);
	// =  > compose the submatrices A, B to C
	// =  > C = (A B)
	// Conditions: A.columns, B.columns <= columns and
	//             A.rows, B.rows <= rows
	//    Version: 2.0
	//

	debug_handler_l(DMESSAGE,
			"compose_h(const sparse_base_matrix< T > &, "
			"const sparse_base_matrix< T > &)",
			DVALUE + 3);

	if (A.rows > this->rows || A.columns > this->columns ||
	    B.rows > this->rows || B.columns > this->columns)
		precondition_error_handler(A.rows, "A.rows", "A.rows <= rows",
				    A.columns, "A.columns", "A.columns <= columns",
				    B.rows, "B.rows", "B.rows <= rows",
				    B.columns, "B.columns", "B.columns <= columns",
				    this->rows, "rows", "",
				    this->columns, "columns", "",
				    "void sparse_base_matrix< T >::"
				    "compose_h(const sparse_base_matrix< T > &A, "
				    "const sparse_base_matrix< T > &B)",
				    DMESSAGE, EMESSAGE[4]);

	this->SS_base_modul.insert_at(*this, 0, 0, A, 0, 0, A.rows, A.columns);
	this->SS_base_modul.insert_at(*this, 0, this->columns - B.columns, B, 0, 0, B.rows, B.columns);
}



template< class T >
void
sparse_base_matrix< T >::compose_h (const T *v, const sparse_base_matrix< T > &A)
{
	//
	//       Task: C.compose_h(v, A);
	// =  > compose the submatrix A and the vector v to matrix C
	// =  > C = (v A)
	// Conditions: A.columns <= columns and
	//             A.rows <= rows and
	//             v != NULL
	//    Version: 2.0
	//

	debug_handler_l(DMESSAGE, "in member - function "
			"compose_h(const T *, const sparse_base_matrix< T > &)",
			DVALUE + 3);

	if (A.rows > this->rows || A.columns > this->columns || v == NULL)
		precondition_error_handler(A.rows, "A.rows", "A.rows <= rows",
				    A.columns, "A.columns", "A.columns <= columns",
				    this->rows, "rows", "",
				    this->columns, "columns", "",
				    PRT, "v", "v != NULL",
				    "void sparse_base_matrix< T >::"
				    "compose_h(const T *v, const sparse_base_matrix< T > &A)",
				    DMESSAGE, EMESSAGE[4]);

	this->S_base_modul.insert_column_at(*this, 0, v, this->rows);
	this->SS_base_modul.insert_at(*this, 0, this->columns - A.columns, A, 0, 0, A.rows, A.columns);
}



template< class T >
void
sparse_base_matrix< T >::compose_h (const base_vector< T > &v, const sparse_base_matrix< T > &A)
{
	//
	//       Task: C.compose_h(v, A);
	// =  > compose the submatrix A and the vector v to matrix C
	// =  > C = (v A)
	// Conditions: A.columns <= columns and
	//             A.rows <= rows and
	//             v.size <= rows
	//    Version: 2.0
	//

	debug_handler_l(DMESSAGE,
			"compose_h(const base_vector< T > &, "
			"const sparse_base_matrix< T > &)", DVALUE + 3);

	if (A.rows > this->rows || A.columns > this->columns || v.size() > this->rows)
		precondition_error_handler(A.rows, "A.rows", "A.rows <= rows",
				    A.columns, "A.columns", "A.columns <= columns",
				    this->rows, "rows", "",
				    this->columns, "columns", "",
				    v.size(), "v.size", "v.size <= rows",
				    "void sparse_base_matrix< T >::"
				    "compose_h(const base_vector< T > &v, "
				    "const sparse_base_matrix< T > &A)",
				    DMESSAGE, EMESSAGE[4]);

	this->S_base_modul.insert_column_at(*this, 0, v.get_data_address(), v.size());
	this->SS_base_modul.insert_at(*this, 0, this->columns - A.columns, A, 0, 0, A.rows, A.columns);
}



template< class T >
void
sparse_base_matrix< T >::compose_h (const sparse_base_matrix< T > &A, const T *v)
{
	//
	//       Task: C.compose_h(A, v);
	// =  > compose the submatrix A and vector v to matrix C
	// =  > C = (A v)
	// Conditions: A.columns <= columns and
	//             A.rows <= rows and
	//             v != NULL
	//    Version: 2.0
	//

	debug_handler_l(DMESSAGE,
			"compose_h(const sparse_base_matrix< T > &, const T *)",
			DVALUE + 3);

	if (A.rows > this->rows || A.columns > this->columns || v == NULL)
		precondition_error_handler(A.rows, "A.rows", "A.rows <= rows",
				    A.columns, "A.columns", "A.columns <= columns",
				    this->rows, "rows", "",
				    this->columns, "columns", "",
				    PRT, "v", "v != NULL",
				    "void sparse_base_matrix< T >::"
				    "compose_h(const sparse_base_matrix< T > &A, const T *v)",
				    DMESSAGE, EMESSAGE[4]);

	this->SS_base_modul.insert_at(*this, 0, 0, A, 0, 0, A.rows, A.columns);
	this->S_base_modul.insert_column_at(*this, this->columns - 1, v, this->rows);
}



template< class T >
void
sparse_base_matrix< T >::compose_h (const sparse_base_matrix< T > &A, const base_vector< T > &v)
{
	//
	//       Task: C.compose_h(A, v);
	// =  > compose the submatrix A and vector v to matrix C
	// =  > C = (A v)
	// Conditions: A.columns <= columns and
	//             A.rows <= rows and
	//             v.size <= rows
	//    Version: 2.0
	//

	debug_handler_l(DMESSAGE,
			"compose_h(const sparse_base_matrix< T > &, base_vector< T > &)",
			DVALUE + 3);

	if (A.rows > this->rows || A.columns > this->columns || v.size() > this->rows)
		precondition_error_handler(A.rows, "A.rows", "A.rows <= rows",
				    A.columns, "A.columns", "A.columns <= columns",
				    this->rows, "rows", "",
				    this->columns, "columns", "",
				    v.size(), "v.size", "v.size <= rows",
				    "void sparse_base_matrix< T >::"
				    "compose_h(const sparse_base_matrix< T > &A, "
				    "const base_vector< T > &v)",
				    DMESSAGE, EMESSAGE[4]);

	this->SS_base_modul.insert_at(*this, 0, 0, A, 0, 0, A.rows, A.columns);
	this->S_base_modul.insert_column_at(*this, this->columns - 1, v.get_data_address(), v.size());
}



template< class T >
void
sparse_base_matrix< T >::compose_v (const sparse_base_matrix< T > &A, const sparse_base_matrix< T > &B)
{
	//
	//       Task: C.compose_v(A, B);
	// =  > compose the submatrices A, B to C
	//                    (A)
	// =  > C = (B)
	// Conditions: A.columns, B.columns <= columns and
	//             A.rows, B.rows <= rows
	//    Version: 2.0
	//

	debug_handler_l(DMESSAGE,
			"compose_v(const sparse_base_matrix< T > &, "
			"const sparse_base_matrix< T > &)",
			DVALUE + 3);

	if (A.rows > this->rows || A.columns > this->columns ||
	    B.rows > this->rows || B.columns > this->columns)
		precondition_error_handler(A.rows, "A.rows", "A.rows <= rows",
				    A.columns, "A.columns", "A.columns <= columns",
				    B.rows, "B.rows", "B.rows <= rows",
				    B.columns, "B.columns", "B.columns <= columns",
				    this->rows, "rows", "",
				    this->columns, "columns", "",
				    "void sparse_base_matrix< T >::"
				    "compose_v(const sparse_base_matrix< T > &A, "
				    "const sparse_base_matrix< T > &B)",
				    DMESSAGE, EMESSAGE[4]);

	this->SS_base_modul.insert_at(*this, 0, 0, A, 0, 0, A.rows, A.columns);
	this->SS_base_modul.insert_at(*this, this->rows - B.rows, 0, B, 0, 0, B.rows, B.columns);
}



template< class T >
void
sparse_base_matrix< T >::compose_v (const T *v, const sparse_base_matrix< T > &A)
{
	//
	//       Task: C.compose_v(v, A);
	// =  > compose the vector v and submatrix A to matrix C
	//                    (v)
	// =  > C = (A)
	// Conditions: A.columns <= columns and
	//             A.rows <= rows and
	//             v != NULL
	//    Version: 2.0
	//

	debug_handler_l(DMESSAGE,
			"compose_v(const T *, const sparse_base_matrix< T > &)",
			DVALUE + 3);

	if (A.rows > this->rows || A.columns > this->columns || v == NULL)
		precondition_error_handler(A.rows, "A.rows", "A.rows <= rows",
				    A.columns, "A.columns", "A.columns <= columns",
				    this->rows, "rows", "",
				    this->columns, "columns", "",
				    PRT, "v", "v != NULL",
				    "void sparse_base_matrix< T >::"
				    "compose_v(const T *v, const sparse_base_matrix< T > &A)",
				    DMESSAGE, EMESSAGE[4]);

	this->S_base_modul.insert_row_at(*this, 0, v, this->columns);
	this->SS_base_modul.insert_at(*this, this->rows - A.rows, 0, A, 0, 0, A.rows, A.columns);
}



template< class T >
void
sparse_base_matrix< T >::compose_v (const base_vector< T > &v, const sparse_base_matrix< T > &A)
{
	//
	//       Task: C.compose_v(v, A);
	// =  > compose the vector v and submatrix A to matrix C
	//                    (v)
	// =  > C = (A)
	// Conditions: A.columns <= columns and
	//             A.rows <= rows and
	//             v.size <= columns
	//    Version: 2.0
	//

	debug_handler_l(DMESSAGE,
			"compose_v(const base_vector< T > &, "
			"const sparse_base_matrix< T > &)",
			DVALUE + 3);

	if (A.rows > this->rows || A.columns > this->columns || v.size() > this->columns)
		precondition_error_handler(A.rows, "A.rows", "A.rows <= rows",
				    A.columns, "A.columns", "A.columns <= columns",
				    this->rows, "rows", "",
				    this->columns, "columns", "",
				    v.size(), "v.size", "v.size <= coluns",
				    "void sparse_base_matrix< T >::"
				    "compose_v(const base_vector< T > &v, "
				    "const sparse_base_matrix< T > &A)",
				    DMESSAGE, EMESSAGE[4]);

	this->S_base_modul.insert_row_at(*this, 0, v.get_data_address(), v.size());
	this->SS_base_modul.insert_at(*this, this->rows - A.rows, 0, A, 0, 0, A.rows, A.columns);
}



template< class T >
void
sparse_base_matrix< T >::compose_v (const sparse_base_matrix< T > &A, const T *v)
{
	//
	//       Task: C.compose_v(A, v);
	// =  > compose the submatrix A and the vector v to matrix C
	//                    (A)
	// =  > C = (v)
	// Conditions: A.columns <= columns and
	//             A.rows <= rows and
	//             v != NULL
	//    Version: 2.0
	//

	debug_handler_l(DMESSAGE,
			"compose_v(const sparse_base_matrix< T > &, const T *)",
			DVALUE + 3);

	if (A.rows > this->rows || A.columns > this->columns || v == NULL)
		precondition_error_handler(A.rows, "A.rows", "A.rows <= rows",
				    A.columns, "A.columns", "A.columns <= columns",
				    this->rows, "rows", "",
				    this->columns, "columns", "",
				    PRT, "v", "v != NULL",
				    "void sparse_base_matrix< T >::"
				    "compose_v(const sparse_base_matrix< T > &A, const T *v)",
				    DMESSAGE, EMESSAGE[4]);

	this->SS_base_modul.insert_at(*this, 0, 0, A, 0, 0, A.rows, A.columns);
	this->S_base_modul.insert_row_at(*this, this->rows - 1, v, this->columns);
}



template< class T >
void
sparse_base_matrix< T >::compose_v (const sparse_base_matrix< T > &A, const base_vector< T > &v)
{
	//
	//       Task: C.compose_v(A, v);
	// =  > compose the submatrix A and the vector v to matrix C
	//                    (A)
	// =  > C = (v)
	// Conditions: A.columns <= columns and
	//             A.rows <= rows and
	//             v.size <= columns
	//    Version: 2.0
	//

	debug_handler_l(DMESSAGE,
			"compose_v(const sparse_base_matrix< T > &, "
			"const base_vector< T > &)",
			DVALUE + 3);

	if (A.rows > this->rows || A.columns > this->columns || v.size() > this->columns)
		precondition_error_handler(A.rows, "A.rows", "A.rows <= rows",
				    A.columns, "A.columns", "A.columns <= columns",
				    this->rows, "rows", "",
				    this->columns, "columns", "",
				    v.size(), "v.size", "v.size <= columns",
				    "void sparse_base_matrix< T >::"
				    "compose_v(const sparse_base_matrix< T > &A, "
				    "const base_vector< T > &v)",
				    DMESSAGE, EMESSAGE[4]);

	this->SS_base_modul.insert_at(*this, 0, 0, A, 0, 0, A.rows, A.columns);
	this->S_base_modul.insert_row_at(*this, this->rows - 1, v.get_data_address(), v.size());
}



//
// exchange functions / swap functions
//

template< class T >
void
sparse_base_matrix< T >::swap (sparse_base_matrix< T > &B)
{
	//
	//    Task: A.swap(A) exchanges matrix A with matrix B.
	// Version: 2.0
	//

	debug_handler_l(DMESSAGE, "swap(sparse_base_matrix< T > &)", DVALUE + 4);

	// swap no_of_rows
	LiDIA::swap(this->rows, B.rows);

	// swap no_of_columns
	LiDIA::swap(this->columns, B.columns);

	// swap no_of_sparse_rows
	LiDIA::swap(this->sparse_rows, B.sparse_rows);

	// swap no_of_sparse_columns
	LiDIA::swap(this->sparse_columns, B.sparse_columns);

	// swap allocated
	register lidia_size_t *tmp2 = this->allocated;
	this->allocated = B.allocated;
	B.allocated = tmp2;

	// swap value_counter
	tmp2 = this->value_counter;
	this->value_counter = B.value_counter;
	B.value_counter = tmp2;

	// swap values
	register T **tmp = this->value;
	this->value = B.value;
	B.value = tmp;

	// swap values
	register lidia_size_t **tmp1 = this->index;
	this->index = B.index;
	B.index = tmp1;

	// swap bitfield
	matrix_flags TMP = this->bitfield;
	this->bitfield = B.bitfield;
	B.bitfield = TMP;
}



template< class T >
void
sparse_base_matrix< T >::swap_columns (lidia_size_t i, lidia_size_t j)
{
	//
	//       Task: A.swap_columns(i, j) exchanges column i with column j
	//             in matrix A.
	// Conditions: 0 <= i, j < A.columns
	//    Version: 2.0
	//

	debug_handler_l(DMESSAGE, "swap_columns(lidia_size_t , lidia_size_t)", DVALUE + 4);

	if (i< 0 || i >= this->columns || j< 0 || j >= this->columns)
		precondition_error_handler(i, "i", "0 <= i < columns",
				    j, "j", "0 <= j < columns",
				    this->columns, "columns", "",
				    "void sparse_base_matrix< T >::"
				    "swap_columns(lidia_size_t i, lidia_size_t j)",
				    DMESSAGE, EMESSAGE[1]);

	this->SS_base_modul.swap_columns(*this, i, j);
}



template< class T >
void
sparse_base_matrix< T >::swap_rows (lidia_size_t i, lidia_size_t j)
{
	//
	//       Task: A.swap_rows(i, j) exchanges row i with row j
	//             in matrix A.
	// Conditions: 0 <= i, j < A.rows
	//    Version: 2.0
	//

	debug_handler_l(DMESSAGE, "swap_rows(lidia_size_t , lidia_size_t)", DVALUE + 4);

	if (i< 0 || i >= this->rows || j< 0 || j >= this->rows)
		precondition_error_handler(i, "i", "0 <= i < rows",
				    j, "j", "0 <= j < rows",
				    this->rows, "rows", "",
				    "void sparse_base_matrix< T >::"
				    "swap_rows(lidia_size_t i, lidia_size_t j)",
				    DMESSAGE, EMESSAGE[1]);

	this->SS_base_modul.swap_rows(*this, i, j);
}



//
// structur functions
//

template< class T >
void
sparse_base_matrix< T >::set_no_of_rows (lidia_size_t r)
{
	//
	//       Task: A.set_no_of_rows(r) sets number of rows of matrix A to r.
	// Conditions: r >= 0
	//    Version: 2.0
	//

	debug_handler_l(DMESSAGE, "set_no_of_rows(lidia_size_t)", DVALUE + 5);

	if (r < 0)
		precondition_error_handler(r, "r", "r >= 0",
				    "void sparse_base_matrix< T >::"
				    "set_no_of_rows(lidia_size_t r)",
				    DMESSAGE, EMESSAGE[1]);

	this->S_base_modul.set_no_of_rows(*this, r);
}



template< class T >
void
sparse_base_matrix< T >::set_no_of_columns (lidia_size_t c)
{
	//
	//       Task: A.set_no_of_columns(c) sets number of columns of matrix A to c.
	// Conditions: 0 <= c
	//    Version: 2.0
	//

	debug_handler_l(DMESSAGE, "set_no_of_columns(lidia_size_t)", DVALUE + 5);

	if (c < 0)
		precondition_error_handler(c, "c", "c >= 0",
				    "void sparse_base_matrix< T >::"
				    "set_no_of_columns(lidia_size_t c)",
				    DMESSAGE, EMESSAGE[1]);

	this->S_base_modul.set_no_of_columns(*this, c);
}



template< class T >
void
sparse_base_matrix< T >::resize (lidia_size_t r, lidia_size_t c)
{
	//
	//       Task: A.resize(r, c)
	// =  > sets number of columns of matrix A to c and
	//                number of rows of matrix A to r.
	// Conditions: 0 <= r, c
	//    Version: 2.0
	//

	debug_handler_l(DMESSAGE,
			"resize(lidia_size_t, lidia_size_t)",
			DVALUE + 5);

	if (c < 0 || r < 0)
		precondition_error_handler(r, "r", "r >= 0",
				    c, "c", "c >= 0",
				    "void sparse_base_matrix< T >::"
				    "resize(lidia_size_t r, lidia_size_t c)",
				    DMESSAGE, EMESSAGE[1]);

	this->SS_base_modul.resize(*this, r, c);
}



template< class T >
void
sparse_base_matrix< T >::kill ()
{
	//
	//    Task: A.kill()
	// =  > sets number of columns of matrix A to 1 and
	//          number of rows of matrix A to 0.
	// Version: 2.0
	//

	debug_handler_l(DMESSAGE, "kill()", DVALUE + 5);

	this->SS_base_modul.kill(*this);
}



//
// assignments
//

template< class T >
void
sparse_base_matrix< T >::assign (const sparse_base_matrix< T > &M)
{
	//
	//    Task: A.assign(B);
	// =  > A.value[x][y] = B.value[x][y],
	//             x = 0, ..., A.rows-1, y = 0, ..., A.columns-1
	// =  > A.rows = B.rows and A.columns = B.columns
	// Version: 2.0
	//

	debug_handler_l(DMESSAGE, "assign(const sparse_base_matrix< T > &)",
			DVALUE + 6);

	this->S_base_modul.assign(*this, M);
}



//
// diagonal function
//

template< class T >
void
sparse_base_matrix< T >::diag (const T &a, const T &b)
{
	//
	//    Task: A.diag(a, b);
	// =  > A.value[i][i] = a, i = 0, ..., min(columns, rows)
	// =  > A.value[i][j] = b, i = 0, ..., rows and j = 0, ..., columns
	//             and i != j
	// Version: 2.0
	//

	debug_handler_l(DMESSAGE, "diag(const T &, const T &)", DVALUE + 7);

	this->SS_base_modul.diag(*this, a, b);
}



//
// transpose function
//

template< class T >
sparse_base_matrix< T > sparse_base_matrix< T >::trans () const
{
	//
	//    Task: AT = A.trans();
	// =  > AT.value[x][y] = A.value[y][x],
	//             x = 0, ..., A.columns-1, y = 0, ..., A.rows-1
	// Version: 2.0
	//

	debug_handler_l(DMESSAGE, "trans()", DVALUE + 7);

	sparse_base_matrix< T > TRANS(this->columns, this->rows);

	this->SS_base_modul.trans(TRANS, *this);
	return TRANS;
}



template< class T >
void
sparse_base_matrix< T >::trans (const sparse_base_matrix< T > &A)
{
	//
	//    Task: AT.trans(A);
	// =  > AT.value[x][y] = A.value[y][x],
	//             x = 0, ..., A.columns-1, y = 0, ..., A.rows-1
	// =  > AT.rows = A.columns and AT.columns = A.rows
	// Version: 2.0
	//

	debug_handler_l(DMESSAGE, "trans(const sparse_base_matrix< T > &)",
			DVALUE + 7);

	this->SS_base_modul.trans(*this, A);
}



//
// stream handling - LIDIA
//

template< class T >
void
sparse_base_matrix< T >::write_to_beauty (std::ostream &out) const
{
	//
	//    Task: A.write_to_stream(out);
	// =  > writes matrix A to stream out in beauty - format.
	// Version: 2.0
	//

	debug_handler_l(DMESSAGE, "write_to_stream(std::ostream &)", DVALUE + 8);

	this->S_base_modul.write_to_beauty(*this, out);
}



template< class T >
void
sparse_base_matrix< T >::write_to_stream (std::ostream &out) const
{
	//
	//    Task: A.write_to_stream(out);
	// =  > writes matrix A to stream out in LiDIA - format.
	// Version: 2.0
	//

	debug_handler_l(DMESSAGE, "write_to_stream(std::ostream &)", DVALUE + 8);

	this->S_base_modul.write_to_stream(*this, out);
}



template< class T >
void
sparse_base_matrix< T >::read_from_stream (std::istream &in)
{
	//
	//    Task: A.read_from_stream(in);
	// =  > reads matrix A in LiDIA - format from stream in .
	// Version: 2.0
	//

	debug_handler_l(DMESSAGE, "read_from_stream(std::istream &)", DVALUE + 8);

	this->S_base_modul.read_from_stream(*this, in);
}



//
// stream handling - MATHEMATICA
//

template< class T >
void
sparse_base_matrix< T >::write_to_mathematica (std::ostream &out) const
{
	//
	//    Task: A.write_to_mathematica(out);
	// =  > writes matrix A to stream out in mathematica format
	// Version: 2.0
	//

	debug_handler_l(DMESSAGE, "write_to_mathematica(std::ostream &)", DVALUE + 8);

	this->S_base_modul.write_to_mathematica(*this, out);
}



template< class T >
void
sparse_base_matrix< T >::read_from_mathematica (std::istream &in)
{
	//
	//    Task: A.read_from_mathematica(dz);
	// =  > reads matrix A from stream dz in mathematica format
	// Version: 2.0
	//

	debug_handler_l(DMESSAGE, "read_from_mathematica(std::istream &)", DVALUE + 8);

	this->S_base_modul.read_from_mathematica(*this, in);
}



//
// stream handling - MAPLE
//

template< class T >
void
sparse_base_matrix< T >::write_to_maple (std::ostream &out) const
{
	//
	//    Task: A.write_to_maple(out);
	// =  > writes matrix A to stream out in maple format
	// Version: 2.0
	//

	debug_handler_l(DMESSAGE, "write_to_maple(std::ostream &)", DVALUE + 8);

	this->S_base_modul.write_to_maple(*this, out);
}



template< class T >
void
sparse_base_matrix< T >::read_from_maple (std::istream &in)
{
	//
	//    Task: A.read_from_maple(dz);
	// =  > reads matrix A from stream dz in maple format
	// Version: 2.0
	//

	debug_handler_l(DMESSAGE, "read_from_maple(std::istream &)", DVALUE + 8);

	this->S_base_modul.read_from_maple(*this, in);
}



//
// stream handling - PARI
//

template< class T >
void
sparse_base_matrix< T >::write_to_gp (std::ostream &out) const
{
	//
	//    Task: A.write_to_gp(out);
	// =  > writes matrix A to stream out in PARI format
	// Version: 2.0
	//

	debug_handler_l(DMESSAGE, "write_to_gp(std::ostream &)", DVALUE + 8);

	this->S_base_modul.write_to_gp(*this, out);
}



template< class T >
void
sparse_base_matrix< T >::read_from_gp (std::istream &in)
{
	//
	//    Task: A.read_from_gp(dz);
	// =  > reads matrix A from stream dz in PARI format
	// Version: 2.0
	//

	debug_handler_l(DMESSAGE, "read_from_gp(std::istream &)", DVALUE + 8);

	this->S_base_modul.read_from_gp(*this, in);
}



//
// stream handling - KANT
//

template< class T >
void
sparse_base_matrix< T >::write_to_kash (std::ostream &out) const
{
	//
	//    Task: A.write_to_kash(out);
	// =  > writes matrix A to stream out in kash format
	// Version: 2.0
	//

	debug_handler_l(DMESSAGE, "write_to_kash(std::ostream &)", DVALUE + 8);

	this->S_base_modul.write_to_kash(*this, out);
}



template< class T >
void
sparse_base_matrix< T >::read_from_kash (std::istream &in)
{
	//
	//    Task: A.read_from_kash(dz);
	// =  > reads matrix A from stream dz in kash format
	// Version: 2.0
	//

	debug_handler_l(DMESSAGE, "read_from_kash(std::istream &)", DVALUE + 8);

	this->S_base_modul.read_from_kash(*this, in);
}



//
// stream handling - LaTeX
//

template< class T >
void
sparse_base_matrix< T >::write_to_latex (std::ostream &out) const
{
	//
	//    Task: A.write_to_latex(out);
	// =  > writes matrix A to stream out in latex format
	// Version: 2.0
	//

	debug_handler_l(DMESSAGE, "write_to_latex(std::ostream &)", DVALUE + 8);

	this->S_base_modul.write_to_latex(*this, out);
}



//
// stream handling - Magma
//

template< class T >
void
sparse_base_matrix< T >::write_to_magma (std::ostream &out) const
{
	//
	//    Task: A.write_to_magma(out);
	// =  > writes matrix A to stream out in magma format
	// Version: 2.0
	//

	debug_handler_l(DMESSAGE, "write_to_magma(std::ostream &)", DVALUE + 8);

	this->S_base_modul.write_to_magma(*this, out);
}



//
// boolean functions
//

template< class T >
bool
sparse_base_matrix< T >::is_column_zero (lidia_size_t c) const
{
	//
	//    Task: A.is_column_zero(c) == true <=  > A.value[x][c] == 0,
	//          x = 0, ..., rows-1
	// Version: 2.0
	//

	debug_handler_l(DMESSAGE, "is_column_zero(lidia_size_t)", DVALUE + 10);

	if (c< 0 || c >= this->columns)
		precondition_error_handler(c, "c", "0 <= c < columns",
				    this->columns , "columns", "",
				    "bool sparse_base_matrix< T >::"
				    "is_column_zero(lidia_size_t c) const",
				    DMESSAGE, EMESSAGE[3]);

	return this->S_base_modul.is_column_zero(*this, c);
}



template< class T >
bool
sparse_base_matrix< T >::is_row_zero (lidia_size_t r) const
{
	//
	//    Task: A.is_row_zero(r) == true <=  > A.value[r][x] == 0,
	//          x = 0, ..., columns-1
	// Version: 2.0
	//

	debug_handler_l(DMESSAGE, "is_row_zero(lidia_size_t)", DVALUE + 10);

	if (r< 0 || r >= this->rows)
		precondition_error_handler(r, "r", "0 <= r < rows",
				    this->rows, "rows", "",
				    "bool sparse_base_matrix< T >::"
				    "is_row_zero(lidia_size_t r) const",
				    DMESSAGE, EMESSAGE[3]);

	return this->S_base_modul.is_row_zero(*this, r);
}



template< class T >
bool
sparse_base_matrix< T >::is_matrix_zero () const
{
	//
	//    Task: A.is_matrix_zero() == true <=  > A.value[x][y] == 0,
	//          x = 0, ..., rows-1 and y = 0, ..., columns-1
	// Version: 2.0
	//

	debug_handler_l(DMESSAGE, "is_matrix_zero()", DVALUE + 10);

	return this->S_base_modul.is_matrix_zero(*this);
}



//
// internal functions
//

template< class T >
void
sparse_base_matrix< T >::change_orientation (unsigned long mode)
{
	debug_handler_l(DMESSAGE, "change_orientation(unsigned long)", DVALUE + 9);

	this->S_base_modul.change_orientation(*this, mode);
}



//
// status report
//

template< class T >
void
sparse_base_matrix< T >::status_report ()
{
	debug_handler_l(DMESSAGE, "status_report()", DVALUE + 9);

	std::cout << "configuration: "
		  << "\n\t rows = " << this->rows
		  << "\n\t columns = " << this->columns
		  << "\n\t sparse_rows = " << this->sparse_rows
		  << "\n\t sparse_columns = " << this->sparse_columns << std::endl;

	for (register lidia_size_t i = 0; i < this->rows; i++) {
		std::cout << "row " << i << ": allocated = " << this->allocated[i]
			  << " value_counter = " << this->value_counter[i] << std::endl;
	}
}



#undef DVALUE
#undef DMESSAGE
#undef EMESSAGE



#ifdef LIDIA_NAMESPACE
# ifndef IN_NAMESPACE_LIDIA
}	// end of namespace LiDIA
# endif
#endif



#endif	// LIDIA_SPARSE_BASE_MATRIX_CC_GUARD_
