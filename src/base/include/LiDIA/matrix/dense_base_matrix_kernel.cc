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


#ifndef LIDIA_DENSE_BASE_MATRIX_KERNEL_CC_GUARD_
#define LIDIA_DENSE_BASE_MATRIX_KERNEL_CC_GUARD_



#ifndef LIDIA_DENSE_BASE_MATRIX_KERNEL_H_GUARD_
# include	"LiDIA/matrix/dense_base_matrix_kernel.h"
#endif



#ifdef LIDIA_NAMESPACE
# ifndef IN_NAMESPACE_LIDIA
namespace LiDIA {
# endif
#endif



#define dense_base_matrix_kernel DBMK

//
// error defines
//


extern const char *PRT;
extern const char *matrix_error_msg[];

#define DMESSAGE "dense_base_matrix_kernel"
#define EMESSAGE matrix_error_msg

//
// constructor kernel
//

template< class T >
inline void dense_base_matrix_kernel< T >::
constructor(MR< T > &A, lidia_size_t r, lidia_size_t c) const
{
	A.rows = r;
	A.sparse_rows = r;
	A.columns = c;
	A.sparse_columns = c;

	A.index = NULL;
	A.allocated = A.value_counter = NULL;

	if (r == 0)
		A.value = NULL;
	else {
		A.value = new T *[r];
		memory_handler(A.value, DMESSAGE, "constructor(MR< T > &, "
			       "lidia_size_t, lidia_size_t) :: "
			       "Error in memory allocation (A.value)");
	}

	for (register lidia_size_t i = 0; i < r; i++) {
		A.value[i] = new T[c];
		memory_handler(A.value[i], DMESSAGE, "constructor(MR< T > &, "
			       "lidia_size_t, lidia_size_t) :: "
			       "Error in memory allocation (A.value[i])");
	}
}



template< class T >
void dense_base_matrix_kernel< T >::
constructor(MR< T > &B, lidia_size_t r, lidia_size_t c, const T **A) const
{
	B.rows = r;
	B.sparse_rows = r;
	B.columns = c;
	B.sparse_columns = c;

	B.index = NULL;
	B.allocated = NULL;
	B.value_counter = NULL;

	if (r == 0)
		B.value = NULL;
	else {
		B.value = new T *[r];
		memory_handler(B.value, DMESSAGE, "constructor(MR< T > &, "
			       "lidia_size_t, lidia_size_t, const T **) :: "
			       "Error in memory allocation (B.value)");
	}

	register T *Btmp;
	register const T *Atmp;
	for (register lidia_size_t i = 0; i < r; i++) {
		Atmp = A[i];

		Btmp = new T[c];
		memory_handler(Btmp, DMESSAGE, "constructor(MR< T > &, "
			       "lidia_size_t, lidia_size_t, const T **) :: "
			       "Error in memory allocation (Btmp)");

		for (register lidia_size_t j = 0; j < c; j++)
			Btmp[j] = Atmp[j];
		B.value[i] = Btmp;
	}
}



template< class T >
void dense_base_matrix_kernel< T >::
constructor(MR< T > &A, const MR< T > &M) const
{
	A.rows = M.rows;
	A.sparse_rows = M.sparse_rows;

	A.columns = M.columns;
	A.sparse_columns = M.sparse_columns;

	A.allocated = NULL;
	A.index = NULL;
	A.value_counter = NULL;

	if (M.rows == 0)
		A.value = NULL;
	else {
		A.value = new T *[A.rows];
		memory_handler(A.value, DMESSAGE,
			       "constructor(MR< T > &, "
			       "const MR< T > &) :: "
			       "Error in memory allocation (A.value)");
	}

	register T *Atmp, *Mtmp;
	for (register lidia_size_t i = 0; i < A.rows; i++) {
		Atmp = new T[A.columns];
		Mtmp = M.value[i];
		memory_handler(Atmp, DMESSAGE,
			       "constructor(MR< T > &, "
			       "const MR< T > &) :: "
			       "Error in memory allocation (Atmp)");

		for (register lidia_size_t j = 0; j < A.columns; j++)
			Atmp[j] = Mtmp[j];
		A.value[i] = Atmp;
	}

	A.bitfield = M.bitfield;
}



//
// destructor
//

template< class T >
inline void dense_base_matrix_kernel< T >::
destructor(MR< T > &A) const
{
	if (A.rows > 0) {
		for (A.rows -= 1; A.rows >= 0; A.rows--)
			delete[] A.value[A.rows];
		delete[] A.value;
	}
}



//
// element access kernel
//

template< class T >
inline void dense_base_matrix_kernel< T >::
sto(MR< T > &A, lidia_size_t x, lidia_size_t y, const T &e) const
{
	A.value[x][y] = e;
}



template< class T >
inline T & dense_base_matrix_kernel< T >::
member(const MR< T > &A, lidia_size_t x, lidia_size_t y) const
{
	return A.value[x][y];
}



//
// column access kernel
//

template< class T >
inline void dense_base_matrix_kernel< T >::
sto_column(MR< T > &A, const T *v, lidia_size_t l,
	   lidia_size_t j, lidia_size_t from) const
{
	for (register lidia_size_t k = from; k < from + l && k < A.rows; k++)
		A.value[k][j] = v[k - from];
}



template< class T >
inline void dense_base_matrix_kernel< T >::
get_column(const MR< T > &A, T *RES, lidia_size_t i) const
{
	for (register lidia_size_t j = A.rows - 1; j >= 0; j--)
		RES[j] = A.value[j][i];
}



//
// row access kernel
//

template< class T >
inline void dense_base_matrix_kernel< T >::
sto_row(MR< T > &A, const T *v, lidia_size_t l,
	lidia_size_t i, lidia_size_t from) const
{
	register T *Atmp = A.value[i];
	for (register lidia_size_t k = from; k < l + from && k < A.columns; k++)
		Atmp[k] = v[k-from];
}



template< class T >
inline void dense_base_matrix_kernel< T >::
get_row(const MR< T > &A, T *RES, lidia_size_t i) const
{
	register T *Atmp = A.value[i];
	for (register lidia_size_t k = 0; k < A.columns; k++)
		RES[k] = Atmp[k];
}



//
// insert and remove kernel
//

template< class T >
void dense_base_matrix_kernel< T >::
insert_columns(MR< T > &A, lidia_size_t *ind, const T **news) const
{
	A.columns += ind[0];
	T **new_value = new T *[A.rows];
	memory_handler(new_value, DMESSAGE, "insert_columns(MR< T > &, "
		       "lidia_size_t *, const T **) :: "
		       "Error in memory allocation (new_value)");

	register lidia_size_t i; // MM, for HP UX 10.20, CC
	for (i = 0; i < A.rows; i++) {
		new_value[i] = new T[A.columns];
		memory_handler(new_value[i], DMESSAGE, "insert_columns(MR< T > &, "
			       "lidia_size_t *, const T **) :: "
			       "Error in memory allocation (new_value[i])");
	}

	lidia_size_t l1 = 0, l2 = 0, l3 = 1;
	for (i = 0; i < A.columns; i++) {
		if (i == ind[l3]) {
			for (register lidia_size_t j = 0; j < A.rows; j++)
				new_value[j][i] = news[j][l2];
			l2++;
			l3++;
		}
		else {
			for (register lidia_size_t j = 0; j < A.rows; j++)
				new_value[j][i] = A.value[j][l1];
			l1++;
		}
	}

	for (i = 0; i < A.rows; i++)
		delete[] A.value[i];
	delete[] A.value;
	A.value = new_value;
}



template< class T >
inline void dense_base_matrix_kernel< T >::
insert_column_at(MR< T > &M, lidia_size_t j, const T *v, lidia_size_t len) const
{
	lidia_size_t i;

	for (i = 0; i < len; i++)
		M.value[i][j] = v[i];
}



template< class T >
inline void dense_base_matrix_kernel< T >::
remove_columns(MR< T > &A, lidia_size_t *rem) const
{
	lidia_size_t i, j, l = 0, l1 = 1;
	for (j = 0; j < A.columns; j++)
		if (rem[l1] != j)
			for (i = 0; i < A.rows; i++) {
				LiDIA::swap(A.value[i][j], A.value[i][l]);
				l++;
			}
		else
			l1++;
	A.columns -= rem[0];
}



template< class T >
void dense_base_matrix_kernel< T >::
insert_rows(MR< T > &A, lidia_size_t *ind, const T **news) const
{
	A.rows += ind[0];
	T **new_value = new T*[A.rows];
	memory_handler(new_value, DMESSAGE, "insert_rows(MR< T > &, "
		       "lidia_size_t *, const T **) :: "
		       "Error in memory allocation (new_value)");

	lidia_size_t l1 = 0, l2 = 1, l3 = 0;
	for (register lidia_size_t i = 0; i <= A.rows; i++) {
		if (i == ind[l2]) {
			new_value[i] = new T[A.columns];
			memory_handler(new_value[i], DMESSAGE, "insert_rows(MR< T > &, "
				       "lidia_size_t *, const T **) :: "
				       "Error in memory allocation (new_value[i])");
			for (register lidia_size_t j = 0; j < A.columns; j++)
				new_value[i][j] = news[l3][j];
			l2++;
			l3++;
		}
		else {
			new_value[i] = A.value[l1];
			l1++;
		}
	}
	delete[] A.value;
	A.value = new_value;
}



template< class T >
inline void dense_base_matrix_kernel< T >::
insert_row_at(MR< T > &M, lidia_size_t j, const T *v, lidia_size_t len) const
{
	lidia_size_t i;

	for (i = 0; i < len; i++)
		M.value[j][i] = v[i];
}



template< class T >
inline void dense_base_matrix_kernel< T >::
remove_rows(MR< T > &A, lidia_size_t *rem) const
{
	lidia_size_t len = A.rows;
	A.rows -= rem[0];

	T **new_value = new T*[A.rows];

	register lidia_size_t i, l = 0, l1 = 1;
	for (i = 0; i <= len; i++)
		if (rem[l1] != i) {
			new_value[l] = A.value[i];
			l++;
		}
		else
			l1++;
	delete[] A.value;

	A.value = new_value;
}



//
// exchange functions / swap functions
//

template< class T >
inline void dense_base_matrix_kernel< T >::
swap_columns(MR< T > &A, lidia_size_t i, lidia_size_t j) const
{
	T *tmp;

	for (register lidia_size_t k = A.rows - 1; k >= 0; k--) {
		tmp = A.value[k];
		LiDIA::swap(tmp[i], tmp[j]);
	}
}



template< class T >
inline void dense_base_matrix_kernel< T >::
swap_rows(MR< T > &A, lidia_size_t i, lidia_size_t j) const
{
	T *tmp = A.value[i];
	A.value[i] = A.value[j];
	A.value[j] = tmp;
}



//
// assignments
//

template< class T >
inline void dense_base_matrix_kernel< T >::
assign(MR< T > &A, const MR< T > &M) const
{
	for (register lidia_size_t len1 = A.rows - 1; len1 >= 0; len1--)
		for (register lidia_size_t len = A.columns - 1; len >= 0; len--)
			A.value[len1][len] = M.value[len1][len];
}



//
// diagonal function
//

template< class T >
inline void dense_base_matrix_kernel< T >::
diag(MR< T > &A, const T &a, const T &b) const
{
	T *tmp;

	for (register lidia_size_t i = 0; i < A.rows; i++) {
		tmp = A.value[i];
		for (register lidia_size_t j = 0; j < A.columns; j++)
			tmp[j] = (i == j) ? a : b;
	}
}



//
// transpose function
//

template< class T >
void dense_base_matrix_kernel< T >::
trans(MR< T > &R, const MR< T > &A) const
{
	register lidia_size_t i, j;
	T *Atmp;

	if (A.value != R.value) {
		if (R.columns != A.rows)
			set_no_of_columns(R, A.rows);
		if (R.rows != A.columns)
			set_no_of_rows(R, A.columns);

		for (i = 0; i < A.rows; i++) {
			Atmp = A.value[i];
			for (j = 0; j < A.columns; j++)
				R.value[j][i] = Atmp[j];
		}
	}
	else {
		register lidia_size_t oldrows = R.rows, oldcolumns = R.columns;
		if (R.rows != R.columns)
			if (R.columns > R.rows) {
				set_no_of_rows(R, oldcolumns);
				for (i = 0; i < R.rows; i++)
					for (j = 0; j < i; j++)
						LiDIA::swap(R.value[i][j], R.value[j][i]);
				set_no_of_columns(R, oldrows);
			}
			else {
				set_no_of_columns(R, oldrows);
				for (i = 0; i < R.rows; i++)
					for (j = 0; j < i; j++)
						LiDIA::swap(R.value[i][j], R.value[j][i]);
				set_no_of_rows(R, oldcolumns);
			}
		else {
			for (i = 0; i < R.rows; i++)
				for (j = 0; j < i; j++)
					LiDIA::swap(R.value[i][j], R.value[j][i]);
		}
	}
}



//
// stream handling - LIDIA
//

template< class T >
void dense_base_matrix_kernel< T >::
write_to_beauty(const MR< T > &A, std::ostream &out) const
{
	register lidia_size_t i, j;
	register T *tmp;
	for (i = 0; i < A.rows; i++) {
		tmp = A.value[i];
		out << std::endl << "(";
		for (j = 0; j < A.columns; j++)
			out << tmp[j] << " ";
		out << ")";
	}
	out << std::endl << std::flush;
}



template< class T >
void dense_base_matrix_kernel< T >::
write_to_stream(const MR< T > &A, std::ostream &out) const
{
	register lidia_size_t i, j, col = A.columns - 1;
	T *tmp;

	out << A.rows << " " << A.columns << std::endl;
	for (i = 0; i < A.rows; i++) {
		tmp = A.value[i];
		for (j = 0; j < col; j++)
			out << tmp[j] << " ";
		out << tmp[j] << std::endl << std::flush;
	}
}



template< class T >
void dense_base_matrix_kernel< T >::
read_from_stream(MR< T > &A, std::istream &in) const
{
    if(in.good()) {
	lidia_size_t i, j;
	T *tmp;
	
	in >> i;
	if(i < 0) {
	    in.setstate(std::ios::failbit);
	}
	if(in.fail()) {
	    return;
	}
	if (i != A.rows)
		set_no_of_rows(A, i);

	in >> j;
	if(j < 0) {
	    in.setstate(std::ios::failbit);
	}
	if(in.fail()) {
	    return;
	}
	if (j != A.columns)
		set_no_of_columns(A, j);

	for (i = 0; i < A.rows; i++) {
		tmp = A.value[i];
		for (j = 0; j < A.columns; j++)
			in >> tmp[j];
	}
    }
}



//
// stream handling - MATHEMATICA
//

template< class T >
void dense_base_matrix_kernel< T >::
write_to_mathematica(const MR< T > &A, std::ostream &out) const
{

	register lidia_size_t i, j, l = A.columns - 1, k = A.rows - 1;
	T *tmp;

	out << "{";
	for (i = 0; i < A.rows; i++) {
		tmp = A.value[i];
		out << "{";
		for (j = 0; j < A.columns; j++) {
			out << tmp[j];
			if (j != l)
				out << ", ";
			else
				if (i != k)
					out << "}, " << std::endl;
		}
	}
	out << "}}" << std::endl << std::flush;
}



template< class T >
void dense_base_matrix_kernel< T >::
read_from_mathematica(MR< T > &A, std::istream &dz) const
{
	char c = 0;

	dz >> std::ws >> c;
	if (c != '{')
		lidia_error_handler("void MR< T >::"
				    "read_from_mathematica(std::istream &dz)",
				    DMESSAGE, EMESSAGE[5]);

	base_vector< T > buffer(0, vector_flags(vector_flags::expand));
	lidia_size_t len = 0, local_rows = 0, local_columns = 0;
	bool END = false;

	dz >> std::ws >> c;
	while (!dz.eof() && END == false) {
		if (c == '{') {
			local_rows++;
			while (c != '}') {
				dz >> buffer[len];
				len++;
				dz >> std::ws >> c;
			}
		}
		dz >> std::ws >> c;
		if (c == '}')
			END = true;
	}

	local_columns = len / local_rows;

	len = 0;
	if (A.rows != local_rows)
		set_no_of_rows(A, local_rows);
	if (A.columns != local_columns)
		set_no_of_columns(A, local_columns);

	for (register lidia_size_t i = 0; i < A.rows; i++) {
		for (register lidia_size_t j = 0; j < A.columns; j++) {
			A.value[i][j] = buffer[len];
			len++;
		}
	}
}



//
// stream handling - MAPLE
//

template< class T >
void dense_base_matrix_kernel< T >::
write_to_maple(const MR< T > &A, std::ostream &out) const
{
	register lidia_size_t i, j, l = A.columns - 1, k = A.rows - 1;
	T *tmp;

	out << "array(1 .. " << A.rows << ", 1 .. " << A.columns << ", [";

	for (i = 0; i < A.rows; i++) {
		tmp = A.value[i];
		for (j = 0; j < A.columns; j++) {
			out << "(" << i+1 << ", " << j+1 << ") = " << tmp[j];
			if (!(j == l && i == k))
				out << ", ";
		}
	}
	out << "]); " << std::endl << std::flush;
}



template< class T >
void dense_base_matrix_kernel< T >::
read_from_maple(MR< T > &A, std::istream &dz) const
{
	char c = 0;
	while (c != '(' && dz.good())
		dz >> std::ws >> c;

	if (dz.fail())
		lidia_error_handler("void MR< T >::"
				    "read_from_maple(std::istream &dz)",
				    DMESSAGE, EMESSAGE[5]);

	register lidia_size_t i, j;
	lidia_size_t startr, startc, endr, endc, x, y;
	T TMP;

	// read row dimension
	dz >> startr >> std::ws >> c >> std::ws >> c >> endr >> std::ws >> c;

	// read column dimension
	dz >> startc >> std::ws >> c >> std::ws >> c >> endc >> std::ws >> c;

	endr = endr - startr + 1;
	if (endr < 0)
		endr = -endr;
	endc = endc - startc + 1;
	if (endc < 0)
		endc = -endc;

	if (A.columns != endc)
		set_no_of_columns(A, endc);
	if (A.rows != endr)
		set_no_of_rows(A, endr);

	dz >> std::ws >> c;

	for (i = 0; i < endc; i++)
		for (j = 0; j < endr; j++) {
			dz >> std::ws >> c >> x >> std::ws >> c >> y;
			dz >> std::ws >> c >> std::ws >> c >> std::ws >> TMP >> std::ws >> c;
			A.value[x - startr][y - startc] = TMP;
		}
	dz >> std::ws >> c >> std::ws >> c;
}



//
// stream handling - PARI
//

template< class T >
void dense_base_matrix_kernel< T >::
write_to_gp(const MR< T > &A, std::ostream &out) const
{
	register lidia_size_t i, j, l = A.columns - 1, k = A.rows - 1;
	T *tmp;

	out << "[";
	for (i = 0; i < A.rows; i++) {
		tmp = A.value[i];
		for (j = 0; j < A.columns; j++) {
			out << tmp[j];
			if (j != l)
				out << ", ";
		}
		if (i != k)
			out << "; ";
	}
	out << "]" << std::flush;
}



template< class T >
void dense_base_matrix_kernel< T >::
read_from_gp(MR< T > &A, std::istream &dz) const
{
	char c = 0;

	dz >> std::ws >> c;
	if (c != '[')
		lidia_error_handler("void MR< T >::"
				    "read_from_gp(std::istream &dz)",
				    DMESSAGE, EMESSAGE[5]);

	base_vector< T > buffer(0, vector_flags(vector_flags::expand));
	lidia_size_t len = 0, local_rows = 0, local_columns = 0;
	bool END = false;

	while (!dz.eof() && END == false) {
		local_rows++;
		do {
			dz >> buffer[len];
			len++;
			dz >> std::ws >> c;
		}
		while (c != ';' && c != ']');

		if (c == ']')
			END = true;
	}

	local_columns = len / local_rows;

	len = 0;
	if (A.rows != local_rows)
		set_no_of_rows(A, local_rows);
	if (A.columns != local_columns)
		set_no_of_columns(A, local_columns);

	for (register lidia_size_t i = 0; i < A.rows; i++) {
		for (register lidia_size_t j = 0; j < A.columns; j++) {
			A.value[i][j] = buffer[len];
			len++;
		}
	}
}



//
// stream handling - KANT
//

template< class T >
void dense_base_matrix_kernel< T >::
write_to_kash(const MR< T > &A, std::ostream &out) const
{
	register lidia_size_t i, j, l = A.columns - 1, k = A.rows - 1;
	T *tmp;

	out << "LIDIA: = Mat(Z, [";
	for (i = 0; i < A.rows; i++) {
		tmp = A.value[i];
		out << "[";
		for (j = 0; j < A.columns; j++) {
			out << tmp[j];
			if (j != l)
				out << ", ";
			else
				if (i != k)
					out << "], " << std::endl;
		}
	}
	out << "]]); " << std::endl << std::flush;
}



template< class T >
void dense_base_matrix_kernel< T >::
read_from_kash(MR< T > &A, std::istream &dz) const
{
	char c = 0;
	do {
		dz >> std::ws >> c;
	} while (c != 'M' || dz.eof());
	if (c != 'M')
		lidia_error_handler("void MR< T >::"
				    "read_from_kash(std::istream &dz)",
				    DMESSAGE, EMESSAGE[5]);

	dz >> std::ws >> c >> std::ws >> c >> std::ws >> c >> std::ws >> c >> std::ws >> c >> std::ws >> c;

	base_vector< T > buffer(0, 0, vector_flags(vector_flags::expand));
	lidia_size_t len = 0, local_rows = 0, local_columns = 0;
	bool END = false;

	dz >> std::ws >> c;
	while (!dz.eof() && END == false) {
		if (c == '[') {
			local_rows++;
			while (c != ']') {
				dz >> buffer[len];
				len++;
				dz >> std::ws >> c;
			}
		}
		dz >> std::ws >> c;
		if (c == ']')
			END = true;
	}

	local_columns = len / local_rows;

	len = 0;
	if (A.rows != local_rows)
		set_no_of_rows(A, local_rows);
	if (A.columns != local_columns)
		set_no_of_columns(A, local_columns);

	for (register lidia_size_t i = 0; i < A.rows; i++) {
		for (register lidia_size_t j = 0; j < A.columns; j++) {
			A.value[i][j] = buffer[len];
			len++;
		}
	}
	dz >> std::ws >> c >> std::ws >> c;
}



//
// stream handling - LATEX
//

template< class T >
void dense_base_matrix_kernel< T >::
write_to_latex(const MR< T > &A, std::ostream &out) const
{
	out << "1" << A.value[0][0] << std::endl;
}



//
// stream handling - MAGMA
//

template< class T >
void dense_base_matrix_kernel< T >::
write_to_magma(const MR< T > &A, std::ostream &out) const
{
	lidia_size_t i, j;
	out << "L : = RMatrixSpace(Integers(), " << A.rows << ", " << A.columns << ") ! [" << std::endl;
	for (i = 0; i < A.rows; i++) {
		for (j = 0; j < A.columns; j++) {
			out << A.value[i][j] << std::flush;
			if (j < A.columns - 1)
				out << ", " << std::flush;
			else if (j == A.columns - 1 && i < A.rows - 1)
				out << ", " << std::endl;
			else if (j == A.columns - 1 && i == A.rows - 1)
				out << "]; " << std::endl;
		}
	}
}



//
// structure functions krenel
//

template< class T >
void dense_base_matrix_kernel< T >::
set_no_of_rows(MR< T > &A, lidia_size_t r) const
{
	if (r == A.rows)
		return;

	lidia_size_t i;
	if (r == 0) {
		for (i = 0; i < A.rows; i++)
			delete[] A.value[i];
		delete[] A.value;
		A.value = NULL;
		A.rows = A.columns = 0;
		return;
	}

	T **tmp = A.value;
	if (r < A.rows) {
		for (i = r; i < A.rows; i++)
			delete[] A.value[i];
		A.value = new T *[r];
		memory_handler(A.value, DMESSAGE, "set_no_of_rows :: "
			       "Error in memory allocation (A.value)");
		for (i = 0; i < r; i++)
			A.value[i] = tmp[i];
	}
	else {
		A.value = new T *[r];
		memory_handler(A.value, DMESSAGE, "set_no_of_rows :: "
			       "Error in memory allocation (A.value)");
		for (i = 0; i < A.rows; i++)
			A.value[i] = tmp[i];
		for (i = A.rows; i < r; i++) {
			A.value[i] = new T[A.columns];
			memory_handler(A.value[i], DMESSAGE, "set_no_of_rows :: "
				       "Error in memory allocation (A.value[i])");
		}
	}
	A.rows = r;
	delete[] tmp;
}



template< class T >
void dense_base_matrix_kernel< T >::
set_no_of_columns(MR< T > &A, lidia_size_t c) const
{
	if (A.columns == c)
		return;

	if (c == 0) {
		for (register lidia_size_t i = 0; i < A.rows; i++)
			delete[] A.value[i];
		delete[] A.value;
		A.value = NULL;
		A.rows = 0;
		A.columns = 0;
		return;
	}

	T *tmp, *tmp1;

	if (c < A.columns) {
		for (register lidia_size_t i = 0; i < A.rows; i++) {
			tmp1 = A.value[i];
			tmp = new T[c];
			memory_handler(tmp, DMESSAGE, "set_no_of_columns :: "
				       "Error in memory allocation (tmp)");
			copy_data(tmp, tmp1, c);
			A.value[i] = tmp;
			delete[] tmp1;
		}
	}
	else {
		for (register lidia_size_t i = 0; i < A.rows; i++) {
			tmp1 = A.value[i];
			tmp = new T[c];
			memory_handler(tmp, DMESSAGE, "set_no_of_columns :: "
				       "Error in memory allocation (tmp)");
			for (register lidia_size_t j = 0; j < A.columns; j++)
				tmp[j] = tmp1[j];
			A.value[i] = tmp;
			delete[] tmp1;
		}
	}
	A.columns = c;
}



template< class T >
inline void dense_base_matrix_kernel< T >::
kill(MR< T > &A) const
{
	for (register lidia_size_t i = 0; i < A.rows; i++)
		delete[] A.value[i];
	delete[] A.value;
	A.value = NULL;
	A.columns = 0;
	A.rows = 0;
}



//
// boolean functions
//

template< class T >
inline bool dense_base_matrix_kernel< T >::
is_column_zero(const MR< T > &RES, lidia_size_t c) const
{
	register lidia_size_t i;

	for (i = 0; i < RES.rows; i++) {
		if (RES.value[i][c] != RES.Zero)
			return false;
	}
	return true;
}



template< class T >
inline bool dense_base_matrix_kernel< T >::
is_row_zero(const MR< T > &RES, lidia_size_t r) const
{
	register lidia_size_t i;
	register T *tmp = RES.value[r];

	for (i = 0; i < RES.columns; i++) {
		if (tmp[i] != RES.Zero)
			return false;
	}
	return true;
}



template< class T >
inline bool dense_base_matrix_kernel< T >::
is_matrix_zero(const MR< T > &RES) const
{
	register lidia_size_t i, j;
	register T *tmp;

	for (i = 0; i < RES.rows; i++) {
		tmp = RES.value[i];
		for (j = 0; j < RES.columns; j++)
			if (tmp[j] != RES.Zero)
				return false;
	}
	return true;
}



//
// change orientation
//

template< class T >
void dense_base_matrix_kernel< T >::
change_orientation(MR< T > &A, unsigned long mode) const
{
	lidia_size_t index1 = 0, index2 = 0;
	if (mode != A.bitfield.get_orientation()) {
		if (mode == matrix_flags::column_oriented) {
			index1 = A.columns;
			index2 = A.rows;
			A.bitfield.set_orientation(matrix_flags::column_oriented);
		}
		if (mode == matrix_flags::row_oriented) {
			index1 = A.rows;
			index2 = A.columns;
			A.bitfield.set_orientation(matrix_flags::row_oriented);
		}

		T **oldvalue = A.value;
		A.value = new T*[index1];
		memory_handler(A.value, DMESSAGE,
			       "change_orientation(MR< T > &, unsigned long) ::"
			       "Error in memory allocation (A.value)");

		for (register lidia_size_t i = 0; i < index1; i++) {
			A.value[i] = new T[A.value_counter[i]];
			memory_handler(A.value[i], DMESSAGE,
				       "change_orientation(MR< T > &, unsigned long) ::"
				       "Error in memory allocation (A.value[i])");

			for (register lidia_size_t j = 0; j < index2; j++) {
				A.value[i][j] = oldvalue[j][i];
			}
			delete[] oldvalue[i];
		}
		delete[] oldvalue;
	}
}



#undef DMESSAGE
#undef EMESSAGE


#undef dense_base_matrix_kernel



#ifdef LIDIA_NAMESPACE
# ifndef IN_NAMESPACE_LIDIA
}	// end of namespace LiDIA
# endif
#endif



#endif	// LIDIA_DENSE_BASE_MATRIX_KERNEL_CC_GUARD_
