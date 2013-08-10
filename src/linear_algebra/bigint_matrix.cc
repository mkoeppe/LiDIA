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


#ifdef HAVE_CONFIG_H
# include	"config.h"
#endif
#include	"LiDIA/timer.h"
#include	"LiDIA/error.h"
#include	"LiDIA/info.h"
#include	"LiDIA/bigint_matrix.h"
#include	"LiDIA/base_power_product.h"
#include	"LiDIA/matrix/dense_bigint_matrix_modules.h"
#include	"LiDIA/matrix/sparse_bigint_matrix_modules.h"
#include	"LiDIA/matrix/hnf_conf.h"
#include	"LiDIA/matrix/hnf_kernel.h"
#include	"LiDIA/matrix/crt_and_prime_handling.h"
#include	"LiDIA/matrix/dense_fp_matrix_kernel.h"
#include	"LiDIA/matrix/sparse_fp_matrix_kernel.h"
#include	"LiDIA/matrix/sparse_fp_matrix_algorithms.h"
#include	"LiDIA/matrix/bigint_matrix_algorithms.h"
#include	"LiDIA/matrix/modular_arithmetic.h"



#ifdef LIDIA_NAMESPACE
namespace LiDIA {
#endif



//
// debug defines / error defines
//

extern const char *PRT;
extern const char *matrix_error_msg[];

#define DVALUE LDBL_MATRIX           // Debug value
#define DMESSAGE "bigint_matrix"     // Debug message
#define EMESSAGE matrix_error_msg    // Error message

//
// debug level
//
//   0 : remainder
//   1 : divide
//   2 : norms and bounds
//   3 : randomize
//   4 : Linear algebra PART 1
//   5 : Linear algebra PART 2
//

#define DRMKex DRMK< bigint >
#define SRMKex SRMK< bigint >

//
// constructor
//

matrix< bigint >::
matrix(const base_matrix< long > &M)
{
	if (M.bitfield.get_representation() == matrix_flags::dense_representation) {
		rows = M.rows;
		sparse_rows = M.sparse_rows;

		columns = M.columns;
		sparse_columns = M.sparse_columns;

		allocated = NULL;
		index = NULL;
		value_counter = NULL;

		if (M.rows == 0)
			value = NULL;
		else {
			value = new bigint *[rows];
			memory_handler(value, DMESSAGE,
				       "constructor(MR< T > &, "
				       "const MR< T > &) :: "
				       "Error in memory allocation (value)");
		}

		register bigint *Atmp;
		register long *Mtmp;
		for (register lidia_size_t i = 0; i < rows; i++) {
			Atmp = new bigint[columns];
			Mtmp = M.value[i];
			memory_handler(Atmp, DMESSAGE,
				       "constructor(MR< T > &, "
				       "const MR< T > &) :: "
				       "Error in memory allocation (Atmp)");

			for (register lidia_size_t j = 0; j < columns; j++)
				Atmp[j] = bigint(Mtmp[j]);
			value[i] = Atmp;
		}
	}
	else {
		rows = M.rows;
		columns = M.columns;

		sparse_rows = M.sparse_rows;
		sparse_columns = M.sparse_columns;

		value = new bigint *[rows];
		memory_handler(value, DMESSAGE,
			       "constructor((MR< T > &, const MR< T > &) ::"
			       "Error in memory allocation (value)");

		index = new lidia_size_t *[rows];
		memory_handler(index, DMESSAGE,
			       "constructor((MR< T > &, const MR< T > &) ::"
			       "Error in memory allocation (index)");

		value_counter = new lidia_size_t[rows];
		memory_handler(value_counter, DMESSAGE,
			       "constructor((MR< T > &, const MR< T > &) ::"
			       "Error in memory allocation (value_counter)");

		allocated = new lidia_size_t[rows];
		memory_handler(allocated, DMESSAGE,
			       "constructor((MR< T > &, const MR< T > &) ::"
			       "Error in memory allocation (allocated)");

		for (register lidia_size_t i = 0; i < rows; i++) {
			register lidia_size_t size = M.allocated[i];
			if (size) {
				index[i] = new lidia_size_t[size];
				memory_handler(index[i], DMESSAGE,
					       "constructor((MR< T > &, const MR< T > &) ::"
					       "Error in memory allocation (index[i])");
				value[i] = new bigint[size];
				memory_handler(value[i], DMESSAGE,
					       "constructor((MR< T > &, const MR< T > &) ::"
					       "Error in memory allocation (value[i])");
				value_counter[i] = M.value_counter[i];
				allocated[i] = size;
				for (register lidia_size_t p = 0; p < value_counter[i]; p++) {
					value[i][p] = bigint(M.value[i][p]);
					index[i][p] = M.index[i][p];
				}
			}
		}
	}

	bitfield = M.bitfield;
}



//
// remainder
//

void
remainder(matrix< bigint > &RES, const matrix< bigint > & M, const bigint & mod)
{
	//
	//    Task: RES.remainder(M,mod);
	//          => RES = M % mod,
	// Version: 2.0
	//

	debug_handler_l(DMESSAGE, "remainder(matrix< bigint >, const matrix< bigint > &, "
			"const bigint &)", DVALUE);

	const bigint_matrix_algorithms< DRMKex, DRMKex, DRMKex > DDD_bigint_modul;
	const bigint_matrix_algorithms< SRMKex, DRMKex, DRMKex > SDD_bigint_modul;
	const bigint_matrix_algorithms< DRMKex, SRMKex, DRMKex > DSD_bigint_modul;
	const bigint_matrix_algorithms< SRMKex, SRMKex, DRMKex > SSD_bigint_modul;

	if (RES.rows != M.rows)
		RES.set_no_of_rows(M.rows);
	if (RES.columns != M.columns)
		RES.set_no_of_columns(M.columns);

	if (M.bitfield.get_representation() == matrix_flags::dense_representation) {
		if (RES.bitfield.get_representation() == matrix_flags::dense_representation)
			DDD_bigint_modul.remainder(RES, M, mod);
		else
			SDD_bigint_modul.remainder(RES, M, mod);
	}
	else {
		if (RES.bitfield.get_representation() == matrix_flags::dense_representation)
			DSD_bigint_modul.remainder(RES, M, mod);
		else
			SSD_bigint_modul.remainder(RES, M, mod);
	}
}



void
remainder(matrix< long > &RES, const matrix< bigint > & M, long mod)
{
	//
	//    Task: RES.remainder(M,mod);
	//          => RES = M % mod,
	// Version: 2.0
	//

	debug_handler_l(DMESSAGE, "remainder(matrix< long > &, const matrix< bigint > &, "
			"const bigint &)", DVALUE);

	const bigint_matrix_algorithms< DRMKex, DRMKex, DRMKex > DDD_bigint_modul;
	const bigint_matrix_algorithms< SRMKex, DRMKex, DRMKex > SDD_bigint_modul;
	const bigint_matrix_algorithms< DRMKex, SRMKex, DRMKex > DSD_bigint_modul;
	const bigint_matrix_algorithms< SRMKex, SRMKex, DRMKex > SSD_bigint_modul;

	if (RES.get_no_of_rows() != M.rows)
		RES.set_no_of_rows(M.rows);
	if (RES.get_no_of_columns() != M.columns)
		RES.set_no_of_columns(M.columns);

	if (RES.bitfield.get_representation() == matrix_flags::dense_representation) {
		if (RES.bitfield.get_representation() == matrix_flags::dense_representation)
			DDD_bigint_modul.remainder(RES, M, mod);
		else
			SDD_bigint_modul.remainder(RES, M, mod);
	}
	else {
		if (RES.bitfield.get_representation() == matrix_flags::dense_representation)
			DSD_bigint_modul.remainder(RES, M, mod);
		else
			SSD_bigint_modul.remainder(RES, M, mod);
	}
}



//
// pseudo-division
//

void matrix< bigint >::
divide(const matrix< bigint > &A, const bigint &k)
{
	//
	//    Task: RES.divide(A, k)
	//          => RES.value[x][y] = A.value[x][y] / k,
	//             x=0,...,A.rows-1, y=0,...,A.columns-1
	// Version: 2.0
	//

	debug_handler_l(DMESSAGE, "divide(const matrix< bigint > &, const T &)", DVALUE + 1);

	const bigint_matrix_algorithms< DRMKex, SRMKex, DRMKex > DSD_bigint_modul;
	const bigint_matrix_algorithms< SRMKex, DRMKex, SRMKex > SDD_bigint_modul;
	const bigint_matrix_algorithms< SRMKex, SRMKex, DRMKex > SSD_bigint_modul;

	if (rows != A.rows)
		set_no_of_rows(A.rows);
	if (columns != A.columns)
		set_no_of_columns(A.columns);

	if (bitfield.get_representation() == matrix_flags::dense_representation) {
		if (A.bitfield.get_representation() == matrix_flags::dense_representation)
			D_bigint_modul.divide(*this, A, k);
		else
			DSD_bigint_modul.divide(*this, A, k);
	}
	else {
		if (A.bitfield.get_representation() == matrix_flags::dense_representation)
			SDD_bigint_modul.divide(*this, A, k);
		else
			SSD_bigint_modul.divide(*this, A, k);
	}
}



void matrix< bigint >::
compwise_divide(const matrix< bigint > &A, const matrix< bigint > &B)
{
	//
	//       Task: RES.compwise_divide(A, B)
	//             => RES.value[x][y] = A.value[x][y] / B.value[x][y],
	//                x=0,...,A.rows-1, y=0,...,A.columns-1
	// Conditions: A.rows == A.rows and
	//             A.columns == B.columns
	//    Version: 2.0
	//

	debug_handler_l(DMESSAGE, "compwise_divide(const matrix< bigint > &, "
			"const matrix< bigint > &)", DVALUE + 1);

	if (A.rows != B.rows || A.columns != B.columns)
		precondition_error_handler(A.rows, "A.rows", "A.rows == B.rows",
				    B.rows, "B.rows", "A.rows == B.rows",
				    A.columns, "A.columns", "A.columns == B.columns",
				    B.columns, "B.columns", "B.columns == B.columns",
				    "void matrix< bigint >::"
				    "compwise_divide(const matrix< bigint > &A, "
				    "const matrix< bigint > &B)",
				    DMESSAGE, EMESSAGE[4]);

	const bigint_matrix_algorithms< DRMKex, DRMKex, SRMKex > DDS_bigint_modul;
	const bigint_matrix_algorithms< DRMKex, SRMKex, DRMKex > DSD_bigint_modul;
	const bigint_matrix_algorithms< SRMKex, DRMKex, SRMKex > SDD_bigint_modul;
	const bigint_matrix_algorithms< DRMKex, SRMKex, SRMKex > DSS_bigint_modul;
	const bigint_matrix_algorithms< SRMKex, DRMKex, SRMKex > SDS_bigint_modul;
	const bigint_matrix_algorithms< SRMKex, SRMKex, DRMKex > SSD_bigint_modul;
	const bigint_matrix_algorithms< SRMKex, SRMKex, SRMKex > SSS_bigint_modul;

	if (rows != A.rows)
		set_no_of_rows(A.rows);
	if (columns != A.columns)
		set_no_of_columns(A.columns);

	if (bitfield.get_representation() == matrix_flags::dense_representation) {
		if (A.bitfield.get_representation() ==
		    matrix_flags::dense_representation)
			if (B.bitfield.get_representation() ==
			    matrix_flags::dense_representation)
				D_bigint_modul.compwise_divide(*this, A, B);
			else
				DDS_bigint_modul.compwise_divide(*this, A, B);
		else
			if (B.bitfield.get_representation() ==
			    matrix_flags::dense_representation)
				DSD_bigint_modul.compwise_divide(*this, A, B);
			else
				DSS_bigint_modul.compwise_divide(*this, A, B);
	}
	else {
		if (A.bitfield.get_representation() ==
		    matrix_flags::dense_representation)
			if (B.bitfield.get_representation() ==
			    matrix_flags::dense_representation)
				SDD_bigint_modul.compwise_divide(*this, A, B);
			else
				SDS_bigint_modul.compwise_divide(*this, A, B);
		else
			if (B.bitfield.get_representation() ==
			    matrix_flags::dense_representation)
				SSD_bigint_modul.compwise_divide(*this, A, B);
			else
				SSS_bigint_modul.compwise_divide(*this, A, B);
	}
}



//
// norms and bounds
//

void matrix< bigint >::
max(bigint &MAX) const
{
	//
	//    Task: A.max(res);
	//          => res = maximum of all members of matrix A
	// Version: 2.0
	//

	debug_handler_l(DMESSAGE, "max(bigint &)", DVALUE + 2);

	if (bitfield.get_representation() == matrix_flags::dense_representation)
		D_bigint_modul.max(*this, MAX);
	else
		S_bigint_modul.max(*this, MAX);
}



void matrix< bigint >::
max_abs(bigint &MAX) const
{
	//
	//    Task: A.max_abs(res);
	//          => res = absolute maximum of all members of matrix A
	// Version: 2.0
	//

	debug_handler_l(DMESSAGE, "max_abs(bigint &)", DVALUE + 2);

	if (bitfield.get_representation() == matrix_flags::dense_representation)
		D_bigint_modul.max_abs(*this, MAX);
	else
		S_bigint_modul.max_abs(*this, MAX);
}



void matrix< bigint >::
max_pos(bigint & MAX, lidia_size_t & x, lidia_size_t & y) const
{
	//
	//    Task: A.max_pos(res, x, y);
	//          => res = maximum of all members of matrix A
	//          => res = A.value[x][y]
	// Version: 2.0
	//

	debug_handler_l(DMESSAGE, "max_pos(bigint &, lidia_size_t &, lidia_size_t &)",
			DVALUE + 2);

	if (bitfield.get_representation() == matrix_flags::dense_representation)
		D_bigint_modul.max_pos(*this, MAX, x, y);
	else
		S_bigint_modul.max_pos(*this, MAX, x, y);
}



void matrix< bigint >::
max_abs_pos(bigint & MAX, lidia_size_t & x, lidia_size_t & y) const
{
	//
	//    Task: A.max_abs_pos(res, x, y);
	//          => res = absolute maximum of all members of matrix A
	//          => res = A.value[x][y]
	// Version: 2.0
	//

	debug_handler_l(DMESSAGE, "max_abs_pos(bigint &, lidia_size_t &, lidia_size_t &)",
			DVALUE + 2);

	if (bitfield.get_representation() == matrix_flags::dense_representation)
		D_bigint_modul.max_abs_pos(*this, MAX, x, y);
	else
		S_bigint_modul.max_abs_pos(*this, MAX, x, y);
}



void matrix< bigint >::
min(bigint &MIN) const
{
	//
	//    Task: A.min(MIN);
	//          => MIN = minimum of all members of matrix A
	// Version: 2.0
	//

	debug_handler_l(DMESSAGE, "min(bigint &)", DVALUE + 2);

	if (bitfield.get_representation() == matrix_flags::dense_representation)
		D_bigint_modul.min(*this, MIN);
	else
		S_bigint_modul.min(*this, MIN);
}



void matrix< bigint >::
min_abs(bigint &MIN) const
{
	//
	//    Task: A.min_abs(MIN);
	//              => MIN = absolute minimum of all members of matrix A
	// Version: 2.0
	//

	debug_handler_l(DMESSAGE, "min_abs(bigint &)", DVALUE + 2);

	if (bitfield.get_representation() == matrix_flags::dense_representation)
		D_bigint_modul.min_abs(*this, MIN);
	else
		S_bigint_modul.min_abs(*this, MIN);
}



void matrix< bigint >::
min_pos(bigint & MIN, lidia_size_t & x, lidia_size_t & y) const
{
	//
	//    Task: A.min_pos(MIN, x, y);
	//          => MIN = minimum of all members of matrix A
	//          => MIN = A.value[x][y]
	// Version: 2.0
	//

	debug_handler_l(DMESSAGE, "min_pos(bigint &, lidia_size_t &, lidia_size_t &)",
			DVALUE + 2);

	if (bitfield.get_representation() == matrix_flags::dense_representation)
		D_bigint_modul.min_pos(*this, MIN, x, y);
	else
		S_bigint_modul.min_pos(*this, MIN, x, y);
}



void matrix< bigint >::
min_abs_pos(bigint & MIN, lidia_size_t & x, lidia_size_t & y) const
{
	//
	//    Task: A.min_abs_pos(MIN, x, y);
	//          => MIN = absolute minimum of all members of matrix A
	//          => MIN = A.value[x][y]
	// Version: 2.0
	//

	debug_handler_l(DMESSAGE, "min_abs_pos(bigint &, lidia_size_t &, lidia_size_t &)",
			DVALUE + 2);

	if (bitfield.get_representation() == matrix_flags::dense_representation)
		D_bigint_modul.min_abs_pos(*this, MIN, x, y);
	else
		S_bigint_modul.min_abs_pos(*this, MIN, x, y);

}



void matrix< bigint >::
hadamard(bigint & H) const
{
	//
	//    Task: A.hadamard(H);
	//          => H = hadamard bound of matrix A
	// Version: 2.0
	//

	debug_handler_l(DMESSAGE, "hadamard(bigint &)", DVALUE + 2);

	if (bitfield.get_representation() == matrix_flags::dense_representation)
		D_bigint_modul.hadamard(*this, H);
	else
		S_bigint_modul.hadamard(*this, H);
}



void matrix< bigint >::
binary_hadamard(lidia_size_t &H) const
{
	//
	//    Task: A.hadamard2(H);
	//          => H = hadamard bound of matrix A
	// Version: 2.0
	//

	debug_handler_l(DMESSAGE, "binary_hadamard(lidia_size_t &)", DVALUE + 2);

	if (bitfield.get_representation() == matrix_flags::dense_representation)
		D_bigint_modul.binary_hadamard(*this, H);
	else
		S_bigint_modul.binary_hadamard(*this, H);
}



void matrix< bigint >::
row_norm(bigint & RES, lidia_size_t pos, long art) const
{
	//
	//    Task: A.row_norm(RES, pos, art);
	//          => RES = A(pos, 0)^art + ..... + A(pos, A.columns - 1)^art
	// Version: 2.0
	//

	debug_handler_l(DMESSAGE, "row_norm(bigint &, lidia_size_t, long)",
			DVALUE + 2);

	if (pos< 0 || pos >= rows)
		precondition_error_handler(pos, "pos", "0 <= pos < rows",
				    rows, "rows", "",
				    "void matrix< bigint >::"
				    "row_norm(bigint & RES, lidia_size_t pos, long art) const",
				    DMESSAGE, EMESSAGE[1]);

	if (bitfield.get_representation() == matrix_flags::dense_representation)
		D_bigint_modul.row_norm(*this, RES, pos, art);
	else
		S_bigint_modul.row_norm(*this, RES, pos, art);
}



void matrix< bigint >::
column_norm(bigint & RES, lidia_size_t pos, long art) const
{
	//
	//    Task: A.column_norm(RES, pos , art);
	//          => RES = A(0, pos)^art + ..... + A(A.rows - 1, pos)^art
	// Version: 2.0
	//

	debug_handler_l(DMESSAGE, "column_norm(bigint &, lidia_size_t, long)",
			DVALUE + 2);

	if (pos< 0 || pos >= columns)
		precondition_error_handler(pos, "pos", "0 <= pos < columns",
				    columns, "columns", "",
				    "void matrix< bigint >::"
				    "column_norm(bigint & RES, lidia_size_t pos, long art) const",
				    DMESSAGE, EMESSAGE[1]);

	if (bitfield.get_representation() == matrix_flags::dense_representation)
		D_bigint_modul.column_norm(*this, RES, pos, art);
	else
		S_bigint_modul.column_norm(*this, RES, pos, art);
}



//
// randomize
//

void matrix< bigint >::
randomize(const bigint & S)
{
	//
	//    Task: RES.randomize(S);
	//          => 0 <= RES.value[i][j] <= S, i=0,...,RES.rows-1,
	//             j=0,...,RES.columns-1; random values
	// Version: 2.0
	//

	debug_handler_l(DMESSAGE, "randomize(const bigint &)", DVALUE + 3);

	if (bitfield.get_representation() == matrix_flags::dense_representation)
		D_bigint_modul.randomize(*this, S);
	else
		S_bigint_modul.randomize(*this, S);
}



void matrix< bigint >::
randomize_with_det(const bigint & S, const bigint & DET)
{
	//
	//    Task: RES.randomize_with_det(S, DET);
	//          => 0 <= RES.value[i][j] <= S, i=0,...,RES.rows-1,
	//             j=0,...,RES.columns-1; random values
	//          => det(RES) == DET
	// Version: 2.0
	//

	debug_handler_l(DMESSAGE, "randomize_with_det(const bigint &, const bigint &)",
			DVALUE + 3);

	if (rows != columns)
		precondition_error_handler(rows, "rows", "rows == columns",
				    columns, "columns", "rows == columns",
				    "void matrix< bigint >::"
				    "randomize_with_det(const bigint & S, const bigint & DET)",
				    DMESSAGE, EMESSAGE[7]);

	if (bitfield.get_representation() == matrix_flags::dense_representation)
		D_bigint_modul.randomize_with_det(*this, S, DET);
	else
		S_bigint_modul.randomize_with_det(*this, S, DET);
}



void matrix< bigint >::
randomize(const bigint & S, const long d)
{
	//
	//    Task: RES.randomize(S, d);
	//          => 0 <= RES.value[i][j] <= S, i=0,...,RES.rows-1,
	//             j=0,...,RES.columns-1; random values (d %)
	// Version: 2.0
	//

	debug_handler_l(DMESSAGE, "randomize(const bigint &, const long)", DVALUE + 3);

	if (bitfield.get_representation() == matrix_flags::dense_representation)
		D_bigint_modul.randomize(*this, S, d);
	else
		S_bigint_modul.randomize(*this, S, d);
}



///////////////////////////
// BEGIN: Linear algebra //
// PART 1                //
///////////////////////////

//
// rank
//

lidia_size_t matrix< bigint >::
rank(const bigint &H) const
{
	//
	//    Task: A.rank(H) = Rank of matrix A with
	//          H = hadamard(A).
	// Version: 2.0
	//

	debug_handler_l(DMESSAGE, "rank()", DVALUE + 4);

	const modular_arithmetic< DRMK< bigint >,
		dense_fp_matrix_kernel< long, MR< long > >,
		dense_fp_matrix_kernel< bigint, MR< bigint > > > Dm_bigint_modul;

	const modular_arithmetic< SRMK< bigint >,
		sparse_fp_matrix_kernel< long, MR< long > >,
		sparse_fp_matrix_kernel< bigint, MR< bigint > > > Sm_bigint_modul;

	if (bitfield.get_representation() == matrix_flags::dense_representation)
		return Dm_bigint_modul.rank(*this, H);
	else
		return Sm_bigint_modul.rank(*this, H);
}



lidia_size_t matrix< bigint >::
rank() const
{
	//
	//    Task: A.rank() = Rank of matrix A.
	// Version: 2.0
	//

	debug_handler_l(DMESSAGE, "rank()", DVALUE + 4);

	bigint H;

	const modular_arithmetic< DRMK< bigint >,
		dense_fp_matrix_kernel< long, MR< long > >,
		dense_fp_matrix_kernel< bigint, MR< bigint > > > Dm_bigint_modul;

	const modular_arithmetic< SRMK< bigint >,
		sparse_fp_matrix_kernel< long, MR< long > >,
		sparse_fp_matrix_kernel< bigint, MR< bigint > > > Sm_bigint_modul;

	if (bitfield.get_representation() == matrix_flags::dense_representation) {
		D_bigint_modul.hadamard(*this, H);
		return Dm_bigint_modul.rank(*this, H);
	}
	else {
		S_bigint_modul.hadamard(*this, H);
		return Sm_bigint_modul.rank(*this, H);
	}
}



//
// rank and linearly independent rows
//

lidia_size_t *matrix< bigint >::
lininr1(const bigint &H) const
{
	//
	//    Task: RES[0] = Rank of matrix (Avalue,r,c).
	//          RES[1],...,RES[RES[0]], such that
	//          row(RES[1]),...,row(RES[RES[0]])
	//          of matrix *this are linearly independent.
	// Version: 2.0
	//

	debug_handler_l(DMESSAGE, "lininr1(const bigint &)", DVALUE + 4);

	const modular_arithmetic< DRMK< bigint >,
		dense_fp_matrix_kernel< long, MR< long > >,
		dense_fp_matrix_kernel< bigint, MR< bigint > > > Dm_bigint_modul;

	const modular_arithmetic< SRMK< bigint >,
		sparse_fp_matrix_kernel< long, MR< long > >,
		sparse_fp_matrix_kernel< bigint, MR< bigint > > > Sm_bigint_modul;

	if (bitfield.get_representation() == matrix_flags::dense_representation)
		return Dm_bigint_modul.lininr1(*this, H);
	else
		return Sm_bigint_modul.lininr1(*this, H);
}



lidia_size_t *matrix< bigint >::
lininr1() const
{
	//
	//    Task: RES[0] = Rank of matrix (Avalue,r,c).
	//          RES[1],...,RES[RES[0]], such that
	//          row(RES[1]),...,row(RES[RES[0]])
	//          of matrix *this are linearly independent.
	// Version: 2.0
	//

	debug_handler_l(DMESSAGE, "lininr1()", DVALUE + 4);

	bigint H;

	const modular_arithmetic< DRMK< bigint >,
		dense_fp_matrix_kernel< long, MR< long > >,
		dense_fp_matrix_kernel< bigint, MR< bigint > > > Dm_bigint_modul;

	const modular_arithmetic< SRMK< bigint >,
		sparse_fp_matrix_kernel< long, MR< long > >,
		sparse_fp_matrix_kernel< bigint, MR< bigint > > > Sm_bigint_modul;

	if (bitfield.get_representation() == matrix_flags::dense_representation) {
		D_bigint_modul.hadamard(*this, H);
		return Dm_bigint_modul.lininr1(*this, H);
	}
	else {
		S_bigint_modul.hadamard(*this, H);
		return Sm_bigint_modul.lininr1(*this, H);
	}
}



void matrix< bigint >::
lininr1(base_vector< lidia_size_t > &RES) const
{
	//
	//    Task: RES[0] = Rank of matrix (Avalue,r,c).
	//          RES[1],...,RES[RES[0]], such that
	//          row(RES[1]),...,row(RES[RES[0]])
	//          of matrix *this are linearly independent.
	// Version: 2.0
	//

	debug_handler_l(DMESSAGE, "lininr1(base_vector< lidia_size_t > &)", DVALUE + 4);

	bigint H;
	lidia_size_t *tmp;

	const modular_arithmetic< DRMK< bigint >,
		dense_fp_matrix_kernel< long, MR< long > >,
		dense_fp_matrix_kernel< bigint, MR< bigint > > > Dm_bigint_modul;

	const modular_arithmetic< SRMK< bigint >,
		sparse_fp_matrix_kernel< long, MR< long > >,
		sparse_fp_matrix_kernel< bigint, MR< bigint > > > Sm_bigint_modul;

	if (bitfield.get_representation() == matrix_flags::dense_representation) {
		D_bigint_modul.hadamard(*this, H);
		tmp = Dm_bigint_modul.lininr1(*this, H);
	}
	else {
		S_bigint_modul.hadamard(*this, H);
		tmp = Sm_bigint_modul.lininr1(*this, H);
	}

	if (RES.capacity() < tmp[0])
		RES.set_capacity(tmp[0]);
	if (RES.size() != tmp[0])
		RES.set_size(tmp[0]);

	for (register lidia_size_t i = 0; i <= tmp[0]; i++)
		RES[i] = tmp[i];

	delete[] tmp;
}



lidia_size_t *matrix< bigint >::
lininr2(const bigint &H) const
{
	//
	//    Task: RES[0] = Rank of matrix (Avalue,r,c).
	//          RES[1],...,RES[RES[0]], such that
	//          row(RES[1]),...,row(RES[RES[0]])
	//          of matrix *this are linearly independent.
	// Version: 2.0
	//

	debug_handler_l(DMESSAGE, "lininr2(const bigint &)", DVALUE + 4);

	const modular_arithmetic< DRMK< bigint >,
		dense_fp_matrix_kernel< long, MR< long > >,
		dense_fp_matrix_kernel< bigint, MR< bigint > > > Dm_bigint_modul;

	const modular_arithmetic< SRMK< bigint >,
		sparse_fp_matrix_kernel< long, MR< long > >,
		sparse_fp_matrix_kernel< bigint, MR< bigint > > > Sm_bigint_modul;

	if (bitfield.get_representation() == matrix_flags::dense_representation)
		return Dm_bigint_modul.lininr2(*this, H);
	else
		return Sm_bigint_modul.lininr2(*this, H);
}



lidia_size_t *matrix< bigint >::
lininr2() const
{
	//
	//    Task: RES[0] = Rank of matrix (Avalue,r,c).
	//          RES[1],...,RES[RES[0]], such that
	//          row(RES[1]),...,row(RES[RES[0]])
	//          of matrix *this are linearly independent.
	// Version: 2.0
	//

	debug_handler_l(DMESSAGE, "lininr2()", DVALUE + 4);

	bigint H;

	const modular_arithmetic< DRMK< bigint >,
		dense_fp_matrix_kernel< long, MR< long > >,
		dense_fp_matrix_kernel< bigint, MR< bigint > > > Dm_bigint_modul;

	const modular_arithmetic< SRMK< bigint >,
		sparse_fp_matrix_kernel< long, MR< long > >,
		sparse_fp_matrix_kernel< bigint, MR< bigint > > > Sm_bigint_modul;

	if (bitfield.get_representation() == matrix_flags::dense_representation) {
		D_bigint_modul.hadamard(*this, H);
		return Dm_bigint_modul.lininr2(*this, H);
	}
	else {
		S_bigint_modul.hadamard(*this, H);
		return Sm_bigint_modul.lininr2(*this, H);
	}
}



void matrix< bigint >::
lininr2(base_vector< lidia_size_t > &RES) const
{
	//
	//    Task: RES[0] = Rank of matrix (Avalue,r,c).
	//          RES[1],...,RES[RES[0]], such that
	//          row(RES[1]),...,row(RES[RES[0]])
	//          of matrix *this are linearly independent.
	// Version: 2.0
	//

	debug_handler_l(DMESSAGE, "lininr2(base_vector< lidia_size_t > &)", DVALUE + 4);

	bigint H;
	lidia_size_t *tmp;

	const modular_arithmetic< DRMK< bigint >,
		dense_fp_matrix_kernel< long, MR< long > >,
		dense_fp_matrix_kernel< bigint, MR< bigint > > > Dm_bigint_modul;

	const modular_arithmetic< SRMK< bigint >,
		sparse_fp_matrix_kernel< long, MR< long > >,
		sparse_fp_matrix_kernel< bigint, MR< bigint > > > Sm_bigint_modul;

	if (bitfield.get_representation() == matrix_flags::dense_representation) {
		D_bigint_modul.hadamard(*this, H);
		tmp = Dm_bigint_modul.lininr2(*this, H);
	}
	else {
		S_bigint_modul.hadamard(*this, H);
		tmp = Sm_bigint_modul.lininr2(*this, H);
	}

	if (RES.capacity() < tmp[0])
		RES.set_capacity(tmp[0]);
	if (RES.size() != tmp[0])
		RES.set_size(tmp[0]);

	for (register lidia_size_t i = 0; i <= tmp[0]; i++)
		RES[i] = tmp[i];

	delete[] tmp;
}



//
// rank linearly independent columns
//

lidia_size_t *matrix< bigint >::
lininc1(const bigint &H) const
{
	//
	//    Task: RES[0] = Rank of matrix (Avalue,r,c).
	//          RES[1],...,RES[RES[0]], such that
	//          column(RES[1]),...,column(RES[RES[0]])
	//          of matrix *this are linearly independent.
	// Version: 2.0
	//

	debug_handler_l(DMESSAGE, "lininc1(const bigint &)", DVALUE + 4);

	const modular_arithmetic< DRMK< bigint >,
		dense_fp_matrix_kernel< long, MR< long > >,
		dense_fp_matrix_kernel< bigint, MR< bigint > > > Dm_bigint_modul;

	const modular_arithmetic< SRMK< bigint >,
		sparse_fp_matrix_kernel< long, MR< long > >,
		sparse_fp_matrix_kernel< bigint, MR< bigint > > > Sm_bigint_modul;

	if (bitfield.get_representation() == matrix_flags::dense_representation)
		return Dm_bigint_modul.lininc1(*this, H);
	else
		return Sm_bigint_modul.lininc1(*this, H);
}



lidia_size_t *matrix< bigint >::
lininc1() const
{
	//
	//    Task: RES[0] = Rank of matrix (Avalue,r,c).
	//          RES[1],...,RES[RES[0]], such that
	//          column(RES[1]),...,column(RES[RES[0]])
	//          of matrix *this are linearly independent.
	// Version: 2.0
	//

	debug_handler_l(DMESSAGE, "lininc1()", DVALUE + 4);

	bigint H;

	const modular_arithmetic< DRMK< bigint >,
		dense_fp_matrix_kernel< long, MR< long > >,
		dense_fp_matrix_kernel< bigint, MR< bigint > > > Dm_bigint_modul;

	const modular_arithmetic< SRMK< bigint >,
		sparse_fp_matrix_kernel< long, MR< long > >,
		sparse_fp_matrix_kernel< bigint, MR< bigint > > > Sm_bigint_modul;

	if (bitfield.get_representation() == matrix_flags::dense_representation) {
		D_bigint_modul.hadamard(*this, H);
		return Dm_bigint_modul.lininc1(*this, H);
	}
	else {
		S_bigint_modul.hadamard(*this, H);
		return Sm_bigint_modul.lininc1(*this, H);
	}
}



void matrix< bigint >::
lininc1(base_vector< lidia_size_t > &RES) const
{
	//
	//    Task: RES[0] = Rank of matrix (Avalue,r,c).
	//          RES[1],...,RES[RES[0]], such that
	//          column(RES[1]),...,column(RES[RES[0]])
	//          of matrix *this are linearly independent.
	// Version: 2.0
	//

	debug_handler_l(DMESSAGE, "lininc1(base_vector< lidia_size_t > &)",
			DVALUE + 4);

	bigint H;
	lidia_size_t *tmp;

	const modular_arithmetic< DRMK< bigint >,
		dense_fp_matrix_kernel< long, MR< long > >,
		dense_fp_matrix_kernel< bigint, MR< bigint > > > Dm_bigint_modul;

	const modular_arithmetic< SRMK< bigint >,
		sparse_fp_matrix_kernel< long, MR< long > >,
		sparse_fp_matrix_kernel< bigint, MR< bigint > > > Sm_bigint_modul;

	if (bitfield.get_representation() == matrix_flags::dense_representation) {
		D_bigint_modul.hadamard(*this, H);
		tmp = Dm_bigint_modul.lininc1(*this, H);
	}
	else {
		S_bigint_modul.hadamard(*this, H);
		tmp = Sm_bigint_modul.lininc1(*this, H);
	}

	if (RES.capacity() < tmp[0])
		RES.set_capacity(tmp[0]);
	if (RES.size() != tmp[0])
		RES.set_size(tmp[0]);

	for (register lidia_size_t i = 0; i <= tmp[0]; i++)
		RES[i] = tmp[i];

	delete[] tmp;
}



lidia_size_t *matrix< bigint >::
lininc2(const bigint &H) const
{
	//
	//    Task: RES[0] = Rank of matrix (Avalue,r,c).
	//          RES[1],...,RES[RES[0]], such that
	//          column(RES[1]),...,column(RES[RES[0]])
	//          of matrix *this are linearly independent.
	// Version: 2.0
	//

	debug_handler_l(DMESSAGE, "lininc2(const bigint &)", DVALUE + 4);

	const modular_arithmetic< DRMK< bigint >,
		dense_fp_matrix_kernel< long, MR< long > >,
		dense_fp_matrix_kernel< bigint, MR< bigint > > > Dm_bigint_modul;

	const modular_arithmetic< SRMK< bigint >,
		sparse_fp_matrix_kernel< long, MR< long > >,
		sparse_fp_matrix_kernel< bigint, MR< bigint > > > Sm_bigint_modul;

	if (bitfield.get_representation() == matrix_flags::dense_representation)
		return Dm_bigint_modul.lininc2(*this, H);
	else
		return Sm_bigint_modul.lininc2(*this, H);
}



lidia_size_t *matrix< bigint >::
lininc2() const
{
	//
	//    Task: RES[0] = Rank of matrix (Avalue,r,c).
	//          RES[1],...,RES[RES[0]], such that
	//          column(RES[1]),...,column(RES[RES[0]])
	//          of matrix *this are linearly independent.
	// Version: 2.0
	//

	debug_handler_l(DMESSAGE, "lininc2()", DVALUE + 4);

	bigint H;

	const modular_arithmetic< DRMK< bigint >,
		dense_fp_matrix_kernel< long, MR< long > >,
		dense_fp_matrix_kernel< bigint, MR< bigint > > > Dm_bigint_modul;

	const modular_arithmetic< SRMK< bigint >,
		sparse_fp_matrix_kernel< long, MR< long > >,
		sparse_fp_matrix_kernel< bigint, MR< bigint > > > Sm_bigint_modul;

	if (bitfield.get_representation() == matrix_flags::dense_representation) {
		D_bigint_modul.hadamard(*this, H);
		return Dm_bigint_modul.lininc2(*this, H);
	}
	else {
		S_bigint_modul.hadamard(*this, H);
		return Sm_bigint_modul.lininc2(*this, H);
	}
}



void matrix< bigint >::
lininc2(base_vector< lidia_size_t > &RES) const
{
	//
	//    Task: RES[0] = Rank of matrix (Avalue,r,c).
	//          RES[1],...,RES[RES[0]], such that
	//          column(RES[1]),...,column(RES[RES[0]])
	//          of matrix *this are linearly independent.
	// Version: 2.0
	//

	debug_handler_l(DMESSAGE, "lininc2(base_vector< lidia_size_t > &)",
			DVALUE + 4);

	bigint H;
	lidia_size_t *tmp;

	const modular_arithmetic< DRMK< bigint >,
		dense_fp_matrix_kernel< long, MR< long > >,
		dense_fp_matrix_kernel< bigint, MR< bigint > > > Dm_bigint_modul;

	const modular_arithmetic< SRMK< bigint >,
		sparse_fp_matrix_kernel< long, MR< long > >,
		sparse_fp_matrix_kernel< bigint, MR< bigint > > > Sm_bigint_modul;

	if (bitfield.get_representation() == matrix_flags::dense_representation) {
		D_bigint_modul.hadamard(*this, H);
		tmp = Dm_bigint_modul.lininc2(*this, H);
	}
	else {
		S_bigint_modul.hadamard(*this, H);
		tmp = Sm_bigint_modul.lininc2(*this, H);
	}

	if (RES.capacity() < tmp[0])
		RES.set_capacity(tmp[0]);
	if (RES.size() != tmp[0])
		RES.set_size(tmp[0]);

	for (register lidia_size_t i = 0; i <= tmp[0]; i++)
		RES[i] = tmp[i];

	delete[] tmp;
}



//
// regular expansion
//

void matrix< bigint >::
regexpansion(const lidia_size_t * v)
{
	//
	//    Task: A.regexpansion(v);
	//          => A = Regular Expansion of the old matrix A relative v.
	// Version: 2.0
	//

	debug_handler_l(DMESSAGE, "regexpansion(const lidia_size_t *)", DVALUE + 4);

	if (bitfield.get_representation() == matrix_flags::dense_representation)
		D_bigint_modul.regexpansion(*this, v);
	else
		S_bigint_modul.regexpansion(*this, v);
}



//
// adjoint matrix
//

void matrix< bigint >::
adj1(const matrix< bigint > & A, const bigint &H, const bigint &DET)
{
	//
	//       Task: B.adj(A);
	//             => B = adjoint matrix of matrix A
	// Conditions: A.columns != A.rows
	//    Version: 2.0
	//

	debug_handler_l(DMESSAGE, "adj1(const matrix< bigint > &, const bigint &, const bigint &)",
			DVALUE + 4);

	if (A.columns != A.rows || DET == 0)
		precondition_error_handler(A.columns, "A.columns", "A.columns == A.rows",
				    A.rows, "A.rows", "A.columns == A.rows",
				    DET, "DET", "DET != 0",
				    "void matrix< bigint >::adj1(const matrix< bigint > & A, "
				    "const bigint &H, const bigint &DET)", DMESSAGE, EMESSAGE[7]);

	const modular_arithmetic< DRMK< bigint >,
		dense_fp_matrix_kernel< long, MR< long > >,
		dense_fp_matrix_kernel< bigint, MR< bigint > > > Dm_bigint_modul;

	const modular_arithmetic< SRMK< bigint >,
		sparse_fp_matrix_kernel< long, MR< long > >,
		sparse_fp_matrix_kernel< bigint, MR< bigint > > > Sm_bigint_modul;

	if (bitfield.get_representation() == matrix_flags::dense_representation)
		Dm_bigint_modul.adj1(*this, A, H, DET);
	else
		Sm_bigint_modul.adj1(*this, A, H, DET);
}



void matrix< bigint >::
adj1(const matrix< bigint > & A, const bigint &H)
{
	//
	//       Task: B.adj(A);
	//             => B = adjoint matrix of matrix A
	// Conditions: A.columns != A.rows
	//    Version: 2.0
	//

	debug_handler_l(DMESSAGE, "adj1(const matrix< bigint > &, const bigint &H)", DVALUE + 4);

	if (A.columns != A.rows)
		precondition_error_handler(A.columns, "A.columns", "A.columns == A.rows",
				    A.rows, "A.rows", "A.columns == A.rows",
				    "void matrix< bigint >::"
				    "adj1(const matrix< bigint > & A, const bigint &H)", DMESSAGE, EMESSAGE[7]);

	const modular_arithmetic< DRMK< bigint >,
		dense_fp_matrix_kernel< long, MR< long > >,
		dense_fp_matrix_kernel< bigint, MR< bigint > > > Dm_bigint_modul;

	const modular_arithmetic< SRMK< bigint >,
		sparse_fp_matrix_kernel< long, MR< long > >,
		sparse_fp_matrix_kernel< bigint, MR< bigint > > > Sm_bigint_modul;

	bigint DET;

	if (bitfield.get_representation() == matrix_flags::dense_representation) {
		Dm_bigint_modul.det(A, DET, H);
		if (DET == 0)
			precondition_error_handler(DET, "det(A)", "det(A) != 0", "void matrix< bigint >::"
					    "adj1(const matrix< bigint > & A, const bigint &H)",
					    DMESSAGE, EMESSAGE[7]);
		Dm_bigint_modul.adj1(*this, A, H, DET);
	}
	else {
		Sm_bigint_modul.det(A, DET, H);
		if (DET == 0)
			precondition_error_handler(DET, "det(A)", "det(A) != 0", "void matrix< bigint >::"
					    "adj1(const matrix< bigint > & A, const bigint &H)",
					    DMESSAGE, EMESSAGE[7]);
		Sm_bigint_modul.adj1(*this, A, H, DET);
	}
}



void matrix< bigint >::
adj1(const matrix< bigint > & A)
{
	//
	//       Task: B.adj(A);
	//             => B = adjoint matrix of matrix A
	// Conditions: A.columns != A.rows
	//    Version: 2.0
	//

	debug_handler_l(DMESSAGE, "adj1(const matrix< bigint > &)", DVALUE + 4);

	if (A.columns != A.rows)
		precondition_error_handler(A.columns, "A.columns", "A.columns == A.rows",
				    A.rows, "A.rows", "A.columns == A.rows",
				    "void matrix< bigint >::"
				    "adj1(const matrix< bigint > & A)", DMESSAGE, EMESSAGE[7]);

	bigint H, DET;

	const modular_arithmetic< DRMK< bigint >,
		dense_fp_matrix_kernel< long, MR< long > >,
		dense_fp_matrix_kernel< bigint, MR< bigint > > > Dm_bigint_modul;

	const modular_arithmetic< SRMK< bigint >,
		sparse_fp_matrix_kernel< long, MR< long > >,
		sparse_fp_matrix_kernel< bigint, MR< bigint > > > Sm_bigint_modul;

	if (bitfield.get_representation() == matrix_flags::dense_representation) {
		D_bigint_modul.hadamard(A, H);
		Dm_bigint_modul.det(A, DET, H);
		if (DET == 0)
			precondition_error_handler(DET, "det(A)", "det(A) != 0", "void matrix< bigint >::"
					    "adj1(const matrix< bigint > & A)", DMESSAGE, EMESSAGE[7]);
		Dm_bigint_modul.adj1(*this, A, H, DET);
	}
	else {
		S_bigint_modul.hadamard(A, H);
		Sm_bigint_modul.det(A, DET, H);
		if (DET == 0)
			precondition_error_handler(DET, "det(A)", "det(A) != 0", "void matrix< bigint >::"
					    "adj1(const matrix< bigint > & A)", DMESSAGE, EMESSAGE[7]);
		Sm_bigint_modul.adj1(*this, A, H, DET);
	}
}



void matrix< bigint >::
adj2(const matrix< bigint > & A, const bigint &H, const bigint &DET)
{
	//
	//       Task: B.adj(A);
	//             => B = adjoint matrix of matrix A
	// Conditions: A.columns != A.rows
	//    Version: 2.0
	//

	debug_handler_l(DMESSAGE, "adj2(const matrix< bigint > &, const bigint &, const bigint &)",
			DVALUE + 4);

	if (A.columns != A.rows || DET == 0)
		precondition_error_handler(A.columns, "A.columns", "A.columns == A.rows",
				    A.rows, "A.rows", "A.columns == A.rows",
				    DET, "DET", "DET != 0",
				    "void matrix< bigint >::adj2(const matrix< bigint > & A, "
				    "const bigint &H, const bigint &DET)", DMESSAGE, EMESSAGE[7]);

	const modular_arithmetic< DRMK< bigint >,
		dense_fp_matrix_kernel< long, MR< long > >,
		dense_fp_matrix_kernel< bigint, MR< bigint > > > Dm_bigint_modul;

	const modular_arithmetic< SRMK< bigint >,
		sparse_fp_matrix_kernel< long, MR< long > >,
		sparse_fp_matrix_kernel< bigint, MR< bigint > > > Sm_bigint_modul;

	if (bitfield.get_representation() == matrix_flags::dense_representation)
		Dm_bigint_modul.adj2(*this, A, H, DET);
	else
		Sm_bigint_modul.adj2(*this, A, H, DET);
}



void matrix< bigint >::
adj2(const matrix< bigint > & A, const bigint &H)
{
	//
	//       Task: B.adj(A);
	//             => B = adjoint matrix of matrix A
	// Conditions: A.columns != A.rows
	//    Version: 2.0
	//

	debug_handler_l(DMESSAGE, "adj2(const matrix< bigint > &, const bigint &H)", DVALUE + 4);

	if (A.columns != A.rows)
		precondition_error_handler(A.columns, "A.columns", "A.columns == A.rows",
				    A.rows, "A.rows", "A.columns == A.rows",
				    "void matrix< bigint >::"
				    "adj2(const matrix< bigint > & A, const bigint &H)", DMESSAGE, EMESSAGE[7]);

	const modular_arithmetic< DRMK< bigint >,
		dense_fp_matrix_kernel< long, MR< long > >,
		dense_fp_matrix_kernel< bigint, MR< bigint > > > Dm_bigint_modul;

	const modular_arithmetic< SRMK< bigint >,
		sparse_fp_matrix_kernel< long, MR< long > >,
		sparse_fp_matrix_kernel< bigint, MR< bigint > > > Sm_bigint_modul;

	bigint DET;

	if (bitfield.get_representation() == matrix_flags::dense_representation) {
		Dm_bigint_modul.det(A, DET, H);
		if (DET == 0)
			precondition_error_handler(DET, "det(A)", "det(A) != 0", "void matrix< bigint >::"
					    "adj2(const matrix< bigint > & A, const bigint &H)",
					    DMESSAGE, EMESSAGE[7]);
		Dm_bigint_modul.adj2(*this, A, H, DET);
	}
	else {
		Sm_bigint_modul.det(A, DET, H);
		if (DET == 0)
			precondition_error_handler(DET, "det(A)", "det(A) != 0", "void matrix< bigint >::"
					    "adj2(const matrix< bigint > & A, const bigint &H)",
					    DMESSAGE, EMESSAGE[7]);
		Sm_bigint_modul.adj2(*this, A, H, DET);
	}
}



void matrix< bigint >::
adj2(const matrix< bigint > & A)
{
	//
	//       Task: B.adj(A);
	//             => B = adjoint matrix of matrix A
	// Conditions: A.columns != A.rows
	//    Version: 2.0
	//

	debug_handler_l(DMESSAGE, "adj2(const matrix< bigint > &)", DVALUE + 4);

	if (A.columns != A.rows)
		precondition_error_handler(A.columns, "A.columns", "A.columns == A.rows",
				    A.rows, "A.rows", "A.columns == A.rows",
				    "void matrix< bigint >::"
				    "adj2(const matrix< bigint > & A)", DMESSAGE, EMESSAGE[7]);

	bigint H, DET;

	const modular_arithmetic< DRMK< bigint >,
		dense_fp_matrix_kernel< long, MR< long > >,
		dense_fp_matrix_kernel< bigint, MR< bigint > > > Dm_bigint_modul;

	const modular_arithmetic< SRMK< bigint >,
		sparse_fp_matrix_kernel< long, MR< long > >,
		sparse_fp_matrix_kernel< bigint, MR< bigint > > > Sm_bigint_modul;

	if (bitfield.get_representation() == matrix_flags::dense_representation) {
		D_bigint_modul.hadamard(A, H);
		Dm_bigint_modul.det(A, DET, H);
		if (DET == 0)
			precondition_error_handler(DET, "det(A)", "det(A) != 0", "void matrix< bigint >::"
					    "adj2(const matrix< bigint > & A)", DMESSAGE, EMESSAGE[7]);
		Dm_bigint_modul.adj2(*this, A, H, DET);
	}
	else {
		S_bigint_modul.hadamard(A, H);
		Sm_bigint_modul.det(A, DET, H);
		if (DET == 0)
			precondition_error_handler(DET, "det(A)", "det(A) != 0", "void matrix< bigint >::"
					    "adj2(const matrix< bigint > & A)", DMESSAGE, EMESSAGE[7]);
		Sm_bigint_modul.adj2(*this, A, H, DET);
	}
}



//
// lattice determinant
//

void matrix< bigint >::
latticedet1(bigint & DET, const bigint &H) const
{
	//
	//    Task: A.latticedet1(DET, H);
	//          => DET = lattice determinant
	//          of the lattice formed by the columns of matrix A
	//          => H = hadamard(A)
	// Version: 2.0
	//

	debug_handler_l(DMESSAGE, "latticedet1(bigint &, const bigint &)", DVALUE + 4);

	const modular_arithmetic< DRMK< bigint >,
		dense_fp_matrix_kernel< long, MR< long > >,
		dense_fp_matrix_kernel< bigint, MR< bigint > > > Dm_bigint_modul;

	const modular_arithmetic< SRMK< bigint >,
		sparse_fp_matrix_kernel< long, MR< long > >,
		sparse_fp_matrix_kernel< bigint, MR< bigint > > > Sm_bigint_modul;

	if (bitfield.get_representation() == matrix_flags::dense_representation)
		Dm_bigint_modul.latticedet1(*this, DET, H);
	else
		Sm_bigint_modul.latticedet1(*this, DET, H);
}



void matrix< bigint >::
latticedet1(bigint & DET) const
{
	//
	//    Task: A.latticedet1(DET);
	//          => DET = lattice determinant
	//          of the lattice formed by the columns of matrix A
	// Version: 2.0
	//

	debug_handler_l(DMESSAGE, "latticedet1(bigint &)", DVALUE + 4);

	bigint H;

	const modular_arithmetic< DRMK< bigint >,
		dense_fp_matrix_kernel< long, MR< long > >,
		dense_fp_matrix_kernel< bigint, MR< bigint > > > Dm_bigint_modul;

	const modular_arithmetic< SRMK< bigint >,
		sparse_fp_matrix_kernel< long, MR< long > >,
		sparse_fp_matrix_kernel< bigint, MR< bigint > > > Sm_bigint_modul;

	if (bitfield.get_representation() == matrix_flags::dense_representation) {
		D_bigint_modul.hadamard(*this, H);
		Dm_bigint_modul.latticedet1(*this, DET, H);
	}
	else {
		S_bigint_modul.hadamard(*this, H);
		Sm_bigint_modul.latticedet1(*this, DET, H);
	}
}



void matrix< bigint >::
latticedet2(bigint & DET, const bigint &H) const
{
	//
	//    Task: A.latticedet2(DET, H);
	//          => DET = lattice determinant
	//          of the lattice formed by the columns of matrix A
	//          => H = hadamard(A)
	// Version: 2.0
	//

	debug_handler_l(DMESSAGE, "latticedet2(bigint &, const bigint &)", DVALUE + 4);

	const modular_arithmetic< DRMK< bigint >,
		dense_fp_matrix_kernel< long, MR< long > >,
		dense_fp_matrix_kernel< bigint, MR< bigint > > > Dm_bigint_modul;

	const modular_arithmetic< SRMK< bigint >,
		sparse_fp_matrix_kernel< long, MR< long > >,
		sparse_fp_matrix_kernel< bigint, MR< bigint > > > Sm_bigint_modul;

	if (bitfield.get_representation() == matrix_flags::dense_representation)
		Dm_bigint_modul.latticedet2(*this, DET, H);
	else
		Sm_bigint_modul.latticedet2(*this, DET, H);
}



void matrix< bigint >::
latticedet2(bigint & DET) const
{
	//
	//    Task: A.latticedet2(DET);
	//          => DET = lattice determinant
	//          of the lattice formed by the columns of matrix A
	// Version: 2.0
	//

	debug_handler_l(DMESSAGE, "latticedet2(bigint &)", DVALUE + 4);

	bigint H;

	const modular_arithmetic< DRMK< bigint >,
		dense_fp_matrix_kernel< long, MR< long > >,
		dense_fp_matrix_kernel< bigint, MR< bigint > > > Dm_bigint_modul;

	const modular_arithmetic< SRMK< bigint >,
		sparse_fp_matrix_kernel< long, MR< long > >,
		sparse_fp_matrix_kernel< bigint, MR< bigint > > > Sm_bigint_modul;

	if (bitfield.get_representation() == matrix_flags::dense_representation) {
		D_bigint_modul.hadamard(*this, H);
		Dm_bigint_modul.latticedet2(*this, DET, H);
	}
	else {
		S_bigint_modul.hadamard(*this, H);
		Sm_bigint_modul.latticedet2(*this, DET, H);
	}
}



void matrix< bigint >::
latticedet3(bigint & DET, const bigint &H) const
{
	//
	//    Task: A.latticedet3(DET, H);
	//          => DET = lattice determinant
	//          of the lattice formed by the columns of matrix A
	//          => H = hadamard(A)
	// Version: 2.0
	//

	debug_handler_l(DMESSAGE, "latticedet3(bigint &, const bigint &)", DVALUE + 4);

	const modular_arithmetic< DRMK< bigint >,
		dense_fp_matrix_kernel< long, MR< long > >,
		dense_fp_matrix_kernel< bigint, MR< bigint > > > Dm_bigint_modul;

	const modular_arithmetic< SRMK< bigint >,
		sparse_fp_matrix_kernel< long, MR< long > >,
		sparse_fp_matrix_kernel< bigint, MR< bigint > > > Sm_bigint_modul;

	if (bitfield.get_representation() == matrix_flags::dense_representation)
		Dm_bigint_modul.latticedet3(*this, DET, H);
	else
		Sm_bigint_modul.latticedet3(*this, DET, H);
}



void matrix< bigint >::
latticedet3(bigint & DET) const
{
	//
	//    Task: A.latticedet3(DET);
	//          => DET = lattice determinant
	//          of the lattice formed by the columns of matrix A
	// Version: 2.0
	//

	debug_handler_l(DMESSAGE, "latticedet3(bigint &)", DVALUE + 4);

	bigint H;

	const modular_arithmetic< DRMK< bigint >,
		dense_fp_matrix_kernel< long, MR< long > >,
		dense_fp_matrix_kernel< bigint, MR< bigint > > > Dm_bigint_modul;

	const modular_arithmetic< SRMK< bigint >,
		sparse_fp_matrix_kernel< long, MR< long > >,
		sparse_fp_matrix_kernel< bigint, MR< bigint > > > Sm_bigint_modul;

	if (bitfield.get_representation() == matrix_flags::dense_representation) {
		D_bigint_modul.hadamard(*this, H);
		Dm_bigint_modul.latticedet3(*this, DET, H);
	}
	else {
		S_bigint_modul.hadamard(*this, H);
		Sm_bigint_modul.latticedet3(*this, DET, H);
	}
}



void matrix< bigint >::
latticedet_special(bigint & DET) const
{
	//
	//    Task: A.latticedet3(DET);
	//          => DET = lattice determinant
	//          of the lattice formed by the columns of matrix A
	// Version: 2.0
	//

	debug_handler_l(DMESSAGE, "latticedet3(bigint &)", DVALUE + 4);

	bigint H;

	const modular_arithmetic< DRMK< bigint >,
		dense_fp_matrix_kernel< long, MR< long > >,
		dense_fp_matrix_kernel< bigint, MR< bigint > > > Dm_bigint_modul;

	const modular_arithmetic< SRMK< bigint >,
		sparse_fp_matrix_kernel< long, MR< long > >,
		sparse_fp_matrix_kernel< bigint, MR< bigint > > > Sm_bigint_modul;

	if (bitfield.get_representation() == matrix_flags::dense_representation) {
		D_bigint_modul.hadamard(*this, H);
		Dm_bigint_modul.latticedet_special(*this, DET, H, 5);
	}
	else {
		S_bigint_modul.hadamard(*this, H);
		Sm_bigint_modul.latticedet_special(*this, DET, H, 5);
	}
}



void matrix< bigint >::
real_latticedet(bigint & DET, const bigint &H) const
{
	//
	//    Task: A.real_latticedet(DET, H);
	//          => DET = lattice determinant
	//          of the lattice formed by the columns of matrix A
	//          => H = hadamard(A)
	// Version: 2.0
	//

	debug_handler_l(DMESSAGE, "real_latticedet(bigint &, const bigint &)",
			DVALUE + 4);

	const modular_arithmetic< DRMK< bigint >,
		dense_fp_matrix_kernel< long, MR< long > >,
		dense_fp_matrix_kernel< bigint, MR< bigint > > > Dm_bigint_modul;

	const modular_arithmetic< SRMK< bigint >,
		sparse_fp_matrix_kernel< long, MR< long > >,
		sparse_fp_matrix_kernel< bigint, MR< bigint > > > Sm_bigint_modul;

	if (bitfield.get_representation() == matrix_flags::dense_representation) {
		matrix< bigint > A(columns, rows);
		A.trans(*this);
		LiDIA::multiply(A, A, *this);
		Dm_bigint_modul.det(A, DET, H*H);
	}
	else {
		sparse_fp_matrix_algorithms< long, ring_matrix< long > > long_modul;
		sparse_fp_matrix_algorithms< bigint, ring_matrix< bigint > > bigint_modul;

		// Step 1,2
		bigint *PRIM = get_primes(bigint(H*H), bigint(1));
		long n, Modlong;
		PRIM[0].longify(n);

		bigint Gmod;

		// Step 3
		bigint *chininput = new bigint[n];
		memory_handler(chininput, DMESSAGE, "det :: "
			       "Error in memory allocation (chininput)");

		lidia_size_t i;
		for (i = 1; i <= n; i++) {
			Gmod.assign(PRIM[i]);
			if (Gmod.bit_length() >= bigint::bits_per_digit()) {
				// bigint part
				base_power_product< ring_matrix< bigint >, lidia_size_t > bpp;
				matrix< bigint > A(rows, columns);
				matrix< bigint > AT(columns, rows);
				remainder(A, *this, Gmod);
				AT.trans(A);

				bpp.append(AT);
				bpp.append(A);

				//chininput[i - 1] = bigint_modul.det(bpp, Gmod);
				chininput[i - 1] = bigint_modul.det(A, Gmod);
			}
			else {
				// long part
				Gmod.longify(Modlong);
				base_power_product< ring_matrix< long >, lidia_size_t > bpp;
				matrix< long > A(rows, columns);
				matrix< long > AT(columns, rows);
				A.set_zero_element(0);
				AT.set_zero_element(0);
				remainder(A, *this, Modlong);
				AT.trans(A);

				bpp.append(AT);
				bpp.append(A);

				//std::cout << "Prime " << i << "/" << n << std::endl;
				//chininput[i - 1] = long_modul.det(bpp, Modlong);
				chininput[i - 1] = long_modul.det(A, Modlong);
			}
		}

		// Step 4
		LiDIA::chinrest(DET, static_cast<const bigint *>(chininput), static_cast<const bigint *>(PRIM));
		delete[] chininput;
	}
}



//
// determinant
//

void matrix< bigint >::
det(bigint & DET, const bigint &H) const
{
	//
	//       Task: A.det(DET, H)
	//             => DET = determinant of matrix A
	//             => H = hadamard(A)
	// Conditions: columns != rows
	//    Version: 2.0
	//

	debug_handler_l(DMESSAGE, "det(bigint &, const bigint &)", DVALUE + 4);

	if (rows != columns)
		precondition_error_handler(rows, "rows", "rows == columns",
				    columns, "columns", "rows == columns",
				    "void matrix< bigint >::"
				    "det(bigint &DET, const bigint &H) const",
				    DMESSAGE, EMESSAGE[7]);

	const modular_arithmetic< DRMK< bigint >,
		dense_fp_matrix_kernel< long, MR< long > >,
		dense_fp_matrix_kernel< bigint, MR< bigint > > > Dm_bigint_modul;

	const modular_arithmetic< SRMK< bigint >,
		sparse_fp_matrix_kernel< long, MR< long > >,
		sparse_fp_matrix_kernel< bigint, MR< bigint > > > Sm_bigint_modul;

	if (bitfield.get_representation() == matrix_flags::dense_representation)
		Dm_bigint_modul.det(*this, DET, H);
	else
		Sm_bigint_modul.det(*this, DET, H);
}



void matrix< bigint >::
det(bigint & DET, const bigint &H, int num) const
{
	//
	//       Task: A.det(DET, H, num)
	//             => DET = determinant of matrix A
	// Conditions: columns != rows
	//    Version: 2.0
	//

	debug_handler_l(DMESSAGE, "det(bigint &)", DVALUE + 4);

	if (rows != columns)
		precondition_error_handler(rows, "rows", "rows == columns",
				    columns, "columns", "rows == columns",
				    "void matrix< bigint >::"
				    "det(bigint & DET) const",
				    DMESSAGE, EMESSAGE[7]);

	const modular_arithmetic< DRMK< bigint >,
		dense_fp_matrix_kernel< long, MR< long > >,
		dense_fp_matrix_kernel< bigint, MR< bigint > > > Dm_bigint_modul;

	const modular_arithmetic< SRMK< bigint >,
		sparse_fp_matrix_kernel< long, MR< long > >,
		sparse_fp_matrix_kernel< bigint, MR< bigint > > > Sm_bigint_modul;

	if (bitfield.get_representation() == matrix_flags::dense_representation)
		Dm_bigint_modul.det(*this, DET, H, num);
	else
		Sm_bigint_modul.det(*this, DET, H, num);
}



void matrix< bigint >::
det(bigint & DET) const
{
	//
	//       Task: A.det(DET)
	//             => DET = determinant of matrix A
	// Conditions: columns != rows
	//    Version: 2.0
	//

	debug_handler_l(DMESSAGE, "det(bigint &)", DVALUE + 4);

	if (rows != columns)
		precondition_error_handler(rows, "rows", "rows == columns",
				    columns, "columns", "rows == columns",
				    "void matrix< bigint >::"
				    "det(bigint & DET) const",
				    DMESSAGE, EMESSAGE[7]);

	bigint H;

	const modular_arithmetic< DRMK< bigint >,
		dense_fp_matrix_kernel< long, MR< long > >,
		dense_fp_matrix_kernel< bigint, MR< bigint > > > Dm_bigint_modul;

	const modular_arithmetic< SRMK< bigint >,
		sparse_fp_matrix_kernel< long, MR< long > >,
		sparse_fp_matrix_kernel< bigint, MR< bigint > > > Sm_bigint_modul;

	if (bitfield.get_representation() == matrix_flags::dense_representation) {
		D_bigint_modul.hadamard(*this, H);
		Dm_bigint_modul.det(*this, DET, H);
	}
	else {
		S_bigint_modul.hadamard(*this, H);
		Sm_bigint_modul.det(*this, DET, H);
	}
}



//
// characteristic polynomial
//

bigint *matrix< bigint >::
charpoly() const
{
	//
	//       Task: RES = A.charpoly();
	//             => RES[0],...,RES[r] are the coefficients of
	//             the characteristic polynomial of matrix A
	// Conditions: rows != columns
	//    Version: 2.0
	//

	debug_handler_l(DMESSAGE, "charpoly()", DVALUE + 4);

	if (columns != rows)
		precondition_error_handler(rows, "rows", "rows == columns",
				    columns, "columns", "rows == columns",
				    "bigint *matrix< bigint >::"
				    "charpoly() const",
				    DMESSAGE, EMESSAGE[7]);

	bigint H;
	bigint *RES = new bigint[columns + 1];
	memory_handler(RES, DMESSAGE, "charpoly :: "
		       "Error in memory allocation (RES)");

	const modular_arithmetic< DRMK< bigint >,
		dense_fp_matrix_kernel< long, MR< long > >,
		dense_fp_matrix_kernel< bigint, MR< bigint > > > Dm_bigint_modul;

	const modular_arithmetic< SRMK< bigint >,
		sparse_fp_matrix_kernel< long, MR< long > >,
		sparse_fp_matrix_kernel< bigint, MR< bigint > > > Sm_bigint_modul;

	if (bitfield.get_representation() == matrix_flags::dense_representation) {
		D_bigint_modul.hadamard(*this, H);
		Dm_bigint_modul.charpoly(*this, RES, H);
	}
	else {
		S_bigint_modul.hadamard(*this, H);
		Sm_bigint_modul.charpoly(*this, RES, H);
	}

	return RES;
}



void matrix< bigint >::
charpoly(base_vector< bigint > &RES) const
{
	//
	//       Task: A.charpoly(RES);
	//             => RES[0],...,RES[r] are the coefficients of
	//             the characteristic polynomial of matrix A
	// Conditions: rows == columns
	//    Version: 2.0
	//

	debug_handler_l(DMESSAGE, "charpoly(base_vector< bigint > &)", DVALUE + 4);

	if (columns != rows)
		precondition_error_handler(rows, "rows", "rows == columns",
				    columns, "columns", "rows == columns",
				    "void matrix< bigint >::"
				    "charpoly(base_vector< bigint > &) const",
				    DMESSAGE, EMESSAGE[7]);

	bigint H;
	RES.set_size(columns + 1);

	const modular_arithmetic< DRMK< bigint >,
		dense_fp_matrix_kernel< long, MR< long > >,
		dense_fp_matrix_kernel< bigint, MR< bigint > > > Dm_bigint_modul;

	const modular_arithmetic< SRMK< bigint >,
		sparse_fp_matrix_kernel< long, MR< long > >,
		sparse_fp_matrix_kernel< bigint, MR< bigint > > > Sm_bigint_modul;

	if (bitfield.get_representation() == matrix_flags::dense_representation) {
		D_bigint_modul.hadamard(*this, H);
		Dm_bigint_modul.charpoly(*this, RES.get_data_address(), H);
	}
	else {
		S_bigint_modul.hadamard(*this, H);
		Sm_bigint_modul.charpoly(*this, RES.get_data_address(), H);
	}
}



/////////////////////////
// END: Linear algebra //
// PART 1              //
/////////////////////////

///////////////////////////
// BEGIN: Linear algebra //
// PART 2                //
///////////////////////////

//
// Hermite normal form
//

void matrix< bigint >::
hnfmod_dkt(const bigint &mod)
{
	//
	//       Task: A.hnfmod_dkt(mod);
	//             => A in Hermite normal form
	//             => h = lattice determinant of lattice formed
	//             by the columns of matrix A
	// Conditions: rank != rows
	//    Version: 2.0
	//

	debug_handler_l(DMESSAGE, "hnfmod_dkt(const bigint &)", DVALUE + 5);

	lidia_size_t r = rank();
	if (r != rows)
		precondition_error_handler(r, "rank", "rank == rows",
				    rows, "rows", "rank == rows",
				    "void matrix< bigint >::"
				    "hnfmod_dkt(const bigint &mod)",
				    DMESSAGE, EMESSAGE[10]);

	if (bitfield.get_representation() == matrix_flags::dense_representation)
		D_bigint_modul.hnfmod_dkt(*this, mod);
	else
                lidia_error_handler("matrix<bigint>", 
                                    "hnfmod_dkt: dense matrices supported only");
}



void matrix< bigint >::
hnfmod_dkt(matrix< bigint > &TR, const bigint &mod)
{
	//
	//       Task: A.hnfmod_dkt(TR, mod);
	//             => A in Hermite normal form
	//             => h = lattice determinant of lattice formed
	//             by the columns of matrix A
	// Conditions: rank != rows
	//    Version: 2.0
	//

	debug_handler_l(DMESSAGE, "hnfmod_dkt(const bigint &)", DVALUE + 5);

	if (rank() != rows)
		precondition_error_handler(rank(), "rank", "rank == rows",
				    rows, "rows", "rank == rows",
				    "void matrix< bigint >::"
				    "hnfmod_dkt(const bigint &mod)",
				    DMESSAGE, EMESSAGE[10]);

	TR.resize(rows + columns, rows + columns);
	TR.diag(1, 0);

	if (bitfield.get_representation() == matrix_flags::dense_representation)
		D_bigint_modul.hnfmod_dkt(*this, TR, mod);
	else
                lidia_error_handler("matrix<bigint>", 
                                    "hnfmod_dkt: dense matrices supported only");
}



void matrix< bigint >::
hnfmod_cohen(const bigint & D)
{
	debug_handler_l(DMESSAGE, "hnfmod_cohen(const bigint &)", DVALUE + 5);

	if (rank() != rows)
		precondition_error_handler(rank(), "rank", "rows == rank",
				    rows, "rows", "rows == rank",
				    "void matrix< bigint >::"
				    "hnfmod_cohen(const bigint & D)",
				    DMESSAGE, EMESSAGE[10]);

	if (bitfield.get_representation() == matrix_flags::dense_representation)
		D_bigint_modul.hnfmod_cohen(*this, D);
	else
                lidia_error_handler("matrix<bigint>", 
                                    "hnfmod_cohen: dense matrices supported only");
}



void matrix< bigint >::
hnfmod_mueller(matrix< bigint > & TRANS)
{
	//
	//       Task: A.hnfmod_mueller(TRANS);
	//             => A in Hermite normal form
	//             => TRANS = transformtion matrix A*TRANS = HNFmod(A)
	// Conditions: TRANS.rows != TRANS.columns or TRANS.rows != columns or
	//             rank != rows
	//    Version: 2.0
	//

	debug_handler_l(DMESSAGE, "hnfmod_mueller(matrix< bigint > &, const bigint &)", DVALUE + 5);

	const modular_arithmetic< DRMK< bigint >, dense_fp_matrix_kernel< long, MR< long > >,
		dense_fp_matrix_kernel< bigint, MR< bigint > > > Dm_bigint_modul;

	lidia_size_t *linuz = lininc();
	if (linuz[0] != rows)
		precondition_error_handler(linuz[0], "rank", "rank == rows",
				    rows, "rows", "rank == rows",
				    "void matrix< bigint >::"
				    "hnfmod_mueller(matrix< bigint > & TRANS)",
				    DMESSAGE, EMESSAGE[10]);

	if (bitfield.get_representation() == matrix_flags::dense_representation) {
		bigint DET, H;
		D_bigint_modul.hadamard(*this, H);
		Dm_bigint_modul.latticedet1(*this, DET, H);
		D_bigint_modul.hnfmod_mueller(*this, TRANS, DET);
	}
	else
                lidia_error_handler("matrix<bigint>", 
                                    "hnfmod_mueller: dense matrices supported only");
}



void matrix< bigint >::
hnf_ref()
{
	if (bitfield.get_representation() == matrix_flags::dense_representation)
		hnf_ref_modul.hnf_havas(*this);
	else
                lidia_error_handler("matrix<bigint>", 
                                    "hnf_ref: dense matrices supported only");
}




//
// hnf_storjohann
//

void matrix< bigint >::
hnf_storjohann()
{
	//
	// Task: HNF Computation
	// Algorithm: Gauss with reduction
	// IMPROVEMENTS: Theory of Havas / best reaminder
	// Version: 2.0
	//

	debug_handler_l(DMESSAGE, "hnf_storjohann()", DVALUE + 5);

	if (bitfield.get_representation() == matrix_flags::dense_representation) {
		D_bigint_modul.hnf_storjohann(*this);
	}
	else {
		S_bigint_modul.hnf_storjohann(*this);
	}
}



void matrix< bigint >::
hnf_storjohann(matrix< bigint > &T, matrix< bigint > &C, matrix< bigint > &Q)
{
	//
	// Task: HNF Computation
	// Algorithm: Gauss with reduction
	// IMPROVEMENTS: Theory of Havas / best reaminder
	// Version: 2.0
	//

	debug_handler_l(DMESSAGE, "hnf_storjohann()", DVALUE + 5);

	lidia_size_t i = rows, j = columns;
	if (T.columns != columns)
		T.set_no_of_columns(columns);
	if (T.rows != columns)
		T.set_no_of_rows(columns);

	if (C.columns != columns)
		C.set_no_of_columns(columns);
	if (C.rows != columns)
		C.set_no_of_rows(columns);

	if (Q.columns != columns)
		Q.set_no_of_columns(columns);
	if (Q.rows != columns)
		Q.set_no_of_rows(columns);

	T.diag(1, 0);
	C.diag(1, 0);
	Q.diag(1, 0);

	if (bitfield.get_representation() == matrix_flags::dense_representation)
		D_bigint_modul.hnf_storjohann(*this, T, C, Q);
	else {

		set_orientation(matrix_flags::column_oriented);
		i = 0;
		j = 0;
		// hnf_smodul2.hnf(*this, i, j);
		set_orientation(matrix_flags::row_oriented);
	}
}



//
// hnf_havas
//

#define schema_havas_call_sequenz(name)                                  \
{if (bitfield.get_representation() == matrix_flags::dense_representation)\
{                                                                        \
  havas_kernel< bigint, RODMM< bigint >, RODMM< bigint >, \
    name< bigint, RODMM< bigint >, RODMM< bigint > > > modul; \
  modul.hnf_Z1(*this, i, j); \
}                                                                        \
else {                                                                   \
  havas_kernel< bigint, COSMM< bigint >, COSMM< bigint >, \
    name< bigint, COSMM< bigint >, COSMM< bigint > > > modul; \
                                                                         \
  set_orientation(matrix_flags::column_oriented); \
  modul.hnf_Z1(*this, i, j); \
  set_orientation(matrix_flags::row_oriented); \
}                                                                        \
break; }

#define schema_havas_stf_call_sequenz(name)                              \
{if (bitfield.get_representation() == matrix_flags::dense_representation)\
{                                                                        \
  havas_kernel< bigint, RODMM< bigint >, RODMM< bigint >, \
    name< bigint, RODMM< bigint >, RODMM< bigint > > > modul; \
  modul.hnf_Z2(*this, i, j); \
  normalization_kernel< bigint, RODMM< bigint >, \
                       RODMM< bigint > > nmodul; \
                                                                         \
  switch(normalizeModul) {                                               \
    case 0:                                                              \
      {                                                                  \
        nmodul.normalize_Std(*this, 0, columns - rows); \
        break; \
      }                                                                  \
    case 1:                                                              \
      {                                                                  \
        nmodul.normalize_ChouCollins(*this, 0, columns - rows); \
        break; \
      }                                                                  \
    case 2:                                                              \
      {                                                                  \
        nmodul.normalizeMod_Std(*this, 0, columns - rows); \
        break; \
      }                                                                  \
    case 3:                                                              \
      {                                                                  \
        nmodul.normalizeMod_ChouCollins(*this, 0, columns - rows); \
        break; \
      }                                                                  \
    case 4:                                                              \
      {                                                                  \
        nmodul.normalizeHybrid_Std(*this, 0, columns - rows); \
        break; \
      }                                                                  \
    case 5:                                                              \
      {                                                                  \
        nmodul.normalizeHybrid_ChouCollins(*this, 0, columns - rows); \
        break; \
      }                                                                  \
    default:                                                             \
      lidia_error_handler("bigint_matrix", "hnf_havas :: "               \
                          "Mode not supported! (normalize "              \
                          "Modul not defined)"); \
    }                                                                    \
}                                                                        \
else {                                                                   \
  havas_kernel< bigint, COSMM< bigint >, COSMM< bigint >, \
    name< bigint, COSMM< bigint >, COSMM< bigint > > > modul; \
                                                                         \
  set_orientation(matrix_flags::column_oriented); \
  modul.hnf_Z2(*this, i, j); \
  normalization_kernel< bigint, COSMM< bigint >, \
                       COSMM< bigint > > nmodul; \
                                                                         \
  switch(normalizeModul) {                                               \
    case 0:                                                              \
      {                                                                  \
        nmodul.normalize_Std(*this, 0, columns - rows); \
        break; \
      }                                                                  \
    case 1:                                                              \
      {                                                                  \
        nmodul.normalize_ChouCollins(*this, 0, columns - rows); \
        break; \
      }                                                                  \
    case 2:                                                              \
      {                                                                  \
        nmodul.normalize_Std(*this, 0, columns - rows); \
        break; \
      }                                                                  \
    case 3:                                                              \
      {                                                                  \
        nmodul.normalize_Std(*this, 0, columns - rows); \
        break; \
      }                                                                  \
    case 4:                                                              \
      {                                                                  \
        nmodul.normalizeHybrid_Std(*this, 0, columns - rows); \
        break; \
      }                                                                  \
    case 5:                                                              \
      {                                                                  \
        nmodul.normalizeHybrid_ChouCollins(*this, 0, columns - rows); \
        break; \
      }                                                                  \
    default:                                                             \
      lidia_error_handler("bigint_matrix", "hnf_havas :: "               \
                          "Mode not supported! (normalize "              \
                          "Modul not defined)"); \
    }                                                                    \
  set_orientation(matrix_flags::row_oriented); \
}                                                                        \
break; }

void matrix< bigint >::
hnf_havas(lidia_size_t KernAlgo, lidia_size_t mgcdModul, lidia_size_t normalizeModul)
{
	//
	// Task: HNF Computation
	// Algorithm: Gauss with reduction
	// IMPROVEMENTS: Theory of Havas / best reaminder
	// Version: 2.0
	//

	debug_handler_l(DMESSAGE, "hnf_havas()", DVALUE + 5);

	lidia_size_t i = rows, j = columns;

	switch(KernAlgo) {
	case 0:
{
		switch(mgcdModul) {
		case 0:
			schema_havas_call_sequenz(hermite)
				case 1:
					schema_havas_call_sequenz(Bradley)
				case 2:
					schema_havas_call_sequenz(Iliopoulos)
				case 3:
					schema_havas_call_sequenz(opt)
				case 4:
					schema_havas_call_sequenz(Blankinship)
				case 5:
					schema_havas_call_sequenz(Best_remainder)
				case 6:
					schema_havas_call_sequenz(havas_best_remainder)
				case 7:
					schema_havas_call_sequenz(havas_euclidean_norm)
				case 8:
					schema_havas_call_sequenz(havas_minimal_norm)
				case 9:
					schema_havas_call_sequenz(havas_sorting_gcd)
				case 10:
					schema_havas_call_sequenz(havas_min_no_of_elements)
				case 11:
					schema_havas_call_sequenz(havas_min_abs_of_row_plus_minimal_euclidean_norm)
				case 12:
					schema_havas_call_sequenz(havas_min_abs_of_row_plus_minimal_norm)
				case 13:
					schema_havas_call_sequenz(havas_min_abs_of_row_plus_min_no_of_elements)
				case 14:
					schema_havas_call_sequenz(havas_minimal_norm_plus_min_abs_of_row)
				case 15:
					schema_havas_call_sequenz(havas_minimal_norm_plus_sorting_gcd)
				case 16:
					schema_havas_call_sequenz(havas_minimal_norm_plus_min_no_of_elements)
				case 17:
					schema_havas_call_sequenz(havas_sorting_gcd_plus_minimal_euclidean_norm)
				case 18:
					schema_havas_call_sequenz(havas_sorting_gcd_plus_minimal_norm)
				case 19:
					schema_havas_call_sequenz(havas_sorting_gcd_plus_min_no_of_elements)
				case 20:
					schema_havas_call_sequenz(havas_min_no_of_elements_plus_min_abs_of_row)
				case 21:
					schema_havas_call_sequenz(havas_min_no_of_elements_plus_minimal_euclidean_norm)
				case 22:
					schema_havas_call_sequenz(havas_min_no_of_elements_plus_minimal_norm)
				case 23:
					schema_havas_call_sequenz(havas_min_no_of_elements_plus_sorting_gcd)
				case 24:
					schema_havas_call_sequenz(Storjohann)
				case 25:
					schema_havas_call_sequenz(Heuristik)
				default:
					lidia_error_handler("bigint_matrix", "hnf_havas :: Mode not supported! (mgcd Modul not defined)");
		}
		break;
	}
	case 1:
	{
		switch(mgcdModul) {
		case 0:
			schema_havas_stf_call_sequenz(hermite)
				case 1:
					schema_havas_stf_call_sequenz(Bradley)
				case 2:
					schema_havas_stf_call_sequenz(Iliopoulos)
				case 3:
					schema_havas_stf_call_sequenz(opt)
				case 4:
					schema_havas_stf_call_sequenz(Blankinship)
				case 5:
					schema_havas_stf_call_sequenz(Best_remainder)
				case 6:
					schema_havas_stf_call_sequenz(havas_best_remainder)
				case 7:
					schema_havas_stf_call_sequenz(havas_euclidean_norm)
				case 8:
					schema_havas_stf_call_sequenz(havas_minimal_norm)
				case 9:
					schema_havas_stf_call_sequenz(havas_sorting_gcd)
				case 10:
					schema_havas_stf_call_sequenz(havas_min_no_of_elements)
				case 11:
					schema_havas_stf_call_sequenz(havas_min_abs_of_row_plus_minimal_euclidean_norm)
				case 12:
					schema_havas_stf_call_sequenz(havas_min_abs_of_row_plus_minimal_norm)
				case 13:
					schema_havas_stf_call_sequenz(havas_min_abs_of_row_plus_min_no_of_elements)
				case 14:
					schema_havas_stf_call_sequenz(havas_minimal_norm_plus_min_abs_of_row)
				case 15:
					schema_havas_stf_call_sequenz(havas_minimal_norm_plus_sorting_gcd)
				case 16:
					schema_havas_stf_call_sequenz(havas_minimal_norm_plus_min_no_of_elements)
				case 17:
					schema_havas_stf_call_sequenz(havas_sorting_gcd_plus_minimal_euclidean_norm)
				case 18:
					schema_havas_stf_call_sequenz(havas_sorting_gcd_plus_minimal_norm)
				case 19:
					schema_havas_stf_call_sequenz(havas_sorting_gcd_plus_min_no_of_elements)
				case 20:
					schema_havas_stf_call_sequenz(havas_min_no_of_elements_plus_min_abs_of_row)
				case 21:
					schema_havas_stf_call_sequenz(havas_min_no_of_elements_plus_minimal_euclidean_norm)
				case 22:
					schema_havas_stf_call_sequenz(havas_min_no_of_elements_plus_minimal_norm)
				case 23:
					schema_havas_stf_call_sequenz(havas_min_no_of_elements_plus_sorting_gcd)
				case 24:
					schema_havas_stf_call_sequenz(Storjohann)
				case 25:
					schema_havas_stf_call_sequenz(Heuristik)
				default:
					lidia_error_handler("bigint_matrix", "hnf_havas :: Mode not supported! (mgcd Modul not defined)");
		}
	}
	}
}



#define schema_havas_call_sequenz_ex(name)                               \
{if (bitfield.get_representation() == matrix_flags::dense_representation)\
{                                                                        \
  havas_kernel< bigint, RODMM< bigint >, RODMM< bigint >, \
    name< bigint, RODMM< bigint >, RODMM< bigint > > > modul; \
  modul.hnf_Z1(*this, T, i, j); \
}                                                                        \
else                                                                     \
{                                                                        \
  havas_kernel< bigint, COSMM< bigint >, COSMM< bigint >, \
    name< bigint, COSMM< bigint >, COSMM< bigint > > > modul; \
                                                                         \
  set_orientation(matrix_flags::column_oriented); \
  modul.hnf_Z1(*this, T, i, j); \
  set_orientation(matrix_flags::row_oriented); \
}                                                                        \
break; }

#define schema_havas_stf_call_sequenz_ex(name)                           \
{if (bitfield.get_representation() == matrix_flags::dense_representation)\
{                                                                        \
  havas_kernel< bigint, RODMM< bigint >, RODMM< bigint >, \
    name< bigint, RODMM< bigint >, RODMM< bigint > > > modul; \
  modul.hnf_Z2(*this, T, i, j); \
  normalization_kernel< bigint, RODMM< bigint >, \
                       RODMM< bigint > > nmodul; \
                                                                         \
  switch(normalizeModul)                                                 \
    {                                                                    \
    case 0:                                                              \
      {                                                                  \
        nmodul.normalize_Std(*this, T, 0, columns - rows); \
        break; \
      }                                                                  \
    case 1:                                                              \
      {                                                                  \
        nmodul.normalize_ChouCollins(*this, T, 0, columns - rows); \
        break; \
      }                                                                  \
    case 2:                                                              \
      {                                                                  \
        nmodul.normalizeHybrid_Std(*this, T, 0, columns - rows); \
        break; \
      }                                                                  \
    case 3:                                                              \
      {                                                                  \
        nmodul.normalizeHybrid_ChouCollins(*this, T, 0, columns - rows); \
        break; \
      }                                                                  \
    default:                                                             \
      lidia_error_handler("bigint_matrix", "hnf_havas :: "               \
                          "Mode not supported! (normalize "              \
                          "Modul not defined)"); \
    }                                                                    \
}                                                                        \
else                                                                     \
{                                                                        \
  havas_kernel< bigint, COSMM< bigint >, COSMM< bigint >, \
    name< bigint, COSMM< bigint >, COSMM< bigint > > > modul; \
                                                                         \
  set_orientation(matrix_flags::column_oriented); \
  modul.hnf_Z2(*this, T, i, j); \
  normalization_kernel< bigint, COSMM< bigint >, \
                       COSMM< bigint > > nmodul; \
                                                                         \
  switch(normalizeModul)                                                 \
    {                                                                    \
    case 0:                                                              \
      {                                                                  \
        nmodul.normalize_Std(*this, T, 0, columns - rows); \
        break; \
      }                                                                  \
    case 1:                                                              \
      {                                                                  \
        nmodul.normalize_ChouCollins(*this, T, 0, columns - rows); \
        break; \
      }                                                                  \
    case 2:                                                              \
      {                                                                  \
        nmodul.normalizeHybrid_Std(*this, T, 0, columns - rows); \
        break; \
      }                                                                  \
    case 3:                                                              \
      {                                                                  \
        nmodul.normalizeHybrid_ChouCollins(*this, T, 0, columns - rows); \
        break; \
      }                                                                  \
    default:                                                             \
      lidia_error_handler("bigint_matrix", "hnf_havas :: "               \
                          "Mode not supported! (normalize "              \
                          "Modul not defined)"); \
    }                                                                    \
  set_orientation(matrix_flags::row_oriented); \
}                                                                        \
break; }

void matrix< bigint >::
hnf_havas(matrix< bigint > & T, lidia_size_t KernAlgo, lidia_size_t mgcdModul, lidia_size_t normalizeModul)
{
	//
	// Task: HNF Computation
	// Algorithm: Gauss with reduction
	// IMPROVEMENTS: Theory of Havas / best reaminder
	// Version: 2.0
	//

	debug_handler_l(DMESSAGE, "hnf_havas(matrix< bigint > &)", DVALUE + 5);

	lidia_size_t i = rows, j = columns;

	if (T.columns != columns)
		T.set_no_of_columns(columns);
	if (T.rows != columns)
		T.set_no_of_rows(columns);

	T.diag(1, 0);

	switch(KernAlgo) {
	case 0:
	{
		switch(mgcdModul) {
		case 0:
			schema_havas_call_sequenz_ex(hermite)
				case 1:
					schema_havas_call_sequenz_ex(Bradley)
				case 2:
					schema_havas_call_sequenz_ex(Iliopoulos)
				case 3:
					schema_havas_call_sequenz_ex(opt)
				case 4:
					schema_havas_call_sequenz_ex(Blankinship)
				case 5:
					schema_havas_call_sequenz_ex(Best_remainder)
				case 6:
					schema_havas_call_sequenz_ex(havas_best_remainder)
				case 7:
					schema_havas_call_sequenz_ex(havas_euclidean_norm)
				case 8:
					schema_havas_call_sequenz_ex(havas_minimal_norm)
				case 9:
					schema_havas_call_sequenz_ex(havas_sorting_gcd)
				case 10:
					schema_havas_call_sequenz_ex(havas_min_no_of_elements)
				case 11:
					schema_havas_call_sequenz_ex(havas_min_abs_of_row_plus_minimal_euclidean_norm)
				case 12:
					schema_havas_call_sequenz_ex(havas_min_abs_of_row_plus_minimal_norm)
				case 13:
					schema_havas_call_sequenz_ex(havas_min_abs_of_row_plus_min_no_of_elements)
				case 14:
					schema_havas_call_sequenz_ex(havas_minimal_norm_plus_min_abs_of_row)
				case 15:
					schema_havas_call_sequenz_ex(havas_minimal_norm_plus_sorting_gcd)
				case 16:
					schema_havas_call_sequenz_ex(havas_minimal_norm_plus_min_no_of_elements)
				case 17:
					schema_havas_call_sequenz_ex(havas_sorting_gcd_plus_minimal_euclidean_norm)
				case 18:
					schema_havas_call_sequenz_ex(havas_sorting_gcd_plus_minimal_norm)
				case 19:
					schema_havas_call_sequenz_ex(havas_sorting_gcd_plus_min_no_of_elements)
				case 20:
					schema_havas_call_sequenz_ex(havas_min_no_of_elements_plus_min_abs_of_row)
				case 21:
					schema_havas_call_sequenz_ex(havas_min_no_of_elements_plus_minimal_euclidean_norm)
				case 22:
					schema_havas_call_sequenz_ex(havas_min_no_of_elements_plus_minimal_norm)
				case 23:
					schema_havas_call_sequenz_ex(havas_min_no_of_elements_plus_sorting_gcd)
				case 24:
					schema_havas_call_sequenz_ex(Storjohann)
				case 25:
					schema_havas_call_sequenz_ex(Heuristik)
				default:
					lidia_error_handler("bigint_matrix", "hnf_havas :: Mode not supported! (mgcd Modul not defined)");
		}
		break;
	}
	case 1:
	{
		switch(mgcdModul) {
		case 0:
			schema_havas_stf_call_sequenz_ex(hermite)
				case 1:
					schema_havas_stf_call_sequenz_ex(Bradley)
				case 2:
					schema_havas_stf_call_sequenz_ex(Iliopoulos)
				case 3:
					schema_havas_stf_call_sequenz_ex(opt)
				case 4:
					schema_havas_stf_call_sequenz_ex(Blankinship)
				case 5:
					schema_havas_stf_call_sequenz_ex(Best_remainder)
				case 6:
					schema_havas_stf_call_sequenz_ex(havas_best_remainder)
				case 7:
					schema_havas_stf_call_sequenz_ex(havas_euclidean_norm)
				case 8:
					schema_havas_stf_call_sequenz_ex(havas_minimal_norm)
				case 9:
					schema_havas_stf_call_sequenz_ex(havas_sorting_gcd)
				case 10:
					schema_havas_stf_call_sequenz_ex(havas_min_no_of_elements)
				case 11:
					schema_havas_stf_call_sequenz_ex(havas_min_abs_of_row_plus_minimal_euclidean_norm)
				case 12:
					schema_havas_stf_call_sequenz_ex(havas_min_abs_of_row_plus_minimal_norm)
				case 13:
					schema_havas_stf_call_sequenz_ex(havas_min_abs_of_row_plus_min_no_of_elements)
				case 14:
					schema_havas_stf_call_sequenz_ex(havas_minimal_norm_plus_min_abs_of_row)
				case 15:
					schema_havas_stf_call_sequenz_ex(havas_minimal_norm_plus_sorting_gcd)
				case 16:
					schema_havas_stf_call_sequenz_ex(havas_minimal_norm_plus_min_no_of_elements)
				case 17:
					schema_havas_stf_call_sequenz_ex(havas_sorting_gcd_plus_minimal_euclidean_norm)
				case 18:
					schema_havas_stf_call_sequenz_ex(havas_sorting_gcd_plus_minimal_norm)
				case 19:
					schema_havas_stf_call_sequenz_ex(havas_sorting_gcd_plus_min_no_of_elements)
				case 20:
					schema_havas_stf_call_sequenz_ex(havas_min_no_of_elements_plus_min_abs_of_row)
				case 21:
					schema_havas_stf_call_sequenz_ex(havas_min_no_of_elements_plus_minimal_euclidean_norm)
				case 22:
					schema_havas_stf_call_sequenz_ex(havas_min_no_of_elements_plus_minimal_norm)
				case 23:
					schema_havas_stf_call_sequenz_ex(havas_min_no_of_elements_plus_sorting_gcd)
				case 24:
					schema_havas_stf_call_sequenz_ex(Storjohann)
				case 25:
					schema_havas_stf_call_sequenz_ex(Heuristik)
				default:
					lidia_error_handler("bigint_matrix", "hnf_havas :: Mode not supported! (mgcd Modul not defined)");
		}
	}
	}
}



//
// mgcd
//

#define schema_mgcd_call_sequenz(name)                                   \
{if (bitfield.get_representation() == matrix_flags::dense_representation)\
{                                                                        \
  havas_kernel< bigint, RODMM< bigint >, RODMM< bigint >, \
    name< bigint, RODMM< bigint >, RODMM< bigint > > > modul; \
  modul.mgcd(*this, i, j); \
}                                                                        \
else                                                                     \
{                                                                        \
  havas_kernel< bigint, COSMM< bigint >, COSMM< bigint >, \
    name< bigint, COSMM< bigint >, COSMM< bigint > > > modul; \
                                                                         \
  set_orientation(matrix_flags::column_oriented); \
  modul.mgcd(*this, i, j); \
  set_orientation(matrix_flags::row_oriented); \
}                                                                        \
break; }



void matrix< bigint >::
mgcd(lidia_size_t mgcdModul)
{
	lidia_size_t i = rows, j = columns;

	switch(mgcdModul) {
	case 0:
		schema_mgcd_call_sequenz(hermite)
			case 1:
				schema_mgcd_call_sequenz(Bradley)
			case 2:
				schema_mgcd_call_sequenz(Iliopoulos)
			case 3:
				schema_mgcd_call_sequenz(opt)
			case 4:
				schema_mgcd_call_sequenz(Blankinship)
			case 5:
				schema_mgcd_call_sequenz(Best_remainder)
			case 6:
				schema_mgcd_call_sequenz(havas_best_remainder)
			case 7:
				schema_mgcd_call_sequenz(havas_euclidean_norm)
			case 8:
				schema_mgcd_call_sequenz(havas_minimal_norm)
			case 9:
				schema_mgcd_call_sequenz(havas_sorting_gcd)
			case 10:
				schema_mgcd_call_sequenz(havas_min_no_of_elements)
			case 11:
				schema_mgcd_call_sequenz(havas_min_abs_of_row_plus_minimal_euclidean_norm)
			case 12:
				schema_mgcd_call_sequenz(havas_min_abs_of_row_plus_minimal_norm)
			case 13:
				schema_mgcd_call_sequenz(havas_min_abs_of_row_plus_min_no_of_elements)
			case 14:
				schema_mgcd_call_sequenz(havas_minimal_norm_plus_min_abs_of_row)
			case 15:
				schema_mgcd_call_sequenz(havas_minimal_norm_plus_sorting_gcd)
			case 16:
				schema_mgcd_call_sequenz(havas_minimal_norm_plus_min_no_of_elements)
			case 17:
				schema_mgcd_call_sequenz(havas_sorting_gcd_plus_minimal_euclidean_norm)
			case 18:
				schema_mgcd_call_sequenz(havas_sorting_gcd_plus_minimal_norm)
			case 19:
				schema_mgcd_call_sequenz(havas_sorting_gcd_plus_min_no_of_elements)
			case 20:
				schema_mgcd_call_sequenz(havas_min_no_of_elements_plus_min_abs_of_row)
			case 21:
				schema_mgcd_call_sequenz(havas_min_no_of_elements_plus_minimal_euclidean_norm)
			case 22:
				schema_mgcd_call_sequenz(havas_min_no_of_elements_plus_minimal_norm)
			case 23:
				schema_mgcd_call_sequenz(havas_min_no_of_elements_plus_sorting_gcd)
			case 24:
				schema_mgcd_call_sequenz(Storjohann)
			case 25:
				schema_mgcd_call_sequenz(Heuristik)
			default:
				lidia_error_handler("bigint_matrix", "mgcd :: Mode not supported! (mgcd Modul not defined)");
	}
}



#define schema_mgcd_call_sequenz_ex(name)                                \
{if (bitfield.get_representation() == matrix_flags::dense_representation)\
{                                                                        \
  havas_kernel< bigint, RODMM< bigint >, RODMM< bigint >, \
    name< bigint, RODMM< bigint >, RODMM< bigint > > > modul; \
  modul.mgcd(*this, T, i, j); \
}                                                                        \
else                                                                     \
{                                                                        \
  havas_kernel< bigint, COSMM< bigint >, COSMM< bigint >, \
    name< bigint, COSMM< bigint >, COSMM< bigint > > > modul; \
                                                                         \
  set_orientation(matrix_flags::column_oriented); \
  modul.mgcd(*this, T, i, j); \
  set_orientation(matrix_flags::row_oriented); \
}                                                                        \
break; }



void matrix< bigint >::
mgcd(matrix< bigint > & T, lidia_size_t mgcdModul)
{
	lidia_size_t i = rows, j = columns;

	if (T.columns != columns)
		T.set_no_of_columns(columns);
	if (T.rows != columns)
		T.set_no_of_rows(columns);

	T.diag(1, 0);

	switch(mgcdModul) {
	case 0:
		schema_mgcd_call_sequenz_ex(hermite)
			case 1:
				schema_mgcd_call_sequenz_ex(Bradley)
			case 2:
				schema_mgcd_call_sequenz_ex(Iliopoulos)
			case 3:
				schema_mgcd_call_sequenz_ex(opt)
			case 4:
				schema_mgcd_call_sequenz_ex(Blankinship)
			case 5:
				schema_mgcd_call_sequenz_ex(Best_remainder)
			case 6:
				schema_mgcd_call_sequenz_ex(havas_best_remainder)
			case 7:
				schema_mgcd_call_sequenz_ex(havas_euclidean_norm)
			case 8:
				schema_mgcd_call_sequenz_ex(havas_minimal_norm)
			case 9:
				schema_mgcd_call_sequenz_ex(havas_sorting_gcd)
			case 10:
				schema_mgcd_call_sequenz_ex(havas_min_no_of_elements)
			case 11:
				schema_mgcd_call_sequenz_ex(havas_min_abs_of_row_plus_minimal_euclidean_norm)
			case 12:
				schema_mgcd_call_sequenz_ex(havas_min_abs_of_row_plus_minimal_norm)
			case 13:
				schema_mgcd_call_sequenz_ex(havas_min_abs_of_row_plus_min_no_of_elements)
			case 14:
				schema_mgcd_call_sequenz_ex(havas_minimal_norm_plus_min_abs_of_row)
			case 15:
				schema_mgcd_call_sequenz_ex(havas_minimal_norm_plus_sorting_gcd)
			case 16:
				schema_mgcd_call_sequenz_ex(havas_minimal_norm_plus_min_no_of_elements)
			case 17:
				schema_mgcd_call_sequenz_ex(havas_sorting_gcd_plus_minimal_euclidean_norm)
			case 18:
				schema_mgcd_call_sequenz_ex(havas_sorting_gcd_plus_minimal_norm)
			case 19:
				schema_mgcd_call_sequenz_ex(havas_sorting_gcd_plus_min_no_of_elements)
			case 20:
				schema_mgcd_call_sequenz_ex(havas_min_no_of_elements_plus_min_abs_of_row)
			case 21:
				schema_mgcd_call_sequenz_ex(havas_min_no_of_elements_plus_minimal_euclidean_norm)
			case 22:
				schema_mgcd_call_sequenz_ex(havas_min_no_of_elements_plus_minimal_norm)
			case 23:
				schema_mgcd_call_sequenz_ex(havas_min_no_of_elements_plus_sorting_gcd)
			case 24:
				schema_mgcd_call_sequenz_ex(Storjohann)
			case 25:
				schema_mgcd_call_sequenz_ex(Heuristik)
			default:
				lidia_error_handler("bigint_matrix", "mgcd :: Mode not supported! (mgcd Modul not defined)");
	}
}



void matrix< bigint >::
normalize(lidia_size_t normalizeModul)
{
	//
	// Task: Normalization
	// Version: 2.0
	//

	debug_handler_l(DMESSAGE, "normalize(lidia_size_t)", DVALUE + 5);

	if (bitfield.get_representation() == matrix_flags::dense_representation) {
		normalization_kernel< bigint, RODMM< bigint >,
			RODMM< bigint > > nmodul;

		switch(normalizeModul) {
		case 0:
		{
			nmodul.normalize_Std(*this, 0, columns - rows);
			break;
		}
		case 1:
		{
			nmodul.normalize_ChouCollins(*this, 0, columns - rows);
			break;
		}
		case 2:
		{
			nmodul.normalizeMod_Std(*this, 0, columns - rows);
			break;
		}
		case 3:
		{
			nmodul.normalizeMod_ChouCollins(*this, 0, columns - rows);
			break;
		}
		case 4:
		{
			nmodul.normalizeHybrid_Std(*this, 0, columns - rows);
			break;
		}
		case 5:
		{
			nmodul.normalizeHybrid_ChouCollins(*this, 0, columns - rows);
			break;
		}
		default:
			lidia_error_handler("bigint_matrix", "normalize :: "
					    "Mode not supported! (normalize "
					    "Modul not defined)");
		}
	}
	else {
		set_orientation(matrix_flags::column_oriented);
		normalization_kernel< bigint, COSMM< bigint >,
			COSMM< bigint > > nmodul;

		switch(normalizeModul) {
		case 0:
		{
			nmodul.normalize_Std(*this, 0, columns - rows);
			break;
		}
		case 1:
		{
			nmodul.normalize_ChouCollins(*this, 0, columns - rows);
			break;
		}
		case 2:
		{
			nmodul.normalizeMod_Std(*this, 0, columns - rows);
			break;
		}
		case 3:
		{
			nmodul.normalizeMod_ChouCollins(*this, 0, columns - rows);
			break;
		}
		case 4:
		{
			nmodul.normalizeHybrid_Std(*this, 0, columns - rows);
			break;
		}
		case 5:
		{
			nmodul.normalizeHybrid_ChouCollins(*this, 0, columns - rows);
			break;
		}
		default:
			lidia_error_handler("bigint_matrix", "normalize :: "
					    "Mode not supported! (normalize "
					    "Modul not defined)");
		}
		set_orientation(matrix_flags::row_oriented);
	}
}



//
// hnf_kannan
//

void matrix< bigint >::
hnf_kannan(lidia_size_t SW)
{
	//
	// Task: HNF Computation
	// Algorithm: Gauss with reduction
	// IMPROVEMENTS: Algorithm of Kannan and Bachem
	// Version: 2.0
	//

	debug_handler_l(DMESSAGE, "hnf_kannan()", DVALUE + 5);

	lidia_size_t i = rows, j = columns;

	switch(SW) {
	case 0:
	{
		if (bitfield.get_representation() == matrix_flags::dense_representation) {
			kannan_kernel< bigint, RODMM< bigint >, RODMM< bigint >,
				Standard_normalization< bigint, RODMM< bigint >, RODMM< bigint > > > modul;

			modul.hnf(*this, i, j);
		}
		else {
			kannan_kernel< bigint, COSMM< bigint >, COSMM< bigint >,
				Standard_normalization< bigint, COSMM< bigint >, COSMM< bigint > > > modul;

			set_orientation(matrix_flags::column_oriented);
			modul.hnf(*this, i, j);
			set_orientation(matrix_flags::row_oriented);
		}
		break;
	}
	case 1:
	{
		if (bitfield.get_representation() == matrix_flags::dense_representation) {
			kannan_kernel< bigint, RODMM< bigint >, RODMM< bigint >,
				ChouCollins_normalization< bigint, RODMM< bigint >, RODMM< bigint > > > modul;

			modul.hnf(*this, i, j);
		}
		else {
			kannan_kernel< bigint, COSMM< bigint >, COSMM< bigint >,
				ChouCollins_normalization< bigint, COSMM< bigint >, COSMM< bigint > > > modul;

			set_orientation(matrix_flags::column_oriented);
			modul.hnf(*this, i, j);
			set_orientation(matrix_flags::row_oriented);
		}
		break;
	}
	case 2:
	{
		if (bitfield.get_representation() == matrix_flags::dense_representation) {
			kannan_kernel< bigint, RODMM< bigint >, RODMM< bigint >,
				Jacobson_normalization< bigint, RODMM< bigint >, RODMM< bigint > > > modul;

			modul.hnf(*this, i, j);
		}
		else {
			kannan_kernel< bigint, COSMM< bigint >, COSMM< bigint >,
				Jacobson_normalization< bigint, COSMM< bigint >, COSMM< bigint > > > modul;

			set_orientation(matrix_flags::column_oriented);
			modul.hnf(*this, i, j);
			set_orientation(matrix_flags::row_oriented);
		}
		break;
	}
	default:
		lidia_error_handler("bigint_matrix", "hnf_kannan :: Mode not supported!");
	}
}



//
// hnf_cg
//

void matrix< bigint >::
hnf_cg(const matrix< long > &B, long BOUND_1, const bigint &BOUND_2, int no)
{
	//
	// modul definitions
	//

	COSMM< long > modul3;

	havas_kernel< long, COSMM< long >, COSMM< bigint >,
		nf_conf3e< long, COSMM< long >, COSMM< bigint > > > hnf_smodul3le;

	normalization_kernel< long, COSMM< long >, COSMM< bigint > > normalize_smodul3le;

	havas_kernel< bigint, RODMM< bigint >, RODMM< bigint >,
		havas_best_remainder_ext< bigint, RODMM< bigint >, RODMM< bigint > > > hnf_dmodul2e;

	havas_kernel< bigint, COSMM< bigint >, COSMM< bigint >,
		havas_best_remainder< bigint, COSMM< bigint >, COSMM< bigint > > > hnf_smodul2;

	normalization_kernel< bigint, COSMM< bigint >, COSMM< bigint > > normalize_smodul2;

	const modular_arithmetic< DRMK< bigint >,
		dense_fp_matrix_kernel< long, MR< long > >,
		dense_fp_matrix_kernel< bigint, MR< bigint > > > Dm_bigint_modul;

	const modular_arithmetic< SRMK< bigint >,
		sparse_fp_matrix_kernel< long, MR< long > >,
		sparse_fp_matrix_kernel< bigint, MR< bigint > > > Sm_bigint_modul;

  //
  // Variables
  //

	timer t;
	bigint H, DET;
	lidia_size_t actual_row = B.rows, actual_column = B.columns;
	lidia_size_t i, j;

  //
  // main algorithm
  //

	lidia_info_handler(t.start_timer(););

	// Change representation of matrix B
	set_storage_mode(matrix_flags::sparse_representation);

	// Change dimensions of member matrix
	resize(B.rows, B.columns);

	//
	// Stage 1: non-modular long
	//

	{
		// Create a copy of the original matrix
		matrix< long > A = B;
		A.set_zero_element(0);

		// change orientation
		A.set_orientation(matrix_flags::column_oriented);

		// Compute the Hadamard Bound
		modul3.hadamard(A, H);
		lidia_info_handler(t.stop_timer();
				   std::cout << std::endl << "Hadamard Bound = " << H << std::endl;
				   std::cout << "Time: " << t << std::endl;
				   t.cont_timer(););

		// STF computation
		if (hnf_smodul3le.hnf_Z2(A, actual_row, actual_column, BOUND_1)) {
			normalize_smodul3le.normalize_ChouCollins(A, 0, actual_column);

			set_orientation(matrix_flags::row_oriented);

			// copy data
			for (i = 0; i < A.columns; i++)
				for (j = 0; j < A.value_counter[i]; j++)
					sto(A.index[i][j], i, bigint(A.value[i][j]));

			lidia_info_handler(t.stop_timer();
					   std::cout << std::endl << "total: " << t << std::endl;);
			return;
		}

		// copy data
		set_orientation(matrix_flags::column_oriented);
		for (i = 0; i < A.columns; i++) {
			lidia_size_t len = A.value_counter[i];
			value[i] = new bigint[len];
			index[i] = new lidia_size_t[len];
			for (j = 0; j < len; j++) {
				value[i][j] = bigint(A.value[i][j]);
				index[i][j] = A.index[i][j];
			}
			value_counter[i] = len;
			allocated[i] = len;
		}
	}

	//
	// Stage 2: non-modular bigint
	//

	lidia_info_handler(t.stop_timer();
			   std::cout << std::endl << "hnf_Z2 (long): " << t << std::endl;
			   t.cont_timer(););

	{
		matrix< bigint > A3(actual_row, actual_column);

		// insert elements
		for (i = 0; i < A3.columns; i++) {
			for (j = 0; j < value_counter[i]; j++)
				A3.value[index[i][j]][i] = value[i][j];
			delete[] value[i];
			delete[] index[i];
			value_counter[i] = 0;
		}

		if (!hnf_dmodul2e.hnf_Z2(A3, actual_row, actual_column, BOUND_2)) {
			lidia_info_handler(t.stop_timer();
					   std::cout << std::endl << "hnf_Z2 (bigint): " << t << std::endl;
					   t.cont_timer(););

			Dm_bigint_modul.latticedet2(A3, DET, H, no);
			D_bigint_modul.hnfmod_dkt_part(A3, DET);
		}

		// insert elements back
		for (i = 0; i < A3.columns; i++) {
			value[i] = new bigint[A3.rows];
			index[i] = new lidia_size_t[A3.rows];
			allocated[i] = A3.rows;

			lidia_size_t k = 0;
			for (j = 0; j < A3.rows; j++) {
				if (A3.value[j][i] != A3.Zero) {
					value[i][k] = A3.value[j][i];
					index[i][k] = j;
					k++;
				}
			}
			value_counter[i] = k;
		}
	}

	// normalize
	normalize_smodul2.normalize_ChouCollins(*this, 0, B.columns - B.rows);

	set_orientation(matrix_flags::row_oriented);
	lidia_info_handler(t.stop_timer();
			   std::cout << std::endl << "total: " << t << std::endl;);
}



void matrix< bigint >::
hnf_cg(const matrix< long > &B, matrix< bigint > &TR, long BOUND_1,
       const bigint &BOUND_2, int no)
{
	//
	// modul definitions
	//

	COSMM< long > modul3;
	COSMM< bigint > modul4;

	const modular_arithmetic< DRMK< bigint >,
		dense_fp_matrix_kernel< long, MR< long > >,
		dense_fp_matrix_kernel< bigint, MR< bigint > > > Dm_bigint_modul;

	const modular_arithmetic< SRMK< bigint >,
		sparse_fp_matrix_kernel< long, MR< long > >,
		sparse_fp_matrix_kernel< bigint, MR< bigint > > > Sm_bigint_modul;



	TR.set_representation(matrix_flags::sparse_representation);

	if (TR.columns != TR.rows || TR.columns != B.columns)
		TR.resize(B.columns, B.columns);

	TR.diag(1, 0);

	timer t;
	lidia_info_handler(t.start_timer(););

	lidia_size_t actual_row, actual_column;

	bigint H;

	set_representation(matrix_flags::sparse_representation);
	resize(B.rows, B.columns);

	lidia_size_t i, j;
	//
	// Stage 1: non-modular long
	//
	{
		matrix< long > A;
		A.set_zero_element(0);

		// Create a copy of the original matrix
		A = B;

		// change orientation
		A.set_orientation(matrix_flags::column_oriented);
		TR.set_orientation(matrix_flags::column_oriented);

		modul3.hadamard(A, H);
		lidia_info_handler(t.stop_timer();
				   std::cout << std::endl << "Hadamard Bound = " << H << std::endl;
				   std::cout << "Time: " << t << std::endl;
				   t.cont_timer(););

		havas_kernel< long, COSMM< long >, COSMM< bigint >,
			nf_conf3e< long, COSMM< long >, COSMM< bigint > > > hnf_smodul3le;

		normalization_kernel< long, COSMM< long >, COSMM< bigint > > normalize_smodul3le;

		if (hnf_smodul3le.hnf_Z2(A, TR, actual_row, actual_column, BOUND_1)) {
			normalize_smodul3le.normalize_Std(A, TR, 0, actual_column);

			set_orientation(matrix_flags::row_oriented);
			TR.set_orientation(matrix_flags::row_oriented);

			// copy data
			for (i = 0; i < A.columns; i++)
				for (j = 0; j < A.value_counter[i]; j++)
					sto(A.index[i][j], i, bigint(A.value[i][j]));

			lidia_info_handler(t.stop_timer();
					   std::cout << std::endl << "total (Phase 1): " << t << std::endl;);
			return;
		}

		// copying matrix
		set_orientation(matrix_flags::column_oriented);

		for (i = 0; i < A.columns; i++) {
			lidia_size_t len = A.value_counter[i];
			value[i] = new bigint[len];
			index[i] = new lidia_size_t[len];
			for (j = 0; j < len; j++) {
				value[i][j] = bigint(A.value[i][j]);
				index[i][j] = A.index[i][j];
			}
			value_counter[i] = len;
			allocated[i] = len;
		}
	}

	//
	// Stage 2: non-modular bigint
	//

	lidia_info_handler(t.stop_timer();
			   std::cout << std::endl << "hnf_Z2 (long): " << t << std::endl;
			   t.cont_timer(););

	{
		matrix< bigint > A3(actual_row, actual_column);

		// insert elements
		for (i = 0; i < A3.columns; i++) {
			for (j = 0; j < value_counter[i]; j++)
				A3.value[index[i][j]][i] = value[i][j];
			delete[] value[i];
			delete[] index[i];
			value_counter[i] = 0;
		}

		havas_kernel< bigint, RODMM< bigint >, RODMM< bigint >,
			havas_best_remainder_ext< bigint, RODMM< bigint >, RODMM< bigint > > > hnf_dmodul2e;

		if (!hnf_dmodul2e.hnf_Z2(A3, actual_row, actual_column, BOUND_2)) {
			lidia_info_handler(t.stop_timer();
					   std::cout << std::endl << "hnf_Z2 (bigint): " << t << std::endl;
					   t.cont_timer(););
			bigint DET;
			Dm_bigint_modul.latticedet2(A3, DET, H, no);
			D_bigint_modul.hnfmod_dkt_part(A3, DET);
		}

		// insert elements back
		for (i = 0; i < A3.columns; i++) {
			value[i] = new bigint[A3.rows];
			index[i] = new lidia_size_t[A3.rows];
			allocated[i] = A3.rows;

			lidia_size_t k = 0;
			for (j = 0; j < A3.rows; j++) {
				if (A3.value[j][i] != A3.Zero) {
					value[i][k] = A3.value[j][i];
					index[i][k] = j;
					k++;
				}
			}
			value_counter[i] = k;
		}
	}

	// normalize
	normalization_kernel< bigint, COSMM< bigint >, COSMM< bigint > > normalize_smodul2;

	normalize_smodul2.normalize_Std(*this, 0, B.columns - B.rows);

	set_orientation(matrix_flags::row_oriented);
	lidia_info_handler(t.stop_timer();
			   std::cout << std::endl << "total: " << t << std::endl;);
}



//
// hnf_gls_solver
//

void  matrix< bigint >::
hnf_gls_solver()
{
	matrix_flags flags;
	flags.set_representation(matrix_flags::sparse_representation);

	//
	// temp variables
	//

	matrix< bigint > C2(columns, columns, flags);
	matrix< bigint > T2(columns, columns, flags);

	bool SW = false;
	bigint DET, DET2;
	bigint H, H1, d;

	matrix< bigint > B, D, B2;
	B.set_representation(matrix_flags::sparse_representation);

	matrix< bigint > CT(rows, rows, flags);
	matrix< bigint > TT(rows, rows, flags);

	lidia_size_t pos = 0;

	math_vector< bigint > b(rows - 1, rows - 1), x(rows, rows);

	//
	// computing the hadamard bound
	//

	hadamard(H);
	std::cout << "hadamard's bound: " << H << "(" << decimal_length(H) << " digits)" << std::endl;

	//
	// computation of the lattice determinante
	//

	det(DET, H);
	std::cout << "determinante: " << DET << "(" << decimal_length(DET) << " digits)" << std::endl;

	DET2 = DET;
	lidia_size_t i, j, l;
	while (abs(DET) != 1) {
		SW = false;
		std::cout << "Rest determinante (" << pos << ") :" << DET << std::endl;

		C2.diag(1, 0);
		T2.diag(1, 0);

		//
		// creating linear system
		//

		b.set_size(rows - pos - 1);
		for (i = pos + 1; i < rows; i++) {
			b[i - pos - 1] = member(i, pos);
			if (b[i - pos - 1] != 0)
				SW = true;
		}

		if (SW) {
			B.resize(rows - pos - 1, columns - pos - 1);
			for (j = pos + 1; j < rows; j++)
				for (i = pos + 1; i < columns; i++)
					B.sto(i - pos - 1, j - pos - 1, member(i, j));

			//
			// new system
			//

			//D.trans(B);

			//B.det(DET2, H);
			//std::cout << "partial determinante (sparse): " << DET2 << std::endl;

			x.set_size(rows - pos);

			//
			// size reduction 1
			//

			//d = b[0];
			// for (i = 1; i < rows - pos - 1; i++)
			//d = gcd(d, b[i]);

			//for (i = 0; i < rows - pos - 1; i++)
			//b[i] /= d;
			// std::cout << "gcd computation (b) " << std::endl;
			// std::cout << "gcd = " << d << std::endl;

			//B = D * B;
			//b = D * b;

			//
			// system solving
			//    PART 1: wiedemann
			//    PART 2: lanczos
			//

			std::cout << "System solving:" << std::endl;
			B.hadamard(H1);

			B.reginvimage(x, b, H1*H1*H, DET2);
			std::cout << "Loesung: " << x << std::endl;

			//
			// size reduction 2
			//

			d = x[0];
			for (i = 1; i < rows - pos; i++)
				d = gcd(d, x[i]);
			std::cout << "gcd computation" << std::endl;
			std::cout << "gcd = " << d << std::endl;

			bigint *x2 = new bigint[rows - pos];
			for (i = 0; i < rows - pos; i++)
				x2[i] = x[i] / d;

			//
			// conditioning matrix
			//

			CT.resize(rows - pos, rows - pos);
			CT.cond_matrix(x2, rows - pos);

			x2 = CT * x2;


			//
			// transformation matrix
			//

			TT.resize(rows - pos, rows - pos);
			TT.diag(1, 0);

			bigint TMP1, TMP2, TMP3;
			TMP3 = xgcd(TMP1, TMP2, x2[0], x2[1]);
			for (i = 0; i < rows - pos; i++)
				TT.sto(i, 0, x2[i]);
			TT.sto(0, 1, -TMP2);
			TT.sto(1, 1, TMP1);
			delete[] x2;

			//
			// update of A and TR
			//

			for (i = 2; i < rows - pos; i++)
				CT.sto(1, i, -CT.member(1, i));
			C2.insert_at(pos, pos, CT, 0, 0, CT.get_no_of_rows(), CT.get_no_of_columns());
			T2.insert_at(pos, pos, TT, 0, 0, TT.get_no_of_rows(), TT.get_no_of_columns());
			multiply(*this, C2);
			multiply(*this, T2);
		}

		std::cout << "OLDdet = " << DET << std::endl;
		if (member(pos, pos) < 0) {
			for (l = 0; l < rows; l++)
				sto(l, pos, -member(l, pos));
		}

		std::cout << "Element: " << member(pos, pos) << std::endl;
		DET = DET / member(pos, pos);

		//
		// final reduction
		//

		std::cout << "NORMALIZE: " << std::endl;
		bigint TMPQ1, TMPQ2;
		for (i = 0; i < pos; i++)
			for (j = i + 1; j < pos; j++) {
				pos_div_rem(TMPQ1, TMPQ2, member(i, j), member(i, i));
				for (l = 0; l < rows; l++)
					sto(l, j, member(l, j) - TMPQ1*member(l, i));
			}
		pos++;
	}
}



void  matrix< bigint >::
hnf_gls_solver(matrix< bigint > &)
{

	matrix_flags flags;
	flags.set_representation(matrix_flags::sparse_representation);

	//
	// temp variables
	//

	matrix< bigint > C(columns, columns, flags);
	matrix< bigint > C2(columns, columns, flags);
	matrix< bigint > T(columns, columns, flags);
	matrix< bigint > T2(columns, columns, flags);

	C.diag(1, 0);
	T.diag(1, 0);

	bool SW = false;
	bigint DET, DET2;
	bigint H, H1, d;

	matrix< bigint > B, D, B2;
	B.set_representation(matrix_flags::sparse_representation);

	matrix< bigint > CT(rows, rows, flags);
	matrix< bigint > TT(rows, rows, flags);

	lidia_size_t pos = 0;

	math_vector< bigint > b(rows - 1, rows - 1), x(rows, rows);

	//
	// computing the hadamard bound
	//

	hadamard(H);
	std::cout << "hadamard's bound: " << H << "(" << decimal_length(H) << " digits)" << std::endl;

	//
	// computation of the lattice determinante
	//

	det(DET, H);
	std::cout << "determinante: " << DET << "(" << decimal_length(DET) << " digits)" << std::endl;

	lidia_size_t i, j, l;
	while (abs(DET) != 1) {
		SW = false;
		std::cout << "rest determinante (" << pos << ") :" << DET << std::endl;

		C2.diag(1, 0);
		T2.diag(1, 0);

		//
		// creating linear system
		//

		b.set_size(rows - pos - 1);
		for (i = pos + 1; i < rows; i++) {
			b[i - pos - 1] = member(i, pos);
			if (b[i - pos - 1] != 0)
				SW = true;
		}

		if (SW) {
			B.resize(rows - pos - 1, columns - pos - 1);
			for (j = pos + 1; j < rows; j++)
				for (i = pos + 1; i < columns; i++)
					B.sto(i - pos - 1, j - pos - 1, member(i, j));

			//
			// new system
			//

			D.trans(B);

			B.det(DET2, H);
			std::cout << "partial determinante (sparse): " << DET2 << std::endl;

			x.set_size(rows - pos);

			//
			// size reduction 1
			//

			d = b[0];
			for (i = 1; i < rows - pos - 1; i++)
				d = gcd(d, b[i]);

			for (i = 0; i < rows - pos - 1; i++)
				b[i] /= d;
			std::cout << "gcd computation (b) " << std::endl;
			std::cout << "gcd = " << d << std::endl;

			B = D * B;
			b = D * b;

			//
			// system solving
			//    PART 1: wiedemann
			//    PART 2: lanczos
			//

			std::cout << "System solving:" << std::endl;
			B.hadamard(H1);

			B.reginvimage(x, b, H1*H1*H, DET2);
			std::cout << "Loesung: " << x << std::endl;

			//
			// size reduction 2
			//

			d = x[0];
			for (i = 1; i < rows - pos; i++)
				d = gcd(d, x[i]);
			std::cout << "gcd computation" << std::endl;
			std::cout << "gcd = " << d << std::endl;

			bigint *x2 = new bigint[rows - pos];
			for (i = 0; i < rows - pos; i++)
				x2[i] = x[i] / d;

			//
			// conditioning matrix
			//

			CT.resize(rows - pos, rows - pos);
			CT.cond_matrix(x2, rows - pos);

			x2 = CT * x2;


			//
			// transformation matrix
			//

			TT.resize(rows - pos, rows - pos);
			TT.diag(1, 0);

			bigint TMP1, TMP2, TMP3;
			TMP3 = xgcd(TMP1, TMP2, x2[0], x2[1]);
			for (i = 0; i < rows - pos; i++)
				TT.sto(i, 0, x2[i]);
			TT.sto(0, 1, -TMP2);
			TT.sto(1, 1, TMP1);
			delete[] x2;

			//
			// update of A and TR
			//

			for (i = 2; i < rows - pos; i++)
				CT.sto(1, i, -CT.member(1, i));
			C2.insert_at(pos, pos, CT, 0, 0, CT.get_no_of_rows(), CT.get_no_of_columns());
			T2.insert_at(pos, pos, TT, 0, 0, TT.get_no_of_rows(), TT.get_no_of_columns());
			multiply(*this, C2);
			multiply(*this, T2);
		}

		std::cout << "OLDdet = " << DET << std::endl;
		if (member(pos, pos) < 0) {
			for (l = 0; l < rows; l++)
				sto(l, pos, -member(l, pos));
		}

		std::cout << "Element: " << member(pos, pos) << std::endl;
		DET = DET / member(pos, pos);

		//
		// final reduction
		//

		std::cout << "NORMALIZE: " << std::endl;
		bigint TMPQ1, TMPQ2;
		for (i = 0; i < pos; i++)
			for (j = i + 1; j < pos; j++) {
				pos_div_rem(TMPQ1, TMPQ2, member(i, j), member(i, i));
				for (l = 0; l < rows; l++)
					sto(l, j, member(l, j) - TMPQ1*member(l, i));
			}
		pos++;
	}
}



//
// Kernel
//

void matrix< bigint >::
kernel1(const matrix< bigint > & A)
{
	//
	// Task: B.kernel1(A);
	//              => The columns of matrix B form a basis
	//              of the kernel of matrix A.
	// Version: 2.0
	//

	debug_handler_l(DMESSAGE, "kernel1(const matrix< bigint > &)", DVALUE + 5);

	const modular_bigint_matrix_algorithms< DRMK< bigint >,
		dense_fp_matrix_kernel< long, MR< long > >,
		dense_fp_matrix_kernel< bigint, MR< bigint > > > Dm_bigint_modul;

	const modular_bigint_matrix_algorithms< SRMK< bigint >,
		sparse_fp_matrix_kernel< long, MR< long > >,
		sparse_fp_matrix_kernel< bigint, MR< bigint > > > Sm_bigint_modul;

	bigint H;

	if (A.bitfield.get_representation() == matrix_flags::dense_representation) {
		D_bigint_modul.hadamard(A, H);
		Dm_bigint_modul.kernel1(*this, A, H);
	}
	else {
		S_bigint_modul.hadamard(A, H);
		Sm_bigint_modul.kernel1(*this, A, H);
	}

}



void matrix< bigint >::
kernel2(const matrix< bigint > & A)
{
	//
	// Task: B.kernel2(A);
	//              => The columns of matrix B form a basis
	//              of the kernel of matrix A.
	// Version: 2.0
	//

	debug_handler_l(DMESSAGE, "kernel2(const matrix< bigint > &)", DVALUE + 5);

	D_bigint_modul.kernel2(*this, A);
}



//
// regular InvImage
//

void matrix< bigint >::
reginvimage(math_vector< bigint > &RES, const math_vector< bigint > & b) const
{
	bigint TMP, H = hadamard();

	sparse_fp_matrix_algorithms< long, MR< long > > long_modul;
	sparse_fp_matrix_algorithms< bigint, MR< bigint > > bigint_modul;

	// Step 2
	long len;
	shift_left(TMP, H, 1);
	bigint *PRIM = get_primes(TMP, bigint(1));
	PRIM[0].longify(len);

	// Step 3
	matrix< bigint > U(columns + 1, static_cast<int>(len));

	bigint MOD;
	long Modlong;

	// Step 3
	lidia_size_t i, j;
	for (i = 1; i <= len; i++) {
		MOD.assign(PRIM[i]);
		if (MOD.bit_length() > bigint::bits_per_digit()) {
			matrix< bigint > B(rows, columns);
			remainder(B, *this, MOD);
			math_vector< bigint > x(columns, columns);

			bigint DET = bigint_modul.det(B, MOD);

			math_vector< bigint > b1(rows, rows);
			for (j = 0; j < rows; j++) {
				LiDIA::remainder(b1[j], b[j], MOD);
				LiDIA::mult_mod(b1[j], b1[j], DET, MOD);
			}

			bigint_modul.lanczos(B, x, b1, MOD);

			U.value[0][i - 1] = -DET;
			for (j = 1; j < columns + 1; j++)
				U.value[j][i - 1].assign(x[j-1]);
		}
		else {
			MOD.longify(Modlong);
			matrix< long > B(rows, columns);
			B.set_zero_element(0);
			remainder(B, *this, Modlong);
			math_vector< long > x(columns, columns);

			long DET = long_modul.det(B, Modlong);
			math_vector< long > b1(rows, rows);
			for (j = 0; j < rows; j++) {
				LiDIA::remainder(b1[j], b[j], Modlong);
				LiDIA::mult_mod(b1[j], b1[j], DET, Modlong);
			}

			long_modul.lanczos(B, x, b1, Modlong);

			U.value[0][i - 1] = -DET;
			for (j = 1; j < columns + 1; j++)
				U.value[j][i - 1].assign(x[j-1]);
		}
	}

	// Step 4,5
	RES.set_size(U.rows);
	for (i = 0; i < columns + 1; i++)
		LiDIA::chinrest(RES[i], U.value[i], PRIM);
}



void matrix< bigint >::
reginvimage(math_vector< bigint > &RES, const math_vector< bigint > & b,
	    const bigint &H, const bigint &DET) const
{
	bigint TMP; //, H = hadamard();

	sparse_fp_matrix_algorithms< long, MR< long > > long_modul;
	sparse_fp_matrix_algorithms< bigint, MR< bigint > > bigint_modul;

	// Step 2
	long len;
	shift_left(TMP, H, 1);
	//bigint *PRIM = get_primes(TMP, bigint(DET));
	bigint *PRIM = get_primes(TMP, bigint(H));
	PRIM[0].longify(len);

	// Step 3
	bigint MOD;
	long Modlong;

	math_vector< bigint > x(columns, columns);
	math_vector< bigint > x2(columns, columns);

	matrix< bigint > Bbigint(rows, columns);
	math_vector< bigint > xbigint(columns, columns);
	math_vector< bigint > b1bigint(rows, rows);


	matrix< long > Blong(rows, columns);
	Blong.set_zero_element(0);
	math_vector< long > xlong(columns, columns);
	math_vector< long > b1long(rows, rows);

	// Step 3
	bigint PROD;
	lidia_size_t i, j;
	for (i = 1; i <= len; i++) {
		MOD.assign(PRIM[i]);
		if (MOD.bit_length() > bigint::bits_per_digit()) {
			remainder(Bbigint, *this, MOD);

			for (j = 0; j < rows; j++)
				LiDIA::mult_mod(b1bigint[j], b[j], DET, MOD);

			lidia_info_handler(std::cout << "- > in lanczos (" << i << ", " << len << ")" << std::endl;);
			bigint_modul.lanczos(Bbigint, xbigint, b1bigint, MOD);

			x2 = xbigint;
		}
		else {
			MOD.longify(Modlong);
			LiDIA::remainder(Blong, *this, Modlong);

			long DETlong;
			best_remainder(DETlong, DET, Modlong);
			for (j = 0; j < rows; j++) {
				LiDIA::best_remainder(b1long[j], b[j], Modlong);
				LiDIA::mult_mod(b1long[j], b1long[j], DETlong, Modlong);
			}

			lidia_info_handler(std::cout << "- > in lanczos (" << i << ", " <<
					   len << ")" << std::endl;);
			long_modul.lanczos(Blong, xlong, b1long, Modlong);

			for (j = 0; j < columns; j++)
				x2[j] = xlong[j];
		}

		if (i == 1) {
			x = x2;
			PROD = MOD;
		}
		else {

			for (j = 0; j < columns; j++) {
				x[j] = chinese_remainder(x[j], PROD, x2[j], MOD);
				best_remainder(x[j], x[j], (PROD*MOD));
			}
			PROD *= MOD;
		}

		if ((*this)*x == (DET*b)) {
			RES[0] = -DET;
			for (j = 1; j <= columns; j++)
				RES[j] = x[j - 1];
			return;
		}
	}
}



void matrix< bigint >::
reginvimage(math_vector< bigint > &RES, const math_vector< bigint > & b,
	    const bigint &H, const bigint &DET, const bigint &MODULUS) const
{
	bigint TMP; //, H = hadamard();

	sparse_fp_matrix_algorithms< long, MR< long > > long_modul;
	sparse_fp_matrix_algorithms< bigint, MR< bigint > > bigint_modul;

	// Step 2
	long len;
	shift_left(TMP, H, 1);
	//bigint *PRIM = get_primes(TMP, bigint(DET));
	bigint *PRIM = get_primes(TMP, bigint(H));
	PRIM[0].longify(len);

	// Step 3
	bigint MOD;
	long Modlong;

	math_vector< bigint > x(columns, columns);
	math_vector< bigint > x2(columns, columns);

	matrix< bigint > Bbigint(rows, columns);
	math_vector< bigint > xbigint(columns, columns);
	math_vector< bigint > b1bigint(rows, rows);


	matrix< long > Blong(rows, columns);
	Blong.set_zero_element(0);
	math_vector< long > xlong(columns, columns);
	math_vector< long > b1long(rows, rows);

	// Step 3
	lidia_size_t i, j;
	bigint PROD;
	for (i = 1; i <= len; i++) {
		MOD.assign(PRIM[i]);
		if (MOD.bit_length() > bigint::bits_per_digit()) {
			remainder(Bbigint, *this, MOD);

			for (j = 0; j < rows; j++)
				LiDIA::mult_mod(b1bigint[j], b[j], DET, MOD);

			lidia_info_handler(std::cout << "- > in lanczos (" << i << ", "
					   << len << ")" << std::endl;);
			bigint_modul.lanczos(Bbigint, xbigint, b1bigint, MOD);

			x2 = xbigint;
		}
		else {
			MOD.longify(Modlong);
			LiDIA::remainder(Blong, *this, Modlong);

			long DETlong;
			best_remainder(DETlong, DET, Modlong);
			for (j = 0; j < rows; j++) {
				LiDIA::best_remainder(b1long[j], b[j], Modlong);
				LiDIA::mult_mod(b1long[j], b1long[j], DETlong, Modlong);
			}

			lidia_info_handler(std::cout << "- > in lanczos (" << i << ", " <<
					   len << ")" << std::endl;);
			long_modul.lanczos(Blong, xlong, b1long, Modlong);

			for (j = 0; j < columns; j++)
				x2[j] = xlong[j];
		}

		if (i == 1) {
			x = x2;
			PROD = MOD;
		}
		else {

			for (j = 0; j < columns; j++) {
				x[j] = chinese_remainder(x[j], PROD, x2[j], MOD);
				best_remainder(x[j], x[j], (PROD*MOD));
			}
			PROD *= MOD;
		}
		bool RET = true;
		math_vector< bigint > c = (*this)*x - DET*b;
		for (j = 0; j < columns; j++)
			if (c[j] != 0)
				RET = false;
		if (RET) {
			RES[0] = -DET;
			for (j = 1; j <= columns; j++)
				RES[j] = x[j - 1];
			return;
		}
	}
}



void matrix< bigint >::
reginvimage_ZmZ(math_vector< bigint > &RES, const math_vector< bigint > & b,
		const bigint &H, bigint &DET) const
{
	bigint TMP; //, H = hadamard();

	sparse_fp_matrix_algorithms< long, MR< long > > long_modul;
	sparse_fp_matrix_algorithms< bigint, MR< bigint > > bigint_modul;

	// Step 2
	long len;
	shift_left(TMP, H, 1);
	bigint *PRIM = get_primes(TMP, bigint(DET));
	PRIM[0].longify(len);

	// Step 3
	bigint MOD, DETtmp;
	long Modlong;

	math_vector< bigint > x(columns, columns);
	math_vector< bigint > x2(columns, columns);

	matrix< bigint > Bbigint(rows, columns);
	math_vector< bigint > xbigint(columns, columns);
	math_vector< bigint > b1bigint(rows, rows);


	matrix< long > Blong(rows, columns);
	Blong.set_zero_element(0);
	math_vector< long > xlong(columns, columns);
	math_vector< long > b1long(rows, rows);

	// Step 3
	lidia_size_t i, j;
	bigint PROD;
	for (i = 1; i <= len; i++) {
		MOD.assign(PRIM[i]);
		if (MOD.bit_length() > bigint::bits_per_digit()) {
			remainder(Bbigint, *this, MOD);

			for (j = 0; j < rows; j++)
				LiDIA::best_remainder(b1bigint[j], b[j], MOD);

			lidia_info_handler(std::cout << "- > in lanczos (" << i << ", " << len << ")" << std::endl;);
			DETtmp = bigint_modul.lanczos_ZmZ(Bbigint, xbigint, b1bigint, MOD);

			std::cout << "Test : " << (Bbigint * xbigint - DETtmp * b1bigint) % MOD << std::endl;
			x2 = xbigint;
		}
		else {
			MOD.longify(Modlong);
			LiDIA::remainder(Blong, *this, Modlong);

			for (j = 0; j < rows; j++)
				LiDIA::best_remainder(b1long[j], b[j], Modlong);

			lidia_info_handler(std::cout << "- > in lanczos (" << i << ", " << len << ")" << std::endl;);
			DETtmp = long_modul.lanczos_ZmZ(Blong, xlong, b1long, Modlong);

			for (j = 0; j < columns; j++)
				x2[j] = xlong[j];
		}

		if (i == 1) {
			x = x2;
			DET = DETtmp;
			PROD = MOD;
		}
		else {

			DET = chinese_remainder(DET, PROD, DETtmp, MOD);
			remainder(DET, DET, PROD*MOD);
			std::cout << "->DET = " << DET << " mod " << PROD*MOD << std::endl;
			for (j = 0; j < columns; j++) {
				x[j] = chinese_remainder(x[j], PROD, x2[j], MOD);
				best_remainder(x[j], x[j], (PROD*MOD));
			}
			PROD *= MOD;
		}

		if ((*this)*x == (DET*b)) {
			RES[0] = -DET;
			for (j = 1; j < columns + 1; j++)
				RES[j] = x[j - 1];
			return;
		}
	}
}



void matrix< bigint >::
reginvimage1(const matrix< bigint > & A, const matrix< bigint > & B)
{
	//
	// Task: C.reginvimage1(A,B);
	//              => A * C.column(j) = g(j)*B.column(j), j=0,...,B.columns
	//              => g(j) minimal
	// Version: 2.0
	//

	debug_handler_l(DMESSAGE, "reginvimage1(const matrix< bigint > &, const matrix< bigint > &",
			DVALUE + 5);

	const modular_bigint_matrix_algorithms< DRMK< bigint >,
		dense_fp_matrix_kernel< long, MR< long > >,
		dense_fp_matrix_kernel< bigint, MR< bigint > > > Dm_bigint_modul;

	const modular_bigint_matrix_algorithms< SRMK< bigint >,
		sparse_fp_matrix_kernel< long, MR< long > >,
		sparse_fp_matrix_kernel< bigint, MR< bigint > > > Sm_bigint_modul;

	bigint H;

	D_bigint_modul.hadamard(A, H);
	Dm_bigint_modul.reginvimage1(*this, A, B, H);
}



void matrix< bigint >::
reginvimage2(const matrix< bigint > & A, const matrix< bigint > & B)
{
	//
	// Task: C.reginvimage2(A,B);
	//              => A * C.column(j) = g(j)*B.column(j), j=0,...,B.columns
	//              => g(j) minimal
	// Version: 2.0
	//

	debug_handler_l(DMESSAGE, "reginvimage2(const matrix< bigint > &, const matrix< bigint > &",
			DVALUE + 5);

	bigint H;

	D_bigint_modul.hadamard(A, H);
	D_bigint_modul.reginvimage2(*this, A, B);
}



//
// Image
//

void matrix< bigint >::
image1(const matrix< bigint > & A)
{
	//
	// Task: B.image1(A);
	//              => The columns of matrix B form a basis
	//                 of the image of matrix A.
	// Version: 2.0
	//

	debug_handler_l(DMESSAGE, "image1(const matrix< bigint > &)", DVALUE + 5);

	const modular_bigint_matrix_algorithms< DRMK< bigint >,
		dense_fp_matrix_kernel< long, MR< long > >,
		dense_fp_matrix_kernel< bigint, MR< bigint > > > Dm_bigint_modul;

	const modular_bigint_matrix_algorithms< SRMK< bigint >,
		sparse_fp_matrix_kernel< long, MR< long > >,
		sparse_fp_matrix_kernel< bigint, MR< bigint > > > Sm_bigint_modul;

	bigint H;

	D_bigint_modul.hadamard(A, H);
	Dm_bigint_modul.image1(*this, A, H);
}



void matrix< bigint >::
image2(const matrix< bigint > & A)
{
	//
	// Task: B.image2(A);
	//              => The columns of matrix B form a basis
	//                 of the image of matrix A.
	// Version: 2.0
	//

	debug_handler_l(DMESSAGE, "image2(const matrix< bigint > &)", DVALUE + 5);

	const modular_bigint_matrix_algorithms< DRMK< bigint >,
		dense_fp_matrix_kernel< long, MR< long > >,
		dense_fp_matrix_kernel< bigint, MR< bigint > > > Dm_bigint_modul;

	const modular_bigint_matrix_algorithms< SRMK< bigint >,
		sparse_fp_matrix_kernel< long, MR< long > >,
		sparse_fp_matrix_kernel< bigint, MR< bigint > > > Sm_bigint_modul;

	Dm_bigint_modul.image2(*this, A);
}



//
// InvImage
//

void matrix< bigint >::
invimage(const matrix< bigint > & B, const bigint * b)
{
	//
	// Task: v = invimage(B,b);
	//              => A is a basis of the solution of B*x=b
	// Version: 2.0
	//

	debug_handler_l(DMESSAGE, "invimage(const matrix< bigint > &, const bigint *)", DVALUE + 5);

	if (b == NULL)
	  precondition_error_handler(PRT, "b", "b != NULL",
				     "void matrix< bigint >::"
				     "invimage(const matrix<bigint>& B,"
				     " const bigint * b)",
				     DMESSAGE, EMESSAGE[1]);

	D_bigint_modul.invimage(*this, B, b);
}



void matrix< bigint >::
invimage(const matrix< bigint > & B, const math_vector< bigint > &b)
{
	//
	// Task: v = invimage(B,b);
	//              => A is a basis of the solution of B*x=b
	// Version: 2.0
	//

	debug_handler_l(DMESSAGE, "invimage(const matrix< bigint > &, const math_vector< bigint > &)", DVALUE + 5);

	if (b.size() != B.rows)
	  precondition_error_handler(b.size(), "b.size", "b.size == B.rows",
				     B.rows, "B.rows", "b.size == B.rows",
				     "void matrix<bigint>::"
				     "invimage(const matrix<bigint>& B,"
				     " const math_vector< bigint > &b)",
				     DMESSAGE, EMESSAGE[1]);

	bigint *tmp = b.get_data_address();
	D_bigint_modul.invimage(*this, B, tmp);
}



//
// Smith normal form
//

void matrix< bigint >::
snf_hartley()
{
	//
	// Task: snf Computation
	// Algorithm: given in Hartley and Hawkes
	// IMPROVEMENTS: Havas, Majewski
	// PAPER: Recognizing badly represented Z-modules, Havas
	// Version: 2.0
	//

	debug_handler_l(DMESSAGE, "snf_hartley()", DVALUE + 5);

	D_bigint_modul.snf_hartley(*this);
}



void matrix< bigint >::
snf_hartley(matrix< bigint > & T1, matrix< bigint > & T2)
{
	//
	// Task: snf Computation
	// Algorithm: given in Hartley and Hawkes
	// PAPER: Recognizing badly represented Z-modules
	// Version: 2.0
	//

	debug_handler_l(DMESSAGE, "snf_hartley(matrix< bigint > &, matrix< bigint > &)", DVALUE + 5);

	if (T1.columns != rows)
		T1.set_no_of_columns(rows);
	if (T1.rows != rows)
		T1.set_no_of_rows(rows);

	if (T2.columns != columns)
		T2.set_no_of_columns(columns);
	if (T2.rows != columns)
		T2.set_no_of_rows(columns);

	T1.diag(1, 0);
	T2.diag(1, 0);

	D_bigint_modul.snf_hartley(*this, T1, T2);
}



void matrix< bigint >::
snf_simple()
{
	//
	// Task: SNF Computation
	// Algorithm: given in Hartley and Hawkes
	// IMPROVEMENT: Havas
	// PAPER: Recognizing badly represented Z-modules
	// Version: 2.0
	//

	debug_handler_l(DMESSAGE, "snf_simple()", DVALUE + 5);

	D_bigint_modul.snf_simple(*this);
}



void matrix< bigint >::
snf_simple(matrix< bigint > & T1, matrix< bigint > & T2)
{
	//
	// Task: SNF Computation
	// Algorithm: given in Hartley and Hawkes
	// IMPROVEMENT: Havas
	// PAPER: Recognizing badly represented Z-modules
	// Version: 2.0
	//

	debug_handler_l(DMESSAGE, "snf_simple(matrix< bigint > &, matrix< bigint > &)", DVALUE + 5);
	if (T1.columns != rows)
		T1.set_no_of_columns(rows);
	if (T1.rows != rows)
		T1.set_no_of_rows(rows);

	if (T2.columns != columns)
		T2.set_no_of_columns(columns);
	if (T2.rows != columns)
		T2.set_no_of_rows(columns);

	T1.diag(1, 0);
	T2.diag(1, 0);

	D_bigint_modul.snf_simple(*this, T1, T2);
}



void matrix< bigint >::
snf_havas()
{
	//
	// Task: snf Computation
	// Algorithm: given in Hartley and Hawkes
	// PAPER: Recognizing badly represented Z-modules
	// IMPROVEMENTS: Havas, best reaminder include mgcd
	// Version: 2.0
	//

	debug_handler_l(DMESSAGE, "snf_havas()", DVALUE + 5);

	D_bigint_modul.snf_havas(*this);
}



void matrix< bigint >::
snf_havas(matrix< bigint > & T1, matrix< bigint > & T2)
{
	//
	// Task: snf Computation
	// Algorithm: given in Hartley and Hawkes
	// PAPER: Recognizing badly represented Z-modules
	// IMPROVEMENTS: Havas, best reaminder include mgcd
	// Version: 2.0
	//

	debug_handler_l(DMESSAGE, "snf_havas(matrix< bigint > &, matrix< bigint > &)", DVALUE + 5);

	if (T1.columns != rows)
		T1.set_no_of_columns(rows);
	if (T1.rows != rows)
		T1.set_no_of_rows(rows);

	if (T2.columns != columns)
		T2.set_no_of_columns(columns);
	if (T2.rows != columns)
		T2.set_no_of_rows(columns);

	T1.diag(1, 0);
	T2.diag(1, 0);

	D_bigint_modul.snf_havas(*this, T1, T2);
}



void matrix< bigint >::
snf_mult(long art)
{
	//
	// Task: snf Computation
	// Algorithm: given in Hartley and Hawkes
	// PAPER: Recognizing badly represented Z-modules + pivot selection
	// Version: 2.0
	//

	debug_handler_l(DMESSAGE, "snf_mult(matrix< bigint > &, matrix< bigint > &)", DVALUE + 5);

	D_bigint_modul.snf_mult(*this, art);
}



void matrix< bigint >::
snf_mult(matrix< bigint > & T1, matrix< bigint > & T2, long art)
{
	//
	// Task: snf Computation
	// Algorithm: given in Hartley and Hawkes
	// PAPER: Recognizing badly represented Z-modules + pivot selection
	// Version: 2.0
	//

	debug_handler_l(DMESSAGE, "snf_mult(matrix< bigint > &, matrix< bigint > &)", DVALUE + 5);

	if (T1.columns != rows)
		T1.set_no_of_columns(rows);
	if (T1.rows != rows)
		T1.set_no_of_rows(rows);

	if (T2.columns != columns)
		T2.set_no_of_columns(columns);
	if (T2.rows != columns)
		T2.set_no_of_rows(columns);

	T1.diag(1, 0);
	T2.diag(1, 0);

	D_bigint_modul.snf_mult(*this, T1, T2, art);
}



void matrix< bigint >::
snf_add(long art)
{
	//
	// Task: snf Computation
	// Algorithm: given in Hartley and Hawkes
	// PAPER: Recognizing badly represented Z-modules + pivot selection
	// Version: 2.0
	//

	debug_handler_l(DMESSAGE, "snf_add(long)", DVALUE + 5);

	D_bigint_modul.snf_add(*this, art);
}



void matrix< bigint >::
snf_add(matrix< bigint > & T1, matrix< bigint > & T2, long art)
{
	//
	// Task: snf Computation
	// Algorithm: given in Hartley and Hawkes
	// PAPER: Recognizing badly represented Z-modules + pivot selection
	// Version: 2.0
	//

	debug_handler_l(DMESSAGE, "snf_add(matrix< bigint > &, matrix< bigint > &, long)", DVALUE + 5);

	if (T1.columns != rows)
		T1.set_no_of_columns(rows);
	if (T1.rows != rows)
		T1.set_no_of_rows(rows);

	if (T2.columns != columns)
		T2.set_no_of_columns(columns);
	if (T2.rows != columns)
		T2.set_no_of_rows(columns);

	T1.diag(1, 0);
	T2.diag(1, 0);

	D_bigint_modul.snf_add(*this, T1, T2, art);
}



void matrix< bigint >::
snf_new(long art)
{
	//
	// Task: snf Computation
	// Algorithm: given in Hartley and Hawkes
	// PAPER: Recognizing badly represented Z-modules + pivot selection
	// Version: 2.0
	//

	debug_handler_l(DMESSAGE, "snf_new(long)", DVALUE + 5);

	D_bigint_modul.snf_new(*this, art);
}



void matrix< bigint >::
snf_new(matrix< bigint > & T1, matrix< bigint > & T2, long art)
{
	//
	// Task: snf Computation
	// Algorithm: given in Hartley and Hawkes
	// PAPER: Recognizing badly represented Z-modules + pivot selection
	// Version: 2.0
	//

	debug_handler_l(DMESSAGE, "snf_new(matrix< bigint > &, matrix< bigint > &, long)", DVALUE + 5);

	if (T1.columns != rows)
		T1.set_no_of_columns(rows);
	if (T1.rows != rows)
		T1.set_no_of_rows(rows);

	if (T2.columns != columns)
		T2.set_no_of_columns(columns);
	if (T2.rows != columns)
		T2.set_no_of_rows(columns);

	T1.diag(1, 0);
	T2.diag(1, 0);

	D_bigint_modul.snf_new(*this, T1, T2, art);
}



void matrix< bigint >::
snfmod_dkt(const bigint &mod)
{
	//
	// Task: A.snfmod_dkt(mod);
	//              => A in Smith normal form
	//              => h = lattice determinant of lattice formed
	//              by the columns of matrix A
	// PAPER: Asymptotically fast triangulazation of matrices over rings
	// IMPROVEMENT: Hafner, McCurly
	// Conditions: rank != rows
	// Version: 2.0
	//

	debug_handler_l(DMESSAGE, "snfmod_dkt(const bigint &)", DVALUE + 5);

	if (rank() != rows)
	  precondition_error_handler(rank(), "rank", "rows == rank",
				     rows, "rows", "rows == rank",
				     "void matrix< bigint >::"
				     "snfmod_dkt(const bigint &mod)",
				     DMESSAGE, EMESSAGE[10]);

	D_bigint_modul.snfmod_dkt(*this, mod);
}



void matrix< bigint >::
snfmod_cohen(const bigint & mod)
{
	//
	// Task: A.snfmod_cohen(mod);
	//              => A in Smith normal form
	//              => mod = multiple of lattice determinant of lattice formed
	//              by the columns of matrix A
	// Conditions: rank != rows
	// Version: 2.0
	//

	debug_handler_l(DMESSAGE, "snfmod_cohen(const bigint &)", DVALUE + 5);

	if (rank() != rows)
	  precondition_error_handler(rank(), "rank", "rank == rows",
				     rows, "rows", "rank == rows",
				     "void matrix< bigint >::"
				     "snfmod_cohen(const bigint & mod)",
				     DMESSAGE, EMESSAGE[10]);

	D_bigint_modul.snfmod_cohen(*this, mod);
}



//
// END: Linear algebra
// PART 2
//

void matrix< bigint >::
gauss()
{
	debug_handler(DMESSAGE, "gauss()");

	D_bigint_modul.gauss(*this);
}



bigint *matrix< bigint >::
mgcd1(const bigint * aconst, lidia_size_t n)
{
	//
	// DESCRIPTION: RES = T.mgcd1(a, n);
	// =  > RES[0] = RES[1]*a[0] + ... + RES[n]*a[n-1]
	// =  > RES[0] = gcd(a[0], ..., a[n-1])
	// =  > T*a = RES
	// ALGORITHM: Blankinship, PIVOT: MINIMUM
	// VERSION: 1.8
	//

	debug_handler("multiple_gcd", "mgcd1(const bigint *, lidia_size_t)");

	register lidia_size_t i, j, index, bound;

	bigint MIN, TMP, q, r, *Ttmp1, *Ttmp2 = NULL;

	if (columns != n)
		set_no_of_columns(n);
	if (rows != n)
		set_no_of_rows(n);

	diag(1, 0);

	bigint *a = new bigint[n + 1];
	memory_handler(a, "multiple_gcd", "mgcd1 :: "
		       "Error in memory allocation (a)");

	for (i = 0; i < n; i++)
		a[i].assign(aconst[i]);

	// Init
	for (index = 0; index < n && a[index].is_zero(); index++);

	if (index == n) {
		delete[] a;
		return new bigint[n];
	}
	else
		bound = index;

	do {
		// Pivot search: MINIMUM
		MIN.assign(a[index]);

		for (i = bound; i < n; i++)
			if ((abs(MIN) > abs(a[i])) && !a[i].is_zero()) {
				MIN.assign(a[i]);
				index = i;
			}

		// first element != 0
		for (i = bound; i < n && (a[i].is_zero() || i == index); i++);
		if (i < n) {
			div_rem(q, r, a[i], MIN);
			a[i].assign(r);
			Ttmp1 = value[i];
			Ttmp2 = value[index];
			for (j = 0; j < n; j++) {
				LiDIA::multiply(TMP, q, Ttmp2[j]);
				LiDIA::subtract(Ttmp1[j], Ttmp1[j], TMP);
			}
		}
	} while (i < n);

	Ttmp2 = value[index];

	// gcd < 0 ?
	if (a[index] < 0) {
		a[index].negate();
		for (i = 0; i < n; i++)
			Ttmp2[i].negate();
	}

	if (index != 0)
		a[0].assign(a[index]);
	for (i = 1; i <= n; i++)
		a[i].assign(Ttmp2[i - 1]);

	return a;
}



bigint *matrix< bigint >::
mgcd2(const bigint * aconst, lidia_size_t n)
{
	//
	// DESCRIPTION: RES = T.mgcd2(a, n);
	// =  > RES[0] = RES[1]*a[0] + ... + RES[n]*a[n-1]
	// =  > RES[0] = gcd(a[0], ..., a[n-1])
	// =  > T*a = RES
	// ALGORITHM: Blankinship
	// IMPROVEMENTS: Havas, Majewski, reduction of all elements, MIN assignments
	// PAPER: Hermite normal form computation for integer matrices, Havas
	// VERSION: 1.8
	//

	debug_handler("multiple_gcd", "mgcd2(const bigint *, lidia_size_t)");

	register lidia_size_t i, j, index, bound, SW;
	bigint MIN, TMP, q, r, *Ttmp1, *Ttmp2 = NULL;

	if (columns != n)
		set_no_of_columns(n);
	if (rows != n)
		set_no_of_rows(n);
	diag(1, 0);

	bigint *a = new bigint[n + 1];
	memory_handler(a, "multiple_gcd", "mgcd2 :: "
		       "Error in memory allocation (a)");

	for (i = 0; i < n; i++)
		a[i].assign(aconst[i]);

  // init
	for (index = 0; index < n && a[index].is_zero(); index++);

	if (index == n) {
		delete[] a;
		return new bigint[n];
	}
	else
		bound = index;

	do {
		MIN.assign(a[index]);

		// Pivot search: MINIMUM
		for (i = bound; i < n; i++)
			if ((abs(MIN) > abs(a[i])) && !a[i].is_zero()) {
				MIN.assign(a[i]);
				index = i;
			}

		// all elements
		SW = 0;
		Ttmp2 = value[index];
		for (i = bound; i < n; i++)
			if ((i != index) && !a[i].is_zero()) {
				SW = 1;
				Ttmp1 = value[i];
				div_rem(q, r, a[i], MIN);
				a[i].assign(r);
				for (j = 0; j < n; j++) {
					LiDIA::multiply(TMP, q, Ttmp2[j]);
					LiDIA::subtract(Ttmp1[j], Ttmp1[j], TMP);
				}
			}
	} while (SW == 1);

	Ttmp2 = value[index];

	// gcd < 0 ?
	if (a[index] < 0) {
		a[index].negate();
		for (i = 0; i < n; i++)
			Ttmp2[i].negate();
	}

	if (index != 0)
		a[0].assign(a[index]);
	for (i = 1; i <= n; i++)
		a[i].assign(Ttmp2[i - 1]);

	return a;
}



void matrix< bigint >::
basis_completion(bigint *v, lidia_size_t n)
{
	//
	//    Task: A.hadamard(H);
	//          => H = hadamard bound of matrix A
	// Version: 2.0
	//

	debug_handler_l(DMESSAGE, "hadamard(bigint &)", DVALUE + 3);

	if (bitfield.get_representation() == matrix_flags::dense_representation)
		D_bigint_modul.basis_completion(*this, v, n);
	else
		S_bigint_modul.basis_completion(*this, v, n);
}



void matrix< bigint >::
simple_basis_completion(bigint *v, lidia_size_t n)
{
	//
	//    Task: A.hadamard(H);
	//          => H = hadamard bound of matrix A
	// Version: 2.0
	//

	debug_handler_l(DMESSAGE, "hadamard(bigint &)", DVALUE + 3);

	if (bitfield.get_representation() == matrix_flags::dense_representation)
		D_bigint_modul.simple_basis_completion(*this, v, n);
	else
		S_bigint_modul.simple_basis_completion(*this, v, n);
}



lidia_size_t matrix< bigint >::
cond_matrix(bigint *v, lidia_size_t n)
{
	//
	//    Task: A.cond_matrix(v, n);
	//
	// Version: 2.0
	//

	debug_handler_l(DMESSAGE, "cond_matrix(bigint *, lidia_size_t)", DVALUE + 3);

	if (bitfield.get_representation() == matrix_flags::dense_representation)
		return D_bigint_modul.cond_matrix(*this, v, n);
	else
		return S_bigint_modul.cond_matrix(*this, v, n);
}



#undef DRMKex
#undef SRMKex

#undef DVALUE
#undef DMESSAGE
#undef EMESSAGE



#ifdef LIDIA_NAMESPACE
}	// end of namespace LiDIA
#endif
