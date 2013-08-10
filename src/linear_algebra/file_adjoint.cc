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
//	Author	: Michael Jacobson (MJ)
//	Changes	: See CVS log
//
//==============================================================================================


#ifdef HAVE_CONFIG_H
# include	"config.h"
#endif
#include	"LiDIA/bigint_matrix.h"
#include	"LiDIA/matrix/dense_fp_matrix_kernel.h"
#include	"LiDIA/matrix/bigint_matrix_algorithms.h"
#include	"LiDIA/matrix/modular_arithmetic.h"
#include	<cstdlib>
#include	<cstdio>
#include	<unistd.h>



#ifdef LIDIA_NAMESPACE
namespace LiDIA {
#endif



file_adjoint::file_adjoint()
{
	static int curr = 0;

	actual_rows = rows = 0;
	actual_rows = columns = 0;
	num_files = 0;
	rows_per_file = 0;
	mods.set_mode(EXPAND);

	index = curr;

	++curr;
}



file_adjoint::~file_adjoint()
{
	kill();
}



void
file_adjoint::kill()
{
	lidia_size_t i;
	char fname[50];

	for (i = 0; i < num_files; ++i) {
		sprintf(fname, "/tmp/%ldQO_ADJ_%i.%i", (long) getpid(), index, i);
		std::remove(fname);
	}

	actual_rows = rows = 0;
	actual_columns = columns = 0;
	num_files = 0;
	rows_per_file = 0;
	mods.set_size(0);
}



void
file_adjoint::init(lidia_size_t nrows, lidia_size_t ncols, lidia_size_t nfiles)
{
	if (nfiles > 0)
		num_files = nfiles;

	if (nrows > 0)
		rows = nrows;

	if (ncols > 0)
		columns = ncols;

	if (rows != columns)
		lidia_error_handler("init", "matrix must be square!");

	rows_per_file = static_cast<lidia_size_t>(std::ceil(static_cast<double>(rows) / static_cast<double>(num_files)));
	actual_rows = rows;
	actual_columns = columns;
}



lidia_size_t
file_adjoint::get_no_of_files()
{
	return num_files;
}



lidia_size_t
file_adjoint::get_no_of_rows()
{
	return actual_rows;
}



lidia_size_t
file_adjoint::get_no_of_columns()
{
	return actual_columns;
}



void
file_adjoint::touch_files()
{
	// touch all temporary files in order to give current date stamps
	register lidia_size_t i;
	char command[100];

	for (i = 0; i < num_files; ++i) {
		sprintf(command, "touch /tmp/%ldQO_ADJ_%i.%i", static_cast<long>(getpid()), index, i);
		std::system(command);
	}

	sprintf(command, "touch /tmp/%ldQO_TRANS*", static_cast<long>(getpid()));
	std::system(command);
}



void
file_adjoint::resize(lidia_size_t nrows, lidia_size_t ncols)
{
	if (nrows != ncols) {
		lidia_error_handler("resize", "must be a square matrix");
		return;
	}

	actual_rows = nrows;
	actual_columns = ncols;
}



void
file_adjoint::add_new_prime(matrix< bigint > & A, bigint & p)
{
	debug_handler("file_adjoint", "add_new_prime");

	register lidia_size_t j, k, curr_rows, nrows;
	lidia_size_t i;
	FILE *out;
	char fname[50];
	matrix< bigint > outmat;
	long val;
	math_vector< bigint > vec;

	vec.set_mode(EXPAND);
	outmat.set_print_mode(LIDIA_MODE);

	// save rows of A to the appropriate files
	curr_rows = 0;
	for (i = 0; i < num_files; ++i) {
		// get size of submatrix
		nrows = rows_per_file;
		if ((curr_rows + nrows) > rows)
			nrows = rows - curr_rows;

		// create submatrix and output
		outmat.resize(nrows, A.columns);
		for (j = 0; j < nrows; ++j) {
			A.get_row_vector(vec, j+curr_rows);
			outmat.sto_row_vector(vec, A.columns, j);
		}

		// get file name
		sprintf(fname, "/tmp/%ldQO_ADJ_%i.%i", static_cast<long>(getpid()), index, i);
		if (mods.size() == 0)
			out = fopen(fname, "wb");
		else
			out = fopen(fname, "ab");
		if (out == NULL) {
			std::cout << "ERROR (add_new_prime)!  ";
			std::cout << "Can't open adjoint file (";
			std::cout << fname << ")!!!" << std::endl;
		}

		// output matrix in binary format
		if (mods.size() == 0) {
			// first matrix - output dimensions
			val = static_cast<long>(outmat.rows);
			fwrite(&val, sizeof(long), 1, out);
			val = static_cast<long>(outmat.columns);
			fwrite(&val, sizeof(long), 1, out);
		}

		for (j = 0; j < outmat.rows; ++j)
			for (k = 0; k < outmat.columns; ++k) {
				outmat.member(j, k).longify(val);
				fwrite(&val, sizeof(long), 1, out);
			}
		fclose(out);

		outmat.reset();
		curr_rows += nrows;
	}

	mods[mods.size()] = p;

	touch_files();
}



bool
file_adjoint::test_adjoint(matrix< bigint > & A, const bigint & DET)
{
	debug_handler("file_adjoint", "test_adjoint");

	register lidia_size_t i, j, k, l, nrows;
	lidia_size_t idx;
	bigint ajk, val;
	matrix< bigint > ADJ;
	bool done;

	// perform test (ADJ*A = det I)
	i = 0;
	idx = 0;
	done = true;
	while ((i < rows) && done) {
		// get next set of rows in adjoint
		combine_rows(ADJ, idx);

		nrows = ADJ.rows;

		for (j = 0; (j < nrows) && done; ++j) {
			for (l = 0; (l < A.columns) && done; ++l) {
				val.assign_zero();
				for (k = 0; k < ADJ.columns; ++k) {
					ajk.assign(ADJ.member(j, k));
					LiDIA::multiply(ajk, ajk, A.member(k, l));
					add(val, val, ajk);
				}
				if ((i+j) == l)
					done = (val == DET);
				else
					done = val.is_zero();
			}
		}

		i += nrows;
		++idx;
	}

	return done;
}



bool
file_adjoint::test_adjoint_diag(matrix< bigint > & A, const bigint & DET)
{
	debug_handler("file_adjoint", "test_adjoint_diag");

	register lidia_size_t i, j, k, l, nrows;
	lidia_size_t idx;
	bigint ajk, val;
	matrix< bigint > ADJ;
	bool done;

	// perform test (ADJ*A = det I)
	i = 0;
	idx = 0;
	done = true;
	while ((i < rows) && done) {
		// get next set of rows in adjoint
		combine_rows(ADJ, idx);

		nrows = ADJ.rows;

		for (j = 0; (j < nrows) && done; ++j) {
			l = i+j;
			val.assign_zero();
			for (k = 0; k < ADJ.columns; ++k) {
				ajk.assign(ADJ.member(j, k));
				LiDIA::multiply(ajk, ajk, A.member(k, l));
				add(val, val, ajk);
			}

			done = (val == DET);

			if (done && (l > 0)) {
				val.assign_zero();
				for (k = 0; k < ADJ.columns; ++k) {
					ajk.assign(ADJ.member(j, k));
					LiDIA::multiply(ajk, ajk, A.member(k, 0));
					add(val, val, ajk);
				}

				done = (val.is_zero());
			}
		}

		i += nrows;
		++idx;
	}

	return done;
}



void
file_adjoint::multiply(math_vector< bigint > & RES, math_vector< bigint > & x)
{
	debug_handler("file_adjoint", "multiply(vector)");

	register lidia_size_t i, j, k, nrows;
	lidia_size_t idx;
	bigint ajk;
	matrix< bigint > ADJ;

	RES.set_mode(EXPAND);

	// perform row by row multiplication
	i = 0;
	idx = 0;
	while (i < rows) {
		// get next set of rows in adjoint
		combine_rows(ADJ, idx);
		ADJ.set_storage_mode(matrix_flags::sparse_representation);
		ADJ.set_orientation(matrix_flags::row_oriented);
		nrows = ADJ.get_no_of_rows();

		for (j = 0; j < nrows; ++j) {
			RES[j+i] = 0;
			for (k = 0; k < ADJ.value_counter[j]; ++k) {
				ajk.assign(ADJ.value[j][k]);
				LiDIA::multiply(ajk, ajk, x[ADJ.index[j][k]]);
				add(RES[j+i], RES[j+i], ajk);
			}
		}

		i += nrows;
		++idx;
	}

	// compute extra entries from unit rows
	for (i = rows; i < actual_rows; ++i)
		RES[i] = x[i];
}



void
file_adjoint::multiply(matrix< bigint > & RES, matrix< bigint > & X)
{
	debug_handler("file_adjoint", "multiply(matrix)");

	register lidia_size_t i, j, k, l, nrows;
	lidia_size_t idx;
	bigint ajk, tmp;
	matrix< bigint > ADJ;

	RES.resize(X.rows, X.columns);

	// perform row by row multiplication
	i = 0;
	idx = 0;
	while (i < rows) {
		// get next set of rows in adjoint
		combine_rows(ADJ, idx);
		ADJ.set_storage_mode(matrix_flags::sparse_representation);
		ADJ.set_orientation(matrix_flags::row_oriented);
		nrows = ADJ.get_no_of_rows();

		for (j = 0; j < nrows; ++j) {
			for (l = 0; l < X.columns; ++l) {
				tmp.assign_zero();
				for (k = 0; k < ADJ.value_counter[j]; ++k) {
					ajk.assign(ADJ.value[j][k]);
					LiDIA::multiply(ajk, ajk, X.member(ADJ.index[j][k], l));
					add(tmp, tmp, ajk);
				}
				RES.sto(j+i, l, tmp);
			}
		}

		i += nrows;
		++idx;
	}

	// compute extra entries from unit rows
	for (i = rows; i < actual_rows; ++i)
		for (j = 0; j < X.columns; ++j)
			RES.sto(i, j, X.member(i, j));
}



void
file_adjoint::combine_rows(matrix< bigint > & RES, lidia_size_t idx)
{
	debug_handler("file_adjoint", "combine_rows");

	register lidia_size_t i, j, k;
	bigint MODULUS, MOD;
	long lrows, lcols, val;
	FILE *in;
	char fname[50];

	const modular_arithmetic< DRMK< bigint >,
		dense_fp_matrix_kernel< long, MR< long > >,
		dense_fp_matrix_kernel< bigint, MR< bigint > > > Dm_bigint_modul;

	matrix< bigint > A;
	A.set_print_mode(LIDIA_MODE);

	// get file name
	sprintf(fname, "/tmp/%ldQO_ADJ_%i.%i", static_cast<long>(getpid()), index, idx);
	if ((in = fopen(fname, "rb")) == NULL) {
		std::cout << "ERROR (combine_rows)!  ";
		std::cout << "Can't open adjoint file (";
		std::cout << fname << ")!!!" << std::endl;
	}

	// input dimensions
	fread(&lrows, sizeof(long), 1, in);
	fread(&lcols, sizeof(long), 1, in);
	A.resize(static_cast<lidia_size_t>(lrows), static_cast<lidia_size_t>(lcols));

	// combine rows for each prime modulus
	MODULUS = 1;
	for (i = 0; i < mods.size(); ++i) {
		MOD.assign(mods[i]);

		// input matrix
		for (j = 0; j < A.rows; ++j)
			for (k = 0; k < A.columns; ++k) {
				fread(&val, sizeof(long), 1, in);
				A.sto(j, k, bigint(val));
			}

		if (i == 0) {
			RES = A;
			MODULUS = MOD;
		}
		else {
			Dm_bigint_modul.chinese_remainder(RES, MODULUS, A, MOD);
			MODULUS *= MOD;
		}
	}

	fclose(in);
}



#ifdef LIDIA_NAMESPACE
}	// end of namespace LiDIA
#endif
