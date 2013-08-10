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
#include	<cstdlib>
#include	<cstdio>
#include	<fstream>
#include	<unistd.h>



#ifdef LIDIA_NAMESPACE
namespace LiDIA {
#endif



trans_matrix::trans_matrix()
{
	static int curr = 0;

	files = false;
	size = 0;
	TR.set_mode(EXPAND);
	index = curr;
	last_cols = 0;
	notone_idx = -1;
	notone_rows = 0;
	diag_element = 1;

	adj_files = false;
	adj_index = -1;

	DET_div.assign_zero();
	det_idx = -1;

	++curr;
}


trans_matrix::~trans_matrix()
{
	kill();
}


//
// trans_matrix::kill()
//
// Deletes any temporary files used by the trans_matrix, and reinitializes
// all members.
//

void
trans_matrix::kill()
{
	if (files) {
		lidia_size_t i;
		char fname[50];

		for (i = 0; i < size; ++i) {
			sprintf(fname, "/tmp/%ldQO_TRANS_%d.%d", (long) getpid(), index, i);
			std::remove(fname);
		}
	}

	if (adj_files)
		ADJ.kill();

	size = 0;
	TR.set_size(0);
	last_cols = 0;
	notone_idx = -1;
	notone_rows = 0;
	diag_element = 1;

	adj_index = -1;
}



//
// trans_matrix::set_mode(int)
//
// Sets the storage mode of the trans_matrix:
//   0 - no files
//   1 - components stored in files
//   2 - components stored in files, adjoint also stored in files
//

void
trans_matrix::set_mode(int mode)
{
	files = adj_files = false;

	if (mode > 0)
		files = true;
	if (mode > 1)
		adj_files = true;
}



//
// trans_matrix::get_adjoint_mode()
//
// Returns true if the adjoint is stored in its own set of files, false
// otherwise.
//

bool
trans_matrix::get_adjoint_mode()
{
	return adj_files;
}



//
// trans_matrix::touch_files
//
// Touches all the temporary files of the trans_matrix in order to give
// them current date stamps (i.e., so they aren't deleted by the system
// before we're finished).
//

void
trans_matrix::touch_files()
{
	// touch all temporary files in order to give current date stamps
	if (files) {
		register lidia_size_t i;
		char command[100];

		for (i = 0; i < size; ++i) {
			if (i == adj_index)
				ADJ.touch_files();
			else {
				sprintf(command, "touch /tmp/%ldQO_TRANS_%i.%i", (long)getpid(), index, i);
				std::system(command);
			}
		}
	}
}



//
// trans_matrix::remove_columns(lidia_dize_t *)
//
// Removes the selected columns from the trans_matrix.
//

void
trans_matrix::remove_columns(lidia_size_t *rcols)
{
	if (files) {
		std::ofstream out;
		std::ifstream in;
		char fname[50];
		matrix< bigint > tran;

		tran.set_storage_mode(matrix_flags::sparse_representation);
		tran.set_orientation(matrix_flags::row_oriented);
		tran.set_print_mode(LIDIA_MODE);

		sprintf(fname, "/tmp/%ldQO_TRANS_%d.%d", (long) getpid(), index, size-1);
		in.open(fname);
		if (in.fail()) {
			std::cout << "ERROR (remove_columns)!  ";
			std::cout << "Can't open transformation matrix file for input!!!" << std::endl;
		}
		in >> tran;
		in.close();

		tran.remove_columns(rcols);

		out.open(fname);
		if (out.fail()) {
			std::cout << "ERROR (remove)!  ";
			std::cout << "Can't open transformation matrix file for output!!!" << std::endl;
		}
		out << tran;
		out.close();

		last_cols = tran.columns;
	}
	else {
		TR[size-1].remove_columns(rcols);
		last_cols = TR[size-1].columns;
	}
}



//
// trans_matrix::set_no_of_columns(lidia_size_t)
//
// Sets the number of columns of the trans_matrix, either by adding or
// deleting.
//

void
trans_matrix::set_no_of_columns(lidia_size_t ncols)
{
	lidia_size_t i, j, k, ii, old_rows;
	bigint newel;

	if (size == 0) {
		cols = ncols;
		return;
	}

	if (files) {
		std::ofstream out;
		std::ifstream in;
		char fname[50];
		matrix< bigint > tran;

		tran.set_storage_mode(matrix_flags::sparse_representation);
		tran.set_orientation(matrix_flags::row_oriented);
		tran.set_print_mode(LIDIA_MODE);

		if (ncols > cols) {
			k = ncols - cols;
			for (i = 0; i < size; ++i) {
				if (i == adj_index) {
					ADJ.resize(ncols, ncols);
					last_cols = ncols;
				}
				else {
					sprintf(fname, "/tmp/%ldQO_TRANS_%d.%d", (long) getpid(), index, i);
					in.open(fname);
					if (in.fail()) {
						std::cout << "ERROR (set_no_of_columns)!  ";
						std::cout << "Can't open trans matrix file for input (1)!!!" << std::endl;
					}
					in >> tran;
					in.close();

					old_rows = tran.get_no_of_rows();
					ii = tran.get_no_of_columns();
					if (i == notone_idx)
						newel = diag_element;
					else
						newel = 1;
					tran.resize(old_rows+k, ii+k);
					for (j = old_rows; j < old_rows+k; ++j) {
						tran.sto(j, ii, newel);
						++ii;
					}

					out.open(fname);
					if (out.fail()) {
						std::cout << "ERROR (set_no_of_columns)!  ";
						std::cout << "Can't open trans matrix file for output (1)!!!" << std::endl;
					}
					out << tran;
					out.close();

					last_cols = tran.columns;
				}
			}
		}
		else {
			sprintf(fname, "/tmp/%ldQO_TRANS_%d.%d", (long) getpid(), index, size-1);
			in.open(fname);
			if (in.fail()) {
				std::cout << "ERROR (set_no_of_columns)!  ";
				std::cout << "Can't open trans matrix file for input (2)!!!" << std::endl;
			}
			in >> tran;
			in.close();

			tran.set_no_of_columns(ncols);

			out.open(fname);
			if (out.fail()) {
				std::cout << "ERROR (set_no_of_columns)!  ";
				std::cout << "Can't open trans matrix file for output (2)!!!" << std::endl;
			}
			out << tran;
			out.close();

			last_cols = tran.columns;
		}
	}
	else {
		if (ncols > cols) {
			k = ncols - cols;
			for (i = 0; i < size; ++i) {
				old_rows = TR[i].get_no_of_rows();
				ii = TR[i].get_no_of_columns();
				if (i == notone_idx)
					newel = diag_element;
				else
					newel = 1;
				TR[i].resize(old_rows+k, ii+k);
				for (j = old_rows; j < old_rows+k; ++j) {
					TR[i].sto(j, ii, newel);
					++ii;
				}
			}
		}
		else
			TR[size-1].set_no_of_columns(ncols);

		last_cols = TR[size-1].columns;
	}

	cols = ncols;
}



//
// trans_matrix::get_no_of_columns()
//
// Returns the number of columns
//

lidia_size_t
trans_matrix::get_no_of_columns()
{
	return cols;
}



//
// trans_matrix::get_size()
//
// Returns the numberer of component matrices
//

lidia_size_t
trans_matrix::get_size()
{
	return size;
}



//
// trans_matrix::set_size()
//
// Sets the numberer of component matrices
//

void
trans_matrix::set_size(lidia_size_t nsize)
{
	if (files) {
		lidia_size_t i;
		std::ofstream out;
		char fname[50];

		if (nsize < size) {
			for (i = nsize; i < size; ++i) {
				sprintf(fname, "/tmp/%ldQO_TRANS_%d.%d", (long) getpid(), index, i);
				std::remove(fname);
			}
		}
		else if (nsize > size) {
			for (i = size; i < nsize; ++i) {
				sprintf(fname, "/tmp/%ldQO_TRANS_%d.%d", (long) getpid(), index, i);
				out.open(fname);
				out.close();
			}
		}
	}
	else
		TR.set_size(nsize);

	size = nsize;
}



//
// trans_matrix::get_matrix(matrix <bigint> &, lidia_size_t)
//
// Loads (from memory or disk) the component matrix with index idx into
// the matrix tran
//

void
trans_matrix::get_matrix(matrix< bigint > & tran, lidia_size_t idx)
{
	if (idx >= size) {
		std::cout << "Only " << size << "matrices (input " << idx << ")" << std::endl;
		return;
	}

	unsigned long osmode = tran.get_storage_mode();
	unsigned long oorient = tran.get_orientation();
	unsigned long opmode = tran.get_print_mode();

	if (files) {
		if (idx == adj_index) {
			std::cout << "Can't get adjoint from file (not implemented (or needed)" << std::endl;
			return;
		}

		std::ifstream in;
		char fname[50];

		tran.set_storage_mode(matrix_flags::sparse_representation);
		tran.set_orientation(matrix_flags::row_oriented);
		tran.set_print_mode(LIDIA_MODE);

		sprintf(fname, "/tmp/%ldQO_TRANS_%d.%d", (long) getpid(), index, idx);
		in.open(fname);
		if (in.fail()) {
			std::cout << "ERROR (get_matrix)!  ";
			std::cout << "Can't open transformation matrix file for input!!!" << std::endl;
		}
		in >> tran;
		in.close();

		tran.set_storage_mode(osmode);
		tran.set_orientation(oorient);
		tran.set_print_mode(opmode);
	}
	else
		tran.assign(TR[idx]);
}



//
// trans_matrix::store_matrix(matrix <bigint> &)
//
// Stores the matrix tran as the next component matrix in the trans_matrix,
// either in memory or to disk as indicated by the mode.
//

void
trans_matrix::store_matrix(matrix< bigint > & tran)
{
	std::ofstream out;
	lidia_size_t i, TRcols = tran.columns, TRrows = tran.rows;
	char fname[50];
	unsigned long osmode = tran.get_storage_mode();
	unsigned long oorient = tran.get_orientation();
	unsigned long opmode = tran.get_print_mode();

	tran.set_storage_mode(matrix_flags::sparse_representation);
	tran.set_orientation(matrix_flags::row_oriented);
	if (size == 0) {
		tran.resize(cols, cols);
		for (i = TRcols; i < cols; ++i)
			tran.sto(i, i, bigint(1));
		last_cols = cols;
	}
	else {
		tran.resize(last_cols, last_cols);
		for (i = TRcols; i < last_cols; ++i)
			tran.sto(i, i, bigint(1));
	}

	if (files) {
		tran.set_print_mode(LIDIA_MODE);
		sprintf(fname, "/tmp/%ldQO_TRANS_%d.%d", (long) getpid(), index, size);
		out.open(fname);
		if (out.fail()) {
			std::cout << "ERROR (store_matrix)!  ";
			std::cout << "Can't open transformation matrix file for output!!!" << std::endl;
		}
		out << tran;
		out.close();
	}
	else
		TR[size].assign(tran);

	++size;

	tran.reset();
	tran.set_storage_mode(osmode);
	tran.set_orientation(oorient);
	tran.set_print_mode(opmode);
	tran.resize(TRrows, TRcols);

	touch_files();
}



//
// trans_matrix::store_matrix(matrix <bigint> &, bigint &)
//
// Stores the matrix tran as the next component matrix in the trans_matrix,
// either in memory or to disk as indicated by the mode.  DET is the multiple
// by which the additional diagonal elements must be multiplied in the case
// that AT = D HNF(A).
//

void
trans_matrix::store_matrix(matrix< bigint > & tran, bigint & DET)
{
	std::ofstream out;
	lidia_size_t i, TRcols = tran.columns, TRrows = tran.rows;
	char fname[50];
	unsigned long osmode = tran.get_storage_mode();
	unsigned long oorient = tran.get_orientation();
	unsigned long opmode = tran.get_print_mode();

	notone_idx = size;
	diag_element.assign(DET);

	tran.set_storage_mode(matrix_flags::sparse_representation);
	tran.set_orientation(matrix_flags::row_oriented);
	if (size == 0) {
		tran.resize(cols, cols);
		for (i = TRcols; i < cols; ++i)
			tran.sto(i, i, DET);
		last_cols = cols;
	}
	else {
		tran.resize(last_cols, last_cols);
		for (i = TRcols; i < last_cols; ++i)
			tran.sto(i, i, DET);
	}


	if (files) {
		tran.set_print_mode(LIDIA_MODE);
		sprintf(fname, "/tmp/%ldQO_TRANS_%d.%d", (long) getpid(), index, size);
		out.open(fname);
		if (out.fail()) {
			std::cout << "ERROR (store_matrix_DET)!  ";
			std::cout << "Can't open transformation matrix file!!! for output" << std::endl;
		}
		out << tran;
		out.close();
	}
	else
		TR[size].assign(tran);

	++size;

	tran.reset();
	tran.set_storage_mode(osmode);
	tran.set_orientation(oorient);
	tran.set_print_mode(opmode);
	tran.resize(TRrows, TRcols);

	touch_files();
}



//
// trans_matrix::store_adjoint
//
// Sets the component index of the adjoint
//

void
trans_matrix::store_adjoint()
{
	ADJ.resize(cols, cols);
	adj_index = size;
	++size;
}



//
// trans_matrix::store_det
//
// Sets the component index of the matrix whose diagonal elements must be
// multiplied by DET
//

void
trans_matrix::store_det(bigint & DET)
{
	DET_div.assign(DET);
	det_idx = size;
}



//
// trans_matrix::get_column_vector
//
// Computes column idx of the trans_matrix
//

void
trans_matrix::get_column_vector(math_vector< bigint > & vec, lidia_size_t idx)
{
	debug_handler("trans_matrix", "make_column");

	lidia_size_t i, j, k, last;
	math_vector< bigint > tmp;
	bigint ajk;

	tmp.set_mode(EXPAND);

	if (files) {
		std::ifstream in;
		char fname[50];
		matrix< bigint > tran;

		tran.set_storage_mode(matrix_flags::sparse_representation);
		tran.set_orientation(matrix_flags::row_oriented);
		tran.set_print_mode(LIDIA_MODE);

		// get last matrix and the desired column
		last = size-1;

		sprintf(fname, "/tmp/%ldQO_TRANS_%d.%d", (long) getpid(), index, last);
		in.open(fname);
		if (in.fail()) {
			std::cout << "ERROR (get_column_vector)!  ";
			std::cout << "Can't open trans matrix file for input (1)!!!" << std::endl;
		}
		in >> tran;
		in.close();

		tran.get_column_vector(vec, idx);
		tran.reset();

		// evaluate product
		for (i = last-1; i >= 0; --i) {
			if (i == adj_index) {
				ADJ.multiply(tmp, vec);
				vec = tmp;
			}
			else {
				// get matrix from file
				sprintf(fname, "/tmp/%ldQO_TRANS_%d.%d", (long) getpid(), index, i);
				in.open(fname);
				if (in.fail()) {
					std::cout << "ERROR (get_column_vector)!  ";
					std::cout << "Can't open trans matrix file for input (2)!!!" << std::endl;
				}
				in >> tran;
				in.close();

				// perform non-trivial multiplication (multiply(vec,TR[i],vec))
				for (j = 0; j < tran.rows; ++j) {
					tmp[j] = 0;
					for (k = 0; k < tran.value_counter[j]; ++k) {
						ajk.assign(tran.value[j][k]);
						LiDIA::multiply(ajk, ajk, vec[tran.index[j][k]]);
						add(tmp[j], tmp[j], ajk);
					}
				}

				vec = tmp;
				tran.reset();
			}

			if (i == det_idx)
				LiDIA::divide(vec, vec, DET_div);
		}
	}
	else {
		last = size-1;
		TR[last].get_column_vector(vec, idx);

		for (i = last-1; i >= 0; --i) {
			// perform non-trivial multiplication (multiply(vec,TR[i],vec))
			for (j = 0; j < TR[i].rows; ++j) {
				tmp[j] = 0;
				for (k = 0; k < TR[i].value_counter[j]; ++k) {
					ajk.assign(TR[i].value[j][k]);
					LiDIA::multiply(ajk, ajk, vec[TR[i].index[j][k]]);
					add(tmp[j], tmp[j], ajk);
				}
			}

			vec = tmp;
			if (i == det_idx)
				LiDIA::divide(vec, vec, DET_div);
		}
	}
}



//
// trans_matrix::get_submatrix
//
// Computes the submatrix containing the given columns of the trans_matrix
//

void
trans_matrix::get_submatrix(matrix< bigint > & V, lidia_size_t *idx)
{
	debug_handler("trans_matrix", "get_submatrix");

	if (!idx)
		return;

	lidia_size_t i, j, k, l, last;
	math_vector< bigint > vec, tmp;
	bigint ajk, val;

	vec.set_mode(EXPAND);
	tmp.set_mode(EXPAND);

	if (files) {
		std::ifstream in;
		char fname[50];
		matrix< bigint > tran;

		tran.set_storage_mode(matrix_flags::sparse_representation);
		tran.set_orientation(matrix_flags::row_oriented);
		tran.set_print_mode(LIDIA_MODE);

		// get last matrix and the desired column
		last = size-1;

		sprintf(fname, "/tmp/%ldQO_TRANS_%d.%d", (long) getpid(), index, last);
		in.open(fname);
		if (in.fail()) {
			std::cout << "ERROR (get_submatrix)!  ";
			std::cout << "Can't open trans matrix file for input (1)!!!" << std::endl;
		}
		in >> tran;
		in.close();

		V.resize(tran.rows, idx[0]);
		for (i = 0; i < idx[0]; ++i) {
			tran.get_column_vector(vec, idx[i+1]);
			V.sto_column_vector(vec, tran.rows, i);
		}
		tran.reset();

		// evaluate product
		for (i = last-1; i >= 0; --i) {
			if (i == adj_index) {
				matrix< bigint > B;
				ADJ.multiply(B, V);
				V = B;
			}
			else {
				// get matrix from file
				sprintf(fname, "/tmp/%ldQO_TRANS_%d.%d", (long) getpid(), index, i);
				in.open(fname);
				if (in.fail()) {
					std::cout << "ERROR (get_submatrix)!  ";
					std::cout << "Can't open trans matrix file for input (2)!!!" << std::endl;
				}
				in >> tran;
				in.close();

				// perform non-trivial multiplication (multiply(vec,TR[i],vec))
				for (l = 0; l < V.columns; ++l) {
					V.get_column_vector(vec, l);

					for (j = 0; j < tran.rows; ++j) {
						tmp[j] = 0;
						for (k = 0; k < tran.value_counter[j]; ++k) {
							ajk.assign(tran.value[j][k]);
							LiDIA::multiply(ajk, ajk, vec[tran.index[j][k]]);
							add(tmp[j], tmp[j], ajk);
						}
					}

					V.sto_column_vector(tmp, V.rows, l);
				}

				tran.reset();
			}

			if (i == det_idx)
				LiDIA::divide(V, V, DET_div);
		}
	}
	else {
		last = size-1;
		V.resize(TR[last].rows, idx[0]);
		for (i = 0; i < idx[0]; ++i) {
			TR[last].get_column_vector(vec, idx[i+1]);
			V.sto_column_vector(vec, TR[last].rows, i);
		}

		for (i = last-1; i >= 0; --i) {
			// perform non-trivial multiplication (multiply(V,TR[i],V))
			for (l = 0; l < V.columns; ++l) {
				V.get_column_vector(vec, l);

				for (j = 0; j < TR[i].rows; ++j) {
					tmp[j] = 0;
					for (k = 0; k < TR[i].value_counter[j]; ++k) {
						ajk.assign(TR[i].value[j][k]);
						LiDIA::multiply(ajk, ajk, vec[TR[i].index[j][k]]);
						add(tmp[j], tmp[j], ajk);
					}
				}

				V.sto_column_vector(tmp, V.rows, l);
			}

			if (i == det_idx)
				LiDIA::divide(V, V, DET_div);
		}
	}
}



//
// trans_matrix::multiply()
//
// Computes the product of the trans_matrix with the vector x, i.e., Tx.
//

math_vector< bigint >
trans_matrix::multiply(math_vector< bigint > & x)
{
	debug_handler("trans_matrix", "multiply");

	lidia_size_t i, j, k, last;
	math_vector< bigint > tmp, vec;
	bigint ajk;

	vec.set_mode(EXPAND);
	tmp.set_mode(EXPAND);
	vec = x;

	if (files) {
		std::ifstream in;
		char fname[50];
		matrix< bigint > tran;

		tran.set_storage_mode(matrix_flags::sparse_representation);
		tran.set_orientation(matrix_flags::row_oriented);
		tran.set_print_mode(LIDIA_MODE);

		// get last matrix
		last = size;

		// evaluate product
		for (i = last-1; i >= 0; --i) {
			if (i == adj_index) {
				ADJ.multiply(tmp, vec);
				vec = tmp;
			}
			else {
				// get matrix from file
				sprintf(fname, "/tmp/%ldQO_TRANS_%d.%d", (long) getpid(), index, i);
				in.open(fname);
				if (in.fail()) {
					std::cout << "ERROR (multiply)!  ";
					std::cout << "Can't open transformation matrix file for input!!!" << std::endl;
				}
				in >> tran;
				in.close();

				// perform non-trivial multiplication (multiply(vec,TR[i],vec))
				for (j = 0; j < tran.rows; ++j) {
					tmp[j] = 0;
					for (k = 0; k < tran.value_counter[j]; ++k) {
						ajk.assign(tran.value[j][k]);
						LiDIA::multiply(ajk, ajk, vec[tran.index[j][k]]);
						add(tmp[j], tmp[j], ajk);
					}
				}

				vec = tmp;
				tran.reset();
			}

			if (i == det_idx)
				LiDIA::divide(vec, vec, DET_div);
		}
	}
	else {
		last = size;

		for (i = last-1; i >= 0; --i) {
			// perform non-trivial multiplication (multiply(vec,TR[i],vec))
			for (j = 0; j < TR[i].rows; ++j) {
				tmp[j] = 0;
				for (k = 0; k < TR[i].value_counter[j]; ++k) {
					ajk.assign(TR[i].value[j][k]);
					LiDIA::multiply(ajk, ajk, vec[TR[i].index[j][k]]);
					add(tmp[j], tmp[j], ajk);
				}
			}

			vec = tmp;

			if (i == det_idx)
				LiDIA::divide(vec, vec, DET_div);
		}
	}

	return vec;
}



#ifdef LIDIA_NAMESPACE
}	// end of namespace LiDIA
#endif
