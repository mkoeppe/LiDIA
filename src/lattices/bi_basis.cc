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
//	Author	: Werner Backes (WB), Thorsten Lauer (TL)
//	Changes	: See CVS log
//
//==============================================================================================


#ifdef HAVE_CONFIG_H
# include	"config.h"
#endif
#include	"LiDIA/lattices/bi_lattice_basis.h"
#include	"LiDIA/bigint_matrix.h"



#ifdef LIDIA_NAMESPACE
namespace LiDIA {
#endif



//
// Algorithms
//
void bigint_lattice_basis::lll(sdigit x_factor)
{
	debug_handler("bigint_lattice_basis", "lll(x_factor)");
	if (!(Tr_check_basis()))
		lidia_error_handler("bigint_lattice_basis", "lll(x_factor) :: "
				    "lattice is no basis  !");
	bigint_lattice_basis A;
	Tr_trans_swap(A);
	A.assign_the_rest(*this);
	if (x_factor < 2)
		A.Tr_lll_dbl();
	else
		A.Tr_lll_bfl(x_factor*8*SIZEOF_DOUBLE);
	A.Tr_trans_swap(*this);
	assign_the_rest(A);
}



void bigint_lattice_basis::lll(const math_matrix< bigint > & B, sdigit x_factor)
{
	debug_handler("bigint_lattice_basis", "lll(B, x_factor)");
	((math_matrix< bigint > *)this)->assign(B);
	lll(x_factor);
}



void bigint_lattice_basis::lll_benne_de_weger()
{
	debug_handler("bigint_lattice_basis", "lll_benne_de_weger()");
	if (!(Tr_check_basis()))
		lidia_error_handler("bigint_lattice_basis", "lll_benne_de_weger() :: "
				    "lattice is no basis  !");
	bigint_lattice_basis A;
	Tr_trans_swap(A);
	A.assign_the_rest(*this);
	A.Tr_lll();
	A.Tr_trans_swap(*this);
	assign_the_rest(A);
}



void bigint_lattice_basis::lll_benne_de_weger(const math_matrix< bigint > & B)
{
	debug_handler("bigint_lattice_basis", "lll_benne_de_weger(B)");
	((math_matrix< bigint > *)this)->assign(B);
	lll_benne_de_weger();
}



void bigint_lattice_basis::lll_trans(math_matrix< bigint > & T, sdigit x_factor)
{
	debug_handler("bigint_lattice_basis", "lll_trans(T, x_factor)");
	if (!(Tr_check_basis()))
		lidia_error_handler("bigint_lattice_basis", "lll_trans(T, x_factor) :: "
				    "lattice is no basis  !");
	bigint_lattice_basis A;
	Tr_trans_swap(A);
	A.assign_the_rest(*this);
	T.set_no_of_rows(columns);
	T.set_no_of_columns(columns);
	if (x_factor < 2)
		A.Tr_lll_trans_dbl(T);
	else
		A.Tr_lll_trans_bfl(T, x_factor*8*SIZEOF_DOUBLE);
	A.Tr_trans_swap(*this);
	assign_the_rest(A);
}



void bigint_lattice_basis::lll_trans_benne_de_weger(math_matrix< bigint > & T)
{
	debug_handler("bigint_lattice_basis", "lll_trans_benne_de_weger(T)");
	if (!(Tr_check_basis()))
		lidia_error_handler("bigint_lattice_basis", "lll_trans_benner_de_weger(T) :: "
				    "lattice is no basis  !");
	bigint_lattice_basis A;
	Tr_trans_swap(A);
	A.assign_the_rest(*this);
	T.set_no_of_rows(columns);
	T.set_no_of_columns(columns);
	A.Tr_lll_trans(T);
	A.Tr_trans_swap(*this);
	assign_the_rest(A);
	if (trans_flag)
		LiDIA::multiply(*this, *this, T);
	else
		LiDIA::multiply(*this, T, *this);
}



void bigint_lattice_basis::gram_schmidt_orth(math_matrix< bigfloat > & my, math_matrix< bigfloat > & G)
{
	debug_handler("bigint_lattice_basis", "gram_schmidt_orth(my, G)");
	bigint_lattice_basis Tr;
	bigfloat **Tr_my;
	bigfloat *Tr_my_del;
	bigfloat **Tr_G;
	bigfloat **my_Addr;
	bigfloat **G_Addr;
	lidia_size_t pseudo_rows;
	lidia_size_t pseudo_columns;

	if (!(Tr_check_basis()))
		lidia_error_handler("bigint_lattice_basis", "gram_schmidt_orth(my, G) :: "
				    "lattice is no basis !");
	//
	// allocating memory for matrices for gram - schmitt - orth
	//
	if (trans_flag) {
		pseudo_rows = columns;
		pseudo_columns = rows;
	}
	else {
		pseudo_rows = rows;
		pseudo_columns = columns;
	}
	Tr_my = new bigfloat*[2*pseudo_rows];
	memory_handler(Tr_my, "bigint_lattice_basis", "gram_schmitt_orth(my, G) :: "
		       "not enough memory !");
	Tr_G = &Tr_my[pseudo_rows];
	Tr_my_del = new bigfloat[pseudo_rows*(pseudo_columns+pseudo_rows)];
	memory_handler(Tr_my[0], "bigint_lattice_basis", "gram_schmitt_orth(my, G) :: "
		       "not enough memory !");

	for (lidia_size_t i = 0; i < pseudo_rows; i++) {
		Tr_my[i] = &Tr_my_del[i*pseudo_rows];
		Tr_G[i] = &Tr_my_del[pseudo_rows*pseudo_rows+i*pseudo_columns];
	}
	//
	// Perform gram - schmitt - orth
	//
	Tr_trans_swap(Tr);
	Tr.assign_the_rest(*this);
	Tr.Tr_gram_schmidt_orth(Tr_my, Tr_G);
	Tr.Tr_trans_swap(*this);
	//
	// Do assignments (transposed !!!)
	//
	my.set_no_of_rows(pseudo_rows);
	my.set_no_of_columns(pseudo_rows);
	G.set_no_of_rows(rows);
	G.set_no_of_columns(columns);

	my_Addr = my.get_data_address();
	G_Addr = G.get_data_address();

	for  (lidia_size_t i = 0; i < rows; i++) {
		for (lidia_size_t j = 0; j < columns; j++)
			if (trans_flag)
				LiDIA::swap(G_Addr[i][j], Tr_G[j][i]);
			else
				LiDIA::swap(G_Addr[i][j], Tr_G[i][j]);
	}
	for  (lidia_size_t i = 0; i < pseudo_rows; i++) {
		for (lidia_size_t j = 0; j < pseudo_rows; j++)
			if (trans_flag)
				LiDIA::swap(my_Addr[i][j], Tr_my[i][j]);
			else
				LiDIA::swap(my_Addr[i][j], Tr_my[j][i]);
	}

	//
	// Free allocated storage
	//
	delete[] Tr_my_del;
	delete[] Tr_my;
}



//
// Tools
//
double bigint_lattice_basis::lll_check_search()
{
	debug_handler("bigint_lattice_basis", "lll_check_search()");
	if (!(Tr_check_basis()))
		lidia_error_handler("bigint_lattice_basis", "lll_check_search() :: "
				    "lattice is no basis !");
	bigint_lattice_basis A;
	sdigit a, b;
	Tr_trans_swap(A);
	A.Tr_lll_check_search(a, b);
	A.Tr_trans_swap(*this);
	return(rint((static_cast<double>(a)/static_cast<double>(b))*1000)/1000.0);
}



void bigint_lattice_basis::lll_check_search(sdigit& a, sdigit& b)
{
	debug_handler("bigint_lattice_basis", "lll_check_search(a, b)");
	if (!(Tr_check_basis()))
		lidia_error_handler("bigint_lattice_basis", "lll_check_search(a, b) :: "
				    "lattice is no basis !");
	bigint_lattice_basis A;
	Tr_trans_swap(A);
	A.Tr_lll_check_search(a, b);
	A.Tr_trans_swap(*this);
}



bool bigint_lattice_basis::lll_check(double y)
{
	debug_handler("bigint_lattice_basis", "lll_check(y)");
	if (!(Tr_check_basis()))
		lidia_error_handler("bigint_lattice_basis", "lll_check(y) :: "
				    "lattice is no basis !");
	bigint_lattice_basis A;
	bool flag;
	Tr_trans_swap(A);
	flag = A.Tr_lll_check(y);
	A.Tr_trans_swap(*this);
	return (flag);
}



bool bigint_lattice_basis::lll_check(sdigit a, sdigit b)
{
	debug_handler("bigint_lattice_basis", "lll_check(a, b)");
	if (!(Tr_check_basis()))
		lidia_error_handler("bigint_lattice_basis", "lll_check(a, b) :: "
				    "lattice is no basis !");
	bigint_lattice_basis A;
	bool flag;
	Tr_trans_swap(A);
	flag = A.Tr_lll_check(static_cast<double>(a)/static_cast<double>(b));
	A.Tr_trans_swap(*this);
	return (flag);
}



//
// Conversion
//
void bigint_lattice_basis::extract_basis(const bigint_lattice_gensys& gsys, lidia_size_t& rank)
{
	debug_handler("bigint_lattice_basis", "extract_basis(gsys, rank)");
	bigint_lattice_basis TrThis;
	//
	// Transpose lattice gsys, then basis - checking
	//
	assign(gsys);
	Tr_trans_swap(TrThis);
	TrThis.assign_the_rest(*this);
	TrThis.Tr_extract_basis(rank);
	TrThis.Tr_trans_swap(*this);
}



bool bigint_lattice_basis::make_basis(const bigint_lattice_gensys& gsys, lidia_size_t& rank, sdigit x_factor)
{
	debug_handler("bigint_lattice_basis", "make_basis(gsys, rank, x_factor)");
	bool flag;
	bigint_lattice_basis TrThis;
	//
	// Transpose lattice gsys, then basis - checking
	//
	assign(gsys);
	Tr_trans_swap(TrThis);
	TrThis.assign_the_rest(*this);
	flag = TrThis.Tr_make_basis(rank, x_factor);
	TrThis.Tr_trans_swap(*this);
	assign_the_rest(TrThis);
	//
	// return answer
	//
	return (flag);
}



//
// real implementations
//

//
// basis checking
//
bool bigint_lattice_basis::Tr_make_basis(lidia_size_t& rank, sdigit x_factor)
{
	bigint determ;
	bigint_matrix TrBiM(*this);
	bigint_matrix TrTrBiM(columns, rows);

	determ.assign_zero();
	if (columns == rows)
		determ = det(TrBiM);
	else {
		if (rows < columns) {
			TrTrBiM.trans(TrBiM);
			determ = det(TrBiM*TrTrBiM);
		}
		else
			determ.assign_zero();
	}

	//
	// Conversion to bigint_matrix, then determinaten checking
	//
	if (determ.is_zero()) {
		if (x_factor < 2)
			Tr_lll_dbl_gensys(rank);
		else
			Tr_lll_bfl_gensys(rank, x_factor*8*SIZEOF_DOUBLE);
		set_no_of_rows(rank);
		return (false);
	}
	else {
		rank = rows;
		return (true);
	}
}



void bigint_lattice_basis::Tr_extract_basis(lidia_size_t& rank)
{
	debug_handler("bigint_lattice_basis", "Tr_extract_basis(rank)");
	bool flag;

	rank = rows;
	for (lidia_size_t i = 0; i < rows; i++) {
		flag = false;
		for (lidia_size_t cx = 0; cx < columns; cx++)
			if (value[i][cx].is_zero() == 0) {
				flag = true;
				break;
			}
		if (!flag) {
			--rank;
			for (lidia_size_t j = i; j < rank; j++)
				bin_swap_bin(value[j], value[j+1]);
		}
	}
	set_no_of_rows(rank);
}



//
// Gram - Schmidt ported from bigfloat version
//
bool bigint_lattice_basis::Tr_gram_schmidt_orth(bigfloat** my, bigfloat** gso)
{
	debug_handler("bigint_lattice_basis", "Tr_gram_schmidt_orth(my, gso)");
	sdigit old_prec;
	bigfloat temp0;
	bigfloat temp1;
	bigfloat *tempvect0;
	bigfloat *tempvect1;

	//
	// Save old precision
	//
	old_prec = bigfloat::get_precision();
	bigfloat::set_precision(compute_precision());
	tempvect0 = new bigfloat[columns+rows];
	tempvect1 = &tempvect0[columns];
	vectsize = columns;
	for (lidia_size_t i = 0; i < rows; i++) {
		for (lidia_size_t cx = 0; cx < columns; cx++)
			gso[i][cx].assign(value[i][cx]);
		for (lidia_size_t j = 0; j < i; j++) {
			temp1.assign_zero();
			for (lidia_size_t cx = 0; cx < columns; cx++) {
				temp0.assign(value[i][cx]);
				LiDIA::multiply(temp0, temp0, gso[j][cx]);
				LiDIA::add(temp1, temp1, temp0);
			}
			if (tempvect1[j].is_approx_zero()) {
				//
				// Free allocated storage
				// and restore precision
				//
				delete[] tempvect0;
				bigfloat::set_precision(old_prec);
				return(false);
			}
			LiDIA::divide(my[j][i], temp1, tempvect1[j]);
		}
		for (lidia_size_t j = 0; j < i; j++) {
			bin_scalmul_bfl(tempvect0, my[j][i], gso[j]);
			bin_subtract_bfl(gso[i], gso[i], tempvect0);
		}
		bin_scalprod_bfl(tempvect1[i], gso[i], gso[i]);
	}
	//
	// Free allocated storage
	// and restore precision
	//
	delete[] tempvect0;
	bigfloat::set_precision(old_prec);
	return(true);
}



bool bigint_lattice_basis::Tr_gram_schmidt_orth_bdw(bigint_lattice_basis& my,
                                        	    bigint_lattice_basis& gso,
                                                    bigint* vect)
{
	debug_handler("bigint_lattice_basis", "Tr_gram_schmidt_orth_bdw(my, gso, vect)");

	bigint* tempvect0;
	bigint* tempvect1;

	if (rows > columns)
		return(false);

	tempvect0 = new bigint[rows+columns+1];
	memory_handler(tempvect0, "bigint_lattice_basis", "Tr_gram_schmidt_orth_bdw(my, gso) :: "
		       "not enough memory !");
	tempvect1 = &tempvect0[rows+1];

	vectsize = columns;

	bin_assign_zero_bin(tempvect0);

	tempvect0[0].assign_one();

	for (lidia_size_t i = 0; i < rows; i++) {
		bin_assign_bin(gso.value[i], value[i]);
		for (lidia_size_t j = 0; j < i; j++) {
			bin_scalprod_bin(my.value[j][i], value[i], gso.value[j]);
			for (lidia_size_t cx = 0; cx < columns; cx++) {
				LiDIA::multiply(gso.value[i][cx], tempvect0[j+1], gso.value[i][cx]);
				LiDIA::multiply(tempvect1[cx], my.value[j][i], gso.value[j][cx]);
				LiDIA::subtract(gso.value[i][cx], gso.value[i][cx], tempvect1[cx]);
			}
			if (tempvect0[j].is_zero()) {
				delete[] tempvect0;
				return(false);
			}
			for (lidia_size_t m = 0; m < columns; m++)
				LiDIA::divide(gso.value[i][m], gso.value[i][m], tempvect0[j]);
		}
		bin_scalquad_bin(tempvect0[i+1], gso.value[i]);
		if (tempvect0[i].is_zero()) {
			delete[] tempvect0;
			return(false);
		}
		else
			LiDIA::divide(tempvect0[i+1], tempvect0[i+1], tempvect0[i]);
	}
	vectsize = rows+1;
	bin_assign_bin(vect, tempvect0);
	delete[] tempvect0;
	return(true);
}



bool bigint_lattice_basis::Tr_lll_check(double param)
{
	debug_handler("bigint_lattice_basis", "Tr_lll_check(param)");
	bigint* tempvect0;
	bigint_lattice_basis my(rows, rows);
	bigint_lattice_basis gso(rows, columns);


	tempvect0 = new bigint[rows+1];
	memory_handler(tempvect0, "bigint_lattice_basis", "Tr_lll_check(param) :: "
		       "not enough memory !");

	if (!Tr_gram_schmidt_orth_bdw(my, gso, tempvect0)) {
		delete[] tempvect0;
		lidia_error_handler("bigint_lattice_basis", "Tr_lll_check(param) :: "
				    "lattice is no basis !");
		return(false);
	}

	for (lidia_size_t i = 0; i < rows; i++)
		for (lidia_size_t j = 0; j < i; j++) {
			tempmz0.assign(abs(my[j][i]));
			tempmz0.multiply_by_2();
			if (tempmz0 > tempvect0[j+1]) {
				delete[] tempvect0;
				return (false);
			}
		}
	vectsize = columns;
	for (lidia_size_t i = 1; i < rows; i++) {
		LiDIA::square(tempmz0, tempvect0[i]);
		LiDIA::multiply(tempmz1, tempmz0, static_cast<sdigit>(param*10000000));

		LiDIA::square(tempmz0, my.value[i-1][i]);
		LiDIA::multiply(tempmz2, tempmz0, 10000000);
		LiDIA::subtract(tempmz0, tempmz1, tempmz2);
		LiDIA::multiply(tempmz1, tempvect0[i-1], tempvect0[i+1]);
		LiDIA::multiply(tempmz2, tempmz1, 10000000);

		if (tempmz2 < tempmz0) {
			delete[] tempvect0;
			return(false);
		}

	}
	delete[] tempvect0;
	return (true);
}



void bigint_lattice_basis::Tr_lll_check_search(sdigit& a, sdigit& b)
{
	debug_handler("bigint_lattice_basis", "Tr_lll_check_search(a, b)");
	bool success;
	sdigit param_nom = 3;
	sdigit param_den = 4;
	sdigit upper_nom = 4;
	sdigit upper_den = 4;
	sdigit downer_nom = 2;
	sdigit downer_den = 4;
	bigint* tempvect0;
	bigint_lattice_basis my(rows, rows);
	bigint_lattice_basis gso(rows, columns);


	tempvect0 = new bigint[rows+1];
	memory_handler(tempvect0, "bigint_lattice_basis", "Tr_lll_check_search() :: "
		       "not enough memory !");

	if (!Tr_gram_schmidt_orth_bdw(my, gso, tempvect0)) {
		delete[] tempvect0;
		lidia_error_handler("bigint_lattice_basis", "Tr_lll_check_search() :: "
				    "lattice is no basis !");
		return;
	}

	a = 0;
	b = 1;
	for (lidia_size_t i = 0; i < rows; i++)
		for (lidia_size_t j = 0; j < i; j++) {
			tempmz0.assign(abs(my[j][i]));
			tempmz0.multiply_by_2();
			if (tempmz0 > tempvect0[j+1]) {
				delete[] tempvect0;
				return;
			}
		}
	vectsize = columns;
	for (lidia_size_t j = 0; j < 20; j++)   // 2^-20
	{
		success = true;
		for (lidia_size_t i = 1; i < rows; i++) {
			LiDIA::square(tempmz0, tempvect0[i]);
			LiDIA::multiply(tempmz1, tempmz0, param_nom);

			LiDIA::square(tempmz0, my.value[i-1][i]);
			LiDIA::multiply(tempmz2, tempmz0, param_den);
			LiDIA::subtract(tempmz0, tempmz1, tempmz2);
			LiDIA::multiply(tempmz1, tempvect0[i-1], tempvect0[i+1]);
			LiDIA::multiply(tempmz2, tempmz1, param_den);

			if (tempmz2 < tempmz0) {
				success = false;
				break;
			}

		}
		if (success) {
			a = param_nom;
			b = param_den;

			downer_nom = param_nom;
			downer_den = param_den << 1;
			upper_den <<= 1;
			param_nom = (downer_nom+upper_nom);
			param_den <<= 1;
			upper_nom <<= 1;
			downer_nom <<= 1;
		}
		else {
			upper_nom = param_nom;
			upper_den = param_den << 1;
			downer_den <<= 1;
			param_nom = (upper_nom+downer_nom);
			param_den <<= 1;
			upper_nom <<= 1;
			downer_nom <<= 1;
		}
	}
	delete[] tempvect0;
}



#ifdef LIDIA_NAMESPACE
}	// end of namespace LiDIA
#endif
