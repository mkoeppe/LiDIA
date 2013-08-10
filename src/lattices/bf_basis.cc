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
#include	"LiDIA/lattices/bf_lattice_basis.h"
#include	"LiDIA/lattices/bi_lattice_basis.h"
#include	"LiDIA/bigint_matrix.h"



#ifdef LIDIA_NAMESPACE
namespace LiDIA {
#endif



//
// Algorithms
//
void bigfloat_lattice_basis::lll(sdigit x_factor)
{
	debug_handler("bigfloat_lattice_basis", "lll(x_factor)");
	if (Tr_check_basis() == false)
		lidia_error_handler("bigfloat_lattice_basis", "lll(x_factor) :: "
				    "lattice is no basis !");
	bigfloat_lattice_basis A;
	Tr_trans_swap(A);
	A.assign_the_rest(*this);
	if (x_factor < 2)
		A.Tr_lll_dbl();
	else
		A.Tr_lll_bfl(x_factor*8*SIZEOF_DOUBLE);
	A.Tr_trans_swap(*this);
	assign_the_rest(A);
}



void bigfloat_lattice_basis::lll(const math_matrix< bigfloat > & B, sdigit x_factor)
{
	debug_handler("bigfloat_lattice_basis", "lll(B, x_factor)");
	((math_matrix< bigfloat > *)this)->assign(B);
	lll(x_factor);
}



void bigfloat_lattice_basis::lll_trans(math_matrix< bigint > & T, sdigit x_factor)
{
	debug_handler("bigfloat_lattice_basis", "lll_trans(T, x_factor)");
	if (!(Tr_check_basis()))
		lidia_error_handler("bigfloat_lattice_basis", "lll_trans(T, x_factor) :: "
				    "lattice is no basis !");
	bigfloat_lattice_basis A;
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



void bigfloat_lattice_basis::lll_trans(math_matrix< bigfloat > & T, sdigit x_factor)
{
	debug_handler("bigfloat_lattice_basis", "lll_trans(T)");
	if (!(Tr_check_basis()))
		lidia_error_handler("bigfloat_lattice_basis", "lll_trans(T) :: "
				    "lattice is no basis !");
	bigfloat_lattice_basis A;
	bigint_lattice_basis TBi;
	Tr_trans_swap(A);
	A.assign_the_rest(*this);
	T.set_no_of_rows(columns);
	T.set_no_of_columns(columns);
	if (x_factor < 2)
		A.Tr_lll_trans_dbl(TBi);
	else
		A.Tr_lll_trans_bfl(TBi, x_factor*8*SIZEOF_DOUBLE);
	A.Tr_trans_swap(*this);
	assign_the_rest(A);
	//
	// Conversion
	//
	T.set_no_of_rows(TBi.rows);
	T.set_no_of_columns(TBi.columns);
	bigfloatify(T, TBi);
}



void bigfloat_lattice_basis::lll_var(sdigit x_factor)
{
	debug_handler("bigfloat_lattice_basis", "lll_var(x_factor)");
	if (!(Tr_check_basis()))
		lidia_error_handler("bigfloat_lattice_basis", "lll_var(x_factor) :: "
				    "lattice is no basis !");
	bigfloat_lattice_basis A;
	bigint_lattice_basis BiA;
	sdigit bit_len;

	Tr_trans_swap(A);
	A.assign_the_rest(*this);
	BiA.bigintify(bit_len, A);
	if (bit_len == 1) {
		BiA.y_par = A.y_par;
		if (x_factor < 2)
			BiA.Tr_lll_dbl();
		else
			BiA.Tr_lll_bfl(x_factor*8*SIZEOF_DOUBLE);
		bigfloatify(A, BiA);
		A.reduction_steps = BiA.reduction_steps;
		A.swapping_steps = BiA.swapping_steps;
		A.correction_steps = BiA.correction_steps;
	}
	else
		if (x_factor < 2)
			A.Tr_lll_var_dbl(BiA, bit_len);
		else
			A.Tr_lll_var_bfl(BiA, bit_len, x_factor*8*SIZEOF_DOUBLE);
	A.Tr_trans_swap(*this);
	assign_the_rest(A);
}



void bigfloat_lattice_basis::lll_var(const math_matrix< bigfloat > & B, sdigit x_factor)
{
	debug_handler("bigfloat_lattice_basis", "lll_var(B, x_factor)");
	((math_matrix< bigfloat > *)this)->assign(B);
	lll_var(x_factor);
}



void bigfloat_lattice_basis::lll_trans_var(math_matrix< bigint > & T, sdigit x_factor)
{
	debug_handler("bigfloat_lattice_basis", "lll_trans_var(T, x_factor)");
	if (!(Tr_check_basis()))
		lidia_error_handler("bigfloat_lattice_basis", "lll_trans_var(T, x_factor) :: "
				    "lattice is no basis !");
	bigfloat_lattice_basis A;
	bigint_lattice_basis BiA;
	sdigit bit_len;

	Tr_trans_swap(A);
	A.assign_the_rest(*this);
	BiA.bigintify(bit_len, A);
	if (bit_len == 1) {
		BiA.y_par = A.y_par;
		if (x_factor < 2)
			BiA.Tr_lll_trans_dbl(T);
		else
			BiA.Tr_lll_trans_bfl(T, x_factor*8*SIZEOF_DOUBLE);
		bigfloatify(A, BiA);
		A.reduction_steps = BiA.reduction_steps;
		A.swapping_steps = BiA.swapping_steps;
		A.correction_steps = BiA.correction_steps;
	}
	else
		if (x_factor < 2)
			A.Tr_lll_trans_var_dbl(BiA, T, bit_len);
		else
			A.Tr_lll_trans_var_bfl(BiA, T, bit_len, x_factor*8*SIZEOF_DOUBLE);
	A.Tr_trans_swap(*this);
	assign_the_rest(A);
}



void bigfloat_lattice_basis::lll_trans_var(math_matrix< bigfloat > & T, sdigit x_factor)
{
	debug_handler("bigfloat_lattice_basis", "lll_trans_var(T, x_factor)");
	if (!(Tr_check_basis()))
		lidia_error_handler("bigfloat_lattice_basis", "lll_trans_var(T, x_factor) :: "
				    "lattice is no basis !");
	bigfloat_lattice_basis A;
	bigint_lattice_basis BiA;
	bigint_lattice_basis TBi;
	sdigit bit_len;

	Tr_trans_swap(A);
	A.assign_the_rest(*this);
	BiA.bigintify(bit_len, A);
	if (bit_len == 1) {
		BiA.y_par = A.y_par;
		if (x_factor < 2)
			BiA.Tr_lll_trans_dbl(TBi);
		else
			BiA.Tr_lll_trans_bfl(TBi, x_factor*8*SIZEOF_DOUBLE);
		bigfloatify(A, BiA);
		A.reduction_steps = BiA.reduction_steps;
		A.swapping_steps = BiA.swapping_steps;
		A.correction_steps = BiA.correction_steps;
	}
	else
		if (x_factor < 2)
			A.Tr_lll_trans_var_dbl(BiA, TBi, bit_len);
		else
			A.Tr_lll_trans_var_bfl(BiA, TBi, bit_len, x_factor*8*SIZEOF_DOUBLE);
	A.Tr_trans_swap(*this);
	assign_the_rest(A);
	//
	// Conversion
	//
	T.set_no_of_rows(TBi.rows);
	T.set_no_of_columns(TBi.columns);
	bigfloatify(T, TBi);
}



void bigfloat_lattice_basis::gram_schmidt_orth(math_matrix< bigfloat > & my, math_matrix< bigfloat > & G)
{
	debug_handler("bigfloat_lattice_basis", "gram_schmidt_orth(my, G)");
	bigfloat_lattice_basis Tr;
	bigfloat **Tr_my;
	bigfloat *Tr_my_del;
	bigfloat **Tr_G;
	bigfloat **G_Addr;
	bigfloat **my_Addr;
	lidia_size_t pseudo_rows;
	lidia_size_t pseudo_columns;


	if (!Tr_check_basis())
		lidia_error_handler("bigfloat_lattice_basis", "gram_schmidt_orth(my, G) :: "
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
	memory_handler(Tr_my, "bigfloat_lattice_basis", "gram_schmidt_orth(my, G) :: "
		       "not enough memory !");
	Tr_G = &Tr_my[pseudo_rows];
	Tr_my_del = new bigfloat[pseudo_rows*(pseudo_columns+pseudo_rows)];
	memory_handler(Tr_my[0], "bigfloat_lattice_basis", "gram_schmidt_orth(my, G) :: "
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
	my.set_no_of_columns(pseudo_columns);
	G.set_no_of_rows(rows);
	G.set_no_of_columns(columns);

	G_Addr = G.get_data_address();
	my_Addr = my.get_data_address();

	for (lidia_size_t i = 0; i < pseudo_rows; i++) {
		for (lidia_size_t j = 0; j < pseudo_rows; j++)
			if (trans_flag)
				LiDIA::swap(my_Addr[i][j], Tr_my[i][j]);
			else
				LiDIA::swap(my_Addr[i][j], Tr_my[j][i]);
	}
	for  (lidia_size_t i = 0; i < rows; i++) {
		for (lidia_size_t j = 0; j < columns; j++)
			if (trans_flag)
				LiDIA::swap(G_Addr[i][j], Tr_G[j][i]);
			else
				LiDIA::swap(G_Addr[i][j], Tr_G[i][j]);
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
double bigfloat_lattice_basis::lll_check_search()
{
	debug_handler("bigfloat_lattice_basis", "lll_check_search()");
	if (!(Tr_check_basis()))
		lidia_error_handler("bigfloat_lattice_basis", "lll_check_search() :: "
				    "lattice is no basis !");
	bigfloat_lattice_basis A(columns, rows);
	bigint_lattice_basis BiA(columns, rows);
	sdigit bit_len;
	sdigit a, b;

	Tr_trans_swap(A);
	BiA.bigintify(bit_len, A);
	BiA.Tr_lll_check_search(a, b);
	A.Tr_trans_swap(*this);
	return(rint((static_cast<double>(a)/static_cast<double>(b))*1000)/1000.0);
}



void bigfloat_lattice_basis::lll_check_search(sdigit& a, sdigit& b)
{
	debug_handler("bigfloat_lattice_basis", "lll_check_search(a, b)");
	if (!(Tr_check_basis()))
		lidia_error_handler("bigfloat_lattice_basis", "lll_check(a, b) :: "
				    "lattice is no basis !");
	bigfloat_lattice_basis A;
	bigint_lattice_basis BiA;
	sdigit bit_len;

	Tr_trans_swap(A);
	BiA.bigintify(bit_len, A);
	BiA.Tr_lll_check_search(a, b);
	A.Tr_trans_swap(*this);
}



bool bigfloat_lattice_basis::lll_check(double y)
{
	debug_handler("bigfloat_lattice_basis", "lll_check(y)");
	if (!(Tr_check_basis()))
		lidia_error_handler("bigfloat_lattice_basis", "lll_check(y) :: "
				    "lattice is no basis !");
	bigfloat_lattice_basis A(columns, rows);
	bigint_lattice_basis BiA(columns, rows);
	sdigit bit_len;
	bool flag;
	Tr_trans_swap(A);
	BiA.bigintify(bit_len, A);
	flag = BiA.Tr_lll_check(y);
	A.Tr_trans_swap(*this);
	return(flag);
}



bool bigfloat_lattice_basis::lll_check(sdigit a, sdigit b)
{
	debug_handler("bigfloat_lattice_basis", "lll_check(a, b)");
	if (!(Tr_check_basis()))
		lidia_error_handler("bigfloat_lattice_basis", "lll_check(a, b) :: "
				    "lattice is no basis !");
	bigfloat_lattice_basis A;
	bigint_lattice_basis BiA;
	sdigit bit_len;
	bool flag;
	Tr_trans_swap(A);
	BiA.bigintify(bit_len, A);
	flag = BiA.Tr_lll_check(static_cast<double>(a)/static_cast<double>(b));
	A.Tr_trans_swap(*this);
	return(flag);
}



//
// Conversion
//
void bigfloat_lattice_basis::extract_basis(const bigfloat_lattice_gensys& gsys, lidia_size_t& rank)
{
	debug_handler("bigfloat_lattice_basis", "extract_basis(gsys, rank)");
	bigfloat_lattice_basis TrThis;
	//
	// Transpose lattice gsys, then basis - checking
	//
	assign(gsys);
	Tr_trans_swap(TrThis);
	TrThis.assign_the_rest(*this);
	TrThis.Tr_extract_basis(rank);
	TrThis.Tr_trans_swap(*this);
}



bool bigfloat_lattice_basis::make_basis(const bigfloat_lattice_gensys& gsys, lidia_size_t& rank, sdigit x_factor)
{
	debug_handler("bigfloat_lattice_basis", "make_basis(gsys, rank)");
	bigfloat_lattice_basis TrThis;
	bool flag;
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
bool bigfloat_lattice_basis::Tr_make_basis(lidia_size_t& rank, sdigit x_factor)
{
	bigint_lattice_basis BiThis(rows, columns);
	bigint_matrix TrBiM(rows, columns);
	bigint_matrix TrTrBiM(columns, rows);
	bigint determ;
	bool flag;
	sdigit bit_len;

	BiThis.bigintify(bit_len, *this);
	if (bit_len == 1) {
		BiThis.trans_flag = trans_flag;
		flag = BiThis.Tr_make_basis(rank, x_factor);
		bigfloatify(*this, BiThis);
		return(flag);
	}
	else {
		TrBiM.assign(BiThis);
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
	}

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



void bigfloat_lattice_basis::Tr_extract_basis(lidia_size_t& rank)
{
	debug_handler("bigfloat_lattice_basis", "Tr_extract_basis(rank)");
	bool flag;

	rank = rows;
	for (lidia_size_t i = 0; i < rows; i++) {
		flag = false;
		for (lidia_size_t cx = 0; cx < columns; cx++)
			if (value[i][cx].is_approx_zero() == 0) {
				flag = true;
				break;
			}
		if (flag == false) {
			rank--;
			for (lidia_size_t j = i; j < rank; j++)
				bfl_swap_bfl(value[j], value[j+1]);
		}
	}
	set_no_of_rows(rank);
}



//
// Gram - Schmidt ported from bigfloat version
//
bool bigfloat_lattice_basis::Tr_gram_schmidt_orth(bigfloat** my, bigfloat** gso)
{
	debug_handler("bigfloat_lattice_basis", "Tr_gram_schmidt_orth(my, gso)");
	sdigit old_prec;
	bigfloat *tempvect0;
	bigfloat *tempvect1;

	old_prec = bigfloat::get_precision();
	bigfloat::set_precision(compute_precision());
	tempvect0 = new bigfloat[columns+rows];
	tempvect1 = &tempvect0[columns];
	vectsize = columns;
	for (lidia_size_t i = 0; i < rows; i++) {
		for (lidia_size_t cx = 0; cx < columns; cx++)
			gso[i][cx].assign(value[i][cx]);
		for (lidia_size_t j = 0; j < i; j++) {
			tempmz1.assign_zero();
			for (lidia_size_t cx = 0; cx < columns; cx++) {
				LiDIA::multiply(tempmz0, value[i][cx], gso[j][cx]);
				LiDIA::add(tempmz1, tempmz1, tempmz0);
			}
			if (tempvect1[j].is_approx_zero()) {
				delete[] tempvect0;
				bigfloat::set_precision(old_prec);
				return(false);
			}
			LiDIA::divide(my[j][i], tempmz1, tempvect1[j]);
		}
		for (lidia_size_t j = 0; j < i; j++) {
			bfl_scalmul_bfl(tempvect0, my[j][i], gso[j]);
			bfl_subtract_bfl(gso[i], gso[i], tempvect0);
		}
		bfl_scalprod_bfl(tempvect1[i], gso[i], gso[i]);
	}
	delete[] tempvect0;
	bigfloat::set_precision(old_prec);
	return(true);
}



bool bigfloat_lattice_basis::Tr_lll_check(double param)
{
	debug_handler("bigfloat_lattice_basis", "Tr_lll_check(param)");
	bigfloat_lattice_basis my(rows, rows);
	bigfloat_lattice_basis gso(rows, columns);
	bigfloat halb(0.5);


	if (Tr_gram_schmidt_orth(my.value, gso.value) == false)
		lidia_error_handler("bigfloat_lattice_basis", "Tr_lll_check(param) :: "
				    "lattice is no basis !");

	for (lidia_size_t i = 0; i < rows; i++)
		for (lidia_size_t j = 0; j < i; j++)
			if (my.value[j][i] > halb)
				return (false);

	for (lidia_size_t i = 1; i < rows; i++) {
		bfl_scalquad_bfl(tempmz2, value[i-1]);
		bfl_scalquad_bfl(tempmz3, value[i]);

		tempmz0.assign(param);
		LiDIA::multiply(tempmz0, tempmz0, tempmz2);
		LiDIA::square(tempmz1, my.value[i-1][i]);
		LiDIA::multiply(tempmz1, tempmz1, tempmz2);
		LiDIA::add(tempmz1, tempmz1, tempmz3);
		if (tempmz0.compare(tempmz1) > 0)
			return(false);
	}
	return (true);
}



void bigfloat_lattice_basis::Tr_lll_check_search(sdigit& a, sdigit& b)
{
	debug_handler("bigfloat_lattice_basis", "Tr_lll_check_search(a, b)");
	bool success;
	sdigit param_nom = 3;
	sdigit param_den = 4;
	sdigit upper_nom = 4;
	sdigit upper_den = 4;
	sdigit downer_nom = 2;
	sdigit downer_den = 4;
	bigfloat halb(0.5);
	bigfloat *vect;
	bigfloat_lattice_basis my(rows, rows);
	bigfloat_lattice_basis gso(rows, columns);


	vect = new bigfloat[rows];
	memory_handler(vect, "bigfloat_lattice_basis", "Tr_lll_check_search() :: "
                       "not enough memory !");

	if (Tr_gram_schmidt_orth(my.value, gso.value) == false)
		lidia_error_handler("bigfloat_lattice_basis", "Tr_lll_check_search() :: "
				    "lattice is no basis !");

	a = 0;
	b = 1;

	for (lidia_size_t i = 0; i < rows; i++) {
		bfl_scalquad_bfl(vect[i], value[i]);
		for (lidia_size_t j = 0; j < i; j++)
			if (my.value[j][i] > halb)
				return;
	}

	vectsize = columns;
	for (lidia_size_t j = 0; j < 20; j++)   // 2^-20
	{
		success = true;
		for (lidia_size_t i = 1; i < rows; i++) {
			LiDIA::divide(tempmz0, param_nom, param_den);
			LiDIA::multiply(tempmz0, tempmz0, vect[i-1]);
			LiDIA::square(tempmz1, my.value[i-1][i]);
			LiDIA::multiply(tempmz1, tempmz1, vect[i-1]);
			LiDIA::add(tempmz1, tempmz1, vect[i]);
			if (tempmz0.compare(tempmz1) > 0) {
				success = false;
				break;
			}

		}
		if (success == true) {
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
}



#ifdef LIDIA_NAMESPACE
}	// end of namespace LiDIA
#endif
