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
#include	"LiDIA/lattices/bf_lattice_gensys.h"



#ifdef LIDIA_NAMESPACE
namespace LiDIA {
#endif



//
// variables for temporary use
//

//
// For storing temporary values
//
bigfloat bigfloat_lattice_gensys::tempmz0;
bigfloat bigfloat_lattice_gensys::tempmz1;
bigfloat bigfloat_lattice_gensys::tempmz2;
bigfloat bigfloat_lattice_gensys::tempmz3;
bigfloat bigfloat_lattice_gensys::tempmz4;
bigfloat bigfloat_lattice_gensys::ergmz;

//
// For storing values during vector - operations
//
double bigfloat_lattice_gensys::vectdblz;
bigint bigfloat_lattice_gensys::vectbinz;
bigfloat bigfloat_lattice_gensys::vectbflz;

//
// vector - size for arithmetic operations
//
lidia_size_t bigfloat_lattice_gensys::vectsize;

//
// Constructor / Destructor
//

//
// Simple constructors
//

bigfloat_lattice_gensys::bigfloat_lattice_gensys():math_matrix< bigfloat > ()
{
	debug_handler("bigfloat_lattice_gensys", "bigfloat_lattice_gensys()");
	init_parameters();
	do_some_other_stuff();
}



// bigfloat_lattice_gensys::bigfloat_lattice_gensys(lidia_size_t n):math_matrix< bigfloat > (n,1)
// {
//   debug_handler("bigfloat_lattice_gensys", "bigfloat_lattice_gensys(n)");
//   init_parameters();
//   do_some_other_stuff();
// }

bigfloat_lattice_gensys::bigfloat_lattice_gensys(lidia_size_t m, lidia_size_t n):math_matrix< bigfloat > (m, n)
{
	debug_handler("bigfloat_lattice_gensys", "bigfloat_lattice_gensys(m, n)");
	init_parameters();
	do_some_other_stuff();
}



//
// Constructor with dimension (m, n)
// rows and columns and values bigfloat **abf
//
bigfloat_lattice_gensys::bigfloat_lattice_gensys(lidia_size_t m, lidia_size_t n, const bigfloat **abf):math_matrix< bigfloat > (m, n, abf)
{
	debug_handler("bigfloat_lattice_gensys", "bigfloat_lattice_gensys(m, n, abf)");
	init_parameters();
	do_some_other_stuff();
}



//
// Copy - constructor
// creating a copy of math_matrix< bigfloat > L
//
bigfloat_lattice_gensys::bigfloat_lattice_gensys(const math_matrix< bigfloat > & L):math_matrix< bigfloat > (L)
{
	debug_handler("bigfloat_lattice_gensys", "bigfloat_lattice_gensys(L)");
	init_parameters();
	do_some_other_stuff();
}



bigfloat_lattice_gensys::bigfloat_lattice_gensys(const bigfloat_lattice_gensys& L):math_matrix< bigfloat > (L)
{
	debug_handler("bigfloat_lattice_gensys", "bigfloat_lattice_gensys(L)");
	assign_the_rest(L);
	do_some_other_stuff();
}



//
// The destructor that frees allocated storage
//
bigfloat_lattice_gensys::~bigfloat_lattice_gensys()
{
	debug_handler("bigfloat_lattice_gensys", "~bigfloat_lattice_gensys()");
}



//
// Initialize Parameters
//
void bigfloat_lattice_gensys::init_parameters()
{
	debug_handler("bigfloat_lattice_gensys", "init_parameters()");
	//
	// set the defaults for reduction quality (same as pari)
	//
	y_par = 0.99;
	y_nom = 99;
	y_denom = 100;
	trans_flag = true;

	//
	// reset algorithm information
	//
	reduction_steps = 0;
	swapping_steps = 0;
	correction_steps = 0;
}



//
// Assignments
//
void bigfloat_lattice_gensys::assign(const bigfloat_lattice_gensys& L)
{
	debug_handler("bigfloat_lattice_gensys", "assign(L)");
	((math_matrix< bigfloat > *)this)->assign(L);
	assign_the_rest(L);
}



bigfloat_lattice_gensys& bigfloat_lattice_gensys::operator = (const bigfloat_lattice_gensys& L)
{
	debug_handler("bigfloat_lattice_gensys", "operator = (L)");
	assign(L);
	return(*this);
}



void assign(bigfloat_lattice_gensys& A, const bigfloat_lattice_gensys& B)
{
	debug_handler("bigfloat_lattice_gensys", "assign(A, B)");
	A.assign(B);
}



void bigfloat_lattice_gensys::assign_the_rest(const bigfloat_lattice_gensys& B)
{
	debug_handler("bigfloat_lattice_gensys", "assign_the_rest()");
	//
	// copy those elements which are not part of math_matrix< bigfloat >
	//
	y_par = B.y_par;
	y_nom = B.y_nom;
	y_denom = B.y_denom;
	trans_flag = B.trans_flag;

	//
	// copy algorithm information
	//
	reduction_steps = B.reduction_steps;
	swapping_steps = B.swapping_steps;
	correction_steps = B.correction_steps;
}



//
// Higher Functions
//
void bigfloat_lattice_gensys::set_reduction_parameter(double par)
{
	debug_handler("bigfloat_lattice_gensys", "set_reduction_parameter(par)");
	//
	// the parameter for reduction has to be in ]0.5 , 1.0]
	// it is not guarateed that the reduction terminates for 1.0
	// it has to be greater than 0.5 due to the schnorr - euchner modifications
	//
	if ((par > 0.5) && (par <= 1.0)) {
		y_par = par;
		y_nom = static_cast<sdigit>(y_par*100);
		y_denom = 100;
		return;
	}
	//
	// not in range -> default settings
	//
	y_par = 0.99;
	y_nom = 99;
	y_denom = 100;
}



void bigfloat_lattice_gensys::set_reduction_parameter(sdigit nom, sdigit denom)
{
	debug_handler("bigfloat_lattice_gensys", "set_reduction_parameter(nom, denom)");
	double par;
	//
	// the parameter nom/denom for reduction has to be in ]0.5 , 1.0]
	// it is not guarateed that the reduction terminates for 1.0
	// it has to be greater than 0.5 due to the schnorr - euchner modifications
	//
	par = static_cast<double>(nom)/static_cast<double>(denom);
	if ((par > 0.5) && (par <= 1.0)) {
		y_nom = nom;
		y_denom = denom;
		y_par = static_cast<double>(nom)/static_cast<double>(denom);
		return;
	}
	//
	// not in range -> default settings
	//
	y_par = 0.99;
	y_nom = 99;
	y_denom = 100;
}



double bigfloat_lattice_gensys::get_reduction_parameter()
{
	debug_handler("bigfloat_lattice_gensys", "get_reduction_parameter()");
	return (y_par);
}



void bigfloat_lattice_gensys::get_reduction_parameter(sdigit& nom, sdigit& denom)
{
	debug_handler("bigfloat_lattice_gensys", "get_reduction_parameter(nom, denom)");
	nom = y_nom;
	denom = y_denom;
}



void bigfloat_lattice_gensys::set_representation_rows()
{
	debug_handler("bigfloat_lattice_gensys", "set_representation_rows()");
	trans_flag = false;
}



void bigfloat_lattice_gensys::set_representation_columns()
{
	debug_handler("bigfloat_lattice_gensys", "set_representation_columns()");
	trans_flag = true;
}



void bigfloat_lattice_gensys::randomize_vectors()
{
	debug_handler("bigfloat_lattice_gensys", "randomize_vectors()");
	bigfloat_lattice_gensys A;
	//
	// First step is to transpose the lattice
	// Second is to generate a permutation of the rows (using rand48 - functions)
	// Third to transpose again
	//
	Tr_trans_swap(A);
	A.Tr_randomize_vectors();
	A.Tr_trans_swap(*this);
}



void bigfloat_lattice_gensys::sort_big_vectors()
{
	debug_handler("bigfloat_lattice_gensys", "sort_big_vectors()");
	bigfloat_lattice_gensys A;
	//
	// First step is to transpose the lattice
	// Second is to sort the rows by length (biggest first)
	// Third to transpose again
	//
	Tr_trans_swap(A);
	A.Tr_sort_vectors(1);
	A.Tr_trans_swap(*this);
}



void bigfloat_lattice_gensys::sort_small_vectors()
{
	debug_handler("bigfloat_lattice_gensys", "sort_small_vectors()");
	bigfloat_lattice_gensys A;
	//
	// First step is to transpose the lattice
	// Second is to sort the rows by length (smallest first)
	// Third to transpose again
	//
	Tr_trans_swap(A);
	A.Tr_sort_vectors(-1);
	A.Tr_trans_swap(*this);
}



void bigfloat_lattice_gensys::sort_vectors(bfl_cmp cpf)
{
	debug_handler("bigfloat_lattice_gensys", "sort_vectors()");
	bigfloat_lattice_gensys A;
	//
	// First step is to transpose the lattice
	// Second is to sort the rows by the given compare functions cpf
	// (See header file for more information)
	// Third to transpose again
	//
	Tr_trans_swap(A);
	A.Tr_sort_vectors(cpf);
	A.Tr_trans_swap(*this);
}



//
// Algorithm`s
//

//
// Buchmann - Kessler Algorthm for
// linaer generating systems
// result: transformation matrix using bigints
//
void bigfloat_lattice_gensys::lin_gen_system(math_matrix< bigint > & L, lidia_size_t& rank)
{
	debug_handler("bigfloat_lattice_gensys", "lin_gen_system(L, rank)");
	bigfloat_lattice_gensys TrThis;
	//
	// Give L the right dimension
	//
	L.set_no_of_rows(columns); //  Give lattice L right size
	L.set_no_of_columns(columns); //  columns x columns

	//
	// Transpose lattice and perform Buchmann - Kessler
	//
	Tr_trans_swap(TrThis);
	TrThis.assign_the_rest(*this);
	TrThis.Tr_lin_gen_system(L, rank);
	TrThis.Tr_trans_swap(*this);
}



//
// Buchmann - Kessler Algorthm for
// linaer generating systems
// result: transformation matrix using bigfloats
//
void bigfloat_lattice_gensys::lin_gen_system(math_matrix< bigfloat > & L, lidia_size_t& rank)
{
	debug_handler("bigfloat_lattice_gensys", "lin_gen_system(L, rank)");
	bigfloat_lattice_gensys TrThis;
	bigint_lattice_gensys TrBi;
	//
	// Transpose lattice and perform Buchmann - Kessler
	//
	Tr_trans_swap(TrThis);
	TrThis.assign_the_rest(*this);
	TrThis.Tr_lin_gen_system(TrBi, rank);
	TrThis.Tr_trans_swap(*this);
	//
	// Give L the right dimension
	//
	L.set_no_of_rows(TrBi.rows); //  Give lattice L right size
	L.set_no_of_columns(TrBi.columns); //  columns x columns

	//
	// convert to bigfloat
	//
	bigfloatify(L, TrBi);
}



//
// Schnorr - Euchner modified lll
//
void bigfloat_lattice_gensys::lll_gensys(lidia_size_t& rank, sdigit x_factor)
{
	debug_handler("bigfloat_lattice_gensys", "lll_gensys(rank, x_factor)");
	bigfloat_lattice_gensys TrThis;
	//
	// Transpose lattice and perform schnorr - euchner modified for
	// linear generating systems using double for approximation
	// return in rank the rank of the lattice
	//
	Tr_trans_swap(TrThis);
	TrThis.assign_the_rest(*this);
	if (x_factor < 2)
		TrThis.Tr_lll_dbl_gensys(rank);
	else
		TrThis.Tr_lll_bfl_gensys(rank, x_factor*8*SIZEOF_DOUBLE);
	TrThis.Tr_trans_swap(*this);
	assign_the_rest(TrThis);
}



void bigfloat_lattice_gensys::lll_gensys(const math_matrix< bigfloat > & A, lidia_size_t& rank,
					 sdigit x_factor)
{
	debug_handler("bigfloat_lattice_gensys", "lll_gensys(A, rank)");
	((math_matrix< bigfloat > *)this)->assign(A);
	lll_gensys(rank, x_factor);
}



//
// Schnorr - Euchner modified lll
//
void bigfloat_lattice_gensys::lll_trans_gensys(math_matrix< bigint > & T, lidia_size_t& rank, sdigit x_factor)
{
	debug_handler("bigfloat_lattice_gensys", "lll_trans_gensys(T, rank, x_factor)");
	bigfloat_lattice_gensys TrThis;

	//
	// Transpose lattice and perform schnorr - euchner modified for
	// linear generating systems using double for approximation
	// return in rank the rank of the lattice
	//
	Tr_trans_swap(TrThis);
	TrThis.assign_the_rest(*this);
	if (x_factor < 2)
		TrThis.Tr_lll_trans_dbl_gensys(T, rank);
	else
		TrThis.Tr_lll_trans_bfl_gensys(T, rank, x_factor*8*SIZEOF_DOUBLE);
	TrThis.Tr_trans_swap(*this);
	assign_the_rest(TrThis);
}



//
// Schnorr - Euchner modified lll
//
void bigfloat_lattice_gensys::lll_trans_gensys(math_matrix< bigfloat > & T, lidia_size_t& rank,
					       sdigit x_factor)
{
	debug_handler("bigfloat_lattice_gensys", "lll_trans_gensys(T, rank)");
	bigfloat_lattice_gensys TrThis;
	bigint_lattice_gensys TrBi;
	//
	// Transpose lattice and perform schnorr - euchner modified for
	// linear generating systems using double for approximation
	// return in rank the rank of the lattice
	//
	Tr_trans_swap(TrThis);
	TrThis.assign_the_rest(*this);
	if (x_factor < 2)
		TrThis.Tr_lll_trans_dbl_gensys(TrBi, rank);
	else
		TrThis.Tr_lll_trans_bfl_gensys(TrBi, rank, x_factor*8*SIZEOF_DOUBLE);
	TrThis.Tr_trans_swap(*this);
	assign_the_rest(TrThis);
	//
	// Convert to bigfloat
	//
	T.set_no_of_columns(TrBi.rows);
	T.set_no_of_rows(TrBi.columns);
	bigfloatify(T, TrBi);
}



void bigfloat_lattice_gensys::lll_var_gensys(lidia_size_t& rank, sdigit x_factor)
{
	debug_handler("bigfloat_lattice_gensys", "lll_var_gensys(rank, x_factor)");
	bigfloat_lattice_gensys A;
	bigint_lattice_gensys BiA;
	sdigit bit_len;

	Tr_trans_swap(A);
	A.assign_the_rest(*this);
	BiA.bigintify(bit_len, A);
	if (bit_len == 1) {
		BiA.y_par = A.y_par;
		if (x_factor < 2)
			BiA.Tr_lll_dbl_gensys(rank);
		else
			BiA.Tr_lll_bfl_gensys(rank, x_factor*8*SIZEOF_DOUBLE);
		A.reduction_steps = BiA.reduction_steps;
		A.swapping_steps = BiA.swapping_steps;
		A.correction_steps = BiA.correction_steps;
		bigfloatify(A, BiA);
	}
	else
		if (x_factor < 2)
			A.Tr_lll_var_dbl_gensys(BiA, bit_len, rank);
		else
			A.Tr_lll_var_bfl_gensys(BiA, bit_len, rank, x_factor*8*SIZEOF_DOUBLE);
	A.Tr_trans_swap(*this);
	assign_the_rest(A);
}



void bigfloat_lattice_gensys::lll_var_gensys(const math_matrix< bigfloat > & B, lidia_size_t& rank,
					     sdigit x_factor)
{
	debug_handler("bigfloat_lattice_gensys", "lll_var_gensys(B, rank, x_factor)");
	((math_matrix< bigfloat > *)this)->assign(B);
	lll_var_gensys(rank, x_factor);
}



void bigfloat_lattice_gensys::lll_trans_var_gensys(math_matrix< bigint > & T, lidia_size_t& rank,
						   sdigit x_factor)
{
	debug_handler("bigfloat_lattice_gensys", "lll_trans_var_gensys(T, rank, x_factor)");
	bigfloat_lattice_gensys A;
	bigint_lattice_gensys BiA;
	sdigit bit_len;

	Tr_trans_swap(A);
	A.assign_the_rest(*this);
	BiA.bigintify(bit_len, A);
	if (bit_len == 1) {
		BiA.y_par = A.y_par;
		if (x_factor < 2)
			BiA.Tr_lll_trans_dbl_gensys(T, rank);
		else
			BiA.Tr_lll_trans_bfl_gensys(T, rank, x_factor*8*SIZEOF_DOUBLE);
		bigfloatify(A, BiA);
		A.reduction_steps = BiA.reduction_steps;
		A.swapping_steps = BiA.swapping_steps;
		A.correction_steps = BiA.correction_steps;
	}
	else
		if (x_factor < 2)
			A.Tr_lll_trans_var_dbl_gensys(BiA, T, bit_len, rank);
		else
			A.Tr_lll_trans_var_bfl_gensys(BiA, T, bit_len, rank, x_factor*8*SIZEOF_DOUBLE);
	A.Tr_trans_swap(*this);
	assign_the_rest(A);
}



void bigfloat_lattice_gensys::lll_trans_var_gensys(math_matrix< bigfloat > & T, lidia_size_t& rank,
						   sdigit x_factor)
{
	debug_handler("bigfloat_lattice_gensys", "lll_trans_var_gensys(T, rank)");
	bigfloat_lattice_gensys A;
	bigint_lattice_gensys BiA;
	bigint_lattice_gensys TBi;
	sdigit bit_len;

	Tr_trans_swap(A);
	A.assign_the_rest(*this);
	BiA.bigintify(bit_len, A);
	if (bit_len == 1) {
		BiA.y_par = A.y_par;
		if (x_factor < 2)
			BiA.Tr_lll_trans_dbl_gensys(TBi, rank);
		else
			BiA.Tr_lll_trans_bfl_gensys(TBi, rank, x_factor*8*SIZEOF_DOUBLE);
		bigfloatify(A, BiA);
		A.reduction_steps = BiA.reduction_steps;
		A.swapping_steps = BiA.swapping_steps;
		A.correction_steps = BiA.correction_steps;
	}
	else
		if (x_factor < 2)
			A.Tr_lll_trans_var_dbl_gensys(BiA, TBi, bit_len, rank);
		else
			A.Tr_lll_trans_var_bfl_gensys(BiA, TBi, bit_len, rank, x_factor*8*SIZEOF_DOUBLE);
	A.Tr_trans_swap(*this);
	assign_the_rest(A);
	//
	// Conversion
	//
	T.set_no_of_rows(TBi.rows);
	T.set_no_of_columns(TBi.columns);
	bigfloatify(T, TBi);
}



//
// Modified lll for genertating systems
// of the form n x n+1
//
// Returns vector of dimension n
// where you can find the relations
//
//
// using bigints
//
void bigfloat_lattice_gensys::mlll(base_vector< bigint > & vect)
{
	debug_handler("bigfloat_lattice_gensys", "mlll(vect)");
	if (Tr_check_mlll() == false)
		lidia_error_handler("bigfloat_lattice_gensys", "mlll(vect) :: "
				    "Illegal size for modified lll (must have n x n+1)");

	bigfloat_lattice_gensys A;
	bigint *v;
	lidia_size_t dim;

	if (trans_flag)
		dim = columns;
	else
		dim = rows;
	v = new bigint[dim];
	memory_handler(v, "bigfloat_lattice_gensys", "mlll(vect) :: "
		       "not enough memory !");
	//
	// Transpose lattice and perform modified lll
	//
	Tr_trans_swap(A);
	A.assign_the_rest(*this);
	A.Tr_mlll_bfl(v);
	A.Tr_trans_swap(*this);
	//
	// return vector which satisfies the conditions
	// discribed in the manual
	// discribtion of the 0 - vector
	//
	assign_the_rest(A);
	base_vector< bigint > temp(v, dim);
	vect = temp;
}



void bigfloat_lattice_gensys::mlll(const math_matrix< bigfloat > & A, base_vector< bigint > & vect)
{
	debug_handler("bigfloat_lattice_gensys", "mlll(A, vect)");
	((math_matrix< bigfloat > *)this)->assign(A);
	mlll(vect);
}



void bigfloat_lattice_gensys::mlll(bigint *&v)
{
	debug_handler("bigfloat_lattice_gensys", "mlll(v)");
	if (!(Tr_check_mlll()))
		lidia_error_handler("bigfloat_lattice_gensys", "mlll(v) :: "
				    "Illegal size for modified lll (must have n x n+1)");

	bigfloat_lattice_gensys A;
	lidia_size_t dim;
	//
	// Generate Vector of bigint
	//

	if (trans_flag)
		dim = columns;
	else
		dim = rows;
	v = new bigint[dim];
	memory_handler(v, "bigfloat_lattice_gensys", "mlll(v) :: "
		       "not enough memory !");
	//
	// Transpose lattice and perform modified lll
	//
	Tr_trans_swap(A);
	A.assign_the_rest(*this);
	A.Tr_mlll_bfl(v);
	A.Tr_trans_swap(*this);
	assign_the_rest(A);
}



void bigfloat_lattice_gensys::mlll(const math_matrix< bigfloat > & A, bigint *&v)
{
	debug_handler("bigfloat_lattice_gensys", "mlll(A, v)");
	((math_matrix< bigfloat > *)this)->assign(A);
	mlll(v);
}



//
// using bigfloats
//
void bigfloat_lattice_gensys::mlll(base_vector< bigfloat > & vect)
{
	debug_handler("bigfloat_lattice_gensys", "mlll(vect)");
	if (!(Tr_check_mlll()))
		lidia_error_handler("bigfloat_lattice_gensys", "mlll(vect) :: "
				    "Illegal size for modified lll (must have n x n+1)");

	bigfloat_lattice_gensys A;
	bigint *v;
	bigfloat *bflv;
	lidia_size_t dim;

	if (trans_flag)
		dim = columns;
	else
		dim = rows;

	v = new bigint[dim];
	memory_handler(v, "bigfloat_lattice_gensys", "mlll(vect) :: "
		       "not enough memory !");
	bflv = new bigfloat[dim];
	memory_handler(bflv, "bigfloat_lattice_gensys", "mlll(vect) :: "
		       "not enough memory !");
	//
	// Transpose lattice and perform modified lll
	//
	Tr_trans_swap(A);
	A.assign_the_rest(*this);
	A.Tr_mlll_bfl(v);
	A.Tr_trans_swap(*this);
	assign_the_rest(A);
	//
	// return vector which satisfies the conditions
	// discribed in the manual
	// discribtion of the 0 - vector
	//
	for (lidia_size_t i = 0; i < dim; i++)
		bflv[i].assign(v[i]);
	base_vector< bigfloat > temp(bflv, dim);
	vect = temp;
}



void bigfloat_lattice_gensys::mlll(const math_matrix< bigfloat > & A, base_vector< bigfloat > & vect)
{
	debug_handler("bigfloat_lattice_gensys", "mlll(A, vect)");
	((math_matrix< bigfloat > *)this)->assign(A);
	mlll(vect);
}



void bigfloat_lattice_gensys::mlll(bigfloat *&v)
{
	debug_handler("bigfloat_lattice_gensys", "mlll(v)");
	if (!(Tr_check_mlll()))
		lidia_error_handler("bigfloat_lattice_gensys", "mlll(v) :: "
				    "Illegal size for modified lll (must have n x n+1)");

	bigfloat_lattice_gensys A;
	bigint *vbin;
	lidia_size_t dim;


	//
	// Generate Vector of bigint
	//
	if (trans_flag)
		dim = columns;
	else
		dim = rows;
	v = new bigfloat[dim];
	memory_handler(v, "bigfloat_lattice_gensys", "mlll(v) :: "
		       "not enough memory !");
	vbin = new bigint[dim];
	memory_handler(vbin, "bigfloat_lattice_gensys", "mlll(v) :: "
                       "not enough memory !");
	//
	// Transpose lattice and perform modified lll
	//
	Tr_trans_swap(A);
	A.assign_the_rest(*this);
	A.Tr_mlll_bfl(vbin);
	A.Tr_trans_swap(*this);
	assign_the_rest(A);
	for (lidia_size_t i = 0; i < dim; i++)
		v[i].assign(vbin[i]);
}



void bigfloat_lattice_gensys::mlll(const math_matrix< bigfloat > & A, bigfloat *&v)
{
	debug_handler("bigfloat_lattice_gensys", "mlll(A, v)");
	((math_matrix< bigfloat > *)this)->assign(A);
	mlll(v);
}



//
// Magie !!!
//

void bigfloat_lattice_gensys::Tr_trans_swap(bigfloat_lattice_gensys& Tr)
{
	debug_handler("bigfloat_lattice_gensys", "Tr_trans_swap(Tr)");
	bigfloat **TrBfl;
	//
	// Swap matrix only
	//
	if (trans_flag == false) {
		Tr.set_no_of_rows(rows);
		Tr.set_no_of_columns(columns);
		TrBfl = Tr.value;
		Tr.value = value;
		value = TrBfl;
		return;
	}
	//
	// Transpose and swap matrix
	//
	Tr.set_no_of_rows(columns);
	Tr.set_no_of_columns(rows);
	TrBfl = Tr.get_data_address();
	for (lidia_size_t i = 0; i < rows; i++)
		for (lidia_size_t j = 0; j < columns; j++)
			LiDIA::swap(value[i][j], Tr.value[j][i]);
}



//
// Real implementation of the algorithms
// They are working on the transposed lattice
//
void bigfloat_lattice_gensys::Tr_randomize_vectors()
{
	debug_handler("bigfloat_lattice_gensys", "Tr_randomize_vectors()");
	random_generator rg;
	char *bitvect;
	sdigit *perm;
	sdigit ran;
	bigfloat **temp;


	//
	// Allocate memory
	//
	bitvect = new char[rows];
	memory_handler(bitvect, "bigfloat_lattice_gensys", "Tr_randomize_vectors() :: "
		       "not enough memory !");
	perm = new sdigit[rows];
	memory_handler(perm, "bigfloat_lattice_gensys", "Tr_randomize_vectors() :: "
		       "not enough memory !");
	temp = new bigfloat*[rows];
	memory_handler(temp, "bigfloat_lattice_gensys", "Tr_randomize_vectors() :: "
		       "not enough memory !");

	//
	// Clear bitvector
	//
	for (lidia_size_t i = 0; i < rows; bitvect[i++] = 0);
	//
	// Generate a permutation
	// Try rows times to find valid value
	//
	for (lidia_size_t i = 0; i < rows; i++) {
		for (lidia_size_t j = 0; j < rows; j++) {
			rg >> ran;
			ran = static_cast<sdigit>(ran%rows);
			if (bitvect[ran] == 0) {
				bitvect[ran] = 1;
				perm[i] = ran;
				break;
			}
			else
				if (j == rows-1) {
					for (lidia_size_t l = 0; l < rows; l++)
						if (bitvect[l] == 0) {
							perm[i] = l;
							bitvect[l] = 1;
							break;
						}
					break;
				}
		}
	}

	//
	// Perform permutation on lattice
	//
	for (lidia_size_t i = 0; i < rows; i++)
		temp[perm[i]] = value[i];
	for (lidia_size_t i = 0; i < rows; i++)
		value[i] = temp[i];

	//
	// Free allocated storage
	//
	delete[] perm;
	delete[] bitvect;
	delete[] temp;
}



void bigfloat_lattice_gensys::Tr_sort_vectors(sdigit vgl)
{
	debug_handler("bigfloat_lattice_gensys", "Tr_sort_vectors()");
	bigfloat *quads;

	//
	// Allocate memory for scalar product of the vectors
	//
	quads = new bigfloat[rows];
	memory_handler(quads, "bigfloat_lattice_gensys", "Tr_sort_vectors() :: "
		       "not enough memory !");

	//
	// Compute scalar products
	//
	vectsize = columns;
	for (lidia_size_t i = 0; i < rows; i++)
		bfl_scalquad_bfl(quads[i], value[i]);
	//
	// sort by scalar products ("length^2")
	// vgl = -1   ->   smallest vector first
	// vgl =  1   ->   biggest vector first
	//
	for (lidia_size_t i = 0; i < rows-1; i++)
		for (lidia_size_t j = i+1; j < rows; j++) {
			if (quads[i].compare(quads[j]) == vgl) {
				bfl_swap_bfl(value[i], value[j]);
				LiDIA::swap(quads[i], quads[j]);
			}
		}
	//
	// Free allocated storage
	//
	delete[] quads;
}



void bigfloat_lattice_gensys::Tr_sort_vectors(bfl_cmp cpf)
{
	debug_handler("bigfloat_lattice_gensys", "Tr_sort_vectors()");
	//
	// Sort vectors by specified compare function cpf
	// Perform bubble sort
	//
	vectsize = columns;
	for (lidia_size_t i = 0; i < rows-1; i++)
		for (lidia_size_t j = i+1; j < rows; j++)
			if (cpf(value[i], value[j], columns) > 0)
				bfl_swap_bfl(value[i], value[j]);
}



//
// Algorithm - Subroutines
//
void bigfloat_lattice_gensys::bigfloatify(math_matrix< bigfloat > & Bf,
                                          const math_matrix< bigint > & Bi)
{
	debug_handler("bigfloat_lattice_gensys", "bigfloatify(Bf, Bi)");
	bigint **BiAddr;
	bigfloat **BfAddr;
	lidia_size_t Bf_rows = 0;
	lidia_size_t Bf_columns = 0;
	Bf.set_no_of_rows(Bf_rows = Bi.get_no_of_rows());
	Bf.set_no_of_columns(Bf_columns = Bi.get_no_of_columns());
	BiAddr = Bi.get_data_address();
	BfAddr = Bf.get_data_address();
	for (lidia_size_t i = 0; i < Bf_rows; i++)
		for (lidia_size_t j = 0; j < Bf_columns; j++)
			BfAddr[i][j].assign(BiAddr[i][j]);
}



sdigit bigfloat_lattice_gensys::compute_read_precision()
{
	debug_handler("bigfloat_lattice_gensys", "compute_read_precision()");
	sdigit new_prec;
	sdigit read_prec = 0;
	for (lidia_size_t x = 0; x < rows; ++x)
		for (lidia_size_t y = 0; y < columns; ++y)
			if (read_prec < (new_prec = static_cast<sdigit>(value[x][y].bit_length()/std::log(10.0)+1)))
				read_prec = new_prec;
	return (read_prec);
}



sdigit bigfloat_lattice_gensys::compute_precision()
{
	debug_handler("bigfloat_lattice_gensys", "compute_precision()");
	bigfloat alpha, zweipotq;
	sdigit read_prec = compute_read_precision();
	sdigit n2 = rows;
	sdigit prec;
	alpha_compute(alpha);
	zwei_pot_q_compute(zweipotq, n2, alpha);
	read_prec = 2*read_prec+rows-1;
	prec = static_cast<long>(2*(((zweipotq.exponent()+zweipotq.bit_length()+1)*
				     std::log(2.0)/std::log(10.0)+1)+read_prec)+(columns+rows)-1);
	return (prec);
}



void bigfloat_lattice_gensys::gamma_compute(bigfloat& g, long l)
{
	debug_handler("bigfloat_lattice_gensys", "gamma_compute(g, l)");
	bigfloat ha[] = {1, 4.0/3.0, 2, 4, 8, 64.0/3.0, 64, 256};
	//
	// Computation of the l-th Hermite - Constant gamma,
	//
	bigfloat lg;
	if ((l > 0) && (l < 9))
		g.assign(ha[l-1]);
	else {
		lg.assign(0.75);
		log(lg, lg);
		g.assign(static_cast<sdigit>(l * (l-1)));
		g.divide_by_2();
		LiDIA::multiply(lg, lg, g);
		exp(g, lg);
	}
}



void bigfloat_lattice_gensys::alpha_compute(bigfloat& alpha)
{
	// alpha = max{vect[1].l2_norm(),...,vect[columns].l2_norm()}
	// norm = vect[i].l2_norm(), 0 <= i < columns
	debug_handler("bigfloat_lattice_gensys", "alpha_compute(alpha)");
	vectsize = columns;
	bfl_l2_norm_bfl(alpha, value[0]);
	for (lidia_size_t i = 1; i < rows; ++i) {
		bfl_l2_norm_bfl(tempmz0, value[i]);
		if (alpha.compare(tempmz0) < 0)
			alpha.assign(tempmz0);
	}
}



void bigfloat_lattice_gensys::zwei_pot_q_compute(bigfloat& zweipotq, long& n2, bigfloat& alpha)
{
	debug_handler("bigfloat_lattice_gensys", "zwei_pot_q_compute(zwpq, n2, alpha)");
	long beta = rows;

	// Computation of beta = min (A.columns, A.rows)
	if (columns < rows)
		beta = columns;

	tempmz0.assign(n2);
	LiDIA::sqrt(tempmz0, tempmz0);
	tempmz0.divide_by_2();
	LiDIA::multiply(tempmz0, rows, tempmz0);
	tempmz1.assign(rows);
	LiDIA::sqrt(tempmz1, tempmz1);
	LiDIA::add(tempmz0, tempmz0, tempmz1);

	log(tempmz1, alpha);
	LiDIA::multiply(tempmz1, tempmz1, beta);
	exp(tempmz1, tempmz1);
	gamma_compute(tempmz2 , beta);

	LiDIA::sqrt(tempmz2, tempmz2);
	LiDIA::multiply(tempmz2, tempmz2, tempmz1);
	LiDIA::multiply(tempmz2, tempmz2, tempmz0);
	tempmz0.assign(n2);
	LiDIA::multiply(tempmz0, tempmz0, rows);
	LiDIA::sqrt(tempmz0, tempmz0);
	LiDIA::add(tempmz0, tempmz0, 2);
	LiDIA::multiply(zweipotq, tempmz0, tempmz2);

	tempmz0.assign(beta + rows);
	tempmz0.divide_by_2();
	LiDIA::subtract(tempmz0, tempmz0, 1);
	tempmz0.multiply_by_2();
	exp(tempmz0, tempmz0);
	LiDIA::multiply(tempmz0, zweipotq, tempmz0);

	tempmz2.assign(beta + 1);
	tempmz2.divide_by_2();
	tempmz1.assign(columns);
	log(tempmz1, tempmz1);
	LiDIA::multiply(tempmz2, tempmz2, tempmz1);
	exp(tempmz2, tempmz2);
	LiDIA::divide(tempmz0, tempmz0, tempmz2);
	ceil(zweipotq, tempmz0);
}



#ifdef LIDIA_NAMESPACE
}	// end of namespace LiDIA
#endif
