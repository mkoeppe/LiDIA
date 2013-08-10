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
#include	"LiDIA/lattices/bi_lattice_gensys.h"



#ifdef LIDIA_NAMESPACE
namespace LiDIA {
#endif



//
// basic funtions for bigint_lattice_gensys
//

//
// variables for temporary use
//

//
// for storing temporary values
//
bigint bigint_lattice_gensys::tempmz0;
bigint bigint_lattice_gensys::tempmz1;
bigint bigint_lattice_gensys::tempmz2;
bigint bigint_lattice_gensys::tempmz3;
bigint bigint_lattice_gensys::tempmz4;
bigint bigint_lattice_gensys::ergmz;

//
// temporary variables used in vector - operations
//
double bigint_lattice_gensys::vectdblz;
bigint bigint_lattice_gensys::vectbinz;
bigfloat bigint_lattice_gensys::vectbflz;

//
// vector - size for arithmetic operations
//
lidia_size_t bigint_lattice_gensys::vectsize;

//
// Constructor / Destructor
//

//
// Simple constructors
//
bigint_lattice_gensys::bigint_lattice_gensys():math_matrix< bigint > ()
{
	debug_handler("bigint_lattice_gensys", "bigint_lattice_gensys()");
	init_parameters();
	do_some_other_stuff();
}



// bigint_lattice_gensys::bigint_lattice_gensys(lidia_size_t n):math_matrix< bigint > (n,1)
// {
//   debug_handler("bigint_lattice_gensys", "bigint_lattice_gensys(n)");
//   init_parameters();
//   do_some_other_stuff();
// }

bigint_lattice_gensys::bigint_lattice_gensys(lidia_size_t m, lidia_size_t n):math_matrix< bigint > (m, n)
{
	debug_handler("bigint_lattice_gensys", "bigint_lattice_gensys(m, n)");
	init_parameters();
	do_some_other_stuff();
}



//
// Constructor with dimension (m, n)
// rows and columns and values bigint **abi
//
bigint_lattice_gensys::bigint_lattice_gensys(lidia_size_t m, lidia_size_t n, const bigint **abi):math_matrix< bigint > (m, n, abi)
{
	debug_handler("bigint_lattice_gensys", "bigint_lattice_gensys(m, n, abi)");
	init_parameters();
	do_some_other_stuff();
}



//
// Copy - constructor
// creating a copy of math_matrix< bigint > L
//
bigint_lattice_gensys::bigint_lattice_gensys(const math_matrix< bigint > & L):math_matrix< bigint > (L)
{
	debug_handler("bigint_lattice_gensys", "bigint_lattice_gensys(L)");
	init_parameters();
	do_some_other_stuff();
}



bigint_lattice_gensys::bigint_lattice_gensys(const bigint_lattice_gensys& L):math_matrix< bigint > (L)
{
	debug_handler("bigint_lattice_gensys", "bigint_lattice_gensys(L)");
	assign_the_rest(L);
}



//
// The destructor that frees allocated storage
//
bigint_lattice_gensys::~bigint_lattice_gensys()
{
	debug_handler("bigint_lattice_gensys", "~bigint_lattice_gensys()");
}



//
// Initialize Parameters
//
void bigint_lattice_gensys::init_parameters()
{
	debug_handler("bigint_lattice_gensys", "init_parameters()");
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
void bigint_lattice_gensys::assign(const bigint_lattice_gensys& L)
{
	debug_handler("bigint_lattice_gensys", "assign(L)");
	((math_matrix< bigint > *)this)->assign(L);
	assign_the_rest(L);
}



bigint_lattice_gensys& bigint_lattice_gensys::operator = (const bigint_lattice_gensys& L)
{
	debug_handler("bigint_lattice_gensys", "operator = (L)");
	assign(L);
	return(*this);
}



void assign(bigint_lattice_gensys& A, const bigint_lattice_gensys& B)
{
	debug_handler("bigint_lattice_gensys", "assign(A, B)");
	A.assign(B);
}



void bigint_lattice_gensys::assign_the_rest(const bigint_lattice_gensys& B)
{
	//
	// copy those elements which are not part of math_matrix< bigint >
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
void bigint_lattice_gensys::set_reduction_parameter(double par)
{
	debug_handler("bigint_lattice_gensys", "set_reduction_parameter(par)");
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



void bigint_lattice_gensys::set_reduction_parameter(sdigit nom, sdigit denom)
{
	debug_handler("bigint_lattice_gensys", "set_reduction_parameter(nom, denom)");
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



double bigint_lattice_gensys::get_reduction_parameter()
{
	debug_handler("bigint_lattice_gensys", "get_reduction_parameter()");
	return (y_par);
}



void bigint_lattice_gensys::get_reduction_parameter(sdigit& nom, sdigit& denom)
{
	debug_handler("bigint_lattice_gensys", "get_reduction_parameter(nom, denom)");
	nom = y_nom;
	denom = y_denom;
}



void bigint_lattice_gensys::set_representation_rows()
{
	debug_handler("bigint_lattice_gensys", "set_representation_rows()");
	trans_flag = false;
}



void bigint_lattice_gensys::set_representation_columns()
{
	debug_handler("bigint_lattice_gensys", "set_representation_columns()");
	trans_flag = true;
}



void bigint_lattice_gensys::randomize_vectors()
{
	debug_handler("bigint_lattice_gensys", "randomize_vectors()");
	bigint_lattice_gensys A;
	//
	// First step is to transpose the lattice
	// Second is to generate a permutation of the rows (using random - functions)
	// Third to transpose again
	//
	Tr_trans_swap(A);
	A.Tr_randomize_vectors();
	A.Tr_trans_swap(*this);
}



void bigint_lattice_gensys::sort_big_vectors()
{
	debug_handler("bigint_lattice_gensys", "sort_big_vectors()");
	bigint_lattice_gensys A;
	//
	// First step is to transpose the lattice
	// Second is to sort the rows by length (biggest first)
	// Third to transpose again
	//
	Tr_trans_swap(A);
	A.Tr_sort_vectors(1);
	A.Tr_trans_swap(*this);
}



void bigint_lattice_gensys::sort_small_vectors()
{
	debug_handler("bigint_lattice_gensys", "sort_small_vectors()");
	bigint_lattice_gensys A;
	//
	// First step is to transpose the lattice
	// Second is to sort the rows by length (smallest first)
	// Third to transpose again
	//
	Tr_trans_swap(A);
	A.Tr_sort_vectors(-1);
	A.Tr_trans_swap(*this);
}



void bigint_lattice_gensys::sort_vectors(bin_cmp cpf)
{
	debug_handler("bigint_lattice_gensys", "sort_vectors()");
	bigint_lattice_gensys A;
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
//
void bigint_lattice_gensys::lin_gen_system(math_matrix< bigint > & L, lidia_size_t& rank)
{
	debug_handler("bigint_lattice_gensys", "lin_gen_system(L, rank)");
	bigint_lattice_gensys TrThis;

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
	//
	// Give Timings and further informations back to *this
	//
	TrThis.Tr_trans_swap(*this);
	assign_the_rest(TrThis);
}



void bigint_lattice_gensys::lll_gensys(lidia_size_t& rank, sdigit x_factor)
{
	debug_handler("bigint_lattice_gensys", "lll_gensys(rank, x_factor)");
	bigint_lattice_gensys TrThis;

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



void bigint_lattice_gensys::lll_gensys(const math_matrix< bigint > & A, lidia_size_t& rank, sdigit x_factor)
{
	debug_handler("bigint_lattice_gensys", "lll_gensys(A, rank, x_factor)");
	((math_matrix< bigint > *)this)->assign(A);
	lll_gensys(rank, x_factor);
}



void bigint_lattice_gensys::lll_trans_gensys(math_matrix< bigint > & Tr, lidia_size_t& rank, sdigit x_factor)
{
	debug_handler("bigint_lattice_gensys", "lll_trans_gensys(Tr, rank, x_factor)");
	bigint_lattice_gensys TrThis;
	//
	// set up to right size
	//
	Tr.set_no_of_rows(columns);
	Tr.set_no_of_columns(columns);

	//
	// Transpose lattice and perform schnorr - euchner modified for
	// linear generating systems ( with transformation matrix `Tr` )
	// using double for approximation
	// return in rank the rank of the lattice
	//
	Tr_trans_swap(TrThis);
	TrThis.assign_the_rest(*this);
	if (x_factor < 2)
		TrThis.Tr_lll_trans_dbl_gensys(Tr, rank);
	else
		TrThis.Tr_lll_trans_bfl_gensys(Tr, rank, x_factor*8*SIZEOF_DOUBLE);
	TrThis.Tr_trans_swap(*this);
	assign_the_rest(TrThis);
}



//
// Modified lll for genertating systems
// of the form n x n+1
//
// Returns vector of dimension n
// where you can find the relations
//
void bigint_lattice_gensys::mlll(base_vector< bigint > & vect)
{
	debug_handler("bigint_lattice_gensys", "mlll(vect)");
	if (!(Tr_check_mlll()))
		lidia_error_handler("bigint_lattice_gensys", "mlll(vect) :: "
				    "Illegal size for modified lll (must have n x n+1)");

	bigint_lattice_gensys A;
	bigint *v;
	lidia_size_t dim;

	if (trans_flag)
		dim = columns;
	else
		dim = rows;
	v = new bigint[dim];
	memory_handler(v, "bigint_lattice_gensys", "mlll(vect) :: "
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
	base_vector< bigint > temp(v, dim);
	vect = temp;
}



void bigint_lattice_gensys::mlll(const math_matrix< bigint > & A, base_vector< bigint > & vect)
{
	debug_handler("bigint_lattice_gensys", "mlll(A, vect)");
	((math_matrix< bigint > *)this)->assign(A);
	mlll(vect);
}



void bigint_lattice_gensys::mlll(bigint *&v)
{
	debug_handler("bigint_lattice_gensys", "mlll(v)");
	if (!(Tr_check_mlll()))
		lidia_error_handler("bigint_lattice_gensys", "mlll(v) :: "
				    "Illegal size for modified lll (must have n x n+1)");

	bigint_lattice_gensys A;
	lidia_size_t dim;

	//
	// Generate Vector of bigints
	//
	if (trans_flag)
		dim = columns;
	else
		dim = rows;
	v = new bigint[dim];
	memory_handler(v, "bigint_lattice_gensys", "mlll(v) :: "
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



void bigint_lattice_gensys::mlll(const math_matrix< bigint > & A, bigint *&v)
{
	debug_handler("bigint_lattice_gensys", "mlll(A, v)");
	((math_matrix< bigint > *)this)->assign(A);
	mlll(v);
}



//
// Swaps two matrices by transposing them
//
//
// Reducing 2 matrix assignments
// 2*rows*columns real bigint (> 4 Bytes) assignments to
// 2*rows*columns*3 pointer assignments (4 Bytes)
//

void bigint_lattice_gensys::Tr_trans_swap(bigint_lattice_gensys& A)
{
	debug_handler("bigint_lattice_gensys", "Tr_trans_swap(A)");
	bigint **address;
	//
	// Swap matrix only
	//
	if (!trans_flag) {
		A.set_no_of_rows(rows);
		A.set_no_of_columns(columns);
		address = A.value;
		A.value = value;
		value = address;
		return;
	}
	//
	// Transpose and swap matrix
	//
	A.set_no_of_rows(columns);
	A.set_no_of_columns(rows);
	address = A.get_data_address();
	for (lidia_size_t i = 0; i < rows; i++)
		for (lidia_size_t j = 0; j < columns; j++)
			LiDIA::swap(value[i][j], address[j][i]);
}



void bigint_lattice_gensys::bigintify(sdigit& bit_len, const bigfloat_matrix & Bfl)
{
	debug_handler("bigint_lattice_gensys", "bigintify(bit_len, Bfl)");
	sdigit min_exp;
	bigint x_factor;
	bigfloat **BflAddr;
	set_no_of_rows(Bfl.get_no_of_rows());
	set_no_of_columns(Bfl.get_no_of_columns());
	BflAddr = Bfl.get_data_address();
	if (BflAddr[0][0].mantissa().is_zero())
		min_exp = 0;
	else
		min_exp = BflAddr[0][0].exponent();
	for (lidia_size_t i = 0; i < rows; i++)
		for (lidia_size_t j = 0; j < columns; j++)
			if (!(BflAddr[i][j].mantissa().is_zero()))
				min_exp = (BflAddr[i][j].exponent() < min_exp)?
					BflAddr[i][j].exponent():min_exp;
	if (min_exp < 0) {
		bit_len = -min_exp;
		for (lidia_size_t i = 0; i < rows; i++)
			for (lidia_size_t j = 0; j < columns; j++) {
				x_factor = BflAddr[i][j].mantissa();
				LiDIA::swap(value[i][j], x_factor);
				if (!(value[i][j].is_zero()))
					shift_left(value[i][j], value[i][j], (BflAddr[i][j].exponent()-min_exp));
			}
	}
	else {
		bit_len = 1;
		for (lidia_size_t i = 0; i < rows; i++)
			for (lidia_size_t j = 0; j < columns; j++) {
				x_factor = BflAddr[i][j].mantissa();
				LiDIA::swap(value[i][j], x_factor);
				if (!(value[i][j].is_zero()))
					shift_left(value[i][j], value[i][j], BflAddr[i][j].exponent());
			}
	}
}



//
// Real implementation of the algorithms
// They are working on the transposed lattice
//
void bigint_lattice_gensys::Tr_randomize_vectors()
{
	debug_handler("bigint_lattice_gensys", "Tr_randomize_vectors()");
	random_generator rg;
	char *bitvect;
	sdigit *perm;
	sdigit ran;
	bigint **temp;


	//
	// Allocate memory
	//
	bitvect = new char[rows];
	memory_handler(bitvect, "bigint_lattice_gensys", "Tr_randomize_vectors() :: "
		       "not enough memory !");
	perm = new sdigit[rows];
	memory_handler(perm, "bigint_lattice_gensys", "Tr_randomize_vectors() :: "
		       "not enough memory !");
	temp = new bigint*[rows];
	memory_handler(temp, "bigint_lattice_gensys", "Tr_randomize_vectors() :: "
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
			ran = static_cast<sdigit>(ran % rows);
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



void bigint_lattice_gensys::Tr_sort_vectors(sdigit vgl)
{
	debug_handler("bigint_lattice_gensys", "Tr_sort_vectors()");
	bigint *quads;

	//
	// Allocate memory for scalar product of the vectors
	//
	quads = new bigint[rows];
	memory_handler(quads, "bigint_lattice_gensys", "Tr_sort_vectors() :: "
		       "not enough memory !");

	//
	// Compute scalar products
	//
	vectsize = columns;
	for (lidia_size_t i = 0; i < rows; i++)
		bin_scalquad_bin(quads[i], value[i]);
	//
	// sort by scalar products ("length^2")
	// vgl = -1   ->   smallest vector first
	// vgl =  1   ->   biggest vector first
	//
	for (lidia_size_t i = 0; i < rows-1; i++)
		for (lidia_size_t j = i+1; j < rows; j++) {
			if (quads[i].compare(quads[j]) == vgl) {
				bin_swap_bin(value[i], value[j]);
				LiDIA::swap(quads[i], quads[j]);
			}
		}
	//
	// Free allocated storage
	//
	delete[] quads;
}



void bigint_lattice_gensys::Tr_sort_vectors(bin_cmp cpf)
{
	debug_handler("bigint_lattice_gensys", "Tr_sort_vectors()");

	//
	// Sort vectors by specified compare function cpf
	// Perform bubble sort (for typical sizes of the
	// lattice quicksort not nessesary)
	//
	vectsize = columns;
	for (lidia_size_t i = 0; i < rows-1; i++)
		for (lidia_size_t j = i+1; j < rows; j++)
			if (cpf(value[i], value[j], columns) > 0)
				bin_swap_bin(value[i], value[j]);
}



sdigit bigint_lattice_gensys::compute_read_precision()
{
	debug_handler("lattice_gensys", "compute_read_precision()");

	//
	// Give the digit size of longest value
	//
	sdigit new_prec;
	sdigit read_prec = 0;
	for (lidia_size_t x = 0; x < rows; ++x)
		for (lidia_size_t y = 0; y < columns; ++y)
			if (read_prec < (new_prec = static_cast<sdigit>(value[x][y].bit_length()/std::log(10.0)+1)))
				read_prec = new_prec;
	return (read_prec);
}



sdigit bigint_lattice_gensys::compute_precision()
{
	debug_handler("bigint_lattice_gensys", "compute_precision()");
	bigfloat alpha, zweipotq;
	sdigit read_prec = compute_read_precision();
	sdigit n2 = rows;
	sdigit prec;
	alpha_compute(alpha);
	zwei_pot_q_compute(zweipotq, n2, alpha);
	read_prec = 2*read_prec+rows-1;
	prec = static_cast<sdigit>(2*(((zweipotq.bit_length()+1)*
				       std::log(2.0)/std::log(10.0)+1)+read_prec)+(columns+rows)-1);
	return (prec);
}



void bigint_lattice_gensys::gamma_compute(bigfloat& g, sdigit l)
{
	debug_handler("bigint_lattice_gensys", "gamma_compute(g, l)");
	bigfloat ha[] = {1, 4.0/3.0, 2, 4, 8, 64.0/3.0, 64, 256};
	//
	// Computation of the l-th Hermite - Constant gamma,
	//
	bigfloat lg;
	if ((l > 0) && (l < 9))
		g.assign(ha[l-1]);
	else {
		lg.assign(0.75);
		LiDIA::log(lg, lg);
		g.assign(static_cast<sdigit>(l * (l-1)));
		g.divide_by_2();
		LiDIA::multiply(lg, lg, g);
		LiDIA::exp(g, lg);
	}
}



void bigint_lattice_gensys::alpha_compute(bigfloat& alpha)
{
	//
	// alpha = max{vect[1].l2_norm(),...,vect[columns].l2_norm()}
	// norm = vect[i].l2_norm(), 0 <= i < columns
	//
	debug_handler("bigint_lattice_gensys", "alpha_compute(alpha)");
	bigint bi_alpha;
	vectsize = columns;
	bin_l2_norm_bin(bi_alpha, value[0]);
	for (lidia_size_t i = 1; i < rows; ++i) {
		bin_l2_norm_bin(tempmz0, value[i]);
		if (bi_alpha.compare(tempmz0) < 0)
			bi_alpha.assign(tempmz0);
	}
	//
	// calculating sqrt of l2_norm
	//
	alpha.assign(bi_alpha);
	LiDIA::sqrt(alpha, alpha);
}



void bigint_lattice_gensys::zwei_pot_q_compute(bigfloat& zweipotq, sdigit& n2, bigfloat& alpha)
{
	debug_handler("bigint_lattice_gensys", "zwei_pot_q_compute(zwpq, n2, alpha)");
	sdigit beta = rows;
	bigfloat temp0;
	bigfloat temp1;
	bigfloat temp2;

	// Computation of beta = min (A.columns, A.rows)
	if (columns < rows)
		beta = columns;

	temp0.assign(n2);
	LiDIA::sqrt(temp0, temp0);
	temp0.divide_by_2();
	LiDIA::multiply(temp0, rows, temp0);
	temp1.assign(rows);
	LiDIA::sqrt(temp1, temp1);
	LiDIA::add(temp0, temp0, temp1);

	LiDIA::log(temp1, alpha);
	LiDIA::multiply(temp1, temp1, beta);
	LiDIA::exp(temp1, temp1);
	gamma_compute(temp2 , beta);

	LiDIA::sqrt(temp2, temp2);
	LiDIA::multiply(temp2, temp2, temp1);
	LiDIA::multiply(temp2, temp2, temp0);
	tempmz0.assign(n2);
	LiDIA::multiply(temp0, temp0, rows);
	LiDIA::sqrt(temp0, tempmz0);
	LiDIA::add(temp0, temp0, 2);
	LiDIA::multiply(zweipotq, temp0, temp2);

	temp0.assign(beta + rows);
	temp0.divide_by_2();
	LiDIA::subtract(temp0, temp0, 1);
	temp0.multiply_by_2();
	LiDIA::exp(temp0, temp0);
	LiDIA::multiply(temp0, zweipotq, temp0);

	temp2.assign(beta + 1);
	temp2.divide_by_2();
	temp1.assign(columns);
	LiDIA::log(temp1, temp1);
	LiDIA::multiply(temp2, temp2, temp1);
	LiDIA::exp(temp2, temp2);
	LiDIA::divide(temp0, temp0, temp2);
	LiDIA::ceil(zweipotq, temp0);
}



#ifdef LIDIA_NAMESPACE
}	// end of namespace LiDIA
#endif
