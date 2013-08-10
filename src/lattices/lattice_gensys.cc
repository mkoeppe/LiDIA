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
#include	"LiDIA/lattice_gensys.h"



#ifdef LIDIA_NAMESPACE
namespace LiDIA {
#endif



//
// Constructors
//
lattice_gensys::lattice_gensys(lidia_size_t n, lidia_size_t m)
{
	debug_handler("lattice_gensys", "lattice_gensys(n, m)");
	if ((n > 0) && (m > 0)) {
		comp_mode = bigfloat_mode;
		io_mode = own_io_mode;
		lbf = new bigfloat_lattice_basis(n, m);
		memory_handler(lbf, "lattice_gensys", "lattice_gensys(n, m) :: "
			       "not enough memory");
		return;
	}
	lidia_error_handler("lattice_gensys", "lattice_gensys(n, m) :: "
			    "illegal lattice size !");
}

//
// Initialising constructors
//
lattice_gensys::lattice_gensys(lidia_size_t n, lidia_size_t m, const double** ado)
{
	debug_handler("lattice_gensys", "lattice_gensys(n, m, ado)");
	if ((n > 0) && (m > 0)) {
		comp_mode = bigfloat_mode;
		io_mode = own_io_mode;
		lbf = new bigfloat_lattice_basis(n, m);
		memory_handler(lbf, "lattice_gensys", "lattice_gensys(n, m, ado) :: "
			       "not enough memory");
		sto(n, m, ado);
		return;
	}
	lidia_error_handler("lattice_gensys", "lattice_gensys(n, m, ado) :: "
			    "illegal lattice size !");
}

lattice_gensys::lattice_gensys(lidia_size_t n, lidia_size_t m, const bigint** abi)
{
	debug_handler("lattice_gensys", "lattice_gensys(n, m, abi)");
	if ((n > 0) && (m > 0)) {
		comp_mode = bigint_mode;
		io_mode = own_io_mode;
		lbi = new bigint_lattice_basis(n, m);
		memory_handler(lbi, "lattice_gensys", "lattice_gensys(n, m, abi) :: "
			       "not enough memory");
		sto(n, m, abi);
		return;
	}
	lidia_error_handler("lattice_gensys", "lattice_gensys(n, m, abi) :: "
			    "illegal lattice size !");
}

lattice_gensys::lattice_gensys(lidia_size_t n, lidia_size_t m, const bigfloat** abf)
{
	debug_handler("lattice_gensys", "lattice_gensys(n, m, abf)");
	if ((n > 0) && (m > 0)) {
		comp_mode = bigfloat_mode;
		io_mode = own_io_mode;
		lbf = new bigfloat_lattice_basis(n, m);
		memory_handler(lbf, "lattice_gensys", "lattice_gensys(n, m, abf) :: "
			       "not enough memory");
		sto(n, m, abf);
		return;
	}
	lidia_error_handler("lattice_gensys", "lattice_gensys(n, m, abf) :: "
			    "illegal lattice size !");
}

//
// Copy - constructor
//
lattice_gensys::lattice_gensys(const lattice_gensys& L)
{
	debug_handler("lattice_gensys", "lattice_gensys(L)");
	comp_mode = L.comp_mode;
	io_mode = L.io_mode;
	if (comp_mode == bigint_mode) {
		lbi = new bigint_lattice_basis(*L.lbi);
		memory_handler(lbi, "lattice_gensys", "lattice_gensys(L) :: "
			       "not enough memory");
	}
	else {
		lbf = new bigfloat_lattice_basis(*L.lbf);
		memory_handler(lbf, "lattice_gensys", "lattice_gensys(L) :: "
			       "not enough memory");
	}
}

//
// Special copy - constructors
//
lattice_gensys::lattice_gensys(const math_matrix< bigint > & A)
{
	debug_handler("lattice_gensys", "lattice_gensys(A)");
	lbi = new bigint_lattice_basis(A);
	memory_handler(lbi, "lattice_gensys", "lattice_gensys(A) :: "
		       "not enough memory");
	comp_mode = bigint_mode;
	io_mode = own_io_mode;
}

lattice_gensys::lattice_gensys(const math_matrix< bigfloat > & A)
{
	debug_handler("lattice_gensys", "lattice_gensys(A)");
	lbf = new bigfloat_lattice_basis(A);
	memory_handler(lbf, "lattice_gensys", "lattice_gensys(A) :: "
		       "not enough memory");
	comp_mode = bigfloat_mode;
	io_mode = own_io_mode;
}

//
// Destructor
//
lattice_gensys::~lattice_gensys()
{
	debug_handler("lattice_gensys", "~lattice_gensys()");
	if (comp_mode == bigint_mode)
		delete lbi;
	else
		delete lbf;
}

//
// Assignments
//
void lattice_gensys::assign(const lattice_gensys& L)
{
	debug_handler("lattice_gensys", "assign(L)");
	set_computing_mode(L.comp_mode);
	io_mode = L.io_mode;
	comp_mode = L.comp_mode;
	if (comp_mode == bigint_mode)
		lbi->assign(*L.lbi);
	else
		lbf->assign(*L.lbf);
}

void lattice_gensys::assign(const math_matrix< bigint > & M)
{
	debug_handler("lattice_gensys", "assign(M)");
	if (comp_mode == bigfloat_mode) {
		delete lbf;
		lbi = new bigint_lattice_basis(M);
		memory_handler(lbi, "lattice_gensys", "assign(M) :: "
			       "not enough memory !");
	}
	else
		lbi->assign(M);
	comp_mode = bigint_mode;
	io_mode = own_io_mode;
}

void lattice_gensys::assign(const math_matrix< bigfloat > & M)
{
	debug_handler("lattice_gensys", "assign(M)");
	if (comp_mode == bigint_mode) {
		delete lbi;
		lbf = new bigfloat_lattice_basis(M);
		memory_handler(lbf, "lattice_gensys", "assign(M) :: "
			       "not enough memory !");
	}
	else
		lbf->assign(M);
	comp_mode = bigfloat_mode;
	io_mode = own_io_mode;
}

lattice_gensys& lattice_gensys::operator = (const lattice_gensys& L)
{
	debug_handler("lattice_gensys", "operator = (L)");
	assign(L);
	return(*this);
}

lattice_gensys& lattice_gensys::operator = (const math_matrix< bigint > &M)
{
	debug_handler("lattice_gensys", "operator = (M)");
	assign(M);
	return(*this);
}

lattice_gensys& lattice_gensys::operator = (const math_matrix< bigfloat > &M)
{
	debug_handler("lattice_gensys", "operator = (M)");
	assign(M);
	return(*this);
}

void assign(lattice_gensys& A, const lattice_gensys& B)
{
	debug_handler("lattice_gensys", "assign(A, B)");
	A.assign(B);
}

void assign(lattice_gensys& A, const math_matrix< bigint > & Bbi)
{
	debug_handler("lattice_gensys", "assign(A, Bbi)");
	A.assign(Bbi);
}

void assign(lattice_gensys& A, const math_matrix< bigfloat > &Bbf)
{
	debug_handler("lattice_gensys", "assign(A, Bbf)");
	A.assign(Bbf);
}

//
// Input / Output
//
std::istream& operator >> (std::istream& in, lattice_gensys& L)
{
	debug_handler("lattice_gensys", "operator >> (in, L)");
	if (L.comp_mode == bigint_mode)
		in >> *(L.lbi);
	else
		in >> *(L.lbf);
	L.io_mode = own_io_mode;
	return (in);
}

std::ostream& operator << (std::ostream& out, const lattice_gensys& L)
{
	debug_handler("lattice_gensys", "operator << (out, L)");
	int mode = 0;
	switch (L.io_mode) {
	case pari_io_mode : {
		mode = GP_MODE;
	} break;
	case maple_io_mode : {
		mode = MAPLE_MODE;
	} break;
	case mathematica_io_mode : {
		mode = MATHEMATICA_MODE;
	} break;
	case own_io_mode : {
		mode = BEAUTY_MODE;
	} break;
	default : {
		lidia_error_handler("lattice_gensys", "operator << (out, L) :: "
				    "illegal IO Mode !");
	} break;
	}
	if (L.comp_mode == bigint_mode) {
		L.lbi->set_print_mode(mode);
		out << *(L.lbi);
	}
	else {
		L.lbf->set_print_mode(mode);
		out << *(L.lbf);
	}
	return (out);
}

//
// Element Operations
//
void lattice_gensys::member(lidia_size_t x, lidia_size_t y, double& d)
{
	debug_handler("lattice_gensys", "member(x, y, d)");
	bigfloat temp;
	if (comp_mode == bigint_mode)
		temp.assign(lbi->member(x, y));
	else
		temp.assign(lbf->member(x, y));
	temp.doublify(d);
}

void lattice_gensys::member(lidia_size_t x, lidia_size_t y, bigint& bi)
{
	debug_handler("lattice_gensys", "member(x, y, bi)");
	if (comp_mode == bigint_mode)
		bi.assign(lbi->member(x, y));
	else
		lbf->member(x, y).bigintify(bi);
}

void lattice_gensys::member(lidia_size_t x, lidia_size_t y, bigfloat& bf)
{
	debug_handler("lattice_gensys", "member(x, y, bf)");
	if (comp_mode == bigint_mode)
		bf.assign(lbi->member(x, y));
	else
		bf.assign(lbf->member(x, y));
}

void lattice_gensys::sto(lidia_size_t x, lidia_size_t y, const double& d)
{
	debug_handler("lattice_gensys", "sto(x, y, d)");
	if (comp_mode == bigint_mode) {
		bigint temp;
		temp.assign(d);
		lbi->sto(x, y, temp);
	}
	else {
		bigfloat temp(d);
		lbf->sto(x, y, temp);
	}
}

void lattice_gensys::sto(lidia_size_t x, lidia_size_t y, const bigint& bi)
{
	debug_handler("lattice_gensys", "sto(x, y, bi)");
	if (comp_mode == bigint_mode)
		lbi->sto(x, y, bi);
	else {
		bigfloat temp(bi);
		lbf->sto(x, y, temp);
	}
}

void lattice_gensys::sto(lidia_size_t x, lidia_size_t y, const bigfloat& bf)
{
	debug_handler("lattice_gensys", "sto(x, y, bf)");
	if (comp_mode == bigint_mode) {
		bigint temp;
		bf.bigintify(temp);
		lbi->sto(x, y, temp);
	}
	else
		lbf->sto(x, y, bf);
}

void lattice_gensys::member(lidia_size_t& n, lidia_size_t& m, double**& ado)
{
	debug_handler("lattice_gensys", "member(n, m, ado)");
	lidia_size_t rows = get_no_of_rows();
	lidia_size_t columns = get_no_of_columns();
	bigfloat temp;

	n = rows;
	m = columns;
	ado = new double*[rows];
	memory_handler(ado, "lattice_gensys", "member(n, m, ado) :: "
		       "not enough memory !");
	ado[0] = new double[rows*columns];
	memory_handler(ado[0], "lattice_gensys", "member(n, m, ado) :: "
		       "not enough memory !");
//
// Set pointers (tricky)
//
	for (lidia_size_t x = 0; x < rows; ado[x] = &ado[0][x++*columns]);

	if (comp_mode == bigint_mode) {
		for (lidia_size_t x = 0; x < rows; ++x)
			for (lidia_size_t y = 0; y < columns; ++y) {
				temp.assign(lbi->member(x, y));
				temp.doublify(ado[x][y]);
			}
	}
	else {
		for (lidia_size_t x = 0; x < rows; ++x)
			for (lidia_size_t y = 0; y < columns; ++y)
				lbf->member(x, y).doublify(ado[x][y]);
	}
}

void lattice_gensys::member(lidia_size_t& n, lidia_size_t& m, bigint**& abi)
{
	debug_handler("lattice_gensys", "member(n, m, abi)");
	lidia_size_t rows = get_no_of_rows();
	lidia_size_t columns = get_no_of_columns();
	bigfloat temp;

	n = rows;
	m = columns;
	abi = new bigint*[rows];
	memory_handler(abi, "lattice_gensys", "member(n, m, abi) :: "
		       "not enough memory !");
	abi[0] = new bigint[rows*columns];
	memory_handler(abi[0], "lattice_gensys", "member(n, m, abi) :: "
		       "not enough memory !");

	for (lidia_size_t x = 0; x < rows; abi[x] = &abi[0][x++*columns]);

	if (comp_mode == bigint_mode) {
		for (lidia_size_t x = 0; x < rows; ++x)
			for (lidia_size_t y = 0; y < columns; ++y)
				abi[x][y].assign(lbi->member(x, y));
	}
	else {
		for (lidia_size_t x = 0; x < rows; ++x)
			for (lidia_size_t y = 0; y < columns; ++y)
				lbf->member(x, y).bigintify(abi[x][y]);
	}
}

void lattice_gensys::member(lidia_size_t& n, lidia_size_t& m, bigfloat**& abf)
{
	debug_handler("lattice_gensys", "member(n, m, abf)");
	lidia_size_t rows = get_no_of_rows();
	lidia_size_t columns = get_no_of_columns();

	n = rows;
	m = columns;
	abf = new bigfloat*[rows];
	memory_handler(abf, "lattice_gensys", "member(n, m, abf) :: "
		       "not enough memory !");
	abf[0] = new bigfloat[rows*columns];
	memory_handler(abf[0], "lattice_gensys", "member(n, m, abf) :: "
		       "not enough memory !");
	for (lidia_size_t x = 0; x < rows; abf[x] = &abf[0][x++*columns]);

	if (comp_mode == bigint_mode) {
		for (lidia_size_t x = 0; x < rows; ++x)
			for (lidia_size_t y = 0; y < columns; ++y)
				abf[x][y].assign(lbi->member(x, y));
	}
	else {
		for (lidia_size_t x = 0; x < rows; ++x)
			for (lidia_size_t y = 0; y < columns; ++y)
				abf[x][y].assign(lbf->member(x, y));
	}
}

void lattice_gensys::member(math_matrix< bigint > & Mbi)
{
	debug_handler("lattice_gensys", "member(Mbi)");
	lidia_size_t rows = get_no_of_rows();
	lidia_size_t columns = get_no_of_columns();
	bigfloat temp;
	bigint temp_bi;

	Mbi.set_no_of_columns(columns);
	Mbi.set_no_of_rows(rows);

	if (comp_mode == bigint_mode) {
		for (lidia_size_t x = 0; x < rows; ++x)
			for (lidia_size_t y = 0; y < columns; ++y)
				Mbi.sto(x, y, lbi->member(x, y));
	}
	else {
		for (lidia_size_t x = 0; x < rows; ++x)
			for (lidia_size_t y = 0; y < columns; ++y) {
				lbf->member(x, y).bigintify(temp_bi);
				Mbi.sto(x, y, temp_bi);
			}
	}
}

void lattice_gensys::member(math_matrix< bigfloat > & Mbf)
{
	debug_handler("lattice_gensys", "member(Mbf)");
	lidia_size_t rows = get_no_of_rows();
	lidia_size_t columns = get_no_of_columns();
	bigfloat temp;

	Mbf.set_no_of_rows(rows);
	Mbf.set_no_of_columns(columns);

	if (comp_mode == bigint_mode) {
		for (lidia_size_t x = 0; x < rows; ++x)
			for (lidia_size_t y = 0; y < columns; ++y) {
				temp.assign(lbf->member(x, y));
				Mbf.sto(x, y, temp);
			}
	}
	else {
		for (lidia_size_t x = 0; x < rows; ++x)
			for (lidia_size_t y = 0; y < columns; ++y)
				Mbf.sto(x, y, lbf->member(x, y));
	}
}

void lattice_gensys::sto(lidia_size_t n, lidia_size_t m, const double** ado)
{
	debug_handler("lattice_gensys", "sto(n, m, ado)");
	lidia_size_t rows = get_no_of_rows();
	lidia_size_t columns = get_no_of_columns();
	bigfloat temp;
	bigint temp_bi;

	if ((n == rows) && (m == columns)) {
		if (comp_mode == bigint_mode) {
			for (lidia_size_t x = 0; x < rows; ++x)
				for (lidia_size_t y = 0; y < columns; ++y) {
					temp.assign(ado[x][y]);
					temp.bigintify(temp_bi);
					lbi->sto(x, y, temp_bi);
				}
		}
		else {
			for (lidia_size_t x = 0; x < rows; ++x)
				for (lidia_size_t y = 0; y < columns; ++y) {
					temp.assign(ado[x][y]);
					lbf->sto(x, y, temp);
				}
		}
		return;
	}
	lidia_error_handler("lattice_gensys", "sto(n, m, ado) :: "
			    "illegal array size !");
}

void lattice_gensys::sto(lidia_size_t n, lidia_size_t m, const bigint** abi)
{
	debug_handler("lattice_gensys", "sto(n, m, abi)");
	lidia_size_t rows = get_no_of_rows();
	lidia_size_t columns = get_no_of_columns();
	bigfloat temp;
	if ((n == rows) && (m == columns)) {
		if (comp_mode == bigint_mode) {
			for (lidia_size_t x = 0; x < rows; ++x)
				for (lidia_size_t y = 0; y < columns; ++y)
					lbi->sto(x, y, abi[x][y]);
		}
		else {
			for (lidia_size_t x = 0; x < rows; ++x)
				for (lidia_size_t y = 0; y < columns; ++y) {
					temp.assign(abi[x][y]);
					lbf->sto(x, y, temp);
				}
		}
		return;
	}
	lidia_error_handler("lattice_gensys", "sto(n, m, abi) :: "
			    "illegal array size !");
}

void lattice_gensys::sto(lidia_size_t n, lidia_size_t m, const bigfloat** abf)
{
	debug_handler("lattice_gensys", "sto(n, m, abf)");
	lidia_size_t rows = get_no_of_rows();
	lidia_size_t columns = get_no_of_columns();
	bigint temp;
	if ((n == rows) && (m == columns)) {
		if (comp_mode == bigint_mode) {
			for (lidia_size_t x = 0; x < rows; ++x)
				for (lidia_size_t y = 0; y < columns; ++y) {
					abf[x][y].bigintify(temp);
					lbi->sto(x, y, temp);
				}
		}
		else {
			for (lidia_size_t x = 0; x < rows; ++x)
				for (lidia_size_t y = 0; y < columns; ++y)
					lbf->sto(x, y, abf[x][y]);
		}
		return;
	}
	lidia_error_handler("lattice_gensys", "sto(n, m, abf) :: "
			    "illegal array size !");
}

lidia_size_t lattice_gensys::get_no_of_rows()
{
	debug_handler("lattice_gensys", "get_no_of_rows()");
	if (comp_mode == bigint_mode)
		return (lbi->get_no_of_rows());
	else
		return (lbf->get_no_of_rows());
}

lidia_size_t lattice_gensys::get_no_of_columns()
{
	debug_handler("lattice_gensys", "get_no_of_columns()");
	if (comp_mode == bigint_mode)
		return (lbi->get_no_of_columns());
	else
		return (lbf->get_no_of_columns());
}

//
// Type Checking
//

// Checks if you can convert to double without losing significant digits
bool lattice_gensys::check_double()
{
	debug_handler("lattice_gensys", "check_double()");
	lidia_size_t rows = get_no_of_rows();
	lidia_size_t columns = get_no_of_columns();
	if (comp_mode == bigfloat_mode) {
		for (lidia_size_t x = 0; x < rows; ++x) {
			for (lidia_size_t y = 0; y < columns; ++y)
				if (is_double(lbf->member(x, y)) == 0)
					return (false);
		}
	}
	else
		return (false);

	return (true);
}

// Checks if you can covert to bigint (no floating point)
bool lattice_gensys::check_bigint()
{
	debug_handler("lattice_gensys", "check_bigint()");
	lidia_size_t rows = get_no_of_rows();
	lidia_size_t columns = get_no_of_columns();
	bigfloat x_factor;
	bigfloat zero;
	if (comp_mode == bigfloat_mode) {
		for (lidia_size_t x = 0; x < rows; ++x) {
			for (lidia_size_t y = 0; y < columns; ++y) {
				x_factor = lbf->member(x, y);
				round(x_factor, zero);
				subtract(zero, zero, x_factor);
				if (zero.is_approx_zero() == 0)
					return(false);
			}
		}
	}
	return (true);
}

void lattice_gensys::set_computing_mode(sdigit m)
{
	debug_handler("lattice_gensys", "set_computing_mode(m)");
	lidia_size_t rows = get_no_of_rows();
	lidia_size_t columns = get_no_of_columns();

	if (m == double_mode)
		m = bigfloat_mode;
	if ((m != bigint_mode) && (m != bigfloat_mode))
		lidia_error_handler("lattice_gensys", "set_computing_mode(m) :: "
				    "illegal computing mode !");

	if (comp_mode == bigfloat_mode) {
		if (m == bigint_mode) {
			bigint temp;
			lbi = new bigint_lattice_basis(rows, columns);
			memory_handler(lbi, "lattice_gensys", "set_computing_mode(m) :: "
				       "not enough memory !");
			for (lidia_size_t x = 0; x < rows; ++x)
				for (lidia_size_t y = 0; y < columns; ++y) {
					lbf->member(x, y).bigintify(temp);
					lbi->sto(x, y, temp);
				}
			comp_mode = bigint_mode;
			lbi->set_reduction_parameter(lbf->get_reduction_parameter());
			delete lbf;
		}
	}
	else {
		if (m == bigfloat_mode) {
			bigfloat temp;
			lbf = new bigfloat_lattice_basis(rows, columns);
			memory_handler(lbf, "lattice_gensys", "set_computing_mode(m) :: "
				       "not enough memory !");
			for (lidia_size_t x = 0; x < rows; ++x)
				for (lidia_size_t y = 0; y < columns; ++y) {
					temp.assign(lbi->member(x, y));
					lbf->sto(x, y, temp);
				}
			comp_mode = bigfloat_mode;
			lbf->set_reduction_parameter(lbi->get_reduction_parameter());
			delete lbi;
		}
	}
}

sdigit lattice_gensys::get_computing_mode()
{
	debug_handler("lattice_gensys", "get_computing_mode()");
	return (comp_mode);
}

void lattice_gensys::set_io_mode(sdigit m)
{
	debug_handler("lattice_gensys", "set_io_mode(m)");
	if ((m == own_io_mode) ||
	    (m == pari_io_mode) ||
	    (m == maple_io_mode) ||
	    (m == mathematica_io_mode)) {
		io_mode = m;
		return;
	}
	lidia_error_handler("lattice_gensys", "set_io_mode(m) :: "
			    "illegal IO mode");
}

sdigit lattice_gensys::get_io_mode()
{
	debug_handler("lattice_gensys", "get_io_mode()");
	return (io_mode);
}


//
// Other Member Functions
//
void lattice_gensys::set_reduction_parameter(double y)
{
	debug_handler("lattice_gensys", "set_reduction_parameter(y)");
	if (comp_mode == bigfloat_mode)
		lbf->set_reduction_parameter(y);
	else
		lbi->set_reduction_parameter(y);
}

void lattice_gensys::set_reduction_parameter(sdigit a, sdigit b)
{
	debug_handler("lattice_gensys", "set_reduction_parameter(a, b)");
	if (comp_mode == bigfloat_mode)
		lbf->set_reduction_parameter(a, b);
	else
		lbi->set_reduction_parameter(a, b);
}

double lattice_gensys::get_reduction_parameter()
{
	debug_handler("lattice_gensys", "get_reduction_parameter()");
	if (comp_mode == bigfloat_mode)
		return(lbf->get_reduction_parameter());
	else
		return(lbi->get_reduction_parameter());
}

void lattice_gensys::get_reduction_parameter(sdigit& a, sdigit& b)
{
	debug_handler("lattice_gensys", "get_reduction_parameter(a, b)");
	if (comp_mode == bigfloat_mode)
		lbf->get_reduction_parameter(a, b);
	else
		lbi->get_reduction_parameter(a, b);
}

sdigit lattice_gensys::get_no_of_reduction_steps()
{
	debug_handler("lattice_gensys", "get_no_of_reduction_steps()");
	if (comp_mode == bigfloat_mode)
		return(lbf->get_no_of_reduction_steps());
	else
		return(lbi->get_no_of_reduction_steps());
}

sdigit lattice_gensys::get_no_of_swaps()
{
	debug_handler("lattice_gensys", "get_no_of_swaps()");
	if (comp_mode == bigfloat_mode)
		return(lbf->get_no_of_swaps());
	else
		return(lbi->get_no_of_swaps());
}

sdigit lattice_gensys::get_no_of_corrections()
{
	debug_handler("lattice_gensys", "get_no_of_corrections()");
	if (comp_mode == bigfloat_mode)
		return(lbf->get_no_of_corrections());
	else
		return(lbi->get_no_of_corrections());
}

//
// Time needed
//
timer& lattice_gensys::get_computing_time()
{
	debug_handler("lattice_gensys", "get_computing_time()");
	return(CT);
}

//
// Algorithms
//
void lattice_gensys::mlll(base_vector< bigint > & v)
{
	debug_handler("lattice_gensys", "mlll(v)");
	CT.start_timer();
	if (comp_mode == bigfloat_mode)
		lbf->mlll(v);
	else
		lbi->mlll(v);
	CT.stop_timer();
}

void lattice_gensys::mlll(const lattice_gensys& B, base_vector< bigint > & v)
{
	debug_handler("lattice_gensys", "lll_gensys(B, rank)");
	assign(B);
	CT.start_timer();
	if (comp_mode == bigfloat_mode)
		lbf->mlll(v);
	else
		lbi->mlll(v);
	CT.stop_timer();
}

void lattice_gensys::mlll(bigint*& v)
{
	debug_handler("lattice_gensys", "mlll(v)");
	CT.start_timer();
	if (comp_mode == bigfloat_mode)
		lbf->mlll(v);
	else
		lbi->mlll(v);
	CT.stop_timer();
}

void lattice_gensys::mlll(const lattice_gensys& B, bigint*& v)
{
	debug_handler("lattice_gensys", "lll_gensys(B, rank)");
	assign(B);
	CT.start_timer();
	if (comp_mode == bigfloat_mode)
		lbf->mlll(v);
	else
		lbi->mlll(v);
	CT.stop_timer();
}

void lattice_gensys::lll_gensys(lidia_size_t& rank)
{
	debug_handler("lattice_gensys", "lll_gensys(rank)");
	CT.start_timer();
	if (comp_mode == bigfloat_mode)
		lbf->lll_var_gensys(rank);
	else
		lbi->lll_gensys(rank);
	CT.stop_timer();
}

void lattice_gensys::lll_gensys(const lattice_gensys& B, lidia_size_t& rank)
{
	debug_handler("lattice_gensys", "lll_gensys(B, rank)");
	assign(B);
	CT.start_timer();
	if (comp_mode == bigfloat_mode)
		lbf->lll_var_gensys(rank);
	else
		lbi->lll_gensys(rank);
	CT.stop_timer();
}

void lattice_gensys::lll_trans_gensys(lattice_gensys& Tr, lidia_size_t& rank)
{
	debug_handler("lattice_gensys", "lll_trans_gensys(Tr, rank)");
	Tr.set_computing_mode(comp_mode);
	CT.start_timer();
	if (comp_mode == bigfloat_mode)
		lbf->lll_trans_var_gensys(*Tr.lbf, rank);
	else
		lbi->lll_trans_gensys(*Tr.lbi, rank);
	CT.stop_timer();
}


// The Buchmann - Kessler algoritm on the identic matrix, result : transformation matrix
void lattice_gensys::lin_gen_system(lattice_gensys& Tr, lidia_size_t& rank)
{
	debug_handler("lattice_gensys", "lin_gen_system(Tr, rank)");
	Tr.set_computing_mode(comp_mode);
	CT.start_timer();
	if (comp_mode == bigfloat_mode)
		lbf->lin_gen_system(*Tr.lbf, rank);
	else
		lbi->lin_gen_system(*Tr.lbi, rank);
	CT.stop_timer();
}

// Multiplies the matrix with the transformation matrix, result : reduced matrix
void lattice_gensys::compute_basis(const lattice_gensys& A, const lattice_gensys& B)
{
	debug_handler("lattice_gensys", "compute_basis(A, B)");
	if (A.comp_mode == B.comp_mode) {
		set_computing_mode(A.comp_mode);
		if (comp_mode == bigfloat_mode)
			lbf->compute_basis(*A.lbf, *B.lbf);
		else
			lbi->compute_basis(*A.lbi, *B.lbi);
		return;
	}
	lidia_error_handler("lattice_gensys", "compute_basis(A, B) :: "
			    "different computing modes !");
}

//
// Friend functions
//

void lll_gensys(math_matrix< bigint > & A,
                const math_matrix< bigint > & B,
                lidia_size_t& rank)
{
	debug_handler("lattice_gensys", "friend lll_gensys(A, B, rank)");
	bigint_lattice_basis temp(B);
	temp.lll_gensys(rank);
	A.assign(temp);
}

void lll_trans_gensys(math_matrix< bigint > & T,
                      math_matrix< bigint > & A,
                      lidia_size_t& rank)
{
	debug_handler("lattice_gensys", "friend lll_trans_gensys(T, A, rank)");
	bigint_lattice_basis temp(A);
	temp.lll_trans_gensys(T, rank);
	A.assign(temp);
}

void lll_gensys(math_matrix< bigfloat > & A,
                const math_matrix< bigfloat > & B,
                lidia_size_t& rank)
{
	debug_handler("lattice_gensys", "friend lll_gensys(A, B, rank)");
	bigfloat_lattice_basis temp(B);
	temp.lll_var_gensys(rank);
	A.assign(temp);
}

void lll_trans_gensys(math_matrix< bigint > & T,
                      math_matrix< bigfloat > & A,
                      lidia_size_t& rank)
{
	debug_handler("lattice_gensys", "friend lll_trans_gensys(T, A, rank)");
	bigfloat_lattice_basis temp(A);
	temp.lll_trans_var_gensys(T, rank);
	A.assign(temp);
}

void lll_trans_gensys(math_matrix< bigfloat > & T,
                      math_matrix< bigfloat > & A,
                      lidia_size_t& rank)
{
	debug_handler("lattice_gensys", "friend lll_trans_gensys(T, A, rank)");
	bigfloat_lattice_basis temp(A);
	temp.lll_trans_var_gensys(T, rank);
	A.assign(temp);
}

void mlll(math_matrix< bigint > & A,
          const math_matrix< bigint > & B,
          base_vector< bigint > & v)
{
	debug_handler("lattice_gensys", "friend mlll(A, B, v)");
	bigint_lattice_basis temp(B);
	temp.mlll(v);
	A.assign(temp);
}

void mlll(math_matrix< bigfloat > & A,
          const math_matrix< bigfloat > & B,
          base_vector< bigint > & v)
{
	debug_handler("lattice_gensys", "friend mlll(A, B, v)");
	bigfloat_lattice_basis temp(B);
	temp.mlll(v);
	A.assign(temp);
}

void  mlll(math_matrix< bigfloat > & A,
           const math_matrix< bigfloat > & B,
           base_vector< bigfloat > & v)
{
	debug_handler("lattice_gensys", "friend lll_gensys(A, B, rank)");
	bigfloat_lattice_basis temp(B);
	temp.mlll(v);
	A.assign(temp);
}

void mlll(math_matrix< bigint > & A,
          const math_matrix< bigint > & B,
          bigint*& v)
{
	debug_handler("lattice_gensys", "friend mlll(A, B, v)");
	bigint_lattice_basis temp(B);
	temp.mlll(v);
	A.assign(temp);
}

void mlll(math_matrix< bigfloat > & A,
          const math_matrix< bigfloat > & B,
          bigint*& v)
{
	debug_handler("lattice_gensys", "friend mlll(A, B, v)");
	bigfloat_lattice_basis temp(B);
	temp.mlll(v);
	A.assign(temp);
}

void mlll(math_matrix< bigfloat > & A,
          const math_matrix< bigfloat > & B,
          bigfloat*& v)
{
	debug_handler("lattice_gensys", "friend mlll(A, B, v)");
	bigfloat_lattice_basis temp(B);
	temp.mlll(v);
	A.assign(temp);
}

void lin_gen_system(math_matrix< bigint > & Tr,
                    const math_matrix< bigint > & A,
                    lidia_size_t& rank)
{
	debug_handler("lattice_gensys", "friend lin_gen_system(Tr, A, rank)");
	bigint_lattice_gensys TrBi(A);
	TrBi.lin_gen_system(Tr, rank);
}

void lin_gen_system(math_matrix< bigint > & Tr,
                    const math_matrix< bigfloat > & A,
                    lidia_size_t& rank)
{
	debug_handler("lattice_gensys", "friend lin_gen_system(Tr, A, rank)");
	bigfloat_lattice_gensys TrBf(A);
	TrBf.lin_gen_system(Tr, rank);
}

void lin_gen_system(math_matrix< bigfloat > & Tr,
                    const math_matrix< bigfloat > & A,
                    lidia_size_t& rank)
{
	debug_handler("lattice_gensys", "friend lin_gen_system(Tr, A, rank)");
	bigfloat_lattice_gensys TrBi(A);
	TrBi.lin_gen_system(Tr, rank);
}



#ifdef LIDIA_NAMESPACE
}	// end of namespace LiDIA
#endif
