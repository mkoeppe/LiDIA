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
//	Author	: Stefan Neis (SN)
//	Changes	: See CVS log
//
//==============================================================================================


#ifdef HAVE_CONFIG_H
# include	"config.h"
#endif
#include	"LiDIA/alg_number.h"
#include	"LiDIA/debug.h"
#include	"LiDIA/base_vector.h"
#include	"LiDIA/Fp_polynomial.h"
#include	"LiDIA/factorization.h"
#include	"LiDIA/rational_factorization.h"
#include	"LiDIA/bigfloat_lattice.h"
#include	<cctype>



#ifdef LIDIA_NAMESPACE
namespace LiDIA {
#endif



#ifdef LIDIA_DEBUG
int number_field::count = 0;
int order::count = 0;
#endif

nf_base * nf_base::dummy_base = new nf_base;
nf_base * nf_base::current_base = nf_base::dummy_base;

// private functions: Computing information, that is needed internally:
lidia_size_t nf_base::degree() const
{
	debug_handler_l("nf_base", "in member - function degree()", -2);
	if (!f.is_zero()) {
		return f.degree();
	}
	else {
		return table.get_no_of_columns();
	}
}



void nf_base::compute_conjugates()
{
	if (f.is_zero()) {
		if (table.get_no_of_columns() == 1) {
			lidia_error_handler("nf_base", "compute_conjugates()::Not initialized!");
			return;
		}
		else
			compute_base();
	}

	lidia_size_t re = 0, img = real_roots;
	conjugates.set_no_of_columns(degree());
	conjugates.set_no_of_rows(degree());

	register lidia_size_t i = 0;
	for (; i < real_roots; i++)
		conjugates.sto(i, 0, 1.0);
	for (; i < degree(); i++) {
		conjugates.sto(i, 0, 1.0);
		conjugates.sto(++i, 0, 0.0);
	}

	bigcomplex * wurz = new bigcomplex [degree()];
	do {
		int cnt = 0;
		cohen(f, wurz, 1, cnt); // the "1" means: polynomial is over R
		// the "cnt" means: store first root at pos. cnt
		for (i = 0; i < degree(); i++) {
			if (wurz[i].imag().is_zero()) {
				conjugates.sto(re++, 1, wurz[i].real());
			}
			else {
				conjugates.sto(img++, 1, wurz[i].real());
				conjugates.sto(img++, 1, wurz[i++].imag());
			}
		}
		if (re != real_roots) {
			bigfloat::set_precision(bigfloat::get_precision() + 1);
			warning_handler("field", "get_conjugate(...)::"
					"increasing precision!");
		}
	} while (re != real_roots);
	delete[] wurz;

	bigfloat tmp_bf, multiplier_bf;
	bigcomplex tmp_bc, multiplier_bc;

	for (i = 0; i < real_roots; i++) {
		tmp_bf.assign(conjugates(i, 1));
		multiplier_bf.assign(tmp_bf);
		for (register lidia_size_t j = 2; j < degree(); j++) {
			multiply(tmp_bf, tmp_bf, multiplier_bf);
			conjugates.sto(i, j, tmp_bf);
		}
	}

	for (; i < degree(); i += 2) {
		tmp_bc.assign(bigcomplex(conjugates(i, 1), conjugates(i+1, 1)));
		multiplier_bc.assign(tmp_bc);
		for (register lidia_size_t j = 2; j < degree(); j++) {
			multiply(tmp_bc, tmp_bc, multiplier_bc);
			conjugates.sto(i, j, tmp_bc.real());
			conjugates.sto(i+1, j, tmp_bc.imag());
		}
	}

	if (base_computed()) {
		math_matrix< bigfloat > TRAFO(degree(), degree());
		for (i = 0; i < degree(); i++)
			for (register lidia_size_t j = 0; j < degree(); j++)
				TRAFO.sto(i, j, bigfloat(base.member(i, j)));
		multiply(conjugates, conjugates, TRAFO);
		divide(conjugates, conjugates, bigfloat(den));
	}
}



void nf_base::compute_table() const
{
	debug_handler("nf_base", "in member - function compute_table()");
	debug_handler_c("", "", 0,
			std::cout << "Computing table for nf_base at" << this << std::endl);

	register lidia_size_t i, j, k = 0;
	polynomial< bigint > p1, p2;
	bigint * tmp;

	table.set_no_of_columns(degree());
	table.set_no_of_rows(static_cast<lidia_size_t>(static_cast<double>(degree())/2*(degree()+1)));

	if (!base_computed()) {
		polynomial< bigint > x; //polynom p(x) = x;
		x.assign_x();
		polynomial< bigint > p1(bigint(1));
		for (i = 0; i < degree(); i++, p1 *= x) {
			polynomial< bigint > p2(bigint(1));
			for (j = 0; j <= i; j++, p2 *= x) {
				debug_handler_c("nf_base", "compute_table", 0,
						std::cout << "computing x^" << j+i << ": " << (p1*p2) << f);
				polynomial< bigint > g = (p1*p2)%f;
				register lidia_size_t l = 0;
				for (; l <= g.degree(); l++)
					table.sto(k, l, g[l]);
				for (; l < degree(); l++)
					table.sto(k, l, 0);
				k++;
			}
		}
		debug_handler_c("nf_base", "compute_table():: Computed table\n", 0,
				std::cout << table << " for nf_base " << (*this) << std::flush);
	}
	else {
		bigint q, r;
		//temporary variable to hold the products of base elements:
		bigint_matrix init(degree(), degree()*(degree()+1)/2);
		for (i = 0; i < degree(); i++) {
			p1.set_data(tmp = base.get_column(i), degree());
			delete[] tmp;
			for (j = 0; j <= i; j++) {
				p2.set_data(tmp = base.get_column(j), degree());
				delete[] tmp;
				polynomial< bigint > g = (p1*p2)%f;
				init.sto_column(tmp = g.get_data(), g.degree()+1, (i*(i+1))/2+j);
				delete[] tmp;
			}
		}
		init.reginvimage(base, init);
		//Evtl. vorher testen, ob letzte Zeile 1 ist!
		init.set_no_of_rows(degree());

		debug_handler_c("nf_base", "in member - function compute_table(...)",
				3, std::cout << "den * MT is" << trans(init) << std::flush);

		for (i = 0; i < init.get_no_of_rows(); i++)
			for (j = 0; j < init.get_no_of_columns(); j++) {
				div_rem(q, r, init.member(i, j), den);
				if (!r.is_zero()) {
					lidia_error_handler_c("nf_base", "compute_table::internal error::"
							      "division unsuccesful!",
							      std::cout << "While dividing" << init.member(i, j);
							      std::cout << " by " << den << std::endl);
					return;
				}
				init.sto(i, j, q);
			}
		table = trans(init);
		debug_handler_c("nf_base", "compute_table():: Computed table\n", 0,
				std::cout << table << " for polynomial " << f;
				std::cout << "with transformation\n" << base << std::flush);
	}
}



void nf_base::compute_base() const
{
	debug_handler_l("nf_base", "in member_function `compute_base()'", 0);
	bigint * tmp = new bigint[degree()];
	bigint tmp2;

	alg_number rho(tmp, 1, this);
	polynomial< bigint > rho_pol(0);
	polynomial< bigint > temp_pol;

	long mu = 0;

	debug_handler_l("nf_base", "compute_base()::initialized", 0);
	while (rho_pol.degree() < degree()) {
		mu++;
		for (register lidia_size_t i = 0; i < degree() &&
			     rho_pol.degree() < degree(); i++) {
			tmp2.assign(tmp[i]);
			tmp[i].assign(mu);
			alg_number temp(tmp, 1, this);

			debug_handler_c("nf_base", "compute_base():", 0,
					std::cout << "\tTesting " << temp << std::endl;
					std::cout << " with matrix " << rep_matrix(temp) << std::flush);
			temp_pol = charpoly(temp);
			divide(temp_pol, temp_pol, gcd(temp_pol, derivative(temp_pol)));

			debug_handler_c("nf_base", "compute_base():", 0,
					std::cout << "\tMin.pol is " << temp_pol << std::endl << std::flush);

			if (temp_pol.degree() > rho_pol.degree()) {
				rho.assign(temp);
				rho_pol = temp_pol;
				tmp2.assign(mu);
			}

			tmp[i].assign(-mu);
			temp = alg_number(tmp, 1, this);

			debug_handler_c("nf_base", "compute_base():", 0,
					std::cout << "\tTesting " << temp << std::endl << std::flush);

			temp_pol = charpoly(temp);
			debug_handler_l("nf_base", "compute_base()::Computed char.pol!", 0);

			divide(temp_pol, temp_pol, gcd(temp_pol, derivative(temp_pol)));

			debug_handler_c("nf_base", "compute_base():", 0,
					std::cout << "\tMin.pol is " << temp_pol << std::endl << std::flush);

			if (temp_pol.degree() > rho_pol.degree()) {
				rho.assign(temp);
				rho_pol = temp_pol;
			}
			else tmp[i].assign(tmp2);
		}
	}
	delete[] tmp;
	f.assign(rho_pol);

	debug_handler_c("nf_base", "compute_base():", 7,
			std::cout << "Found polynomial " << f << " for ";
			std::cout << rho << std::endl << std::flush);

	real_roots = no_of_real_roots(f);

	// Update the nf_base:
	bigint_matrix T(degree(), degree());
	T.sto_column_vector(get_one(), degree(), 0);
	alg_number temp(rho);
	for (register lidia_size_t i = 1; i < degree()-1; i++) {
		T.sto_column_vector(temp.coeff_vector(), degree(), i);
		multiply(temp, temp, rho);
	}
	T.sto_column_vector(temp.coeff_vector(), degree(), degree()-1);

	// Now T expresses the powers of rho by the base elements of the nf_base.
	// To express the elements of the nf_base by powers of rho, we now
	// simply have to invert "T":

	debug_handler_c("nf_base", "compute_base():", 7,
			std::cout << "Compute base" << std::endl);

	base.resize(degree(), degree());
	base.adj(T);
	T.det(tmp2);
	den.assign(tmp2);

	for (register lidia_size_t k = 0; (k < degree()) && (!tmp2.is_one()); k++)
		for (register lidia_size_t j = 0; j < degree() && (!tmp2.is_one()); j++)
			tmp2.assign(gcd (base.member(k, j), tmp2));
	if (!tmp2.is_one()) {
		den /= tmp2;
		base /= tmp2;
	}
	debug_handler_c("nf_base", "compute_base():", 7,
			std::cout << "Trafo is " << base << " / " << den << std::endl);
}



const math_vector< bigint > & nf_base::get_one() const
{
	debug_handler("nf_base", "in member - function get_one()");
	register lidia_size_t i, j;

	if (One.size() > 0) {
		debug_handler_l("nf_base", "exit from get_one()", 4);
		return One;
	}
	One.set_capacity(degree());
	debug_handler("nf_base", "get_one()::capacity set, decide what's used.");

	// Linear equation system for 1: One* MT(nik)_{i,k)=e_n
	// instead of n every number less than n also works, but for n
	//		getting all the columns is easier/faster.
	// construct the linear equation system using MT

	if (using_necessary()) {
		// If necessary construct the multiplication table:
		if (!table_computed()) {
			compute_table();
		}
		debug_handler("get_one", "Constructing LES");
		bigint_matrix A(degree(), degree());
		base_vector< bigint > tmp2;
		for (i = 0, j = static_cast<lidia_size_t>(static_cast<double>(degree()-1)/2.0*degree());
		     i < degree(); i++, j++) {
			debug_handler("get_one", "Constructing LES");
			table.get_row_vector(tmp2, j);
			A.sto_column_vector(tmp2, degree(), i);
		}
		// solve LES
		debug_handler("get_one", "construct right side");
		bigint * tmp = new bigint[degree()];
		tmp[degree() - 1].assign_one();
		debug_handler_c("nf_base", "get_one(..)", 0,
				std::cout << "Solving " << A;
				std::cout << base_vector< bigint > (tmp, degree()));

		bigint_matrix solution = solve(A, tmp);
		if (solution.get_no_of_columns() != 2) {
			lidia_error_handler_c("nf_base", "get_one::nf_base does not "
					      "contain a unique representation of 1",
					      std::cout << "While constructing 1 in the following ";
					      std::cout << "nf_base:" << (*this) << "Solving " << A;
					      std::cout << base_vector< bigint > (tmp, degree());
					      std::cout << "Got solution:" << solution);
#ifdef LIDIA_NO_ERROR_ON_EXIT
			return math_vector< bigint > ();
#endif
		}
		delete[] tmp;
		solution.get_column_vector(One, 1);
	}
	else {
		// If we are in the simple case, 1 is just a polynomial of degree 0
		debug_handler("nf_base", "get_one():: just return [ 1 0 ... 0]");
		bigint * tmp = new bigint[degree()];
		tmp[0].assign_one();
		debug_handler("nf_base", "get_one():: set vector");
		One.set_data(tmp, degree());
		delete[] tmp;
	}
	debug_handler_l("nf_base", "other exit from get_one()", 4);
	return One;
}



const bigfloat & nf_base::get_conjugate(lidia_size_t i, lidia_size_t j)
	// get the j-th conjugate of w_i
{
	static long prec = 0;
	if (i< 0 || i >= degree()) {
		lidia_error_handler("nf_base",
				    "get_conjugate(lidia_size_t i, lidia_size_t j)"
				    "::i out of range (must be between 0 and n-1)");
#ifdef LIDIA_NO_ERROR_ON_EXIT
		return conjugates.member(0, 0);
#endif
	}

	if (j< 1 || j > degree()) {
		lidia_error_handler("nf_base",
				    "get_conjugate(lidia_size_t i, lidia_size_t j)"
				    "::j out of range (must be between 1 and n)");
#ifdef LIDIA_NO_ERROR_ON_EXIT
		return conjugates.member(0, 0);
#endif
	}

	if (conjugates.get_no_of_columns() <= 1 ||
	    prec < bigfloat::get_precision()) {
		compute_conjugates();
		prec = bigfloat::get_precision();
	}
	return conjugates.member(--j, i);
}



const bigfloat_matrix & nf_base::get_conjugates()
{
	static long prec = 0;

	if (conjugates.get_no_of_columns() <= 1 ||
	    prec < bigfloat::get_precision()) {
		compute_conjugates();
		prec = bigfloat::get_precision();
	}
	return conjugates;
}



const polynomial< bigint > & nf_base::which_polynomial() const
{
	if (f.is_zero()) {
		if (table.get_no_of_columns() == 1)
			lidia_error_handler("nf_base", "which_polynomial()::applied on "
					    "non-initialized number_field!");
		else
			compute_base();
	}
	return f;
}



// Copying:

void nf_base::assign(const nf_base &a)
{
	references = 0;
	f.assign(a.f);
	real_roots = a.real_roots;
	base.assign(a.base);
	den.assign(a.den);
	table.assign(a.table);
	One.assign(a.One);
	conjugates.assign(a.conjugates);
}



// Input/Output:

std::ostream& operator << (std::ostream &s, const nf_base &output)
{
	debug_handler("nf_base", "in operator << (std::ostream &, const nf_base &)");
	if (output.f.is_zero()) {
		output.compute_base();
	}
	if (!output.base_computed())
		s << output.f << std::flush;
	else {
		output.base.set_print_mode(LIDIA_MODE);
		s << '(' << output.f << " " << output.base << std::flush;
		if (!output.den.is_one())
			s << '/' << output.den;
		s << ')' << std::flush;
	}
	return s;
	//  s << output.table;
}



std::istream& operator >> (std::istream &s, nf_base &input)
{
	// Input works as follows: either input a polynomial, a matrix (MT or TRAFO),
	// or a pair (polynomial matrix). The last case is easily recognized, since
	// matrizes and vectors never start with '('.
	// Polynomials are only recognized, if they are in verbose form and
	// starting with a letter. This BUG (!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!)
	//                              ====================================
	// simplifies programming a lot.
	// Everything else is read as matrix. If the matrix has as many rows as
	// columns, it is a TRAFO with respect to current_base (especially if it has
	// just one row and column, it is assumed to be the identity matrix),
	// otherwise it is either an MT or a wrong input.

	debug_handler("nf_base", "in operator >> (std::istream &, nf_base &)");

	char c;
	s >> std::ws >> c;
	if (c == '(') {
		// Read Poly
		s >> input.f;
		input.real_roots = no_of_real_roots(input.f);
		// Read Trafo
		s >> input.base;
		if (input.base.get_no_of_columns() != 1 &&
		    (input.base.get_no_of_columns() != input.f.degree() ||
		     input.base.get_no_of_columns() != input.base.get_no_of_rows())) {
			lidia_error_handler ("nf_base", "operator >>::The given matrix "
					     "is not a trafo for the given polynomial!");
			return s;
		}
		s >> std::ws >> c;
		if (c == '/') {
			s >> input.den;
			s >> std::ws >> c;
		}
		else input.den = 1;
		if (c != ')') {
			lidia_error_handler("nf_base", "operator >>::Your input does "
					    "not look like an nf_base to me. Sorry!");
			return s;
		}
		input.table = bigint_matrix(1, 1);
		input.One = math_vector< bigint > ();
		input.conjugates = matrix< bigfloat > (1, 1);
	}
	else if ((c >= 'A' && c <= 'Z') || (c >= 'a' && c <= 'z')) {
		s.putback(c);
		// Read Poly;
		s >> input.f;
		input.real_roots = no_of_real_roots(input.f);
		input.den = 1;
		input.base = (input.table = bigint_matrix(1, 1));
		input.One = math_vector< bigint > ();
		input.conjugates = matrix< bigfloat > (1, 1);
	}
	else {
		s.putback(c);
		s >> input.table;
		if ((input.table.get_no_of_columns() == 1)
		    || (input.table.get_no_of_columns() == input.table.get_no_of_rows())) {
			if (nf_base::current_base == nf_base::dummy_base) {
				lidia_error_handler ("nf_base", "operator >>::You didn't give a mult. "
						     "table and I don't know with respect to what "
						     "your matrix is to be considered as "
						     "transformation matrix!");
				return s;
			}
			if (input.table.get_no_of_columns() != 1 &&
			    input.table.get_no_of_columns() != input.f.degree()) {
				lidia_error_handler ("nf_base", "operator >>::The given matrix "
						     "does not match with the current nf_base!");
				return s;
			}
			input.f.assign(nf_base::current_base->f);
			input.real_roots = nf_base::current_base->real_roots;
			if (nf_base::current_base->base.get_no_of_columns() == 1) {
				swap(input.base, input.table); // tricky & fast!!!
				input.table = bigint_matrix(1, 1);
				input.den.assign_one();
			}
			else {
				if (input.table.get_no_of_columns() == 1)
					input.base = nf_base::current_base->base;
				else
					input.base = nf_base::current_base->base * input.table;
				input.den.assign(nf_base::current_base->den);
				input.table = bigint_matrix(1, 1);
			}
			if (input.table.get_no_of_columns() == 1) {
				input.One.assign(nf_base::current_base->One);
				input.conjugates.assign(nf_base::current_base->conjugates);
			}
			else {
				input.One = math_vector< bigint > ();
				input.conjugates = matrix< bigfloat > (1, 1);
			}
			char c;
			do {
				s.get(c);
			} while (isspace(c));
			if (c == '/') {
				bigint deno;
				s >> deno;
				input.den *= deno;
			}
			else {
				if (c != '\n' && c != '\r')
					s.putback(c);
			}
		}
		else {
			if (static_cast<double>(input.table.get_no_of_columns()+1)/2.0 *
			    input.table.get_no_of_columns() != input.table.get_no_of_rows()) {
				lidia_error_handler ("nf_base", "operator >>::Your input does not "
						     "look like an nf_base to me. Sorry!");
				return s;
			}
			input.f.assign_zero();
			input.base = bigint_matrix(1, 1);
			input.den = 1;
			input.One = math_vector< bigint > ();
			input.conjugates = matrix< bigfloat > (1, 1);
		}
	}
	nf_base::current_base = & input;
	return(s);
}



// **************************************************************************

	// Constructors & destructor:
number_field::number_field(): base(nf_base::current_base)
{
	debug_handler_c("number_field", "in number_field"
			"()", 1, count++;
			std::cout << "\nNow we have " << count << " number_fields!\n");
	nf_base::current_base->inc_ref();
}



number_field::number_field(const polynomial< bigint > &p)
{
	debug_handler_c("number_field", "in number_field"
			"(const polynomial< bigint > &)", 1, count++;
			std::cout << "\nNow we have " << count << " number_fields!\n");
	base = new nf_base;
	base->f = p;
	base->real_roots = no_of_real_roots(p);
	base->references = 1;
	nf_base::current_base = base;
}



number_field::number_field(const bigint *v, lidia_size_t deg)
{
	debug_handler_c("number_field", "in number_field"
			"(const bigint *, lidia_size_t)", 1, count++;
			std::cout << "\nNow we have " << count << " number_fields!\n");
	polynomial< bigint > f(v, deg);
	base = new nf_base;
	base->f.assign(f);
	base->real_roots = no_of_real_roots(f);
	base->references = 1;
	nf_base::current_base = base;
}



number_field::number_field(const order & O): base(O.which_base())
{
	debug_handler_c("number_field", "in number_field(order *)", 1, count++;
			std::cout << "\nNow we have " << count << " number_fields!\n");
	base->inc_ref();
	if (base != nf_base::dummy_base)
		nf_base::current_base = base;
}



number_field::number_field(const number_field & F): base(F.base)
{
	debug_handler_c("number_field", "in number_field(const number_field &)", 1,
			count++; std::cout << "\nNow we have " << count << " number_fields!\n");
	base->inc_ref();
	if (base != nf_base::dummy_base)
		nf_base::current_base = base;
}



number_field::~number_field()
{
	debug_handler_c("number_field", "in ~number_field(const number_field &)", 1,
			count--;
			std::cout << "\nNow we have only " << count << " number_fields!\n");
	base->dec_ref();
}



// member functions:
const bigfloat & number_field::get_conjugate(lidia_size_t i, lidia_size_t j)
{
	// get the j-th conjugate of \alpha^i with precision q.
	return base->get_conjugate(i, j);
}



const bigfloat_matrix & number_field::get_conjugates()
{
	return base->get_conjugates();
}



number_field overfield(const number_field &A, const number_field &B)
{
	lidia_error_handler("field", "Construction of common overfields "
			    "is not yet supported.");
	return A;
}



void number_field::assign(const number_field & F)
{
	if (base != F.base) {
		base->dec_ref();
		base = F.base;
		base->inc_ref();
	}
}



// Comparision:
bool operator == (const number_field & A, const number_field & B)
{
	// Fields are considered to be equal, if generated by the same polynomial;
	if (A.base == B.base) return true;
	if (A.base->f.is_zero()) A.base->compute_base();
	if (B.base->f.is_zero()) B.base->compute_base();
	return (A.base->f == B.base->f && A.base->base == B.base->base
		&& A.base->den == B.base->den);
}



bool operator <= (const number_field & A, const number_field & B)
{
	if (A == B) return true;
	std::cout << "WARNING: number_field::operator <= \n";
	std::cout << "WARNING: Solving the subfield problem is difficult...\n";
	std::cout << "WARNING: So by now, we always return false!!\n";
	return false;
}



lidia_size_t number_field::no_of_real_embeddings() const
{
	if (base->f.is_zero()) {
		if (base->table.get_no_of_columns() == 1)
			lidia_error_handler("number_field", "no_of_real_embeddings()::applied on "
					    "non-initialized number_field!");
		else
			base->compute_base();
	}
	return base->real_roots;
}



// In-/Output:
std::ostream& operator << (std::ostream &s, const number_field  &K)
{
	debug_handler("number_field",
		      "in operator << (std::ostream &, const number_field &)");
	s << *(K.base) << std::flush;
	return s;
}



std::istream& operator >> (std::istream &s, number_field &K)
{
	K.base->dec_ref();
	K.base = new nf_base;
	K.base->inc_ref();
	s >> *(K.base);
	nf_base::current_base = K.base;
	return s;
}



// **************************************************************************

	// Constructors & destructor:
order::order(const nf_base * base1)
	: discriminant(),
	  base(const_cast<nf_base *>(base1))
{
	debug_handler_c("order", "in order()", 1,
			count++; std::cout << "\nNow we have " << count << " orders!\n");
	base->inc_ref();
}



order::order(const matrix< bigint > &in_table)
	: discriminant()
{
	debug_handler_c("order", "in order(const bigint_matrix &)", 1, count++;
			std::cout << "\nNow we have " << count << " orders!\n");
	base = new nf_base;
	base->table.assign(in_table);
	base->references = 1;
	nf_base::current_base = base;
}



order::order(const polynomial< bigint > & p,
	     const matrix< bigint > & coefficients,
	     const bigint & d)
	: discriminant()
{
	debug_handler_c("order", "in order (const polynomial< bigint > &, "
			" const bigint_matrix &, const bigint &)", 1, count++;
			std::cout << "\nNow we have " << count << " orders!\n");
	base = new nf_base;
	base->f.assign(p);
	base->real_roots = no_of_real_roots(p);
	base->base.assign(coefficients);
	base->den.assign(d);
	base->references = 1;
	nf_base::current_base = base;
}



order::order(const matrix< bigint > & coefficients,
	     const bigint & d, const nf_base * base1)
	: discriminant()
{
	debug_handler_c("order", "in order (const polynomial< bigint > &, "
			" const bigint_matrix &, const bigint &)", 1, count++;
			std::cout << "\nNow we have " << count << " orders!\n");
	base = new nf_base;
	base->f.assign(base1->f);
	base->real_roots = base1->real_roots;
	if (base1->base.get_no_of_columns() > 1) {
		LiDIA::multiply(base->base, base1->base, coefficients);
		multiply(base->den, base1->den, d);
	}
	else {
		base->base.assign(coefficients);
		base->den.assign(d);
	}
	base->references = 1;
	nf_base::current_base = base;
}



order::order(const order &O)
	: discriminant(O.discriminant),
	  base(O.base)
{
	debug_handler_c("order", "in order(const order &)", 1, count++;
			std::cout << "\nNow we have " << count << " orders!\n");
	base->inc_ref();
	nf_base::current_base = base;
}



order::order(const number_field & F)
	: discriminant(),
	  base(F)
{
	debug_handler_c("order", "in order(const number_field &)", 1, count++;
			std::cout << "\nNow we have " << count << " orders!\n");
	base->inc_ref();
	assign(maximize());
	nf_base::current_base = base;
}



order::~order() {
	debug_handler_c("order", "in ~order()", 1, count--;
			std::cout << "\nNow we have only " << count << " orders!\n");
	base->dec_ref();
}



// member functions:
const bigint & disc(const order & O)
{
	debug_handler("order", "in member - function disc()");
	if (!O.discriminant.is_zero()) return O.discriminant;

	bigint_matrix trace_matrix(O.degree(), O.degree());
	debug_handler("order", "in member - function disc() -- computing");
	bigint * tmp;
	bigrational Tr;
	// If necessary construct the multiplication table:
	if (!O.base->table_computed()) {
		O.base->compute_table();
	}
	for (register lidia_size_t i = 0; i < O.degree(); i++) {
		for (register lidia_size_t j = 0; j < i; j++) {
			tmp = O.base->table.get_row(static_cast<lidia_size_t>((static_cast<double>(i)/2.0)*(i+1))+j);
			Tr = trace(alg_number(tmp, 1, O.base));
			if (!Tr.denominator().is_one()) {
				debug_handler_c("order", "in member - function "
						"disc()", 0,
						std::cout << "Found an element with denominator ";
						std::cout << Tr.denominator() << std::endl);
				lidia_error_handler("order", "disc::internal error::"
						    "denominator is not 1");
#ifdef LIDIA_NO_EXIT_ON_ERROR
				return bigint();
#endif
			}
			trace_matrix.sto(i, j, Tr.numerator());
			trace_matrix.sto(j, i, Tr.numerator());
			delete[] tmp;
		}
		tmp = O.base->table.get_row(static_cast<lidia_size_t>(static_cast<double>(i)/2.0*(i+3)));
		Tr = trace(alg_number(tmp, 1, O.base));
		trace_matrix.sto(i, i, Tr.numerator());
		if (!Tr.denominator().is_one()) {
			debug_handler_c("order", "in member - function "
					"disc()", 0,
					std::cout << "Found an element with denominator ";
					std::cout << Tr.denominator() << std::endl);
			lidia_error_handler("order", "disc::internal error::"
					    "denominator is not 1");
#ifdef LIDIA_NO_ERROR_ON_EXIT
			return bigint();
#endif
		}
		delete[] tmp;
	}
	O.discriminant = trace_matrix.det();
	return O.discriminant;
}



const bigint & order::MT(lidia_size_t i, lidia_size_t j, lidia_size_t k)
{
	register lidia_size_t n = degree();
	if (!base->table_computed()) {
		if (base->f.is_zero()) {
			lidia_error_handler("order", "MT::applied on non-initialized order");
			return base->table.member(0, 0);
		}
		base->compute_table();
	}
	if ((i< 1) || (i > n)) {
		lidia_error_handler("order", "MT::first index out of range\n");
		return base->table.member(0, 0);
	}
	if ((j< 1) || (j > n)) {
		lidia_error_handler("order", "MT::second index out of range\n");
		return base->table.member(0, 0);
	}
	if ((k< 1) || (k > n)) {
		lidia_error_handler("order", "MT::third index out of range\n");
		return base->table.member(0, 0);
	}
	if (i < j) {
		n = i; i = j; j = n;
	}
	if (i%2)
		return base->table.member(i*((i-1) >> 1)+j-1, k-1);
	else
		return base->table.member((i >> 1)*(i-1)+j-1, k-1);
}



const bigfloat & order::get_conjugate(lidia_size_t i, lidia_size_t j)
{
	// get the j-th conjugate of w_i
	return base->get_conjugate(i, j);
}



const bigfloat_matrix & order::get_conjugates()
{
	return base->get_conjugates();
}



void order::assign(const order & O)
{
	if (base != O.base) {
		base->dec_ref();
		base = O.base;
		base->inc_ref();
	}
	discriminant.assign(O.discriminant);
}



// Cast to module
order::operator alg_ideal() const
{
	debug_handler("order", "in operator alg_ideal()");
	bigmod_matrix A(degree(), 1, 1);
	alg_ideal M(A, 1, base);
	debug_handler("order", "leaving operator alg_ideal()");
	return M;
}



// Comparision:
// We are interested in comparing orders only if they are over the same field!

// One function to compute the always needed transformation matrix:
void order::compare(const order & O, bool & this_in_O, bool & O_in_this) const
{
	debug_handler("order", "in member-function compare");
	if (base == O.base) {
		this_in_O = true;
		O_in_this = true;
		return;
	}
	if (base->f.is_zero() || !base->base_computed()) {
		if (base->table.get_no_of_columns() == 1) {
			lidia_error_handler("order", "compare(...)::comparing non-initialized"
					    "order");
			return;
		}
		nf_base * O1 = static_cast<nf_base *>(base);
		O1->compute_base();
	}
	if (O.base->f.is_zero() || !O.base->base_computed()) {
		if (O.base->table.get_no_of_columns() == 1) {
			lidia_error_handler("order", "compare(...)::comparing non-initialized"
					    "order");
			return;
		}
		nf_base * O1 = static_cast<nf_base *>(O.base);
		O1->compute_base();
	}
	if (O.base->f != base->f) {
		lidia_error_handler("order", "compare::Comparing orders with "
				    "different polynomials is not supported!");
		return;
	}
	bigint_matrix T;
	T.reginvimage(base->base, O.base->base);
	O_in_this = true;
	for (lidia_size_t i = 0; O_in_this && (i < degree()); i++)
		if (T.member(degree(), i) != bigint(1)) O_in_this = false;
	T.set_no_of_rows(degree());
	this_in_O = (abs(T.det()) == bigint(1));
}



bool operator == (const order & O1, const order & O2)
{
	debug_handler("order", "in operator == ");
	bool O1_in_O2;
	bool O2_in_O1;
	O1.compare(O2, O1_in_O2, O2_in_O1);
	return O1_in_O2 && O2_in_O1;
}



bool operator != (const order & O1, const order & O2)
{
	debug_handler("order", "in operator != ");
	bool O1_in_O2;
	bool O2_in_O1;
	O1.compare(O2, O1_in_O2, O2_in_O1);
	return !(O1_in_O2 && O2_in_O1);

}



bool operator <= (const order & O1, const order & O2)
{
	debug_handler("order", "in operator <= ");
	bool O1_in_O2;
	bool O2_in_O1;
	O1.compare(O2, O1_in_O2, O2_in_O1);
	return O1_in_O2;
}



bool operator < (const order & O1, const order & O2)
{
	debug_handler("order", "in operator < ");
	bool O1_in_O2;
	bool O2_in_O1;
	O1.compare(O2, O1_in_O2, O2_in_O1);
	return O1_in_O2 && !O2_in_O1;
}



bool operator >= (const order & O1, const order & O2)
{
	debug_handler("order", "in operator >= ");
	bool O1_in_O2;
	bool O2_in_O1;
	O1.compare(O2, O1_in_O2, O2_in_O1);
	return O2_in_O1;
}



bool operator > (const order & O1, const order & O2)
{
	debug_handler("order", "in operator >");
	bool O1_in_O2;
	bool O2_in_O1;
	O1.compare(O2, O1_in_O2, O2_in_O1);
	return O2_in_O1 && !O1_in_O2;
}



// friend functions:
void swap(order & O, order & Q) {
	swap (O.discriminant, Q.discriminant);

	nf_base * help = O.base;
	O.base = Q.base;
	Q.base = help;
}



// Number-theoretic functions:
// Dedekind test:
bool order::dedekind(const bigint & p, polynomial< bigint > & h2) const
	// returns true, if p might be an index divisor
	// and false otherwise.
	// More precisely: true is returned if either p _is_ an index divisor or
	//                 if no polynomial is given!
{
	debug_handler("order", "in member - function dedekind(const bigint &, polynomial< bigint > &)");

	if (base->f.is_zero()) return true;
	//Attention: h not initialized

	register lidia_size_t j;
	Fp_polynomial f(base->f, p);
	factorization< Fp_polynomial > fa;

	// Compute the factors fa[i] of f with fa[i]^{e_i}| f, e_1<e_2<....
	// We assume that fa will not contain any prime component !!!
	debug_handler_c("order", "dedekind(...)", 4,
			std::cout << "call squarefree_factor " << std::flush);
	square_free_decomp(fa, f);
	debug_handler_c("order", "dedekind(...)", 4,
			std::cout << "Factorisation of f = " << f << " is:" << fa << std::endl);

	polynomial< bigint > h(bigint(1)), h1(bigint(1)), f_factor_lifted;

	for (j = fa.no_of_composite_components() - 1; j >= 0; j--) {
		lidia_size_t i, n;

		const Fp_polynomial & f_factor = fa.composite_base(j).base();
		n = f_factor.degree();
		f_factor_lifted.set_degree(n);
		for (i = 0; i <= n; i++)
			f_factor_lifted[i] = f_factor[i];
		if (fa.composite_exponent(j) >= 2)
			h1 *= f_factor_lifted;
		debug_handler_c("order", "dedekind(...)", 4,
				std::cout << "(" << f_factor_lifted << ") ^ "
				<< fa.composite_exponent(j) << " = ");

		power(f_factor_lifted, f_factor_lifted, fa.composite_exponent(j));
		debug_handler_c("order", "dedekind(...)", 4,
				std::cout << f_factor_lifted << std::endl);
		h *= f_factor_lifted;
		debug_handler_c("order", "dedekind(...)", 4,
				std::cout << "intermediate h is " << h << std::endl << std::flush);
	}
	debug_handler_c("order", "dedekind(...)", 4,
			std::cout << "h is " << h << std::endl << std::flush);
	polynomial< bigint > f1 = base->f-h;

	debug_handler_c("order", "dedekind(...)", 4,
			std::cout << "f1 is " << f1 << std::endl << std::flush);
	f1 /= bigint(p); // h(t) is f1
	debug_handler_c("order", "dedekind(...)", 4,
			std::cout << "call gcd for " << h1 << " and ";
			std::cout << f1 << " mod " << p << std::endl << std::flush);

	Fp_polynomial h2_mod, h1_mod (h1, p), f1_mod(f1, p);
	gcd(h2_mod, h1_mod, f1_mod);

	register lidia_size_t i = h2_mod.degree();
	h2.set_degree(i);
	for (; i >= 0; i--)
		h2[i] = h2_mod[i];

	if (h2.degree()) {
		debug_handler_c("order", "dedekind(...)", 5,
				std::cout << p << " is an index divisor." << std::endl << std::flush);
		div_rem(h2, h1, base->f, h2);
		return true;
	}

	debug_handler_c("order", "dedekind(...)", 5,
			std::cout << p << " is not an index divisor." << std::endl << std::flush);
	return false;
}



// Computation of the pseudo radical

module order::pseudo_radical(const bigint & p) const
	// computes Ip + pO, the pseudo-radical of O/pO
	// p is assumed to be either a (composite) number with gcd(p, degree)= 1 or
	// a prime number.
	// Alle Rechnungen sollten intern modulo p geschehen!!?!?!
{
	debug_handler("order", "in member - function "
		      "pseudo_radical(const bigint &, bigint &)");
	register lidia_size_t j, k, n = degree();
	alg_number e(base), im_e(base);
	//initialised with zero!
	bigint q = p;
	bigmod_matrix A(n, n, p);
	bigint *tmp = new bigint[n];

	if ((p < 1000) && is_prime(p, 10)) {
		debug_handler_l("order", "pseudo_radical: prime case", 0);
		j = 1;
		// compute j such that the j-th power of p exceeds the dimension
		while (q < n) {
			j++;
			q = q*p;
		}
		for (k = 0; k < n; tmp[k++] = 0) {
			tmp[k] = 1; //tmp = e_k
			e = alg_number(tmp, 1, base);
			// e=e_k
			debug_handler_c("order", "pseudo_radical", 4,
					std::cout << "compute" << e << "^" << p << std::flush);
			power_mod_p(im_e, e, q, p); // Should be done mod p
			debug_handler_c("order", "pseudo_radical", 4,
					std::cout << "It's" << im_e << std::flush);
			A.sto_column_vector(im_e.coeff_vector(), degree(), k);
		}
	}
	else {
		// we know gcd(p,n!)=1 (otherwise we would have factored p)
		// Compute A as matrix of the map:
		//      (O/pO -> (O/pO->Z/pZ)
		//	( a   -> (b   ->Tr(ab)))
		// For this, compute (b_j -> Tr(b_j*b_k)) for each pair of
		// base_vectors (b_j,b_k) of O/pO and store this as A_{j,k}.
		debug_handler_l("order", "pseudo_radical: non-prime case", 0);
		bigrational Trace;
		for (j = 0; j < n; tmp[n-1] = 0, j++) {
			tmp[j] = 1; //tmp = e_j;
			e = alg_number(tmp, 1, base);
			// e = e_j;
			tmp[j] = 0;
			for (k = 0; k < n; tmp[k++] = 0) {
				tmp[k] = 1; //tmp = e_k;
				im_e = alg_number(tmp, 1, base);
				// im_e = e_k
				if ((Trace = trace(e*im_e)).denominator() != 1)
					A.sto(j, k, Trace.numerator());
			}
		}
	}
	//gemeinsamer Code:
	debug_handler_l("order", "pseudo_radical::common part of code", 4);
	A = kernel(A);
	debug_handler_l("order", "pseudo_radical::kernel computed part of code", 4);
	bigint_matrix COPY(n, k = A.get_no_of_columns());
	delete[] tmp;
	base_vector< bigint > tmp2;
	for (j = 0; j < n; j++) {
		A.get_row_vector(tmp2, j);
		COPY.sto_row_vector(tmp2, k, j);
	}
	debug_handler_l("order", "pseudo_radical: construct Ip", 0);
	module Ip(module(COPY, 1, base) +
		  module(alg_number(bigint(p), base)));
	debug_handler_c("order", "in member - function pseudo_radical(...)", 4,
			std::cout << "Radical is " << Ip << std::flush);
	return Ip;
}



// Locally maximize the order:

order order::maximize(const bigint & p) const
{
	debug_handler("order", "in member - function "
		      "maximize(const bigint &)");
	debug_handler_l("order", "in member - function maximize(...)::"
			"call pseudo_radical", 4);
	module Ip = pseudo_radical(p);
	debug_handler_l("order", "in member - function maximize(...)::"
			"call ring_of_multipliers", 4);
	order multipliers = Ip.ring_of_multipliers(p);
	order old_multipliers = *this;
	debug_handler("order", "in member - function "
		      "maximize(const bigint &, bigint &)::start loop");

	while ((multipliers.base->f.is_zero())?
	       (multipliers.base->table != old_multipliers.base->table):
	       (multipliers.base->base != old_multipliers.base->base ||
		multipliers.base->den != old_multipliers.base->den)) {
		// if true superset:
		debug_handler("order", "in member - function "
			      "maximize(const bigint &, bigint &)::in loop");
		debug_handler_l("order", "in member - function maximize(...)::"
				"call pseudo_radical", 4);
		Ip = multipliers.pseudo_radical(p);
		swap(old_multipliers, multipliers);
		debug_handler_l("order", "in member - function maximize(...)::"
				"call ring_of_multipliers", 4);
		multipliers.assign(Ip.ring_of_multipliers(p));
	}
	debug_handler("order", "in member - function "
		      "maximize(const bigint &, bigint &)::after loop");
	return multipliers;
}



// Maximizing the order:

order order::maximize() const
{
	debug_handler("order", "in member - function maximize()");
	rational_factorization a = disc(*this); // initialize a with discriminant
	order maximized(*this); // initialize maximal order
	bigint factor;
	polynomial< bigint > h;
	bigint* tmp;

	a.factor(); // find all factors (and smaller)
	debug_handler_c("order", "in member - function maximize()", 4,
			std::cout << "Factorization of Diskriminant is " << a << std::endl);
	for (register lidia_size_t i = 0; i < a.no_of_comp(); i++)
		if (a.exponent(i) > 1) {
			if (a.is_prime_factor(i)) {
				debug_handler_c("order", "maximize()", 3,
						std::cout << "calling dedekind for " << a.base(i);
						std::cout << "\n on\n" << maximized << std::flush);
				if (!base->f.is_zero()) {
					if (!maximized.dedekind(a.base(i), h))
						continue;
					// i.e. if the dedekind test could't maximize anything
					// start the next loop iteration, otherwise do a first
					//maximizing step:
					bigint_matrix A(degree(), 1);
					bigint deno(1);

					A.sto_column(tmp = h.get_data(), h.degree()+1, 0);
					delete[] tmp;
					if (maximized.base->base_computed()) {
						A.reginvimage(maximized.base->base, A);
						deno = A.member(degree(), 0);
						A.set_no_of_rows(degree());
						A *= base->den;
					}
					alg_number x(tmp = A.get_column(0), deno, maximized.base);
					delete[] tmp;
					debug_handler_c("order", "maximize()", 3,
							std::cout << "Polynomial is " << h;
							std::cout << "x is " << x;
							std::cout << "p is " << alg_number(bigint(a.base(i)),
											   maximized.base);
							std::cout << "Radical is ";
							std::cout << module(x, alg_number(a.base(i), maximized.base))
							<< std::flush);
					maximized.assign(module(x, alg_number(a.base(i), maximized.base)).ring_of_multipliers(a.base(i)));
					debug_handler_c("order", "maximize()", 3,
							std::cout << "maximized is " << maximized << std::flush);
					// if the discriminant will shrink sufficiently,
					// we can start the next loop iteration:
					if (a.exponent(i) < 2*(degree()-h.degree()+1))
						continue;
				}
			}
			// and use the Buchmann/Lenstra - algorithm for further
			// maximization, if necessary
			debug_handler_c("order", "in member - function maximize()::"
					"next  call to maximize(p)", 4,
					std::cout << "with p = " << a.base(i) << " exponent is "
					<< a.exponent(i));
			maximized.assign(maximized.maximize(a.base(i)));
			// i.e. the order is maximized locally at a.base(i)
			debug_handler("order", "in member - function maximize()::"
				      "returned from maximize(p)");
		}
	return maximized;
}



order order::maximize2() const
{
	debug_handler("order", "in member - function maximize()");
	rational_factorization a = disc(*this);
	// initialize a with discriminant
	order maximized(*this); // initialize maximal order
	bigint product(1);
	polynomial< bigint > h;
	bigint* tmp;

	a.factor(); // find factors
	debug_handler_c("order", "in member - function maximize()", 4,
			std::cout << "Factorization of Diskriminant is " << a << std::endl);
	for (register lidia_size_t i = 0; i < a.no_of_comp(); i++)
		if (a.exponent(i) > 1) {
			if (a.is_prime_factor(i)) {
				debug_handler_c("order", "maximize()", 3,
						std::cout << "calling dedekind for " << a.base(i);
						std::cout << "\n on\n" << maximized << std::flush);
				if (!base->f.is_zero()) {
					if (!maximized.dedekind(a.base(i), h))
						continue;
					// i.e. if the dedekind test could't maximize anything
					// start the next loop iteration, otherwise do a first
					//maximizing step:
					bigint_matrix A(degree(), 1);
					bigint deno(1);

					A.sto_column(tmp = h.get_data(), h.degree()+1, 0);
					delete[] tmp;
					if (maximized.base->base_computed()) {
						A.reginvimage(maximized.base->base, A);
						deno = A.member(degree(), 0);
						A.set_no_of_rows(degree());
						A *= base->den;
					}
					alg_number x(tmp = A.get_column(0), deno, maximized.base);
					delete[] tmp;
					debug_handler_c("order", "maximize()", 3,
							std::cout << "Polynomial is " << h << std::endl;
							std::cout << "x is " << x << std::endl;
							std::cout << "p is " << alg_number(bigint(a.base(i)),
											   maximized.base) << std::endl;
							std::cout << "Radical is ";
							std::cout << module(x, alg_number(a.base(i),
											  maximized.base)) << std::endl);
					maximized.assign(module(x, alg_number(a.base(i), maximized.base)).ring_of_multipliers(a.base(i)));
					debug_handler_c("order", "maximize()", 3,
							std::cout << "maximized is " << maximized << std::flush);
					// if the discriminant will shrink sufficiently,
					// we can start the next loop iteration:
					if (a.exponent(i) < 2*(degree()-h.degree()+1)) {
						continue;
					}
				}
			}
			// and use the Buchmann/Lenstra - algorithm for further
			// maximization, if necessary
			if (a.base(i) <= degree()) {
				debug_handler_c("order", "in member - function maximize()::"
						"next  call to maximize(p)", 4,
						std::cout << "with p = " << a.base(i) << " exponent is "
						<< a.exponent(i));
				maximized.assign(maximized.maximize(a.base(i)));
				// i.e. the order is maximized locally at a.base(i)
				debug_handler("order", "in member - function maximize()::"
					      "returned from maximize(p)");
			}
			else {
				multiply(product, product, a.base(i));
			}
		}
	if (!product.is_one()) {
		maximized.assign(maximized.maximize(product));
	}
	return maximized;
}



lidia_size_t order::no_of_real_embeddings() const
{
	if (base->f.is_zero()) {
		if (base->table.get_no_of_columns() == 1)
			lidia_error_handler("number_field", "no_of_real_embeddings()::applied on "
					    "non-initialized number_field!");
		else
			base->compute_base();
	}
	return base->real_roots;
}



lidia_size_t order::no_of_roots_of_unity()
{
	if (no_of_real_embeddings() > 0)
		return 2;

	lidia_size_t n = degree();
	lidia_size_t prec = 500;
	bigint_lattice lat(n, n);
	bigint tmp;
	lidia_size_t i, j;
	for (i = 0; i < n; i++)
		for (j = 0; j < n; j++) {
			(get_conjugate(i, j+1) << prec).bigintify(tmp);
			lat.sto(j, i, tmp);
		}
	lat.set_basis_flag();

	matrix< bigint > result;
	bigrational limit = 16;
	limit.invert();
	add(limit, limit, n/2);
	shift_left(limit, limit, 2*prec);
	lat.short_vector_coeff(result, limit, 20);

	int w = 0;
	for (i = 0; i < result.get_no_of_rows(); i++) {
		bigrational N = norm(alg_number(result.get_row_vector(i), 1, base));
		if (N.is_one() || (-N).is_one()) {
			w++;
		}
	}
	return w*2;
}



// In-/Output:
std::ostream& operator << (std::ostream &s, const order &O) {
	debug_handler("order", "in operator << (std::ostream &, const order &)");
	s << *(O.base) << std::flush;
	return s;
}



std::istream& operator >> (std::istream &s, order &input) {
	debug_handler("order", "in operator >> (std::istream &, order &)");
	input.base->dec_ref();
	input.base = new nf_base;
	input.base->inc_ref();
	s >> *(input.base);
	nf_base::current_base = input.base;
	return s;
}



#ifdef LIDIA_NAMESPACE
}	// end of namespace LiDIA
#endif
