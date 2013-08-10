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
#include	<cctype>



#ifdef LIDIA_NAMESPACE
namespace LiDIA {
#endif



#ifdef LIDIA_IMPLICIT_CAST_EXPLICIT
#define alg_ideal_cast(O) alg_ideal(O)
#else
#define alg_ideal_cast(O) O
#endif

#ifdef LIDIA_DEBUG
int module::count = 0;
#endif

// Constructors & destructor:
// WILL BE REMOVED IN NEXT RELEASE:
module::module(const alg_number & a, const alg_number & b)
	: base(a.degree(), 1),
	  den(a.den),
	  O(a.which_base()),
	  is_exp(true)
{
	debug_handler_c("module", "in module(const a_n, const a_n)", 1,
			count++;
			std::cout << "\nNow we have " << count << " modules!\n");
	O->inc_ref();
	if (b.is_zero())
		base.sto_column_vector(a.coeff, a.degree(), 0);
	else if (a.is_zero()) {
		base.sto_column_vector(b.coeff, b.degree(), 0);
		den.assign(b.den);
	}
	else {
		bigint d = gcd(a.den, b.den);
		bigint e = b.den / d;
		den *= e;
		base.set_no_of_columns(2);
		base.sto_column_vector((a.coeff) * e, a.degree(), 0);
		base.sto_column_vector((b.coeff) * (a.den / d), a.degree(), 1);
	}
	debug_handler_c("module", "in module(const a_n, const a_n)", 3,
			std::cout << " Now the module is " << (*this));
	multiply(*this, *this, alg_ideal_cast(order(static_cast<const nf_base *>(O))));
	debug_handler_c("module", "in module(const a_n, const a_n)", 3,
			std::cout << " So the ideal is " << (*this));
}



module::module(const matrix< bigint > & B, const bigint & d,
               const nf_base * O1)
	: base(B),
	  den(d),
	  O(const_cast<nf_base *>(O1)),
	  is_exp(false)
{
	debug_handler_c("module",
			"in module(const matrix< bigint > &, const bigint &, "
			"const nf_base *)", 1,
			count++;
			std::cout << "\nNow we have " << count << " modules!\n");
	O->inc_ref();
	if (B.get_no_of_rows() != O1->degree()) {
		lidia_error_handler("module", "module(const matrix< bigint > &, ....):: "
				    "Dimension of matrix does not match degree of order");
		assign(module());
	}
	base.image(base);
	normalize();
	// WILL BE REMOVED IN NEXT RELEASE:
	if (B.get_no_of_columns() != O1->degree())
		multiply(*this, *this, alg_ideal_cast(order(static_cast<const nf_base *>(O))));
}



module::module(const bigmod_matrix & B, const bigint & d,
               const nf_base * O1)
	: base(B),
	  den(d),
	  O(const_cast<nf_base *>(O1)),
	  is_exp(false)
{
	debug_handler_c("module",
			"in module(const bigmod_matrix &, const bigint &, "
			"const nf_base *)", 1,
			count++;
			std::cout << "\nNow we have " << count << " modules!\n");
	O->inc_ref();
	if (B.get_no_of_rows() != O1->degree()) {
		lidia_error_handler("module", "module(const bigmod_matrix &, ....):: "
				    "Dimension of matrix does not match degree of order");
		assign(module());
	}
	base.image(base);
	normalize();
}



module::module(const nf_base * O1)
	: base(O1->degree(), 1),
	  den(1),
	  O(const_cast<nf_base *>(O1)),
	  is_exp(true)
{
	debug_handler_c("module", "in module(const nf_base *)", 1,
			count++;
			std::cout << "\nNow we have " << count << " modules!\n");
	O->inc_ref();
}



module::module(const module & M)
	: base(M.base),
	  den(M.den),
	  O(M.O),
	  is_exp(M.is_exp)
{
	debug_handler_c("module", "in module(const module &)", 1,
			count++;
			std::cout << "\nNow we have " << count << " modules!\n");
	O->inc_ref();
}



module::~module()
{
	debug_handler_c("module", "in ~module()", 1,
			count--;
			std::cout << "\nNow we have only " << count << " modules!\n");
	O->dec_ref();
}



// member-functions
bool module::is_zero() const
{
	//   if (base.is_matrix_zero() && base.get_no_of_columns() > 1)
	//     lidia_error_handler_c("module","is_zero()",
	//                           std::cout <<"\nNon standard representation of zer "
	//                           "found for "<<base<<std::endl);
	return (base.get_modulus().is_zero()
		&& base.get_no_of_columns() == 1 && base.is_matrix_zero());
}



bool module::is_whole_order() const
{
	if (!(den.is_one()))
		return false;
	if (base.get_modulus().is_one())
		return true;
	lidia_size_t n = base.get_no_of_rows();
	if (base.get_no_of_columns() != n)
		return false;
	bigint_matrix A(n, n);
	A.diag(1, 0);
	return (A == base);
}



void module::normalize()
{
	if (den.is_negative())
		den.negate();
	bigint d = gcd(den, base.get_modulus());
	register lidia_size_t i, j;
	for (i = 0; !(d.is_one()) && i < degree(); i++)
		for (j = 0; !(d.is_one()) && j < base.get_no_of_columns(); j++)
			d = gcd (base.member(i, j), d);
	if (!(d.is_one())) {
		den /= d;
		base /= d;
		base.reduce(base.get_modulus()/d);
	}
}



void module::invert()
{
	debug_handler("module", "in member - function invert()");
	lidia_size_t n = base.get_no_of_columns();
	lidia_size_t m = base.get_no_of_rows();

	if (base.get_modulus().is_zero() && n != m) {
		lidia_error_handler("module", "invert():: Sorry, I can only "
				    "invert modules of full rank");
		return;
	}
	bigint d = den;
	den.assign_one(); //Multiply *this by den !!

	//Now we have an integral Z -- module
	// Compute and save the norm
	bigint N = exponent(*this).numerator(); //it's anyway an integer!!
	// Compute: O/(exp(M)) ---> Hom( M/(exp(M)) --> O/(exp(M)) )
	//              a    ---> (       b     -->    a b   )

	alg_number a, b, y;
	base.reduce(N);

	n = base.get_no_of_columns();
	bigmod_matrix Map(n*m, m, N);
	bigint * tmp = new bigint[m], rem;

	register lidia_size_t i;

	base_vector< bigint > tmp_vec(a.degree(), a.degree());
	for (i = 0; i < m; tmp[i++].assign_zero()) {
		tmp[i].assign_one(); //tmp = e_i;
		// compute the images of the basis of O
		a = alg_number(tmp, 1, O);
		//a=e_i
		for (register lidia_size_t j = 0; j < n; j++) {
                        // compute the images of the basis of M
                        // under the MULT-with-a-homomorphism
			base.get_column_vector(tmp_vec, j);
			b = alg_number(tmp_vec, 1, O);
			multiply(y, a, b); // the image under the MULT-with-a-homomorphism
			tmp_vec = y.coeff_vector();
			Map.sto_column_vector(y.coeff_vector(), m, i, j*m);
		}
	}
	delete[] tmp;

	// Computing the kernel of this matrix and selecting the right rows
	debug_handler_c("module", "invert", 2,
			std::cout << "Computing kernel of " << Map << std::endl);
	base.kernel(Map);
	debug_handler_c("module", "invert", 2,
			std::cout << "Kernel is " << base << std::endl);
	den.assign(N);
	base.image(base);
	if (!(d.is_one())) multiply(*this, *this, d);
	normalize();
	is_exp = false;
	debug_handler_c("module", "invert", 2,
			std::cout << "so the inverse is " << (*this) << std::endl);
}



void module::assign_zero()
{
	base.set_no_of_columns(1);
	base.set_modulus(0);
	for (lidia_size_t i = 0; i < degree(); i++)
		base.sto(i, 0, bigint(0));
	den.assign_one();
	is_exp = true;
}



// WILL BE REMOVED IN NEXT RELEASE:
void module::assign(const bigint & a)
{
	debug_handler("module", "in member - function assign(const bigint &)");
	base.set_no_of_rows(degree());
	base.set_no_of_columns(1);
	for (lidia_size_t i = 0; i < degree(); i++)
		base.sto(i, 0, bigint(0));
	base.set_modulus(abs(a));
	den.assign_one();
	is_exp = true;
}



// WILL BE REMOVED IN NEXT RELEASE:
void module::assign(const alg_number & a)
{
	debug_handler("module", "in member - function assign(a_n &)");
	base.set_no_of_columns(1);
	base.set_no_of_rows(a.degree());
	base.set_modulus(0);
	if (O != a.O) {
		O->dec_ref();
		O = a.O;
		O->inc_ref();
	}
	base.sto_column_vector(a.coeff_vector(), degree(), 0);
	den = a.den;
	is_exp = false;
	multiply(*this, *this, alg_ideal_cast(order(static_cast<const nf_base *>(O))));
}



void module::assign(const module & a)
{
	debug_handler("module", "in member - function assign(module &)");
	if (O != a.O) {
		O->dec_ref();
		O = a.O;
		O->inc_ref();
	}
	base = a.base;
	den = a.den;
	is_exp = a.is_exp;
}



// Procedural versions:
void add(module & c, const module & a, const module & b)
{
	debug_handler("module", "in function "
		      "add(module &, const module &, const module &)");
	if (a.O != b.O) {
		lidia_error_handler("module", "add(...):: Addition of modules from "
				    "different orders is not yet implemented");
		return;
	}

	if (c.O != a.O) {
		c.O->dec_ref();
		c.O = a.O;
		c.O->inc_ref();
	}
	c.is_exp = false;

	bigint d = gcd(a.den, b.den);
	bigint e = (a.den.is_zero())? bigint(1) : a.den / d;
	bigint f = (b.den.is_zero())? bigint(1) : b.den / d;

	debug_handler_c("module", "add", 0,
			std::cout << "Adding " << a << "(" << a.base << ") and ";
			std::cout << b << "(" << b.base << ")" << std::flush);

	bigmod_matrix
		multhelp1(a.base.get_no_of_rows(), a.base.get_no_of_columns(),
			  a.base.get_modulus() * f),
		multhelp2(a.base.get_no_of_rows(), b.base.get_no_of_columns(),
			  b.base.get_modulus() * e);

	debug_handler_c("module", "add", 5,
			std::cout << "multiply a by " << f << " and b by " << e);

	multiply(multhelp1, a.base, f);
	multiply(multhelp2, b.base, e);

	debug_handler_c("module", "add", 5,
			std::cout << "results are " << multhelp1 << " and " << multhelp2);


	d.assign(gcd(multhelp1.get_modulus(), multhelp2.get_modulus()));
	multhelp1.reduce(d);
	multhelp2.reduce(d);

	debug_handler_c("module", "add", 5,
			std::cout << "reduced results are " << multhelp1 << " and " << multhelp2);
	debug_handler_c("module", "in add(...) -- tmp should be a ", 0,
			std::cout << a.degree() << "x" << multhelp1.get_no_of_columns()
			+ multhelp2.get_no_of_columns() << "matrix");
	bigmod_matrix tmp(a.degree(), multhelp1.get_no_of_columns()
			  + multhelp2.get_no_of_columns(), d);
	debug_handler_c("module", "in add(...) -- tmp is", 0, std::cout << tmp);
	tmp.compose_h(multhelp1, multhelp2);
	//  i.e. tmp.compose_h(a.base * f, b.base * e);
	debug_handler_c("module", "in add(...) -- tmp is", 0, std::cout << tmp);

	c.base.image(tmp);
	debug_handler_c("module", "in add(...) -- reduced tmp is", 0, std::cout << c.base);

	multiply(c.den, e, b.den);
	c.normalize();
}



void intersect(module & c, const module & a, const module & b)
{
	debug_handler("module", "in function intersect(module, "
		      "const & module, const & module)");
	if (a.O != b.O) {
		lidia_error_handler("module", "intersect(...):: "
				    "Intersection of modules from "
				    "different orders is not yet implemented");
		return;
	}

	debug_handler_c("module", "intersect", 5,
			std::cout << "process" << a << " and " << b;
			std::cout << "represented by " << a.base << " and " << b.base);

	if (c.O != a.O) {
		c.O->dec_ref();
		c.O = a.O;
		c.O->inc_ref();
	}
	c.is_exp = false;

	bigint d = gcd(a.den, b.den);
	bigint e = a.den / d;
	bigint f = b.den / d;
	c.base.set_no_of_rows(a.degree());

	bigmod_matrix kern,
		tmp1(a.base.get_no_of_rows(), a.base.get_no_of_columns(),
		     a.base.get_modulus() * f),
		tmp2(a.base.get_no_of_rows(), b.base.get_no_of_columns(),
		     b.base.get_modulus() * e);

	debug_handler_c("module", "intersect", 5,
			std::cout << "multiply a by " << f << " and b by " << -e);

	multiply(tmp1, a.base, f);
	multiply(tmp2, b.base, -e);

	debug_handler_c("module", "intersect", 5,
			std::cout << "results are " << tmp1 << " and " << tmp2);

	d = lcm(a.base.get_modulus()*f, b.base.get_modulus()*e);
	tmp1.lift(d);
	tmp2.lift(d);

	debug_handler_c("module", "intersect", 5,
			std::cout << "Lifting mod " << d
			<< " results in " << tmp1 << " and " << tmp2);

	bigmod_matrix tmp(a.degree(), tmp1.get_no_of_columns()
			  + tmp2.get_no_of_columns(), d);
	debug_handler_c("module", "intersect", 5,
                        std::cout << "compose " << tmp1 << " and " << tmp2);
	tmp.compose_h(tmp1, tmp2);
	debug_handler_c("module", "intersect", 5,
                        std::cout << "compute kernel of " << tmp);
	kern.kernel(tmp);
	debug_handler_c("module", "intersect", 5,
                        std::cout << "kernel is " << kern);
	tmp.set_no_of_columns(kern.get_no_of_columns());

	math_vector< bigint > tmp_vec(a.degree(), a.degree());
	for (lidia_size_t i = 0; i < kern.get_no_of_columns(); i++) {
		kern.get_column_vector(tmp_vec, i);
		tmp.sto_column_vector((tmp1 * tmp_vec), a.degree(), i);
	}
	c.base.image(tmp);
	debug_handler_c("module", "intersect", 5,
                        std::cout << "base is " << c.base);
	multiply(c.den, e, b.den);
	c.normalize();
}



#if 0
void multiply2(module &c, const module & a, const module & b)
{
	debug_handler("module","in function multiply2 (module, "
		      "const & module, const & module)");
	c.O->dec_ref();
	c.O = a.O;
	c.O->inc_ref();

	c.base.set_no_of_rows(a.degree());

	bigmod_matrix tmp1;
	tmp1.special_multiply(a.base, b.base, (a.O)->table);
	tmp1.image(tmp1);
	register lidia_size_t i = 0;
	for (register lidia_size_t j = 0; j<tmp1.get_no_of_columns(); j++)
		if (!tmp1.is_column_zero(j)) i++;
	if (i==0){
		c.base.set_no_of_columns(1);
	}
	else {
		c.base.set_no_of_columns(i);
	}
	tmp1.split_h(tmp1, c.base);
	c.den = a.den * b.den;
	c.normalize();
}
#endif



void multiply(module & c, const module & a, const module & b)
{
	debug_handler("module", "in function multiply(module, "
		      "const & module, const & module)");
	if (a.O != b.O) {
		lidia_error_handler("module", "multiply(...):: Multiplication of modules "
				    "from different orders is not yet implemented");
		return;
	}

	if (c.O != a.O) {
		c.O->dec_ref();
		c.O = a.O;
		c.O->inc_ref();
	}
	if (a.is_zero() || b.is_zero()) {
		c.assign_zero();
		return;
	}
	c.is_exp = false;

	c.base.set_no_of_rows(a.degree());


	// Obvious improvement: First lift the original moduls, if they don't have
	// no_of_columns == rank and avoid appending b*a.modulus()*I,a*b.modulus()*I.

	alg_number n1, n2;
	bigint d;
	multiply(d, a.base.get_modulus(), b.base.get_modulus());
	lidia_size_t column_no1
		= a.base.get_no_of_columns()*b.base.get_no_of_columns();
	lidia_size_t column_no = column_no1;
	lidia_size_t column_no2 = column_no1;
	if (!(a.base.get_modulus().is_zero()))
		column_no2 =
			(column_no += a.base.get_no_of_rows() * b.base.get_no_of_columns());
	if (!(b.base.get_modulus().is_zero()))
		column_no += b.base.get_no_of_rows() * a.base.get_no_of_columns();

	bigmod_matrix tmp1(a.degree(), column_no, d);

	register lidia_size_t i;

	debug_handler_c("module", "multiply(...)", 1,
			std::cout << a << " and " << b;
			std::cout << "represented by" << a.base << " and " << b.base);

	math_vector< bigint > tmp_vec(a.degree(), a.degree());
	for (i = 0; i < a.base.get_no_of_columns(); i++) {
		a.base.get_column_vector(tmp_vec, i);
		alg_number n1(tmp_vec, 1, a.which_base());
		for (register lidia_size_t j = 0; j < b.base.get_no_of_columns(); j++) {
			b.base.get_column_vector(tmp_vec, j);
			alg_number n2(tmp_vec, 1, b.which_base());
			debug_handler_c("module", "multiply::multiply elements", 0,
					std::cout << i << " " << j << " " << i*b.base.get_no_of_columns()+j);
			tmp1.sto_column_vector((n1*n2).coeff_vector(),
					       b.degree(), i*b.base.get_no_of_columns()+j);
		}
	}

	// append b*a.modulus()*I
	bigint * tmp2 = new bigint[b.base.get_no_of_rows()];
	if (!(a.base.get_modulus().is_zero())) {
		d.assign(a.base.get_modulus());
		for (i = 0; i < b.base.get_no_of_columns(); i++) {
			b.base.get_column_vector(tmp_vec, i);
			alg_number n1(tmp_vec, 1, b.which_base());
			for (register lidia_size_t j = 0; j < a.base.get_no_of_rows();
			     tmp2[j++].assign_zero()) {
				tmp2[j].assign(d);
				alg_number n2(tmp2, 1, a.which_base());
				debug_handler_c("module", "multiply::multiply elements", 0,
						std::cout << n1 << " " << n2 << " " << column_no1+i*a.base.get_no_of_rows()+j);
				tmp1.sto_column_vector((n1*n2).coeff_vector(), b.degree(),
						       column_no1 + i * a.base.get_no_of_rows() + j);
			}
		}
	}

	// append a*b.modulus()*I
	if (!(b.base.get_modulus().is_zero())) {
		d.assign(b.base.get_modulus());
		for (i = 0; i < a.base.get_no_of_columns(); i++) {
			a.base.get_column_vector(tmp_vec, i);
			alg_number n1(tmp_vec, 1, a.which_base());
			for (register lidia_size_t j = 0; j < b.base.get_no_of_rows();
			     tmp2[j++].assign_zero()) {
				tmp2[j].assign(d);
				alg_number n2(tmp2, 1, b.which_base());
				debug_handler_c("module", "multiply::multiply elements", 0,
						std::cout << n1 << " " << n2 << " " << column_no2 + i*b.base.get_no_of_rows()+j);
				tmp1.sto_column_vector((n1*n2).coeff_vector(), b.degree(),
						       column_no2 + i * b.base.get_no_of_rows() + j);
			}
		}
	}
	delete[] tmp2;

	c.base.image(tmp1);
	c.den = a.den * b.den;
	c.normalize();
	debug_handler_c("module", "multiply::end of call", 1,
			std::cout << "Product is " << c << "represented by " << c.base << std::endl);
}



void multiply(module &c, const module &a, const bigint & b)
{
	bigint d = gcd(a.den, b);
	if (c.O != a.O) {
		c.O->dec_ref();
		c.O = a.O;
		c.O->inc_ref();
	}

	if (!(b.is_zero())) {
		bigint e;
		divide(e, b, d);
		e.absolute_value();
		c.is_exp = a.is_exp;
		c.den = a.den / d; // Might change value of b, if &b == &c.den !!
		c.base.set_modulus(a.base.get_modulus()*e);
		multiply(c.base, a.base, e);
	}
	else c.assign_zero();
}



void multiply(module &c, const bigint &b, const module &a)
{
	bigint d = gcd(a.den, b);

	if (c.O != a.O) {
		c.O->dec_ref();
		c.O = a.O;
		c.O->inc_ref();
	}
	if (!(b.is_zero())) {
		bigint e;
		divide(e, b, d);
		e.absolute_value();
		c.is_exp = a.is_exp;
		c.den = a.den / d; // Might change value of b, if &b == &c.den !!
		c.base.set_modulus(a.base.get_modulus()*e);
		multiply(c.base, a.base, e);
	}
	else c.assign_zero();
}



void divide(module &c, const module &a, const bigint & b)
{
	debug_handler("alg_numbers", "in function divide(module &, "
		      "const module &, const bigint &)");
	if (c.O != a.O) {
		c.O->dec_ref();
		c.O = a.O;
		c.O->inc_ref();
	}
	c.is_exp = a.is_exp;

	c.base = a.base;

	if (&c != &a) {
		c.den = b;
		c.normalize();
		c.den *= a.den;
	}
	else {
		bigint tmp = a.den;
		c.den = b;
		c.normalize();
		c.den *= tmp;
	}
}



void divide(module & c, const module & a, const module & bb)
{
	debug_handler_c("module", "in function divide(module &, const module &, "
			"const module &", 2, std::cout << "Divide " << a << " by " << bb);
	if (a.O != bb.O) {
		lidia_error_handler("module", "divide(...):: Division of modules from "
				    "different orders is not yet implemented");
		return;
	}
	if (bb.base.get_modulus().is_zero() &&
	    bb.base.get_no_of_columns() != bb.base.get_no_of_rows()) {
		lidia_error_handler("module", "divide():: Sorry, I can only "
				    "divide by modules of full rank");
		return;
	}

	if (c.O != a.O) {
		c.O->dec_ref();
		c.O = a.O;
		c.O->inc_ref();
	}

	if (a.is_zero()) {
		c.base.set_no_of_rows(a.base.get_no_of_rows());
		c.assign_zero();
		return;
	}
	c.is_exp = false;

	module b(bb);

	bigint d = b.den;
	b.den.assign_one(); // Multiply b by den !!

	//Now 'b' is an integral Z -- modules.
	// Compute and save the exponent of B (which is an integer !!)
	bigint N = exponent(b).numerator();
	b.base.lift(0);
	bigmod_matrix tmp(a.base);
	tmp.lift(0);

	lidia_size_t m = tmp.get_no_of_columns();
	lidia_size_t n = b.base.get_no_of_columns();

	//Now we consider two integral Z -- modules, where 'b' has full rank.
	// Compute: A ---> Hom( B --> O / (exp(B) A) )
	//            a ---> (    b -->      a b       )

	alg_number x, y, z;

	bigint_matrix Map(a.degree() * n, m * (n+1));

	register lidia_size_t i;

	base_vector< bigint > tmp_vec(a.degree(), a.degree());
	for (i = 0; i < m; i++) {
		// compute the images of the basis of A
		tmp.get_column_vector(tmp_vec, i);
		x = alg_number(tmp_vec, 1, a.O);
		for (register lidia_size_t j = 0; j < n; j++) {
			// compute the images of the basis of B
			// under the MULT-with-x-homomorphism
			b.base.get_column_vector(tmp_vec, j);
			y = alg_number(tmp_vec, 1, b.O);
			multiply(z, x, y); // the image under the MULT-with-a-homomorphism
			tmp_vec = z.coeff_vector();
			//  Put the image in the matrix of the map
			Map.sto_column_vector(tmp_vec, a.degree(), i, j*a.degree());
		}
	}

	// Complete initialization of map by use of a
	// (n*degree())x(m*n) diagonal matrix,
	// containing N(B)* A ... on the diagonal.
	if (N.is_negative()) {
		lidia_error_handler_c("module", "divide",
				      std::cout << "Norm of B < 0 for B = " << b << ": " << N << std::endl);
		return;
	}
	bigint_matrix NtimesA;
	multiply(NtimesA, *(static_cast<bigint_matrix *>(&tmp)), N);

	register lidia_size_t j = m;
	for (i = 0; i < n * m; i++, j++)
		for (register lidia_size_t k = 0; k < a.degree(); k++)
			Map.sto(k + (i/n) * n, j, NtimesA.member(k, i % m));

	// Computing the kernel of this matrix and selecting the right rows
	debug_handler_c("module", "divide", 2,
			std::cout << "Computing kernel of " << Map << std::endl);
	bigint_matrix tmp2;
	tmp2.kernel(Map);
	debug_handler_c("module", "divide", 2,
			std::cout << "Verify kernel: " << Map*tmp2 << std::endl);
	tmp2.set_no_of_rows(m);
	debug_handler_c("module", "divide", 2,
			std::cout << "Kernel is " << tmp2 << std::endl);
	bigint reduction_para;
	multiply(reduction_para, a.base.get_modulus(), N);
	c.base.set_modulus(0);
	multiply(*(static_cast<bigint_matrix *>(&c.base)),
		 *(static_cast<bigint_matrix *>(&tmp)), tmp2); //bigint_matrix - Multiplikation !!
	if (reduction_para.is_zero())
		c.base.image(c.base);
	else
		c.base.reduce(reduction_para);
	debug_handler_c("module", "divide", 2,
			std::cout << "so c.base is " << c.base << std::endl);
	if (!(a.den.is_one()))
		multiply(c.den, a.den, N);
	else
		c.den.assign(N);
	if (!(d.is_one())) multiply(c, c, d);
	c.normalize();
	debug_handler_c("module", "divide", 2,
			std::cout << "so c is " << c << std::endl);
	debug_handler_c("module", "divide", 2,
			std::cout << "Verify: " << bb*c << a << std::endl);
}



void remainder(module &c, const module &a, const bigint &p)
{
	if (!(a.den.is_one())) {
		lidia_error_handler("module", "remainder:: Reducing mod p for a "
				    "fractional ideal ??");
		return;
	}
	c.assign(a);
	bigint new_mod = a.base.get_modulus() %p;
	if (new_mod.is_zero()) new_mod.assign(p);
	c.base.set_modulus(new_mod << 1);
	c.base.reduce(new_mod);
	c.base.image(c.base);
}



void power(module & c, const module & a, const bigint & b)
{
	debug_handler("module", "in function power(module &, "
		      "const module &, const bigint &)");
	bigint expo;
	nf_base * O = a.which_base();
	module multiplier(O);

	if (c.O != a.O) {
		c.O->dec_ref();
		c.O = a.O;
		c.O->inc_ref();
	}
	if (b.is_negative())
		power(c, inverse(a), -b);
	else if (b.is_zero() || a.is_one())
		c.assign_one();
	else {
		expo.assign(b);
		multiplier.assign(a);
		c.assign_one();
		while (expo.is_gt_zero()) {
			if (!(expo.is_even()))
				multiply(c, c, multiplier);
			square(multiplier, multiplier);
			expo.divide_by_2();
		}
	}
}



// Comparision:
// By now, only comparision of modules over the same order is implemented.
// One function to compute the needed information:

void module::compare(const module & a,
                     bool & equal,
                     bool & a_in_this) const
{
	debug_handler("module", "in member - function compare(...)");
	if (equal = (this == &a))return;
	if (O != (a.O)) {       //we should use *O != *(a.O)
		lidia_error_handler("module", "compare::You tried to compare modules over "
				    "different bases!");
		return;
	}

	debug_handler_c("module", "compare", 9,
			std::cout << "comparing" << (*this) << " and " << a;
			std::cout << "represented by " << base << " and " << a.base);
	debug_handler_c("module", "compare", 9,
			std::cout << " and these are true exponents: ";
			std::cout << is_exp << " " << a.is_exp << std::endl);

	if (den != a.den) {
		bigint d = gcd(den, a.den);
		((*this)*(a.den/d)*den).compare(a*(den/d)*a.den, equal, a_in_this);
		return;
	}

	// First compute exponents
	bigint expo = base.get_modulus(),
		aexpo = a.base.get_modulus();

	if (!is_exp) {
		base.exponent(expo);
		base.reduce(expo);
		is_exp = true;
	}
	if (!a.is_exp) {	// Wouldn't change, if &a = this, so check
		a.base.exponent(aexpo); // at beginning is necessary!!!
		a.base.reduce(aexpo);
		a.is_exp = true;
	}

	// Check for equality - also useful, in case you're not checking
	// for equality since either fast or sufficient
	debug_handler_c("module", "compare", 9,
			std::cout << "Reduced to " << (*this) << " and " << a;
			std::cout << "represented by " << base << " and " << a.base);
	bigmod_matrix cmp1(base), cmp2(a.base);
	if (aexpo == expo) {
		if (equal = (den == a.den && base == a.base)) return;
		if (den == a.den) {
			cmp1.unique_image(cmp1);
			cmp2.unique_image(cmp2);
			if (equal = (cmp1 == cmp2)) return;
		}
	}
	else equal = false;

	//Check for inclusion - useless, if you check for (in-)equality!!!
	bigint rem;
	if (expo.is_zero()) {
		if (!aexpo.is_zero()) {
			a_in_this = false;
			return;
		}
	}
	else {
		remainder(rem, aexpo, expo);
		if (!(rem.is_zero())) {
			a_in_this = false;
			return;
		}
	}
	cmp1.lift(0);
	bigint_matrix T;
	debug_handler_c("module", "compare", 9,
			std::cout << "Trying to solve " << cmp1 << a.base << std::endl);
	if (cmp1.get_no_of_columns() == cmp1.get_no_of_rows()) {
		T.reginvimage(cmp1, a.base);
		debug_handler_c("module", "compare", 9,
				std::cout << "Solution is " << T << std::endl);
		a_in_this = true;
		for (lidia_size_t i = 0; a_in_this && (i < T.get_no_of_columns()); i++)
			if (!(T.member(degree(), i).is_one())) a_in_this = false;
	}
	else {
		a_in_this = true;
		for (lidia_size_t i = 0; a_in_this && (i < a.base.get_no_of_columns()); i++) {
			debug_handler_c("module", "compare", 9,
					std::cout << " Calling solve for " << cmp1 << " and ";
					std::cout << a.base.get_column_vector(i) << std::endl);
			T.solve(cmp1, a.base.get_column_vector(i));
			debug_handler_c("module", "compare", 9,
					std::cout << "Solution for column " << i << " is:" << T << std::endl);
			if (T.get_no_of_columns() == 1) a_in_this = false;
		}
	}
}



bool operator == (const module & a, const module & b)
{
	bool equal, b_in_a;
	debug_handler_l("module", "operator == ", 2);
	a.compare(b, equal, b_in_a);
	return equal;
}



bool operator != (const module & a, const module & b)
{
	bool equal, b_in_a;
	a.compare(b, equal, b_in_a);
	return !equal;
}



bool operator <= (const module & a, const module & b) // Same as divisibility!!!!
{
	bool equal, a_in_b;
	b.compare(a, equal, a_in_b);
	return (equal || a_in_b);
}



bool operator < (const module & a, const module & b)
{
	bool equal, a_in_b;
	b.compare(a, equal, a_in_b);
	return (!equal && a_in_b);
}



bool operator >= (const module & a, const module & b)
{
	bool equal, b_in_a;
	a.compare(b, equal, b_in_a);
	return (equal || b_in_a);
}



bool operator > (const module & a, const module & b)
{
	bool equal, b_in_a;
	a.compare(b, equal, b_in_a);
	return (!equal && b_in_a);
}



// Some number-theoretic functions:
bigrational norm(const module & a)      // Norm
{
	if (a.base.get_modulus().is_zero() &&
	    a.base.get_no_of_columns() != a.degree()) {
		return 0;
	}
	bigmod_matrix tmp (a.base);
	tmp.lift(0);
	bigint normdenominator;
	power(normdenominator, a.den, bigint(a.degree()));
	return bigrational((bigint_matrix(tmp)).det(), normdenominator);
}



bigrational exponent(const module & a)
{

	debug_handler_c("module",
			"in function exp(const & module)", 3,
			std::cout << "Computing exponent of " << a << std::endl);
	if (!a.is_exp) {
		bigint e(a.base.exponent());
		a.base.reduce(e);
		a.is_exp = true;
	}
	return bigrational(a.base.get_modulus(), a.den);
}



order module::ring_of_multipliers(const bigint &p) const
{
	// We are here in a member function of `Ip'
	// Consider the kernel C of the map:
	// O/pO -> End(Ip/pIp)
	//   a  -> (b -> ab)
	// then the result O' fullfills: O' =1/p * C
	// quite similar to `pseudo-radical', but more complicated, since
	// now we have matrices for each (b -> ab).

	debug_handler("module", "in member - function ring_of_multipliers(const bigint &, bigint &)");
	register lidia_size_t i, j, k;
	register lidia_size_t n = base.get_no_of_columns();
	register lidia_size_t m = degree();
	alg_number a, b, y;

	bigint_matrix init(m, n*m);
	bigmod_matrix B(base);
	bigmod_matrix C;
	bigmod_matrix VW(n * m, m, p);
	bigint_matrix v;
	bigint * tmp = new bigint[m];
	bigint rem;

	if (base.get_modulus().is_zero() && n != base.get_no_of_rows()) {
		lidia_error_handler("module", "ring_of_multipliers(const bigint &):: "
				    "Sorry, I can compute this only for "
				    "modules of full rank");
		delete [] tmp;
		return order();
	}

	math_vector< bigint > tmp_vec(a.degree(), a.degree());
	for (i = 0; i < m; tmp[i++].assign_zero()) {
		tmp[i].assign_one(); //tmp = e_i;
		// compute the images of the basis of O/pO
		a = alg_number(tmp, 1, O);
		//a=e_i
		for (j = 0; j < n; j++) {
			// compute the images of the basis of A+pO
			// under the MULT-with-a-homomorphism
			B.get_column_vector(tmp_vec, j);
			b = alg_number(tmp_vec, 1, O);
			y = a*b; // the image under the MULT-with-a-homomorphism
			init.sto_column_vector(y.coeff_vector(), m, i*n+j);
		}
	}
	delete[] tmp;

	// the image must be written relative to base;
	B.lift(0);
	debug_handler_c("module", "in member - function ring_of_multipliers(...)",
			1, std::cout << "solve " << B << init << std::endl << std::flush);
	v.reginvimage(B, init);
	debug_handler_c("module", "in member - function ring_of_multipliers(...)",
			1, std::cout << "solution is " << v << std::endl << std::flush);
	// move the result to the Matrix VW
	for (i = 0; i < m; i++) {
		for (j = 0; j < n; j++) {
			for (k = 0; k < m; k++) {
				VW.sto(k+j*m, i, v.member(k, i*n+j));
			}
		}
	}
	debug_handler_c("module", "in member - function ring_of_multipliers(...)",
			1, std::cout << "Compute kernel of" << VW);
	C.kernel(VW);
	debug_handler_c("module", "in member - function ring_of_multipliers(...)",
			1, std::cout << "Kernel is " << C << std::flush);

	// interpret C as a module M and lift it into the order!!
	init.set_no_of_columns(k = C.get_no_of_columns());
	for (j = 0; j < n; j++) {
		C.get_row_vector(tmp_vec, j);
		init.sto_row_vector(tmp_vec, k, j);
	}
	module M(init, 1, O);
	M += p * order(static_cast<const nf_base *>(O)); // hopefully interpreted as module-mult.
	debug_handler_c("module", "in member - function ring_of_multipliers(...)",
			4, std::cout << "Module is" << M << std::flush);

	//  Instead of creating the multiplication table, we create the base:
	if (!(O->f.is_zero())) {
		M.den *= p;
		M.normalize();
		debug_handler_c("module", "ring_of_multipliers(...)", 4,
				std::cout << "Transformation is  1/" << M.den << " * ";
				std::cout << M.base << std::endl << std::flush);

		order result(M.z_basis(), M.den, O);
		debug_handler_c("module", "ring_of_multipliers(...)", 4,
				std::cout << "So the new order has the following table";
				std::cout << result << std::endl << std::flush;
				std::cout << "and the Discriminant is" << disc(result);
				std::cout << std::endl << std::flush);
		return result;
	}
	else {
		// If we can't create the base (because there is no field)
		// create the MT instead.
		base_vector< bigint > tmp;
		n = M.base.get_no_of_rows();
		bigint_matrix help = M.z_basis();
		init.set_no_of_columns((n*(n+1))/2);
		for (i = 0; i < n; i++) {
			help.get_column_vector(tmp, i);
			a = alg_number(tmp, 1, O);
			for (j = 0; j <= i; j++) {
				help.get_column_vector(tmp, j);
				b = alg_number(tmp, 1, O);
				init.sto_column_vector((a*b).coeff_vector(), n, (i*(i+1))/2+j);
			}
		}
		debug_handler_c("module",
				"in member - function ring_of_multipliers(...)", 3,
				std::cout << "solve " << help << init << std::endl << std::flush);

		init.reginvimage(help, init);
		init.set_no_of_rows(n);

		debug_handler_c("module",
				"in member - function ring_of_multipliers(...)", 3,
				std::cout << "p * MT is" << trans(init) << std::flush);

		for (i = 0; i < init.get_no_of_rows(); i++)
			for (j = 0; j < init.get_no_of_columns(); j++) {
				bigint q, r;
				div_rem(q, r, init.member(i, j), p);
				if (!(r.is_zero())) {
					lidia_error_handler("module", "ring_of multipliers::internal error::"
							    "division by p unsuccesful");
					return order();
				}
				init.sto(i, j, q);
			}
		debug_handler_c("module",
				"in member - function ring_of_multipliers(...)", 3,
				std::cout << "MT is" << trans(init) << std::flush);
		return order(trans(init));
	}
}



// Other functions:
void invert(module &c, const module & a)
{
	c.assign(a);
	c.invert();
}



module inverse(const module & a)
{
	module c = a;
	c.invert();
	return c;
}



void square(module & a, const module & b)
{
	debug_handler("modules",
		      "in function square(a_n &, const a_n &)");
	multiply(a, b, b);
}



void swap(module & a, module & b)
{
	swap (a.den, b.den);
	swap(a.base, b.base);
	nf_base * O = a.O;
	a.O = b.O;
	b.O = O;
	bool help = a.is_exp;
	a.is_exp = b.is_exp;
	b.is_exp = help;
}



// random modules
void module::randomize(const bigint & b)
{
	debug_handler("module", "in member-function "
		      "randomize(const bigint &)");
	den.assign_one();
	base.set_no_of_columns(degree());
	base.set_modulus(0);
	bigint tmp;
	for (register lidia_size_t i = 0; i < degree(); i++) {
		tmp.assign(LiDIA::randomize(b));
		base.sto(i, i, tmp);
		if (!tmp.is_zero())
			for (register lidia_size_t j = i + 1; j < degree(); j++)
				base.sto(i, j, LiDIA::randomize(tmp));
		else
			for (register lidia_size_t j = i + 1; j < degree(); j++)
				base.sto(i, j, LiDIA::randomize(b));
	}
	add(den, LiDIA::randomize(b), 1);
	is_exp = false;
	base.image(base);
	normalize();
}



// In-/Output:
std::ostream& operator << (std::ostream & s, const module & a)
{
	debug_handler("module", "in function "
		      "std::ostream& operator << (std::ostream &, const module &)");
	s << a.z_basis();
	if (!(a.is_zero() || a.den.is_one()))
		s << " / " << a.den;
	return s;
}



std::istream& operator >> (std::istream & s, module & a)
{
	debug_handler("module", "in function "
		      "std::istream& operator >> (std::istream &, module &)");
	bigint_matrix tmp;
	if (nf_base::current_base == nf_base::dummy_base) {
		lidia_error_handler ("module", "operator >>::No number_field or order!");
		return s;
	}
	if (a.O != nf_base::current_base) {
		a.O->dec_ref();
		a.O = nf_base::current_base;
		a.O->inc_ref();
	}
	s >> tmp;
	a.base.assign(bigmod_matrix(tmp, 0));
	a.is_exp = false;
	a.base.image(a.base);
	char c;
	do {
		s.get(c);
	} while (isspace(c) && c != '\n');
	if (c == '/') {
		s >> a.den;
	}
	else {
		a.den.assign_one();
		if (c != '\n' && c != '\r')
			s.putback(c);
	}
	a.normalize();
	return s;
}



#ifdef LIDIA_NAMESPACE
}	// end of namespace LiDIA
#endif
