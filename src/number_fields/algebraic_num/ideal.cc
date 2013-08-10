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
#include	"LiDIA/bigfloat_lattice.h"



#ifdef LIDIA_NAMESPACE
namespace LiDIA {
#endif



#ifdef LIDIA_IMPLICIT_CAST_EXPLICIT
#define alg_ideal_cast(O) alg_ideal(O)
#else
#define alg_ideal_cast(O) O
#endif



// Constructors & destructor:
alg_ideal::alg_ideal(const bigint & a, const alg_number & b)
	: module(rep_matrix(b), b.denominator(), b.which_base())
{
	bigint m(gcd(base.get_modulus(), a));
	base.reduce(m);
	base.set_modulus(m);
	debug_handler_c("alg_ideal", "in alg_ideal(const bigint &, const a_n)", 3,
			std::cout << " The ideal is " << (*this));
}



alg_ideal::alg_ideal(const alg_number & a, const alg_number & b)
	: module(rep_matrix(a), a.denominator(), a.which_base())
{
	debug_handler_c("alg_ideal", "in alg_ideal(const a_n, const a_n)", 3,
			std::cout << " The ideal is " << (*this));
	if (!b.is_zero()) {
		if (a.which_base() != b.which_base())
			lidia_error_handler("alg_ideal",
					    "alg_ideal (alg_number, alg_number):: "
					    " numbers must be given relative to the same base.");
		bigmod_matrix tmp(degree(), base.get_no_of_columns()+degree(), 0);
		bigint d = gcd(a.denominator(), b.denominator());
		bigint e = b.denominator() / d;
		den *= e;
		multiply(base, base, e);
		tmp.compose_h(base, rep_matrix(b)*(a.denominator()/d));
		is_exp = false;
		base.image(tmp);
		normalize();
	}
}



alg_ideal::alg_ideal(const matrix< bigint > & B,
		     const bigint & d,
		     const nf_base * O1): module(B, d, O1)
{
}



alg_ideal::alg_ideal(const bigmod_matrix & B,
		     const bigint & d,
		     const nf_base * O1)
	: module(B, d, O1)
{
	if (B.get_no_of_columns() != O1->degree() && B.get_modulus().is_zero())
		multiply(*this, *(static_cast<module *>(this)),
			 alg_ideal_cast(order(const_cast<nf_base *>(O))));
}



alg_ideal::alg_ideal(const nf_base * O1)
	: module(O1)
{
}



alg_ideal::alg_ideal(const alg_ideal & I)
	: module(I)
{
}



alg_ideal::~alg_ideal()
{
}



// Procedural versions:
void multiply(alg_ideal & c, const alg_ideal & a, const alg_ideal & b)
{
	debug_handler("alg_ideal", "in function multiply("
		      "alg_ideal, const & alg_ideal, const & alg_ideal)");
	if (a.O != b.O) {
		lidia_error_handler("alg_ideal", "multiply(...):: Multiplication of ideals "
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

	alg_number n1, n2;
	bigint d;
	multiply(d, a.base.get_modulus(), b.base.get_modulus());

	bigmod_matrix tmp1(a.degree(), (a.base.get_no_of_columns()+1)
			   * (b.base.get_no_of_columns()+1) - 1, d);

	register lidia_size_t i, col_no = 0;

	debug_handler_c("alg_ideal", "multiply(...)", 1,
			std::cout << a << " and " << b;
			std::cout << "represented by" << a.base << " and " << b.base);

	// Using get_column in a loop is a bad idea, since it requires continuously
	// allocating and freeing memory. Better allocate _one_ vector and use it
	// all the time !

	math_vector< bigint > tmp_vec(a.degree(), a.degree());
	for (i = 0; i < a.base.get_no_of_columns(); i++) {
		a.base.get_column_vector(tmp_vec, i);
		alg_number n1(tmp_vec, 1, a.which_base());
		for (register lidia_size_t j = 0; j < b.base.get_no_of_columns(); j++) {
			b.base.get_column_vector(tmp_vec, j);
			alg_number n2(tmp_vec, 1, b.which_base());
			debug_handler_c("alg_ideal", "multiply::multiply elements", 1,
					std::cout << col_no);
			tmp1.sto_column_vector((n1*n2).coeff_vector(), b.degree(), col_no++);
		}
	}

	// append b*a.modulus()
	if (!(a.base.get_modulus().is_zero())) {
		d.assign(a.base.get_modulus());
		for (i = 0; i < b.base.get_no_of_columns(); i++) {
			b.base.get_column_vector(tmp_vec, i);
			debug_handler_c("alg_ideal", "multiply::multiply elements", 0,
					std::cout << tmp_vec << " " << d << " " << col_no);
			tmp1.sto_column_vector((tmp_vec * d), b.degree(), col_no++);
		}
	}

	// append a*b.modulus()
	if (!(b.base.get_modulus().is_zero())) {
		d.assign(b.base.get_modulus());
		for (i = 0; i < a.base.get_no_of_columns(); i++) {
			a.base.get_column_vector(tmp_vec, i);
			debug_handler_c("alg_ideal", "multiply::multiply elements", 0,
					std::cout << tmp_vec << " " << d << " " << col_no);
			tmp1.sto_column_vector((tmp_vec * d), a.degree(), col_no++);
		}
	}
	debug_handler_c("alg_ideal", "multiply::Compute image of", 0,
			std::cout << tmp1);
	c.base.image(tmp1);
	c.den = a.den * b.den;
	c.normalize();
	debug_handler_c("alg_ideal", "multiply::end of call", 1,
			std::cout << "Product is " << c << "represented by " << c.base << std::endl);
}



void square(alg_ideal & c, const alg_ideal & a)
{
	debug_handler("alg_ideal", "in function square("
		      "alg_ideal, const & alg_ideal)");

	if (c.O != a.O) {
		c.O->dec_ref();
		c.O = a.O;
		c.O->inc_ref();
	}
	if (a.is_zero()) {
		c.assign_zero();
		return;
	}
	c.is_exp = false;
	c.base.set_no_of_rows(a.degree());

	alg_number n1, n2;
	bigint d;
	square(d, a.base.get_modulus());

	bigmod_matrix tmp1(1, 1, d);
	if (!(a.base.get_modulus().is_zero()))
		tmp1.set_no_of_columns(a.base.get_no_of_columns() +
				       a.base.get_no_of_columns()*(a.base.get_no_of_columns()+1)/2);
	else
		tmp1.set_no_of_columns(a.base.get_no_of_columns()*(a.base.get_no_of_columns()+1)/2);

	tmp1.set_no_of_rows(a.degree());

	register lidia_size_t i, col_no = 0;

	debug_handler_c("alg_ideal", "square(...)", 1,
			std::cout << a << "represented by" << a.base);

	math_vector< bigint > tmp_vec(a.degree(), a.degree());
	alg_number erg;
	for (i = 0; i < a.base.get_no_of_columns(); i++) {
		a.base.get_column_vector(tmp_vec, i);
		alg_number n1(tmp_vec, 1, a.which_base());
		debug_handler_c("alg_ideal", "square::square element", 1,
				std::cout << col_no);
		square(erg, n1);
		tmp1.sto_column_vector(erg.coeff_vector(), a.degree(), col_no++);
		for (lidia_size_t j = i+1; j < a.base.get_no_of_columns(); j++) {
			a.base.get_column_vector(tmp_vec, j);
			alg_number n2(tmp_vec, 1, a.which_base());
			debug_handler_c("alg_ideal", "square::multiply element", 1,
					std::cout << col_no);
			multiply(erg, n1, n2);
			tmp1.sto_column_vector(erg.coeff_vector(), a.degree(), col_no++);
		}
	}

	// append a*a.modulus()
	if (!(a.base.get_modulus().is_zero())) {
		d.assign(a.base.get_modulus());
		for (i = 0; i < a.base.get_no_of_columns(); i++) {
			a.base.get_column_vector(tmp_vec, i);
			debug_handler_c("alg_ideal", "multiply::multiply elements", 0,
					std::cout << tmp_vec << " " << d << " " << col_no);
			tmp1.sto_column_vector((tmp_vec * d), a.degree(), col_no++);
		}
	}
	debug_handler_c("alg_ideal", "square::Compute image of", 0,
			std::cout << tmp1);
	c.base.image(tmp1);
	square(c.den, a.den);
	c.normalize();
	debug_handler_c("alg_ideal", "square::end of call", 1,
			std::cout << "Square is " << c << "represented by " << c.base << std::endl);
}



void multiply(alg_ideal &c, const alg_ideal &a, const alg_number & b)
{
	debug_handler("alg_ideal", "in function "
		      "multiply(alg_ideal, const & alg_ideal, const & alg_number)");
	if (a.O != b.O) {
		lidia_error_handler("alg_ideal", "multiply(...)::alg_number and alg_ideal"
				    "which are to be multiplied must be over the"
				    "same order!");
		return;
	}
	if (c.O != a.O) {
		c.O->dec_ref();
		c.O = a.O;
		c.O->inc_ref();
	}
	c.is_exp = false;
	lidia_size_t n = a.base.get_no_of_columns();
	bigint d = a.base.get_modulus();
	bigint_matrix multiplier = rep_matrix(b);

	multiply(c.den, a.den, b.den);

	if (!d.is_zero()) {
		c.base.assign(a.base);
		bigint N = abs(det(multiplier));
		lidia_size_t rows = c.base.get_no_of_rows();
		if (rows > n) {
			lidia_size_t columns = rows;
			if (n > 1 || !c.base.is_column_zero(0))
				columns += n;
			c.base.set_no_of_columns(columns);
			c.base.set_modulus(d*N);
			for (register lidia_size_t i = rows; i; i--)
				c.base.sto(rows - i, columns-i, d);
		}
		else
			c.base.set_modulus(d*N);
		bigmod_matrix mod_multiplier(multiplier, d*N);
		multiply(c.base, mod_multiplier, c.base);
	}
	else {
		c.base.set_modulus(0);
		multiply(*(bigint_matrix*)&c.base, multiplier, *(bigint_matrix*)&a.base);
	}
	debug_handler_c("alg_ideal", "in multiply()", 3,
			std::cout << "compute image of " << c.base << std::endl);
	c.base.image(c.base);
	c.normalize();
}



void multiply(alg_ideal &c, const alg_number & b, const alg_ideal &a)
{
	debug_handler("alg_ideal", "in function "
		      "multiply(alg_ideal, const & alg_ideal, const & alg_number)");
	if (a.O != b.O) {
		lidia_error_handler("alg_ideal", "multiply(...)::alg_number and alg_ideal"
				    "which are to be multiplied must be over the"
				    "same order!");
		return;
	}
	if (c.O != a.O) {
		c.O->dec_ref();
		c.O = a.O;
		c.O->inc_ref();
	}
	c.is_exp = false;
	lidia_size_t n = a.base.get_no_of_columns();
	bigint d = a.base.get_modulus();
	bigint_matrix multiplier = rep_matrix(b);

	multiply(c.den, a.den, b.den);

	if (!d.is_zero()) {
		c.base.assign(a.base);
		bigint N = abs(det(multiplier));
		lidia_size_t rows = c.base.get_no_of_rows();
		if (rows > n) {
			lidia_size_t columns = rows;
			if (n > 1 || !c.base.is_column_zero(0))
				columns += n;
			c.base.set_no_of_columns(columns);
			c.base.set_modulus(d*N);
			for (register lidia_size_t i = rows; i; i--)
				c.base.sto(rows - i, columns-i, d);
		}
		else
			c.base.set_modulus(d*N);
		bigmod_matrix mod_multiplier(multiplier, d*N);
		multiply(c.base, mod_multiplier, c.base);
	}
	else {
		c.base.set_modulus(0);
		multiply(*(bigint_matrix*)&c.base, multiplier, *(bigint_matrix*)&a.base);
	}
	c.base.image(c.base);
	c.normalize();
}



void divide(alg_ideal &c, const alg_ideal &a, const alg_number & b)
{
	debug_handler("alg_ideal", "in function "
		      "divide(alg_ideal, const & alg_ideal, const & alg_number)");
	if (a.O != b.which_base()) {
		lidia_error_handler("alg_ideal", "divide(...)::alg_number and alg_ideal"
				    "must be over the same order!");
		return;
	}
	if (c.O != a.O) {
		c.O->dec_ref();
		c.O = a.O;
		c.O->inc_ref();
	}

	bigmod_matrix abase(a.base);
	abase.lift(0);

	c.is_exp = false;
	c.base.set_modulus(0); // Do this after lifting a.base,
                                // since &a==&c is possible.


	bigint d = gcd(a.den, b.den);
	divide(c.den, a.den, d);
	bigint multiplier;
	divide(multiplier, b.den, d);
	bigint old_multiplier(multiplier);

	alg_number divisor;
	multiply(divisor, b, b.den); // forget about denominator

	c.base.set_no_of_columns(abase.get_no_of_columns());
	alg_number * results;
	results = new alg_number[c.base.get_no_of_columns()];
	alg_number tmp(bigint(1), c.O);

	register lidia_size_t i;
	for (i = 0; i < abase.get_no_of_columns(); i++) {
		abase.get_column_vector(tmp.coeff, i);
		divide(results[i], tmp, divisor);
		multiplier = lcm(multiplier, results[i].den);
	}
	for (i = 0; i < abase.get_no_of_columns(); i++) {
		multiply(results[i], results[i], multiplier);
		c.base.sto_column_vector(results[i].coeff, c.degree(), i);
		// Only works without reducing since c.base.modulus is set to zero !!!!
	}
	delete[] results;

	divide(multiplier, multiplier, old_multiplier);
	multiply(c.den, c.den, multiplier);
	c.base.image(c.base);
	c.normalize();
}



void divide(alg_ideal & c, const alg_ideal & a, const module & bb)
{
	debug_handler("alg_ideal", "in function divide(alg_ideal &, const alg_ideal &, "
		      "const module &");
	if (a.O != bb.O) {
		lidia_error_handler("alg_ideal", "divide(...):: Division of ideals from "
				    "different orders is not yet implemented");
		return;
	}

	debug_handler_c("alg_ideal", "divide(...)", 2,
			std::cout << a << " by " << bb;
			std::cout << "represented by" << a.base << " and " << bb.base);

	if (bb.base.get_modulus().is_zero() &&
	    bb.base.get_no_of_columns() != bb.base.get_no_of_rows()) {
		lidia_error_handler("alg_ideal", "divide():: Sorry, I can only "
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
		c.is_exp = true;
		return;
	}
	c.is_exp = false;

	module b(bb);

	bigint d = b.den;
	b.den.assign_one(); // Multiply b by den !!

	//Now we have two integral Z -- ideals of full rank.
	// Compute and save the norm of B (which is an integer !!)
	bigint N = exponent(b).numerator();
	//    bigint Nsicher(N);
	//    power(N,a.den, a.degree());
	//    N = (N * norm(a)).numerator();
	//    multiply(N,N,Nsicher);

	// Compute: A ---> Hom( B --> A/(N(B) A) )
	//          a ---> (    b --> a b + (N(B) A) )

	b.base.lift(0);
	lidia_size_t n = b.base.get_no_of_columns();

	alg_number x, y, z;
	bigint_matrix LGS(n, n*n);

	register lidia_size_t i;

	bigmod_matrix tmp(a.base);
	tmp.lift(0);


	base_vector< bigint > tmp_vec(a.degree(), a.degree());
	for (i = 0; i < n; i++) {
		// compute the images of the basis of A
		tmp.get_column_vector(tmp_vec, i);
		x = alg_number(tmp_vec, 1, a.O);
		for (register lidia_size_t j = 0; j < n; j++) {
			// compute the images of the basis of B
			// under the MULT-with-x-homomorphism
			b.base.get_column_vector(tmp_vec, j);
			y = alg_number(tmp_vec, 1, b.O);
			multiply(z, x, y); // the image under the MULT-with-a-homomorphism
			// Put the image in the matrix of the LGS to compute
			// the representation relative to A
			LGS.sto_column_vector(z.coeff_vector(), n, i*n+j);
		}
	}

	// represent the columns of LGS relative to A.
	debug_handler_c("alg_ideal", "divide", 2,
			std::cout << "Solve " << tmp << LGS << std::endl);

	LGS.reginvimage(tmp, LGS);
	debug_handler_c("alg_ideal", "divide", 2,
			std::cout << "Solution is " << LGS << std::endl);

	// since B now is an integral ideal, the denominators of the solutions
	// must all be 1. Therefore just ignore the last row.

	// Put the columns of LGS into Map.
	bigmod_matrix Map(n*n, n, N); // set ideal to 'N'.
	for (i = 0; i < n * n; i++) {
		for (register lidia_size_t j = 0; j < n; j++) {
			Map.sto(j + (i % n)*n, i / n, LGS.member(j, i));
		}
	}

	// Compute the kernel of this matrix mod 'N'.
	debug_handler_c("alg_ideal", "divide", 2,
			std::cout << "Computing kernel of " << Map << std::endl);
	bigmod_matrix kern(n, n, N);
	kern.kernel(Map);
	debug_handler_c("alg_ideal", "divide", 2,
			std::cout << "Kernel is " << kern;
			std::cout << "Verify kernel: " << Map*kern << std::endl);
	kern.set_no_of_rows(n);

	bigint reduction_para;
	multiply(reduction_para, a.base.get_modulus(), N);
	kern.lift(reduction_para);
	tmp.reduce(reduction_para);

	debug_handler_c("alg_ideal", "divide", 2,
			std::cout << "multiply " << tmp << "\nby" << kern << std::endl);
	multiply(c.base, tmp, kern);
	debug_handler_c("alg_ideal", "divide", 2,
			std::cout << "result is " << c.base << std::endl);
	debug_handler_c("alg_ideal", "divide", 2,
			std::cout << "so c.base is " << c.base << std::endl);

	if (!(a.den.is_one()))
		multiply(c.den, a.den, N);
	else
		c.den.assign(N);
	if (!(d.is_one())) multiply(c, *(module *)&c, d);
	c.normalize();
	debug_handler_c("alg_ideal", "divide", 2,
			std::cout << "so c is " << c << std::endl);
	debug_handler_c("alg_ideal", "divide", 2,
			std::cout << "Verify: " << bb*c << a << std::endl);
}



#if 0
void remainder(alg_ideal &c, const alg_ideal &a, const bigint &p)
{
	bigint_matrix garbage;
	c.O = a.which_base();
	if (!a.den.is_one())
		error_handler("alg_ideal", "remainder:: Reducing mod p for a fractional "
			      "alg_ideal ??");
	c.den.assign_one();

	remainder(garbage, a.base, p);
	garbage.hnf();
	register long i=garbage.rank();
	if (i==0){
		c.base.set_no_of_columns(1);
	}
	else {
		c.base.set_no_of_columns(i);
	}
	garbage.split_h(garbage, c.base);
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

	debug_handler("module","in member - function ring_of_multipliers(const bigint &, bigint &)");
	register lidia_size_t i,j,k;
	register lidia_size_t n=degree();
	alg_number a,b,y;

	bigint_matrix init(n,n*n);
	bigint_matrix B(base);
	bigmod_matrix C;
	bigmod_matrix VW(n * n, n, p);
	bigint_matrix v;
	bigint * tmp = new bigint[n];
	bigint * tmp2;
	bigint rem;

	if (base.get_no_of_columns()!= n)
		B.assign((*this * order(O)).base);

	for (i=0;i<n;tmp[i++] = 0){
		tmp[i]=1;               //tmp=e_i;
		// compute the images of the basis of O/pO
		a= alg_number(tmp,1,O);
		//a=e_i
		for (j=0; j < n; j++){
			// compute the images of the basis of A+pO
			// under the MULT-with-a-homomorphism
			b=alg_number(tmp2 = B.get_column(j), 1, O);
			delete[] tmp2;
			y=a*b;   // the image under the MULT-with-a-homomorphism
			init.sto_column_vector(y.coeff_vector(),n,i*n+j);
		}
	}
	delete[] tmp;

	// the image must be written relative to base;
	debug_handler_c("module", "in member - function ring_of_multipliers(...)",
			1, std::cout <<"solve "<<B<<init<<std::endl<<std::flush );
	v.reginvimage(B, init);
	debug_handler_c("module", "in member - function ring_of_multipliers(...)",
			1, std::cout <<"solution is "<<v<<std::endl<<std::flush);
	// move the result to the Matrix VW
	for (i=0;i<n;i++){
		for (j=0; j < n; j++){
			for (k=0;k<n;k++){
				VW.sto(k+j*n,i,v.member(k,i*n+j));
			}
		}
	}
	debug_handler_c("module", "in member - function ring_of_multipliers(...)",
			1, std::cout << "Compute kernel of"<<VW);
	C.kernel(VW);
	debug_handler_c("module", "in member - function ring_of_multipliers(...)",
			1, std::cout << "Kernel is " << C << std::flush);

	// interpret C as a module M and lift it into the order!!
	init.set_no_of_columns(k=C.get_no_of_columns());
	base_vector <bigint> new_tmp;
	for (j = 0; j < n; j++){
		C.row_vector(new_tmp,j);
		init.sto_row_vector(new_tmp,k,j);
	}
	module M(init,1,O);
	M+= p * order(O);  // hopefully interpreted as module-mult.
	debug_handler_c("module","in member - function ring_of_multipliers(...)",
			4,std::cout <<"Module is"<<M<<std::flush);

	//  Instead of creating the multiplication table, we create the base:
	if (!O->f.is_zero()){
		M.den *= p;
		if (O->base_computed()){
			debug_handler_c("module",
					"in member - function ring_of_multipliers(...)",
					4,std::cout <<"Multiply module by "<<O->base<<std::flush);
			multiply(M.base, O->base, M.base);
			M.den *= O->den;
		}
		M.normalize();
		debug_handler_c("module","ring_of_multipliers(...)",4,
				std::cout << "Transformation is  1/" << M.den << " * ";
				std::cout << M.base << std::endl << std::flush);

		order result(M.base, M.den, O);
		debug_handler_c("module","ring_of_multipliers(...)",4,
				std::cout << "So the new order has the following table";
				std::cout << result << std::endl << std::flush;
				std::cout << "and the Discriminant is" <<disc(result);
				std::cout << std::endl << std::flush);
		return result;
	}
	else {
		// If we can't create the base (because there is no field)
		// create the MT instead.
		base_vector <bigint> tmp;
		init.set_no_of_columns((n*(n+1))/2);
		for (i = 0; i < n; i++){
			M.base.get_column_vector(tmp,i);
			a=alg_number(tmp,1,O);
			for (j = 0; j <=i; j++){
				M.base.get_column_vector(tmp,j);
				b=alg_number(tmp,1,O);
				init.sto_column_vector((a*b).coeff_vector(), n, (i*(i+1))/2+j);
			}
		}
		debug_handler_c("module",
				"in member - function ring_of_multipliers(...)",3,
				std::cout <<"solve "<<M.base<<init<<std::endl<<std::flush );

		init.reginvimage(M.base, init);
		init.set_no_of_rows(n);

		debug_handler_c("module",
				"in member - function ring_of_multipliers(...)",3,
				std::cout << "p * MT is" << trans(init) << std::flush);

		for (i = 0; i < init.get_no_of_rows(); i++)
			for (j = 0; j < init.get_no_of_columns(); j++){
				bigint q,r;
				div_rem(q,r,init.member(i,j),p);
				if (!r.is_zero())
					error_handler("module", "ring_of multipliers::internal error::"
						      "division by p unsuccesful");
				init.sto(i,j,q);
			}
		debug_handler_c("module",,
				"in member - function ring_of_multipliers(...)",3,
				std::cout << "MT is" << trans(init) << std::flush);
		return order(trans(init));
	}
}
#endif



// Reduction:
void alg_ideal::reduce(alg_number & divisor)
{
	lidia_size_t i, n = degree();
	bigfloat_lattice lat(n, n);
	lat.set_basis_flag();

	bigmod_matrix tmp(base);
	tmp.lift(0);

	int prec = (tmp.max_abs()).length();
	long store = bigfloat::get_precision();
	bigfloat::set_precision(static_cast<long>((prec-3)*L2B10*bits_per_base_digit+4));
	for (i = 0; i < n; i++)
		for (register lidia_size_t j = 0; j < n; j++)
			lat.sto(i, j, bigfloat(tmp.member(i, j)));

	//   std::cout <<"Zeros are " << ((nf_base *)O)->get_conjugates()<<std::flush;
	//   std::cout << "Trafo is "<< lat<<std::flush;

	multiply(lat, (static_cast<nf_base *>(O))->get_conjugates(), lat);

	bigfloat multiplier (std::sqrt(2.0));
	bigfloat_matrix mult(degree(), degree());
	for (i = 0; i < O->real_roots; i++)
		mult.sto(i, i, 1);
	for (; i < n; i++)
		mult.sto(i, i, multiplier);
	multiply(lat, mult, lat);

	//   std::cout << "Lattice to reduce is "<<lat<<std::flush;
	//   math_vector<bigfloat> a,b;
	//   lat.row_vector(a, 0);
	//   lat.row_vector(b, 1);
	//   std::cout << "Difference is "<<a-b<<std::flush;
	bigint_matrix TRAFO(n, n);
	lat.lll(TRAFO, 1, 4);
	//   std::cout << "Reduced lattice is "<<lat<<std::flush;
	// Search the smallest column of lat:
	// ==================================
	// MM, not used: lidia_size_t min_col;
	bigfloat minimum, tmp1, tmp2;
	// MM, not used: bigfloat * lat_col;

	// compute the alg_number corresponding to column 0:
	bigint * tmp_trafo, * prod = new bigint[n];
	tmp_trafo = TRAFO.get_column(0);
	multiply(prod, tmp, tmp_trafo);
	delete[] tmp_trafo;
	{
		alg_number tmp_alg_number(prod, 1, O);
		swap(divisor, tmp_alg_number);
	}
	delete[] prod;

	//  std::cout << "Reducing by " <<divisor<<std::flush;
	// divide by this number and guarantee an integral ideal:
	divide(*this, *this, divisor);
	//  std::cout << " leads to denominator = "<<den<<std::endl<<std::flush;
	divide(divisor, divisor, den);
	den = 1;
	bigfloat::set_precision(store);
}



// Other functions:
alg_ideal inverse(const alg_ideal & a)
{
	alg_ideal c = a;
	c.invert();
	return c;
}



// random (???) ideales
void alg_ideal::randomize(const bigint & b)
{
	debug_handler("alg_ideal", "in member-function "
		      "randomize(const bigint &)");
	den.assign_one();
	base.set_no_of_columns(1);
	base.set_modulus(LiDIA::randomize(b)+1);
	for (register lidia_size_t j = 0; j < degree(); j++)
		base.sto(j, 0, LiDIA::randomize(b));
	is_exp = false;
	multiply(*this, *(static_cast<module *>(this)), alg_ideal_cast(order(O)));
	den = LiDIA::randomize(b)+1;
	normalize();
}



#ifdef LIDIA_NAMESPACE
}	// end of namespace LiDIA
#endif
