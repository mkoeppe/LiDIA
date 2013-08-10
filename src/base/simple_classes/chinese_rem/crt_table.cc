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
//	Author	: Frank J. Lehmann (FL), Thomas Pfahler (TPF)
//	Changes	: See CVS log
//
//==============================================================================================


#ifdef HAVE_CONFIG_H
# include	"config.h"
#endif
#include	"LiDIA/bigint.h"
#include	"LiDIA/bigmod.h"
#include	"LiDIA/base_vector.h"
#include	"LiDIA/crt_table.h"



#ifdef LIDIA_NAMESPACE
namespace LiDIA {
#endif



//***************************************************************
//
//				  class crt_table
//
//***************************************************************


// FIXME: not portable!
#define SDIGIT_MAX 0x7fffffff


const char crt_table::INTEGER = 1;
const char crt_table::MODULAR = 2;


sdigit
crt_table::next_prime (const sdigit & udq)
{
	sdigit rc;
	bigint q = udq;

	dec(q);

	if (q.is_even())
		dec(q);

	while (q > 0) {
		if (is_prime(q, 8)) {
			q.longify (rc);
			return rc;
		}
		q -= 2;
	}

	lidia_error_handler("crt_table", "next_prime(sdigit)::out of prime numbers");

	return static_cast<sdigit>(0);
}



void
crt_table::gen_primes (bigint B , const sdigit & bound)
{
	debug_handler("crt_table", "gen_primes (const bigint& , const sdigit&)");

	lidia_size_t    n;
	sdigit help;
	base_vector< sdigit > temp (100, vector_flags(vector_flags::expand));

	B.multiply_by_2();

	M.assign_one();
	n = 0;
	help = bound;

	while (B >= M) {
		help = NEXT_PRIME (help);

		if (help < 0)
			help = -help;

		temp[n] = help;
		multiply (M, M, help);
		n++;
	}

	if (n == 0)
		lidia_error_handler("crt_table", "crt_table (bigint , const sdigit &):: upper bound too low");

	temp.set_size (n);
	swap (primes, temp);
	num_primes = n;

	add (M_over_2, M, 1UL);
	shift_right (M_over_2, M_over_2, 1);
}



void
crt_table::set_primes (const base_vector< sdigit > & Prime_vec)
{
	debug_handler("crt_table", "set_primes (const base_vector< sdigit > &)");

	lidia_size_t i;
	const base_vector< sdigit > & Primes = Prime_vec;

	M.assign_one();

	primes = Primes;
	num_primes = primes.size();
	primes.set_capacity (num_primes);


	for (i = 0; i < num_primes; i++) {
		primes[i] = (Primes[i] > 0) ? (Primes[i]) : ((-1) * Primes[i]);
		multiply (M, M, primes[i]);
	}

	add (M_over_2, M, 1UL);
	shift_right (M_over_2, M_over_2, 1);
}



void
crt_table::build_tables()
{
	debug_handler("crt_table", "build_tables ()");

	sdigit tu;

	bigint q, t, t1;
	bigint M1, B, gcdrc;
	bigint minus_M;

	lidia_size_t  i;


	coeff_mod_p.set_capacity(num_primes);
	x.set_capacity (num_primes);

	negate (minus_M, M);

	if (mode == MODULAR) //case bigmod: minus_M_mod_p = -M % p
	{
		remainder (minus_M_mod_p, minus_M, p);
	}
	else {
		minus_M_mod_p = minus_M;
	}

	for (i = 0; i < num_primes; i++) {
		q = bigint(primes[i]);

		divide (M1, M, q);
		remainder (t1, M1, q);
		gcdrc = xgcd_left (t, t1, q);
		if (t.is_negative()) t += q;

		if (! gcdrc.is_one())
			lidia_error_handler ("crt_table", "build_tables()::moduli with non-trivial gcd");

		multiply (M1, M1, t);

		if (mode == MODULAR)
			remainder(coeff_mod_p[i] , M1, p);
		else
			coeff_mod_p[i] = M1;

		t.longify (tu);
		x[i] = static_cast<double>(tu)/static_cast<double>(primes[i]);
	}
}



// * * * Initialization functions  * * *

void
crt_table::init (const bigint &B, const bigint & P)
{
	debug_handler("crt_table", "init (const bigint &, const bigint)");

	if (reference_counter > 0)
		lidia_error_handler ("crt_table", "init(const bigint &, const bigint)::table still in use");

	if (! B.is_positive())
		lidia_error_handler ("crt_table", "init(const bigint &, const bigint)::non-positive bound for absolute values");

	if (! P.is_zero())
		mode = MODULAR;
	else
		mode = INTEGER;

	p = P;

	gen_primes (B , SDIGIT_MAX);
	build_tables ();
}



void
crt_table::init (const bigint &B, const sdigit & prime_bound, const bigint & P)
{
	debug_handler("crt_table", "init (const bigint &, const sdigit &, const bigint)");

	if (reference_counter > 0)
		lidia_error_handler ("crt_table", "init(const bigint &, const sdigit &, const bigint)::table still in use");

	if (! B.is_positive())
		lidia_error_handler ("crt_table", "init(const bigint &, const sdigit &, const bigint)::non-positive bound for absolute values");

	if (prime_bound < 0)
		lidia_error_handler ("crt_table", "init(const bigint &, const sdigit &, const bigint)::negative bound for positive primes");

	if (! P.is_zero())
		mode = MODULAR;
	else
		mode = INTEGER;

	p = P;

	gen_primes (B , prime_bound);
	build_tables ();
}



void
crt_table::init (const base_vector< sdigit > & Primes, const bigint & P)
{
	debug_handler("crt_table", "init (const base_vector< sdigit > & , const bigint)");

	if (reference_counter > 0)
		lidia_error_handler ("crt_table", "init(const base_vector< sdigit > & , const bigint)::table still in use");

	if (! P.is_zero())
		mode = MODULAR;
	else
		mode = INTEGER;

	p = P;

	set_primes (Primes);
	build_tables ();
}



void
crt_table::clear ()
{
	debug_handler("crt_table", "clear ()");

	if (reference_counter > 0)
		lidia_error_handler ("crt_table", "clear()::table still in use");

	num_primes = 0;

	primes.kill();
	coeff_mod_p.kill();
	x.kill();

	mode = INTEGER;

	p = 0;
	minus_M_mod_p = M = 1;

	NEXT_PRIME = &next_prime;
}



// * * *  constructors  * * *

crt_table::crt_table ()
{
	debug_handler("crt_table", "crt_table ()");

	reference_counter = 0;
	num_primes = 0;

	NEXT_PRIME = &next_prime;
}



crt_table::crt_table (const bigint &B, const bigint & P)
{
	debug_handler("crt_table", "crt_table (const bigint &, const bigint&)");

	reference_counter = 0;

	NEXT_PRIME = &next_prime;
	init (B, P);
}



crt_table::crt_table (const bigint &B)
{
	debug_handler("crt_table", "crt_table (const bigint &)");

	reference_counter = 0;

	NEXT_PRIME = &next_prime;
	init (B);
}



crt_table::crt_table (const bigint &B, const sdigit & prime_bound, const bigint & P)
{
	debug_handler("crt_table", "crt_table (const bigint &, const sdigit &, const bigint&)");

	reference_counter = 0;

	NEXT_PRIME = &next_prime;
	init (B, prime_bound, P);
}



crt_table::crt_table (const bigint &B, const sdigit & prime_bound)
{
	debug_handler("crt_table", "crt_table (const bigint &, const sdigit &)");

	reference_counter = 0;

	NEXT_PRIME = &next_prime;
	init (B, prime_bound);
}



crt_table::crt_table (const base_vector< sdigit > & Primes, const bigint & P)
{
	debug_handler("crt_table", "crt_table (const base_vector< sdigit > &, const bigint&)");

	reference_counter = 0;

	NEXT_PRIME = &next_prime;
	init (Primes, P);
}



crt_table::crt_table (const base_vector< sdigit > & Primes)
{
	debug_handler("crt_table", "crt_table (const base_vector< sdigit > &)");

	reference_counter = 0;

	NEXT_PRIME = &next_prime;
	init (Primes);
}



crt_table::~crt_table()
{
	debug_handler("crt_table", "~crt_table()");
}



sdigit
crt_table::get_prime (lidia_size_t index) const
{
	debug_handler("crt_table", "get_prime(lidia_size_t)");

	if (index< 0 || index >= num_primes)
		lidia_error_handler("crt_table", "get_prime(lidia_size_t)::index out of range");

	return (primes[index]);
}



lidia_size_t
crt_table::how_much_to_use (bigint  B) const
{
	debug_handler("crt_table", "how_much_to_use (bigint)");

	lidia_size_t    rc;
	bigint N;

	if (! B.is_positive())
		lidia_error_handler("crt_table", "how_much_to_use(bigint)::non-positive bound for absolute values");

	B.multiply_by_2();

	if (num_primes == 0 || B >= M) {
		rc = -1;
	}
	else if (mode == MODULAR) {
		// you have to use every prime in MODULAR mode

		rc = num_primes;
	}
	else {
		N.assign_one();
		rc = 0;

		while (N <= B) {
			multiply (N , N , primes[rc]);
			rc++;
		}
	}

	return rc;
}



void
crt_table::set_prime_generator (sdigit (*np) (const sdigit & q))
{
	debug_handler("crt_table", "set_prime_generator(sdigit (*np) (const sdigit & q))");

	if (np != NULL)
		NEXT_PRIME = np;
	else
		NEXT_PRIME = &next_prime;
}



void
crt_table::info () const
{
	std::cout << "==================crt_table==================\n" << std::endl;

	std::cout << "     Mode        :  " << ((mode == INTEGER) ? "INTEGER" : "MODULAR") << std::endl;
	std::cout << "\n";
	std::cout << "     References  :  " << reference_counter << std::endl;
	std::cout << "\n";
	std::cout << "     Primes      :  ";
	if (num_primes < 20)
		std::cout << primes << std::endl;
	else
		std::cout << num_primes << std::endl;
	std::cout << "\n";
	if (mode == MODULAR)
		std::cout << "     Modul       :  " << p << std::endl;
	std::cout << "     Produkt     :  " << M << std::endl;
	std::cout << "     Negative    :  " << minus_M_mod_p << std::endl;
	std::cout << "\n";
	std::cout << "     Current Sto.:  \n" << std::endl;
	std::cout << "           " << coeff_mod_p << std::endl;
	std::cout << "           " << x << std::endl;
	std::cout << "\n";
	std::cout << "==================*********==================\n" << std::endl;
}



#ifdef LIDIA_NAMESPACE
}	// end of namespace LiDIA
#endif
