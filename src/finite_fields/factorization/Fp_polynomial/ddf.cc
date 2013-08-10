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
//	Author	: Victor Shoup (VS), Thomas Pfahler (TPf)
//	Changes	: See CVS log
//
//==============================================================================================


#ifdef HAVE_CONFIG_H
# include	"config.h"
#endif
#include	"LiDIA/Fp_polynomial.h"
#include	"LiDIA/finite_fields/Fp_polynomial_util.h"
#include	"LiDIA/Fp_poly_modulus.h"
#include	"LiDIA/lidia_signal.h"
#include	"LiDIA/lidia_file.h"
#include	"LiDIA/single_factor.h"
#include	"LiDIA/factorization.h"

#include	<fstream>
#include	<cstdio>
#include	<unistd.h>



#ifdef LIDIA_NAMESPACE
namespace LiDIA {
#endif



int DDF_GCD_BLOCKING_FACTOR = 4;

#define BABY 0
#define GIANT 1


//
// Temporary files will be located in the current directory or in /tmp
// and will have the following names:
// 	ddf1234b.001		(baby steps)
//	ddf1234b.002
//	...
//	ddf1234g.001		(giant steps)
//	ddf1234g.002
//	...
// where 1234 is the pid of the current process (modulo 10^5).
//


//
// Task:        tries to create the file "ddf....b.001", where "...." is
//		replaced by the pid of the current process (mod 10^5),
//		 a) in the current directory
//		 b) if that fails, in the directory /tmp
//		an error will occur, if the directories are not writeable or
//		 if a file with the same name already exists in both dirs.
//		finally, stem="[successful dir]/ddf[pid]"
//
// Conditions:	enough space for stem must be allocated
//

static
void check_stem(char* stem)
{
	memory_handler(tmp, "Fp_polynomial",
		       "create_stem(void)::Error in memory allocation");

	int num = static_cast<int>(getpid()) % 10000; // 4 digits are enough
	bool ok = false;

	sprintf(stem, "./ddf%db.001", num); // try current directory...

	// <MM>
	// Replaced by use of lidia_file,
	// s.is_open() is not portable
	//
	//s.open(stem, std::ios::out|std::ios::noreplace); // test if pid is already used
	//if (s.is_open()) {	// PT
	//	sprintf(stem, "./ddf%d", num);
	//	ok = true; // file was successfully created
	//}
	//s.close();

	if (!lidia_file::file_exists(stem)) {
		sprintf(stem, "./ddf%d", num);
		ok = true; // filename was successfully created
	}
	// </MM>


	if (!ok) {
		// try "/tmp"-directory...
		sprintf(stem, "/tmp/ddf%db.001", num);

		//<MM>
		//
		//s.open(stem, std::ios::out|std::ios::noreplace);
		//if (s.is_open()) {	// PT
		//	sprintf(stem, "/tmp/ddf%d", num);
		//	ok = true;
		//}
		//s.close();

		if (!lidia_file::file_exists(stem)) {
			sprintf(stem, "./ddf%d", num);
			ok = true; // filename was successfully created
		}
		//</MM>
	}

	if (!ok)
		lidia_error_handler("Fp_polynomial", "create_stem(...)::cannot write "
				    "temporary files into current nor into /tmp-directory");
}



//
// Task:	returns the string "[stem][ext].[d]", where d is filled up
//		with zeros if necessary to get a 3-digit number
//
// Conditions:	the whole string must have less than 32 characters
//

static
char * create_file_name(const char *stem, int ext, lidia_size_t d)
{
	debug_handler("Fp_polynomial", "create_file_name (char*, char*, lidia_size_t)");

	static char sbuf[32];

	strcpy(sbuf, stem);
	switch(ext) {
	case(BABY):  strcat(sbuf, "b.");
		break;
	case(GIANT): strcat(sbuf, "g.");
		break;
	default: lidia_error_handler("Fp_polynomial", "create_file_name(...)"
				     "::wrong ext");
	}

	char dbuf[4];
	dbuf[3] = '\0';
	lidia_size_t i, dd;
	dd = d;

	for (i = 2; i >= 0; i--) {
		dbuf[i] = static_cast<char>((dd % 10) + '0');
		dd = dd / 10;
	}

	strcat(sbuf, dbuf);
	return sbuf;
}



//
// Task:	delete all temporary files
//

static
void file_cleanup(lidia_size_t k, lidia_size_t l, const char *new_ddf_stem)
{
	debug_handler("Fp_polynomial", "file_cleanup(lidia_size_t, lidia_size_t, char*)");
	lidia_size_t i;

	for (i = 1; i <= k; i++)
		std::remove(create_file_name(new_ddf_stem, BABY, i));

	for (i = 1; i <= l; i++)
		std::remove(create_file_name(new_ddf_stem, GIANT, i));
}




//
// Task :	remove all temporary files in case of an interruption
//

LIDIA_SIGNAL_FUNCTION(tidy_up)
{
	int num = static_cast<int>(getpid()) % 10000;
	char stem[32];

	// delete all temporary files (if any) in directory "./"
	sprintf(stem, "./ddf%d", num);
	int i = 1;
	while (! std::remove(create_file_name(stem, BABY, i)))
		i++;
	i = 1;
	while (! std::remove(create_file_name(stem, GIANT, i)))
		i++;

	// delete all temporary files (if any) in directory "/tmp/"
	sprintf(stem, "/tmp/ddf%d", num);
	i = 1;
	while (! std::remove(create_file_name(stem, BABY, i)))
		i++;
	i = 1;
	while (! std::remove(create_file_name(stem, GIANT, i)))
		i++;

	lidia_error_handler("Fp_polynomial", "ddf : the program has been "
			    "interrupted");
}



//
// Task:	reads giant step #gs into g
//

static
void fetch_giant_step(Fp_polynomial & g, lidia_size_t gs,
		      const Fp_poly_modulus & F, const char *new_ddf_stem)
{
	debug_handler("Fp_polynomial", "fetch_giant_step(Fp_polynomial&, lidia_size_t, Fp_poly_modulus&, char*)");

	std::ifstream s;

	s.open(create_file_name(new_ddf_stem, GIANT, gs), std::ios::in);
	if (!s) {
		lidia_error_handler_c("Fp_polynomial",
				      "fetch_giant_step(...)::open read error",
				      std::cout << "open read error: " << create_file_name(new_ddf_stem,
											   GIANT, gs) << "\n";);
		return; // LC
	}

	s >> g;
	s.close();
	remainder(g, g, F);
}



//
// Task:	reads baby steps #1 to #k into v
//

static
void fetch_baby_steps(base_vector< Fp_polynomial > &v, lidia_size_t k,
		      const char *new_ddf_stem, const bigint & p)
{
	debug_handler("Fp_polynomial", "fetch_baby_steps(base_vector< Fp_polynomial > &, lidia_size_t, char*, bigint&)");

	std::ifstream s;

	if (v.capacity() < k)
		v.set_capacity(k);
	v.set_size(k);


	lidia_size_t i;
	for (i = 1; i <= k - 1; i++) {
		s.open(create_file_name(new_ddf_stem, BABY, i), std::ios::in);
		if (!s) {
			lidia_error_handler_c("Fp_polynomial",
					      "fetch_baby_steps(...):::open read error",
					      std::cout << "open read error: "
					      << create_file_name(new_ddf_stem, BABY, i) << "\n";);
			return; // LC
		}

		s >> v[i];
		s.close();
	}

	v[0].set_modulus(p);
	v[0].assign_x();
}



//
// Task:	writes baby steps to disk
//

static
void generate_baby_steps(Fp_polynomial & h1, const Fp_polynomial & f,
			 const Fp_polynomial & h, lidia_size_t k, const char *new_ddf_stem)
{
	debug_handler("Fp_polynomial", "generate_baby_steps(Fp_polynomial&, Fp_polynomial&, Fp_polynomial&, lidia_size_t, char*)");

	f.comp_modulus(h, "generate_baby_steps");

	my_timer t;
	std::ofstream s;
	bool verbose = single_factor< Fp_polynomial >::verbose();

	if (verbose) t.start("generating baby steps...");

	Fp_poly_modulus F(f);

	poly_argument H;
	H.build(h, F, 2 * square_root(f.degree()));


	h1.assign(h);

	lidia_size_t i;

	for (i = 1; i <= k - 1; i++) {
		s.open(create_file_name(new_ddf_stem, BABY, i), std::ios::out);
		if (!s) {
			lidia_error_handler_c("Fp_polynomial",
					      "generate_baby_steps(...)::open write error",
					      std::cout << "open write error: "
					      << create_file_name(new_ddf_stem, BABY, i) << "\n";);
			return; // LC
		}

		s << h1 << "\n";
		s.close();

		H.compose(h1, h1, F);
		if (verbose) std::cerr << "+";
	}

	if (verbose) t.stop();
}


//
// Task:	writes giant steps to disk
//

static
void generate_giant_steps(const Fp_polynomial & f, const Fp_polynomial & h,
			  lidia_size_t l, const char *new_ddf_stem)
{
	debug_handler("Fp_polynomial", "generate_giant_steps(Fp_polynomial&, Fp_polynomial&, lidia_size_t, char*)");

	f.comp_modulus(h, "generate_giant_steps");

	my_timer t;
	std::ofstream s;
	bool verbose = single_factor< Fp_polynomial >::verbose();

	if (verbose) t.start("generating giant steps...");

	Fp_poly_modulus F(f);

	poly_argument H;
	H.build(h, F, 2 * square_root(f.degree()));

	Fp_polynomial h1;

	h1.assign(h);

	lidia_size_t i;

	for (i = 1; i <= l; i++) {
		s.open(create_file_name(new_ddf_stem, GIANT, i), std::ios::out);
		if (!s) {
			lidia_error_handler_c("Fp_polynomial",
					      "generate_giant_steps(...)::open write error",
					      std::cout << "open write error: "
					      << create_file_name(new_ddf_stem, GIANT, i) << "\n";);
			return; // LC
		}

		s << h1 << "\n";
		s.close();

		if (i != l)
			H.compose(h1, h1, F);

		if (verbose) std::cerr << "+";
	}

	if (verbose) t.stop();
}



//
// Task:	appends g, a product of irreducible polynomials of degree m,
//		to u
//		(the exponents of u are the degrees of the irreducible
//		polynomials !)
//

static
void new_add_factor(factorization< Fp_polynomial > &u,
		    const Fp_polynomial & g, lidia_size_t m)
{
	debug_handler("Fp_polynomial", "new_add_factor(factorization< Fp_polynomial > &, Fp_polynomial&, lidia_size_t)");

	u.append(g, m);

	if (single_factor< Fp_polynomial >::verbose())
		std::cerr << "split " << m << " " << g.degree() << std::endl;
}



//
// Task:	splits f by computing gcd(f, buf[i])
//
// Algorithm:	instead of computing gcd(f,buf[i]) for all i, we compute
//		g = product of all buf[i]. If gcd(f,g) = 1, we do not have to
//		compute any further gcd and can return immediately
//

static
void new_process_table(factorization< Fp_polynomial > &u, Fp_polynomial & f,
		       const Fp_poly_modulus & F, base_vector< Fp_polynomial > &buf,
		       lidia_size_t size, lidia_size_t StartInterval,
		       lidia_size_t IntervalLength)
{
	debug_handler("Fp_polynomial", "new_process_table(factorization< Fp_polynomial > &, Fp_polynomial&, Fp_poly_modulus&, base_vector< Fp_polynomial > &, lidia_size_t, lidia_size_t, lidia_size_t)");

	if (size == 0)
		return;

	Fp_polynomial & g = buf[size - 1];

	lidia_size_t i;

	for (i = 0; i < size - 1; i++)
		multiply(g, g, buf[i], F);

	gcd(g, f, g);
	if (g.degree() == 0)
		return;

	divide(f, f, g);

	lidia_size_t d = (StartInterval - 1) * IntervalLength + 1;
	i = 0;
	lidia_size_t interval = StartInterval;

	//next, we 'refine' our gcd computations
	while (i < size - 1 && 2 * d <= g.degree()) {
		gcd(buf[i], buf[i], g);
		if (buf[i].degree() > 0) {
			new_add_factor(u, buf[i], interval);
			divide(g, g, buf[i]);
		}
		i++;
		interval++;
		d += IntervalLength;
	}

	if (g.degree() > 0) {
		if (i == size - 1)
			new_add_factor(u, g, interval);
		else
			new_add_factor(u, g, (g.degree()+IntervalLength-1)/IntervalLength);
	}
}



static
void giant_refine(factorization< Fp_polynomial > &u, const Fp_polynomial & ff,
		  lidia_size_t k, const char *new_ddf_stem)
{
	debug_handler("Fp_polynomial", "giant_refine(factorization< Fp_polynomial > &, Fp_polynomial&, lidia_size_t, char*)");

	my_timer t;
	bool verbose = single_factor< Fp_polynomial >::verbose();

	if (verbose) {
		std::cerr << "giant refine...\n";
		t.start("giant refine time: ");
	}

	u.kill();

	base_vector< Fp_polynomial > BabyStep;

	fetch_baby_steps(BabyStep, k, new_ddf_stem, ff.modulus());

	base_vector< Fp_polynomial >
		buf(static_cast<lidia_size_t>(DDF_GCD_BLOCKING_FACTOR),
		    static_cast<lidia_size_t>(DDF_GCD_BLOCKING_FACTOR));

	Fp_polynomial f;
	f.assign(ff);

	Fp_poly_modulus F;
	F.build(f);

	Fp_polynomial g;
	Fp_polynomial h;

	lidia_size_t size = 0;
	lidia_size_t first_gs = 0; //initialized only because of compiler warnings

	lidia_size_t d = 1;

	while (2 * d <= f.degree()) {

		lidia_size_t old_n = f.degree();

		lidia_size_t gs = (d + k - 1) / k;
		lidia_size_t bs = gs * k - d;

		if (bs == k - 1) {
			size++;
			if (size == 1)
				first_gs = gs;
			fetch_giant_step(g, gs, F, new_ddf_stem);
			subtract(buf[size - 1], g, BabyStep[bs]);
		}
		else {
			subtract(h, g, BabyStep[bs]);
			multiply(buf[size - 1], buf[size - 1], h, F); //Fp_poly_modulus

		}

		if (verbose && bs == 0)
			std::cerr << "+";

		if (size == DDF_GCD_BLOCKING_FACTOR && bs == 0) {
			new_process_table(u, f, F, buf, size, first_gs, k);
			if (verbose) std::cerr << "*";
			size = 0;
		}

		d++;

		if (2 * d <= f.degree() && f.degree() < old_n) {
			F.build(f);

			lidia_size_t i;
			for (i = 1; i <= k - 1; i++)
				remainder(BabyStep[i], BabyStep[i], F);
		}
	}

	if (size > 0) {
		new_process_table(u, f, F, buf, size, first_gs, k);
		if (verbose) std::cerr << "*";
	}

	if (f.degree() > 0)
		new_add_factor(u, f, -1);

	if (verbose) t.stop();
}



static
void interval_refine(factorization< Fp_polynomial > &factors,
		     const Fp_polynomial &ff, lidia_size_t k, lidia_size_t gs,
		     const base_vector< Fp_polynomial > &BabyStep, const char *new_ddf_stem)
{
	debug_handler("Fp_polynomial", "interval_refine(factorization< Fp_polynomial > &, Fp_polynomial&, lidia_size_t, lidia_size_t, base_vector< Fp_polynomial > &, char*)");

	if (BabyStep.size() != 0)
		ff.comp_modulus(BabyStep[0], "interval_refine");

	base_vector< Fp_polynomial >
		buf(static_cast<lidia_size_t>(DDF_GCD_BLOCKING_FACTOR),
		    static_cast<lidia_size_t>(DDF_GCD_BLOCKING_FACTOR));

	Fp_polynomial f(ff);

	Fp_poly_modulus F;
	F.build(f);

	Fp_polynomial g;

	fetch_giant_step(g, gs, F, new_ddf_stem);

	lidia_size_t size = 0;

	lidia_size_t first_d = 0; //initialized only because of compiler warnings

	lidia_size_t d = (gs - 1) * k + 1;
	lidia_size_t bs = k - 1;

	while (2 * d <= f.degree()) {
		lidia_size_t old_n = f.degree();

		if (size == 0)
			first_d = d;
		remainder(buf[size], BabyStep[bs], F);
		subtract(buf[size], buf[size], g);
		size++;

		if (size == DDF_GCD_BLOCKING_FACTOR) {
			new_process_table(factors, f, F, buf, size, first_d, 1);
			size = 0;
		}

		d++;
		bs--;

		if (bs < 0) {
			//exit
			d = f.degree() + 1;
		}

		if (2 * d <= f.degree() && f.degree() < old_n) {
			F.build(f);
			remainder(g, g, F);
		}
	}

	new_process_table(factors, f, F, buf, size, first_d, 1);

	if (f.degree() > 0)
		new_add_factor(factors, f, f.degree());
}



static
void baby_refine(factorization< Fp_polynomial > &factors,
		 const factorization< Fp_polynomial > &u,
		 lidia_size_t k, const char *new_ddf_stem)
{
	debug_handler("Fp_polynomial", "baby_refine(factorization< Fp_polynomial > &, factorization< Fp_polynomial > &, lidia_size_t, char*)");

	my_timer t;
	bool verbose = single_factor< Fp_polynomial >::verbose();

	if (verbose) {
		std::cerr << "baby refine...\n";
		t.start("baby refine time: ");
	}

	factors.kill();

	base_vector< Fp_polynomial > BabyStep;

	lidia_size_t i;

	for (i = 0; i < u.no_of_composite_components(); i++) {
		const Fp_polynomial & g = u.composite_base(i).base();
		lidia_size_t gs = u.composite_exponent(i);

		if (gs == -1 || 2 * ((gs - 1) * k + 1) > g.degree())
			new_add_factor(factors, g, g.degree());
		else {
			if (BabyStep.size() == 0)
				fetch_baby_steps(BabyStep, k, new_ddf_stem,
						 u.composite_base(0).base().modulus());
			interval_refine(factors, g, k, gs, BabyStep, new_ddf_stem);
		}
	}

	if (verbose) t.stop();
}



//*************************************************************************
//
//			the main routine
//
//*************************************************************************

//
// Literature:	V. Shoup
//		A New Polynomial Factorization Algorithm and its Implementation
//		Universitaet des Saarlandes, 1994
//		Shoup, J. Symbolic Comp. 20:363-397, 1995
//
// Task:        Performs the distinct degree factorization of f = F.modulus()
//		(i.e. splits f in products of irreducibles of the same
//		degree)
//		when finished, the exponents of 'factors' are the degrees of
//		the irreducible factors of 'f' !!!
//
// Conditions:  f is square-free, monic, deg(f) > 0 and f(0) != 0 (?).
//              h = x^p mod f
//
// Algorithm:   DDF (Shoup)
//		some intermediate results are written to disk, they are
//		deleted at the end of the routine
//		baby step/giant step approach
//

void
ddf(factorization< Fp_polynomial > &factors,
    const Fp_polynomial & f, const Fp_polynomial & h)
{
	debug_handler("Fp_polynomial", "ddf(factorization< Fp_polynomial > &, Fp_polynomial&, Fp_polynomial&)");

	f.comp_modulus(h, "ddf");

	if (f.degree() <= 0) {
		lidia_error_handler("Fp_polynomial",
				    "ddf(...)::polynomial has degree <= 0");
		return; // LC
	}

	if (f.degree() == 1) {
		single_factor< Fp_polynomial > tmp(f);
		factors.assign(tmp);
		return;
	}

	char *new_ddf_stem = new char[32];
	memory_handler(new_ddf_stem, "Fp_polynomial",
		       "ddf::Error in memory allocation");
	check_stem(new_ddf_stem);

	lidia_signal sig1 (SIGTERM, tidy_up);
	lidia_signal sig2 (SIGINT, tidy_up);
	//PT lidia_signal sig3 (SIGHUP, tidy_up);

	lidia_size_t B = f.degree() / 2;
	lidia_size_t k = square_root(B);
	lidia_size_t l = (B + k - 1) / k;

	Fp_polynomial h1;
	generate_baby_steps(h1, f, h, k, new_ddf_stem);
	generate_giant_steps(f, h1, l, new_ddf_stem);

	factorization< Fp_polynomial > u;
	giant_refine(u, f, k, new_ddf_stem);
	baby_refine(factors, u, k, new_ddf_stem);

	file_cleanup(k, l, new_ddf_stem);
	delete[] new_ddf_stem;
}



//*************************************************************************
//
//			    "interface"
//
//*************************************************************************

factorization< Fp_polynomial > ddf(const Fp_polynomial& f)
{
	debug_handler("Fp_polynomial", "ddf(Fp_polynomial&)");

	factorization< Fp_polynomial > factors;
	Fp_poly_modulus F(f);
	Fp_polynomial b;
	power_x(b, f.modulus(), F);
	ddf(factors, f, b);
	return factors;
}



factorization< Fp_polynomial > single_factor< Fp_polynomial >::ddf() const
{
	debug_handler("single_factor< Fp_polynomial >", "ddf()");
	return LiDIA::ddf(rep);
}



#ifdef LIDIA_NAMESPACE
}	// end of namespace LiDIA
#endif
