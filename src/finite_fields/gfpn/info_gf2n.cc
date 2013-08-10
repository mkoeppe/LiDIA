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
//	Author	: Franz-Dieter Berger (FDB), Patric Kirsch (PK),
//                Volker Mueller (VM), Thomas Pfahler (TPf)
//	Changes	: See CVS log
//
//==============================================================================================


#ifdef HAVE_CONFIG_H
# include	"config.h"
#endif
#include	"LiDIA/finite_fields/info_gf2n.h"
#include        "LiDIA/osstream.h"
#include	"LiDIA/path.h"
#include	<cstdlib>
#include	<cstring>
#include	<fstream>
#include        <string>



#ifdef LIDIA_NAMESPACE
namespace LiDIA {
#endif



udigit *info_gf2n::B;
udigit *info_gf2n::C;
udigit *info_gf2n::F;
udigit *info_gf2n::G;
unsigned int info_gf2n::max_anzBI;



//*********************************************************************

static inline unsigned int max(unsigned int a, unsigned int b)
{
	if (a >= b)
		return a;
	else return b;
}



static inline unsigned int min(unsigned int a, unsigned int b)
{
	if (a <= b)
		return a;
	else return b;
}



//*********************************************************************

static
char* read_modulus_from_database (unsigned int in_degree)
{
	std::string name_database = "";
	std::string line;
	unsigned int d, i;
	char ch;
	bool found = false;

	if (getenv("LIDIA_GF2N") != NULL) {
		name_database = getenv("LIDIA_GF2N");
	}
	if(name_database.length() == 0) {
		name_database = LIDIA_GF2N_DB;
	}
	std::ifstream database(name_database.c_str());
	if (database.fail()) {
	        std::string msg= "Can't open file ";
		msg += name_database;
		lidia_error_handler("info_gf2n::read_modulus_from_database",
				    msg.c_str());
		return NULL;
	}

	while (database.good() && !found) {
		database >> d;
		if (database.fail()) {
		    break;
		}

		if (d == in_degree) {
		        osstream oss;
			oss << d;
			i = 1;
			do
			{
				i++;
				database >> d;
				oss << ' ' << d;
			}
			while (d != 0 && !database.fail());
			line = extractString(oss);
			found = !database.fail();
		}
		else {
			do
			{
				database.get(ch);
			}
			while (ch != '\n' && !database.fail());
		}
	}

	if (!found) { // database.fail() == true!
	        std::string msg;
	        if(database.eof()) {
		    msg = "Modulus not found in database ";
		}
		else {
		    msg = "Error while reading database ";
		}
		msg += name_database;
		lidia_error_handler("info_gf2n::read_modulus_from_database",
				    msg.c_str());		
		return NULL;
	}

	database.close();
	char* p = new char[line.length()+1];
	strcpy(p, line.c_str());

	return p;
}



//*******************************************************************
// the two possible initialization routines
//*******************************************************************

info_gf2n::info_gf2n() :
	degree(0), anzBI(0), exponents(0), exp1(0), exp2(0), exp3(0),
	anz_exponents(0), mul_p(0), invert_p(0), partial_reduce1_p(0),
	partial_reduce2_p(0)
{
}



void info_gf2n::init(unsigned int degree)
{
	char* poly = read_modulus_from_database(degree);
	init(poly, degree);
	delete[] poly;
}



void info_gf2n::init(char* poly, unsigned int deg)
{
	char *str = poly;
	unsigned int i, h;
	static bool already_initialized = false;
	bool constant_term = false;

	//------- compute number of exponents and degree  ----------

	i = 1;
	sscanf(str, "%ud", &h);
	degree = h;

	if (h == 0)
		constant_term = true;

	do {
		str = strchr(str, ' ');
		if (str == NULL)
			break;
		else
			str++;
		sscanf(str, "%ud", &h);
		if (h > degree) degree = h;
		if (h == 0)       constant_term = true;
		i++;
	} while (true);

	if (constant_term == false) {
		lidia_error_handler("info_gf2n",
				    "gf2n_init: input polynomial not irreducible");
		return;
	}
	if (deg != degree) {
		lidia_error_handler("info_gf2n",
				    "init: degrees don't match");
		return;
	}

	anzBI = (deg / BITS_PER_LONG) + 1;

	if (anzBI > FIXMUL)
		mul_p = gf2nmul[FIXMUL]; // init. mult. routine
	else
                mul_p = gf2nmul[anzBI];

	anz_exponents = i-2;

	//delete[] exponents;
	exponents = new unsigned int[anz_exponents];

	str = poly;
	for (i = 0; i < anz_exponents; i++) {
		str = strchr(str, ' ');
		str++;
		sscanf(str, "%ud", & h);
		while (h == 0 || h == degree) {
			str = strchr(str, ' ');
			str++;
			sscanf(str, "%ud", & h);
		}
		exponents[i] = h;
	}

	//---- set invsel and initialize variables for inversion and reduction

	if ((anz_exponents == 1) &&
	    ((degree - exponents[0]) >= BITS_PER_LONG) &&
	    (exponents[0] >= BITS_PER_LONG)) {
		invert_p = &info_gf2n::tri_invert;
		partial_reduce1_p = &info_gf2n::tri_partial_reduce1;
		partial_reduce2_p = &info_gf2n::tri_partial_reduce2;
		exp1 = exponents[0];
	}
	else
		if (anz_exponents == 3) {
			unsigned int ma, mi;

			ma = max(exponents[0], max(exponents[1], exponents[2]));
			mi = min(exponents[0], min(exponents[1], exponents[2]));
			if ((degree - ma) >= BITS_PER_LONG && mi >= BITS_PER_LONG) {
				invert_p = &info_gf2n::pent_invert;
				partial_reduce1_p = &info_gf2n::pent_partial_reduce1;
				partial_reduce2_p = &info_gf2n::pent_partial_reduce2;
				exp1 = exponents[0];
				exp2 = exponents[1];
				exp3 = exponents[2];
			}
			else
				invert_p = &info_gf2n::general_invert;
			partial_reduce1_p = &info_gf2n::general_partial_reduce1;
			partial_reduce2_p = &info_gf2n::general_partial_reduce2;
		}
		else {
			invert_p = &info_gf2n::general_invert;
			partial_reduce1_p = &info_gf2n::general_partial_reduce1;
			partial_reduce2_p = &info_gf2n::general_partial_reduce2;
		}


	// ---- INITIIALIZE OPERATIONS  --------------------------

	if (!already_initialized) {
		gen_tables();
		already_initialized = true;
	}

	if (anzBI > info_gf2n::max_anzBI) {
		delete[] info_gf2n::B;
		delete[] info_gf2n::C;
		delete[] info_gf2n::F;
		delete[] info_gf2n::G;
		info_gf2n::B = new udigit [2*anzBI];
		info_gf2n::C = new udigit [2*anzBI];
		info_gf2n::F = new udigit [2*anzBI];
		info_gf2n::G = new udigit [2*anzBI];
		memory_handler(info_gf2n::B, "info_gf2n", "init: out of memory");
		memory_handler(info_gf2n::C, "info_gf2n", "init: out of memory");
		memory_handler(info_gf2n::F, "info_gf2n", "init: out of memory");
		memory_handler(info_gf2n::G, "info_gf2n", "init: out of memory");
		info_gf2n::max_anzBI = anzBI;
	}

	// ---------- Initialisation of  I/O ----------------

#if 0
	gf2nIO::setbase (mode);
	gf2nIO::setprefix (mode);

	std::cout << "anzBI = " << anzBI << std::endl
		  << "degree = " << degree << std::endl;
#endif
}



#ifdef LIDIA_NAMESPACE
}	// end of namespace LiDIA
#endif
