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
//		  Volker Mueller (VM)
//	Changes	: See CVS log
//
//==============================================================================================


#ifdef HAVE_CONFIG_H
# include	"config.h"
#endif
#include	"LiDIA/gf2n.h"
#include	"LiDIA/finite_fields/gf2nIO.h"
#include        "LiDIA/osstream.h"
#include	"LiDIA/path.h"

#include	<fstream>
#include        <string>
#include	<cstdlib>
#include	<cstring>



#ifdef LIDIA_NAMESPACE
namespace LiDIA {
#endif



//*********************************************************************

inline unsigned int max(unsigned int a, unsigned int b)
{
	if (a >= b)
		return a;
	else return b;
}



inline unsigned int min(unsigned int a, unsigned int b)
{
	if (a <= b)
		return a;
	else return b;
}



//*********************************************************************

char* gf2n::read_modulus_from_database (unsigned int in_degree)
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
		lidia_error_handler("gf2n_init::read_modulus_from_database",
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
		lidia_error_handler("gf2n_init::read_modulus_from_database",
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

void gf2n_init (unsigned int degree, gf2nIO::base base)
{
	char* poly;

	if (degree <= 1) {
		lidia_error_handler("gf2n_init", "degree <= 1 not allowed");
		return;
	}

	poly = gf2n::read_modulus_from_database(degree);
	gf2n_init (poly, base);
	delete[] poly;
}



//**********************************************************************

void gf2n_init (char * poly, gf2nIO::base mode)
{
	char *str = poly;
	unsigned int i, h;
	static bool already_initialized = false;
	bool constant_term = false;

	//------- compute number of exponents and degree  ----------

	i = 1;
	sscanf(str, "%ud", &h);
	gf2n::degree = h;

	if (h == 0)
		constant_term = true;

	do {
		str = strchr(str, ' ');
		if (str == NULL)
			break;
		else str++;
		sscanf(str, "%ud", &h);
		if (h > gf2n::degree)
			gf2n::degree = h;
		if (h == 0)
			constant_term = true;
		i++;
	} while (true);

	if (constant_term == false) {
		lidia_error_handler("gf2n", "gf2n_init::input polynomial not irreducible");
		return;
	}

	gf2n::anzBI = (gf2n::degree / (8*sizeof(gf2n_word))) + 1;

	if (gf2n::anzBI > FIXMUL)
		gf2n::mulsel = FIXMUL;
	else
		gf2n::mulsel = gf2n::anzBI;

	gf2n::anz_exponents = i-2;

	delete[] gf2n::exponents;
	gf2n::exponents = new unsigned int[gf2n::anz_exponents];

	gf2n::ord = rational_factorization(1);
	gf2n::ord_is_fact = false;

	str = poly;
	for (i = 0; i < gf2n::anz_exponents; i++) {
		str = strchr(str, ' ');
		str++;
		sscanf(str, "%ud", & h);
		while (h == 0 || h == gf2n::degree) {
			str = strchr(str, ' ');
			str++;
			sscanf(str, "%ud", & h);
		}
		gf2n::exponents[i] = h;
	}

	//---- set invsel and initialize variables for inversion and reduction

	if ((gf2n::anz_exponents == 1) &&
	    ((gf2n::degree - gf2n::exponents[0]) >= 8*sizeof(gf2n_word)) &&
	    (gf2n::exponents[0] >= 8*sizeof(gf2n_word))) {
		gf2n::invsel = 0;
		gf2n::exp1 = gf2n::exponents[0];
	}
	else
		if (gf2n::anz_exponents == 3) {
			unsigned int ma, mi;

			ma = max(gf2n::exponents[0], max(gf2n::exponents[1], gf2n::exponents[2]));
			mi = min(gf2n::exponents[0], min(gf2n::exponents[1], gf2n::exponents[2]));
			if ((gf2n::degree - ma) >= GF2N_WORDSIZE && mi >= GF2N_WORDSIZE) {
				gf2n::invsel = 1;
				gf2n::exp1 = gf2n::exponents[0];
				gf2n::exp2 = gf2n::exponents[1];
				gf2n::exp3 = gf2n::exponents[2];
			}
			else
				gf2n::invsel = 2;
		}
		else
			gf2n::invsel = 2;

	// ---- INITIIALIZE OPERATIONS  --------------------------

	if (!already_initialized) {
		gf2n::generate_mul_table ();
		gf2n::generate_square_table ();
	}

	delete[] gf2n::B; delete[] gf2n::C;
	delete[] gf2n::F; delete[] gf2n::G;

	gf2n::B = new gf2n_word [2*gf2n::anzBI];
	gf2n::C = new gf2n_word [2*gf2n::anzBI];
	gf2n::F = new gf2n_word [2*gf2n::anzBI];
	gf2n::G = new gf2n_word [2*gf2n::anzBI];

	if ((gf2n::B == NULL) || (gf2n::C == NULL) ||
	    (gf2n::F == NULL) || (gf2n::G == NULL)) {
		lidia_error_handler("gf2n", "gf2n_init::can't allocate memory");
		return;
	}

	// ---------- Initialisation of  I/O ----------------

	gf2nIO::setbase (mode);
	gf2nIO::setprefix (mode);
}



#ifdef LIDIA_NAMESPACE
}	// end of namespace LiDIA
#endif
