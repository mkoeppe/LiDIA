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
//	Author	: Frank Lehmann, Markus Maurer, Volker Mueller
//	Changes	: See CVS log
//
//==============================================================================================



#ifdef HAVE_CONFIG_H
# include	"config.h"
#endif
#include	"LiDIA/path.h"
#include	"LiDIA/meq_prime.h"
#include	"LiDIA/bigmod.h"
#include        "LiDIA/osstream.h"
#include	<cstdio>
#include	<cstdlib>
#include        <string>
#include        <vector>


#ifdef LIDIA_NAMESPACE
namespace LiDIA {
#endif



//  *****  definition of static members of class meq_prime  *****

const char meq_prime::ASCII  = 2;
const char meq_prime::BINARY = 4;
const char meq_prime::GZIP   = 1;
const char meq_prime::BIN_GZIP  = 3;
const char meq_prime::BZIP   = 5;
const char meq_prime::BIN_BZIP = 6;

const udigit meq_prime::MAX_MEQ_PRIME = 1019; // VM, has to be changed
const unsigned int BUFSIZE = 512 * 1024;


static void copy_file(char const file_in[], char const file_out[])
{
  FILE *fin, *fout;

  fin = fopen(file_in, "r");
  fout = fopen(file_out, "w");

  if (fin == NULL)
    lidia_error_handler("meq_prime", "can't open file for copy read");

  if (fout == NULL)
    lidia_error_handler("meq_prime", "can't open file for copy write");

  char buffer[BUFSIZE+1];
  size_t i;

  rewind(fout);
  rewind(fin);

  while ((i = fread((char*) buffer, sizeof(char), BUFSIZE, fin)) == BUFSIZE) {
    fwrite((char *) buffer, sizeof(char), BUFSIZE, fout);
  }
  
  fwrite((char *) buffer, sizeof(char), i, fout);
  fflush(fout);
  fclose(fin);
  fclose(fout);
}

//--------------------------------------------------------------

meq_prime::meq_prime ()
{
	debug_handler ("meq_prime", "meq_prime (void)");

	l = 0;
	v = 0;
	s = 0;

	type = mode = 0;
	hecke_op1 = hecke_op2 = 0;
	file_format = ASCII;
	fp_me  = NULL;
	me_row =    0;
	row.set_mode(EXPAND);
}



meq_prime::meq_prime (const meq_prime & p)
{
	debug_handler ("meq_prime", "meq_prime (meq_prime)");

	l = p.l;
	v = p.v;
	s = p.s;

	type = p.type;
	mode = p.mode;
	hecke_op1 = p.hecke_op1;
	hecke_op2 = p.hecke_op2;

	file_format = p.file_format;
	if (fp_me != NULL) {
		fclose(fp_me);
		if (file_format == GZIP || file_format == BIN_GZIP ||
		    file_format == BZIP || file_format == BIN_BZIP)
			std::remove(temp_file);
	}

	fp_me  = p.fp_me;
	me_row = p.me_row;
	row.set_mode(EXPAND);
	row.set_capacity(l+2);
}



meq_prime::~meq_prime ()
{
	debug_handler ("meq_prime", "destructor");

	if (fp_me != NULL) {
		fclose (fp_me);
		if (file_format == GZIP || file_format == BIN_GZIP ||
		    file_format == BZIP || file_format == BIN_BZIP)
			std::remove(temp_file);
		fp_me = NULL;
	}
}



//
// Returns true, if the prime is found in the database
// (by examine_prime) and the meq file exists on disk.
//

bool meq_prime
::set_prime (meq_prime_t p)
{
	debug_handler ("meq_prime", "set_prime (meq_prime_t)");

	if (!is_prime(p) || p <= 2)
		return false;

	if (p != l) {
		if (fp_me != NULL) {
			fclose (fp_me);
			if (file_format == GZIP || file_format == BIN_GZIP ||
			    file_format == BZIP || file_format == BIN_BZIP) {
				std::remove(temp_file);
			}
			fp_me = NULL;
		}

		l = p;
		v = 0;
		s = 0;

		type = mode = 0;
		hecke_op1 = hecke_op2 = 0;

		if (examine_prime() == true) {
			file_format = ASCII;
			fp_me  = NULL;
			me_row =    0;
			row.kill();
			row.set_capacity(l+2);

			return meq_file_is_available();
		}
		else
			return false;
	}
	else {
		if (fp_me != NULL) {
			fseek (fp_me, 0, SEEK_SET);
			me_row = 0;
		}
		return true;
	}
}



meq_prime & meq_prime
::operator = (meq_prime_t p)
{
	debug_handler ("meq_prime", "operator = (meq_prime_t)");

	if (!is_prime(p) || p <= 2)
		lidia_error_handler("meq_prime", "operator=::p is no prime");

	if (p != l) {
		if (fp_me != NULL) {
			fclose (fp_me);
			if (file_format == GZIP || file_format == BIN_GZIP ||
			    file_format == BZIP || file_format == BIN_BZIP)
				std::remove(temp_file);
		}

		l = p;
		v = 0;
		s = 0;

		type = mode = 0;
		hecke_op1 = hecke_op2 = 0;

		examine_prime ();
		file_format = ASCII;

		fp_me  = NULL;
		me_row =    0;
		row.kill();
		row.set_capacity(l+2);
	}
	else {
		if (fp_me != NULL) {
			fseek (fp_me, 0, SEEK_SET);
			me_row = 0;
		}
	}

	return *this;
}



meq_prime &
meq_prime::operator = (const meq_prime & p)
{
	debug_handler ("meq_prime", "operator = (meq_prime)");

	if (fp_me != NULL) {
		fclose (fp_me);
		if (file_format == GZIP || file_format == BIN_GZIP ||
		    file_format == BZIP || file_format == BIN_BZIP)
			std::remove(temp_file);
	}

	l = p.l;
	v = p.v;
	s = p.s;

	type = p.type;
	mode = p.mode;
	hecke_op1 = p.hecke_op1;
	hecke_op2 = p.hecke_op2;

	file_format = p.file_format;
	fp_me  = p.fp_me;
	me_row = p.me_row;
	row.set_capacity(l+2);

	return *this;
}



meq_prime_t
meq_prime::valence ()
{
	debug_handler ("meq_prime", "valence (void)");

	if (type == 1)
		return v;
	else
		return (v / 2);
}



int
meq_prime::trans_type ()
{
	debug_handler ("meq_prime", "trans_type (void)");
	return type;
}



int meq_prime
::reset ()
{
	debug_handler ("meq_prime", "reset (void)");

	int rc = 0;

	if (fp_me == NULL) {
		rc = ropen_meq_file ();
	}

	if (!rc) {
		fseek (fp_me, 0, SEEK_SET);
	}
	else
		lidia_error_handler ("meq_prime", "reset()::file open error");

	return rc;
}


std::string const meq_prime::bin_suffix(".bin");
std::string const meq_prime::gz_suffix(".gz");
std::string const meq_prime::bin_gz_suffix(".bin.gz");
std::string const meq_prime::bz2_suffix(".bz2");
std::string const meq_prime::bin_bz2_suffix(".bin.bz2");


void meq_prime
::get_filename_vector(std::vector<std::string>& fn_vec) {
	debug_handler ("meq_prime", "get_filename_vector()");

	osstream oss;
	char const* name_env = getenv("LIDIA_ECO");
	if (name_env != NULL && name_env[0] != '\0') {
	    oss << name_env;
	}
	else {
	    oss << LIDIA_ECO_DB;
	}
	oss << MEQ_PREFIX << l;
	std::string name_database = extractString(oss);

	if (type == 2) {
	    name_database += "_II";
	}
	else if(type != 1) {
	    lidia_error_handler ("meq_prime::meq_file_is_available()",
				 "Unknown type.");
	}

	fn_vec.resize(0);
	fn_vec.reserve(6);
	fn_vec.push_back(name_database);
	fn_vec.push_back(name_database + meq_prime::bin_suffix);
	fn_vec.push_back(name_database + meq_prime::gz_suffix);
	fn_vec.push_back(name_database + meq_prime::bin_gz_suffix);
	fn_vec.push_back(name_database + meq_prime::bz2_suffix);
	fn_vec.push_back(name_database + meq_prime::bin_bz2_suffix);
}


bool meq_prime
::meq_file_is_available()
{
	debug_handler ("meq_prime", "meq_file_is_available()");

	typedef std::vector<std::string> SVector;
	SVector fname;
	this->get_filename_vector(fname);
	for (SVector::const_iterator iter = fname.begin(),
		 endIter = fname.end();
	     iter != endIter; ++iter) {
	    if (meq_prime::fexist(iter->c_str()))
		return true;
	}

	return false;
}



int meq_prime
::ropen_meq_file ()     // open file for reading
{
  debug_handler ("meq_prime", "ropen_meq_file()");
  
  int    i, rc = 0;
  bool   found;
  char   name_database[2048], name[2048], command[2048];

  std::vector<std::string> fname;
  get_filename_vector(fname);
  
  if (fp_me != NULL) {
    fclose(fp_me);
    fp_me = NULL;
  }
  
  found = false;
  for (i = 0; i < 6 && !found; i++) {
    if (fexist (fname[i].c_str())) {
      switch (i) {
      case 0 :
	fp_me = fopen (fname[i].c_str(), "r");
	file_format = ASCII;
	break;
	
      case 1 :
	fp_me = fopen (fname[i].c_str(), "rb");
	file_format = BINARY;
	break;
	
      case 2 :
#ifdef HAVE_MKSTEMP
	strcpy(temp_file,"meq_tempXXXXXX");
	if (mkstemp(temp_file) == -1)
	  lidia_error_handler("meq_prime", "r_open_meq_file"
			      "::can't generate temporary"
			      " file");
#else
	if (tmpnam(temp_file)  == NULL)
	  lidia_error_handler("meq_prime", "r_open_meq_file"
			      "::can't generate temporary"
			      " file");
#endif
#ifdef HAVE_GUNZIP
	strcpy(name, temp_file);
	strcat(name, ".gz");
	copy_file(fname[i].c_str(), name);
	sprintf(command, "gunzip %s", name);
	std::system(command);
	fp_me = fopen (temp_file, "r");
	file_format = GZIP;
#else
	lidia_error_handler("meq_prime","gunzip seems to be not installed"
			    " on your computer, format not usable");
#endif
	break;
	
      case 3 :
#ifdef HAVE_MKSTEMP
	strcpy(temp_file,"meq_tempXXXXXX");
	if (mkstemp(temp_file) == -1)
	  lidia_error_handler("meq_prime", "r_open_meq_file"
			      "::can't generate temporary"
			      " file");
#else
	if (tmpnam(temp_file)  == NULL)
	  lidia_error_handler("meq_prime", "r_open_meq_file"
			      "::can't generate temporary"
			      " file");
#endif
#ifdef HAVE_GUNZIP
	strcpy(name, temp_file);
	strcat(name, ".bin.gz");
	copy_file(fname[i].c_str(), name);
	sprintf(command, "gunzip %s", name);
	std::system(command);
	fp_me = fopen (temp_file, "rb");
	file_format = BIN_GZIP;
#else
	lidia_error_handler("meq_prime","gunzip seems to be not installed"
			    " on your computer, format not usable");
#endif
	break;
	
      case 4 :
#ifdef HAVE_MKSTEMP
	strcpy(temp_file,"meq_tempXXXXXX");
	if (mkstemp(temp_file) == -1)
	  lidia_error_handler("meq_prime", "r_open_meq_file"
			      "::can't generate temporary"
			      " file");
#else
	if (tmpnam(temp_file)  == NULL)
	  lidia_error_handler("meq_prime", "r_open_meq_file"
			      "::can't generate temporary"
			      " file");
#endif
#ifdef HAVE_BUNZIP2
	strcpy(name, temp_file);
	strcat(name, ".bz2");
	copy_file(fname[i].c_str(), name);
	sprintf(command, "bunzip2 %s", name);
	std::system(command);
	fp_me = fopen (temp_file, "r");
	file_format = BZIP;
#else
	lidia_error_handler("meq_prime","bunzip2 seems to be not installed"
			    " on your computer, format not usable");
#endif
	break;
	
      case 5 :
#ifdef HAVE_MKSTEMP
	strcpy(temp_file,"meq_tempXXXXXX");
	if (mkstemp(temp_file) == -1)
	  lidia_error_handler("meq_prime", "r_open_meq_file"
			      "::can't generate temporary"
			      " file");
#else
	if (tmpnam(temp_file)  == NULL)
	  lidia_error_handler("meq_prime", "r_open_meq_file"
			      "::can't generate temporary"
			      " file");
#endif
#ifdef HAVE_BUNZIP2
	strcpy(name, temp_file);
	strcat(name, ".bin.bz2");
	copy_file(fname[i].c_str(), name);
	sprintf(command, "bunzip2 %s", name);
	std::system(command);
	fp_me = fopen (temp_file, "rb");
	file_format = BIN_BZIP;
#else
	lidia_error_handler("meq_prime","bunzip2 seems to be not installed"
			    " on your computer, format not usable");
#endif
	break;
	
      default :
	lidia_error_handler("meq_prime", "ropen_meq_file::format "
			    "not recognized");
	break;
      }
      me_row = 0; // number of first line in file
      rc = 0;
      found  = true;
    }
  }
  if ((!found) || (fp_me == NULL)) {
    char s[1000];
    sprintf(s, "ropen_meq_file::can't open meq file for l = %ld", l);
    lidia_error_handler("meq_prime", s);
    rc = 1;
  }
  
  return rc;
}



int
meq_prime::close_meq_file (void)
{
	debug_handler ("meq_prime", "close_meq_file (void)");

	if (fp_me != NULL) {
		fclose (fp_me);
		if (file_format == GZIP || file_format == BIN_GZIP ||
		    file_format == BZIP || file_format == BIN_BZIP)
			std::remove(temp_file);
		fp_me = NULL;
	}

	me_row = 0;
	return 0;
}



int
meq_prime::fexist (const char *fname)
{
	FILE *fp;

	if ((fp = fopen (fname, "r")) != NULL) {
		fclose (fp);
		return 1;
	}
	else
		return 0;
}



int
meq_prime::read_row_ascii (base_vector< bigint > & vec)
{
	debug_handler ("meq_prime", "read_row_ascii (base_vector &)");

	int  i;

	if (fp_me == NULL) {
		lidia_error_handler ("meq_prime", "read_row_ascii(base_vector &)::file not yet open");
		return 0;
	}

	if (me_row == static_cast<int>(l+2)) {
		// after last row
		fseek (fp_me, 0, SEEK_SET);
		me_row = 0;
	}

	vec.set_capacity (v + 1);

	for (i = 0; i <= static_cast<int>(v); i++)
		vec[i].scan_from_file (fp_me);

	me_row ++;
	return 0;
}



int
meq_prime::read_row_ascii (base_vector< bigint > & vec, const bigint & p)
{
	debug_handler ("meq_prime", "read_row_ascii (base_vector &, const bigint &)");

	int    i, rc;

	rc = read_row_ascii(vec);

	for (i = 0; i <= static_cast<int>(v); i++) {
		remainder (vec[i], vec[i], p);
		if (vec[i].is_negative())
			add (vec[i], vec[i], p);
	}
	return rc;
}



int
meq_prime::read_row_bin   (base_vector< bigint > & vec)
{
	debug_handler ("meq_prime", "read_row_bin (base_vector &)");

	int  i;

	if (fp_me == NULL) {
		lidia_error_handler ("meq_prime", "read_row_bin(base_vector &)::file not yet opened");
		return 0;
	}

	if (me_row == static_cast<int>(l+2)) {
		fseek (fp_me, 0, SEEK_SET);
		me_row = 0;
	}

	vec.set_capacity (v + 1);

	for (i = 0; i <= static_cast<int>(v); i++) {
		vec[i].read_from_file (fp_me);
	}

	me_row ++;
	return 0;
}



int
meq_prime::read_row_bin (base_vector< bigint > & vec, const bigint & p)
{
	debug_handler ("meq_prime", "read_row_bin (base_vector &, const bigint &)");

	int    i, rc;

	rc = read_row_bin(vec);

	for (i = 0; i <= static_cast<int>(v); i++) {
		remainder (vec[i], vec[i], p);
		if (vec[i].is_negative())
			add (vec[i], vec[i], p);
	}
	return rc;
}



int
meq_prime::read_row ()
{
	debug_handler ("meq_prime", "read_row ()");
	int rc = 0;

	if (fp_me == NULL) {
		rc = ropen_meq_file ();
	}

	if (!rc) {
		if (file_format == ASCII)
			rc = read_row_ascii (row);
		else if (file_format == BINARY)
			rc = read_row_bin (row);
		else {
			lidia_error_handler ("meq_prime", "read_row::invalid file format");
			rc = 1;
		}
	}

	return rc;
}



int
meq_prime::read_row   (const bigint & p)
{
	debug_handler ("meq_prime", "read_row (const bigint &)");

	int  rc = 0;

	if (fp_me == NULL)
		rc = ropen_meq_file ();


	if (! rc) {
		if (file_format == ASCII  || file_format == GZIP || file_format == BZIP)
			rc = read_row_ascii (row, p);
		else
			if (file_format == BINARY || file_format == BIN_GZIP ||
			    file_format == BIN_BZIP)
				rc = read_row_bin (row, p);
			else {
				lidia_error_handler ("meq_prime", "read_row (const bigint &)::invalid file format");
				rc = 1;
			}
	}
	return rc;
}



void
meq_prime::build_poly_in_X (Fp_polynomial & p, const bigmod & y)
{
	debug_handler ("meq_prime", "build_poly_in_X ()");

	bigint m;
	bigmod c;
	bigmod ypow;
	int i, j;

	m = bigmod::modulus ();
	p.set_modulus(m);
	p.assign_zero();
	meq_prime::reset();

	if (y.is_zero()) {
#ifdef AIX_DEBUG
		std::cout << "meq_prime::build_poly_in_X:: y is zero." << std::endl;
#endif
		for (i = 0; i <= static_cast<int>(l+1); i++) {
			read_row (m);
			p.set_coefficient (row[0], i);
		}
	}
	else {
#ifdef AIX_DEBUG
		std::cout << "meq_prime::build_poly_in_X:: y is NOT zero." << std::endl;
		std::cout << "meq_prime::build_poly_in_X:: v is " << v << std::endl;
#endif
		for (i = 0; i <= static_cast<int>(l+1); i++) {
			read_row (m);
			c.assign(row[v]);

			for (j = static_cast<int>(v-1); j >= 0; j--) {
				multiply(c, c, y);
				add(c, c, row[j]);
			}
			p.set_coefficient (c.mantissa(), i);
		}
	}
}



void
meq_prime::build_poly_in_Y (Fp_polynomial & p, const bigmod & x)
{
	debug_handler ("meq_prime", "build_poly_in_Y()");

	bigint m(bigmod::modulus());
	bigmod c;
	bigmod xpow;
	int i, j;

	meq_prime::reset ();

	p.set_modulus(m);
	p.assign_zero();

	if (x.is_zero()) {
		read_row (m);

		for (j = 0; j <= static_cast<int>(v); j++)
			p.set_coefficient (row[j], j);
	}
	else {
		read_row (m);

		for (j = 0; j <= static_cast<int>(v); j++)
			p.set_coefficient (row[j], j);

		xpow.assign(x);

		for (i = 1; i <= static_cast<int>(l+1); i++) {
			read_row (m);

			for (j = 0; j <= static_cast<int>(v); j++) {
				add(c, static_cast<bigint>(p[j]), row[j] * xpow);
				p.set_coefficient (c.mantissa(), j);
			}
			multiply(xpow, xpow, x);
		}
	}
}



#ifdef LIDIA_NAMESPACE
}	// end of namespace LiDIA
#endif
