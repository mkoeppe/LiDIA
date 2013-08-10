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
//	Author	:
//	Changes	: See CVS log
//
//==============================================================================================


#ifdef HAVE_CONFIG_H
# include	"config.h"
#endif
#include	"LiDIA/multi_bigmod.h"
#include	<cstring>



#ifdef LIDIA_NAMESPACE
namespace LiDIA {
#endif



//
// class variables
//

residue_class< bigint > * multi_bigmod::zero = NULL;



//
// modulus modifiers
//

void
multi_bigmod::set_modulus (const bigint & m)
{
	if (m.is_zero())
		lidia_error_handler ("multi_bigmod", "set_modulus::zero modulus");
	else {
		if (Mp != NULL)
			base_bigmod::L->clear(Mp);

		if (m.is_positive())
			Mp = base_bigmod::L->insert(m);
		else
			Mp = base_bigmod::L->insert(-m);
	}
}



void
multi_bigmod::initialize ()
{
	if (base_bigmod::L == NULL)
		base_bigmod::L = new residue_class_list< bigint >;

	if (zero == NULL) {
		bigint z;
		z.assign_zero();
		zero = base_bigmod::L->insert(z);
	}
}



//
// modifiers
//

void
multi_bigmod::randomize (const bigint & m)
{
	set_modulus(m);
	I.randomize(Mp->get_mod());
	if (I == Mp->get_mod())
		I.assign_zero();
}



void
multi_bigmod::randomize (const multi_bigmod & a)
{
	set_modulus(a.modulus());
	I.randomize(a.I);
	if (I == modulus())
		I.assign_zero();
}



//
// arithmetical procedures
//

void
add (multi_bigmod & c, const multi_bigmod & a, const multi_bigmod & b)
{
	if (a.Mp == b.Mp) {
		add (c, a, b, a.modulus());
		c.assign_modulus (a.Mp);
	}
	else
		lidia_error_handler ("multi_bigmod", "add::different moduli");
}



void
add (multi_bigmod & c, const multi_bigmod & a, const bigmod & b)
{
	if (bigmod::residue_class() == a.Mp) {
		add (c, a, b, bigmod::modulus());
		c.assign_modulus (bigmod::residue_class());
	}
	else
		lidia_error_handler ("multi_bigmod", "add::different moduli");
}



void
subtract (multi_bigmod & c, const multi_bigmod & a, const multi_bigmod & b)
{
	if (a.Mp == b.Mp) {
		subtract (c, a, b, a.modulus());
		c.assign_modulus (a.Mp);
	}
	else
		lidia_error_handler ("multi_bigmod", "subtract::different moduli");
}



void
subtract (multi_bigmod & c, const multi_bigmod & a, const bigmod & b)
{
	if (bigmod::residue_class() == a.Mp) {
		subtract (c, a, b, bigmod::modulus());
		c.assign_modulus (bigmod::residue_class());
	}
	else
		lidia_error_handler ("multi_bigmod", "subtract::different moduli");
}



void
multiply (multi_bigmod & c, const multi_bigmod & a, const multi_bigmod & b)
{
	if (a.Mp == b.Mp) {
		multiply (c, a, b, a.modulus());
		c.assign_modulus (a.Mp);
	}
	else
		lidia_error_handler ("multi_bigmod", "multiply::different moduli");
}



void
multiply (multi_bigmod & c, const multi_bigmod & a, const bigmod & b)
{
	if (bigmod::residue_class() == a.Mp) {
		multiply (c, a, b, bigmod::modulus());
		c.assign_modulus (bigmod::residue_class());
	}
	else
		lidia_error_handler ("multi_bigmod", "multiply::different moduli");
}



void
square (multi_bigmod & a, const multi_bigmod & b)
{
	square(a.I, b.I);
	remainder(a.I, a.I, (b.Mp)->get_mod());
	a.assign_modulus(b.Mp);
}



void
divide (multi_bigmod & c, const multi_bigmod & a, const multi_bigmod & b)
{
	if (a.Mp == b.Mp) {
		divide (c, a, b, a.modulus());
		c.assign_modulus (a.Mp);
	}
	else
		lidia_error_handler ("multi_bigmod", "divide::different moduli");
}



void
divide (multi_bigmod & c, const multi_bigmod & a, const bigmod & b)
{
	if (bigmod::residue_class() == a.Mp) {
		divide (c, a, b, bigmod::modulus());
		c.assign_modulus (bigmod::residue_class());
	}
	else
		lidia_error_handler ("multi_bigmod", "add::different moduli");
}



//
// I/O
//

void
multi_bigmod::read (std::istream & in)
{
	bigint m;
	char c;


	// read white spaces

	in >> c;
	while (c == ' ') in >> c;

	// read '('

	if (c != '(') {
		in.putback(c);
		in >> I;
		this->normalize();
	}
	else {
		// read mantissa
		in >> I;

		// read white spaces

		in >> c;
		while (c == ' ')
			in >> c;

		// read ','

		if (c != ',') {
			lidia_error_handler ("multi_bigmod::read(std::istream & in)",
					     "',' expected.");
		}
		else {
			// read modulus
			in >> m;
			base_bigmod::L->clear(Mp);
			Mp = base_bigmod::L->insert(m);

			this->normalize();

			// read white spaces

			in >> c;
			while (c == ' ')
				in >> c;

			// read ')'

			if (c != ')') {
				lidia_error_handler ("multi_bigmod::read(std::istream & in)",
						     "')' expected.");
			}
		}
	}
}



void
multi_bigmod::print (std::ostream & out) const
{
	out << "(" << I << ", " << Mp->get_mod() << ")";
}



std::istream &
operator >> (std::istream & in, multi_bigmod & a)
{
	a.read(in);
	return(in);
}



std::ostream &
operator << (std::ostream & out, const multi_bigmod & a)
{
	a.print(out);
	return out;
}



int
string_to_multi_bigmod (char *s, multi_bigmod & a)
{
	long l = strlen(s);
	long  n, i;
	char *h;
	bigint m;


	if (l == 0) {
		lidia_error_handler ("multi_bigmod::string_to_multi_bigmod(...)",
				     "Length of input string is zero.");
	}
	else if (s[0] != '(') {
		string_to_bigint(s, a.I);
		a.normalize();
	}
	else {
		// count the number n of ',' in s

		n = 0;
		i = 0;
		while (i < l) {
			if (s[i] == ',')
				n++;
			i++;
		}

		if (n == 0 || !(n&1) || s[l-1] != ')')
			lidia_error_handler ("multi_bigmod::string_to_multi_bigmod(...)",
					     "Invalid multi_bigmod format.");
		else {
			n = (n >> 1)+1;

			// find s[i], the n-th ',' in s
			i = 0;
			while (n > 0) {
				if (s[i] == ',')
					n--;
				i++;
			}
			i--;

			// extract mantissa

			s[i] = '\0';
			h = &(s[1]);
			string_to_bigint(h, a.I);
			s[i] = ',';

			// extract modulus

			s[l-1] = '\0';
			h = &(s[i+1]);
			string_to_bigint(h, m);
			s[l-1] = ')';

			a.set_modulus(m);
			a.normalize();
		}
	}

	return (strlen(s));
}



int
multi_bigmod_to_string (const multi_bigmod & a, char * &s)
{
	char *man, *mod;
	long  l;

	man = new char[a.mantissa().bit_length()/3 + 10];
	mod = new char[a.modulus().bit_length()/3 + 10];

	bigint_to_string(a.mantissa(), man);
	bigint_to_string(a.modulus(), mod);

	if (s)
		delete [] s;

	s = new char [strlen(man) + strlen(mod) + 4];

	s[0] = '(';
	s[1] = '\0';
	s = strcat (s, man);
	l = strlen(s);
	s[l] = ',';
	s[l+1] = '\0';
	s = strcat (s, mod);
	l = strlen(s);
	s[l] = ')';
	s[l+1] = '\0';

	delete[] man; delete[] mod;

	return (strlen(s));
}



#ifdef C_STDIO

//
// using fread/fwrite
//

void
multi_bigmod::read_from_file (FILE * fp)
{
	bigint m;

	I.read_from_file(fp);
	if (feof(fp)) {
		lidia_error_handler ("multi_bigmod::read_from_file(FILE*)",
				     "Invalid multi_bigmod format.");
	}
	else {
		m.read_from_file(fp);
		this->set_modulus(m);
		this->normalize();
	}
}



void
multi_bigmod::write_to_file (FILE * fp)
{
	I.write_to_file(fp);
	(Mp->get_mod()).write_to_file(fp);
}



//
// using fscanf/fprintf
//

void
multi_bigmod::scan_from_file (FILE * fp)
{
	char c;
	bigint m;

	fscanf (fp, "%c", &c);
	if (feof(fp)) {
		I.assign_zero();
	}
	else if (c != '(') {
		lidia_error_handler ("multi_bigmod::scan_from_file(FILE * fp)",
				     "'(' expected.");
	}
	else {
		I.scan_from_file(fp);
		fscanf (fp, "%c", &c);
		if (c != ',')
			lidia_error_handler ("multi_bigmod::scan_from_file(FILE * fp)",
					     "',' expected.");
		else {
			m.scan_from_file(fp);
			this->set_modulus(m);
			this->normalize();

			fscanf (fp, "%c", &c);
			if (c != ')')
				lidia_error_handler ("multi_bigmod::scan_from_file(FILE * fp)",
						     "')' expected.");
		}
	}
}



void
multi_bigmod::print_to_file (FILE * fp)
{
	fprintf (fp, "(");
	I.print_to_file(fp);
	fprintf (fp, ",");
	(Mp->get_mod()).print_to_file(fp);
	fprintf (fp, ")");
}

#endif	// C_STDIO



#ifdef LIDIA_NAMESPACE
}	// end of namespace LiDIA
#endif
