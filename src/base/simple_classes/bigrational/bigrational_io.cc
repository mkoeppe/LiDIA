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
#include	"LiDIA/bigrational.h"
#include	<cctype>



#ifdef LIDIA_NAMESPACE
namespace LiDIA {
#endif



//
// input / output
//

std::istream &
operator >> (std::istream & in, bigrational & a)
{
	char s[10000];
	char *p = s;
	char c;

	a.num.assign_zero();
	a.den.assign_one();

	do {
		in.get(c);
	} while (isspace(c));
	if ((c == '+') || (c == '-')) {
		*p++ = c;
		do {
			in.get(c);
		} while (isspace(c));
	}
	if (!isdigit(c))
		std::cerr << "digit expected";
	while (isdigit(c)) {
		*p++ = c;
		in.get(c);
	}
	*p = '\0';
	string_to_bigint(s, a.num);

	if (c != '\n') {
		*s = '\0';
		p = s;
		while (isspace(c)) {
			in.get(c);
		}

		if (c == '/') {
			do {
				in.get(c);
			} while (isspace(c));
			if ((c == '+') || (c == '-')) {
				*p++ = c;
				do {
					in.get(c);
				} while (isspace(c));
			}
			if (!isdigit(c))
				std::cerr << "digit expected";
			while (isdigit(c)) {
				*p++ = c;
				in.get(c);
			}
			*p = '\0';
			string_to_bigint(s, a.den);
		}
	}
	if (c != '\n' && c != '\r')
		in.putback(c);
	a.normalize();
	return in;
}



std::ostream &
operator << (std::ostream & out, const bigrational & a)
{
	bigint n = a.num;
	bigint d = a.den;

	if (n.is_zero() || d.is_one())
		out << n;
	else
		out << n << "/" << d;
	return out;
}



int
string_to_bigrational (char *s, char *t, bigrational & a)
{
	long i = string_to_bigint(s, a.num);
	long j = string_to_bigint(t, a.den);
	a.normalize();
	return (i + j);
}



int
bigrational_to_string (const bigrational & a, char *s, char *t)
{
	int i = bigint_to_string(a.num, s);
	int j = bigint_to_string(a.den, t);
	return (i + j);
}



#ifdef C_STDIO

//
// using fread/fwrite
//

void
bigrational::read_from_file (FILE * fp)
{
	num.read_from_file(fp);
	den.read_from_file(fp);
	this->normalize();
}



void
bigrational::write_to_file (FILE * fp)
{
	num.write_to_file(fp);
	den.write_to_file(fp);
}



//
// using fscanf/fprintf
//

void
bigrational::scan_from_file (FILE * fp)
{
	char s[10000];
	char *p = s;
	char c;

	num.assign_zero();
	den.assign_one();

	do {
		c = getc(fp);
	} while (isspace(c));
	if ((c == '+') || (c == '-')) {
		*p++ = c;
		do {
			c = getc(fp);
		} while (isspace(c));
	}
	else {
		while (isspace(c)) {
			c = getc(fp);
		}
	}
	if (!isdigit(c))
		std::cerr << "digit expected";
	while (isdigit(c)) {
		*p++ = c;
		c = getc(fp);
	}
	*p = '\0';
	string_to_bigint(s, num);

	if (c != '\n') {
		*s = '\0';
		p = s;
		while (isspace(c)) {
			c = getc(fp);
		}

		if (c == '/') {
			do {
				c = getc(fp);
			} while (isspace(c));
			if ((c == '+') || (c == '-')) {
				*p++ = c;
				do {
					c = getc(fp);
				} while (isspace(c));
			}
			else {
				while (isspace(c)) {
					c = getc(fp);
				}
			}
			if (!isdigit(c))
				std::cerr << "digit expected";
			while (isdigit(c)) {
				*p++ = c;
				c = getc(fp);
			}
			*p = '\0';
			string_to_bigint(s, den);
		}
	}
	if (c != '\n' && c != '\r')
		ungetc(c, fp);
	this->normalize();
}



void
bigrational::print_to_file (FILE * fp)
{
	long l, k;
	char *s;

	if (num.is_zero() || den.is_one()) {
		l = num.bit_length();
		s = new char[l / 3 + 20];
		bigint_to_string(num, s);
	}
	else {
		l = num.bit_length();
		k = den.bit_length();
		if (l < k)
			l = k;
		s = new char[l / 3 + 20];
		bigint_to_string(num, s);
		fputs(s, fp);
		fputs("/", fp);
		bigint_to_string(den, s);
	}
	fputs(s, fp);
	delete[] s;
}

#endif	// C_STDIO



#ifdef LIDIA_NAMESPACE
}	// end of namespace LiDIA
#endif
