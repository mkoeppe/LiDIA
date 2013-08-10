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
#include	"LiDIA/bigfloat.h"
#include	<cstring>
#include	<cctype>



#ifdef LIDIA_NAMESPACE
namespace LiDIA {
#endif



#ifdef C_STDIO

int
bigfloat::print_to_file (FILE * fp)
{
	long sl = m.bit_length() + 1;
	if (e > 0)
		sl += e;
	sl = sl / 3 + 20;

	char *s = new char[sl];
	if (s == 0)
		lidia_error_handler("bigfloat", "print_to_file()::no memory available");

	int count = bigfloat_to_string(*this, s);

	for (long i = 0; i <= count; i++)
		putc(s[i], fp);

	delete[] s;

	return count;
}



int
bigfloat::scan_from_file (FILE * fp)
{
	static char s[100000];
	int count = 0;
	char *p = s;
	char c;

	do {
		c = getc(fp);
		count++;
	} while (isspace(c) || iscntrl(c));
	if ((c == '+') || (c == '-')) {
		*p++ = c;
		c = getc(fp);
		count++;
	}
	while (c == ' ') {
		c = getc(fp);
		count++;
	}
	if (!isdigit(c) && c != '.')
		lidia_error_handler("bigfloat", "scan_from_file()::digit/point expected");
	while (isdigit(c)) {
		*p++ = c;
		c = getc(fp);
		count++;
	}
	if (c == '.') {
		*p++ = c;
		c = getc(fp);
		count++;
		while (isdigit(c)) {
			*p++ = c;
			c = getc(fp);
			count++;
		}
	}
	while (c == ' ') {
		c = getc(fp);
		count++;
	}
	if (c == 'E' || c == 'e') {
		*p++ = c;
		c = getc(fp);
		count++;
		while (c == ' ') {
			c = getc(fp);
			count++;
		}
		if ((c == '+') || (c == '-')) {
			*p++ = c;
			c = getc(fp);
			count++;
		}
		while (c == ' ') {
			c = getc(fp);
			count++;
		}
		if (!isdigit(c))
			lidia_error_handler("bigfloat", "scan_from_file()::digit expected");
		while (isdigit(c)) {
			*p++ = c;
			c = getc(fp);
			count++;
		}
	}
	ungetc(c, fp);
	count--;
	*p = '\0';
	return static_cast<int>(string_to_bigfloat(s, *this));
}

#endif	// C_STDIO



void
bigfloat::print (char *msg)
{
	long l = m.bit_length() + 1;
	if (e > 0)
		l += e;
	if (l < bigfloat::binary_precision)
		l = bigfloat::binary_precision + 1;
	char *s = new char[l / 3 + 20];
	bigfloat_to_string(*this, s);
	std::cout << msg << s << std::flush;
	delete[] s;
}



std::ostream &
operator << (std::ostream & out, const bigfloat & a)
{
	long l = a.m.bit_length() + 1;
	if (a.e > 0)
		l += a.e;
	if (l < bigfloat::binary_precision)
		l = bigfloat::binary_precision + 1;
	char *s = new char[l / 3 + 20];
	bigfloat_to_string(a, s);
	out << s;
	delete[] s;
	return out;
}



std::istream &
operator >> (std::istream & in, bigfloat & a)
{
	char s[100000];
	char *p = s;
	char c;

	do {
		in.get(c);
	} while (isspace(c));
	if ((c == '+') || (c == '-')) {
		*p++ = c;
		do {
			in.get(c);
		} while (isspace(c));
	}
	else {
		while (isspace(c)) {
			in.get(c);
		}
	}

	if (!isdigit(c) && c != '.')
		lidia_error_handler("bigfloat", "cin::digit/point expected");
	while (isdigit(c)) {
		*p++ = c;
		in.get(c);
	}

	if (c == '.') {
		*p++ = c;
		in.get(c);
		while (isdigit(c)) {
			*p++ = c;
			in.get(c);
		}
	}
	while (isspace(c) && c != '\n') {
		in.get(c);
	}

	if (c == 'E' || c == 'e') {
		*p++ = c;
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
			lidia_error_handler("bigfloat", "cin::digit expected");
		while (isdigit(c)) {
			*p++ = c;
			in.get(c);
		}
	}
	in.putback(c);
	*p = '\0';
	string_to_bigfloat(s, a);
	return in;
}



#ifdef LIDIA_NAMESPACE
}	// end of namespace LiDIA
#endif
