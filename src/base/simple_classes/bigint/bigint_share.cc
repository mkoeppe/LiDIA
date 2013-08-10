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
//	Author	: Markus Maurer (MM)
//	Changes	: See CVS log
//
//==============================================================================================


#ifdef HAVE_CONFIG_H
# include	"config.h"
#endif
#include	"LiDIA/bigint.h"
#include	<cctype>
#include	<cstring>
#include        <vector>

#ifdef LIDIA_NAMESPACE
namespace LiDIA {
#endif



bool bigint::is_seeded = false;



//
// input / output
//


// number of characters in one line
// when printing a bigint


int bigint::chars_per_line = LIDIA_CHARS_PER_LINE;



//
// s has size old_size and will be resized to new_size.
//

void
bigint::allocate (char * &s, int old_size, int new_size)
{
	if (old_size > new_size) {
		old_size = new_size;
	}

	if (s == NULL) {
		old_size = 0;
	}

	char *t = new char[new_size];
	memory_handler(t, "bigint", "allocate::out of memory error");

	int i;
	for (i = 0; i < old_size; i++) {
		t[i] = s[i];
	}

	delete[] s;
	s = t;
}



//
// appends c to s at position pos.
// sz is the size of s
//

void
bigint::append_char (char * &s, int & sz, int pos, char c)
{
	if (pos > sz) {
		lidia_error_handler("bigint", "append_char::invalid argument");
	}

	if (pos == sz) {
		bigint::allocate(s, sz, 2*sz);
		sz *= 2;
	}

	s[pos] = c;
}



//
// skips '\\' followed by '\n'
//

int
bigint::skip_backslash_new_line (std::istream & in)
{
	int c = in.get();

	while (c == '\\') {
		c = in.get();

		if (c == '\n') {
			c = in.get();
		}
		else {
   		        in.setstate(std::ios::failbit);
//			lidia_error_handler("bigint::operator >>",
//					    "\\ must be immediately followed by new line.");
			return EOF;
		}
	}

	return c;
}



//
// reads in a bigint in ASCII-format from in
//

void
bigint::scan (std::istream & in)
{
        // g++ 2.95.2 bails out if isdigit(), isxdigit() are fully
        // qualified with std::
        using namespace std;

        if(!in.good()) {
	        return;
	}

	typedef std::vector<char> CharVector;
	CharVector buffer;
	buffer.reserve(1000);

	bool stop;
	int c;

	// skip blanks and line breaks
	do {
		c = in.get();
	}
	while (isspace(c));

	in.putback(c);


	// skip blanks, line breaks, and backslashs
	// followed by line break
	do {
		c = bigint::skip_backslash_new_line(in);

		// end of stream, complain

		if (c == EOF) {
   		        in.setstate(std::ios::failbit|std::ios::eofbit);
// 			lidia_error_handler("bigint::operator >>",
// 					    "bigint expected.");
			return;
		}
	} while (isspace(c));


	// handle sign
	if (c == '+' || c == '-') {
		buffer.push_back(c);
		c = bigint::skip_backslash_new_line(in);
	}

	std::ios::fmtflags streamflags = in.flags();
	int base = ((streamflags & std::ios::hex) ? 16 :
		    (streamflags & std::ios::oct) ? 8 :
		    (streamflags & std::ios::dec) ? 10 : 0);
	if(base == 0) { // input indicates basis
	    if(c != '0') {  // assume decimal input
		base = 10;
	    }
	    else {
		c = bigint::skip_backslash_new_line(in);
		if(c == 'x' || c == 'X') {
		    base = 16;
		    c = bigint::skip_backslash_new_line(in);
		}
		else if(isdigit(c) && c != '8' && c != '9') {
		    base = 8;
		}
		else {
		    // the relevant input is simply "0"
		    in.putback(c);
		    this->assign_zero();
		    return;
		}
	    }
	}
	if(base == 8) {
	        buffer.push_back('0');
	}
	else if(base == 16) {
	        buffer.push_back('0');
	        buffer.push_back('x');
	}
	
	// require digit now
	if ((base == 16 && isxdigit(c)) ||
	    (base == 10 && isdigit(c)) ||
	    (isdigit(c) && c != '8' && c != '9')) {
	        buffer.push_back(std::tolower(c));
		stop = false;
	}
	else {
		stop = true;
		in.setstate(std::ios::failbit);
		return;
	}


	while (!stop) {
		c = in.get();

		// store digit and continue
		if ((base == 16 && isxdigit(c)) ||
		    (base == 10 && isdigit(c)) ||
		    (isdigit(c) && c != '8' && c != '9')) {
		        buffer.push_back(std::tolower(c));
		}

		// if line break, stop
		else if (c == '\n') {
			stop = true;
		        in.putback(c);
		}

		// skip "\\\n".
		// stop if "\\c" and putback "\\c"
		else if (c == '\\') {
			c = in.get();

			if (c != '\n') {
				stop = true;
				in.putback(c);
				in.putback('\\');
			}
		}

		// stop if EOF
		else if (c == EOF) {
			stop = true;
                        // clear the failbit; EOF is no error at this point
                        in.clear(in.rdstate() & ~std::ios::failbit);
		}

		// stop otherwise and putback c
		else {
			stop = true;
		        in.putback(c);
		}
	}

	buffer.push_back('\0');
	string_to_bigint(&buffer[0], *this);
	return;
}



//
// Prints a bigint (represented by the string s)
// to out, where at most bigint::chars_per_line
// characters are printed in one line.
//

void
bigint::print (std::ostream & out, char const *s) const
{
	if (bigint::chars_per_line <= 1) {
		out << s;
	}
	else {
		int l = strlen(s);
		std::vector<char> v(s, s + l + 1); // keep '\0'
		
		std::vector<char>::size_type start, end;
		char c;

		for (start = 0; start < l; start += bigint::chars_per_line) {
			if (start + bigint::chars_per_line >= l) {
				out << &v[start];
			}
			else {
				end = start + bigint::chars_per_line;

				c = v[end];
				v[end] = '\0';

				out << &v[start];
				out << '\\';
				out << '\n';
				v[end] = c;
			}
		}
	}
}



#ifdef LIDIA_NAMESPACE
}	// end of namespace LiDIA
#endif
