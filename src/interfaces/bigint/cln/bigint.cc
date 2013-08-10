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
//	$Id: bigint.cc,v 2.10 2004/06/15 10:19:48 lidiaadm Exp $
//
//	Author	: Thomas Papanikolaou (TP), Bruno Haible (HB)
//	Changes	: See CVS log
//
//==============================================================================================


#ifdef HAVE_CONFIG_H
# include	"config.h"
#endif
#include	"LiDIA/bigint.h"
#include	"LiDIA/random_generator.h"
#include        "LiDIA/osstream.h"
#include        <vector>
#include	<cctype>
#include        <cstring>
#include        <cassert>



#ifdef LIDIA_NAMESPACE
namespace LiDIA {
#endif



#if !INLINE_INTERFACE
# define inline
# include	"LiDIA/kernel/bigint_interface.h"
# undef inline
#endif



void
power (bigint & c, const bigint & a, const bigint & b)
{
	if ((a.I == 1) || (a.I == -1)) {
		if (b.is_odd()) {
			c.I = a.I;
		}
		else {
			c.assign_one();
		}
	}
	else if (b.is_negative()) {
		c.assign_zero();
	}
	else if (b.is_zero()) {
		c.assign_one();
	}
	else {
		c.I = expt_pos(a.I, b.I);
	}
}



void
power (bigint & c, const bigint & a, long b)
{
	if ((a.I == 1) || (a.I == -1)) {
		if (b & 1) {
			c.I = a.I;
		}
		else {
			c.assign_one();
		}
	}
	else if (b < 0) {
		c.assign_zero();
	}
	else if (b == 0) {
		c.assign_one();
	}
	else {
		c.I = expt_pos(a.I, cln::cl_I(b));
	}
}



bigint
xgcd (bigint & u, bigint & v, const bigint & a, const bigint & b)
{
	if (zerop(a.I)) {
		u.I = 0;
		v = b.sign();
		return abs(b);
	}
	if (zerop(b.I)) {
		v.I = 0;
		u = a.sign();
		return abs(a);
	}
	if (abs(a.I) == abs(b.I)) {
		u.I = 0;
		v = b.sign();
		return abs(b);
	}

	cln::cl_I g = xgcd(a.I, b.I, &u.I, &v.I);

	if (minusp(g)) {
		g = -g;
		u.I = - u.I;
		v.I = - v.I;
	}

	if (abs(2*u.I*g) > abs(b.I)) {
		if (u.sign() == b.sign()) {
			u.I = u.I - exquo(b.I, g);
			v.I = v.I + exquo(a.I, g);
		}
		else {
			u.I = u.I + exquo(b.I, g);
			v.I = v.I - exquo(a.I, g);
		}
	}
	else if (abs(2*v.I*g) > abs(a.I)) {
		if (v.sign() == a.sign()) {
			u.I = u.I + exquo(b.I, g);
			v.I = v.I - exquo(a.I, g);
		}
		else {
			u.I = u.I - exquo(b.I, g);
			v.I = v.I + exquo(a.I, g);
		}
	}

	return bigint(g);
}



//
// input / output
//


std::istream &
operator >> (std::istream & in, bigint & a)
{
	a.scan (in);
	return in;
}



std::ostream &
operator << (std::ostream & out, const bigint & a)
{
    if(out.good()) 
    {
	std::ios::fmtflags streamflags = out.flags();
	int base = ((streamflags & std::ios::hex) ? 16 :
		    (streamflags & std::ios::oct) ? 8 : 10);
	bool printbase = streamflags & std::ios::showbase;

	osstream oss;
	if(a.is_negative()) {
	    oss << '-';
	}
	if( base == 8) {
	    if(printbase) {
		oss << '0';
	    }
	    fprintoctal(oss, abs(a.I));
	}
	else if(base == 16) {
	    if(printbase) {
		oss << "0x";
	    }
	    fprinthexadecimal(oss, abs(a.I));
	}
	else {
	    fprintdecimal(oss, abs(a.I));
	}
	std::string output = extractString(oss);

	a.print(out, output.c_str());
// Changed by G.A. - we must return in any case, otherwise gcc 4.1.2
// misunderstands this trying to optimize and the program crashes!
	//return out;
    }
    return out;
}


int
string_to_bigint (const char *s, bigint & a)
{
    assert(s);

    char const* ss = s;
    bool negative = false;
    if(ss[0] == '+' || ss[0] == '-') {
	negative = (ss[0] == '-');
	++ss;
    }

    if(ss[0] == '0' && ss[1] == '\0') {
	a.assign_zero();
	return (ss - s) + 1;
    }
    // now any leading '0' indicates a base != 10

    int base = 10;
    if(ss[0] == '0') {
	if(ss[1] == 'x' || ss[1] == 'X') {
	    base = 16;
	    ss += 2;
	}
	else {
	    base = 8;
	    ++ss;
	}
    }

    if(ss[0] == '+' || ss[0] == '-') {
	negative = (negative != (ss[0] == '-'));
	++ss;
    }
    
    std::vector<char> v;
    v.reserve(1000);

    if(base == 8) {
	v.push_back('#');
	v.push_back('o');
    }
    else if(base == 16) {
	v.push_back('#');
	v.push_back('x');
    }

    if(negative) {
	v.push_back('-');
    }

    bool finished = false;
    while(ss[0] != '\0' && !finished) {
	switch(ss[0]) {
	    case '0': case '1': case '2': case '3': case '4':
	    case '5': case '6': case '7':
		v.push_back(ss[0]);
		++ss;
		break;
	    case '8': case '9':
		if(base != 8) {
		    v.push_back(ss[0]);
		    ++ss;
		}
		else {
		    finished = true;
		}
		break;
	    case 'A': case 'B': case 'C': case 'D': case 'E': case 'F':
	    case 'a': case 'b': case 'c': case 'd': case 'e': case 'f':
		if(base == 16) {
		    v.push_back(std::tolower(ss[0]));
		    ++ss;
		}
		else {
		    finished = true;
		}
		break;
	    default:
		finished = true;
	}
    }
    v.push_back('\0');

    a = bigint(cln::cl_I(&v[0]));
    return ss - s;
}



int
bigint_to_string (const bigint & a, char *s)
{
	char* r = print_integer_to_string(10, a.I);

	cl_decimal_string(a.I);
	strcpy(s, r);
	cln::free_hook(r);
	return strlen(s);
}



#ifdef C_STDIO

//
// using fread/fwrite
//


void
bigint::read_from_file (FILE * fp)
{
	scan_from_file(fp);
	int c = getc(fp);
	if (c != EOF && c != '\n') {
		ungetc(c, fp);
	}
}



void
bigint::write_to_file (FILE * fp)
{
	print_to_file(fp);
	putc('\n', fp);
}



//
// using fscanf/fprintf
//


void
bigint::scan_from_file (FILE * fp)
{
	int alloc = 32; // must be > 0
	char* buffer = 0;
	bigint::allocate(buffer, static_cast<int>(0), alloc);
	int index = 0;

	int c;
	while ((c = getc(fp)) != EOF)
		if (!isspace(c))
			break;

	switch(c) {
	case EOF:
		lidia_error_handler("bigint", "scan_from_file(...)::end of file");
	case '-':
		bigint::append_char(buffer, alloc, index++, c);
	case '+':
		break;
	default :
		ungetc(c, fp);
		break;
	}

	while ((c = getc(fp)) != EOF) {
		if (isdigit(c))
			bigint::append_char(buffer, alloc, index++, c);
		else {
			ungetc(c, fp);
			break;
		}
	}

	bigint::append_char(buffer, alloc, index, '\0');
	string_to_bigint(buffer, *this);

	delete[] buffer;
}



void
bigint::print_to_file (FILE * fp)
{
	char* s = print_integer_to_string(10, I);

	fprintf(fp, "%s", s);
	cln::free_hook(s);
}



#endif	// C_STDIO



//
// random numbers
//

void
bigint::seed ()
{
	unsigned long	hi, lo;
	random_generator	rg;

	rg >> hi >> lo;

	cln::default_random_state.seed.lo = hi;
	cln::default_random_state.seed.hi = lo;

	is_seeded = true;
}



void
bigint::seed (const bigint & a)
{
	cln::default_random_state.seed.lo = cl_I_to_UL(ldb(a.I,
							 cln::cl_byte(32, 0)));
	cln::default_random_state.seed.hi = cl_I_to_UL(ldb(a.I,
							 cln::cl_byte(32, 32)));

	is_seeded = true;
}



void
bigint::randomize (const bigint & a)
{
	if (!is_seeded) {
		seed();
	}
	if (minusp(a.I)) {
		*this = bigint(-random_I(-a.I));
	}
	else {
		*this = bigint(random_I(a.I));
	}
}



//
// Error handler
//


void
cl_abort (void)
{
	lidia_error_handler("cln", "cl_abort() called");
}



#ifdef LIDIA_NAMESPACE
}	// end of namespace LiDIA
#endif
