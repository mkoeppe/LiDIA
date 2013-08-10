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
#include	"LiDIA/gf2n.h"
#include	"LiDIA/finite_fields/gf2nIO.h"

#include	<cstring>
#include	<cctype>

const int STRLEN = 1024;


#ifdef LIDIA_NAMESPACE
namespace LiDIA {
#endif



void decin(gf2n & a, char* input_str);
void hexin(gf2n & a, char* input_str);

// ================================================================

std::istream & operator >> (std::istream & in_strm, gf2n & a)
{
	char line[STRLEN], c, *str;
	int i, length = 0;

	str = line;
        do
         {  
	  c = in_strm.get();
	  if (c == (char) EOF) {
		in_strm.setstate(std::ios::failbit|std::ios::eofbit);
		return in_strm;
        	}
         }
        while (isspace(c));

	while (isalpha(c) || isdigit(c) || c == ':') {
		length++;
		if (isalpha(c))
			*str++ = toupper(c);
		else
			*str++ = c;
		c = in_strm.get();
		if (c == (char) EOF)
			*str = '\0';
                if (length >= STRLEN)
                  lidia_error_handler("gf2nin","input buffer overflow");
	}
	line[length] = '\0';

	i = 0;
	while (line[i] != ':' && i < length)
		i++;

	if (i == length) {
		// no format given, take decimal input
		decin(a, line);
		return (in_strm);
	}
	else {
		str = line+i+1;
		line[i] = (char) NULL;

		if (strcmp(line, "HEX") == 0) {
			hexin(a, str);
			return (in_strm);
		}

		if (strcmp(line, "DEC") == 0) {
			decin(a, str);
			return (in_strm);
		}

		lidia_error_handler ("gf2n", " >> (...)::wrong base");
	}
	return (in_strm);
}



// ================================================================
// Function hexin                                                  
// ================================================================

void hexin(gf2n & a, char* input_str)
{
	unsigned int i;
	int j;
	bigint h(0);
	char c;

	for (i = 0; i < strlen(input_str); i++) {
		c = input_str[i];

		multiply(h, h, 16);
		j = c - '0';
		if (j >= 0 && j < 10) {
			add(h, h, j);
			continue;
		}
		c = toupper(c);

		j = c - 'A';
		if (j >= 0 && j < 6) {
			add(h, h, j+10);
			continue;
		}
		else {
			lidia_error_handler("gf2n", "hexin::wrong input character ");
			return;
		}
	}

	a.assign(h);
}



// ================================================================
// Function decin                                                  
// ================================================================

void decin(gf2n & a, char* input_str)
{
	bigint h(0);
	char c;
	unsigned int i;
	int  j;

	for (i = 0; i < strlen(input_str); i++) {
		c = input_str[i];
		if (c == ' ')
			continue;

		multiply(h, h, 10);
		j = c - '0';
		if (j >= 0 && j < 10) {
			add(h, h, j);
			continue;
		}
		else {
			lidia_error_handler("gf2n", "decin()::wrong input character ");
			return;
		}
	}

	a.assign(h);
}



#ifdef LIDIA_NAMESPACE
}	// end of namespace LiDIA
#endif
