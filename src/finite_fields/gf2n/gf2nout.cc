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



#ifdef LIDIA_NAMESPACE
namespace LiDIA {
#endif



std::ostream & operator << (std::ostream & out_strm, const gf2n & a)
{
	bigint h(0);
	char *s;
	int i;
	long k;

	if (a.is_reduced() == false)
		partial_reduce2[gf2n::invsel](a.element);

	for (i = gf2n::anzBI-1; i >= 0; i--) {
		// convert gf2n element into bigint

		// this hack was only introduced to test whether the 64bit patches
		// work (to test them on a machine with 32 bit register length).
		// no harm will be done if you remove it.
		// -rpw
#if GF2N_WORDSIZE == 64 && SIZEOF_LONG != 8
		shift_left(h, h, 32);
		add(h, h, bigint(static_cast<long>(a.element[i] >> 32)));

		shift_left(h, h, 32);
		add(h, h, bigint(static_cast<long>(a.element[i] & 0xffffffffL)));
#else
		shift_left(h, h, CHAR_BIT*sizeof(gf2n_word));
		add(h, h, bigint(a.element[i]));
#endif
	}

	if (gf2nIO::ioprefix != NULL)
		out_strm << gf2nIO::ioprefix << ':';

	s = new char[h.bit_length()];

	switch (gf2nIO::iobase) {
	case gf2nIO::Hex :

		if (h.is_zero())
			out_strm << '0';
		else {
			i = 0;
			while (!h.is_zero()) {
				remainder(k, h, 16);
				shift_right(h, h, 4);
				if (k < 10)
					s[i++] = static_cast<char>(static_cast<long>('0') + k);
				else
					s[i++] = static_cast<char>(static_cast<long>('a') + k - 10);
			}
			for (k = i-1; k >= 0; k--)
				out_strm << s[k];
		}
		break;

	case gf2nIO::Dec :
		if (h.is_zero())
			out_strm << '0';
		else {
			i = 0;
			while (!h.is_zero()) {
				remainder(k, h, 10);
				divide(h, h, 10);
				s[i++] = static_cast<char>(static_cast<long>('0') + k);
			}
			for (k = i-1; k >= 0; k--)
				out_strm << s[k];
		}
		break;

	default:
		lidia_error_handler("gf2n", " << (...)::wrong base");
		break;
	}
	delete[] s;
	return out_strm;
}



#ifdef LIDIA_NAMESPACE
}	// end of namespace LiDIA
#endif
