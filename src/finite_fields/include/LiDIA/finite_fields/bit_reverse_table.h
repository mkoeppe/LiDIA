// -*- C++ -*-
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
//	Author	: Thorsten Rottschaefer (TR)
//                Victor Shoup (VS) and Thomas Pfahler (TPf)
//	Changes	:
//
//==============================================================================================


#ifndef LIDIA_BIT_REVERSE_TABLE_H_GUARD_
#define LIDIA_BIT_REVERSE_TABLE_H_GUARD_


#ifndef LIDIA_UDIGIT_H_GUARD_
# include	"LiDIA/udigit.h"
#endif



#ifdef LIDIA_NAMESPACE
namespace LiDIA {
# define IN_NAMESPACE_LIDIA
#endif



class bit_reverse_table
{
private:
        lidia_size_t** mem;
        lidia_size_t allocated;

        lidia_size_t rev_inc(lidia_size_t a, lidia_size_t k);
	//increases 'a' in "bitreverse order" (size: 'k' bits)
        lidia_size_t* table(lidia_size_t k);
	//returns mem[k], initializes if mem[k] == 0 or k > allocated
public:

	bit_reverse_table() : mem(0), allocated(0)
	{
		debug_handler("fft_table::bit_reverse_table", "bit_reverse_table()");
	};

        ~bit_reverse_table();

        void copy(udigit* A, const udigit* a, lidia_size_t k);
	//copies a in "bitreverse order" into A
};



#ifdef LIDIA_NAMESPACE
}	// end of namespace LiDIA
# undef IN_NAMESPACE_LIDIA
#endif



#endif	// LIDIA_BIT_REVERSE_TABLE_H_GUARD_
