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
//	Author	: Patrick Theobald (PT)
//	Changes	: See CVS log
//
//==============================================================================================


#ifndef LIDIA_VECTOR_FLAGS_H_GUARD_
#define LIDIA_VECTOR_FLAGS_H_GUARD_



#ifndef LIDIA_ARITH_INL_GUARD_
# include	"LiDIA/arith.inl"
#endif



#ifdef LIDIA_NAMESPACE
namespace LiDIA {
# define IN_NAMESPACE_LIDIA
#endif



class vector_flags
{

public:

	unsigned long print_mode;
	unsigned long info_mode;
	unsigned long sort_mode;
	float dyn_exp_ratio;

public:

	//
	// PRINT_MODE SETTINGS
	//

	enum {
		beauty_mode = 0,
		lidia_mode = 1,
		gp_mode = 2,
		maple_mode = 3,
		mathematica_mode = 4,
		kash_mode = 5,
		latex_mode = 6,
		default_print_mode = beauty_mode
	};

	//
	// INFO MODE SETTINGS
	//

	enum info_settings {
		fixed = 1,
		expand = 2,
		default_info_mode = fixed
	};

	//
	// SORT MODE SETTINGS
	//

	// MM, for AIX start with 1
	enum {
		sort_vector_def = 1,
		sort_vector_up = 2,
		sort_vector_down = 3,
		sort_vector_cmp = 4
	};

	//
	// constructors
	//

	vector_flags()
	{
		debug_handler ("vector_flags", "vector_flags()");
		print_mode = vector_flags::default_print_mode;
		info_mode = vector_flags::default_info_mode;
		sort_mode = 1;
		dyn_exp_ratio = 2.0;
	}

	vector_flags(info_settings md)
	{
		debug_handler ("vector_flags", "vector_flags(info_settings)");
		print_mode = vector_flags::default_print_mode;
		info_mode = md;
		sort_mode = 1;
		dyn_exp_ratio = 2.0;
	}

	vector_flags(unsigned int md)
	{
		debug_handler ("vector_flags", "vector_flags(unsigned int)");
		print_mode = vector_flags::default_print_mode;
		if (md == 'E')
			info_mode = 2;
		if (md == 'F')
			info_mode = 1;
		sort_mode = 1;
		dyn_exp_ratio = 2.0;
	}

	//
	// destructor
	//

	~vector_flags()
	{
		debug_handler ("vector_flags", "~vector_flags()");
	}

	//
	// set info_mode / get info_mode
	//

public:

	void set_info_mode(unsigned long md)
	{
		debug_handler ("vector_flags", "set_info_mode(unsigned long)");

		if (md == vector_flags::fixed || md == vector_flags::expand)
			info_mode = md;
	}

	unsigned long get_info_mode() const
	{
		debug_handler ("vector_flags", "get_info_mode()");
		return info_mode;
	}

	//
	// set print_mode / get print_mode
	//

public:

	void set_print_mode(unsigned long md)
	{
		debug_handler ("vector_flags", "set_print_mode(unsigned long)");
		print_mode = (md != 0 && md != 1) ? print_mode : md;
	}

	unsigned long get_print_mode()
	{
		debug_handler ("vector_flags", "get_print_mode()");
		return print_mode;
	}

	//
	// assignment
	//

	vector_flags &
	operator = (const vector_flags &b)
	{
		debug_handler ("vector_flags", "operator = (const vector_flags&)");

		print_mode = b.print_mode;
		info_mode = b.info_mode;
		sort_mode = b.sort_mode;
		dyn_exp_ratio = b.dyn_exp_ratio;
		return *this;
	}

	//
	// swap function
	//

public:

	void swap(vector_flags &v)
	{
		debug_handler ("vector_flags", "swap(vector_flags&)");

		LiDIA::swap(print_mode, v.print_mode);
		LiDIA::swap(info_mode, v.info_mode);
		LiDIA::swap(dyn_exp_ratio, v.dyn_exp_ratio);
	}

};



#ifdef LIDIA_NAMESPACE
}	// end of namespace LiDIA
# undef IN_NAMESPACE_LIDIA
#endif



#endif	// LIDIA_VECTOR_FLAGS_H_GUARD_
