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
//	Author	: Volker Mueller (VM), Markus Maurer (MM)
//	Changes	: See CVS log
//
//==============================================================================================



#ifdef HAVE_CONFIG_H
# include	"config.h"
#endif
#include	"LiDIA/trace_mod.h"



#ifdef LIDIA_NAMESPACE
namespace LiDIA {
#endif



std::ostream & operator << (std::ostream & o, const trace_mod & t)
{
	o <<  t.cpos << " " << t.modulus << std::flush;
	return o;
}



std::istream & operator >> (std::istream & istr, trace_mod & t)
{
	istr >> t.cpos;
	istr >> t.modulus;
	t.inv_density = static_cast<float>(t.cpos.get_size()) / static_cast<float>(t.modulus);
	return istr;
}



trace_mod::trace_mod ()
{
	modulus = 0;
	inv_density = 0;
	cpos.set_mode(EXPAND);
}



trace_mod::~trace_mod ()
{
}



trace_mod & trace_mod::operator = (const trace_mod & tm)
{
	modulus = tm.modulus;
	cpos = tm.cpos;
	inv_density = tm.inv_density;
	return *this;
}



void trace_mod::set_first_element (udigit p_l, udigit c)
{
	modulus = p_l;
	cpos.set_capacity(1);
	cpos[0] = c;
	inv_density = 1.0 / static_cast<float>(p_l);
}



void trace_mod::set_vector (udigit p_l, const base_vector< udigit > & p_cpos)
{
	modulus = p_l;
	cpos = p_cpos;
	inv_density = static_cast<float>(cpos.size()) / static_cast<float>(p_l);
}



void swap (trace_mod & t1, trace_mod & t2)
{
	trace_mod h = t1;
	t1 = t2; t2 = h;
}



#ifdef LIDIA_NAMESPACE
}	// end of namespace LiDIA
#endif
