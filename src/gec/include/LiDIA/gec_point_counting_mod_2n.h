// -*- C++ -*-
//==============================================================================================
//
//	This file is part of LiDIA --- a library for computational number theory
//
//	Copyright (c) 1994--2002 the LiDIA Group.  All rights reserved.
//
//	See http://www.informatik.tu-darmstadt.de/TI/LiDIA/
//
//----------------------------------------------------------------------------------------------
//
//
//  File    : gec_point_counting_mod_2n.h
//  Author  : Harald Baier (HB)
//	Changes	: See CVS log
//
//==============================================================================================

#ifndef LIDIA_GEC_POINT_COUNTING_MOD_2N_H_GUARD_
#define LIDIA_GEC_POINT_COUNTING_MOD_2N_H_GUARD_


#ifndef LIDIA_ECO_GF2N_H_GUARD_
# include "LiDIA/eco_gf2n.h"
#endif
#ifndef LIDIA_GEC_H_GUARD_
# include "LiDIA/gec.h"
#endif
#ifndef LIDIA_TRACE_LIST_H_GUARD_
# include "LiDIA/trace_list.h"
#endif


#ifdef LIDIA_NAMESPACE
namespace LiDIA {
# define IN_NAMESPACE_LIDIA
#endif


class 
gec_point_counting_mod_2n :
public gec
{
 //
 // degree 
 //
 private: unsigned int degree;

 //
 // Attribute to store number of curve initializations
 //
 private: int tries;

 //
 // Attributes to store information on failures
 //
 private: int both_abort;
 private: int no_abort;

 /*************************************************
 *
 * Class methods 
 *
 **************************************************/

 //
 // constructor and destructor
 //
 public: gec_point_counting_mod_2n();

 public: ~gec_point_counting_mod_2n();

 //
 // Mutators
 //
 public: void set_degree( const unsigned int );

 //
 // High level functions
 // 
 public: void generate();
 public: void get_twist_coeff( gf_element & , gf_element & );

};


#ifdef LIDIA_NAMESPACE
}	// end of namespace LiDIA
# undef IN_NAMESPACE_LIDIA
#endif

#endif	// LIDIA_GEC_POINT_COUNTING_MOD_2N_H_GUARD_
