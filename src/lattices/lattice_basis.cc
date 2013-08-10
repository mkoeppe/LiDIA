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
//	Author	: Werner Backes (WB), Thorsten Lauer (TL)
//	Changes	: See CVS log
//
//==============================================================================================


#ifndef HEADBANGER
#ifdef HAVE_CONFIG_H
# include	"config.h"
#endif
#include	"LiDIA/lattice_gensys.h"
#include	"LiDIA/lattice_basis.h"
#endif



#ifdef LIDIA_NAMESPACE
namespace LiDIA {
#endif



// MM, The include of gensys is a hack, because the cpp of CC
//     on SunOs4.1.4 has a maximum include nesting depth of 12.

//
// Algorithms
//

void lattice_basis::lll()
{
	debug_handler("lattice_basis", "lll()");
	CT.start_timer();
	if (comp_mode == bigfloat_mode)
		lbf->lll_var();
	else
		lbi->lll();
	CT.stop_timer();
}

void lattice_basis::lll(const lattice_basis& B)
{
	debug_handler("lattice_basis", "lll(B)");
	assign(B);
	CT.start_timer();
	if (comp_mode == bigfloat_mode)
		lbf->lll_var();
	else
		lbi->lll();
	CT.stop_timer();
}

void lattice_basis::lll_trans(lattice_basis& Tr)
{
	debug_handler("lattice_basis", "lll_trans(Tr)");
	Tr.set_computing_mode(comp_mode);
	CT.start_timer();
	if (comp_mode == bigfloat_mode)
		lbf->lll_trans_var(*Tr.lbf);
	else
		lbi->lll_trans(*Tr.lbi);
	CT.stop_timer();
}

//
// Conversion
//

void lattice_basis::extract_basis(const lattice_gensys& L, lidia_size_t& rank)
{
	debug_handler("lattice_basis", "extract_basis(L, rank)");
	set_computing_mode(L.comp_mode);
	if (comp_mode == bigint_mode)
		lbi->extract_basis(*L.lbi, rank);
	else
		lbf->extract_basis(*L.lbf, rank);
}

bool lattice_basis::make_basis(const lattice_gensys& L, lidia_size_t& rank)
{
	debug_handler("lattice_basis", "make_basis(L, rank)");
	set_computing_mode(L.comp_mode);
	if (comp_mode == bigint_mode)
		return(lbi->make_basis(*L.lbi, rank));
	else
		return(lbf->make_basis(*L.lbf, rank));
}

//
// lll - checking
//
bool lattice_basis::lll_check(double y)
{
	debug_handler("lattice_basis", "lll_check(y)");
	if (comp_mode == bigint_mode)
		return(lbi->lll_check(y));
	else
		return(lbf->lll_check(y));
}

bool lattice_basis::lll_check(sdigit a, sdigit b)
{
	debug_handler("lattice_basis", "lll_check(a, b)");
	if (comp_mode == bigint_mode)
		return(lbi->lll_check(a, b));
	else
		return(lbf->lll_check(a, b));
}

double lattice_basis::lll_check_search()
{
	debug_handler("lattice_basis", "lll_check_search()");
	if (comp_mode == bigint_mode)
		return(lbi->lll_check_search());
	else
		return(lbf->lll_check_search());
}

void lattice_basis::lll_check_search(sdigit& a, sdigit& b)
{
	debug_handler("lattice_basis", "lll_check_search(a, b)");
	if (comp_mode == bigint_mode)
		lbi->lll_check_search(a, b);
	else
		lbf->lll_check_search(a, b);
}

void lattice_basis::gram_schmidt_orth(lattice_basis& Mue, lattice_basis& Gso)
{
	debug_handler("lattice_basis", "gram_schmidt_orth(Mue, Gso)");
	Mue.set_computing_mode(bigfloat_mode);
	Gso.set_computing_mode(bigfloat_mode);
	if (comp_mode == bigint_mode)
		lbi->gram_schmidt_orth(*Mue.lbf, *Gso.lbf);
	else
		lbf->gram_schmidt_orth(*Mue.lbf, *Gso.lbf);
}

void lll(math_matrix< bigint > & A, const math_matrix< bigint > & B)
{
	debug_handler("lattice_basis", "friend lll(A, B)");
	bigint_lattice_basis temp(B);
	temp.lll();
	A.assign(temp);
}

void lll_trans(math_matrix< bigint > & T, math_matrix< bigint > & A)
{
	debug_handler("lattice_basis", "friend lll_trans(T, A)");
	bigint_lattice_basis temp(A);
	temp.lll_trans(T);
	A.assign(temp);
}

void lll(math_matrix< bigfloat > & A, const math_matrix< bigfloat > & B)
{
	debug_handler("lattice_basis", "friend lll(A, B)");
	bigfloat_lattice_basis temp(B);
	temp.lll_var();
	A.assign(temp);
}

void lll_trans(math_matrix< bigint > & T, math_matrix< bigfloat > & A)
{
	debug_handler("lattice_basis", "friend lll_trans(T, A)");
	bigfloat_lattice_basis temp(A);
	temp.lll_trans_var(T);
	A.assign(temp);
}

void lll_trans(math_matrix< bigfloat > & T, math_matrix< bigfloat > & A)
{
	debug_handler("lattice_basis", "friend lll_trans(T, A)");
	bigfloat_lattice_basis temp(A);
	temp.lll_trans_var(T);
	A.assign(temp);
}



#ifdef LIDIA_NAMESPACE
}	// end of namespace LiDIA
#endif
