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
//	Author	: Safuat Hamdy (SH)
//	Changes	: See CVS log
//
//==============================================================================================


#include	"LiDIA/bigint.h"
#include	"LiDIA/xdouble.h"


#define E_TYPE bigint
#define A_TYPE xdouble

#define VAR Normal
#define MODULES basis_modules
#define VECTOR vector_op_SP


#include	"LiDIA/instantiate/lattices.cc"
