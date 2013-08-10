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
//	$Id: crossover.tbl,v 2.1 2001/05/15 12:05:30 hamdy Exp $
//
//	Author	:
//	Changes	: See CVS log
//
//==============================================================================================

//
// These are default values
// Smarter values can be computed by 'make optimize'
//

static int x_val[CROV_NUM_VALUES] = { 0, 10, 20, 40, 80, 160, 320, 640, 1280, 2560 };
static int fftmul_val[CROV_NUM_VALUES] = { 64, 64, 64, 64, 64, 64, 64, 64, 64, 64 };
static int fftdiv_val[CROV_NUM_VALUES] = { 64, 64, 64, 64, 64, 64, 64, 64, 64, 64 };
static int fftrem_val[CROV_NUM_VALUES] = { 64, 64, 64, 64, 64, 64, 64, 64, 64, 64 };
static int inv_val[CROV_NUM_VALUES] = { 32, 32, 32, 32, 32, 32, 32, 32, 32, 32 };
static int gcd_val[CROV_NUM_VALUES] = { 700, 700, 700, 700, 700, 700, 700, 700, 700, 700 };
static int xgcd_val[CROV_NUM_VALUES] = { 200, 200, 200, 200, 200, 200, 200, 200, 200, 200 };
int halfgcd_val = 16;
int log2_newton_val = 5;
