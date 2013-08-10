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
//	$Id: bigfloat_config.h,v 2.4 2002/06/24 09:41:23 lidiaadm Exp $
//
//	Author	: Thomas Papanikolaou (TP)
//	Changes	: See CVS log
//
//==============================================================================================


#ifndef LIDIA_BIGFLOAT_CONFIG_H_GUARD_
#define LIDIA_BIGFLOAT_CONFIG_H_GUARD_



#if gcos
#define bits_per_char   9
#else
#define bits_per_char   8
#endif

#define bits_per_(type) (bits_per_char * static_cast<int>(sizeof(type)))
#define bits_per_short  bits_per_(short)
#define bits_per_int    bits_per_(int)
#define bits_per_double bits_per_(double)

#ifndef base_digit
#define base_digit unsigned long
#endif

#define base_digit_radix             (bigint::radix())
#define bits_per_base_digit          (bigint::bits_per_digit())
#define bits_per_base_digit_minus_1  (bits_per_base_digit - 1)

#define high_bit_short (1 << (bits_per_(short) - 1))
#define high_bit_int   (1 << (bits_per_(int) - 1))
#define high_bit_long  (1L << (bits_per_(long) - 1))

#define max_short        (~high_bit_short)
#define max_int          (~high_bit_int)
#define max_long         (~high_bit_long)

#define rounding_bits   (3 * bits_per_base_digit)
#define digits_per_double sizeof(double)/sizeof(base_digit)

#define max_base_digit       ((1UL) << bits_per_base_digit_minus_1)
#define max_base_digit_2     (2.0*static_cast<double>(max_base_digit))

#if SIZEOF_LONG == 4
# define log_base_digit_bits  5
# define log_max_pow_10       9
# define max_pow_10           1000000000UL
#elif SIZEOF_LONG == 8
# define log_base_digit_bits  6
# define log_max_pow_10       18
# define max_pow_10           1000000000000000000UL
#else
# define log_base_digit_bits  static_cast<base_digit>(std::log(static_cast<double>(bits_per_base_digit))/std::log(2.0))
# define log_max_pow_10       static_cast<base_digit>(bits_per_base_digit*std::log(2.0)/std::log(10.0))
# define max_pow_10           static_cast<base_digit>(std::exp(log_max_pow_10*std::log(10.0))+0.5)
#endif

// Constants

#define L2B10           0.30102999566398119521	// log(2)/log(10)
#define LOG2            0.69314718055994530942	// log(2)
#define LOG10           2.30258509299404568401	// log(10)
#define INVLOG2         1.44269504088896340735	// 1/log(2)
#define INVLOG10        0.43429448190325182765	// 1/log(10)


#endif	// LIDIA_BIGFLOAT_CONFIG_H_GUARD_
