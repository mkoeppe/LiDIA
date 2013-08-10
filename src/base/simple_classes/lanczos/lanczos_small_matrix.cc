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
#include	"LiDIA/lanczos.h"


#ifdef LIDIA_NAMESPACE
namespace LiDIA {
#endif


lanczos_small_matrix *
lanczos_small_matrix::mult_right(const lanczos_small_matrix& m) const
{
  lanczos_small_matrix *result;

  result = new lanczos_small_matrix();
  mult_right_to(m, *result);
  return result;
}



lanczos_small_matrix *
lanczos_small_matrix::mult_right_transpose(const lanczos_small_matrix& m) const
{
  lanczos_small_matrix *result;
  
  result = new lanczos_small_matrix();
  mult_right_transpose_to(m, *result);
  return result;
}



void
lanczos_small_matrix::mult_right_to(const lanczos_small_matrix& m,
				    lanczos_small_matrix& result) const
{
  size_type i, j;
  value_type s;
  
  for (i = 0; i < WordSize; i++)
    {
      s = 0;
      for (j = 0; j < WordSize; j++)
	if (get_row(i) & Bit_mask(j))
	  s = s ^ m.get_row(j);
      
      result.put_row(i, s);
    }
}



void
lanczos_small_matrix::mult_right_transpose_to(const lanczos_small_matrix& m,
					      lanczos_small_matrix& result) const
{
  size_type i, j;
  value_type s, maske;
  
  for (i = 0; i < WordSize; i++)
    {
      s = 0;
      maske = Bit_mask(i);
      
      for (j = 0; j < WordSize; j++)
	if (get_row(j) & maske)
	  s = s ^ m.get_row(j);
      
      result.put_row(i, s);
    }
}


void
lanczos_small_matrix::print() const
{

  for (size_type i = 0; i < WordSize; i++)
    bin_out(rows[i]);
  std::cout << "-------------------------------" << std::endl;
}



#ifdef LIDIA_NAMESPACE
}	// end of namespace LiDIA
#endif
