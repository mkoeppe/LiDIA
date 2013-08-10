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
#include	<fstream>


#ifdef LIDIA_NAMESPACE
namespace LiDIA {
#endif


void
lanczos_vector_block::read(result_vector_type const& res)
{
  clear();
  size_type col = 0;
  
  while (res[col].size() > 0 && (col < WordSize)) 
    {
      for (size_type i = 1; i <= res[col][0]; i++)
	{
	  size_type value = res[col][i];
	  put_row(value-1, get_row(value-1)|Bit_mask(col));
	}
      col++;
    }
}

  

lanczos_vector_block::result_vector_type
lanczos_vector_block::result() const
{
  size_type j, i, k;
  size_type eintr;
  
  result_vector_type result(WordSize + 1);
  
  for (j = 0, k = 0; j < WordSize; j++, k++) 
    {
      eintr = 0;
      for (i = 0; i < length; i++)  // count ones
	if (rows[i] & Bit_mask(j)) 
	    eintr ++;

      if (eintr == 0)  // Null solution
	{
	  k --;
	  continue;
	}
      
      result[k].resize(eintr + 1);
      result[k][0] = eintr;
      eintr = 1;
      
      for (i = 0; i < length; i++) 
	if (rows[i] & Bit_mask(j)) 
	  result[k][eintr++] = i+1;
    }
  
  return result;
}


lanczos_small_matrix *
lanczos_vector_block::mult(const lanczos_vector_block& v) const
{
  lanczos_small_matrix *result;

  result = new lanczos_small_matrix();
  mult_to(v, *result);
  return result;
}



void
lanczos_vector_block::mult_to(const lanczos_vector_block& v,
			      lanczos_small_matrix& result) const
{
  size_type i, k;
  value_type C[SIZEOF_LONG][256];  

  memset(C, 0, sizeof(C));

  for (i = 0; i < length; i++)
    {
     value_type row_val = v.get_row(i);

     for(k = 0; k < SIZEOF_LONG; k++)
       {
         C[SIZEOF_LONG-1-k][(rows[i] >> (k*8)) & 0xFF] ^= row_val;
       }
    }
  
  result.clear();
  
  i = 0;
  for (k = 0; k < SIZEOF_LONG; k++) 
    {
      for(size_type fj = 128; fj; fj >>= 1)
	{
	  for (size_type j = fj; j < 256; j += fj)
	    {
	      for (size_type m = 0; m < fj; m++)
		{
		  result.put_row(i, result.get_row(i)^C[k][j]);
		  j++;
		}
	    }
	  i++;
	}
    }
}


lanczos_vector_block *
lanczos_vector_block::mult_small (const lanczos_small_matrix& m) const
{
  lanczos_vector_block *result;
  
  result = new lanczos_vector_block(length);
  mult_small_to(m, *result);
  return result;
}



void
lanczos_vector_block::mult_small_to (const lanczos_small_matrix& m, 
				     lanczos_vector_block& result) const
{
  for (size_type i = 0; i < length; i++)
    {
      // s holds row of a * b
      value_type s = 0;
      // trying all columns
      for (size_type j = 0; j < WordSize; j++)
	if (rows[i]  &  Bit_mask(j))
	  s = s ^ m.get_row(j);
      
      result.put_row(i, s);
    }
}

void
lanczos_vector_block::print() const
{
  std::cout << "-Vectorblock--------------------";
  for (size_type i = 0; i < length; i++)
    bin_out(rows[i]);
  std::cout << "-------------------------------" << std::endl;
}




#ifdef LIDIA_NAMESPACE
}	// end of namespace LiDIA
#endif
