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
#include	"LiDIA/random_generator.h"
#include	<cstring>


#ifdef LIDIA_NAMESPACE
namespace LiDIA {
#endif


lanczos_sparse_vector::lanczos_sparse_vector()
{
  length = 0;
  number_of_entries = 0;
  entries = NULL;
}
  
lanczos_sparse_vector::
lanczos_sparse_vector(const lanczos_sparse_vector &vector)
{
  length = vector.get_length();
  number_of_entries = vector.get_number_of_entries();
  entries = new value_type [number_of_entries];
  memcpy (entries, vector.entries, number_of_entries * sizeof(value_type));
}



lanczos_sparse_vector::lanczos_sparse_vector(const size_type len,
					     const size_type number)
{
  length = len;
  number_of_entries = number;
  if( number ) {
    entries = new value_type [number];
    memset(entries, 0, number * sizeof(value_type));
  }
  else entries = NULL;
}

lanczos_sparse_vector::~lanczos_sparse_vector()
{
    delete[] entries;
}


void
lanczos_sparse_vector::print() const
{
  std::cout << "Sparse Vector ";
  std::cout << "(Length: " << length;
  std::cout << ", Non Zero Entries: " << number_of_entries;
  std::cout << ", Indexlist: ";
  
  for (size_type i = 0; i < number_of_entries; i++)
    std::cout << i << ": " << get_entry(i) << std::endl;
}



lanczos_sparse_vector &
lanczos_sparse_vector::fill_random(long part, int ratio)
{
  value_type index;
  size_type  count;
  random_generator rg;
  size_type number_dense;
  size_type dense_bound;
  size_type i, j;
  
  number_dense = (number_of_entries * ratio) / 100;
  dense_bound = (length * part)/100;
  
  for (i = 0; i < number_dense; i++) 
    {
      rg >> index;
      index %= dense_bound-1;
      count = 0;
      while ((i > count) &&
	     (index < get_entry(i-count-1)))
	count++;
      if (index != get_entry(i-count))
	for (j = i; j > i-count; j--)
	  put_entry(j, get_entry(j-1));
      put_entry(i-count, index);
    }
  
  for (i = number_dense; i < number_of_entries; i++)
    {
      rg >> index;
      index = index % (length-dense_bound-1);
      index = index + dense_bound;
      count = 0;
      while ((i > count) &&
	     (index < get_entry(i-count-1)))
	count++;
      if (index != get_entry(i-count))
	for (j = i; j > i-count; j--)
	  put_entry(j, get_entry(j-1));
      put_entry(i-count, index);
    }
  
  return *this;
}



#ifdef LIDIA_NAMESPACE
}	// end of namespace LiDIA
#endif
