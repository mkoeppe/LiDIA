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
#include	"config.h"
#endif
#include	"LiDIA/lanczos.h"
#include        "LiDIA/random_generator.h"


#ifdef LIDIA_NAMESPACE
namespace LiDIA {
#endif



struct entrie_list {
	size_t number;
	unsigned long* entries;
	struct entrie_list *next;
};


lanczos_sparse_matrix::lanczos_sparse_matrix(const size_type row,
                                             const size_type col)
{
  rows = row;
  columns = col;
  entries = new lanczos_sparse_vector [columns];
  row_weight = new index_list(rows);
  row_weight->clear();
}

lanczos_sparse_matrix::lanczos_sparse_matrix(const char *filename)
{
  // size_type total = 0;
  value_type  M = 0, N = 0, sn, enr, dummy, max_val = 0, wert, index;
  size_type i, j;
  FILE *fp;
  
  struct entrie_list *entrie_list;
  struct entrie_list *next_entrie;
  
  entrie_list = new (struct entrie_list);
  entrie_list->entries = NULL;
  entrie_list->next = NULL;
  next_entrie = entrie_list;

  if ((fp = fopen(filename, "r")) == NULL) 
    lidia_error_handler("lanczos","Can't open file!");
  
  // now the reading starts
  
  while (fscanf(fp, "%ld %ld %ld ", &sn, &enr, &dummy) != EOF) 
    {
      next_entrie->next = new (struct entrie_list);
      next_entrie = next_entrie->next;
      next_entrie->next = NULL;
      if (enr & 1)
      {
	lidia_error_handler("lanczos","lanczos_sparse_matrix::Format of"
			    " input matrix wrong");
      }
      else
	enr >>= 1;
      next_entrie->entries = new unsigned long [enr];
      next_entrie->number = enr;
      for (j = 0; j < enr; ++j) 
	{
	  fscanf(fp, "%ld %ld", &wert, &index);
	  if (wert & 1) 
	    {
	      // ++total;
	      if (index > max_val)
		max_val = index;
	      next_entrie->entries[j] = index;
	    }
          else
            {
              next_entrie->entries[j] = 0;
	    }
	}
      fscanf(fp, "%ld", &dummy);
      
      M = sn;
      N = max_val+1;
    }
  fclose(fp);

  rows = N;
  columns = M;

  entries = new lanczos_sparse_vector [columns];

  next_entrie = entrie_list->next;
  for (i = 0; i < columns; ++i) 
    {
      entries[i].set_entries(next_entrie->entries, next_entrie->number, 
			     rows);
      next_entrie = next_entrie->next;
    }
  row_weight = new index_list(rows);

  // killing the help construction
  
  while (entrie_list->next) 
    {
      next_entrie = entrie_list->next;
      entrie_list->next = next_entrie->next;
      delete(next_entrie);
    }
  delete(entrie_list);
}



void
lanczos_sparse_matrix::write(const char *filename) const
{
  FILE *fp;
  size_type number, sn, i;
  
  if ((fp = fopen(filename, "w")) == NULL) {
    lidia_error_handler ("lanczos_sparse_matrix::write(char *filename)",
			 "Can't open file for writing.");
  }
  else 
    {
      for (sn = 0; sn < number_of_columns(); sn++) {
	number = entries[sn].get_number_of_entries();
	fprintf(fp, "%ld %ld %d ", sn+1, number*2, 0);
	for (i = 0; i < number; i++)
	  fprintf(fp, "1 %ld ", entries[sn].get_entry(i));
	fprintf(fp, "0\n");
      }
      fclose(fp);
    }
}


lanczos_sparse_matrix::~lanczos_sparse_matrix()
{
  delete[] entries;
  delete row_weight;
}


void
lanczos_sparse_matrix::delete_vector(const size_type pos)
{
  if (pos >= columns)
    lidia_error_handler("lanczos","delete_vector::pos >= columns");
  //  entries[pos].~lanczos_sparse_vector();
  entries[pos] = lanczos_sparse_vector(); // change with empty vector
}

void
lanczos_sparse_matrix::delete_rows(const index_list& row_list)
{

	size_type index, count, count_max;
	size_type i;
	size_type j;

#ifdef DEBUG
	assert(row_list != NULL);
#endif

	for (j = 0; j < columns; j++) {
	  count = 0;
	  for (i = 0; i < entries[j].get_number_of_entries(); i++) {
	    index = entries[j].get_entry(i);
	    if (count < row_list.number())
	      while ((count < row_list.number()) &&
		     (index > row_list.get(count)))
		count++;
	    if (count < row_list.number())
	      count_max = count;
	    else
	      count_max = count-1;
	    if (index == row_list.get(count_max))
	      {
		lidia_error_handler("lanczos_sparse_matrix",
				    "deleting none zero");
#ifdef DEBUG
		assert (index != row_list.get(count_max));
#endif
	      }
	    // only 0-rows are allowed 
	    entries[j].put_entry(i, index - count);
	  }
	  entries[j].length = entries[j].length - row_list.number();
	  
	  for (i = 0; i < entries[j].get_number_of_entries(); i++) 
	    {
	      if (entries[j].get_entry(i) >= entries[j].length)
		lidia_error_handler("lanczos_sparse_matrix",
				    "Wrong entry in sparse_vector");
	    }
	}
	rows = rows - row_list.number();
}



void
lanczos_sparse_matrix::print() const
{
  for (size_type i = 0; i < columns; i++)
    {
      printf("%ld: ", i);
      entries[i].print();
  }
}



void
lanczos_sparse_matrix::fill_random(const size_type number, const long part, const int ratio)
{

  size_type i;
  
  for (i = 0; i < columns; i++) {
    if (number < rows)
      put_vector(i, lanczos_sparse_vector(rows, number).fill_random(part, ratio));
  }
  
#if 0
  lanczos_sparse_vector vector = lanczos_sparse_vector(rows, rows/10);
  
  for (i = 0; i < columns; i++) {
    vector = new lanczos_sparse_vector(rows, rows/10);
    vector->fill_random(seed+i);
    put_vector(i, *vector);
  }
#endif
}



void
lanczos_sparse_matrix::mult_vectorblock_to(const lanczos_vector_block& vector_block,
					   lanczos_vector_block& result) const
{
	//lanczos_sparse_vector *pdw; // Spalte
	lanczos_vector_block dummy1_vec(rows);

#ifdef DEBUG
	assert (vector_block.get_length() == columns);
#endif

	// bilde zuerst B * a  gehe alle Eintraege durch
	for(size_type i=0; i < columns; i++) {
		size_type h = get_vector(i).get_number_of_entries();
		// Wieviel Eintraeg gibt es in der i-ten Spalte ?
		// solange Spalteneintraege gehe alle Spalten einer Zeile durch
		for(size_type j=0; j < h; j++) {
		  size_type g = get_vector(i).get_entry(j);
 		  dummy1_vec.put_row(g, vector_block.get_row(i) ^ dummy1_vec.get_row(g));
		}
	}


	// dummy1_vec = B * a nun folgt Anwendung von 0B_(transponiert)gehe alle Eintraege durch
	for(size_type i=0;i < columns; i++) {
		value_type s = 0;
		size_type h = get_vector(i).get_number_of_entries();

		// solange Spalteneintraege gehe alle Spalten einer Zeile durch
		for(size_type j=0; j < h; j++) {
			s = s ^ (dummy1_vec.get_row(get_vector(i).get_entry(j)));
		}
		result.put_row(i, s);
	}
}



lanczos_vector_block *
lanczos_sparse_matrix::mult_vectorblock(const lanczos_vector_block& vector_block) const
{
	lanczos_vector_block *result;

	result = new lanczos_vector_block(columns);
	result->clear();
	mult_vectorblock_to(vector_block, *result);

	return result;
}



long
lanczos_sparse_matrix::calculate_weight() const
{
	long counter = 0; // counts the 1 weights

	row_weight->clear();

	for (size_type i = 0; i < columns; i++)
		for (size_type j = 0; j < entries[i].get_number_of_entries(); j++) {
			size_type index = entries[i].get_entry(j);

			if (row_weight->get(index) == 0) {
				row_weight->put(index, -static_cast<long>(i)-1);
				counter++;
			}
			else {
				if (row_weight->get(index) < 0) {
					row_weight->put(index, 2);
					counter--;
				}
				else
					row_weight->put(index, row_weight->get(index)+1);
			}
		}
	return counter;
}



#ifdef LIDIA_NAMESPACE
}	// end of namespace LiDIA
#endif
