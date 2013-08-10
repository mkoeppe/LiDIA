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


// Pre and Postprocessing of the matrix


#ifdef HAVE_CONFIG_H
# include	"config.h"
#endif
#include	"LiDIA/lanczos.h"
#include    <cassert>


#ifdef LIDIA_NAMESPACE
namespace LiDIA {
#endif



std::auto_ptr<index_list>
preprocess::process (lanczos_sparse_matrix& matrix) const
{
	size_type correction_list_size;
	size_type count;
	index_list *zero_rows = NULL;

	//
	// preprocess Matrix
	//

#ifdef DEBUG
	printf("preprocess: rows %ld, columns %ld\n", matrix.number_of_rows(),
	       matrix.number_of_columns());
	printf("Number of singeltons: %ld\n", matrix.calculate_weight());
#endif
	long weight_counter = matrix.calculate_weight();

	while (weight_counter > 0) {

		for (size_type i = 0; i < matrix.number_of_rows(); i++) {
			long index = matrix.row_weight->get(i);
			if (index < 0) {
				// merke dir diese Spalte
				matrix.get_vector(-index-1).put_number_of_entries(0);
				// korigiere noch die Zeileneintraege
			}
		}
#ifdef DEBUG
		printf("Number of singeltons: %ld\n", matrix.calculate_weight());
#endif
		weight_counter = matrix.calculate_weight();
	}

	// removing all 0-columns

	count = 0;
	for (size_type i = 0; i < matrix.number_of_rows(); i++)
		if (matrix.row_weight->get(i) == 0)
			count++;
	zero_rows = new index_list(count);
	if (!zero_rows) {
		printf("Error: Not enough Memory\n");
		return std::auto_ptr<index_list>();
	}
#ifdef DEBUG
	printf("Killing Zero rows: %ld\n", count);

#endif

	count = 0;
	for (size_type i = 0; i < matrix.number_of_rows(); i++) {
		if (matrix.row_weight->get(i) == 0) {
			zero_rows->put(count, i);
			count++;
		}
	}

	if (count > 0)
		matrix.delete_rows(*zero_rows);

#ifdef DEBUG
	weight_counter = matrix.calculate_weight();
	for (size_type i = 0; i < matrix.number_of_rows(); i++)
		if (matrix.row_weight->get(i) == 0)
			printf("Warning after killing zero rows:\n  Row %ld is zero\n", i);
#endif


	// counting all 0-vectos in matrix B

	correction_list_size = 0;

	for (size_type i = 0; i < matrix.number_of_columns(); i++) {
		if (matrix.get_vector(i).get_number_of_entries() == 0)
			correction_list_size++;
	}
#ifdef DEBUG
	printf("Number of zero columns: %ld\n", correction_list_size);
#endif

	std::auto_ptr<index_list> correction_list;
	if (correction_list_size > 0) {
		correction_list.reset(new index_list (correction_list_size));
        }

	long counter = 0;

	if (correction_list.get()) {

#ifdef DEBUG
		printf("Packing matrix...\n");
#endif

		for (size_type i = 0; i < matrix.number_of_columns(); i++) {


			if (matrix.get_vector(i).get_number_of_entries() == 0) {

				assert( counter < correction_list_size );

				correction_list->put(counter, i);
				counter++;
				matrix.delete_vector(i);
			}
			else if (counter > 0) {
				//kopieren des Vektors
				matrix.put_vector(i-counter, matrix.get_vector(i));

			}
		}

		// Loesche die letzten Vektoren
		for (size_type i = 0; i < counter; i++)
			matrix.delete_vector(matrix.number_of_columns()-i-1);
		// Setze Columns neu
		matrix.put_columns(matrix.number_of_columns()-counter);
#ifdef DEBUG
		size_type temp_count = 0;
		for (size_type i = 0; i < matrix.number_of_columns(); i++) {
			if (matrix.get_vector(i).get_number_of_entries() == 0)
				temp_count++;
		}


		printf("After preprocessing: rows %ld, columns %ld, zero-columns %ld\n",
		       matrix.number_of_rows(), matrix.number_of_columns(), temp_count);
		weight_counter = matrix.calculate_weight();
		for (size_type i = 0; i < matrix.number_of_rows(); i++) {
			if (matrix.row_weight->get(i) == 0)
				printf("Warning::Row %ld is a zero-row\n", i);
			if (matrix.row_weight->get(i) < 0)
				printf("Warning::Row %ld is a singelton at column %ld\n",
				       i, -matrix.row_weight->get(i)-1);
		}
#endif
	}
	if (zero_rows)
		delete(zero_rows);

        // postprocess cannot handle a null list from preprocess. So either 
        // arguments and content of postprocess are changed, or return
        // an empty list.

        if (!correction_list.get())
        {
	    correction_list.reset(new index_list (0));
        }
	return correction_list;
}



std::auto_ptr<lanczos_vector_block>
postprocess::process(const lanczos_vector_block& vector,
		     const index_list& list) const
{
	std::auto_ptr<lanczos_vector_block> result;

	//
	// creating a vector
	//

	size_type count = 0;

	result.reset(new lanczos_vector_block(vector.get_length() +
                                              list.number()));

	for (size_type i = 0; i < vector.get_length(); ++i) {

		while ((count < list.number()) &&
		       (i+count >= list.get(count)))
			count++;

		result->put_row(i+count, vector.get_row(i));
	}
	return result;
}



#ifdef LIDIA_NAMESPACE
}	// end of namespace LiDIA
#endif
