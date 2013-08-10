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
//	Author	: Patrick Theobald (PT)
//	Changes	: See CVS log
//
//==============================================================================================


#ifndef LIDIA_FILE_VECTOR_CC_GUARD_
#define LIDIA_FILE_VECTOR_CC_GUARD_



#ifndef LIDIA_FILE_VECTOR_H_GUARD_
# include	"LiDIA/file_vector.h"
#endif
#ifndef LIDIA_ERROR_H_GUARD_
# include	"LiDIA/error.h"
#endif
#include "LiDIA/precondition_error.h"


#ifdef LIDIA_NAMESPACE
# ifndef IN_NAMESPACE_LIDIA
namespace LiDIA {
# endif
#endif



//
// debug defines / error defines
//

extern const char *PRT;
extern const char *vector_error_msg[];

#define DV_FV LDBL_VECTOR + 10   // Debug value
#define DM_FV "file_vector"      // Debug message / Error message
#define LFV_ERROR vector_error_msg

//
// debug level
//
//   0 : writing to file in binary format
//   1 : reading from file in binary format
//   2 : appending to file in binary format
//   3 : printing to file in ascii format
//   4 : scanning form file in ascii format
//   5 : appending to file in ascii format
//

//
// writing to file in binary-format
//

template< class T >
void
file_vector< T >::write_to_file (FILE *fp) const
{
	debug_handler_l(DM_FV, "in member - function "
			"write_to_file(FILE *)", DV_FV);

	lidia_size_t i;
	fwrite (&this->length, sizeof(lidia_size_t), 1, fp);
	for (i = 0; i < this->length; i++)
		file_io_class< T >::write_to_file(this->value[i], fp);
}



//
// reading from file in binary-format
//

template< class T >
void
file_vector< T >::read_from_file (FILE *fp)
{
	debug_handler_l(DM_FV, "in member - function "
			"read_from_file(FILE *)", DV_FV + 1);

	lidia_size_t l, i;
	fread(&l, sizeof (lidia_size_t), 1, fp);
	if (l > 0) {
		this->set_capacity(l);
		this->set_size(l);
		for (i = 0; i < l && !feof(fp); i++)
			file_io_class< T >::read_from_file (this->value[i], fp);
		if (i < l)
			lidia_error_handler("void file_vector< T >::"
					    "read_from_file(FILE *fp)",
					    DM_FV, LFV_ERROR[2]);
	}
	else
		this->set_size(0);
}



//
// appending one element to file in binary-format
//
// assumes a vector to be the last element in the file;
// if there is no vector to append to, a new vector containing the
// (n+1)-st element will be appended to the file.
//



template< class T >
void
file_vector< T >::append_b (FILE *fp, lidia_size_t n) const
{
	debug_handler_l(DM_FV, "in member - function "
			"append_b(FILE *, lidia_size_t)", DV_FV + 2);

	lidia_size_t l;
	long old_fpos;

	if (n >= this->length)
		precondition_error_handler(n, "n", "n < length",
				    "void file_vector< T >::"
				    "append_b(FILE *fp, lidia_size_t n) const",
				    DM_FV, LFV_ERROR[3]);
	else {
		old_fpos = ftell(fp);
		fread(&l, sizeof(lidia_size_t), 1, fp);

		// no vector was found
		if (feof(fp)) {
			l = 1;
			fwrite(&l, sizeof(lidia_size_t), 1, fp);
			file_io_class< T >::write_to_file(this->value[n], fp);
		}

		// append to existing vector
		else {
			// increment length
			l += 1;

			// store new length
			fseek(fp, -(static_cast<long>(sizeof(lidia_size_t))), SEEK_CUR);
			fwrite(&l, sizeof(lidia_size_t), 1 , fp);

			// append the element to the end of the vector
			fseek(fp, 0, SEEK_END);
			file_io_class< T >::write_to_file(this->value[n], fp);
		}

		// reset file-pointer
		fseek(fp, old_fpos, SEEK_SET);
	}
}



//
// appending the vector to file in binary-format
//
// assumes a vector to be the last element in the file;
// if there is no vector to append to, the appending-operation
// works like a write_to_file()-operation.
//

template< class T >
void
file_vector< T >::append_b (FILE *fp) const
{
	debug_handler_l(DM_FV, "in member - function "
			"append_b(FILE *)", DV_FV + 2);

	lidia_size_t l, i;
	long old_fpos = ftell(fp);

	fread(&l, sizeof(lidia_size_t), 1, fp);

	// no vector to append to
	if (feof(fp)) {
		fwrite(&this->length, sizeof(lidia_size_t), 1, fp);
		for (i = 0; i < this->length; i++)
			file_io_class< T >::write_to_file(this->value[i], fp);
	}

	// append to the existing vector
	else {
		// increment the length
		l += this->length;

		// store the new length
		fseek(fp, -(static_cast<long>(sizeof(lidia_size_t))), SEEK_CUR);
		fwrite(&l, sizeof(lidia_size_t), 1, fp);

		// append the vector to the end of the vector in the file
		fseek(fp, 0, SEEK_END);
		for (i = 0; i < this->length; i++)
			file_io_class< T >::write_to_file(this->value[i], fp);
	}

	// reset file-pointer
	fseek(fp, old_fpos, SEEK_SET);
}



//
// printing to file in ASCII-format
//

template< class T >
void
file_vector< T >::print_to_file (FILE *fp) const
{
	debug_handler_l(DM_FV, "in member - function "
			"print_to_file(FILE *)", DV_FV + 3);

	lidia_size_t i;

	fprintf(fp, "[ ");
	for (i = 0; i < this->length; i++) {
		file_io_class< T >::print_to_file(this->value[i], fp);
		fprintf(fp, " ");
	}
	fprintf(fp, "]");
}



//
// scanning from file in ASCII-format
//

template< class T >
void
file_vector< T >::scan_from_file (FILE *fp)
{
	debug_handler_l(DM_FV, "in member - function "
			"scan_from_file(FILE*)", DV_FV + 4);

	lidia_size_t i;
	int end;
	long fpos = 0;
	char c = ' ';

	base_vector< T > v(0, vector_flags(vector_flags::expand));

	while (!feof(fp) && c != '[')
		fscanf(fp, "%c", &c);

	if (feof(fp))
		return;
	else {
		i = 0;
		end = 0;

		while (!end) {
			c = ' ';

			while (!feof(fp) && c == ' ') {
				fpos = ftell(fp);
				fscanf(fp, "%c", &c);
			}

			if (!feof(fp) && c != ']') {
				fseek(fp, fpos, SEEK_SET);
				file_io_class< T >::scan_from_file(v[i], fp);
				i++;
			}
			else
				end = 1;
		}

		if (c == ']')
			*this = v;
		else
			lidia_error_handler("void file_vector< T >::"
					    "scan_from_file(FILE *fp) ",
					    DM_FV, LFV_ERROR[5]);
	}
}



//
// appending one element to file in ASCII-format
//
// assumes a vector to be the last element in the file;
// if there is no vector to append to, a new vector containing the
// (n+1)-st element will be appended to the file.
//

template< class T >
void
file_vector< T >::append_a (FILE *fp, lidia_size_t n) const
{
	debug_handler_l(DM_FV, "in member - function "
			"append_a(FILE *, lidia_size_t)", DV_FV + 5);

	char c;
	int end;
	long fpos = 0;

	T tmp;

	// value[n] not accessible
	if (n >= this->length)
		precondition_error_handler(n, "n", "n < length",
				    "void file_vector< T >::"
				    "append_a(FILE *fp, lidia_size_t n) const",
				    DM_FV, LFV_ERROR[1]);

	long old_fpos = ftell(fp);
	fscanf(fp, "%c", &c);

	// no vector was found
	if (feof(fp)) {
		fprintf(fp, "[ ");
		file_io_class< T >::print_to_file(this->value[n], fp);
		fprintf(fp, " ]");

		// reset file-pointer
		fseek(fp, old_fpos, SEEK_SET);
		return;
	}

	// append to existing vector
	// read blanks
	while (!feof(fp) && c != '[')
		fscanf(fp, "%c", &c);

	if (feof(fp))
		lidia_error_handler("void file_vector< T >::"
				    "append_a(FILE *fp, lidia_size_t n) const",
				    DM_FV, LFV_ERROR[5]);

	// read vector
	else {
		end = 0;

		while (!end) {
			c = ' ';

			while (!feof(fp) && c == ' ') {
				fpos = ftell(fp);
				fscanf(fp, "%c", &c);
			}

			if (!feof(fp) && c != ']') {
				fseek(fp, fpos, SEEK_SET);
				file_io_class< T >::scan_from_file(tmp, fp);
			}
			else
				end = 1;
		}

		if (c == ']') {
			fseek(fp, fpos, SEEK_SET);
			file_io_class< T >::print_to_file(this->value[n], fp);
			fprintf(fp, " ]");
		}
		else
			lidia_error_handler("void file_vector< T >::"
					    "append_a(FILE *fp, lidia_size_t n) const",
					    DM_FV, LFV_ERROR[5]);
	}

	// reset file-pointer
	fseek(fp, old_fpos, SEEK_SET);
}



//
// appending the vector to file in binary-format
//
// assumes a vector to be the last element in the file;
// if there is no vector to append to, the appending-operation
// works like a print_to_file()-operation.
//

template< class T >
void
file_vector< T >::append_a (FILE *fp) const
{
	debug_handler_l(DM_FV, "in member - function :: "
			"append_a(FILE *)", DV_FV + 5);

	char c;
	lidia_size_t i;
	int end;
	long fpos = 0;

	T tmp;

	long old_fpos = ftell (fp);
	fscanf (fp, "%c", &c);

	// no vector was found
	if (feof(fp)) {
		fprintf(fp, "[ ");
		for (i = 0; i < this->length; i++) {
			file_io_class< T >::print_to_file(this->value[i], fp);
			fprintf(fp, " ");
		}
		fprintf(fp, "]");

		// reset file-pointer
		fseek(fp, old_fpos, SEEK_SET);
		return;
	}

	// append to existing vector
	// read blanks
	while (!feof(fp) && c != '[')
		fscanf (fp, "%c", &c);

	if (feof(fp))
		lidia_error_handler("void file_vector< T >::"
				    "append_a(FILE *fp) const",
				    DM_FV, LFV_ERROR[5]);
	// read vector
	else {
		end = 0;
		while (!end) {
			c = ' ';
			while (!feof(fp) && c == ' ') {
				fpos = ftell(fp);
				fscanf(fp, "%c", &c);
			}

			if (!feof(fp) && c != ']') {
				fseek(fp, fpos, SEEK_SET);
				file_io_class< T >::scan_from_file(tmp, fp);
			}
			else
				end = 1;
		}

		if (c == ']') {
			fseek(fp, fpos, SEEK_SET);
			for (i = 0; i < this->length; i++) {
				file_io_class< T >::print_to_file(this->value[i], fp);
				fprintf(fp, " ");
			}

			fprintf(fp, "]");
		}
		else
			lidia_error_handler("void file_vector< T >::"
					    "append_a(FILE *fp) const",
					    DM_FV, LFV_ERROR[5]);
	}

	// reset file-pointer
	fseek(fp, old_fpos, SEEK_SET);
}



#undef DV_FV
#undef DM_FV
#undef LFV_ERROR



#ifdef LIDIA_NAMESPACE
# ifndef IN_NAMESPACE_LIDIA
}	// end of namespace LiDIA
# endif
#endif



#endif	// LIDIA_FILE_VECTOR_CC_GUARD_
