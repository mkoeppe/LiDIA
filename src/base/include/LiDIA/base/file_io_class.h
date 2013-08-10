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
//	$Id: file_io_class.h,v 2.4 2002/06/24 09:41:36 lidiaadm Exp $
//
//	Author	: Frank Lehmann (FL), Markus Maurer (MM),
//                Thomas Papanikolaou (TP), Patrick Theobald (PT)
//	Changes	: See CVS log
//
//==============================================================================================


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
// *
// *    File        :  file_io_class < T >
// *
// *    Description :  -) definition of template class 'file_io_class < T > '
// *			  to provide functions to write or read elements
// *			  of any type T to or from files, respectively
// *
// *		       -) to use this class, a none built-in type T
// *			  must support the functions
// *			       - print_to_file()
// *			       - scan_from_file()
// *			       - write_to_file()
// *			       - read_from_file()
// *
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *


#ifndef LIDIA_FILE_IO_CLASS_H_GUARD_
#define LIDIA_FILE_IO_CLASS_H_GUARD_



#include	<cstdio>



#ifdef LIDIA_NAMESPACE
namespace LiDIA {
# define IN_NAMESPACE_LIDIA
#endif



template< class T > class  file_io_class
{
public:

	static void print_to_file(T & t , FILE * fp)
	{
		t.print_to_file(fp);
	}

	static void scan_from_file (T & t , FILE * fp)
	{
		t.scan_from_file (fp);
	}

	static void write_to_file (T & t , FILE * fp)
	{
		t.write_to_file (fp);
	}

	static void read_from_file (T & t , FILE * fp)
	{
		t.read_from_file (fp);
	}
};



template<>
class file_io_class< char >
{
public:

	static void print_to_file (char & t , FILE * fp)
	{
		fprintf (fp , "%c" , t);
	}

	static void scan_from_file (char & t , FILE * fp)
	{
		fscanf (fp , "%c" , &t);
	}

	static void write_to_file (char & t , FILE * fp)
	{
		fwrite (&t , sizeof (char) , 1 , fp);
	}

	static void read_from_file (char & t , FILE * fp)
	{
		fread (&t , sizeof (char) , 1 , fp);
	}
};



template<>
class file_io_class< unsigned char >
{
public:

	static void print_to_file (unsigned char & t , FILE * fp)
	{
		fprintf (fp , "%c" , t);
	}

	static void scan_from_file (unsigned char & t , FILE * fp)
	{
		fscanf (fp , "%c" , &t);
	}

	static void write_to_file (unsigned char & t , FILE * fp)
	{
		fwrite (&t , sizeof (unsigned char) , 1 , fp);
	}

	static void read_from_file (unsigned char & t , FILE * fp)
	{
		fread (&t , sizeof (unsigned char) , 1 , fp);
	}
};



template<>
class file_io_class< double >
{
public:

	static void print_to_file (double & t , FILE * fp)
	{
		fprintf (fp , "%f" , t);
	}

	static void scan_from_file (double & t , FILE * fp)
	{
		fscanf (fp , "%lf" , &t);
	}

	static void write_to_file (double & t , FILE * fp)
	{
		fwrite (&t , sizeof (double) , 1 , fp);
	}

	static void read_from_file (double & t , FILE * fp)
	{
		fread (&t , sizeof (double) , 1 , fp);
	}
};



template<>
class file_io_class< float >
{
public:

	static void print_to_file (float & t , FILE * fp)
	{
		fprintf (fp , "%f" , t);
	}

	static void scan_from_file (float & t , FILE * fp)
	{
		fscanf (fp , "%f" , &t);
	}

	static void write_to_file (float & t , FILE * fp)
	{
		fwrite (&t , sizeof (float) , 1 , fp);
	}

	static void read_from_file (float & t , FILE * fp)
	{
		fread (&t , sizeof (float) , 1 , fp);
	}
};



template<>
class file_io_class< int >
{
public:

	static void print_to_file (int & t , FILE * fp)
	{
		fprintf (fp , "%d" , t);
	}

	static void scan_from_file (int & t , FILE * fp)
	{
		fscanf (fp , "%d" , &t);
	}

	static void write_to_file (int & t , FILE * fp)
	{
		fwrite (&t , sizeof (int) , 1 , fp);
	}

	static void read_from_file (int & t , FILE * fp)
	{
		fread (&t , sizeof (int) , 1 , fp);
	}
};



template<>
class file_io_class< unsigned int >
{
public:

	static void print_to_file (unsigned int & t , FILE * fp)
	{
		fprintf (fp , "%u" , t);
	}

	static void scan_from_file (unsigned int & t , FILE * fp)
	{
		fscanf (fp , "%u" , &t);
	}

	static void write_to_file (unsigned int & t , FILE * fp)
	{
		fwrite (&t , sizeof (unsigned int) , 1 , fp);
	}

	static void read_from_file (unsigned int & t , FILE * fp)
	{
		fread (&t , sizeof (unsigned int) , 1 , fp);
	}
};



template<>
class file_io_class< long >
{
public:

	static void print_to_file (long & t , FILE * fp)
	{
		fprintf (fp , "%ld" , t);
	}

	static void scan_from_file (long & t , FILE * fp)
	{
		fscanf (fp , "%ld" , &t);
	}

	static void write_to_file (long & t , FILE * fp)
	{
		fwrite (&t , sizeof (long) , 1 , fp);
	}

	static void read_from_file (long & t , FILE * fp)
	{
		fread (&t , sizeof (long) , 1 , fp);
	}
};



template<>
class file_io_class< unsigned long >
{
public:

	static void print_to_file (unsigned long & t , FILE * fp)
	{
		fprintf (fp , "%lu" , t);
	}

	static void scan_from_file (unsigned long & t , FILE * fp)
	{
		fscanf (fp , "%lu" , &t);
	}

	static void write_to_file (unsigned long & t , FILE * fp)
	{
		fwrite (&t , sizeof (unsigned long) , 1 , fp);
	}

	static void read_from_file (unsigned long & t , FILE * fp)
	{
		fread (&t , sizeof (unsigned long) , 1 , fp);
	}
};



#ifdef LIDIA_NAMESPACE
}	// end of namespace LiDIA
# undef IN_NAMESPACE_LIDIA
#endif



#endif	// LIDIA_FILE_IO_CLASS_H_GUARD_
