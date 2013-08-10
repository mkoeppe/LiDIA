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
//	$Id$
//
//	Author	: Patrick Theobald (PT)
//	Changes	: See CVS log
//
//==============================================================================================


#ifndef LIDIA_MATRIX_FLAGS_H_GUARD_
#define LIDIA_MATRIX_FLAGS_H_GUARD_



#ifndef LIDIA_ARITH_INL_GUARD_
# include	"LiDIA/arith.inl"
#endif



#ifdef LIDIA_NAMESPACE
namespace LiDIA {
# define IN_NAMESPACE_LIDIA
#endif



class matrix_flags
{

	unsigned long print_mode;
	unsigned long storage_mode;

public: // should be removed in LiDIA 1.4

	unsigned long structure_mode;
	unsigned long info_mode;
	unsigned long lattice_mode;

public:

	//
	// PRINT MODE SETTINGS
	//

	enum {
		beauty_mode = 0,
		lidia_mode = 1,
		gp_mode = 2,
		maple_mode = 3,
		mathematica_mode = 4,
		kash_mode = 5,
		latex_mode = 6,
		magma_mode = 7,

		default_print_mode = beauty_mode
	};

	//
	// STORAGE MODE SETTINGS
	//

	enum {
		representation = 3,
		dense_representation = 1,
		sparse_representation = 2,
		mixed_representation = 3,

		orientation = 4,
		row_oriented = 0,
		column_oriented = 4,

		default_storage_mode = dense_representation | row_oriented
	};

	//
	// STRUCTURE MODE SETTINGS
	//

	enum {
		diag = 1,
		upper_diag = 2,
		lower_diag = 4,
		upper_tria = 8,
		lower_tria = 16,
		columns_linind = 32,
		rows_linind = 64,

		default_structure_mode = 0
	};

	//
	// INFO MODE SETTINGS (Position der Diagonalen)
	//

	enum {
		diag_up = 1,
		diag_right = 4,
		diag_ld_to_ru = 8,

		default_info_mode = 0
	};

	//
	// LATTICE MODE SETTINGS
	//

	enum {
		default_lattice_mode = 0
	};

	//
	// constructor
	//

	matrix_flags()
	{
		info_mode = default_info_mode;
		print_mode = default_print_mode;
		structure_mode = default_structure_mode;
		storage_mode = default_storage_mode;
		lattice_mode = default_lattice_mode;
	}

	//
	// destructor
	//

	~matrix_flags() {}

	//
	// assignment
	//

	matrix_flags & operator = (const matrix_flags &v)
	{
		info_mode = v.info_mode;
		print_mode = v.print_mode;
		structure_mode = v.structure_mode;
		storage_mode = v.storage_mode;
		lattice_mode = v.lattice_mode;

		return *this;
	}

	//
	// swap
	//

	friend void swap(matrix_flags &, matrix_flags &);

	//
	// print_mode
	//

public:

	void set_print_mode(unsigned long art)
	{
		print_mode = art;
	}

	unsigned long get_print_mode() const
	{
		return print_mode;
	}

	//
	// set storage_mode / get storage_mode
	//

public:

	void set_storage_mode(unsigned long art)
	{
		storage_mode = art;
	}

	unsigned long get_storage_mode() const
	{
		return storage_mode;
	}

	//
	// set structure_mode / get structure_mode
	//

public:

	void set_structure_mode(unsigned long art)
	{
		structure_mode = art;
	}

	unsigned long get_structure_mode() const
	{
		return structure_mode;
	}

	//
	// set info_mode / get info_mode
	//

public:

	void set_info_mode(unsigned long art)
	{
		info_mode = art;
	}

	unsigned long get_info_mode() const
	{
		return info_mode;
	}

	//
	// set lattice_mode / get lattice_mode
	//

public:

	void set_lattice_mode(unsigned long art)
	{
		lattice_mode = art;
	}

	unsigned long get_lattice_mode() const
	{
		return lattice_mode;
	}


	//
	// additional functions for changing storage_mode value
	//

	//
	// get_representation
	//

public:

	unsigned long get_representation() const
	{
		return (storage_mode & representation);
	}

	void set_representation(unsigned long art)
	{
		storage_mode = storage_mode - (storage_mode & representation) + art;
	}

	//
	// get_orientation // set_orientation
	//

public:

	unsigned long get_orientation() const
	{
		return (storage_mode & orientation);
	}

	void set_orientation(unsigned long art)
	{
		storage_mode = storage_mode - (storage_mode & orientation) + art;
	}
};



inline void
swap(matrix_flags &w, matrix_flags &v)
{
	LiDIA::swap(w.info_mode, v.info_mode);
	LiDIA::swap(w.print_mode, v.print_mode);
	LiDIA::swap(w.structure_mode, v.structure_mode);
	LiDIA::swap(w.storage_mode, v.storage_mode);
	LiDIA::swap(w.lattice_mode, v.lattice_mode);
}



#ifdef LIDIA_NAMESPACE
}	// end of namespace LiDIA
# undef IN_NAMESPACE_LIDIA
#endif



#endif	// LIDIA_MATRIX_FLAGS_H_GUARD_
