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


#ifndef LIDIA_LIDIA_DEFINES_H_GUARD_
#define LIDIA_LIDIA_DEFINES_H_GUARD_


#ifndef LIDIA_MATRIX_FLAGS_H_GUARD_
# include	"LiDIA/matrix_flags.h"
#endif
#ifndef LIDIA_VECTOR_FLAGS_H_GUARD_
# include	"LiDIA/vector_flags.h"
#endif



//
// PRINT MODE SETTINGS
//

#define BEAUTY_MODE             matrix_flags::beauty_mode
#define LIDIA_MODE              matrix_flags::lidia_mode
#define GP_MODE                 matrix_flags::gp_mode
#define MAPLE_MODE              matrix_flags::maple_mode
#define MATHEMATICA_MODE        matrix_flags::mathematica_mode
#define KASH_MODE               matrix_flags::kash_mode
#define LATEX_MODE              matrix_flags::latex_mode
#define MAGMA_MODE              matrix_flags::magma_mode
#define DEFAULT_PRINT_MODE      matrix_flags::default_print_mode

//
// STORAGE MODE SETTINGS
//

#define REPRESENTATION		matrix_flags::representation
#define DENSE_REPRESENTATION	matrix_flags::dense_representation
#define SPARSE_REPRESENTATION   matrix_flags::sparse_representation
#define	MIXED_REPRESENTATION    matrix_flags::mixed_representation

#define ORIENTATION		matrix_flags::orientation
#define ROW_ORIENTED		matrix_flags::row_oriented
#define COLUMN_ORIENTED		matrix_flags::column_oriented

#define DEFAULT_STORAGE_MODE	matrix_flags::default_storage_mode

//
// STRUCTURE MODE SETTINGS
//

#define DIAG			matrix_flags::diag
#define UPPER_DIAG		matrix_flags::upper_diag
#define	LOWER_DIAG		matrix_flags::lower_diag
#define UPPER_TRIA		matrix_flags::upper_tria
#define LOWER_TRIA		matrix_flags::lower_tria;
#define COLUMNS_LININD		matrix_flags::columns_linind
#define ROWS_LININD		matrix_flags::rows_linind

#define DEFAULT_STRUCTURE_MODE	matrix_flags::default_structure_mode

//
// INFO MODE SETTINGS (Position der Diagonalen)
//

#define DIAG_UP			matrix_flags::diag_up
#define DIAG_RIGHT		matrix_flags::diag_right
#define DIAG_LD_TO_RU		matrix_flags::diag_ld_to_ru

#define DEFAULT_INFO_MODE	matrix_flags::default_info_mode

//
// LATTICE MODE SETTINGS
//

#define DEFAULT_LATTICE_MODE    matrix_flags::default_lattice_mode

//
// VECTOR SETTINGS
//

#define EXPAND                  vector_flags(vector_flags::expand)
#define FIXED                   vector_flags(vector_flags::fixed)
#define SORT_VECTOR_DOWN        vector_flags::sort_vector_down
#define SORT_VECTOR_UP          vector_flags::sort_vector_up



#endif	// LIDIA_LIDIA_DEFINES_H_GUARD_
