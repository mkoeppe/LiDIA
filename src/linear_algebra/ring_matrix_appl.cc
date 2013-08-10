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


#include	"LiDIA/ring_matrix.h"
#include	<fstream>

#ifndef MATRIX
# if defined SPARSE
#  define MATRIX sparse_ring_matrix
# elif defined DENSE
#  define MATRIX dense_ring_matrix
# else
#  define MATRIX ring_matrix
# endif
#endif

#ifndef TYPE
# define TYPE bigint
# define INCLUDE "LiDIA/bigint.h"
#endif

#define IN_NAME argv[1]
#define OUT_NAME argv[2]


#include	INCLUDE



#ifdef LIDIA_NAMESPACE
using namespace LiDIA;
#endif



int main_LiDIA (int argc, char* argv[])
{
	if (argc != 3) {
		std::cerr << "usage: " << argv[0] << " input-file output-file" << std::endl;
		return 4;
	}

	std::ifstream in(IN_NAME);
	std::ofstream dz(OUT_NAME);

	register int i, j;

	// pointer
	TYPE *pointer_1 = NULL;
	TYPE *pointer_2 = NULL;
	TYPE *pointer_3 = new TYPE[4];
	TYPE *pointer_4 = new TYPE[6];
	TYPE *pointer_5 = NULL;
	TYPE *pointer_6 = NULL;

	// vectors
	math_vector< TYPE > v1, v2;
	math_vector< TYPE > v3(4, 4);
	math_vector< TYPE > v4(6, 6);

	// scalar
	TYPE s1, s2, s3;
	in >> s2 >> s3;
	lidia_size_t c, r;

	// array
	TYPE **wertearray = new TYPE *[2];
	for (i = 0; i < 2; i++) {
		wertearray[i] = new TYPE[3];
		for (j = 0; j < 3; j++) {
			in >> wertearray[i][j];
			//std::cout << wertearray[i][j] << std::flush;
		}
	}

	std::cout << "**********************************************************" << std::endl;
	std::cout << "***             Test for class ring_matrix             ***" << std::endl;
	std::cout << "***                     Version 3.0                    ***" << std::endl;
	std::cout << "**********************************************************" << std::endl;
	std::cout.flush();

	std::cout << std::endl << "testing constructors" << std::flush;
	dz << "testing constructor" << std::endl << std::endl;

	MATRIX< TYPE > A;
	A.set_zero_element(0);
	MATRIX< TYPE > B(2, 3);
	B.set_zero_element(0);
	MATRIX< TYPE > C(2, 3, const_cast<const TYPE **>(wertearray));
	C.set_zero_element(0);
	MATRIX< TYPE > D(C);
	D.set_zero_element(0);
	MATRIX< TYPE > E = C;
	E.set_zero_element(0);
	MATRIX< TYPE > F = C;
	F.set_zero_element(0);
	MATRIX< TYPE > G(v3);
	G.set_zero_element(0);

	dz << A << B << C << D << E << F << G << std::flush;
	std::cout << ".............................completed" << std::endl;

	A.read_from_stream(in);
	B.read_from_stream(in);
	C.read_from_stream(in);

	std::cout << "testing set_print_mode" << std::flush;
	dz << "testing set_print_mode" << std::endl;
	dz << " set A BEAUTY_MODE" << std::endl;
	A.set_print_mode(BEAUTY_MODE);
	dz << " set B LIDIA_MODE" << std::endl;
	B.set_print_mode(LIDIA_MODE);
	dz << " set C GP_MODE" << std::endl;
	C.set_print_mode(GP_MODE);
	dz << " set D MAPLE_MODE" << std::endl;
	D.set_print_mode(MAPLE_MODE);
	dz << " set E MATHEMATICA_MODE" << std::endl;
	E.set_print_mode(MATHEMATICA_MODE);
	dz << " set F KASH_MODE" << std::endl;
	F.set_print_mode(KASH_MODE);
	dz << " set G LATEX_MODE" << std::endl;
	G.set_print_mode(LATEX_MODE);
	std::cout << "...........................completed" << std::endl;

	std::cout << "testing get_print_mode" << std::flush;
	dz << "testing get_print_mode" << std::endl;
	dz << A.get_print_mode() << std::endl;
	dz << B.get_print_mode() << std::endl;
	dz << C.get_print_mode() << std::endl;
	dz << D.get_print_mode() << std::endl;
	dz << E.get_print_mode() << std::endl;
	dz << F.get_print_mode() << std::endl;
	dz << G.get_print_mode() << std::endl;
	std::cout << "...........................completed" << std::endl;

	dz << D << std::endl;
	std::cout << "testing write_to_lidia" << std::flush;
	dz << "testing write_to_lidia" << std::endl;
	std::ofstream out0("LIDIA_I");
	D.write_to_stream(out0);
	out0.close();
	std::cout << "...........................completed" << std::endl;

	std::cout << "testing write_to_stream" << std::flush;
	dz << "testing write_to_stream" << std::endl;
	std::ofstream out1("LIDIA_II");
	D.write_to_stream(out1);
	out1.close();
	std::cout << "..........................completed" << std::endl;

	std::cout << "testing write_to_gp" << std::flush;
	dz << "testing write_to_gp" << std::endl;
	std::ofstream out2("GP");
	D.write_to_gp(out2);
	out2.close();
	std::cout << "..............................completed" << std::endl;

	std::cout << "testing write_to_maple" << std::flush;
	dz << "testing write_to_maple" << std::endl;
	std::ofstream out3("MAPLE");
	D.write_to_maple(out3);
	out3.close();
	std::cout << "...........................completed" << std::endl;

	std::cout << "testing write_to_mathematica" << std::flush;
	dz << "testing write_to_mathematica" << std::endl;
	std::ofstream out4("MATHEMATICA");
	D.write_to_mathematica(out4);
	out4.close();
	std::cout << ".....................completed" << std::endl;

	std::cout << "testing write_to_kash" << std::flush;
	dz << "testing write_to_kash" << std::endl;
	std::ofstream out5("KASH");
	D.write_to_kash(out5);
	out5.close();
	std::cout << "............................completed" << std::endl;

	A.set_print_mode(BEAUTY_MODE);
	B.set_print_mode(BEAUTY_MODE);
	C.set_print_mode(BEAUTY_MODE);
	D.set_print_mode(BEAUTY_MODE);
	E.set_print_mode(BEAUTY_MODE);
	F.set_print_mode(BEAUTY_MODE);
	G.set_print_mode(BEAUTY_MODE);

	std::cout << "testing read_from_lidia" << std::flush;
	dz << "testing read_from_lidia" << std::endl;
	std::ifstream in0("LIDIA_I");
	D.read_from_lidia(in0);
	dz << D << std::flush;
	in0.close();
	std::cout << "..........................completed" << std::endl;

	std::cout << "testing read_from_stream" << std::flush;
	dz << "testing read_from_stream" << std::endl;
	std::ifstream in1("LIDIA_II");
	D.read_from_stream(in1);
	dz << D << std::flush;
	in1.close();
	std::cout << ".........................completed" << std::endl;

	std::cout << "testing read_from_gp" << std::flush;
	dz << "testing read_from_gp" << std::endl;
	std::ifstream in2("GP");
	D.read_from_gp(in2);
	dz << D << std::flush;
	in2.close();
	std::cout << ".............................completed" << std::endl;

	std::cout << "testing read_from_maple" << std::flush;
	dz << "testing read_from_maple" << std::endl;
	std::ifstream in3("MAPLE");
	D.read_from_maple(in3);
	dz << D << std::flush;
	in3.close();
	std::cout << "..........................completed" << std::endl;

	std::cout << "testing read_from_mathematica" << std::flush;
	dz << "testing read_from_mathematica" << std::endl;
	std::ifstream in4("MATHEMATICA");
	D.read_from_mathematica(in4);
	dz << D << std::flush;
	in4.close();
	std::cout << "....................completed" << std::endl;

	std::cout << "testing read_from_kash" << std::flush;
	dz << "testing read_from_kash" << std::endl;
	std::ifstream in5("KASH");
	D.read_from_kash(in5);
	dz << D << std::flush;
	in5.close();
	std::cout << "...........................completed" << std::endl;

	dz << A << std::endl;

	std::cout << "testing member" << std::flush;
	dz << "testing member" << std::endl;
	s1 = A.member(2, 1);
	dz << s1 << std::endl;
	s1 = A(2, 1);
	dz << s1 << std::endl;
	std::cout << "...................................completed" << std::endl;

	std::cout << "testing column" << std::flush;
	dz << "testing column" << std::endl;
	v1 = A(0);
	dz << v1 << std::endl;
	pointer_1 = A.get_column(1);
	for (i = 0; i < 4; i++)
		dz << pointer_1[i] << " " << std::flush;
	dz << std::endl;
	A.get_column(pointer_1, 2);
	for (i = 0; i < 4; i++)
		dz << pointer_1[i] << " " << std::flush;
	dz << std::endl;
	A.get_column_vector(v1, 2);
	dz << v1 << std::endl;
	v1 = A.get_column_vector(1);
	dz << v1 << std::endl;
	std::cout << "...................................completed" << std::endl;

	std::cout << "testing row" << std::flush;
	dz << "testing row" << std::endl;
	v2 = A[0];
	dz << v2 << std::endl;
	dz << std::endl;
	pointer_2 = A.get_row(1);
	for (i = 0; i < 6; i++)
		dz << pointer_2[i] << " " << std::flush;
	dz << std::endl;
	A.get_row(pointer_2, 2);
	for (i = 0; i < 6; i++)
		dz << pointer_2[i] << " " << std::flush;
	dz << std::endl;
	A.get_row_vector(v2, 2);
	dz << v2 << std::endl;
	v2 = A.get_row_vector(1);
	dz << v2 << std::endl;
	std::cout << "......................................completed" << std::endl;

	std::cout << "testing sto" << std::flush;
	E = A;
	dz << "testing sto" << std::endl;
	E.sto(2, 0, s2);
	dz << E << std::flush;
	std::cout << "......................................completed" << std::endl;

	std::cout << "testing sto_column" << std::flush;
	dz << "testing sto_column" << std::endl;
	E = A;
	E.sto_column(pointer_1, 3, 0);
	dz << E << std::flush;
	E.sto_column(pointer_1, 2, 1, 1);
	dz << E << std::flush;
	E.sto_column_vector(v1, 3, 2);
	dz << E << std::flush;
	E.sto_column_vector(v1, 2, 3, 1);
	dz << E << std::flush;
	std::cout << "...............................completed" << std::endl;

	std::cout << "testing sto_row" << std::flush;
	dz << "testing sto_row" << std::endl;
	E = A;
	E.sto_row(pointer_2, 4, 0);
	dz << E << std::flush;
	E.sto_row(pointer_2, 3, 1, 1);
	dz << E << std::flush;
	E.sto_row_vector(v2, 4, 2);
	dz << E << std::flush;
	E.sto_row_vector(v2, 3, 3, 1);
	dz << E << std::flush;
	std::cout << "..................................completed" << std::endl;

	std::cout << "testing get_data" << std::flush;
	dz << "testing get_data" << std::endl;
	r = E.get_no_of_rows();
	c = E.get_no_of_columns();
	TYPE **value = E.get_data();
	for (i = 0; i < r; i++)
		for (j = 0; j < c; j++)
			dz << value[i][j] << " ";
	dz << std::endl;
	std::cout << ".................................completed" << std::endl;

	std::cout << "testing set_data" << std::flush;
	dz << "testing set_data" << std::endl;
	E = A;
	E.set_data(const_cast<const TYPE **>(value), r, c);
	dz << E << std::endl;
	for (i = 0; i < r; i++)
		delete[] value[i];
	delete[] value;
	std::cout << ".................................completed" << std::endl;

	std::cout << "testing swap" << std::flush;
	dz << "testing swap" << std::endl;
	dz << E << A << std::flush;
	swap(E, A);
	dz << E << A << std::flush;
	swap(E, A);
	std::cout << ".....................................completed" << std::endl;

	std::cout << "testing swap_columns" << std::flush;
	dz << "testing swap_columns" << std::endl;
	E = A;
	E.swap_columns(0, 5);
	dz << E << std::flush;
	std::cout << ".............................completed" << std::endl;

	std::cout << "testing swap_rows" << std::flush;
	dz << "testing swap_rows" << std::endl;
	E = A;
	E.swap_rows(3, 0);
	dz << E << std::flush;
	std::cout << "................................completed" << std::endl;

	std::cout << "testing split_t" << std::flush;
	dz << "testing split_t" << std::endl;
	MATRIX< TYPE > Part1(2, 4), Part2(3, 2), Part3(2, 4), Part4(1, 2);
	A.split_t(Part1, Part2, Part3, Part4);
	dz << Part1 << Part2 << Part3 << Part4 << std::flush;
	std::cout << "..................................completed" << std::endl;

	std::cout << "testing split_h" << std::flush;
	dz << "testing split_h" << std::endl;
	MATRIX< TYPE > Part5(4, 4), Part6(4, 2);
	MATRIX< TYPE > PartA(4, 5), PartB(4, 5), PartC(4, 5), PartD(4, 5);
	A.split_h(Part5, Part6);
	dz << Part5 << Part6 << std::flush;

	A.split_h(pointer_1, PartA);
	dz << PartA;
	for (i = 0; i < 4; i++)
		dz << pointer_1[i] << " ";
	dz << std::endl;
	A.split_h(v1, PartB);
	dz << PartB << v1 << std::endl;

	A.split_h(PartC, pointer_3);
	dz << PartC;
	for (i = 0; i < 4; i++)
		dz << pointer_3[i] << " ";
	dz << std::endl;
	A.split_h(PartD, v3);
	dz << PartD << v3 << std::endl;
	std::cout << "..................................completed" << std::endl;

	std::cout << "testing split_v" << std::flush;
	dz << "testing split_v" << std::endl;
	MATRIX< TYPE > Part7(1, 6), Part8(3, 6);
	MATRIX< TYPE > PartE(3, 6), PartF(3, 6), PartG(3, 6), PartH(3, 6);
	A.split_v(Part7, Part8);
	dz << Part7 << Part8 << std::flush;

	A.split_v(pointer_2, PartE);
	dz << PartE;
	for (i = 0; i < 6; i++)
		dz << pointer_2[i] << " ";
	dz << std::endl;
	A.split_v(v2, PartF);
	dz << PartF << v2 << std::endl;

	A.split_v(PartG, pointer_4);
	dz << PartG;
	for (i = 0; i < 6; i++)
		dz << pointer_4[i] << " ";
	dz << std::endl;
	A.split_v(PartH, v4);
	dz << PartH << v4 << std::endl;
	std::cout << "..................................completed" << std::endl;

	std::cout << "testing compose_t" << std::flush;
	dz << "testing compose_t" << std::endl;
	{
		MATRIX< TYPE > TMP(4, 6);
		TMP.compose_t(Part1, Part2, Part3, Part4);
		dz << TMP << std::flush;
	}
	std::cout << "................................completed" << std::endl;

	std::cout << "testing compose_h" << std::flush;
	dz << "testing compose_h" << std::endl;
	{
		MATRIX< TYPE > TMP(4, 6);
		TMP.compose_h(Part5, Part6);
		dz << TMP << std::flush;

		TMP.compose_h(pointer_1, PartA);
		dz << TMP << std::flush;

		TMP.compose_h(v1, PartB);
		dz << TMP << std::flush;

		TMP.compose_h(PartC, pointer_3);
		dz << TMP << std::flush;

		TMP.compose_h(PartD, v3);
		dz << TMP << std::flush;
	}
	std::cout << "................................completed" << std::endl;

	std::cout << "testing compose_v" << std::flush;
	dz << "testing compose_v" << std::endl;
	{
		MATRIX< TYPE > TMP(4, 6);
		TMP.compose_v(Part7, Part8);
		dz << TMP << std::flush;

		TMP.compose_v(pointer_2, PartE);
		dz << TMP << std::flush;

		TMP.compose_v(v2, PartF);
		dz << TMP << std::flush;

		TMP.compose_v(PartG, pointer_4);
		dz << TMP << std::flush;

		TMP.compose_v(PartH, v4);
		dz << TMP << std::flush;
	}
	std::cout << "................................completed" << std::endl;

	std::cout << "testing get_no_of_columns" << std::flush;
	dz << "testing get_no_of_columns" << std::endl;
	s1 = A.get_no_of_columns();
	dz << s1 << std::endl;
	std::cout << "........................completed" << std::endl;

	std::cout << "testing get_no_of_rows" << std::flush;
	dz << "testing get_no_of_rows" << std::endl;
	s1 = A.get_no_of_rows();
	dz << s1 << std::endl;
	std::cout << "...........................completed" << std::endl;

	std::cout << "testing set_no_of_columns" << std::flush;
	dz << "set_no_of_columns" << std::endl << std::flush;
	E = A;
	E.set_no_of_columns(3);
	dz << E << std::flush;
	std::cout << "........................completed" << std::endl;

	std::cout << "testing set_no_of_rows" << std::flush;
	dz << "set_no_of_rows" << std::endl;
	E.set_no_of_rows(3);
	dz << E << std::flush;
	std::cout << "...........................completed" << std::endl;

	std::cout << "testing resize" << std::flush;
	dz << "resize" << std::endl;
	E = A;
	E.resize(3, 3);
	dz << E << std::flush;
	std::cout << "...................................completed" << std::endl;

	std::cout << "testing kill / reset" << std::flush;
	dz << "kill / reset" << std::endl;
	E = A;
	E.kill();
	dz << E << std::flush;
	std::cout << ".............................completed" << std::endl;

	std::cout << "testing assign" << std::flush;
	dz << "testing assign" << std::endl;
	E = A;
	dz << E << std::flush;
	E.kill();
	assign(E, A);
	dz << E << std::flush;
	std::cout << "...................................completed" << std::endl;

	std::cout << "testing trans" << std::flush;
	dz << "testing trans" << std::endl;
	E = trans(A);
	dz << E << std::flush;
	E = E.trans();
	dz << E << std::flush;
	E.trans(E);
	dz << E << std::flush;
	trans(E, E);
	dz << E << std::flush;
	std::cout << "....................................completed" << std::endl;

	std::cout << "testing diag" << std::flush;
	dz << "diag" << std::endl;
	E.diag(s1, s2);
	dz << E << std::flush;
	diag(E, s2, s3);
	dz << E << std::flush;
	std::cout << ".....................................completed" << std::endl;
	G.diag(0, 0);
	dz << E << G << std::endl;

	std::cout << "testing is_column_zero" << std::flush;
	dz << "is_column_zero" << std::endl;
	dz << E.is_column_zero(1) << std::endl;
	dz << G.is_column_zero(0) << std::endl;
	std::cout << "...........................completed" << std::endl;

	std::cout << "testing is_row_zero" << std::flush;
	dz << "is_row_zero" << std::endl;
	dz << E.is_row_zero(1) << std::endl;
	dz << G.is_row_zero(0) << std::endl;
	std::cout << "..............................completed" << std::endl;

	std::cout << "testing is_matrix_zero" << std::flush;
	dz << "is_matrix_zero" << std::endl;
	dz << E.is_matrix_zero() << std::endl;
	dz << G.is_matrix_zero() << std::endl;
	std::cout << "...........................completed" << std::endl;

	//***************************** END   OF PART: BASE_MATIRX ****************
	//***************************** BEGIN Of PART: MATH_MATRIX ****************

	dz << A << B << C << D << E << std::flush;
	dz << v1 << v2 << v3 << v4 << std::flush;

	std::cout << "testing addition" << std::flush;
	dz << "testing addition" << std::endl << std::flush;
	A += B;
	dz << A << std::flush;
	add(A, A, B);
	dz << A << std::flush;
	F = A + B;
	dz << F << std::flush;
	add(F, A, B);
	dz << F << std::flush;
	std::cout << ".................................completed" << std::endl << std::flush;

	std::cout << "testing subtraction" << std::flush;
	dz << "testing subtraction" << std::endl << std::flush;
	A -= B;
	dz << A << std::flush;
	subtract(A, A, B);
	dz << A << std::flush;
	F = A - B;
	dz << F << std::flush;
	subtract(F, A, B);
	dz << F << std::flush;
	std::cout << "..............................completed" << std::endl << std::flush;

	std::cout << "testing negate" << std::flush;
	dz << "testing negate" << std::endl << std::flush;
	F = -A;
	dz << F << std::flush;
	negate(F, A);
	dz << F << std::flush;
	std::cout << "...................................completed" << std::endl << std::flush;

	std::cout << "testing multiplication" << std::flush;
	dz << "testing multiplication" << std::endl << std::flush;
	E = A, F = A;
	E *= C;
	dz << E << std::flush;
	E = A * C;
	dz << E << std::flush;
	multiply(F, A, C);
	dz << F << std::flush;
	std::cout << "...........................completed" << std::endl << std::flush;

	std::cout << "testing addition with scalar" << std::flush;
	dz << "testing addition with scalar" << std::endl << std::flush;
	A += s2;
	dz << A << std::flush;
	F = A + s2;
	dz << F << std::flush;
	add(F, A, s2);
	dz << F << std::flush;
	std::cout << ".....................completed" << std::endl << std::flush;

	std::cout << "testing subtraction with scalar" << std::flush;
	dz << "testing subtraction with scalar" << std::endl << std::flush;
	A -= s2;
	dz << A << std::flush;
	F = A - s2;
	dz << F << std::flush;
	subtract(F, A, s2);
	dz << F << std::flush;
	std::cout << "..................completed" << std::endl << std::flush;

	std::cout << "testing multiplication with scalar" << std::flush;
	dz << "testing multiplication with scalar" << std::endl << std::flush;
	F = A;
	F *= s3;
	dz << F << std::flush;
	F = A * s3;
	dz << F << std::flush;
	multiply(F, A, s3);
	dz << F << std::flush;
	std::cout << "...............completed" << std::endl << std::flush;

	std::cout << "testing multiplication with array" << std::flush;
	dz << "testing multiplication with array" << std::endl << std::flush;
	pointer_3 = A * pointer_2;
	for (i = 0; i < 4; i++)
		dz << pointer_3[i] << " ";
	dz << std::endl << std::flush;

	multiply(pointer_5, A, pointer_2);
	for (i = 0; i < 4; i++)
		dz << pointer_5[i] << " ";
	dz << std::endl << std::flush;

	pointer_4 = pointer_1 * A;
	for (i = 0; i < 6; i++)
		dz << pointer_4[i] << " ";
	dz << std::endl << std::flush;
	multiply(pointer_6, pointer_1, A);
	for (i = 0; i < 6; i++)
		dz << pointer_6[i] << " ";
	dz << std::endl << std::flush;
	std::cout << "................completed" << std::endl << std::flush;

	std::cout << "testing multiplication with vector" << std::flush;
	dz << "testing multiplication with vector" << std::endl << std::flush;
	v3 = A * v2;
	dz << v3 << std::endl << std::flush;
	multiply(v3, A, v2);
	dz << v3 << std::endl << std::flush;

	v4 = v1 * A;
	dz << v4 << std::endl << std::flush;
	multiply(v4, v1, A);
	dz << v4 << std::endl << std::flush;
	std::cout << "...............completed" << std::endl << std::flush;

	std::cout << "testing trace" << std::flush;
	dz << "testing trace" << std::endl << std::flush;
	F = A * (MATRIX< TYPE > )trans(A);
	s1 = F.trace();
	dz << s1 << std::endl << std::flush;
	F.trace(s1);
	dz << s1 << std::endl << std::flush;
	s1 = trace(F);
	dz << s1 << std::endl << std::flush;
	std::cout << "....................................completed" << std::endl << std::flush;

	dz.close();
	std::cout << "\nPlease use now the 'diff' or 'cmp' command to verify" << std::endl;
	std::cout << "the equality of {dense_, sparse_}ring_matrix_appl_< type > .dat";
	std::cout << " and {dense_, sparse_}ring_matrix_appl_< type > .out." << std::endl;

	return 0;
}


int main(int argc, char** argv) {

#if defined(LIDIA_EXCEPTIONS)
    try {
#endif

	main_LiDIA(argc, argv);
	
#if defined(LIDIA_EXCEPTIONS)
    }
    catch(basic_error const& ex) {
	ex.traditional_error_handler();
	return 1;
    }
    catch(std::exception const& ex) {
	std::cerr << "unexpected exception: " << ex.what() << "\n";
	return 1;
    }
#endif
    
}
