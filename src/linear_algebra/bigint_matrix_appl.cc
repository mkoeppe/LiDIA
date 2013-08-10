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

// Application Test Routine     


#include	"LiDIA/bigint.h"
#include	"LiDIA/bigint_matrix.h"
#include	<fstream>

#define IN_NAME argv[1]
#define OUT_NAME argv[2]



#ifdef LIDIA_NAMESPACE
using namespace LiDIA;
#endif



int main_LiDIA(int argc, char** argv)
{
	if (argc != 3) {
		std::cerr << "usage: " << argv[0] << " input-file output-file" << std::endl;
		return 4;
	}

	std::ifstream in(IN_NAME);
	std::ofstream dz(OUT_NAME);

	int i, j;

	// pointer
	bigint *pointer_1 = NULL;
	bigint *pointer_2 = NULL;
	bigint *pointer_3 = new bigint[4];
	bigint *pointer_4 = new bigint[6];
	bigint *pointer_5 = new bigint[4];
	bigint *pointer_6 = new bigint[6];

	// vectors
	math_vector< bigint > v1, v2;
	math_vector< bigint > v3(4, 4);
	math_vector< bigint > v4(6, 6);

	// scalar
	bigint scalar_1, scalar_2, scalar_3;
	in >> scalar_2 >> scalar_3;
	//std::cout << scalar_2 << scalar_3 << std::flush;
	lidia_size_t c, r;

	// array
	bigint **wertearray = new bigint *[2];
	for (i = 0; i < 2; i++) {
		wertearray[i] = new bigint[3];
		for (j = 0; j < 3; j++) {
			in >> wertearray[i][j];
			//std::cout << wertearray[i][j] << std::flush;
		}
	}

	std::cout << "**********************************************************" << std::endl;
	std::cout << "***             Test for class bigint_matrix           ***" << std::endl;
	std::cout << "***                     Version 2.0                    ***" << std::endl;
	std::cout << "**********************************************************" << std::endl;
	std::cout.flush();

	std::cout << "\ntesting constructors" << std::flush;
	dz << "testing constructor" << std::endl << std::flush;
	bigint_matrix  A(3, 1);
	bigint_matrix  B(2, 3);
	bigint_matrix  C(2, 3, const_cast<const bigint **>(wertearray));
	bigint_matrix  D(C);
	bigint_matrix  E = C;
	bigint_matrix  F = C;
	dz << A << B << C << D << E << F << std::flush;
	std::cout << ".............................completed\n" << std::flush;

	A.read_from_stream(in);
	B.read_from_stream(in);
	C.read_from_stream(in);

	std::cout << "testing set_print_mode" << std::flush;
	dz << "testing set_print_mode\n" << std::flush;
	dz << " set A BEAUTY_MODE\n" << std::flush;
	A.set_print_mode(BEAUTY_MODE);
	dz << " set B LIDIA_MODE\n" << std::flush;
	B.set_print_mode(LIDIA_MODE);
	dz << " set C GP_MODE\n" << std::flush;
	C.set_print_mode(GP_MODE);
	dz << " set D MAPLE_MODE\n" << std::flush;
	D.set_print_mode(MAPLE_MODE);
	dz << " set E MATHEMATICA_MODE\n" << std::flush;
	E.set_print_mode(MATHEMATICA_MODE);
	dz << " set E KASH_MODE\n" << std::flush;
	E.set_print_mode(KASH_MODE);
	std::cout << "...........................completed\n" << std::flush;

	std::cout << "testing get_print_mode" << std::flush;
	dz << "testing get_print_mode\n" << std::flush;
	dz << A.get_print_mode() << std::endl << std::flush;
	dz << B.get_print_mode() << std::endl << std::flush;
	dz << C.get_print_mode() << std::endl << std::flush;
	dz << D.get_print_mode() << std::endl << std::flush;
	dz << E.get_print_mode() << std::endl << std::flush;
	dz << F.get_print_mode() << std::endl << std::flush;
	std::cout << "...........................completed\n" << std::flush;

	std::cout << "testing write_to_stream" << std::flush;
	dz << "testing write_to_stream\n" << std::flush;
	std::ofstream out1("LIDIA");
	D.write_to_stream(out1);
	out1.close();
	std::cout << "..........................completed\n" << std::flush;

	std::cout << "testing write_to_gp" << std::flush;
	dz << "testing write_to_gp\n" << std::flush;
	std::ofstream out2("GP");
	D.write_to_gp(out2);
	out2.close();
	std::cout << "..............................completed\n" << std::flush;

	std::cout << "testing write_to_maple" << std::flush;
	dz << "testing write_to_maple\n" << std::flush;
	std::ofstream out3("MAPLE");
	D.write_to_maple(out3);
	out3.close();
	std::cout << "...........................completed\n" << std::flush;

	std::cout << "testing write_to_mathematica" << std::flush;
	dz << "testing write_to_mathematica\n" << std::flush;
	std::ofstream out4("MATHEMATICA");
	D.write_to_mathematica(out4);
	out4.close();
	std::cout << ".....................completed\n" << std::flush;

	std::cout << "testing write_to_kash" << std::flush;
	dz << "testing write_to_kash\n" << std::flush;
	std::ofstream out5("KASH");
	D.write_to_kash(out5);
	out5.close();
	std::cout << "............................completed\n" << std::flush;

	A.set_print_mode(BEAUTY_MODE);
	B.set_print_mode(BEAUTY_MODE);
	C.set_print_mode(BEAUTY_MODE);
	D.set_print_mode(BEAUTY_MODE);
	E.set_print_mode(BEAUTY_MODE);
	F.set_print_mode(BEAUTY_MODE);

	std::cout << "testing read_from_stream" << std::flush;
	dz << "testing read_from_stream\n" << std::flush;
	std::ifstream in1("LIDIA");
	D.read_from_stream(in1);
	dz << D << std::flush;
	in1.close();
	std::cout << ".........................completed\n" << std::flush;

	std::cout << "testing read_from_gp" << std::flush;
	dz << "testing read_from_gp\n" << std::flush;
	std::ifstream in2("GP");
	D.read_from_gp(in2);
	dz << D << std::flush;
	in2.close();
	std::cout << ".............................completed\n" << std::flush;

	std::cout << "testing read_from_maple" << std::flush;
	dz << "testing read_from_maple\n" << std::flush;
	std::ifstream in3("MAPLE");
	D.read_from_maple(in3);
	dz << D << std::flush;
	in3.close();
	std::cout << "..........................completed\n" << std::flush;

	std::cout << "testing read_from_mathematica" << std::flush;
	dz << "testing read_from_mathematica\n" << std::flush;
	std::ifstream in4("MATHEMATICA");
	D.read_from_mathematica(in4);
	dz << D << std::flush;
	in4.close();
	std::cout << "....................completed\n" << std::flush;

	std::cout << "testing read_from_kash" << std::flush;
	dz << "testing read_from_kash\n" << std::flush;
	std::ifstream in5("KASH");
	D.read_from_kash(in5);
	in5.close();
	std::cout << "...........................completed\n" << std::flush;

	std::cout << "testing member" << std::flush;
	dz << "testing member" << std::endl << std::flush;
	scalar_1 = A.member(2, 1);
	dz << scalar_1 << std::endl;
	scalar_1 = A(2, 1);
	dz << scalar_1 << std::endl << std::flush;
	std::cout << "...................................completed\n" << std::flush;

	std::cout << "testing column" << std::flush;
	dz << "testing column" << std::endl << std::flush;
	v1 = A(0);
	dz << v1 << std::endl << std::flush;
	dz << "\n" << std::flush;
	pointer_1 = A.get_column(1);
	for (i = 0; i < 4; i++)
		dz << pointer_1[i] << " " << std::flush;
	dz << "\n" << std::flush;
	A.get_column(pointer_1, 2);
	for (i = 0; i < 4; i++)
		dz << pointer_1[i] << " " << std::flush;
	dz << "\n" << std::flush;
	A.get_column_vector(v1, 2);
	dz << v1 << std::endl << std::flush;
	std::cout << "...................................completed\n" << std::flush;

	std::cout << "testing row" << std::flush;
	dz << "testing row" << std::endl << std::flush;
	v2 = A[0];
	dz << v2 << std::endl << std::flush;
	dz << "\n" << std::flush;
	pointer_2 = A.get_row(1);
	for (i = 0; i < 6; i++)
		dz << pointer_2[i] << " " << std::flush;
	dz << "\n" << std::flush;
	A.get_row(pointer_2, 2);
	for (i = 0; i < 6; i++)
		dz << pointer_2[i] << " " << std::flush;
	dz << "\n" << std::flush;
	A.get_row_vector(v2, 2);
	dz << v2 << std::endl << std::flush;
	std::cout << "......................................completed\n" << std::flush;

	std::cout << "testing sto" << std::flush;
	E = A;
	dz << "testing sto" << std::endl << std::flush;
	E.sto(2, 0, 23);
	dz << E << std::flush;
	std::cout << "......................................completed\n" << std::flush;

	std::cout << "testing sto_column" << std::flush;
	dz << "testing sto_column" << std::endl << std::flush;
	E = A;
	E.sto_column(pointer_1, 3, 0);
	dz << E << std::flush;
	E.sto_column_vector(v1, 3, 1);
	dz << E << std::flush;
	std::cout << "...............................completed\n" << std::flush;

	std::cout << "testing sto_row" << std::flush;
	dz << "testing sto_row" << std::endl << std::flush;
	E = A;
	E.sto_row(pointer_2, 4, 0);
	dz << E << std::flush;
	E.sto_row_vector(v2, 4, 1);
	dz << E << std::flush;
	std::cout << "..................................completed\n" << std::flush;

	std::cout << "testing get_data" << std::flush;
	dz << "testing get_data" << std::endl << std::flush;
	r = E.get_no_of_rows();
	c = E.get_no_of_columns();
	bigint **value = E.get_data();
	for (i = 0; i < r; i++)
		for (j = 0; j < c; j++)
			dz << value[i][j] << " ";
	dz << std::endl << std::flush;
	for (i = 0; i < r; i++)
		delete[] value[i];
	delete[] value;
	std::cout << ".................................completed\n" << std::flush;

	std::cout << "testing swap" << std::flush;
	dz << "testing swap" << std::endl << std::flush;
	dz << E << A << std::flush;
	swap(E, A);
	dz << E << A << std::flush;
	swap(E, A);
	std::cout << ".....................................completed\n" << std::flush;

	std::cout << "testing swap_columns" << std::flush;
	dz << "testing swap_columns" << std::endl << std::flush;
	E = A;
	E.swap_columns(0, 5);
	dz << E << std::flush;
	std::cout << ".............................completed\n" << std::flush;

	std::cout << "testing swap_rows" << std::flush;
	dz << "testing swap_rows" << std::endl << std::flush;
	E = A;
	E.swap_rows(3, 0);
	dz << E << std::flush;
	std::cout << "................................completed\n" << std::flush;

	std::cout << "testing split_t" << std::flush;
	dz << "testing split_t" << std::endl << std::flush;
	bigint_matrix  Part1(2, 4), Part2(3, 2), Part3(2, 4), Part4(1, 2);
	A.split_t(Part1, Part2, Part3, Part4);
	dz << Part1 << Part2 << Part3 << Part4 << std::flush;
	std::cout << "..................................completed\n" << std::flush;

	std::cout << "testing split_h" << std::flush;
	dz << "testing split_h" << std::endl << std::flush;
	bigint_matrix  Part5(4, 4), Part6(4, 2);
	bigint_matrix  PartA(4, 5), PartB(4, 5), PartC(4, 5), PartD(4, 5);
	A.split_h(Part5, Part6);
	dz << Part5 << Part6 << std::flush;

	A.split_h(pointer_1, PartA);
	dz << PartA;
	for (i = 0; i < 4; i++)
		dz << pointer_1[i] << " ";
	dz << std::endl << std::flush;
	A.split_h(v1, PartB);
	dz << PartB << v1 << std::endl << std::flush;

	A.split_h(PartC, pointer_3);
	dz << PartC;
	for (i = 0; i < 4; i++)
		dz << pointer_3[i] << " ";
	dz << std::endl << std::flush;
	A.split_h(PartD, v3);
	dz << PartD << v3 << std::endl << std::flush;
	std::cout << "..................................completed\n" << std::flush;

	std::cout << "testing split_v" << std::flush;
	dz << "testing split_v" << std::endl << std::flush;
	bigint_matrix  Part7(1, 6), Part8(3, 6);
	bigint_matrix  PartE(3, 6), PartF(3, 6), PartG(3, 6), PartH(3, 6);
	A.split_v(Part7, Part8);
	dz << Part7 << Part8 << std::flush;

	A.split_v(pointer_2, PartE);
	dz << PartE;
	for (i = 0; i < 6; i++)
		dz << pointer_2[i] << " ";
	dz << std::endl << std::flush;
	A.split_v(v2, PartF);
	dz << PartF << v2 << std::endl << std::flush;

	A.split_v(PartG, pointer_4);
	dz << PartG;
	for (i = 0; i < 6; i++)
		dz << pointer_4[i] << " ";
	dz << std::endl << std::flush;
	A.split_v(PartH, v4);
	dz << PartH << v4 << std::endl << std::flush;
	std::cout << "..................................completed\n" << std::flush;

	std::cout << "testing compose_t" << std::flush;
	dz << "testing compose_t\n" << std::flush;
	{
		bigint_matrix  TMP(4, 6);
		TMP.compose_t(Part1, Part2, Part3, Part4);
		dz << TMP << std::flush;
	}
	std::cout << "................................completed\n" << std::flush;

	std::cout << "testing compose_h" << std::flush;
	dz << "testing compose_h" << std::endl << std::flush;
	{
		bigint_matrix  TMP(4, 6);
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
	std::cout << "................................completed\n" << std::flush;

	std::cout << "testing compose_v" << std::flush;
	dz << "testing compose_v" << std::endl << std::flush;
	{
		bigint_matrix TMP(4, 6);
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
	std::cout << "................................completed\n" << std::flush;

	std::cout << "testing get_no_of_columns" << std::flush;
	dz << "testing get_no_of_columns" << std::endl << std::flush;
	scalar_1 = A.get_no_of_columns();
	dz << scalar_1 << std::endl << std::flush;
	std::cout << "........................completed\n" << std::flush;

	std::cout << "testing get_no_of_rows" << std::flush;
	dz << "testing get_no_of_rows" << std::endl << std::flush;
	scalar_1 = A.get_no_of_rows();
	dz << scalar_1 << std::endl << std::flush;
	std::cout << "...........................completed\n" << std::flush;

	std::cout << "testing set_no_of_columns" << std::flush;
	dz << "set_no_of_columns" << std::endl << std::flush;
	E = A;
	E.set_no_of_columns(3);
	dz << E << std::flush;
	std::cout << "........................completed\n" << std::flush;

	std::cout << "testing set_no_of_rows" << std::flush;
	dz << "set_no_of_rows" << std::endl << std::flush;
	E.set_no_of_rows(3);
	dz << E << std::flush;
	std::cout << "...........................completed\n" << std::flush;

	std::cout << "testing assign" << std::flush;
	dz << "testing assign" << std::endl << std::flush;
	E = A;
	dz << E << std::flush;
	E.assign(A);
	dz << E << std::flush;
	assign(E, A);
	dz << E << std::flush;
	std::cout << "...................................completed\n" << std::flush;

	std::cout << "testing trans" << std::flush;
	dz << "testing trans" << std::endl << std::flush;
	E = trans(A);
	dz << E << std::flush;
	E = E.trans();
	dz << E << std::flush;
	E.trans(E);
	dz << E << std::flush;
	trans(E, E);
	dz << E << std::flush;
	std::cout << "....................................completed\n" << std::flush;

	std::cout << "testing diag" << std::flush;
	dz << "diag" << std::endl << std::flush;
	E.diag(scalar_1, scalar_2);
	dz << E << std::flush;
	diag(E, scalar_2, scalar_3);
	dz << E << std::flush;
	std::cout << ".....................................completed\n" << std::flush;

	//***************************** END   OF PART: BASE_MATIRX ***************
	//**************************** BEGIN Of PART: MATH_MATRIX ****************

	dz << A << B << C << D << E << std::flush;
	dz << v1 << v2 << v3 << v4 << std::endl << std::flush;

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
	std::cout << ".................................completed\n" << std::flush;

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
	std::cout << "..............................completed\n" << std::flush;

	std::cout << "testing negate" << std::flush;
	dz << "testing negate" << std::endl << std::flush;
	F = -A;
	dz << F << std::flush;
	negate(F, A);
	dz << F << std::flush;
	std::cout << "...................................completed\n" << std::flush;

	std::cout << "testing multiplication" << std::flush;
	dz << "testing multiplication" << std::endl << std::flush;
	E = A, F = A;
	E *= C;
	dz << E << std::flush;
	E = A * C;
	dz << E << std::flush;
	multiply(F, A, C);
	dz << F << std::flush;
	std::cout << "...........................completed\n" << std::flush;

	std::cout << "testing addition with scalar" << std::flush;
	dz << "testing addition with scalar" << std::endl << std::flush;
	A += scalar_2;
	dz << A << std::flush;
	F = A + scalar_2;
	dz << F << std::flush;
	add(F, A, scalar_2);
	dz << F << std::flush;
	std::cout << ".....................completed\n" << std::flush;

	std::cout << "testing subtraction with scalar" << std::flush;
	dz << "testing subtraction with scalar" << std::endl << std::flush;
	A -= scalar_2;
	dz << A << std::flush;
	F = A - scalar_2;
	dz << F << std::flush;
	subtract(F, A, scalar_2);
	dz << F << std::flush;
	std::cout << "..................completed\n" << std::flush;

	std::cout << "testing multiplication with scalar" << std::flush;
	dz << "testing multiplication with scalar" << std::endl << std::flush;
	F = A;
	F *= scalar_3;
	dz << F << std::flush;
	F = A * scalar_3;
	dz << F << std::flush;
	multiply(F, A, scalar_3);
	dz << F << std::flush;
	std::cout << "...............completed\n" << std::flush;

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
	std::cout << "................completed\n" << std::flush;

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
	std::cout << "...............completed\n" << std::flush;

	std::cout << "testing trace" << std::flush;
	dz << "testing trace" << std::endl << std::flush;
	F = A * trans(A);
	scalar_1 = F.trace();
	dz << scalar_1 << std::endl << std::flush;
	F.trace(scalar_1);
	dz << scalar_1 << std::endl << std::flush;
	scalar_1 = trace(F);
	dz << scalar_1 << std::endl << std::flush;
	std::cout << "....................................completed\n" << std::flush;

	//***************************** END   OF PART: MATH_MATIRX ***************
	//**************************** BEGIN Of PART: BIGINT_MATRIX **************

	std::cout << "testing remainder" << std::flush;
	dz << "testing remainder" << std::endl << std::flush;
	F = A;
	F %= 13;
	dz << F << std::flush;
	F = A % 13;
	dz << F << std::flush;
	remainder(F, A, bigint(13));
	dz << F << std::flush;
	std::cout << "................................completed\n" << std::flush;

	std::cout << "testing equal / unequal" << std::flush;
	bigint_matrix G = F;
	dz << "testing equal/unequal" << std::endl << std::flush;
	dz << (F == G) << std::endl << std::flush;
	dz << F.equal(G) << std::endl << std::flush;
	dz << equal(F, G) << std::endl << std::flush;
	dz << (F != G) << std::endl << std::flush;
	dz << F.unequal(G) << std::endl << std::flush;
	dz << unequal(F, G) << std::endl << std::flush;
	std::cout << "..........................completed\n" << std::flush;

	std::cout << "testing is_column_zero" << std::flush;
	D.diag(bigint(0), bigint(0));
	dz << "testing is_column_zero" << std::endl << std::flush;

	//dz << D.is_column_zero(0) << std::endl << std::flush;
	if (D.is_column_zero(0) == true)
		dz << 1 << std::endl;
	else
		dz << 0 << std::endl;

	if (F.is_column_zero(0) == true)
		dz << 1 << std::endl;
	else
		dz << 0 << std::endl;
	std::cout << "...........................completed\n" << std::flush;

	std::cout << "testing is_row_zero" << std::flush;
	dz << "testing is_row_zero" << std::endl << std::flush;
	if (D.is_row_zero(0) == true)
		dz << 1 << std::endl;
	else
		dz << 0 << std::endl;

	if (F.is_row_zero(0) == true)
		dz << 1 << std::endl;
	else
		dz << 0 << std::endl;
	std::cout << "..............................completed\n" << std::flush;

	std::cout << "testing is_matrix_zero" << std::flush;
	dz << "testing is_matrix_zero" << std::endl << std::flush;
	if (D.is_matrix_zero() == true)
		dz << 1 << std::endl;
	else
		dz << 0 << std::endl;

	if (F.is_matrix_zero() == true)
		dz << 1 << std::endl;
	else
		dz << 0 << std::endl;
	std::cout << "...........................completed\n" << std::flush;

	std::cout << "testing max / min " << std::flush;
	dz << "testing max / min" << std::endl << std::flush;
	A.max(scalar_1);
	dz << scalar_1 << std::endl << std::flush;
	dz << A.max() << std::endl << std::flush;
	dz << max(A) << std::endl << std::flush;
	A.min(scalar_1);
	dz << scalar_1 << std::endl << std::flush;
	dz << A.min() << std::endl << std::flush;
	dz << min(A) << std::endl << std::flush;
	std::cout << "...............................completed\n" << std::flush;

	std::cout << "testing max_pos / min_pos " << std::flush;
	dz << "testing min_pos / max_pos" << std::endl << std::flush;
	A.max_pos(scalar_1, i, j);
	dz << scalar_1 << " " << i << " " << j << std::endl << std::flush;
	dz << A.max_pos(i, j) << " " << i << " " << j << std::endl << std::flush;
	dz << max_pos(A, i, j) << " " << i << " " << j << std::endl << std::flush;
	A.min_pos(scalar_1, i, j);
	dz << scalar_1 << " " << i << " " << j << std::endl << std::flush;
	dz << A.min_pos(i, j) << " " << i << " " << j << std::endl << std::flush;
	dz << min_pos(A, i, j) << " " << i << " " << j << std::endl << std::flush;
	std::cout << ".......................completed\n" << std::flush;

	std::cout << "testing max_abs / min_abs " << std::flush;
	dz << "testing min_abs / max_abs" << std::endl << std::flush;
	A.max_abs(scalar_1);
	dz << scalar_1 << std::endl << std::flush;
	dz << A.max_abs() << std::endl << std::flush;
	dz << max_abs(A) << std::endl << std::flush;
	A.min_abs(scalar_1);
	dz << scalar_1 << std::endl << std::flush;
	dz << A.min_abs() << std::endl << std::flush;
	dz << min_abs(A) << std::endl << std::flush;
	std::cout << ".......................completed\n" << std::flush;

	std::cout << "testing max_abs_pos / min_abs_pos " << std::flush;
	dz << "testing min_abs_pos / max_abs_pos" << std::endl << std::flush;
	A.max_abs_pos(scalar_1, i, j);
	dz << scalar_1 << " " << i << " " << j << std::endl << std::flush;
	dz << A.max_abs_pos(i, j) << " " << i << " " << j << std::endl << std::flush;
	dz << max_abs_pos(A, i, j) << " " << i << " " << j << std::endl << std::flush;
	A.min_abs_pos(scalar_1, i, j);
	dz << scalar_1 << " " << i << " " << j << std::endl << std::flush;
	dz << A.min_abs_pos(i, j) << " " << i << " " << j << std::endl << std::flush;
	dz << min_abs_pos(A, i, j) << " " << i << " " << j << std::endl << std::flush;
	std::cout << "...............completed\n" << std::flush;

	std::cout << "testing hadamard " << std::flush;
	dz << "testing hadamard" << std::endl << std::flush;
	A.hadamard(scalar_1);
	dz << scalar_1 << std::endl << std::flush;
	dz << A.hadamard() << std::endl << std::flush;
	dz << hadamard(A) << std::endl << std::flush;
	std::cout << "................................completed\n" << std::flush;

	std::cout << "testing row_norm / column_norm " << std::flush;
	dz << "testing column_norm / row_norm" << std::endl << std::flush;
	A.row_norm(scalar_1, 0, 2);
	dz << scalar_1 << std::endl << std::flush;
	dz << A.row_norm(0, 2) << std::endl << std::flush;
	dz << row_norm(A, 0, 2) << std::endl << std::flush;
	A.column_norm(scalar_1, 0, 2);
	dz << scalar_1 << std::endl << std::flush;
	dz << A.column_norm(0, 2) << std::endl << std::flush;
	dz << column_norm(A, 0, 2) << std::endl << std::flush;
	std::cout << "..................completed\n" << std::flush;

	//std::cout << A << B << C << D << E << F << std::flush;

	std::cout << "testing rank" << std::flush;
	dz << "testing rank" << std::endl << std::flush;
	dz << rank(A) << std::endl << std::flush;
	dz << A.rank() << std::endl << std::flush;
	std::cout << ".....................................completed\n" << std::flush;

	std::cout << "testing lininr" << std::flush;
	dz << "testing lininr" << std::endl << std::flush;
	lidia_size_t *q = lininr(F);
	for (i = 0; i <= q[0]; i++)
		dz << q[i] << " ";
	dz << std::endl << std::flush;
	delete[] q;

	q = F.lininr();
	for (i = 0; i <= q[0]; i++)
		dz << q[i] << " ";
	dz << std::endl << std::flush;
	delete[] q;
	std::cout << "...................................completed\n" << std::flush;

	std::cout << "testing lininc" << std::flush;
	dz << "testing lininc" << std::endl << std::flush;
	lidia_size_t *q1 = lininc(F);
	for (i = 0; i <= q1[0]; i++)
		dz << q1[i] << " ";
	dz << std::endl << std::flush;
	delete[] q1;

	q1 = F.lininc();
	for (i = 0; i <= q1[0]; i++)
		dz << q1[i] << " ";
	dz << std::endl << std::flush;
	delete[] q1;
	std::cout << "...................................completed\n" << std::flush;

	std::cout << "testing RegExpansion" << std::flush;
	dz << "testing RegExpansion" << std::endl << std::flush;
	F = trans(A);
	q = lininr(F);
	D = A;
	E = A;
	regexpansion(D, q);
	dz << D << std::flush;
	E.regexpansion(q);
	dz << E << std::flush;
	std::cout << ".............................completed\n" << std::flush;

	std::cout << "testing adj" << std::flush;
	dz << "testing adj" << std::endl << std::flush;
#ifdef _MSC_VER
	F = A * bigint_matrix(trans(A));
#else
	F = A * trans(A);
#endif
	bigint_matrix ADJ = adj(F);
	dz << ADJ << ADJ * F << std::flush;
	ADJ.adj(F);
	dz << ADJ << ADJ * F << std::flush;
	std::cout << "......................................completed\n" << std::flush;

	std::cout << "testing latticedet / latticedet2" << std::flush;
	dz << "testing latticedet / latticedet2" << std::endl << std::flush;
	A.latticedet(scalar_1);
	dz << scalar_1 << std::endl << std::flush;
	scalar_1 = A.latticedet();
	dz << scalar_1 << std::endl << std::flush;
	scalar_1 = latticedet(A);
	dz << scalar_1 << std::endl << std::flush;

	A.latticedet2(scalar_1);
	dz << scalar_1 << std::endl << std::flush;
	A.latticedet2(scalar_1);
	dz << scalar_1 << std::endl << std::flush;
	A.latticedet2(scalar_1);
	dz << scalar_1 << std::endl << std::flush;
	std::cout << ".................completed\n" << std::flush;

	std::cout << "testing det" << std::flush;
	dz << "testing det" << std::endl << std::flush;
	F.det(scalar_1);
	dz << scalar_1 << std::endl << std::flush;
	scalar_1 = F.det();
	dz << scalar_1 << std::endl << std::flush;
	scalar_1 = det(F);
	dz << scalar_1 << std::endl << std::flush;
	std::cout << "......................................completed\n" << std::flush;

	std::cout << "testing charpoly" << std::flush;
	dz << "charpoly" << std::endl << std::flush;
	c = F.get_no_of_columns();
	pointer_5 = charpoly(F);
	for (i = 0; i < c + 1; i++)
		dz << pointer_5[i] << " ";
	dz << std::endl << std::flush;

	pointer_5 = F.charpoly();
	for (i = 0; i < c + 1; i++)
		dz << pointer_5[i] << " ";
	dz << std::endl << std::flush;
	std::cout << ".................................completed\n" << std::flush;

	std::cout << "testing hnf_havas(0) " << std::flush;
	dz << "testing hnf_havas(0) " << std::endl << std::flush;
	F = A;
	F.hnf_havas(0);
	dz << F << std::flush;
	F = A;
	F.hnf_havas(D, 0);
	dz << F << A*D << std::flush;

	F = A;
	F.hnf_havas(0);
	dz << F << std::flush;
	F = A;
	F.hnf_havas(D, 0);
	dz << F << A*D << std::flush;
	std::cout << "...............................completed\n" << std::flush;

	std::cout << "testing hnf_havas" << std::flush;
	dz << "testing hnf_havas" << std::endl << std::flush;
	F = A;
	F.hnf_havas(1);
	dz << F << std::flush;
	F = A;
	F.hnf_havas(D);
	dz << F << A*D << std::flush;

	F = A;
	F.hnf_havas(1);
	dz << F << std::flush;
	F = A;
	F.hnf_havas(D, 1);
	dz << F << A*D << std::flush;
	std::cout << "................................completed\n" << std::flush;

	std::cout << "testing hnfmod_dkt" << std::flush;
	dz << "testing hnfmod_dkt" << std::endl << std::flush;
	scalar_1 = latticedet(A);
	F = A;
	F.hnfmod_dkt(scalar_1);
	dz << F << std::flush;
	F = A;
	F.hnfmod_dkt();
	dz << F << std::flush;
	F = A;
	F.hnfmod_dkt();
	dz << F << std::flush;
	F = A;
	F.hnfmod_dkt(scalar_1);
	dz << F << std::flush;
	std::cout << "...............................completed\n" << std::flush;

	std::cout << "testing hnfmod_cohen" << std::flush;
	dz << "testing hnfmod_cohen" << std::endl << std::flush;
	scalar_1 = latticedet(A);
	F = A;
	F.hnfmod_cohen(scalar_1);
	dz << F << std::flush;
	F = A;
	F.hnfmod_cohen();
	dz << F << std::flush;
	F = A;
	F.hnfmod_cohen();
	dz << F << std::flush;
	F = A;
	F.hnfmod_cohen(scalar_1);
	dz << F << std::flush;
	std::cout << ".............................completed\n" << std::flush;

	std::cout << "testing hnfmod_mueller" << std::flush;
	dz << "testing hnfmod_mueller" << std::endl << std::flush;
	F = A;
	F.hnfmod_mueller(D);
	dz << F << A*D << std::flush;

	F = A;
	F.hnfmod_mueller(D);
	dz << F << A*D << std::flush;
	std::cout << "...........................completed\n" << std::flush;

	std::cout << "testing kernel / kernel2" << std::flush;
	dz << "testing kernel / kernel2" << std::endl << std::flush;
	F = kernel(A);
	dz << F << std::flush;
	D.kernel(A);
	dz << D << std::flush;

	F.kernel(A);
	dz << F << std::flush;
	D.kernel(A);
	dz << D << std::flush;
	std::cout << ".........................completed\n" << std::flush;

	std::cout << "testing reginvimage / reginvimage2" << std::flush;
	dz << "testing reginvimage / reginvimage2" << std::endl << std::flush;
	C = A * trans(A);
	F = reginvimage(C, B);
	dz << F << std::flush;
	D.reginvimage(C, B);
	dz << D << std::flush;

	F.reginvimage(C, B);
	dz << F << std::flush;
	D.reginvimage(C, B);
	dz << D << std::flush;
	std::cout << "...............completed\n" << std::flush;

	std::cout << "testing image" << std::flush;
	dz << "testing image" << std::endl << std::flush;
	F = image(A);
	dz << F << std::flush;
	D.image(A);
	dz << D << std::flush;

	F.image2(A);
	dz << F << std::flush;
	D.image2(A);
	dz << D << std::flush;
	std::cout << "....................................completed\n" << std::flush;

	std::cout << "testing invimage" << std::flush;
	dz << "testing invimage" << std::endl << std::flush;
	pointer_1 = new bigint[4];
	pointer_1[0] = 30;
	pointer_1[1] = 40;
	pointer_1[2] = 50;
	pointer_1[3] = 60;
	F = invimage(A, pointer_1);
	dz << F << std::flush;
	D.invimage(A, pointer_1);
	dz << D << std::flush;
	std::cout << ".................................completed\n" << std::flush;

	std::cout << "testing solve" << std::flush;
	dz << "testing solve" << std::endl << std::flush;
	pointer_1[0] = 30;
	pointer_1[1] = 40;
	pointer_1[2] = 50;
	pointer_1[3] = 60;
	F = solve(A, pointer_1);
	dz << F << std::flush;
	D.solve(A, pointer_1);
	dz << D << std::flush;
	std::cout << "....................................completed\n" << std::flush;

	std::cout << "testing snf_hartley" << std::flush;
	dz << "testing snf_hartley" << std::endl << std::flush;
	F = A;
	F.snf_hartley();
	dz << F << std::flush;
	F = A;
	F.snf_hartley(D, E);
	dz << F << D * A * E << std::flush;
	F = A;
	F.snf_hartley();
	dz << F << std::flush;
	F = A;
	F.snf_hartley(D, E);
	dz << F << D * A * E << std::flush;
	std::cout << "..............................completed\n" << std::flush;

	std::cout << "testing snf_simple" << std::flush;
	dz << "testing snf_simple" << std::endl << std::flush;
	F = A;
	F.snf_simple();
	dz << F << std::flush;
	F = A;
	F.snf_simple(D, E);
	dz << F << D * A * E << std::flush;
	F = A;
	F.snf_simple();
	dz << F << std::flush;
	F = A;
	F.snf_simple(D, E);
	dz << F << D * A * E << std::flush;
	std::cout << "...............................completed\n" << std::flush;

	std::cout << "testing snf_havas" << std::flush;
	dz << "testing snf_havas" << std::endl << std::flush;
	F = A;
	F.snf_havas();
	dz << F << std::flush;
	F = A;
	F.snf_havas(D, E);
	dz << F << D * A * E << std::flush;
	F = A;
	F.snf_havas();
	dz << F << std::flush;
	F = A;
	F.snf_havas(D, E);
	dz << F << D * A * E << std::flush;
	std::cout << "................................completed\n" << std::flush;

	std::cout << "testing snf_mult" << std::flush;
	dz << "testing snf_mult" << std::endl << std::flush;
	F = A;
	F.snf_mult();
	dz << F << std::flush;
	F = A;
	F.snf_mult(D, E);
	dz << F << D * A * E << std::flush;
	F = A;
	F.snf_mult();
	dz << F << std::flush;
	F = A;
	F.snf_mult(D, E);
	dz << F << D * A * E << std::flush;
	std::cout << ".................................completed\n" << std::flush;

	std::cout << "testing snf_add" << std::flush;
	dz << "testing snf_add" << std::endl << std::flush;
	F = A;
	F.snf_add();
	dz << F << std::flush;
	F = A;
	F.snf_add(D, E);
	dz << F << D * A * E << std::flush;
	F = A;
	F.snf_add();
	dz << F << std::flush;
	F = A;
	F.snf_add(D, E);
	dz << F << D * A * E << std::flush;
	std::cout << "..................................completed\n" << std::flush;

	std::cout << "testing snf_new" << std::flush;
	dz << "testing snf_new" << std::endl << std::flush;
	F = A;
	F.snf_new();
	dz << F << std::flush;
	F = A;
	F.snf_new(D, E);
	dz << F << D * A * E << std::flush;
	F = A;
	F.snf_new();
	dz << F << std::flush;
	F = A;
	F.snf_new(D, E);
	dz << F << D * A * E << std::flush;
	std::cout << "..................................completed\n" << std::flush;

	std::cout << "testing snfmod_dkt" << std::flush;
	dz << "testing snfmod_dkt" << std::endl << std::flush;
	scalar_1 = latticedet(A);
	F = A;
	F.snfmod_dkt(scalar_1);
	dz << F << std::flush;
	F = A;
	F.snfmod_dkt();
	dz << F << std::flush;
	F = A;
	F.snfmod_dkt();
	dz << F << std::flush;
	F = A;
	F.snfmod_dkt(scalar_1);
	dz << F << std::flush;
	std::cout << "...............................completed\n" << std::flush;

	std::cout << "testing snfmod_cohen" << std::flush;
	dz << "testing snfmod_cohen" << std::endl << std::flush;
	scalar_1 = latticedet(A);
	F = A;
	F.snfmod_cohen(scalar_1);
	dz << F << std::flush;
	F = A;
	F.snfmod_cohen();
	dz << F << std::flush;
	F = A;
	F.snfmod_cohen();
	dz << F << std::flush;
	F = A;
	F.snfmod_cohen(scalar_1);
	dz << F << std::flush;
	std::cout << ".............................completed\n" << std::flush;

	std::cout << "\nThe functions ChinRest, get_primes, mgcd, mgcd1, mgcd2 and mgcd3 are\n"
		" tested by implication." << std::flush;

	dz.close();
	std::cout << "\nPlease use now the 'diff' or 'cmp' command to verify\n";
	std::cout << "the equality of bigint_matrix_appl.dat";
	std::cout << " and bigint_matrix_appl.out.\n\n" << std::flush;

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
