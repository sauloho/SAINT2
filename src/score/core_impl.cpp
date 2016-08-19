#include <cstdlib>
#include <iostream>
#include <cstdio>
#include <cassert>
#include <cmath>

#include "peptide.h"
#include "residue.h"
#include "atom.h"
#include "amino.h"
#include "c_file.h"
#include "scorer_combined.h"
#include "core.h"
#include "core_impl.h"


#ifndef M_SQRT1_2
#define M_SQRT1_2 0.70710678119
#endif

CORE_impl::CORE_impl()
	: m_data_loaded(false)
{
	m_data = NULL;
}

CORE_impl::~CORE_impl()
{
	if (m_data != NULL)
	{
		for (int b1 = 0;b1 < NUM_RAPDF_IDS;b1++)
		{
			for (int b2 = 0;b2 <= b1;b2++)
			{
				delete [] m_data[b1][b2];
			}

			delete [] m_data[b1];
		}

		delete [] m_data;
	}
}

CORE_impl::CORE_Params CORE_impl::m_lj[Num_Backbone][Num_Backbone];

CORE_impl::CORE_Params::CORE_Params()
{
}

CORE_impl::CORE_Params::CORE_Params(double c12_val, double c6_val)
	: c12(c12_val), c6(c6_val)
{
	// sigma is the distance at which the potential equals zero
	double sigma = pow(c12 / c6, 1.0 / 6.0);

	// recommended cutoff distance
	// (http://en.wikipedia.org/wiki/Lennard-Jones_potential)
	//
	// This causes a small discontinuity that is probably not worth
	// fixing for our purposes (the jump is about 1.5 percent of the
	// well depth).
	this->max_dist = 2.5 * sigma;
}

// return a single letter representing an atom type
// (backbone atoms only)

inline char atom_char1(Atom_Id a)
{
	switch (a)
	{
		case Atom_C:
		case Atom_CA:
		case Atom_CB:
			return 'C';

		case Atom_N:
			return 'N';

		case Atom_O:
			return 'O';

		default:
			std::cerr << "Error: Illegal value in atom_char1() in "
				<< __FILE__ << "\n";
			exit(1);
	}
}

class CORE_Static_Init
{
public:
	// Values taken from:
	//
	// http://web.njit.edu/all_topics/Prog_Lang_Docs/html/autodock/
	//		AD3.a.0UserGuide.html#29613
	//
	// (Table 1, "Self-consistent Lennard-Jones 12-6 parameters")

	CORE_Static_Init()
	{
		// parameters taken from
		// http://dasher.wustl.edu/tinker/distribution/params/oplsaa.prm

		const double cc_sigma = 3.75;
		const double cc_epsilon = 0.105;

		const double oo_sigma = 2.96;
		const double oo_epsilon = 0.21;

		const double nn_sigma = 3.25;
		const double nn_epsilon = 0.17;

		for (int a1 = 0;a1 < Num_Backbone;a1++)
		{
			for (int a2 = 0;a2 < Num_Backbone;a2++)
			{
				char c1 = atom_char1((Atom_Id) a1);
				char c2 = atom_char1((Atom_Id) a2);
				double e, s;

				if (c1 > c2)
				{
					swap(c1, c2);
				}

				if (c1 == 'C' && c2 == 'C')
				{
					e = cc_epsilon;
					s = cc_sigma;
				}
				else
				if (c1 == 'C' && c2 == 'N')
				{
					e = sqrt(cc_epsilon * nn_epsilon);
					s = 0.5 * (cc_sigma + nn_sigma);
				}
				else
				if (c1 == 'C' && c2 == 'O')
				{
					e = sqrt(cc_epsilon * oo_epsilon);
					s = 0.5 * (cc_sigma + oo_sigma);
				}
				else
				if (c1 == 'N' && c2 == 'N')
				{
					e = nn_epsilon;
					s = nn_sigma;
				}
				else
				if (c1 == 'N' && c2 == 'O')
				{
					e = sqrt(nn_epsilon * oo_epsilon);
					s = 0.5 * (nn_sigma + oo_sigma);
				}
				else
				if (c1 == 'O' && c2 == 'O')
				{
					e = oo_epsilon;
					s = oo_sigma;
				}
				else
				{
					std::cerr << "Error in CORE_Static_Init()\n";
					exit(1);
				}

				double c12 = e * pow(s, 12.0);
				double c6 = 2.0 * e * pow(s, 6.0);

				CORE_impl::m_lj[a1][a2] = CORE_impl::CORE_Params(c12, c6);

				// check if minimum is where it should be: at (s, -e)
				if (fabs(c12 / pow(s, 12.0) - c6 / pow(s, 6.0) + e) > 0.001)
				{
					std::cerr << "Incorrect values in CORE_Static_Init()\n";
					exit(1);
				}
			}
		}
	}
};

namespace
{
	CORE_Static_Init core_static_init;
}

void CORE_impl::set_data_file(const std::string &filename)
{
	m_filename = filename;
}

void CORE_impl::set_data_file_ori(const std::string &filename)
{
	m_filename_ori = filename;
}

void CORE_impl::load_data()
{
	if (m_data_loaded)
	{
		return;
	}

	if (m_filename.empty())
	{
		std::cerr << "Error: RAPDF data file undefined\n";
		exit(1);
	}

//std::cout << "LOADING RAPDF DATA\n";

	// exits with error message if file does not exist
	C_File file(m_filename, "r", "RAPDF data file");

	// std::cout << "Loading RAPDF data" << std::endl;

	static const int Max_Len = 1000;
	char buffer[Max_Len];

	if (!file.next_line(buffer, Max_Len) ||
		sscanf(buffer, "%d %d", &m_first_bin, &m_top_bin) != 2)
	{
		std::cerr << "Error: expected RAPDF data file "
			<< m_filename
			<< " to start with two numbers (first bin, top bin)\n";
		exit(1);
	}

	int b1, b2, d;

	m_data = new double** [NUM_RAPDF_IDS];

	for (b1 = 0;b1 < NUM_RAPDF_IDS;b1++)
	{
		m_data[b1] = new double* [NUM_RAPDF_IDS];

		for (b2 = 0;b2 < NUM_RAPDF_IDS;b2++)
		{
			m_data[b1][b2] = new double[m_top_bin];
		}
	}

	m_rapdf_ids = new int* [100];

	for (b1 = 0;b1 < 100;b1++)
		m_rapdf_ids[b1] = new int [NUM_RAPDF_IDS];

	load_rapdf_ids();

	for (b1 = 0;b1 < NUM_RAPDF_IDS;b1++)
	{
		std::string b1_aa_str = Amino::rapdf_amino(b1).abbr();
		std::string b1_a_str = Amino::rapdf_atom(b1).name();

		for (b2 = 0;b2 <= b1;b2++)
		{
			std::string b2_aa_str = Amino::rapdf_amino(b2).abbr();
			std::string b2_a_str = Amino::rapdf_atom(b2).name();

			char aa1_str[100], a1_str[100];
			char aa2_str[100], a2_str[100];

			if (!file.next_line(buffer, Max_Len) ||
				sscanf(buffer, "%s %s %s %s",
				aa1_str, a1_str, aa2_str, a2_str) != 4 ||
				std::string(aa1_str) != b1_aa_str ||
				std::string(a1_str) != b1_a_str ||
				std::string(aa2_str) != b2_aa_str ||
				std::string(a2_str) != b2_a_str)
			{
				std::cerr << "Error on line " << file.line_num()
					<< " of " << m_filename
					<< ": expected \""
					<< b1_aa_str << " " << b1_a_str << "  "
					<< b2_aa_str << " " << b2_a_str << "\"\n";
				exit(1);
			}

			int dist;
			double val;

			for (d = m_first_bin;d < m_top_bin;d++)
			{
				if (!file.next_line(buffer, Max_Len) ||
					sscanf(buffer, "%d %lf", &dist, &val) != 2 ||
					dist != d)
				{
					std::cerr << "Error on line " << file.line_num()
						<< " of " << m_filename
						<< ": expected \""
						<< d
						<< "\" followed by value\n";
					exit(1);
				}

				m_data[b1][b2][d] = m_data[b2][b1][d] = val;
			}
			for( d = 0 ; d < m_first_bin; d++)
			{
				m_data[b1][b2][d]=m_data[b1][b2][m_first_bin];
				m_data[b2][b1][d]=m_data[b2][b1][m_first_bin];
			}
		}
	}

	if (file.next_line(buffer, Max_Len))
	{
		std::cerr << "Warning: ignoring extra data on line "
			<< file.line_num()
			<< " of "
			<< m_filename
			<< "\n";
	}

	if (m_filename.empty())
	{
		std::cerr << "Error: Orientation data file undefined\n";
		exit(1);
	}

	//std::cout << "LOADING ORIENTATION DATA\n";
/*
	// exits with error message if the file does not exist
	C_File file_ori(m_filename_ori, "r", "Orientation data file");

	double val;

	for (int a1 = 0;a1 < Amino::Num;a1++)
	{
		for (int a2 = a1;a2 < Amino::Num;a2++)
		{
			if (!file_ori.next_line(buffer, Max_Len) || std::string(buffer, 3) != Amino(a1).abbr() || std::string(buffer + 4, 3) != Amino(a2).abbr())
			{
				std::cerr << "Error: expected \"" << Amino(a1).abbr() << ' ' << Amino(a2).abbr() << "\" on line " << file_ori.line_num()
					<< " of orientation data file " << m_filename_ori << "\n"; exit(1);
			}

			for (int angle = 0;angle < ORIENT_ANGLES;angle++)
			{
				for (int dist = 0;dist < ORIENT_DISTS;dist++)
				{
					int d;

					if (!file_ori.next_line(buffer, Max_Len) || sscanf(buffer, "%d %lf", &d, &val) != 2 || d != dist + 3)
					{
						std::cerr << "Error: expected \"" << dist + 3 << "\" followed by value on line " << file_ori.line_num() << " of orientation data file "
							<< m_filename_ori << "\n"; exit(1);
					}

//					std::cout << ": " << dist << " " << angle << " " << a1 << " " << a2 << " " << " = " << val << "\n";

					m_data_ori[dist][angle][a1][a2] = val;
					m_data_ori[dist][angle][a2][a1] = val;
				}
			}
		}
	}

	if (file.next_line(buffer, Max_Len))
		std::cerr << "Warning: ignoring extra data on line " << file_ori.line_num() << " of orientation data file " << m_filename_ori << "\n";
*/

	m_data_loaded = true;
}


void CORE_impl::load_rapdf_ids()
{
	// Cysteine
	
	m_rapdf_ids[4][Atom_N] = 0;	// 0
	m_rapdf_ids[4][Atom_CA] = 1;
	m_rapdf_ids[4][Atom_C] = 2;
	m_rapdf_ids[4][Atom_O] = 3;
	m_rapdf_ids[4][Atom_CB] = 4;
	m_rapdf_ids[4][Atom_SG] = 5;
	
	// Glutamine
	
	m_rapdf_ids[6][Atom_N] = 6;	// 6
	m_rapdf_ids[6][Atom_CA] = 7;
	m_rapdf_ids[6][Atom_C] = 8;
	m_rapdf_ids[6][Atom_O] = 9;
	m_rapdf_ids[6][Atom_CB] = 10;
	m_rapdf_ids[6][Atom_CG] = 11;
	m_rapdf_ids[6][Atom_CD] = 12;
	m_rapdf_ids[6][Atom_OE1] = 13;
	m_rapdf_ids[6][Atom_NE2] = 14;
	
	// Aspartic acid
	
	m_rapdf_ids[3][Atom_N] = 15;	// 15
	m_rapdf_ids[3][Atom_CA] = 16;
	m_rapdf_ids[3][Atom_C] = 17;
	m_rapdf_ids[3][Atom_O] = 18;
	m_rapdf_ids[3][Atom_CB] = 19;
	m_rapdf_ids[3][Atom_CG] = 20;
	m_rapdf_ids[3][Atom_OD1] = 21;
	m_rapdf_ids[3][Atom_OD2] = 22;
	
	// Serine
	
	m_rapdf_ids[15][Atom_N] = 23;	// 23
	m_rapdf_ids[15][Atom_CA] = 24;
	m_rapdf_ids[15][Atom_C] = 25;
	m_rapdf_ids[15][Atom_O] = 26;
	m_rapdf_ids[15][Atom_CB] = 27;
	m_rapdf_ids[15][Atom_OG] = 28;
	
	// Valine
	
	m_rapdf_ids[19][Atom_N] = 29;	// 29
	m_rapdf_ids[19][Atom_CA] = 30;
	m_rapdf_ids[19][Atom_C] = 31;
	m_rapdf_ids[19][Atom_O] = 32;
	m_rapdf_ids[19][Atom_CB] = 33;
	m_rapdf_ids[19][Atom_CG1] = 34;
	m_rapdf_ids[19][Atom_CG2] = 35;
	
	// Methionine
	
	m_rapdf_ids[12][Atom_N] = 36;	// 36
	m_rapdf_ids[12][Atom_CA] = 37;
	m_rapdf_ids[12][Atom_C] = 38;
	m_rapdf_ids[12][Atom_O] = 39;
	m_rapdf_ids[12][Atom_CB] = 40;
	m_rapdf_ids[12][Atom_CG] = 41;
	m_rapdf_ids[12][Atom_SD] = 42;
	m_rapdf_ids[12][Atom_CE] = 43;
	
	// Proline
	
	m_rapdf_ids[14][Atom_N] = 44;	// 44
	m_rapdf_ids[14][Atom_CA] = 45;
	m_rapdf_ids[14][Atom_C] = 46;
	m_rapdf_ids[14][Atom_O] = 47;
	m_rapdf_ids[14][Atom_CB] = 48;
	m_rapdf_ids[14][Atom_CG] = 49;
	m_rapdf_ids[14][Atom_CD] = 50;
	
	// Lysine
	
	m_rapdf_ids[11][Atom_N] = 51;	// 51
	m_rapdf_ids[11][Atom_CA] = 52;
	m_rapdf_ids[11][Atom_C] = 53;
	m_rapdf_ids[11][Atom_O] = 54;
	m_rapdf_ids[11][Atom_CB] = 55;
	m_rapdf_ids[11][Atom_CG] = 56;
	m_rapdf_ids[11][Atom_CD] = 57;
	m_rapdf_ids[11][Atom_CE] = 58;
	m_rapdf_ids[11][Atom_NZ] = 59;
	
	// Threonine
	
	m_rapdf_ids[16][Atom_N] = 60;	// 60
	m_rapdf_ids[16][Atom_CA] = 61;
	m_rapdf_ids[16][Atom_C] = 62;
	m_rapdf_ids[16][Atom_O] = 63;
	m_rapdf_ids[16][Atom_CB] = 64;
	m_rapdf_ids[16][Atom_OG1] = 65;
	m_rapdf_ids[16][Atom_CG2] = 66;
	
	// Phenylalanine
	
	m_rapdf_ids[13][Atom_N] = 67;	// 67
	m_rapdf_ids[13][Atom_CA] = 68;
	m_rapdf_ids[13][Atom_C] = 69;
	m_rapdf_ids[13][Atom_O] = 70;
	m_rapdf_ids[13][Atom_CB] = 71;
	m_rapdf_ids[13][Atom_CG] = 72;
	m_rapdf_ids[13][Atom_CD1] = 73;
	m_rapdf_ids[13][Atom_CD2] = 74;
	m_rapdf_ids[13][Atom_CE1] = 75;
	m_rapdf_ids[13][Atom_CE2] = 76;
	m_rapdf_ids[13][Atom_CZ] = 77;
	
	// Alanine
	
	m_rapdf_ids[0][Atom_N] = 78;	// 78
	m_rapdf_ids[0][Atom_CA] = 79;
	m_rapdf_ids[0][Atom_C] = 80;
	m_rapdf_ids[0][Atom_O] = 81;
	m_rapdf_ids[0][Atom_CB] = 82;
	
      // Histidine
        
      m_rapdf_ids[8][Atom_N] = 83;        // 83
      m_rapdf_ids[8][Atom_CA] = 84;
      m_rapdf_ids[8][Atom_C] = 85;
      m_rapdf_ids[8][Atom_O] = 86;
      m_rapdf_ids[8][Atom_CB] = 87;
      m_rapdf_ids[8][Atom_CG] = 88;
      m_rapdf_ids[8][Atom_ND1] = 89;
      m_rapdf_ids[8][Atom_CD2] = 90;
      m_rapdf_ids[8][Atom_CE1] = 91;
      m_rapdf_ids[8][Atom_NE2] = 92;


	// Glycine
	
	m_rapdf_ids[7][Atom_N] = 93;	// 93
	m_rapdf_ids[7][Atom_CA] = 94;
	m_rapdf_ids[7][Atom_C] = 95;
	m_rapdf_ids[7][Atom_O] = 96;
	
	// Isoleucine
	
	m_rapdf_ids[9][Atom_N] = 97;	// 97
	m_rapdf_ids[9][Atom_CA] = 98;
	m_rapdf_ids[9][Atom_C] = 99;
	m_rapdf_ids[9][Atom_O] = 100;
	m_rapdf_ids[9][Atom_CB] = 101;
	m_rapdf_ids[9][Atom_CG1] = 102;
	m_rapdf_ids[9][Atom_CG2] = 103;
	m_rapdf_ids[9][Atom_CD1] = 104;
	
	// Glutamic acid
	
	m_rapdf_ids[5][Atom_N] = 105;	// 105
	m_rapdf_ids[5][Atom_CA] = 106;
	m_rapdf_ids[5][Atom_C] = 107;
	m_rapdf_ids[5][Atom_O] = 108;
	m_rapdf_ids[5][Atom_CB] = 109;
	m_rapdf_ids[5][Atom_CG] = 110;
	m_rapdf_ids[5][Atom_CD] = 111;
	m_rapdf_ids[5][Atom_OE1] = 112;
	m_rapdf_ids[5][Atom_OE2] = 113;
	
	// Leucine
	
	m_rapdf_ids[10][Atom_N] = 114;	// 114
	m_rapdf_ids[10][Atom_CA] = 115;
	m_rapdf_ids[10][Atom_C] = 116;
	m_rapdf_ids[10][Atom_O] = 117;
	m_rapdf_ids[10][Atom_CB] = 118;
	m_rapdf_ids[10][Atom_CG] = 119;
	m_rapdf_ids[10][Atom_CD1] = 120;
	m_rapdf_ids[10][Atom_CD2] = 121;
	
	// Arginine
	
	m_rapdf_ids[1][Atom_N] = 122;	// 122
	m_rapdf_ids[1][Atom_CA] = 123;
	m_rapdf_ids[1][Atom_C] = 124;
	m_rapdf_ids[1][Atom_O] = 125;
	m_rapdf_ids[1][Atom_CB] = 126;
	m_rapdf_ids[1][Atom_CG] = 127;
	m_rapdf_ids[1][Atom_CD] = 128;
	m_rapdf_ids[1][Atom_NE] = 129;
	m_rapdf_ids[1][Atom_CZ] = 130;
	m_rapdf_ids[1][Atom_NH1] = 131;
	m_rapdf_ids[1][Atom_NH2] = 132;
	
	// Tryptophan
	
	m_rapdf_ids[17][Atom_N] = 133;	// 133
	m_rapdf_ids[17][Atom_CA] = 134;
	m_rapdf_ids[17][Atom_C] = 135;
	m_rapdf_ids[17][Atom_O] = 136;
	m_rapdf_ids[17][Atom_CB] = 137;
	m_rapdf_ids[17][Atom_CG] = 138;
	m_rapdf_ids[17][Atom_CD1] = 139;
	m_rapdf_ids[17][Atom_CD2] = 140;
	m_rapdf_ids[17][Atom_NE1] = 141;
	m_rapdf_ids[17][Atom_CE2] = 142;
	m_rapdf_ids[17][Atom_CE3] = 143;
	m_rapdf_ids[17][Atom_CZ2] = 144;
	m_rapdf_ids[17][Atom_CZ3] = 145;
	m_rapdf_ids[17][Atom_CH2] = 146;
	
	// Asparagine
	
	m_rapdf_ids[2][Atom_N] = 147;	// 147
	m_rapdf_ids[2][Atom_CA] = 148;
	m_rapdf_ids[2][Atom_C] = 149;
	m_rapdf_ids[2][Atom_O] = 150;
	m_rapdf_ids[2][Atom_CB] = 151;
	m_rapdf_ids[2][Atom_CG] = 152;
	m_rapdf_ids[2][Atom_OD1] = 153;
	m_rapdf_ids[2][Atom_ND2] = 154;
	
	// Tyrosine
	
	m_rapdf_ids[18][Atom_N] = 155;	// 155
	m_rapdf_ids[18][Atom_CA] = 156;
	m_rapdf_ids[18][Atom_C] = 157;
	m_rapdf_ids[18][Atom_O] = 158;
	m_rapdf_ids[18][Atom_CB] = 159;
	m_rapdf_ids[18][Atom_CG] = 160;
	m_rapdf_ids[18][Atom_CD1] = 161;
	m_rapdf_ids[18][Atom_CD2] = 162;
	m_rapdf_ids[18][Atom_CE1] = 163;
	m_rapdf_ids[18][Atom_CE2] = 164;
	m_rapdf_ids[18][Atom_CZ] = 165;
	m_rapdf_ids[18][Atom_OH] = 166;

	// end marker

	m_rapdf_ids[26][Atom_C] = 167;	// 167
}

double CORE_impl::score(const Peptide &p,double w_LJ, double w_RAPDF, bool verbose)
{
	// (does nothing if already loaded)
	load_data();

	double total_RAPDF = 0.0;
	double total_LJ = 0.0, d6;
	double d_dist;//, frac, weight;
	int n1,n2,a1,a2;
	int id_1,id_2;
//	int dist, other_bin;

	Atom_Id t1,t2;
	Point pos1, pos2;

	for (n1 = p.start() + 2;n1 <= p.end();n1++)
	{
		const Residue &res1 = p.res(n1);

		for (a1 = 0;a1 < res1.num_atoms();a1++)
		{
			if (!p.atom_exists(n1, (Atom_Id) a1)) continue;

			t1 = res1.m_atom[a1].m_type.m_type;
//			if (t1 == Atom_Undef)				continue;

			id_1 = m_rapdf_ids[res1.m_amino.m_val][t1]; //res1.amino().rapdf_id(t1);
			pos1 = p.m_conf.m_backbone[n1 * Num_Backbone + t1];//atom_pos(n1, t1);
			
			for (n2 = p.start();n2 < n1 - 1; n2++)
			{
				const Residue &res2 = p.res(n2);

				for (a2 = 0;a2 < res2.num_atoms();a2++)
				{
					if (!p.atom_exists(n2, (Atom_Id) a2)) continue;
					t2 = res2.m_atom[a2].m_type.m_type;
//					if (t2 == Atom_Undef)	continue;

					id_2 = m_rapdf_ids[res2.m_amino.m_val][t2];   // res2.amino().rapdf_id(t2);	
					pos2 = p.m_conf.m_backbone[n2 * Num_Backbone + t2]; // pos2(n2, t2);

					const CORE_Params &lj = m_lj[a1][a2];
//					p2 = p.atom_pos2(n2, (Atom_Id) a2);

					d_dist = pos1.distance(pos2);

					/* RAPDF Potential */
					if ( m_top_bin - d_dist > 0 )
					{
//						dist = (int) d_dist;
//						frac = d_dist - (double) dist;
//						weight=fabs(frac-0.5);
//						other_bin = (frac > 0.5 ? dist + 1 : dist - 1);
//						if (other_bin < m_first_bin || other_bin >= m_top_bin)
						total_RAPDF += m_data[id_1][id_2][(int)d_dist];
//						else
//							total_RAPDF += (m_data[id_1][id_2][dist] * 0.5 + m_data[id_1][id_2][other_bin] * weight) / (weight+0.5); 
					}
					else
						total_RAPDF += m_data[id_1][id_2][m_top_bin - 1];

					/* Lennard-Jones Potential */
					if (fabs(pos1.x - pos2.x) < lj.max_dist && fabs(pos1.y - pos2.y) < lj.max_dist && fabs(pos1.z - pos2.z) < lj.max_dist )
					{
						if (d_dist < lj.max_dist)
						{
							if (d_dist < 1.0)
								total_LJ += lj.c12 - lj.c6;
							else
							{
								d6 = square(d_dist * d_dist * d_dist);
								total_LJ += (lj.c12 / square(d6)) - (lj.c6 / d6);
							}
						}
					}
				}
			}
		}
	}

#ifndef RAW_SCORE
	int len = p.length();
	// normalise so that all score types have approximately the same range
	total_LJ /= (len + 200.0);
	total_LJ *= 350.0;

	if (len <= SHORT_PEPTIDE)
	{
		total_RAPDF /= (double) len;
		// put on the same scale as the other scoring terms
		total_RAPDF *= 7.0;
	}
	else
	{
		total_RAPDF /= sqrt((double) len);
		// put on the same scale as the other scoring terms
		total_RAPDF *= 0.75;
	}
#endif //RAW_SCORE
//	std::cout << "LJ\t" << total_LJ << "\tRAPDF\t" << total_RAPDF << "\n";

	return w_LJ*total_LJ +w_RAPDF*total_RAPDF;
}

