#include <cstdlib> // PG added this
#include <iostream>
#include <cstring>
#include <cassert>
#include "config.h"
#include "common.h"
#include "c_file.h"
#include "stream_printf.h"
#include "codon.h"
#include "runner.h"
#include "sequence.h"
#include "extender.h"
#include "extender_codon.h"

// static data members
const char *Extender_Codon::c_type = "codon";
const char *Extender_Codon::c_param_file = "codonfile";
const char *Extender_Codon::c_param_tunnel_length = "tunnel_length";

Extender_Codon::Extender_Codon()
	: m_tunnel_length(0)
{
	for (int n = 0;n < Codon::Num;n++)
	{
		m_codon_factor[n] = 1.0;
	}
}

Extender_Codon::~Extender_Codon()
{
}

void Extender_Codon::set_filename(const std::string &filename)
{
	m_filename = filename;
}

bool Extender_Codon::parse_parameter(const std::string &name,
	const std::string &value)
{
	std::string full_name = Extender::config_section();
	full_name += " ";
	full_name += name;

	if (name == c_param_file)
	{
		set_filename(value);
		return true;
	}
	else
	if (name == c_param_tunnel_length)
	{
		m_tunnel_length = parse_integer(value, full_name, 0);
		return true;
	}

	return false;
}

void Extender_Codon::verify_parameters()
{
	if (m_filename.empty())
	{
		Config::missing_parameter(Extender::config_section(),
			c_type, c_param_file);
	}

	// Extender_Codon requires one residue to be extruded at a time
	if (residues_to_extrude() != 1)
	{
		std::cerr << Config::cmd()
			<< ": warning: "
			<< Extender::config_section()
			<< " type "
			<< c_type
			<< ": changing \""
			<< Extender::extrude_res_param()
			<< "\" parameter from "
			<< residues_to_extrude()
			<< " to 1\n";

		set_extrude_res(1);
	}

	read_file();
}

void Extender_Codon::read_file()
{
	// file format:
	//
	// 1,ATG,1
	// 2,CAG,0.33
	// 3,TTT,0.10975
	// ...

	int n;
	for (n = 0;n < Codon::Num;n++)
	{
		// default value if it doesn't appear in the file
		m_codon_factor[n] = 1.0;
	}

	C_File file(m_filename.c_str(), "r", "Codon file");

	char buffer[1000];

	while (file.next_line(buffer, 1000))
	{
		bool ok = false;
		char *p = strchr(buffer, ',');

		if (p != NULL)
		{
			Codon c(p[1], p[2], p[3]);
			double val;

			if (!c.illegal() && p[4] == ',' &&
				sscanf(p + 5, "%lf", &val) == 1 &&
				val > 0.0 && val <= 1.0)
			{
				m_codon_factor[c.num()] = 1.0 / val;
				ok = true;
			}
		}

		if (!ok)
		{
			std::cerr << "Error on line "
				<< file.line_num()
				<< " of codon file "
				<< m_filename
				<< " (expected format like: 10,AGA,0.15)\n";
			exit(1);
		}
	}
}

void Extender_Codon::calculate_num_moves(const Sequence &seq)
{
	assert(m_extrude_res == 1);

	int full = seq.length();
	m_moves.resize(full);

	for (int x = 0;x < full;x++)
	{
		m_moves[x] = -1;
	}

	double total_factor = 0.0;

	int i;
	for (i = m_initial_res;i < full;i++)
	{
		int j = i + m_tunnel_length;

		if (reverseSaint)
		{
			j = seq.length() - 1 - j;
		}

		double f = 1.0;

		// extrusions which happen after the peptide is released
		// from the ribosome are assigned a factor of 1.0, ie.
		// the fastest possible codon speed (smallest #moves)

		if (j < full && j >= 0)
		{
			int c = seq.codon(j).num();
			f = m_codon_factor[c];
// std::cout << "CALC MOVES: " << i
// << " (" << j << ") = " << 1.0/m_codon_factor[c] << "\n";
		}
// else { std::cout << "CALC MOVES: " << i << " = default\n"; }

		total_factor += f;
	}

	int total_moves = 0;

	for (i = m_initial_res;i < full;i++)
	{
		int j = i + m_tunnel_length;

		if (reverseSaint)
		{
			j = seq.length() - 1 - j;
		}

		double f = 1.0;

		if (j < full && j >= 0)
		{
			int c = seq.codon(i).num();
			f = m_codon_factor[c];
		}

		m_moves[i] = (int) (m_growth_moves * f / total_factor + 0.5);
		total_moves += m_moves[i];
	}

//std::cout << "CALC MOVES: initial total = " << total_moves << "\n";

	// adjust the total number of moves to exactly m_growth_moves

	while (total_moves != m_growth_moves)
	{
		if (total_moves < m_growth_moves)
		{
			// add moves from the end down

			for (i = full - 1;i >= m_initial_res;i--)
			{
//std::cout << "ADDING 1 TO " << i << "\n";
				++m_moves[i];
				if (++total_moves >= m_growth_moves) { break; }
			}
		}
		else   // total_moves > m_growth_moves
		{
			// remove moves from the start up (but not below 0)

			for (i = m_initial_res;i < full;i++)
			{
				if (m_moves[i] > 0)
				{
//std::cout << "TAKING 1 FROM " << i << "\n";
					--m_moves[i];
					if (--total_moves <= m_growth_moves) { break; }
				}
			}
		}
	}
}

void Extender_Codon::print_template(std::ostream &out,
	bool commented /*= true*/)
{
	const char *c = (commented ? "#" : "");

	out << c << "type = " << c_type << "\n"
		<< c << c_param_file << " = ...\t# TAI profile csv file\n"
		//<< c << c_param_multiplier << " = "
		//	<< Printf("%.1f", c_default_multiplier)
		//	<< "\t# multiplication factor (codon speed "
		//		"to #moves)\n"
		<< c << c_param_tunnel_length << " = 0"
			<< "\t# number of residues in tunnel (delay)\n"
		<< "\n";
}

