#include <cstdlib> // PG added this
#include <iostream>
#include "config.h"
#include "sequence.h"
#include "extender.h"
#include "extender_fixed.h"

// static data members
const char *Extender_Fixed::c_type = "fixed";
const char *Extender_Fixed::c_param_move_distrib = "move_distribution";
const char *Extender_Fixed::c_default_move_distrib = "fixed";

Extender_Fixed::Extender_Fixed()
	: m_move_distrib(Distrib_Fixed)
{
}

bool Extender_Fixed::parse_parameter(const std::string &name,
	const std::string &value)
{
	if (name == c_param_move_distrib)
	{
		if (value == "fixed")
		{
			m_move_distrib = Distrib_Fixed;
		}
		else if (value == "linear")
		{
			m_move_distrib = Distrib_Linear;
		}
		else
		{
			std::cerr << Config::cmd()
				<< ": illegal value \"" << value
				<< "\" for " << Extender::c_config_section << " parameter \""
				<< name << "\", type " << c_type << "\n";
			exit(1);
		}

		return true;
	}

	return false;
}

void Extender_Fixed::verify_parameters()
{
}

Extender_Fixed::~Extender_Fixed()
{
}

void Extender_Fixed::calculate_num_moves(const Sequence &seq)
{
	int full = seq.length();
	m_moves.resize(full);

	// set everything to -1 (undefined)

	for (int x = 0;x < full;x++)
	{
		m_moves[x] = -1;
	}

	int n = 0;	// total number of extrusions

	// slow but sure way to get the total number of
	// extrusions (only done once)

	for (int j = m_initial_res;j < full;j += m_extrude_res)
	{
		n++;
	}

	// number of extra moves to add at the end to bring
	// the total to exactly m_growth_moves

	int extra = 0;

	int num = m_growth_moves / n;
	extra = m_growth_moves - (num * n);

	for (int i = m_initial_res;i < full;i += m_extrude_res)
	{
		m_moves[i] = num;
	}
	
	if (m_move_distrib == Distrib_Linear)
	{
		int n1 = m_initial_res;
		int n2 = full - 1;
		while (m_moves[n2] == -1) { n2--; }

		int start = n1;
		int end = n2;
		int mid = (start + end) / 2;

// std::cout << "MID = " << mid << ", START/END = "
// << start << "/" << end << "\n";

		while (n1 < mid)
		{
			double frac = (double) (n1 - start) / (double) (mid - start);
			int val = (int) (frac * m_moves[n1] + 0.5);

			if (val < 1) { val = 1; }

			int diff = m_moves[n1] - val;

			m_moves[n1] = val;
			m_moves[n2] += diff;

			for (n1++;m_moves[n1] == -1;n1++)
			{ }

			for (n2--;m_moves[n2] == -1;n2--)
			{ }
		}
	}

	// add the left over moves to the later extrusions
// std::cout << "ADDING " << extra << " EXTRA\n";

	for (int k = full - 1;k > 0 && extra > 0;k--)
	{
		if (m_moves[k] != -1)
		{
			++m_moves[k];
			--extra;
		}
	}

/*
	std::cout << "MOVES:\n";
	int total = 0;
	for (int z = 0;z < full;z++)
	{
		if (m_moves[z] != -1)
		{
			std::cout << z << ") " << m_moves[z] << "\n";
			total += m_moves[z];
		}
	}
	std::cout << "TOTAL " << total << " = " << m_growth_moves << "\n";
	exit(0);
*/
}

void Extender_Fixed::print_template(std::ostream &out,
	bool commented /*= true*/)
{
	const char *c = (commented ? "#" : "");

	out << c << "type = " << c_type << "\n"
		<< c_param_move_distrib << " = " << c_default_move_distrib
			<< "\t\t# distribution of max moves between extrusions "
			   "(\"fixed\" or \"linear\")\n"
		<< "\n";
}

