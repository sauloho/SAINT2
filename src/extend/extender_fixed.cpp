#include <cstdlib> // PG added this
#include <iostream>
#include <math.h>
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

	if (m_move_distrib == Distrib_Fixed)
	{
		for (int i = m_initial_res;i < full;i += m_extrude_res)
		{
			m_moves[i] = num;
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
	}
	
	if (m_move_distrib == Distrib_Linear)
	{
		// sum of all lengths at which moves are performed
		// up to and including this length
		std::vector<int> cumulative_lengths;
		cumulative_lengths.resize(full);

		// rounded number of moves which should have been
		// performed before start of next extrusion
		std::vector<int> cumulative_moves;
		cumulative_moves.resize(full);

		std::cout << "total moves: " << m_growth_moves << "\n";
		std::cout << "i/Cumulative lengths/Cumulative moves/moves:\n";
		for (int i = 0; i < m_initial_res; i++)
		{
			cumulative_lengths[i] = 0;
			cumulative_moves[i] = 0;
			std::cout << i << "\t";
			std::cout << cumulative_lengths[i] << "\t";
			std::cout << cumulative_moves[i] << "\t";
			std::cout << m_moves[i] << "\n";
		}

		for (int i = m_initial_res; i < full; i++)
		{
			cumulative_lengths[i] = cumulative_lengths[i-1] + i;
		}


		for (int i = m_initial_res; i < full; i++)
		{
			cumulative_moves[i] = round( m_growth_moves * cumulative_lengths[i] / (float) cumulative_lengths[full - 1] );
			m_moves[i] = cumulative_moves[i] - cumulative_moves[i - 1];
			std::cout << i << "\t";
			std::cout << cumulative_lengths[i] << "\t";
			std::cout << cumulative_moves[i] << "\t";
			std::cout << m_moves[i] << "\n";
		}

	}

	// print out to check move distibution
	std::cout << "Linear move distribution:\n";
	for (int i = 0;i < full;i++)
	{
		std::cout << i << ": " << m_moves[i] << "\n";
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

