
#include <iostream>
#include <cassert>
#include "common.h"
#include "strategy_strict.h"

// static data members
const char *Strategy_Strict::c_type = "strict";
const char *Strategy_Strict::c_param_candidates = "number";
const int Strategy_Strict::c_default_candidates = 1;

Strategy_Strict::Strategy_Strict()
	: m_candidates(c_default_candidates)
{
}

Strategy_Strict::~Strategy_Strict()
{
}

bool Strategy_Strict::parse_parameter(const std::string &name,
    const std::string &value)
{
	std::string full_name = Strategy::config_section();
	full_name += " ";
	full_name += c_type;

	if (name == c_param_candidates)
	{
		set_num_candidates(parse_integer(value, full_name, 1));
		return true;
	}

	return false;
}

void Strategy_Strict::verify_parameters()
{
}

void Strategy_Strict::set_num_candidates(int num)
{
	assert(num >= 1);
	m_candidates = 1;
}

int Strategy_Strict::select(double old_score, const Double_Vec &new_score)
{
	double best_score = old_score;
	int best = -1;

	for (unsigned int n = 0;n < new_score.size();n++)
	{
		if (new_score[n] < best_score)
		{
			best_score = new_score[n];
			best = n;
		}
	}

	return best;
}

void Strategy_Strict::start_run(Runner *runner)
{
}

void Strategy_Strict::end_run(Runner *runner)
{
}

bool Strategy_Strict::stop()
{
	// stopping is handled by the Runner
	return false;
}

void Strategy_Strict::print_template(std::ostream &out,
	bool commented /*= true*/)
{
	const char *c = (commented ? "#" : "");

    out << c << "type = " << c_type << "\n"
		<< c << c_param_candidates << " = " << c_default_candidates << "\n"
		<< "\n";
}

