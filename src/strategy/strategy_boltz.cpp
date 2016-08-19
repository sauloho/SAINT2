
#include <cassert>
#include <cmath>
#include "config.h"
#include "random.h"
#include "stream_printf.h"
#include "strategy_boltz.h"

// static data members
const char *Strategy_Boltz::c_type = "boltz";
const char *Strategy_Boltz::c_param_temp = "temperature";
const char *Strategy_Boltz::c_param_num = "number";
const double Strategy_Boltz::c_default_temp = 1.0;
const int Strategy_Boltz::c_default_num = 10;

Strategy_Boltz::Strategy_Boltz()
	: m_temp(c_default_temp)
{
	set_number(c_default_num);
}

Strategy_Boltz::~Strategy_Boltz()
{
}

bool Strategy_Boltz::parse_parameter(const std::string &name,
    const std::string &value)
{
	std::string full_name = Strategy::config_section();
	full_name += " ";
	full_name += c_type;
	full_name += " ";
	full_name += name;

	if (name == c_param_temp)
	{
		set_temperature(parse_double(value, full_name, 0.01));
		return true;
	}
	else
	if (name == c_param_num)
	{
		set_number(parse_integer(value, full_name, 1));
		return true;
	}

	return false;
}

void Strategy_Boltz::verify_parameters()
{
}

void Strategy_Boltz::set_temperature(double val)
{
	m_temp = val;
}

void Strategy_Boltz::set_number(int num)
{
	assert(num > 0);
	m_num = num;
	m_prob.resize(m_num);
}

int Strategy_Boltz::select(double old_score, const Double_Vec &new_score)
{
	assert(new_score.size() == (unsigned) m_num);
	assert(m_prob.size() == (unsigned) m_num);

	double prob_total = 0.0;
	int n;

//std::cout << "\n";

//double highest = 0.0;
//int highest_n = -1;
//double highest_change = 0;

	for (n = 0;n < m_num;n++)
	{
		double change = new_score[n] - old_score;
		m_prob[n] = exp(-change / m_temp);

//std::cout << n << ") old " << old_score << " new " << new_score[n]
//	<< " diff " << change << " prob " << m_prob[n] << "\n";

		// these numbers can be extremely small, so to find the total
		// accurately they should be added from smallest to largest -
		// but for these purposes it should work OK to be approximate

		prob_total += m_prob[n];

//if (m_prob[n] > highest) { highest = m_prob[n]; highest_n = n; highest_change = change; }
	}

//std::cout << "STRATEGY: highest prob = " << highest
//<< " (change " << highest_change << ")"
//<< " = #" << highest_n << ", prob total = " << prob_total << "\n";
//double orig_prob_total = prob_total;

//std::cout << "\nProb total " << prob_total << "\n";

	// check if the current structure is selected
	// (probability value = exp(0) = 1.0)

	if (Random::rnd(prob_total + 1.0) < 1.0)
	{
//std::cout << "SELECTED: -1\n";
		// -1 means do nothing
		return -1;
	}

	for (n = m_num - 1;n >= 1;n--)
	{
		if (Random::rnd(prob_total) < m_prob[n])
		{
//std::cout << "SELECTED: " << n << " with prob "
//<< m_prob[n] << "\n";

//std::cout << "Selected " << n << "\n";
			return n;
		}

		prob_total -= m_prob[n];
	}
//std::cout << "Selected last (0)\n";

	// remaining probability = 100%
	return 0;
}

void Strategy_Boltz::start_run(Runner *runner)
{
}

void Strategy_Boltz::end_run(Runner *runner)
{
}

bool Strategy_Boltz::stop()
{
	return false;
}

void Strategy_Boltz::print_template(std::ostream &out,
	bool commented /*= true*/)
{
	const char *c = (commented ? "#" : "");

    out << c << "type = " << c_type << "\n"
        << c << c_param_temp << " = "
			<< Printf("%.1f", c_default_temp) << "\n"
		<< c << c_param_num << " = "
			<< c_default_num << "\n"
		<< "\n";
}

