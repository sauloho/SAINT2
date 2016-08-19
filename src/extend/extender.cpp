#include <cstdlib>
#include "extender.h"
#include "extender_fixed.h"
#include "extender_codon.h"
#include "common.h"
#include "config.h"
#include "runner.h"
#include "peptide.h"

// macro to call another macro for each sublass of Extender

#define for_each_Extender_subclass(MACRO_NAME) \
do { \
	MACRO_NAME(Extender_Fixed); \
	MACRO_NAME(Extender_Codon); \
} while (0)

// static data members

const char *Extender::c_config_section = "extension";

const char *Extender::c_param_type = "type";
const char *Extender::c_param_initial = "initial";
const char *Extender::c_param_extrude = "extrude";
const char *Extender::c_param_growth_moves = "growth_moves";

const char *Extender::c_default_type = Extender_Fixed::type();
const int Extender::c_default_initial = 9;
const int Extender::c_default_extrude = 1;
const long Extender::c_default_growth_moves = 10000;

Extender::Extender()
    : m_initial_res(c_default_initial),
	  m_extrude_res(c_default_extrude),
	  m_growth_moves(c_default_growth_moves)
{
}

Extender::~Extender()
{
}

void Extender::set_initial_res(int val)
{
	m_initial_res = val;
}

void Extender::set_extrude_res(int val)
{
	m_extrude_res = val;
}

void Extender::set_growth_moves(long val)
{
	m_growth_moves = val;
}

int Extender::curr_length_move_limit(const Peptide &p)
{
	assert(p.length() < (int) m_moves.size());
	return m_moves[p.length()];
}

int Extender::check_extend(Runner *runner)
{
	const Peptide &p = runner->peptide();

	if (p.full_grown())
	{
		return 0;
	}

	assert(p.length() < (int) m_moves.size());
	assert(m_moves[p.length()] != -1);

	if (this->must_extend(runner->peptide(), runner) ||
		runner->curr_length_moves() >= m_moves[p.length()])
	{
		if (p.length() + m_extrude_res > p.full_length())
		{
			// extrude the rest of the residues
			return p.full_length() - p.length();
		}

		return m_extrude_res;
	}

	return 0;
}

void Extender::start_run(const Sequence &seq)
{
	calculate_num_moves(seq);
}

void Extender::after_extend(const Peptide &p, Runner *runner)
{
}

bool Extender::must_extend(const Peptide &p, Runner *runner)
{
	return false;
}

Extender *Extender::create(const Param_List &params)
{
	std::string type = find(params, c_param_type, c_default_type);
	Extender *e = NULL;

#define CREATE_OBJECT(c) if (type == c::type()) { e = new c; break; }
	for_each_Extender_subclass(CREATE_OBJECT);
#undef CREATE_OBJECT

	if (e == NULL)
	{
        std::cerr << Config::cmd()
            << ": unknown " << c_config_section << " type \""
            << type << "\"\n";
        exit(1);
	}

	// parse the parameters

    Param_List::const_iterator i = params.begin();
    Param_List::const_iterator end = params.end();

    for ( ;i != end;++i)
    {
        std::string full_name = c_config_section;
        full_name += " ";
        full_name += i->name;

		if (i->name == c_param_type)
		{
			// already handled
			continue;
		}
		else
        if (i->name == c_param_initial)
        {
            e->set_initial_res(parse_integer(i->value, full_name, 1));
        }
        else
        if (i->name == c_param_extrude)
        {
            e->set_extrude_res(parse_integer(i->value, full_name, 1));
        }
        else
        if (i->name == c_param_growth_moves)
        {
            e->set_growth_moves(parse_integer(i->value, full_name, 0));
        }
		/*	
		else
        if (i->name == c_param_move_distrib)
		{
			if (i->value == "fixed")
			{
				e->m_move_distrib = Distrib_Fixed;
			}
			else if (i->value == "linear")
			{
				e->m_move_distrib = Distrib_Linear;
			}
			else
			{
				std::cerr << Config::cmd()
					<< ": illegal value \"" << i->value
					<< "\" for " << c_config_section << " parameter \""
	                << i->name << "\", type " << type << "\n";
				exit(1);
			}
		}
		*/
        else
		if (!e->parse_parameter(i->name, i->value))
        {
            std::cerr << Config::cmd()
                << ": unknown " << c_config_section << " parameter \""
                << i->name << "\" for type " << type << "\n";
            exit(1);
        }
	}

	e->verify_parameters();
	return e;
}

void Extender::print_template(std::ostream &out)
{
	out << "[" << capitalise(c_config_section) << "]\n\n";

	std::string def = c_default_type;

	out << c_param_initial << " = " << c_default_initial
			<< "\t\t\t# initial number of residues extruded\n"
		<< c_param_extrude << " = " << c_default_extrude
			<< "\t\t\t# number of residues to extrude at a time\n"
		<< c_param_growth_moves << " = " << c_default_growth_moves
			<< "\t\t# total number of moves during growth\n"
	<< "\n";

#define PRINT_TEMPLATE(c) c::print_template(out, (def != c::type()))
	for_each_Extender_subclass(PRINT_TEMPLATE);
#undef PRINT_TEMPLATE

}

