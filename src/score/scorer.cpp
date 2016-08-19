
#include <cstdlib>
#include <iostream>
#include "scorer.h"
#include "scorer_combined.h"
#include "common.h"
#include "config.h"

// macro to call another macro for each sublass of Scorer

#define for_each_Scorer_subclass(MACRO_NAME) \
do { \
	MACRO_NAME(Scorer_Combined); \
} while (0)

// static data members

const char *Scorer::c_config_section = "scoring";
const char *Scorer::c_param_type = "type";
const char *Scorer::c_default_type = Scorer_Combined::type();

Scorer::Scorer() :
	m_score_info_on(false)
{
}

Scorer::~Scorer()
{
}

void Scorer::print_info_when_scoring(bool info_on)
{
	m_score_info_on = info_on;
}

void Scorer::set_verbose(bool val)
{
	m_verbose = val;
}

Scorer *Scorer::create(const Param_List &params)
{
	std::string type = find(params, c_param_type, c_default_type);
	Scorer *s = NULL;
#define CREATE_OBJECT(c) if (type == c::type()) { s = new c; break; }
    for_each_Scorer_subclass(CREATE_OBJECT);
#undef CREATE_OBJECT
	if (s == NULL)
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
		if (i->name == c_param_type)
		{
			// already handled
			continue;
		}
		else
		if (!s->parse_parameter(i->name, i->value))
		{
			std::cerr << Config::cmd()
				<< ": unknown " << c_config_section << " parameter \""
				<< i->name << "\" for type " << type << "\n";
			exit(1);
		}
	}
	s->verify_parameters();
	return s;
}

void Scorer::print_template(std::ostream &out)
{
	out << "[" << capitalise(c_config_section) << "]\n\n";

    std::string def = c_default_type;

#define PRINT_TEMPLATE(c) c::print_template(out, (def != c::type()))
	for_each_Scorer_subclass(PRINT_TEMPLATE);
#undef PRINT_TEMPLATE
}

