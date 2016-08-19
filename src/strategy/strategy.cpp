
#include <cstdlib>
#include <iostream>
#include "config.h"
#include "stream_printf.h"
#include "strategy.h"
#include "common.h"
#include "strategy_strict.h"
#include "strategy_monte.h"
#include "strategy_boltz.h"
#include "strategy_always.h"

#define for_each_Strategy_subclass(MACRO_NAME) \
do { \
	MACRO_NAME(Strategy_Strict); \
	MACRO_NAME(Strategy_Monte); \
	MACRO_NAME(Strategy_Boltz); \
	MACRO_NAME(Strategy_Always); \
} while (0)

// static data members

const char *Strategy::c_config_section = "strategy";
const char *Strategy::c_param_type = "type";
const char *Strategy::c_default_type = Strategy_Monte::type();

Strategy::Strategy()
{
}

Strategy::~Strategy()
{
}

Strategy *Strategy::create(const Param_List &params)
{
	std::string type = find(params, c_param_type, c_default_type);
	Strategy *s = NULL;

#define CREATE_OBJECT(c) if (type == c::type()) { s = new c; break; }
	for_each_Strategy_subclass(CREATE_OBJECT);
#undef CREATE_OBJECT

	if (s == NULL)
	{
		std::cerr << Config::cmd()
			<< ": unknown "
			<< c_config_section
			<< " type \""
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
				<< i->name << "\"\n";
			exit(1);
		}
	}

	s->verify_parameters();
    return s;
}

void Strategy::print_template(std::ostream &out)
{
	out << "[" << capitalise(c_config_section) << "]\n\n";
	
	std::string def = c_default_type;

#define PRINT_TEMPLATE(c) c::print_template(out, (def != c::type()))
    for_each_Strategy_subclass(PRINT_TEMPLATE);
#undef PRINT_TEMPLATE
}

