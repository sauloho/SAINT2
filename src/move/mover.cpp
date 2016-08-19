
#include <cstdlib>
#include <iostream>
#include "config.h"
#include "common.h"
#include "mover.h"
#include "mover_fragment.h"
#include "mover_fragment_fwd.h"
#include "mover_fragment_rev.h"

// macro to call another macro for each sublass of Mover

#define for_each_Mover_subclass(MACRO_NAME) \
do { \
      if (reverseSaint) { MACRO_NAME(Mover_Fragment_Rev); } \
	  else { MACRO_NAME(Mover_Fragment_Fwd); } \
} while (0)

// static data members

const char *Mover::c_config_section = "movement";
const char *Mover::c_param_type = "type";
const char *Mover::c_default_type =	Mover_Fragment::type();

Mover::Mover()
{
}

Mover::~Mover()
{
}

Mover *Mover::create(const Param_List &params)
{
	std::string type = find(params, c_param_type, c_default_type);
	Mover *m = NULL;

#define CREATE_OBJECT(c) if (type == c::type()) { m = new c; break; }
    for_each_Mover_subclass(CREATE_OBJECT);
#undef CREATE_OBJECT

	if (m == NULL)
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
		if (!m->parse_parameter(i->name, i->value))
		{
			std::cerr << Config::cmd()
				<< ": unknown " << c_config_section << " parameter \""
				<< i->name << "\" for type " << type << "\n";
			exit(1);
		}
	}

	m->verify_parameters();
	return m;
}

void Mover::print_template(std::ostream &out)
{
	out << "[" << capitalise(c_config_section) << "]\n\n";

	std::string def = c_default_type;

#define PRINT_TEMPLATE(c) c::print_template(out, (def != c::type()))
    for_each_Mover_subclass(PRINT_TEMPLATE);
#undef PRINT_TEMPLATE
}

