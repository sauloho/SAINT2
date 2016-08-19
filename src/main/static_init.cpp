
#include "static_init.h"
#include "amino.h"
#include "codon.h"
#include "atom_type.h"

static Static_Init the_static_initialiser;

Static_Init::Static_Init()
{
	Amino::do_static_init();
	Codon::do_static_init();
	Atom_Type::do_static_init();
}

Static_Init::~Static_Init()
{
	Codon::do_static_delete();
	Amino::do_static_delete();
	Atom_Type::do_static_delete();
}

