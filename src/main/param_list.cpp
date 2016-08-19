
#include "config.h"
#include "param_list.h"

void insert(Param_List &nv, const std::string &name,
	const std::string &value, const std::string &section_name)
{
	Param_List::iterator i = nv.begin();
	Param_List::iterator end = nv.end();

	for ( ;i != end;++i)
	{
		if (i->name == name)
		{
			if (i->value == value)
			{
				return;
			}

			std::cerr
				<< Config::cmd()
				<< ": warning: multiple values for "
				<< section_name << " parameter \""
				<< name << "\"; using last\n";

			nv.erase(i);
			break;
		}
	}

	nv.push_back(Name_Value(name, value));
}

void append_to_value(Param_List &nv, const std::string &name,
	const std::string &value)
{
	std::string *v = find(nv, name);

	if (v != NULL)
	{
		*v += value;
	}
	else
	{
		nv.push_back(Name_Value(name, value));
	}
}

std::string *find(Param_List &nv, const std::string &name)
{
	Param_List::iterator i = nv.begin();
	Param_List::iterator end = nv.end();

	for ( ;i != end;++i)
	{
		if (i->name == name)
		{
			return &(i->value);
		}
	}

	return NULL;
}

const std::string *find(const Param_List &nv, const std::string &name)
{
	return find(const_cast<Param_List&>(nv), name);
}

std::string find(const Param_List &nv, const std::string &name,
	const std::string &default_val)
{
	const std::string *v = find(nv, name);
	return (v == NULL ? default_val : *v);
}

