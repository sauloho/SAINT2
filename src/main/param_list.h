#ifndef PARAM_LIST_H_INCLUDED
#define PARAM_LIST_H_INCLUDED

/// @file Functions to do with configuration parameters (name/value pairs).
///
/// This is a set of functions rather than a class because
/// std::list cannot be subclassed; it could all be wrapped inside
/// a class, but it doesn't seem worth the effort.

#include <string>
#include <list>

/// @brief A name/value pair.

struct Name_Value
{
	std::string name;
	std::string value;

	Name_Value()
	{ }

	Name_Value(const std::string &n, const std::string &v)
		: name(n), value(v)
	{ }
};

/// @brief (A std::list is used instead of a std::map so that the
/// order is maintained; later values have priority over earlier ones)
typedef std::list<Name_Value> Param_List;

/// @brief Insert a name/value into a parameter list
/// (replaces the existing value, if any)
void insert(Param_List &nv, const std::string &name,
	const std::string &value, const std::string &section_name);

/// @brief Insert a name/value into parameter list (appends to existing
/// value, if any)
void append_to_value(Param_List &nv, const std::string &name,
	const std::string &value);

/// @brief Retrieve a value from a parameter list.
/// @return Pointer to the value string (NULL if name is not found).
std::string *find(Param_List &nv, const std::string &name);

/// @brief As above (const version).
const std::string *find(const Param_List &nv, const std::string &name);

/// @brief Retrieve a value from a parameter list.
/// @return The value string, or \a default_val if not found.
std::string find(const Param_List &nv, const std::string &name,
	const std::string &default_val);

#endif // PARAM_LIST_H_INCLUDED

