#ifndef COMMON_H_INCLUDED
#define COMMON_H_INCLUDED

/// @file Common types and convenience functions.

#include <string>
#include <vector>

extern bool reverseSaint;

/* types */

/// @brief A vector of doubles.
typedef std::vector<double> Double_Vec;

/* string related functions */

/// @brief Get a lower case copy of a string.
/// @param str Source string.
/// @return Result string.
std::string tolower(const std::string &str);

/// @brief Get a copy of a string with the first character capitalised.
/// @param str Source string.
/// @return Result string.
std::string capitalise(const std::string &str);

/// @brief Convert a boolean value to a string.
/// @param val Value.
/// @return The string "true" or "false".
const char *bool_str(bool val);

/// @brief Check if two values are approximately equal.
bool approx_equal(double x, double y, double tolerance = 0.01);

/// @brief Check if two angles (in radians) are approximately equal.
bool approx_equal_angle(double x, double y);

/* parsing functions */

/// @brief Parse an integer. Exits with an error message on failure.
/// @param Source string.
/// @param name Description (used in error messages).
/// @return The number.
int parse_integer(const std::string &str, const std::string &name);

/// @brief As above, with an optional minimum value.
int parse_integer(const std::string &str,
	const std::string &name, int min_val, bool no_min_val = false);

/// @brief Parse a floating point number.
/// Exits with an error message on failure.
/// @param Source string.
/// @param name Description (used in error messages).
/// @return The number.
double parse_double(const std::string &str, const std::string &name);

/// @brief As above, with an optional minimum value.
double parse_double(const std::string &str,
	const std::string &name, double min_val, bool no_min_val = false);

/// @brief Convert a string to a boolean value (eg. true/false, t/f,
/// yes/no or y/n). Exits with an error message on failure
/// @param Source string.
/// @param name Description (used in error messages).
/// @return The value.
bool parse_bool(const std::string &str, const std::string &name);

/// @brief Get the current date and time as a string
std::string get_date();

/* miscellaneous functions */

/// @brief Replacement for std::min() in <algorithm> (doesn't seem to exist).
template <typename T> T min_val(const T &val1, const T &val2)
{
	return (val1 < val2 ? val1 : val2);
}

/// @brief Replacement for std::max() in <algorithm> (doesn't seem to exist).
template <typename T> T max_val(const T &val1, const T &val2)
{
	return (val1 > val2 ? val1 : val2);
}

template <typename T> void swap(T &val1, T &val2)
{
	T temp = val1;
	val1 = val2;
	val2 = temp;
}

template <typename T> T square(const T &val)
{
	return val * val;
}

#endif // COMMON_H_INCLUDED
