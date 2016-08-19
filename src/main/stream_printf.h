#ifndef STREAM_PRINTF_H_INCLUDED
#define STREAM_PRINTF_H_INCLUDED

/// Class for conveniently using printf() with a stream.
///
/// The format string is checked before printing (some rare formats
/// may not be handled).
///
/// Example:
/// <pre>
/// std::cout << Printf("%.3f", val) << "\n";
/// </pre>

#include <string>
#include <iostream>

class Printf
{
	friend std::ostream &operator << (std::ostream &out, const Printf &p);

public:
	// Constructors for various types.
	Printf(const char *fmt, double val);
	Printf(const char *fmt, short val);
	Printf(const char *fmt, int val);
	Printf(const char *fmt, long val);
	Printf(const char *fmt, char val);
	Printf(const char *fmt, const char *val);

private:
	/// Make sure the format string is correct. Exits with a message on error.
	/// @param fmt Format string (eg. "%d").
	/// @param type Expected format char(s) (eg. "duxX" for int).
	/// @param type_prefix Permitted modifiers (eg. "l", "h").
	void verify_format(const char *fmt, const char *type,
		const char *type_prefix = NULL);

	/// Called if "*" is found in the format string (exits with message).
	void star_illegal(const char *fmt);

	/// Called if there are too many "%"s in a format (exits with message).
	void too_many_error(const char *fmt);

	/// Called if the wrong "%" type is found in a format (exits with message).
	void wrong_type_error(const char *fmt, const char *type,
		const char *type_prefix);

	/// Called if no "%" is found in a format (exits with message).
	void not_found_error(const char *fmt);

private:
	/// String to print.
	std::string m_str;
};

std::ostream &operator << (std::ostream &out, const Printf &p);

#endif // STREAM_PRINTF_H_INCLUDED
