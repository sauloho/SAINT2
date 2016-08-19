#include <cstring> //PG added this
#include <cstdio>
#include <cctype>
#include <cstdlib>
#include "config.h"
#include "stream_printf.h"

Printf::Printf(const char *fmt, double val)
{
	verify_format(fmt, "eEfFgG");
	char buffer[1000];
	snprintf(buffer, 1000, fmt, val);
	m_str = buffer;
}

Printf::Printf(const char *fmt, short val)
{
	verify_format(fmt, "duxX", "h");
	char buffer[1000];
	snprintf(buffer, 1000, fmt, val);
	m_str = buffer;
}

Printf::Printf(const char *fmt, int val)
{
	verify_format(fmt, "duxX");
	char buffer[1000];
	snprintf(buffer, 1000, fmt, val);
	m_str = buffer;
}

Printf::Printf(const char *fmt, long val)
{
	verify_format(fmt, "duxX", "l");
	char buffer[1000];
	snprintf(buffer, 1000, fmt, val);
	m_str = buffer;
}

Printf::Printf(const char *fmt, char val)
{
	verify_format(fmt, "c");
	char buffer[1000];
	snprintf(buffer, 1000, fmt, val);
	m_str = buffer;
}

Printf::Printf(const char *fmt, const char *val)
{
	verify_format(fmt, "s");
	char buffer[1000];
	snprintf(buffer, 1000, fmt, val);
	m_str = buffer;
}

void Printf::verify_format(const char *fmt, const char *type,
	const char *type_prefix /* = NULL*/)
{
	const char *orig_fmt = fmt;
	bool found = false;

	for ( ;*fmt;fmt++)
	{
		if (*fmt == '\\')
		{
			// (doesn't check for \x.. etc.)
			if (*(++fmt) == '\0')
			{
				break;
			}
		}
		else
		if (*fmt == '%')
		{
			if (*(++fmt) == '\0')
			{
				break;
			}

			if (*fmt == '%')
			{
				continue;
			}

			if (*fmt == '+' || *fmt == '-' || *fmt == '#' || *fmt == ' ' ||
				*fmt == '0' || *fmt == '\'' || *fmt == 'I')
			{
				if (*(++fmt) == '\0')
				{
					break;
				}
			}
			else
			if (*fmt == '*')
			{
				star_illegal(orig_fmt);
			}

			for ( ;isdigit(*fmt);fmt++)
			{ }

			if (*fmt == '.')
			{
				if (*(++fmt) == '\0')
				{
					break;
				}

				if (*fmt == '*')
				{
					star_illegal(orig_fmt);
				}

				for ( ;isdigit(*fmt);fmt++)
				{ }
			}

			if (found)
			{
				too_many_error(orig_fmt);
			}
			else
			{
				if (type_prefix != NULL && *type_prefix)
				{
					if (strchr(type_prefix, *fmt) != NULL)
					{
						if (*(++fmt) == '\0')
						{
							break;
						}
					}
					else
					{
						wrong_type_error(orig_fmt, type, type_prefix);
					}
				}

				if (strchr(type, *fmt) != NULL)
				{
					found = true;
				}
				else
				{
					wrong_type_error(orig_fmt, type, type_prefix);
				}
			}
		}
	}

	if (!found)
	{
		not_found_error(orig_fmt);
	}
}

void Printf::star_illegal(const char *fmt)
{
	std::cerr << Config::cmd()
		<< ": '*' in Printf format illegal: \"" << fmt << "\"\n";
	exit(1);
}

void Printf::too_many_error(const char *fmt)
{
	std::cerr << Config::cmd()
		<< ": too many values in Printf format: \"" << fmt << "\"\n";
	exit(1);
}

void Printf::wrong_type_error(const char *fmt, const char *type,
	const char *type_prefix)
{
	std::cerr << Config::cmd()
		<< ": wrong type in Printf format (expected ";

	if (type_prefix != NULL)
	{
		std::cerr << "prefix \"" << type_prefix << "\" then ";
	}

	if (strlen(type) > 1)
	{
		std::cerr << "one of ";
	}

	std::cerr << '\"' << type << "\"): \"" << fmt << "\"\n";
	exit(1);
}

void Printf::not_found_error(const char *fmt)
{
	std::cerr << Config::cmd()
		<< ": value missing in Printf format: \"" << fmt << "\"\n";
	exit(1);
}

std::ostream &operator << (std::ostream &out, const Printf &p)
{
	return out << p.m_str;
}

