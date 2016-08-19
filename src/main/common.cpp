
#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <string>
#include <cctype>
#include <ctime>
#include <cmath>
#include "config.h"
#include "common.h"

std::string tolower(const std::string &str)
{
	std::string s = str;

	for (unsigned int n = 0;n < s.length();n++)
	{
		s[n] = tolower(s[n]);
	}

	return s;
}

std::string capitalise(const std::string &str)
{
	std::string s = str;
	s[0] = toupper(s[0]);
	return s;
}

bool approx_equal(double x, double y, double tolerance /*= 0.01*/)
{
	return (fabs(x - y) < tolerance);
}

bool approx_equal_angle(double x, double y)
{
	static const double pi2 =  M_PI * 2.0;

	x = fmod(x, pi2);
	y = fmod(y, pi2);

	if (x < 0.0) { x += pi2; }
	if (y < 0.0) { y += pi2; }

	if (fabs(x - pi2) < 0.1 || fabs(y - pi2) < 0.1)
	{
		x = fmod(x + M_PI, pi2);
		y = fmod(y + M_PI, pi2);
	}

	return (fabs(x - y) < M_PI / 180.0);	// pi/180 is 1 degree
}

const char *bool_str(bool val)
{
	return (val ? "true" : "false");
}

int parse_integer(const std::string &str,
	const std::string &name, int min_val, bool no_min_val /*= false*/)
{
	int val, pos;
	if (sscanf(str.c_str(), "%d%n", &val, &pos) != 1 ||
		pos != (int) str.length())
	{
		std::cerr << Config::cmd() << ": " << name << " must be an integer\n";
		exit(1);
	}

	if (!no_min_val && val < min_val)
	{
		std::cerr << Config::cmd() << ": " << name
			<< " must be at least " << min_val << "\n";
		exit(1);
	}

	return val;
}

int parse_integer(const std::string &str, const std::string &name)
{
	return parse_integer(str, name, 0, true);
}

double parse_double(const std::string &str,
	const std::string &name, double min_val, bool no_min_val /*= false*/)
{
	int pos;
	double val;
	if (sscanf(str.c_str(), "%lf%n", &val, &pos) != 1 ||
		pos != (int) str.length())
	{
		std::cerr << Config::cmd() << ": " << name << " must be a number\n";
		exit(1);
	}

	if (!no_min_val && val < min_val)
	{
		std::cerr << Config::cmd() << ": " << name
			<< " must be at least " << min_val << "\n";
		exit(1);
	}

	return val;
}

double parse_double(const std::string &str, const std::string &name)
{
	return parse_double(str, name, 0.0, true);
}

bool parse_bool(const std::string &str, const std::string &name)
{
	std::string s = tolower(str);

	if (s == "true" || s == "t" || s == "yes" || s == "y")
	{
		return true;
	}
	if (s == "false" || s == "f" || s == "no" || s == "n")
	{
		return false;
	}

	std::cerr << Config::cmd() << ": " << name
		<< " can only be true or false\n";
	exit(1);
}

std::string get_date()
{
    char date_str[100];
    time_t t = time(0);
    strftime(date_str, 100, "%a %b %d %H:%M:%S", localtime(&t));
    return std::string(date_str);
}

