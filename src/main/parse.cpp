
#include <cstdlib>
#include <cstdio>
#include <string>
#include <iostream>
#include "parse.h"

// extract an integer value from a string of length len
// (exits with a message on error)

int parse_integer(const char *s, int len, const char *filename,
	int nline, const char *desc)
{
	int val;
	std::string p(s, len);

	if (sscanf(p.c_str(), "%d", &val) != 1)
	{
		std::cerr << "Error parsing " << desc
			<< " on line " << nline
			<< " of " << filename
			<< "\n";
		exit(1);
	}

	return val;
}

double parse_double(const char *s, int len, const char *filename,
	int nline, const char *desc)
{
	double val;
	std::string p(s, len);

	if (sscanf(p.c_str(), "%lf", &val) != 1)
	{
		std::cerr << "Error parsing " << desc
			<< " on line " << nline
			<< " of " << filename
			<< "\n";
		exit(1);
	}

	return val;
}

