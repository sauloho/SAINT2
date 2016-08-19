#ifndef PARSE_H_INCLUDED
#define PARSE_H_INCLUDED

/// @file Miscellaneous parsing functions.

/// @brief Convert a string to an integer. Exits with an error message
/// on error.
///
/// @param s Source string.
/// @param len Maximum length (ie. only s[0 .. (len - 1)] is parsed).
/// @param filename File being parsed.
/// @param nline Current line number in file.
/// @param desc Description of file (used in error messages).
/// @return Integer value.
int parse_integer(const char *s, int len, const char *filename,
	int nline, const char *desc);

/// @brief Convert a string to a double. Exits with an error message
/// on error.
///
/// @param s Source string.
/// @param len Maximum length (ie. only s[0 .. (len - 1)] is parsed).
/// @param filename File being parsed.
/// @param nline Current line number in file.
/// @param desc Description of file (used in error messages).
/// @return Double value.
double parse_double(const char *s, int len, const char *filename,
	int nline, const char *desc);

#endif // PARSE_H_INCLUDED
