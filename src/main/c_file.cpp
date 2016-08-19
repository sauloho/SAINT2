
#include <cstdlib>
#include <cassert>
#include "c_file.h"
#include "config.h"

C_File::C_File(const char *filename, const char *mode, const char *desc)
{
	init(filename, mode, desc);
}

C_File::C_File(const std::string &filename, const char *mode,
	const char *desc)
{
	init(filename, mode, desc);
}

void C_File::init(const std::string &filename, const char *mode,
	const char *desc)
{
	if (filename.empty())
	{
		m_std_file = true;

		if (*mode == 'r')
		{
			m_filename = "stdin";
			m_file = stdin;
		}
		else
		if (*mode == 'w' || *mode == 'a')
		{
			m_filename = "stdout";
			m_file = stdout;
		}
		else
		{
			std::cerr << Config::cmd()
				<< ": unknown mode '" << mode
				<< "' in C_File::init()\n";
			exit(1);
		}
	}
	else
	{
		m_std_file = false;
		m_filename = filename;
		m_file = fopen(m_filename.c_str(), mode);

		if (m_file == NULL)
		{
			std::cerr << Config::cmd()
				<< ": cannot open " << desc
				<< " \"" << filename
				<< "\"\n";
			exit(1);
		}
	}

	m_line_num = 0;
}

C_File::~C_File()
{
	if (m_file != NULL && !m_std_file)
	{
		fclose(m_file);
	}
}

C_File::operator FILE* ()
{
	assert(m_file != NULL);
	return m_file;
}

bool C_File::next_line(char *line, int max_len)
{
	assert(m_file != NULL);

	for ( ; ; )
	{
		if (fgets(line, max_len, m_file) == NULL)
		{
			return false;
		}

		m_line_num++;

		char *p;
		for (p = line;*p && isspace(*p);p++)
		{
		}

		if (*p && *p != '#')
		{
			return true;
		}
	}
}

bool file_exists(const char *filename)
{
	// fstat() isn't portable, so just try opening the file

	FILE *f = fopen(filename, "r");

	if (f != NULL)
	{
		fclose(f);
		return true;
	}

	return false;
}


