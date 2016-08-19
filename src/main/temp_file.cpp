#include <cstdlib>
#include <cstdio>
#include <unistd.h>
#include <list>
#include <iostream>
#include "temp_file.h"

#include <algorithm> // PG added this
#include <vector> // PG added this

// static member variables
int Temp_File::m_seq = 0;

// local functions

// (The purpose of these functions is to delete any temporary
// files created by class Temp_File if the program terminates by calling
// exit(), which unfortunately prevents any destructors from being called).

namespace
{
	typedef std::list<std::string> Filename_List;
	Filename_List files_to_delete;

	void delete_temp_files_on_exit()
	{
		Filename_List::iterator i = files_to_delete.begin();

		for ( ;i != files_to_delete.end();++i)
		{
			std::cout << "At exit: removing " << *i << "\n";
			remove(i->c_str());
		}
	}

	void add_to_delete_on_exit_list(const std::string &f)
	{
		// check if it is already in the list

		Filename_List::iterator i = std::find(
			files_to_delete.begin(), files_to_delete.end(), f);

		if (i == files_to_delete.end())
		{
			files_to_delete.push_back(f);
		}
	}

	void remove_from_delete_on_exit_list(const std::string &f)
	{
		Filename_List::iterator i = std::find(
			files_to_delete.begin(), files_to_delete.end(), f);

		if (i != files_to_delete.end())
		{
			files_to_delete.erase(i);
		}
	}

	// This class exists so that we can have an object whose
	// constructor is called at program initialisation to register
	// an "atexit" function.

	class Register_atexit
	{
	public:
		Register_atexit()
		{
			atexit(delete_temp_files_on_exit);
		}
	};

	Register_atexit reg_atexit;
}

///////////////////////////////////////////////////////////////////

Temp_File::Temp_File(const char *prefix /*= "temp" */)
{
	char buffer[1000];

	if (prefix == NULL)
	{
		prefix = "";
	}

	const char *job_id = getenv("LSB_JOBID");

	if (job_id == NULL)
	{
		job_id = "";
	}

	sprintf(buffer, "%s%ld_%s_%d.tmp",
		prefix, (long) getpid(), job_id, m_seq++);
	m_name = buffer;
	add_to_delete_on_exit_list(m_name);
}

Temp_File::~Temp_File()
{
	remove(m_name.c_str());
	remove_from_delete_on_exit_list(m_name);
}

const char *Temp_File::name() const
{
	return m_name.c_str();
}

