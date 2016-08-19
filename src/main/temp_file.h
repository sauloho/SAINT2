#ifndef TEMP_FILE_H_INCLUDED
#define TEMP_FILE_H_INCLUDED

#include <string>

/// @brief Temporary filename.
///
/// A temporary filename, automatically deleted on destruction
/// or on program exit.

class Temp_File
{
public:
	/// Constructor.
	
	/// A prefix of "abc" will generate a name like "abc123_1.tmp",
	/// where 123 is the process id and 1 is a sequence number
	/// (incremented each time a Temp_File is created).
	Temp_File(const char *prefix = "temp");

	/// Destructor (deletes the file).
	~Temp_File();

	/// Get the name of the file.
	const char *name() const;

private:
	/// Sequence number (used in filename, incremented each time a
	/// Temp_File object is created).
	static int m_seq;

	/// The filename.
	std::string m_name;
};

#endif // TEMP_FILE_H_INCLUDED
