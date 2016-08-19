#ifndef C_FILE_H_INCLUDED
#define C_FILE_H_INCLUDED

#include <cstdio>
#include <string>

/// @brief Check if a file exists and is readable.
bool file_exists(const char *filename);

/// @brief Wrapper for C FILE* type.
///
/// Can be more convenient to use than a std::iostream.
/// Closes the file automatically on destruction.

class C_File
{
public:
	/// @brief Constructor.
	///
	/// Opens a file in a particular mode (as used by fopen(), eg. "r", "w").
	/// If the filename is empty, either stdin or stdout is used,
	/// depending on the mode.
	/// An object of this class is implicitly converted to type FILE*,
	/// so it can be used as an argument to functions such as fprintf()
	/// and getc().
	/// Exits with an error message if the file cannot be opened.
	///
	/// Example:
	/// <pre>
	/// C_File outfile("fred", "w", "results file");
	/// fprintf(outfile, "Hello\n");
	/// </pre>
	/// @param filename Name of file; empty string means stdin or stdout
	/// (depending on whether the mode starts with "r" or "w"/"a")
	/// @param mode File mode (as in fopen())
	/// @param desc Brief description of the file (used in error messages)
	C_File(const char *filename, const char *mode,
		const char *desc);

	/// @brief Constructor (same as above using a std::string).
	C_File(const std::string &filename, const char *mode,
		const char *desc);

	/// @brief Destructor. Closes the file if it is not stdin or stdout.
	~C_File();

	/// @brief Get the name of the file.
	/// @return The filename (or the string "stdin" or "stdout" if the
	/// original filename was empty).
	const std::string &name() const
	{ return m_filename; }

	/// @brief (Implicit) conversion to a FILE*.
	/// @return The file pointer.
	operator FILE* ();

	/// @brief Get the next line, ignoring blank lines and comments ("# ...")
	/// @return False at end of file.
	bool next_line(char *line, int max_len);

	/// @breief Get the current line number (if next_line() used)
	/// @return The number of the line just read, or 0 if none have been read.
	int line_num() const
	{ return m_line_num; }

private:
	// disallow copy and assignment by making them private.
	C_File(const C_File &);
	C_File &operator = (const C_File &);

	/// Open the file in the mode specified (called by constructors).
	void init(const std::string &filename, const char *mode,
		const char *desc);

private:
	/// @brief Name of file (set to "stdin" or "stdout" if original file name
	/// was empty).
	std::string m_filename;

	/// File pointer.
	FILE *m_file;

	/// Set to true if the file is stdin or stdout.
	bool m_std_file;

	/// Number of lines read by next_line()
	int m_line_num;
};


#endif // C_FILE_H_INCLUDED

