#ifndef CONFIG_H_INCLUDED
#define CONFIG_H_INCLUDED

#include <iostream>
#include <string>
#include "param_list.h"

/// @breief Program version number.
#define SAINT2_VERSION "1.1"

// forward declarations
class Runner;
class Sequence;

extern bool reverseSaint;

/// @brief Class for managing program configuration information
/// (from a configuration file, the command line or a mixture of both).
///
/// The general format of a configuration file is:
/// <pre> 
///   [section]
///   name = value
///   ...
///   # comment ...
/// </pre>
/// Long lines may be continued on the next line by ending them with "+".
///
/// Only one Config object can exist (this is enforced using the
/// m_instantiated variable).
///
/// Example:
///
/// <pre>
/// int main(int argc, const char **argv)
/// {
///     Config(argc, argv);
///     Runner runner(config);
///     ...
/// }
/// </pre>

class Config
{
public:
	/// @brief Constructor.
	///
	/// @param argc Command line argument count (argument to main()).
	/// @param argv Command line arguments (argument to main()).
	Config(
		int argc,				// in - number of command line arguments
		const char *argv[]		// in - command line arguments
	);

	/// @brief Destructor.
	~Config();

	/// @brief Get the number of runs to perform ("-n" command line option).
	/// @return Number of runs.
	int num_runs() const
	{ return m_num_runs; }

	/// @brief Get the random number seed to use.
	/// @return The seed (if not specified on the command line with "-s",
	/// the seed is based on the current time (note: resolution is seconds).
	long random_seed() const
	{ return m_rnd_seed; }

	/// @brief Get the name of the output file.
	/// @return The filename: if no "-o" option was specified on the
	/// command line, returns an empty string; otherwise, if just one
	/// run is being performed, returns the "-o" argument; otherwise
	/// returns the "-o" argument with the run number appended.
	const std::string outfilename(int run_number) const;

	/// @brief Initialise a Runner object with the appropriate type of
	/// scorer, strategy, etc.
	///
	/// @param r [out] Runner object.
	void init_runner(Runner *r) const;

	/// @brief Initialise an amino acid (or codon) sequence
	/// with the configuration values (ie. either from a FASTA file
	/// or a literal sequence in the configuration file).
	///
	/// @param s [out] Amino acid sequence.
	void create_sequence(Sequence *s) const;

	/// @brief Print parameter values in configuration file format.
	void print_params() const;

	/// @brief Get the PDB filename specified on the command line ("--" option).
	/// @return Filename.
	const std::string &pdb_filename() const
	{ return m_pdb_filename; }

	/// @brief Check if the "-t" flag was included on the command line
	/// @return True if flag was included.
	bool show_torsion_angles() const
	{ return m_show_torsion_angles; }

	/// @brief Check if the "-b" flag was included on the command line
	/// @return True if flag was included.
	bool backbone_only() const
	{ return m_backbone_only; }

	/// @brief Check if the "-v" flag was included on the command line
	/// @return True if flag was included.
	bool verbose() const
	{ return m_verbose; }

	/// @brief Get the chain id specified on the command line (following
	/// the PDB filename).
	/// @return Chain id (' ' if not specified on the command line).
	char pdb_chain() const
	{ return m_pdb_chain; }

	/// @brief Get the executable name (argv[0] from the command line).
	/// @return Executable name.
	static const std::string &cmd()
	{ return m_cmd; }

	/// @brief Print an error message about a missing configuration file
	/// parameter and exit.
	static void missing_parameter(
		const char *section,		// config file section
		const char *type,			// "type" value (may be empty / NULL)
		const char *param_name		// parameter name
	);

private:
	/// @brief Sections of the configuration file
	/// (each section has a corresponding class).
	enum Section
	{
		S_Sequence, S_Runner, S_Extend, S_Scoring, S_Strategy, S_Movement,
		Num_Sections,
		S_None
	};

private:
	// disallow copy and assignment by declaring private versions.
	Config(const Config&);
	Config& operator = (const Config&);

	// @brief Print usage message
	void show_usage();

	/// @brief Parse command line arguments.
	///
	/// @param argc Command line argument count.
	/// @param argv Command line arguments.
	void parse(int argc, const char *argv[]);

	/// @brief Parse a single command line option.
	///
	/// @param arg [in,out] Index of current argument in \a argv.
	/// @param argc Command line argument count. 
	/// @param argv Command line arguments.
	void parse_option(int *arg, int argc, const char *argv[]);

	/// @brief Check if the end of the argument list has been reached
	/// prematurely (after the specified command line option).
	///
	/// @param argc Command line argument count. 
	/// @param argv Command line arguments.
	/// @param option Current option (without "-", eg. "o").
	void check_end_of_args(int arg, int argc, const std::string &option);

	/// @brief Check if a command line argument is a (unique) abbreviation
	/// of a configuration file section name; if so, parse the
	/// parameters that follow.
	///
	/// @param arg [in,out] Index of current argument in \a argv.
	/// @param argc Command line argument count. 
	/// @param argv Command line arguments.
	bool parse_section_args(int *arg, int argc, const char *argv[]);

	/// @brief Parse a list of configuration paramaters from the command line
	/// (ie. arguments containing "=").
	///
	/// @param sec Current section.
	/// @param arg [in,out] Index of current argument in \a argv.
	/// @param argc Command line argument count. 
	/// @param argv Command line arguments.
	void parse_param_args(Section sec, int *arg, int argc,
	    const char *argv[]);

	/// @brief Read a configuration file (the filename was specified
	/// on the command line, and is stored in m_config_filename).
	void read_config_file();

	/// @brief Convert a section name to the corresponding enum
	/// (\a str should point to the character after the initial '[').
	/// Exits with a message on error.
	///
	/// @param str Source string.
	/// @return The section.
	Section parse_section_name(const char *str);

	/// @brief Parse a parameter inside a config file
	/// (leading white space is assumed to have already been skipped,
	/// and the line is assumed not to be a comment).
	/// Exits with a message on error.
	///
	/// @param sec Current section in config file.
	/// @param name [in,out] (Last) parameter name; if non-empty, the last
	/// line may be continued with the character "+".
	/// @param str Current line.
	void parse_config_param(Section sec, std::string &name,
		const char *str);

	/// @brief Check if there were any incompatible parameter values
	/// in different sections of the configuration file.
	void verify_interdependent_params();

	/// @brief Transform all escape sequences in a string (such as \")
	/// to single characters.
	///
	/// @param str [in,out] String to transform.
	void transform_escape_sequences(std::string &str);

	/// @brief Print an error message about a syntax error on the current
	/// line of the config file and exit.
	void file_syntax_error();

	/// @brief Output a template configuration file.
	void print_config_template();

private:
	/// Line continuation character (used at the start of continued lines).
	static const char c_continue_char;

	/// Default number of runs to perform.
	static const int c_default_runs;

	/// Parameters for each section of the configuration file.
	Param_List m_param[Num_Sections];

	/// Name of each section.
	std::string m_section_name[Num_Sections];

	/// (Used for making sure only one object of this class ever exists).
	static bool m_instantiated;

	/// Command name (argv[0] from command line).
	static std::string m_cmd;

	/// Name of configuration file (if any).
	std::string m_config_filename;

	/// Current line number in config file.
	int m_config_line_num;

	/// Number of runs to perform.
	int m_num_runs;

	/// Output file (if any).
	std::string m_outfile;

	/// Random number seed.
	long m_rnd_seed;

	/// PDB file name (after "--" on command line).
	std::string m_pdb_filename;

	/// "-t" flag.
	bool m_show_torsion_angles;

	/// "-b" flag.
	bool m_backbone_only;

	/// "-v" flag.
	bool m_verbose;

	/// Chain id following PDF filename (if any).
	char m_pdb_chain;
};

#endif // CONFIG_H_INCLUDED

