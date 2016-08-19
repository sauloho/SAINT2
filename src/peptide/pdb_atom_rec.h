#ifndef PDB_ATOM_REC_H_INCLUDED
#define PDB_ATOM_REC_H_INCLUDED

#include <string>
#include "point.h"
#include "amino.h"
#include "atom_type.h"

// Data from a PDB "ATOM" record, eg.
// ATOM   1666  CG1 ILE A 222      25.647  35.896  17.376  1.00 55.51

struct PDB_Atom_Rec
{
	std::string record;			// usually "ATOM"
	int			serial;			// unique id for this atom
	std::string name;			// four letter atom type, eg. " CG1"
	char		alt_loc;		// (usually ' ')
	std::string res_name;		// three letter amino acid code
	char		chain_id;		// chain (may be ' ')
	int			res_seq;		// unique residue id (along with i_code)
	char		i_code;			// (usually ' ')
	Point		pos;			// position
	double		occupancy;		// occupancy fraction (usually 1.0)
	double		temp_factor;	// temperature factor
	std::string seg_id;			// (optional)
	std::string element;		// (optional)
	std::string charge;			// (optional)

	// parse an "ATOM" line
	// Returns false on error (with *err_msg containing the error; may be NULL)
	bool parse(const char *line, std::string *err_msg);

	// set the main fields (the other fields are set to reasonable defaults).
	// Call this function before using operator << to print a record
	void set_values(int serial_num, Atom_Type a_type, int res_seq_num,
		char icode, Amino amino, const Point &position, char chain = 'A');

	// print a PDB "ATOM" line
	// (the buffer should be at least 80 characters long)
	void format(char *buffer) const;
};

// print a PDB "ATOM" line
std::ostream& operator << (std::ostream &out, const PDB_Atom_Rec &a);

#endif // PDB_ATOM_REC_H_INCLUDED
