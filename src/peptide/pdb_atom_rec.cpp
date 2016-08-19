#include <stdio.h>
#include <cstring>
#include <cctype>
#include <iostream>
#include "pdb_atom_rec.h"

bool PDB_Atom_Rec::parse(const char *line, std::string *err_msg)
{
// std::cout << "LINE: '" << line << "'\n";
	int line_len = strlen(line);

	if (line_len != 0)
	{
		// skip spaces at end
		for (line_len--;line_len >= 0 && isspace(line[line_len]);line_len--)
		{ }

		line_len++;
	}

	// if (line_len < 66)
	if (line_len < 47)
	{
		if (err_msg != NULL) { *err_msg = "line too short"; }
		return false;
	}

	record = std::string(line, 6);

	if (sscanf(std::string(line + 6, 5).c_str(), "%d", &serial) != 1)
	{
		if (err_msg != NULL) { *err_msg = "illegal serial value"; }
		return false;
	}

	name = std::string(line + 12, 4);
	alt_loc = line[16];
	res_name = std::string(line + 17, 3);
	chain_id = line[21];

	if (sscanf(std::string(line + 22, 4).c_str(), "%d", &res_seq) != 1)
	{
		if (err_msg != NULL) { *err_msg = "illegal residue sequence value"; }
		return false;
	}

	i_code = line[26];

	if (sscanf(std::string(line + 30, 8).c_str(), "%lf", &pos.x) != 1 ||
		sscanf(std::string(line + 38, 8).c_str(), "%lf", &pos.y) != 1 ||
		sscanf(std::string(line + 46, 8).c_str(), "%lf", &pos.z) != 1)
	{
		if (err_msg != NULL) { *err_msg = "illegal coordinate"; }
		return false;
	}

	if (line_len >= 55)
	{
		if (sscanf(std::string(line + 54, 6).c_str(), "%lf", &occupancy) != 1)
		{
			if (err_msg != NULL) { *err_msg = "illegal occupancy value"; }
			return false;
		}
	}
	else
	{
		occupancy = 1.0;
	}

	if (line_len >= 61)
	{
		if (sscanf(std::string(line + 60, 6).c_str(), "%lf", &temp_factor) != 1)
		{
			// if (err_msg != NULL) { *err_msg = "illegal temperature factor"; }
			// return false;
	
			temp_factor = 50.0;
		}
	}
	else
	{
		temp_factor = 50.0;
	}

	// seg_id, element and charge are all optional

	if (line_len >= 73)
	{
		seg_id = std::string(line + 72, 4);
	}
	else
	{
		seg_id.clear();
	}

	if (line_len >= 77)
	{
		element = std::string(line + 76, 2);
	}
	else
	{
		element.clear();
	}

	if (line_len >= 79)
	{
		charge = std::string(line + 78, 2);
	}
	else
	{
		charge.clear();
	}

	return true;
}

void PDB_Atom_Rec::set_values(int serial_num, Atom_Type a_type, int res_seq_num,
	char icode, Amino amino, const Point &position, char chain /*= 'A'*/)
{
	record = "ATOM";
	serial = serial_num;
	name = a_type.pdb_name();
	alt_loc = ' ';
	res_name = amino.abbr();
	chain_id = chain;
	res_seq = res_seq_num;
	i_code = icode;
	pos = position;
	occupancy = 1.0;
	temp_factor = 20.0;
	seg_id.clear();
	element.clear();
	charge.clear();
}

void PDB_Atom_Rec::format(char *buffer) const
{
	// format:
	//
	// RRRRRRSSSSS NNNNARRR CSSSSI   XXXXXXXXYYYYYYYYZZZZZZZZOOOOOOTTTTTT
	// ATOM   1666  CG1 ILE A 222      25.647  35.896  17.376  1.00 55.51
	//
	// + optional seg_id, element, charge

	sprintf(buffer,
		"%-6.6s"	// record (eg. "ATOM")
		"%5d"		// serial
		" "
		"%4.4s"		// atom name (eg. " CG1")
		"%c"		// alt loc
		"%3.3s"		// residue name (amino acid)
		" "
		"%c"		// chain
		"%4d"		// residue sequence number
		"%c"		// i code
		"   "
		"%8.3f"		// x coordinate
		"%8.3f"		// y coordinate
		"%8.3f"		// z coordinate
		"%6.2f"		// occupancy
		"%6.2f" 	// temperature factor
		"     "
		"%4.4s"		// seg_id (optional)
		"%2.2s" 	// element (optional)
		"%2.2s",	// charge (optional)
		this->record.c_str(),
		this->serial,
		this->name.c_str(),
		this->alt_loc,
		this->res_name.c_str(),
		this->chain_id,
		this->res_seq,
		this->i_code,
		this->pos.x,
		this->pos.y,
		this->pos.z,
		this->occupancy,
		this->temp_factor,
		this->seg_id.c_str(),
		this->element.c_str(),
		this->charge.c_str());
}

std::ostream& operator << (std::ostream &out, const PDB_Atom_Rec &a)
{
	char buffer[100];

	a.format(buffer);
	return out << buffer << '\n';
}

