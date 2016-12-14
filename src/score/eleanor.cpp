// -*- mode: c++; c-file-style: "ellemtel"; encoding: utf-8; -*-
//
// Copyright (C) 2012 Saulo H. P. de Oliveira
//
// Saulo H. P. de Oliveira <saulo.deoliveira@pmb.ox.ac.uk>
// 2012-08-06
//
// Predicted Contact Scoring Method.
//
// 


#include <string>
#include <cstdlib>
#include <cstdio>
#include <cstring>
#include <cmath>
#include <iostream>
#include <fstream>
#include "amino.h"
#include "residue.h"
#include "peptide.h"
#include "atom.h"
#include "scorer_combined.h"
#include "eleanor.h"

Eleanor::Eleanor()
{
	m_data_loaded=false;
}

Eleanor::~Eleanor()
{
}

void Eleanor::set_short_data_file(const std::string &filename)
{
	m_filename = filename;
}

void Eleanor::set_long_data_file(const std::string &filename)
{
	m_filename  = filename;
}

void Eleanor::load_data(const Peptide& p)
{
	FILE *input_file; 
	int i,j,span,tot_resis;
	int doubled_midpoints[5000];

	if (m_data_loaded)
	{
		return;
	}

	if (m_filename.empty())
	{
		std::cerr << "Error:  Predicted Contacts File not defined!\n";
		exit(1);
	}

	input_file = fopen(m_filename.c_str(),"r");

	if (input_file == NULL)
	{
		std::cerr << "Predicted Contacts File not found !\n";
		exit(1);
	}

	m_con[0] = (int *)malloc(sizeof(int)*5000);
	m_con[1] = (int *)malloc(sizeof(int)*5000);

	//Rosetta-generated spanfile from SpanningTopology object
	//4 208 <- this is a mystery, should be 132 for 1orsC I think 
	//antiparallel
	//n2c
	//	7	27  
	//	35	56  
	//	81	88  
	//	98	122 
	
	fscanf(input_file,"Rosetta-generated spanfile from SpanningTopology object\n");
	fscanf(input_file,"%d %d\n",&m_num_spans,&tot_resis);
	std::cout << "number of spans: " << m_num_spans << "\n";
	// check that tot_resis == peptide full length
	fscanf(input_file,"antiparallel\n");
	fscanf(input_file,"n2c\n");

	for (span=0;fscanf(input_file,"\t%d\t%d\n",&i,&j) != EOF ; span++)
	{
		m_con[0][span] = i;
		m_con[1][span] = j;
		doubled_midpoints[span] = i + j;
	}
	
	m_num_spans=span;
	num_con=span;
	std::cout << "number of spans: " << span << "\n";
	
	// initialise array values to zeros
	for (i=0; i<p.full_length(); i++)
	{
		m_layer[i] = 0;
	}

	// go through and allocate 1 or 2 to central residues in each helix
	//  midpoint-¬
	//           v 
	// odd:  1,2,2,2,1
	// even:  2,2,2,2
	for (i=0; i<span; i++)
	{
		if (doubled_midpoints[i] % 2 == 1)
		{
			m_layer[(doubled_midpoints[i]-1)/2 - 1] = 2;
			std::cout << "Efirst " << (doubled_midpoints[i]-1)/2 - 1 << "\n";
			m_layer[(doubled_midpoints[i]-1)/2 + 0] = 2;
			m_layer[(doubled_midpoints[i]-1)/2 + 1] = 2;
			m_layer[(doubled_midpoints[i]-1)/2 + 2] = 2;
			std::cout << "Elast  " << (doubled_midpoints[i]-1)/2 + 2 << "\n";
		}
		else
		{
			m_layer[doubled_midpoints[i]/2 - 2] = 1;
			std::cout << "Ofirst " << doubled_midpoints[i]/2 - 2 << "\n";
			m_layer[doubled_midpoints[i]/2 - 1] = 2;
			m_layer[doubled_midpoints[i]/2 + 0] = 2;
			m_layer[doubled_midpoints[i]/2 + 1] = 2;
			m_layer[doubled_midpoints[i]/2 + 2] = 1;
			std::cout << "Olast  " << doubled_midpoints[i]/2 + 2 << "\n";
		}
	}
	//for (i=0; i<p.full_length(); i++)
	//{
	//	m_layer[i] = 0;
	//	for (span=0;span<m_num_spans;span++)
	//	{
	//		if (i>=m_con[0][span] && i<=m_con[1][span])
	//		{
	//			m_layer[i] = 1;
	//		}
	//	}
	//	std::cout << m_layer[i];
	//}
	//std::cout << "\n";


    m_previous_len = 0;
    // p.alloc_satisfied_con (num_con);
	m_data_loaded = true;
	fclose(input_file);	
}

double Eleanor::score(const Peptide& p, bool verbose)
{
	load_data(p);

	int len = p.length();
	int k = 0;

	double total=0.0;
	Point cb_i, cb_j;
	float tolerance = 10;

	if(len != m_previous_len)
	{
		m_previous_len=len;
	}

	for (k=0; k<len; k++)
	{
		if (fabs(p.atom_pos(k, Atom_CA).z) > tolerance) // not within 5A of middle of membrane
		{
			total += (fabs(p.atom_pos(k, Atom_CA).z) - tolerance) * (double) m_layer[k];
			//std::cout << k << " CA z-coord: " << p.atom_pos(k, Atom_CA).z << "  m_layer[" << k << "]: " <<  m_layer[k] <<  "\n";
			//std::cout << "total: " << total << "\n";
		}
	}

	
#ifndef RAW_SCORE

	// Alter the score so that it has about the same distribution
	// for all lengths
	//total /= sqrt((double) len);
	// normalise so that all score types have approximately the same range
	total = total * 28.0;
	
#endif // RAW_SCORE
//	std::cout << "Total = " << total << " !!\n";	
	return total;
}

