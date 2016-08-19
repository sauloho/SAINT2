// -*- mode: c++; c-file-style: "ellemtel"; encoding: utf-8; -*-
//
// Copyright (C) 2012 Saulo H. P. de Oliveira
//
// Saulo H. P. de Oliveira <saulo.deoliveira@pmb.ox.ac.uk>
// 2012-08-06
//
// Saulo's scoring method.
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
#include "predtor.h"

PredTor::PredTor()
{
	m_data_loaded=false;
	m_maxtor=0;
}

PredTor::~PredTor()
{
}

void PredTor::set_short_data_file(const std::string &filename)
{
	m_filename = filename;
}

void PredTor::set_long_data_file(const std::string &filename)
{
	m_filename  = filename;
}

void PredTor::load_data(const Peptide& p)
{
	FILE *input_file; 
	int i,j;
	char c,Res[5],SS[5];

	if (m_data_loaded)
	{
		return;
	}

	if (m_filename.empty())
	{
		std::cerr << "Error:  Predicted Torsion Angle File not defined!\n";
		exit(1);
	}

	input_file = fopen(m_filename.c_str(),"r");

	if (input_file == NULL)
	{
		std::cerr << "Predicted Torsion Angle File not found !\n";
		exit(1);
	}

        m_phipsi[0] = (double *) malloc(sizeof(double) * (p.full_length()+1));
        m_phipsi[1] = (double *) malloc(sizeof(double) * (p.full_length()+1));


        /* Remove Header from SPINE-X output file */
        for(c=fgetc(input_file); c!='\n'; c=fgetc(input_file) );

        for(i=0; fscanf(input_file,"%d %s %s %lf %lf",&j,Res,SS,&m_phipsi[0][i], &m_phipsi[1][i])!=EOF && i < p.full_length(); i++)
	{
	        for(c=fgetc(input_file);c!='\n' ; c=fgetc(input_file));
//		std::cerr << "Predicted Torsion Angles: " << m_phipsi[0][i] << " " << m_phipsi[1][i] << "\n";
	}

	m_data_loaded = true;
	fclose(input_file);	
}

double PredTor::score(const Peptide& p, bool verbose)
{
	load_data(p);
	int len = p.length();
	int i;
	double phi2,psi2,total=0.0;

        for (i = p.start() ; i <= p.end()-2 ; i++)
        {
                phi2 = p.conf().phi(i+1);
                psi2 = p.conf().psi(i+1);

                phi2 = fmod(phi2*57.29578,360);
                psi2 = fmod(psi2*57.29578,360);
		total += abs(phi2 - m_phipsi[0][i+1]) + abs(psi2 - m_phipsi[1][i+1]);	
        }



#ifndef RAW_SCORE
	// Alter the score so that it has about the same distribution
	// for all lengths
	total /= (double) len;
	// normalise so that all score types have approximately the same range
//	total = total*28;
	
#endif // RAW_SCORE
	//std::cerr << "Total = " << total << " !!\n";	
	return total;
}

