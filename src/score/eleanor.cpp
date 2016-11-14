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
	double scr;

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
		m_con[0][span]=i;
		m_con[1][span]=j;
	}
	
	m_num_spans=span;
	num_con=span;
	std::cout << "number of spans: " << span << "\n";
	
	for (i=0; i<p.full_length(); i++)
	{
		m_layer[i] = 0;
		for (span=0;span<m_num_spans;span++)
		{
			if (i>=m_con[0][span] && i<=m_con[1][span])
			{
				m_layer[i] = 1;
			}
		}
		std::cout << m_layer[i];
	}
	std::cout << "\n";


    m_satisfied_con = (int *)malloc(sizeof(int)*(num_con+1));
    m_previous_len = 0;
    // p.alloc_satisfied_con (num_con);
	m_data_loaded = true;
	fclose(input_file);	
}

double Eleanor::score(const Peptide& p, bool verbose)
{
	load_data(p);

	int len = p.length();
	int i,j,k,cont=0;
    int possible_con=0;

	double total=0.0;
	Point cb_i, cb_j;

    if(len != m_previous_len)
    {
        if((p.length()-1) % (p.full_length()/10) == 0)
        {
		    std::ofstream myfile;
	        char myString[150];
		    std::strcpy(myString,"scmatrix_perc");
		    std::strcat(myString,p.get_filename());
			myfile.open(myString,std::ios_base::app);
          	  myfile<<len<<" ";
           	 for(k=0;k<num_con+1;k++)
               		 myfile<<  m_satisfied_con[k]<<" ";
           	 myfile<<"\n";
           	 m_previous_len=len;
			myfile.close();
        }

        if((p.length()-1) % 25 == 0 )
        {
		    std::ofstream myfile;
		    char myString[150];
		    std::strcpy(myString,"scmatrix_part");
		    std::strcat(myString,p.get_filename());
            if((p.length()-1) % (p.full_length()/10) == 0) /* This means that p.get_filename() retrieved the perc file instead of the part file. */
            {
                int c;
                for(c=std::strlen(myString);myString[c]!='p';c--);
                myString[c+1]='a';
                myString[c+3]='t';
		    }
            myfile.open(myString,std::ios_base::app);
            myfile<<len<<" ";
            for(k=0;k<num_con+1;k++)
                myfile<< m_satisfied_con[k] << " ";
            myfile<<"\n";
		    myfile.close();
        }
        m_previous_len=len;
    }

    for(k=0;k<num_con;k++) m_satisfied_con[k]=-1;

	for (k=0;k<num_con;k++)
	{
		i = m_con[0][k]; 	j=m_con[1][k];

		if( i-1 >= p.start() &&  i-1 <= p.end() && j-1 <= p.end() && j-1 >= p.start() )
		{
			possible_con++;
			if(std::strcmp(p.res(i-1).amino().abbr(),"GLY") == 0)
			{
				cb_i = p.atom_pos(i-1, Atom_CA);				
			}
			else
				cb_i = p.atom_pos(i-1, Atom_CB);				

			if(!std::strcmp(p.res(j-1).amino().abbr(),"GLY"))
			{
				cb_j = p.atom_pos(j-1, Atom_CA);				
			}
			else
				cb_j = p.atom_pos(j-1, Atom_CB);				

			if ( cb_i.distance(cb_j) > 8.0) /* They are predicted to be contacts, but are far away in the model! */
			{
				total += cb_i.distance(cb_j) - 8.0;
				cont++;
                //std::cout << "0 ";
                m_satisfied_con[k]=0;
			}
            else
                //std::cout << "1 ";
                m_satisfied_con[k]=1;
		}
	}
	
    if(possible_con)
	    m_satisfied_con[num_con]=100*(possible_con - cont)/possible_con;
    else
        m_satisfied_con[num_con]=-1;
#ifndef RAW_SCORE

	// Alter the score so that it has about the same distribution
	// for all lengths
	total /= sqrt((double) len);
	// normalise so that all score types have approximately the same range
	total = total * 28.0;
	
#endif // RAW_SCORE
//	std::cout << "Total = " << total << " !!\n";	
	return total;
}

