
#include <cstdlib>
#include <stdio.h>
#include "runner.h"
#include "run_observer.h"
#include "config.h"
#include "scorer.h"
#include "scorer_combined.h"
#include "strategy.h"
#include "extender.h"
#include "peptide.h"
#include "conformation.h"
#include "mover.h"

const char *Runner::m_config_section =		"general";
//int template_count = 0;

const char *Runner::c_param_sequential =	"sequential";
//const char *Runner::c_param_coil =			"coil";
const char *Runner::c_param_move_limit =	"moves";
const char *Runner::c_param_no_sel_limit =	"none selected";
//const char *Runner::c_param_scwrl_executable = "scwrl_executable";
const char *Runner::c_param_reverse =		"reverse";
const char *Runner::c_param_start_struct =	"start_structure";
const char *Runner::c_param_native_struct =	"native_structure";

// default parameter values
const bool Runner::c_default_sequential		= true;
//const bool Runner::c_default_coil			= false;
const long Runner::c_default_move_limit		= 10000;
const long Runner::c_default_no_sel_limit	= 0;

Runner::Runner(Config &config) :
	m_num_runs(0),
	m_scorer(NULL), m_strategy(NULL), m_extender(NULL), m_mover(NULL),
	m_run(0), m_curr_score(0.0), m_prev_score(0.0),
	m_move_failed(false),
	m_curr_length_moves(0), m_no_sel_count(0),
	m_sequential(c_default_sequential),
	//m_coil(c_default_coil),
	m_move_limit(c_default_move_limit),
	m_no_sel_limit(c_default_no_sel_limit),
	m_native_known(false)
{
	// create the scorer, strategy, etc.
	config.init_runner(this);
}

Runner::~Runner()
{
}

void Runner::set_scorer(Scorer *s)
{
	delete m_scorer;
	m_scorer = s;
}

void Runner::set_strategy(Strategy *s)
{
	delete m_strategy;
	m_strategy = s;
}

void Runner::set_extender(Extender *e)
{
	delete m_extender;
	m_extender = e;
}

void Runner::set_mover(Mover *m)
{
	delete m_mover;
	m_mover = m;
}

void Runner::set_num_runs(int num_runs)
{
	m_num_runs = num_runs;
}

void Runner::set_output(std::string outfile)
{
	m_outfile = outfile;
}

void Runner::parse_params(const Param_List &params)
{
	Param_List::const_iterator i = params.begin();
	Param_List::const_iterator end = params.end();

	for ( ;i != end;++i)
	{
		std::string full_name = m_config_section;
		full_name += " ";
		full_name += i->name;

        if (i->name == c_param_sequential)
        {
            m_sequential = parse_bool(i->value, full_name);
        }
		else
		/*
        if (i->name == c_param_coil)
        {
            m_coil = parse_bool(i->value, full_name);
        }
		else
		*/
		if (i->name == c_param_move_limit)
		{
			m_move_limit = parse_integer(i->value, full_name, 0);
		}
		else
		if (i->name == c_param_no_sel_limit)
		{
			m_no_sel_limit = parse_integer(i->value, full_name, 0);
		}
		else
		/*
		if (i->name == c_param_scwrl_executable)
		{
			Peptide::set_scwrl_exec(i->value);
		}
		else
		*/
        if (i->name == c_param_reverse)
        {
			// (global variable)
            reverseSaint = parse_bool(i->value, full_name);
        }
		else
		if (i->name == c_param_start_struct)
		{
			m_start_struct = i->value;
		}
		else
		if (i->name == c_param_native_struct)
		{
			m_native_struct = i->value;
		}
		else
		{
			std::cerr << Config::cmd()
				<< ": unknown " << m_config_section << " parameter \""
				<< i->name << "\"\n";
			exit(1);
		}
	}

	//if (!m_start_struct.empty() && m_sequential)
	//{
	//	std::cerr <<
	//		"Warning: using start structure for sequential folding - may not work!\n";
	//	//m_start_struct = "";
	//}

	// no other parameter verification needed
}

void Runner::do_runs(Sequence &seq, Run_Observer &observer)
{
	// check if the strategy has a specialised do_runs() function
	// (eg. strategies that require an ensemble of peptides
	// instead of just one at a time)

	if (m_strategy->do_runs(*this, seq, observer))
	{
		return;
	}

	// set up the candidate vectors

	const int num_candidates = m_strategy->num_candidates();
	Conf_Vec candidate;
	Double_Vec candidate_score, candidate_progress1_score;
	candidate.resize(num_candidates);
	candidate_score.resize(num_candidates);
	candidate_progress1_score.resize(num_candidates);

	observer.before_start(this);

	if (!m_native_struct.empty())
	{
		if (!m_native_peptide.read_pdb(m_native_struct.c_str()))
		{
			std::cerr << "Error while reading native structure "
				<< m_native_struct << "\n";
			exit(1);
		}

		m_native_peptide.conf().calc_torsion_angles();
		m_native_known = true;
	}

	for (m_run = 0;m_run < m_num_runs;m_run++)
	{
		std::cout << "Run #" << m_run << "\n";
		m_peptide.create_from_sequence(seq);
		std::cout << "created peptide from sequence\n";
		std::cout << "peptide full_length: " << m_peptide.full_length() << "\n";
// if (m_run % 2 == 1)
// { std::cout << "Switching to in vitro\n"; m_sequential = false; }
// else { std::cout << "Switching to cotrans\n"; m_sequential = true; }

		if (m_sequential)
		{
			if (!m_start_struct.empty())
			{
				if (!m_peptide.read_segment_from_pdb(m_start_struct.c_str(), m_extender->initial_residues() + 1))
				{
					std::cerr << "Error while reading start structure "
						<< m_start_struct << "\n";
					exit(1);
				}

				std::cout << "read start structure\n";
				std::cout << "peptide full_length: " << m_peptide.full_length() << "\n";
				std::cout << "peptide length: " << m_peptide.length() << "\n";
				m_peptide.remove_non_backbone_atoms();
				std::cout << "removed non-backbone atoms\n";
				m_peptide.conf().calc_torsion_angles();
				std::cout << "calculated torsion angles\n";
				m_peptide.idealise_bond_lengths();

				char pdb_out[150];
				sprintf(pdb_out,"%s_part%dtest1",m_outfile.c_str(),m_peptide.length());
				m_peptide.write_pdb(pdb_out);

				std::cout << "starting init sequential from pdb segment\n";
				m_mover->init_sequential_from_segment(m_peptide,
						m_extender->initial_residues() + 1, &observer);

				sprintf(pdb_out,"%s_part%dtest2",m_outfile.c_str(),m_peptide.length());
				m_peptide.write_pdb(pdb_out);

			}
			else
			{
				std::cout << "starting init sequential from fragment\n";
				m_mover->init_sequential(m_peptide,
						m_extender->initial_residues(), &observer);
			}
		}
		else
		{
			if (!m_start_struct.empty())
			{
				if (!m_peptide.read_pdb(m_start_struct.c_str()))
				{
					std::cerr << "Error while reading start structure "
						<< m_start_struct << "\n";
					exit(1);
				}

				m_peptide.remove_non_backbone_atoms();
				m_peptide.conf().calc_torsion_angles();
				m_peptide.idealise_bond_lengths();
				m_mover->init_from_peptide(m_peptide, &observer);
			}
			else
			{
				//m_mover->init_non_sequential(m_peptide, m_coil, &observer);
				m_mover->init_non_sequential(m_peptide, false, &observer);
			}
		}

		m_curr_score = m_scorer->score(m_peptide, 0.0);
		m_prev_score = m_curr_score;
		m_move_failed = false;

		double best_score = 9e99;	// best score found (after full grown)
		Conformation best_conf;		// best scoring conformation

		m_curr_length_moves = 0;
		m_no_sel_count = 0;

		m_strategy->start_run(this);
		m_extender->start_run(seq);
		observer.start_run(this);

		if (m_sequential)
		{
			observer.after_extend(this, m_peptide.length());
		}

		for ( ; ; )
		{
			// check if it is time to extend
			while (!m_peptide.full_grown())
			{
				int num_res = m_extender->check_extend(this);

				if (num_res == 0)
				{
					break;
				}

                                /* EDIT: Saulo on March, 18th */
                                /* Print decoys once 25,50,75,100... residues have been extruded */
                                if(m_peptide.length() % 1 == 0)
                                {
                                        char pdb_out[150];
                                        sprintf(pdb_out,"%s_part%d",m_outfile.c_str(),m_peptide.length());
                                        m_peptide.write_pdb(pdb_out);
                                }
                                /* Print decoys once ~10%,20%,30%,... of the residues have been extruded */
              //                  if(m_peptide.length() % (m_peptide.full_length()/10) == 0)
                //                {
                  //                      char pdb_out[150];
                    //                    sprintf(pdb_out,"%s_perc%d",m_outfile.c_str(),m_peptide.length());
                      //                  m_peptide.write_pdb(pdb_out);
                        //        }



				/* (incorrect -- progress based weight is wrong)
				// use the best scoring structure found at this length
				if (best_score < 9e99)
				{
					m_peptide.conf().swap(best_conf);
					m_curr_score = best_score;
					best_score = 9e99;
				}
				*/

				// TO DO: should make sure this cast is possible first
				bool ribosome_wall =
					((Scorer_Combined*) m_scorer)->ribosome_wall();

				m_mover->extend(m_peptide, num_res, ribosome_wall, &observer);
				m_extender->after_extend(m_peptide, this);

				m_prev_score = m_curr_score;
				m_curr_score = m_scorer->score(m_peptide, 1.0);
				m_move_failed = false;
				m_no_sel_count = 0;

				observer.after_extend(this, num_res);
	//			char pdb_out[150];
	//			sprintf(pdb_out,"%s_part%d",m_outfile.c_str(),m_peptide.length());
	//			m_peptide.write_pdb(pdb_out);

				m_curr_length_moves = 0;
				m_no_sel_count = 0;

				best_score = m_curr_score;
				best_conf = m_peptide.conf();
			}

			// check if it is time to stop
			if (m_peptide.full_grown())
			{
				if (m_curr_length_moves >= m_move_limit ||
					(m_no_sel_limit > 0 && m_no_sel_count >= m_no_sel_limit) ||
					m_strategy->stop())
				{
					break;
				}
			}

			bool exhaustive_for_pos = false;

			m_mover->do_random_move(m_peptide, num_candidates,
				exhaustive_for_pos, candidate, &observer);


			double progress = m_curr_length_moves /
				(double) (m_peptide.full_grown() ? m_move_limit :
					m_extender->curr_length_move_limit(m_peptide));

			for (int m = 0;m < num_candidates;m++)
			{
				m_peptide.conf().swap(candidate[m]);
				candidate_score[m] = m_scorer->score(m_peptide, progress,
					&candidate_progress1_score[m]);
				m_peptide.conf().swap(candidate[m]);
			}

			m_prev_score = m_curr_score;

			// select one of the candidates
			int choice = m_strategy->select(m_curr_score, candidate_score);
			bool is_best = false;

			if (choice != -1)
			{
				m_peptide.conf().swap(candidate[choice]);
				m_curr_score = candidate_score[choice];
				m_move_failed = false;
				m_no_sel_count = 0;
				//sprintf(pdb_out,"Test%5d.pdb",template_count);
				//if(template_count % 2 == 0)
				//{
				//	m_peptide.write_pdb(pdb_out);
				//}
                                //template_count++;

				double p1_score = candidate_progress1_score[choice];

				// if (m_peptide.full_grown() && p1_score < best_score)
				if (p1_score < best_score)
				{
					best_score = p1_score;
					best_conf = m_peptide.conf();
					is_best = true;
				}
			}
			else   // no candidate selected
			{
				m_move_failed = true;
				m_no_sel_count++;
			}

			m_curr_length_moves++;
			observer.after_move(this, candidate, candidate_score, choice,
				is_best, best_score);
			// added to get printout after every proposed move
			//char pdb_out[150];
			//sprintf(pdb_out,"%s_part%d_%ld",m_outfile.c_str(),m_peptide.length(),m_curr_length_moves);
			//m_peptide.write_pdb(pdb_out);
		}

		if (best_score < m_curr_score)
		{
			m_peptide.conf().swap(best_conf);
			m_curr_score = best_score;
		}

		m_strategy->end_run(this);
		observer.end_run(this);
	}
	observer.after_end(this);
}

const char *Runner::config_section()
{
	return m_config_section;
}

void Runner::print_template(std::ostream &out)
{
	out << "[" << capitalise(m_config_section) << "]\n"
		"\n"
        << c_param_sequential << " = " << bool_str(c_default_sequential)
            << "\t# whether is sequential or non-sequential\n"
        //<< c_param_coil << " = "
        //    << bool_str(c_default_coil)
        //    << "\t\t# non-sequential start state (coil vs fully extended)\n"
		<< c_param_move_limit << " = " << c_default_move_limit
			<< "\t\t# move limit (excluding growth phase)\n"
		<< "#" << c_param_no_sel_limit << " = " << c_default_no_sel_limit
			<< "\t# stop after none selected this many "
				"times in a row\n"
		//<< c_param_scwrl_executable << " = ...\t\tPath of SCWRL executable\n"
		<< c_param_reverse << " = false"
			<< "\t\t# whether to extrude in reverse order (C first)\n"
		<< "#" << c_param_start_struct
			<< " = ...\t\t# Starting structure (PDB file)\n"
		<< "#" << c_param_native_struct
			<< " = ...\t\t# Native structure (PDB file)\n"
		<< "\n";
}

