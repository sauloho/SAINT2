
EXEC=saint2
#EXEC=saint2_raw

#CPPFLAGS=-Wall -pedantic -g $(INCLUDES)
CPPFLAGS=-Wall -pedantic -O3 $(INCLUDES) -lstatic

#CPPFLAGS=-Wall -pedantic -O3 -DNDEBUG $(INCLUDES)
#CPPFLAGS=-Wall -pedantic -O3 -DNDEBUG $(INCLUDES) -DRAW_SCORE

#CPPFLAGS=-Wall -pedantic -g $(INCLUDES) -DRAW_SCORE

INCLUDES=-I main -I move -I score -I strategy -I peptide -I extend

MAIN=main/static_init.cpp main/common.cpp main/random.cpp main/c_file.cpp main/config.cpp main/param_list.cpp main/reporter.cpp main/runner.cpp main/stream_printf.cpp main/point.cpp main/matrix.cpp main/transform.cpp main/parse.cpp main/temp_file.cpp main/distribution.cpp main/geom.cpp main/rmsd.cpp
MOVE=move/mover.cpp move/mover_fragment.cpp move/mover_fragment_fwd.cpp move/mover_fragment_rev.cpp move/fragment.cpp
SCORE=score/scorer.cpp score/scorer_combined.cpp score/rapdf.cpp score/rapdf_impl.cpp score/solvation.cpp score/solvation_impl.cpp score/torsion.cpp score/torsion_impl.cpp score/hbond.cpp score/predtor.cpp score/saulo.cpp score/eleanor.cpp score/core.cpp score/core_impl.cpp score/predss.cpp score/rgyr.cpp score/contact.cpp score/crowding.cpp score/randomscr.cpp score/orientation.cpp score/orientation_impl.cpp score/lennard_jones.cpp score/ribosome.cpp 
STRATEGY=strategy/strategy.cpp strategy/strategy_strict.cpp strategy/strategy_monte.cpp strategy/strategy_boltz.cpp strategy/strategy_always.cpp
PEPTIDE=peptide/peptide.cpp peptide/residue.cpp peptide/atom.cpp peptide/sequence.cpp peptide/amino.cpp peptide/codon.cpp peptide/atom_type.cpp peptide/pdb_atom_rec.cpp peptide/conformation.cpp
EXTEND=extend/extender.cpp extend/extender_fixed.cpp extend/extender_codon.cpp

SRCS=$(MAIN) $(MOVE) $(STRATEGY) $(SCORE) $(PEPTIDE) $(EXTEND)
SAINT_SRCS=main/main.cpp $(SRCS)
MKDATA_SRCS=mkdata/main.cpp $(SRCS)
FRAG_SRCS=mkdata/convert_fragments.cpp $(SRCS)
CHECK_CHAINS_SRCS=mkdata/check_chains.cpp $(SRCS)
WRITE_CHAINS_SRCS=mkdata/write_chains.cpp $(SRCS)
GET_FASTA_SRCS=mkdata/get_fasta.cpp $(SRCS)
TRIPLETS_SRCS=mkdata/triplets.cpp $(SRCS)
DGEN_SRCS=decoygen/main.cpp $(SRCS)
DGEN_CREATE_DATA_SRCS=decoygen/create_dgen_data.cpp $(SRCS)
RMSD_SRCS=utils/get_pdb_rmsds.cpp $(SRCS)
FINDMIN_SRCS=findmin/findmin.cpp $(SRCS)
FILTER_LEN_SRCS=mkdata/filter_len_11.cpp $(SRCS)
WRITE_PDB_SRCS=write/main.cpp $(SRCS)
CLUSTER_SRCS=clustering/cluster.cpp $(SRCS)
TIMING_SRCS=timing/timing.cpp $(SRCS)
CIM_SRCS=timing/calc_invitro_moves.cpp $(SRCS)
BEND_SRCS=bend/bend.cpp $(SRCS)

SAINT_OBJS=$(SAINT_SRCS:.cpp=.o)
MKDATA_OBJS=$(MKDATA_SRCS:.cpp=.o)
FRAG_OBJS=$(FRAG_SRCS:.cpp=.o)
CHECK_CHAINS_OBJS=$(CHECK_CHAINS_SRCS:.cpp=.o)
WRITE_CHAINS_OBJS=$(WRITE_CHAINS_SRCS:.cpp=.o)
GET_FASTA_OBJS=$(GET_FASTA_SRCS:.cpp=.o)
TRIPLETS_OBJS=$(TRIPLETS_SRCS:.cpp=.o)
DGEN_OBJS=$(DGEN_SRCS:.cpp=.o)
DGEN_CREATE_DATA_OBJS=$(DGEN_CREATE_DATA_SRCS:.cpp=.o)
RMSD_OBJS=$(RMSD_SRCS:.cpp=.o)
FINDMIN_OBJS=$(FINDMIN_SRCS:.cpp=.o)
FILTER_LEN_OBJS=$(FILTER_LEN_SRCS:.cpp=.o)
WRITE_PDB_OBJS=$(WRITE_PDB_SRCS:.cpp=.o)
CLUSTER_OBJS=$(CLUSTER_SRCS:.cpp=.o)
TIMING_OBJS=$(TIMING_SRCS:.cpp=.o)
CIM_OBJS=$(CIM_SRCS:.cpp=.o)
BEND_OBJS=$(BEND_SRCS:.cpp=.o)
LIBS=-lstdc++ -lm

# standard include files to ignore in "make depend" output
STDINC="iostream"|"algorithm"|"cassert"|"string"|"cstring"|"cctype"|"ctime"|"cmath"|"cstdlib"|"cstdio"|"map"|"set"|"vector"|"list"|"sstream"|"csignal"

$(EXEC): $(SAINT_OBJS)
	gcc -o $(EXEC) $(SAINT_OBJS) $(LIBS)

mkd: $(MKDATA_OBJS)
	gcc -o mkd $(MKDATA_OBJS) $(LIBS)

convert_fragments: $(FRAG_OBJS)
	gcc -o convert_fragments $(FRAG_OBJS) $(LIBS)

check_chains: $(CHECK_CHAINS_OBJS)
	gcc -o check_chains $(CHECK_CHAINS_OBJS) $(LIBS)

write_chains: $(WRITE_CHAINS_OBJS)
	gcc -o write_chains $(WRITE_CHAINS_OBJS) $(LIBS)

triplets: $(TRIPLETS_OBJS)
	gcc -o triplets $(TRIPLETS_OBJS) $(LIBS)

dgen: $(DGEN_OBJS)
	gcc -o dgen $(DGEN_OBJS) $(LIBS)

create_dgen_data: $(DGEN_CREATE_DATA_OBJS)
	gcc -o create_dgen_data $(DGEN_CREATE_DATA_OBJS) $(LIBS)

get_pdb_rmsds: $(RMSD_OBJS)
	gcc -o get_pdb_rmsds $(RMSD_OBJS) $(LIBS)

cluster: $(CLUSTER_OBJS)
	gcc -o cluster $(CLUSTER_OBJS) $(LIBS)

time_moves: $(TIMING_OBJS)
	gcc -o time_moves $(TIMING_OBJS) $(LIBS)

calc_invitro_moves: $(CIM_OBJS)
	gcc -o calc_invitro_moves $(CIM_OBJS) $(LIBS)

find_min: $(FINDMIN_OBJS)
	gcc -o find_min $(FINDMIN_OBJS) $(LIBS)

write_pdb: $(WRITE_PDB_OBJS)
	gcc -o write_pdb $(WRITE_PDB_OBJS) $(LIBS)

get_fasta: $(GET_FASTA_OBJS)
	gcc -o get_fasta $(GET_FASTA_OBJS) $(LIBS)

filter_len_11: $(FILTER_LEN_OBJS)
	gcc -o filter_len_11 $(FILTER_LEN_OBJS) $(LIBS)

bend_pdb: $(BEND_OBJS)
	gcc -o bend_pdb $(BEND_OBJS) $(LIBS)

depend:
	makedepend -- $(CPPFLAGS) -- main/main.cpp mkdata/main.cpp decoygen/main.cpp decoygen/create_dgen_data.cpp mkdata/filter_len_11.cpp clustering/cluster.cpp timing/timing.cpp $(SRCS) $(CLUSTER_SRCS) 2>&1 | egrep -v 'not in |$(STDINC)'; cat Makefile | awk '/^# DO NOT/ { p = 1; next; } p == 1 { print; }'; echo

clean:
	rm -f $(EXEC) `find . -name '*.o' -print`

remake: clean $(EXEC)

MATRIX_TEST_OBJS=main/matrix.o main/point.o main/matrix_test.o
matrix_test: $(MATRIX_TEST_OBJS)
	gcc -o matrix_test  $(MATRIX_TEST_OBJS) $(LIBS)

# DO NOT DELETE THIS LINE -- make depend depends on it.

main/main.o: main/config.h main/param_list.h peptide/sequence.h
main/main.o: peptide/amino.h peptide/atom_id.h peptide/codon.h main/c_file.h
main/main.o: main/runner.h peptide/peptide.h peptide/residue.h peptide/atom.h
main/main.o: peptide/atom_type.h peptide/pdb_atom_rec.h main/point.h
main/main.o: main/transform.h main/matrix.h peptide/conformation.h
main/main.o: main/common.h score/scorer.h score/lennard_jones.h
main/main.o: main/reporter.h main/run_observer.h main/stream_printf.h
main/main.o: main/geom.h main/static_init.h
mkdata/main.o: main/config.h main/param_list.h peptide/sequence.h
mkdata/main.o: peptide/amino.h peptide/atom_id.h peptide/codon.h
mkdata/main.o: main/c_file.h main/runner.h peptide/peptide.h
mkdata/main.o: peptide/residue.h peptide/atom.h peptide/atom_type.h
mkdata/main.o: peptide/pdb_atom_rec.h main/point.h main/transform.h
mkdata/main.o: main/matrix.h peptide/conformation.h main/common.h
mkdata/main.o: score/scorer.h main/reporter.h main/run_observer.h
mkdata/main.o: main/static_init.h score/orientation.h score/torsion.h
mkdata/main.o: main/geom.h
decoygen/main.o: /usr/include/unistd.h /usr/include/features.h
decoygen/main.o: /usr/include/stdc-predef.h /usr/include/sys/cdefs.h
decoygen/main.o: /usr/include/bits/wordsize.h /usr/include/gnu/stubs.h
decoygen/main.o: /usr/include/bits/posix_opt.h
decoygen/main.o: /usr/include/bits/environments.h /usr/include/bits/types.h
decoygen/main.o: /usr/include/bits/typesizes.h /usr/include/bits/confname.h
decoygen/main.o: /usr/include/getopt.h score/torsion.h peptide/amino.h
decoygen/main.o: peptide/atom_id.h peptide/peptide.h peptide/residue.h
decoygen/main.o: peptide/codon.h peptide/atom.h peptide/atom_type.h
decoygen/main.o: peptide/pdb_atom_rec.h main/point.h main/transform.h
decoygen/main.o: main/matrix.h peptide/conformation.h main/common.h
decoygen/main.o: main/stream_printf.h main/c_file.h main/geom.h
decoygen/main.o: score/lennard_jones.h
decoygen/create_dgen_data.o: peptide/peptide.h peptide/residue.h
decoygen/create_dgen_data.o: peptide/amino.h peptide/atom_id.h
decoygen/create_dgen_data.o: peptide/codon.h peptide/atom.h
decoygen/create_dgen_data.o: peptide/atom_type.h peptide/pdb_atom_rec.h
decoygen/create_dgen_data.o: main/point.h main/transform.h main/matrix.h
decoygen/create_dgen_data.o: peptide/conformation.h main/common.h
decoygen/create_dgen_data.o: score/torsion.h main/c_file.h
mkdata/filter_len_11.o: main/c_file.h
clustering/cluster.o: main/c_file.h clustering/cluster_template.h
timing/timing.o: /usr/include/sys/times.h /usr/include/features.h
timing/timing.o: /usr/include/stdc-predef.h /usr/include/sys/cdefs.h
timing/timing.o: /usr/include/bits/wordsize.h /usr/include/gnu/stubs.h
timing/timing.o: /usr/include/time.h /usr/include/bits/types.h
timing/timing.o: /usr/include/bits/typesizes.h main/config.h
timing/timing.o: main/param_list.h main/runner.h peptide/peptide.h
timing/timing.o: peptide/residue.h peptide/amino.h peptide/atom_id.h
timing/timing.o: peptide/codon.h peptide/atom.h peptide/atom_type.h
timing/timing.o: peptide/pdb_atom_rec.h main/point.h main/transform.h
timing/timing.o: main/matrix.h peptide/conformation.h main/common.h
timing/timing.o: main/run_observer.h peptide/sequence.h main/c_file.h
main/static_init.o: main/static_init.h peptide/amino.h peptide/atom_id.h
main/static_init.o: peptide/codon.h peptide/atom_type.h
main/common.o: main/config.h main/param_list.h main/common.h
main/random.o: main/random.h
main/c_file.o: main/c_file.h main/config.h main/param_list.h
main/config.o: main/config.h main/param_list.h main/random.h
main/config.o: peptide/sequence.h peptide/amino.h peptide/atom_id.h
main/config.o: peptide/codon.h main/c_file.h main/runner.h peptide/peptide.h
main/config.o: peptide/residue.h peptide/atom.h peptide/atom_type.h
main/config.o: peptide/pdb_atom_rec.h main/point.h main/transform.h
main/config.o: main/matrix.h peptide/conformation.h main/common.h
main/config.o: extend/extender.h score/scorer.h strategy/strategy.h
main/config.o: main/stream_printf.h move/mover.h
main/param_list.o: main/config.h main/param_list.h
main/reporter.o: main/geom.h main/point.h main/runner.h main/param_list.h
main/reporter.o: peptide/peptide.h peptide/residue.h peptide/amino.h
main/reporter.o: peptide/atom_id.h peptide/codon.h peptide/atom.h
main/reporter.o: peptide/atom_type.h peptide/pdb_atom_rec.h main/transform.h
main/reporter.o: main/matrix.h peptide/conformation.h main/common.h
main/reporter.o: main/config.h main/reporter.h main/run_observer.h
main/runner.o: /usr/include/stdio.h /usr/include/features.h
main/runner.o: /usr/include/stdc-predef.h /usr/include/sys/cdefs.h
main/runner.o: /usr/include/bits/wordsize.h /usr/include/gnu/stubs.h
main/runner.o: /usr/include/bits/types.h /usr/include/bits/typesizes.h
main/runner.o: /usr/include/libio.h /usr/include/_G_config.h
main/runner.o: /usr/include/wchar.h /usr/include/bits/stdio_lim.h
main/runner.o: /usr/include/bits/sys_errlist.h main/runner.h
main/runner.o: main/param_list.h peptide/peptide.h peptide/residue.h
main/runner.o: peptide/amino.h peptide/atom_id.h peptide/codon.h
main/runner.o: peptide/atom.h peptide/atom_type.h peptide/pdb_atom_rec.h
main/runner.o: main/point.h main/transform.h main/matrix.h
main/runner.o: peptide/conformation.h main/common.h main/run_observer.h
main/runner.o: main/config.h score/scorer.h score/scorer_combined.h
main/runner.o: strategy/strategy.h extend/extender.h move/mover.h
main/stream_printf.o: main/config.h main/param_list.h main/stream_printf.h
main/point.o: main/point.h
main/matrix.o: main/matrix.h main/point.h
main/transform.o: main/matrix.h main/point.h main/transform.h
main/parse.o: main/parse.h
main/temp_file.o: /usr/include/unistd.h /usr/include/features.h
main/temp_file.o: /usr/include/stdc-predef.h /usr/include/sys/cdefs.h
main/temp_file.o: /usr/include/bits/wordsize.h /usr/include/gnu/stubs.h
main/temp_file.o: /usr/include/bits/posix_opt.h
main/temp_file.o: /usr/include/bits/environments.h /usr/include/bits/types.h
main/temp_file.o: /usr/include/bits/typesizes.h /usr/include/bits/confname.h
main/temp_file.o: /usr/include/getopt.h main/temp_file.h
main/distribution.o: main/random.h main/distribution.h
main/geom.o: main/matrix.h main/point.h main/geom.h
main/rmsd.o: main/rmsd.h
move/mover.o: main/config.h main/param_list.h main/common.h move/mover.h
move/mover.o: peptide/peptide.h peptide/residue.h peptide/amino.h
move/mover.o: peptide/atom_id.h peptide/codon.h peptide/atom.h
move/mover.o: peptide/atom_type.h peptide/pdb_atom_rec.h main/point.h
move/mover.o: main/transform.h main/matrix.h peptide/conformation.h
move/mover.o: move/mover_fragment.h move/fragment.h /usr/include/assert.h
move/mover.o: /usr/include/features.h /usr/include/stdc-predef.h
move/mover.o: /usr/include/sys/cdefs.h /usr/include/bits/wordsize.h
move/mover.o: /usr/include/gnu/stubs.h move/mover_fragment_fwd.h
move/mover.o: main/distribution.h move/mover_fragment_rev.h
move/mover_fragment.o: move/mover.h peptide/peptide.h peptide/residue.h
move/mover_fragment.o: peptide/amino.h peptide/atom_id.h peptide/codon.h
move/mover_fragment.o: peptide/atom.h peptide/atom_type.h
move/mover_fragment.o: peptide/pdb_atom_rec.h main/point.h main/transform.h
move/mover_fragment.o: main/matrix.h peptide/conformation.h main/common.h
move/mover_fragment.o: main/param_list.h move/mover_fragment.h
move/mover_fragment.o: move/fragment.h /usr/include/assert.h
move/mover_fragment.o: /usr/include/features.h /usr/include/stdc-predef.h
move/mover_fragment.o: /usr/include/sys/cdefs.h /usr/include/bits/wordsize.h
move/mover_fragment.o: /usr/include/gnu/stubs.h main/run_observer.h
move/mover_fragment.o: main/c_file.h main/config.h main/runner.h
move/mover_fragment.o: main/random.h main/geom.h main/stream_printf.h
move/mover_fragment_fwd.o: move/mover.h peptide/peptide.h peptide/residue.h
move/mover_fragment_fwd.o: peptide/amino.h peptide/atom_id.h peptide/codon.h
move/mover_fragment_fwd.o: peptide/atom.h peptide/atom_type.h
move/mover_fragment_fwd.o: peptide/pdb_atom_rec.h main/point.h
move/mover_fragment_fwd.o: main/transform.h main/matrix.h
move/mover_fragment_fwd.o: peptide/conformation.h main/common.h
move/mover_fragment_fwd.o: main/param_list.h move/mover_fragment_fwd.h
move/mover_fragment_fwd.o: move/mover_fragment.h move/fragment.h
move/mover_fragment_fwd.o: /usr/include/assert.h /usr/include/features.h
move/mover_fragment_fwd.o: /usr/include/stdc-predef.h
move/mover_fragment_fwd.o: /usr/include/sys/cdefs.h
move/mover_fragment_fwd.o: /usr/include/bits/wordsize.h
move/mover_fragment_fwd.o: /usr/include/gnu/stubs.h main/distribution.h
move/mover_fragment_fwd.o: main/run_observer.h main/c_file.h main/config.h
move/mover_fragment_fwd.o: main/runner.h main/random.h main/geom.h
move/mover_fragment_fwd.o: main/stream_printf.h
move/mover_fragment_rev.o: move/mover.h peptide/peptide.h peptide/residue.h
move/mover_fragment_rev.o: peptide/amino.h peptide/atom_id.h peptide/codon.h
move/mover_fragment_rev.o: peptide/atom.h peptide/atom_type.h
move/mover_fragment_rev.o: peptide/pdb_atom_rec.h main/point.h
move/mover_fragment_rev.o: main/transform.h main/matrix.h
move/mover_fragment_rev.o: peptide/conformation.h main/common.h
move/mover_fragment_rev.o: main/param_list.h move/mover_fragment_rev.h
move/mover_fragment_rev.o: move/mover_fragment.h move/fragment.h
move/mover_fragment_rev.o: /usr/include/assert.h /usr/include/features.h
move/mover_fragment_rev.o: /usr/include/stdc-predef.h
move/mover_fragment_rev.o: /usr/include/sys/cdefs.h
move/mover_fragment_rev.o: /usr/include/bits/wordsize.h
move/mover_fragment_rev.o: /usr/include/gnu/stubs.h main/distribution.h
move/mover_fragment_rev.o: main/run_observer.h main/c_file.h main/config.h
move/mover_fragment_rev.o: main/runner.h main/random.h main/geom.h
move/mover_fragment_rev.o: main/stream_printf.h
move/fragment.o: move/fragment.h /usr/include/assert.h
move/fragment.o: /usr/include/features.h /usr/include/stdc-predef.h
move/fragment.o: /usr/include/sys/cdefs.h /usr/include/bits/wordsize.h
move/fragment.o: /usr/include/gnu/stubs.h
strategy/strategy.o: main/config.h main/param_list.h main/stream_printf.h
strategy/strategy.o: strategy/strategy.h main/common.h
strategy/strategy.o: strategy/strategy_strict.h strategy/strategy_monte.h
strategy/strategy.o: strategy/strategy_boltz.h strategy/strategy_always.h
strategy/strategy_strict.o: main/common.h strategy/strategy_strict.h
strategy/strategy_strict.o: strategy/strategy.h main/param_list.h
strategy/strategy_monte.o: main/config.h main/param_list.h
strategy/strategy_monte.o: strategy/strategy_monte.h main/common.h
strategy/strategy_monte.o: strategy/strategy.h main/random.h
strategy/strategy_monte.o: main/stream_printf.h
strategy/strategy_boltz.o: main/config.h main/param_list.h main/random.h
strategy/strategy_boltz.o: main/stream_printf.h strategy/strategy_boltz.h
strategy/strategy_boltz.o: main/common.h strategy/strategy.h
strategy/strategy_always.o: main/common.h strategy/strategy_always.h
strategy/strategy_always.o: strategy/strategy.h main/param_list.h
score/scorer.o: score/scorer.h main/param_list.h score/scorer_combined.h
score/scorer.o: main/common.h main/config.h
score/scorer_combined.o: score/scorer_combined.h score/scorer.h
score/scorer_combined.o: main/param_list.h main/common.h peptide/atom.h
score/scorer_combined.o: peptide/amino.h peptide/atom_id.h
score/scorer_combined.o: peptide/atom_type.h peptide/pdb_atom_rec.h
score/scorer_combined.o: main/point.h peptide/peptide.h peptide/residue.h
score/scorer_combined.o: peptide/codon.h main/transform.h main/matrix.h
score/scorer_combined.o: peptide/conformation.h main/c_file.h
score/scorer_combined.o: main/stream_printf.h score/rapdf.h score/solvation.h
score/scorer_combined.o: score/torsion.h score/hbond.h score/saulo.h score/eleanor.h score/predss.h
score/scorer_combined.o: score/rgyr.h score/contact.h score/crowding.h score/predtor.h
score/scorer_combined.o: score/randomscr.h score/orientation.h score/core.h
score/scorer_combined.o: score/lennard_jones.h score/ribosome.h
score/rapdf.o: score/rapdf.h peptide/amino.h peptide/atom_id.h
score/rapdf.o: score/rapdf_impl.h peptide/peptide.h peptide/residue.h
score/rapdf.o: peptide/codon.h peptide/atom.h peptide/atom_type.h
score/rapdf.o: peptide/pdb_atom_rec.h main/point.h main/transform.h
score/rapdf.o: main/matrix.h peptide/conformation.h main/common.h
score/rapdf.o: score/scorer_combined.h score/scorer.h main/param_list.h
score/rapdf_impl.o: peptide/atom.h peptide/amino.h peptide/atom_id.h
score/rapdf_impl.o: peptide/atom_type.h peptide/pdb_atom_rec.h main/point.h
score/rapdf_impl.o: peptide/peptide.h peptide/residue.h peptide/codon.h
score/rapdf_impl.o: main/transform.h main/matrix.h peptide/conformation.h
score/rapdf_impl.o: main/common.h main/c_file.h score/scorer_combined.h
score/rapdf_impl.o: score/scorer.h main/param_list.h score/rapdf_impl.h
score/rapdf_impl.o: score/rapdf.h
score/solvation.o: score/solvation.h score/solvation_impl.h peptide/peptide.h
score/solvation.o: peptide/residue.h peptide/amino.h peptide/atom_id.h
score/solvation.o: peptide/codon.h peptide/atom.h peptide/atom_type.h
score/solvation.o: peptide/pdb_atom_rec.h main/point.h main/transform.h
score/solvation.o: main/matrix.h peptide/conformation.h main/common.h
score/solvation.o: score/scorer_combined.h score/scorer.h main/param_list.h
score/solvation_impl.o: peptide/peptide.h peptide/residue.h peptide/amino.h
score/solvation_impl.o: peptide/atom_id.h peptide/codon.h peptide/atom.h
score/solvation_impl.o: peptide/atom_type.h peptide/pdb_atom_rec.h
score/solvation_impl.o: main/point.h main/transform.h main/matrix.h
score/solvation_impl.o: peptide/conformation.h main/common.h main/c_file.h
score/solvation_impl.o: score/scorer_combined.h score/scorer.h
score/solvation_impl.o: main/param_list.h score/solvation.h
score/solvation_impl.o: score/solvation_impl.h
score/torsion.o: score/torsion.h score/torsion_impl.h peptide/amino.h
score/torsion.o: peptide/atom_id.h peptide/peptide.h peptide/residue.h
score/torsion.o: peptide/codon.h peptide/atom.h peptide/atom_type.h
score/torsion.o: peptide/pdb_atom_rec.h main/point.h main/transform.h
score/torsion.o: main/matrix.h peptide/conformation.h main/common.h
score/torsion.o: main/geom.h score/scorer_combined.h score/scorer.h
score/torsion.o: main/param_list.h
score/torsion_impl.o: peptide/peptide.h peptide/residue.h peptide/amino.h
score/torsion_impl.o: peptide/atom_id.h peptide/codon.h peptide/atom.h
score/torsion_impl.o: peptide/atom_type.h peptide/pdb_atom_rec.h main/point.h
score/torsion_impl.o: main/transform.h main/matrix.h peptide/conformation.h
score/torsion_impl.o: main/common.h main/geom.h main/c_file.h
score/torsion_impl.o: main/stream_printf.h score/scorer_combined.h
score/torsion_impl.o: score/scorer.h main/param_list.h score/torsion.h
score/torsion_impl.o: score/torsion_impl.h
score/hbond.o: peptide/peptide.h peptide/residue.h peptide/amino.h
score/hbond.o: peptide/atom_id.h peptide/codon.h peptide/atom.h
score/hbond.o: peptide/atom_type.h peptide/pdb_atom_rec.h main/point.h
score/hbond.o: main/transform.h main/matrix.h peptide/conformation.h
score/hbond.o: main/common.h score/scorer_combined.h score/scorer.h
score/hbond.o: main/param_list.h score/hbond.h
score/saulo.o: peptide/amino.h peptide/atom_id.h peptide/residue.h
score/saulo.o: peptide/codon.h peptide/atom.h peptide/atom_type.h
score/saulo.o: peptide/pdb_atom_rec.h main/point.h peptide/peptide.h
score/saulo.o: main/transform.h main/matrix.h peptide/conformation.h
score/saulo.o: main/common.h score/scorer_combined.h score/scorer.h
score/eleanor.o: main/param_list.h score/saulo.h
score/eleanor.o: peptide/amino.h peptide/atom_id.h peptide/residue.h
score/eleanor.o: peptide/codon.h peptide/atom.h peptide/atom_type.h
score/eleanor.o: peptide/pdb_atom_rec.h main/point.h peptide/peptide.h
score/eleanor.o: main/transform.h main/matrix.h peptide/conformation.h
score/eleanor.o: main/common.h score/scorer_combined.h score/scorer.h
score/eleanor.o: main/param_list.h score/eleanor.h
score/predtor.o: peptide/amino.h peptide/atom_id.h peptide/residue.h
score/predtor.o: peptide/codon.h peptide/atom.h peptide/atom_type.h
score/predtor.o: peptide/pdb_atom_rec.h main/point.h peptide/peptide.h
score/predtor.o: main/transform.h main/matrix.h peptide/conformation.h
score/predtor.o: main/common.h score/scorer_combined.h score/scorer.h
score/predtor.o: main/param_list.h score/predtor.h
score/core.o: peptide/amino.h peptide/atom_id.h peptide/residue.h
score/core.o: peptide/codon.h peptide/atom.h peptide/atom_type.h
score/core.o: peptide/pdb_atom_rec.h main/point.h peptide/peptide.h
score/core.o: main/transform.h main/matrix.h peptide/conformation.h
score/core.o: main/common.h score/scorer_combined.h score/scorer.h
score/core.o: main/param_list.h score/core.h
score/core_impl.o: peptide/peptide.h peptide/residue.h peptide/amino.h
score/core_impl.o: peptide/atom_id.h peptide/codon.h peptide/atom.h
score/core_impl.o: peptide/atom_type.h peptide/pdb_atom_rec.h main/point.h
score/core_impl.o: main/transform.h main/matrix.h peptide/conformation.h
score/core_impl.o: main/common.h main/geom.h main/c_file.h
score/core_impl.o: main/stream_printf.h score/scorer_combined.h
score/core_impl.o: score/scorer.h main/param_list.h score/torsion.h
score/core_impl.o: score/torsion_impl.h
score/predss.o: peptide/amino.h peptide/atom_id.h peptide/residue.h
score/predss.o: peptide/codon.h peptide/atom.h peptide/atom_type.h
score/predss.o: peptide/pdb_atom_rec.h main/point.h peptide/peptide.h
score/predss.o: main/transform.h main/matrix.h peptide/conformation.h
score/predss.o: main/common.h score/scorer_combined.h score/scorer.h
score/predss.o: main/param_list.h score/predss.h
score/rgyr.o: peptide/peptide.h peptide/residue.h peptide/amino.h
score/rgyr.o: peptide/atom_id.h peptide/codon.h peptide/atom.h
score/rgyr.o: peptide/atom_type.h peptide/pdb_atom_rec.h main/point.h
score/rgyr.o: main/transform.h main/matrix.h peptide/conformation.h
score/rgyr.o: main/common.h score/scorer_combined.h score/scorer.h
score/rgyr.o: main/param_list.h score/rgyr.h
score/contact.o: peptide/peptide.h peptide/residue.h peptide/amino.h
score/contact.o: peptide/atom_id.h peptide/codon.h peptide/atom.h
score/contact.o: peptide/atom_type.h peptide/pdb_atom_rec.h main/point.h
score/contact.o: main/transform.h main/matrix.h peptide/conformation.h
score/contact.o: main/common.h score/scorer_combined.h score/scorer.h
score/contact.o: main/param_list.h score/contact.h
score/crowding.o: peptide/peptide.h peptide/residue.h peptide/amino.h
score/crowding.o: peptide/atom_id.h peptide/codon.h peptide/atom.h
score/crowding.o: peptide/atom_type.h peptide/pdb_atom_rec.h main/point.h
score/crowding.o: main/transform.h main/matrix.h peptide/conformation.h
score/crowding.o: main/common.h score/scorer_combined.h score/scorer.h
score/crowding.o: main/param_list.h score/crowding.h
score/randomscr.o: peptide/peptide.h peptide/residue.h peptide/amino.h
score/randomscr.o: peptide/atom_id.h peptide/codon.h peptide/atom.h
score/randomscr.o: peptide/atom_type.h peptide/pdb_atom_rec.h main/point.h
score/randomscr.o: main/transform.h main/matrix.h peptide/conformation.h
score/randomscr.o: main/common.h score/scorer_combined.h score/scorer.h
score/randomscr.o: main/param_list.h score/randomscr.h
score/orientation.o: score/orientation.h peptide/amino.h peptide/atom_id.h
score/orientation.o: score/orientation_impl.h peptide/peptide.h
score/orientation.o: peptide/residue.h peptide/codon.h peptide/atom.h
score/orientation.o: peptide/atom_type.h peptide/pdb_atom_rec.h main/point.h
score/orientation.o: main/transform.h main/matrix.h peptide/conformation.h
score/orientation.o: main/common.h main/geom.h score/scorer_combined.h
score/orientation.o: score/scorer.h main/param_list.h
score/orientation_impl.o: peptide/peptide.h peptide/residue.h peptide/amino.h
score/orientation_impl.o: peptide/atom_id.h peptide/codon.h peptide/atom.h
score/orientation_impl.o: peptide/atom_type.h peptide/pdb_atom_rec.h
score/orientation_impl.o: main/point.h main/transform.h main/matrix.h
score/orientation_impl.o: peptide/conformation.h main/common.h main/geom.h
score/orientation_impl.o: main/c_file.h score/scorer_combined.h
score/orientation_impl.o: score/scorer.h main/param_list.h
score/orientation_impl.o: score/orientation.h score/orientation_impl.h
score/lennard_jones.o: peptide/peptide.h peptide/residue.h peptide/amino.h
score/lennard_jones.o: peptide/atom_id.h peptide/codon.h peptide/atom.h
score/lennard_jones.o: peptide/atom_type.h peptide/pdb_atom_rec.h
score/lennard_jones.o: main/point.h main/transform.h main/matrix.h
score/lennard_jones.o: peptide/conformation.h main/common.h
score/lennard_jones.o: score/scorer_combined.h score/scorer.h
score/lennard_jones.o: main/param_list.h score/lennard_jones.h
score/ribosome.o: main/common.h peptide/peptide.h peptide/residue.h
score/ribosome.o: peptide/amino.h peptide/atom_id.h peptide/codon.h
score/ribosome.o: peptide/atom.h peptide/atom_type.h peptide/pdb_atom_rec.h
score/ribosome.o: main/point.h main/transform.h main/matrix.h
score/ribosome.o: peptide/conformation.h score/ribosome.h
peptide/peptide.o: /usr/include/sys/stat.h /usr/include/features.h
peptide/peptide.o: /usr/include/stdc-predef.h /usr/include/sys/cdefs.h
peptide/peptide.o: /usr/include/bits/wordsize.h /usr/include/gnu/stubs.h
peptide/peptide.o: /usr/include/bits/types.h /usr/include/bits/typesizes.h
peptide/peptide.o: /usr/include/time.h /usr/include/bits/stat.h
peptide/peptide.o: /usr/include/sys/types.h /usr/include/endian.h
peptide/peptide.o: /usr/include/bits/endian.h /usr/include/bits/byteswap.h
peptide/peptide.o: /usr/include/bits/byteswap-16.h /usr/include/sys/select.h
peptide/peptide.o: /usr/include/bits/select.h /usr/include/bits/sigset.h
peptide/peptide.o: /usr/include/bits/time.h /usr/include/sys/sysmacros.h
peptide/peptide.o: /usr/include/bits/pthreadtypes.h peptide/peptide.h
peptide/peptide.o: peptide/residue.h peptide/amino.h peptide/atom_id.h
peptide/peptide.o: peptide/codon.h peptide/atom.h peptide/atom_type.h
peptide/peptide.o: peptide/pdb_atom_rec.h main/point.h main/transform.h
peptide/peptide.o: main/matrix.h peptide/conformation.h main/common.h
peptide/peptide.o: peptide/sequence.h main/param_list.h main/c_file.h
peptide/peptide.o: main/config.h main/temp_file.h main/stream_printf.h
peptide/peptide.o: main/parse.h main/rmsd.h main/geom.h
peptide/residue.o: peptide/residue.h peptide/amino.h peptide/atom_id.h
peptide/residue.o: peptide/codon.h peptide/atom.h peptide/atom_type.h
peptide/residue.o: peptide/pdb_atom_rec.h main/point.h
peptide/atom.o: peptide/atom.h peptide/amino.h peptide/atom_id.h
peptide/atom.o: peptide/atom_type.h peptide/pdb_atom_rec.h main/point.h
peptide/sequence.o: main/common.h main/c_file.h peptide/peptide.h
peptide/sequence.o: peptide/residue.h peptide/amino.h peptide/atom_id.h
peptide/sequence.o: peptide/codon.h peptide/atom.h peptide/atom_type.h
peptide/sequence.o: peptide/pdb_atom_rec.h main/point.h main/transform.h
peptide/sequence.o: main/matrix.h peptide/conformation.h peptide/sequence.h
peptide/sequence.o: main/param_list.h main/config.h
peptide/amino.o: /usr/include/strings.h /usr/include/features.h
peptide/amino.o: /usr/include/stdc-predef.h /usr/include/sys/cdefs.h
peptide/amino.o: /usr/include/bits/wordsize.h /usr/include/gnu/stubs.h
peptide/amino.o: /usr/include/xlocale.h peptide/amino.h peptide/atom_id.h
peptide/amino.o: peptide/atom_type.h
peptide/codon.o: peptide/amino.h peptide/atom_id.h peptide/codon.h
peptide/atom_type.o: peptide/atom_type.h peptide/atom_id.h
peptide/pdb_atom_rec.o: /usr/include/stdio.h /usr/include/features.h
peptide/pdb_atom_rec.o: /usr/include/stdc-predef.h /usr/include/sys/cdefs.h
peptide/pdb_atom_rec.o: /usr/include/bits/wordsize.h /usr/include/gnu/stubs.h
peptide/pdb_atom_rec.o: /usr/include/bits/types.h
peptide/pdb_atom_rec.o: /usr/include/bits/typesizes.h /usr/include/libio.h
peptide/pdb_atom_rec.o: /usr/include/_G_config.h /usr/include/wchar.h
peptide/pdb_atom_rec.o: /usr/include/bits/stdio_lim.h
peptide/pdb_atom_rec.o: /usr/include/bits/sys_errlist.h
peptide/pdb_atom_rec.o: peptide/pdb_atom_rec.h main/point.h peptide/amino.h
peptide/pdb_atom_rec.o: peptide/atom_id.h peptide/atom_type.h
peptide/conformation.o: peptide/peptide.h peptide/residue.h peptide/amino.h
peptide/conformation.o: peptide/atom_id.h peptide/codon.h peptide/atom.h
peptide/conformation.o: peptide/atom_type.h peptide/pdb_atom_rec.h
peptide/conformation.o: main/point.h main/transform.h main/matrix.h
peptide/conformation.o: peptide/conformation.h main/common.h main/geom.h
extend/extender.o: extend/extender.h main/param_list.h
extend/extender.o: extend/extender_fixed.h extend/extender_codon.h
extend/extender.o: peptide/codon.h peptide/amino.h peptide/atom_id.h
extend/extender.o: main/common.h main/config.h main/runner.h
extend/extender.o: peptide/peptide.h peptide/residue.h peptide/atom.h
extend/extender.o: peptide/atom_type.h peptide/pdb_atom_rec.h main/point.h
extend/extender.o: main/transform.h main/matrix.h peptide/conformation.h
extend/extender_fixed.o: main/config.h main/param_list.h peptide/sequence.h
extend/extender_fixed.o: peptide/amino.h peptide/atom_id.h peptide/codon.h
extend/extender_fixed.o: main/c_file.h extend/extender.h
extend/extender_fixed.o: extend/extender_fixed.h
extend/extender_codon.o: main/config.h main/param_list.h main/common.h
extend/extender_codon.o: main/c_file.h main/stream_printf.h peptide/codon.h
extend/extender_codon.o: peptide/amino.h peptide/atom_id.h main/runner.h
extend/extender_codon.o: peptide/peptide.h peptide/residue.h peptide/atom.h
extend/extender_codon.o: peptide/atom_type.h peptide/pdb_atom_rec.h
extend/extender_codon.o: main/point.h main/transform.h main/matrix.h
extend/extender_codon.o: peptide/conformation.h peptide/sequence.h
extend/extender_codon.o: extend/extender.h extend/extender_codon.h
clustering/cluster.o: main/c_file.h clustering/cluster_template.h
main/static_init.o: main/static_init.h peptide/amino.h peptide/atom_id.h
main/static_init.o: peptide/codon.h peptide/atom_type.h
main/common.o: main/config.h main/param_list.h main/common.h
main/random.o: main/random.h
main/c_file.o: main/c_file.h main/config.h main/param_list.h
main/config.o: main/config.h main/param_list.h main/random.h
main/config.o: peptide/sequence.h peptide/amino.h peptide/atom_id.h
main/config.o: peptide/codon.h main/c_file.h main/runner.h peptide/peptide.h
main/config.o: peptide/residue.h peptide/atom.h peptide/atom_type.h
main/config.o: peptide/pdb_atom_rec.h main/point.h main/transform.h
main/config.o: main/matrix.h peptide/conformation.h main/common.h
main/config.o: extend/extender.h score/scorer.h strategy/strategy.h
main/config.o: main/stream_printf.h move/mover.h
main/param_list.o: main/config.h main/param_list.h
main/reporter.o: main/geom.h main/point.h main/runner.h main/param_list.h
main/reporter.o: peptide/peptide.h peptide/residue.h peptide/amino.h
main/reporter.o: peptide/atom_id.h peptide/codon.h peptide/atom.h
main/reporter.o: peptide/atom_type.h peptide/pdb_atom_rec.h main/transform.h
main/reporter.o: main/matrix.h peptide/conformation.h main/common.h
main/reporter.o: main/config.h main/reporter.h main/run_observer.h
main/runner.o: /usr/include/stdio.h /usr/include/features.h
main/runner.o: /usr/include/stdc-predef.h /usr/include/sys/cdefs.h
main/runner.o: /usr/include/bits/wordsize.h /usr/include/gnu/stubs.h
main/runner.o: /usr/include/bits/types.h /usr/include/bits/typesizes.h
main/runner.o: /usr/include/libio.h /usr/include/_G_config.h
main/runner.o: /usr/include/wchar.h /usr/include/bits/stdio_lim.h
main/runner.o: /usr/include/bits/sys_errlist.h main/runner.h
main/runner.o: main/param_list.h peptide/peptide.h peptide/residue.h
main/runner.o: peptide/amino.h peptide/atom_id.h peptide/codon.h
main/runner.o: peptide/atom.h peptide/atom_type.h peptide/pdb_atom_rec.h
main/runner.o: main/point.h main/transform.h main/matrix.h
main/runner.o: peptide/conformation.h main/common.h main/run_observer.h
main/runner.o: main/config.h score/scorer.h score/scorer_combined.h
main/runner.o: strategy/strategy.h extend/extender.h move/mover.h
main/stream_printf.o: main/config.h main/param_list.h main/stream_printf.h
main/point.o: main/point.h
main/matrix.o: main/matrix.h main/point.h
main/transform.o: main/matrix.h main/point.h main/transform.h
main/parse.o: main/parse.h
main/temp_file.o: /usr/include/unistd.h /usr/include/features.h
main/temp_file.o: /usr/include/stdc-predef.h /usr/include/sys/cdefs.h
main/temp_file.o: /usr/include/bits/wordsize.h /usr/include/gnu/stubs.h
main/temp_file.o: /usr/include/bits/posix_opt.h
main/temp_file.o: /usr/include/bits/environments.h /usr/include/bits/types.h
main/temp_file.o: /usr/include/bits/typesizes.h /usr/include/bits/confname.h
main/temp_file.o: /usr/include/getopt.h main/temp_file.h
main/distribution.o: main/random.h main/distribution.h
main/geom.o: main/matrix.h main/point.h main/geom.h
main/rmsd.o: main/rmsd.h
move/mover.o: main/config.h main/param_list.h main/common.h move/mover.h
move/mover.o: peptide/peptide.h peptide/residue.h peptide/amino.h
move/mover.o: peptide/atom_id.h peptide/codon.h peptide/atom.h
move/mover.o: peptide/atom_type.h peptide/pdb_atom_rec.h main/point.h
move/mover.o: main/transform.h main/matrix.h peptide/conformation.h
move/mover.o: move/mover_fragment.h move/fragment.h /usr/include/assert.h
move/mover.o: /usr/include/features.h /usr/include/stdc-predef.h
move/mover.o: /usr/include/sys/cdefs.h /usr/include/bits/wordsize.h
move/mover.o: /usr/include/gnu/stubs.h move/mover_fragment_fwd.h
move/mover.o: main/distribution.h move/mover_fragment_rev.h
move/mover_fragment.o: move/mover.h peptide/peptide.h peptide/residue.h
move/mover_fragment.o: peptide/amino.h peptide/atom_id.h peptide/codon.h
move/mover_fragment.o: peptide/atom.h peptide/atom_type.h
move/mover_fragment.o: peptide/pdb_atom_rec.h main/point.h main/transform.h
move/mover_fragment.o: main/matrix.h peptide/conformation.h main/common.h
move/mover_fragment.o: main/param_list.h move/mover_fragment.h
move/mover_fragment.o: move/fragment.h /usr/include/assert.h
move/mover_fragment.o: /usr/include/features.h /usr/include/stdc-predef.h
move/mover_fragment.o: /usr/include/sys/cdefs.h /usr/include/bits/wordsize.h
move/mover_fragment.o: /usr/include/gnu/stubs.h main/run_observer.h
move/mover_fragment.o: main/c_file.h main/config.h main/runner.h
move/mover_fragment.o: main/random.h main/geom.h main/stream_printf.h
move/mover_fragment_fwd.o: move/mover.h peptide/peptide.h peptide/residue.h
move/mover_fragment_fwd.o: peptide/amino.h peptide/atom_id.h peptide/codon.h
move/mover_fragment_fwd.o: peptide/atom.h peptide/atom_type.h
move/mover_fragment_fwd.o: peptide/pdb_atom_rec.h main/point.h
move/mover_fragment_fwd.o: main/transform.h main/matrix.h
move/mover_fragment_fwd.o: peptide/conformation.h main/common.h
move/mover_fragment_fwd.o: main/param_list.h move/mover_fragment_fwd.h
move/mover_fragment_fwd.o: move/mover_fragment.h move/fragment.h
move/mover_fragment_fwd.o: /usr/include/assert.h /usr/include/features.h
move/mover_fragment_fwd.o: /usr/include/stdc-predef.h
move/mover_fragment_fwd.o: /usr/include/sys/cdefs.h
move/mover_fragment_fwd.o: /usr/include/bits/wordsize.h
move/mover_fragment_fwd.o: /usr/include/gnu/stubs.h main/distribution.h
move/mover_fragment_fwd.o: main/run_observer.h main/c_file.h main/config.h
move/mover_fragment_fwd.o: main/runner.h main/random.h main/geom.h
move/mover_fragment_fwd.o: main/stream_printf.h
move/mover_fragment_rev.o: move/mover.h peptide/peptide.h peptide/residue.h
move/mover_fragment_rev.o: peptide/amino.h peptide/atom_id.h peptide/codon.h
move/mover_fragment_rev.o: peptide/atom.h peptide/atom_type.h
move/mover_fragment_rev.o: peptide/pdb_atom_rec.h main/point.h
move/mover_fragment_rev.o: main/transform.h main/matrix.h
move/mover_fragment_rev.o: peptide/conformation.h main/common.h
move/mover_fragment_rev.o: main/param_list.h move/mover_fragment_rev.h
move/mover_fragment_rev.o: move/mover_fragment.h move/fragment.h
move/mover_fragment_rev.o: /usr/include/assert.h /usr/include/features.h
move/mover_fragment_rev.o: /usr/include/stdc-predef.h
move/mover_fragment_rev.o: /usr/include/sys/cdefs.h
move/mover_fragment_rev.o: /usr/include/bits/wordsize.h
move/mover_fragment_rev.o: /usr/include/gnu/stubs.h main/distribution.h
move/mover_fragment_rev.o: main/run_observer.h main/c_file.h main/config.h
move/mover_fragment_rev.o: main/runner.h main/random.h main/geom.h
move/mover_fragment_rev.o: main/stream_printf.h
move/fragment.o: move/fragment.h /usr/include/assert.h
move/fragment.o: /usr/include/features.h /usr/include/stdc-predef.h
move/fragment.o: /usr/include/sys/cdefs.h /usr/include/bits/wordsize.h
move/fragment.o: /usr/include/gnu/stubs.h
strategy/strategy.o: main/config.h main/param_list.h main/stream_printf.h
strategy/strategy.o: strategy/strategy.h main/common.h
strategy/strategy.o: strategy/strategy_strict.h strategy/strategy_monte.h
strategy/strategy.o: strategy/strategy_boltz.h strategy/strategy_always.h
strategy/strategy_strict.o: main/common.h strategy/strategy_strict.h
strategy/strategy_strict.o: strategy/strategy.h main/param_list.h
strategy/strategy_monte.o: main/config.h main/param_list.h
strategy/strategy_monte.o: strategy/strategy_monte.h main/common.h
strategy/strategy_monte.o: strategy/strategy.h main/random.h
strategy/strategy_monte.o: main/stream_printf.h
strategy/strategy_boltz.o: main/config.h main/param_list.h main/random.h
strategy/strategy_boltz.o: main/stream_printf.h strategy/strategy_boltz.h
strategy/strategy_boltz.o: main/common.h strategy/strategy.h
strategy/strategy_always.o: main/common.h strategy/strategy_always.h
strategy/strategy_always.o: strategy/strategy.h main/param_list.h
score/scorer.o: score/scorer.h main/param_list.h score/scorer_combined.h
score/scorer.o: main/common.h main/config.h
score/scorer_combined.o: score/scorer_combined.h score/scorer.h
score/scorer_combined.o: main/param_list.h main/common.h peptide/atom.h
score/scorer_combined.o: peptide/amino.h peptide/atom_id.h
score/scorer_combined.o: peptide/atom_type.h peptide/pdb_atom_rec.h
score/scorer_combined.o: main/point.h peptide/peptide.h peptide/residue.h
score/scorer_combined.o: peptide/codon.h main/transform.h main/matrix.h
score/scorer_combined.o: peptide/conformation.h main/c_file.h
score/scorer_combined.o: main/stream_printf.h score/rapdf.h score/solvation.h
score/scorer_combined.o: score/torsion.h score/hbond.h score/saulo.h score/eleanor.h score/predss.h
score/scorer_combined.o: score/rgyr.h score/contact.h score/crowding.h score/predtor.h
score/scorer_combined.o: score/randomscr.h score/orientation.h score/core.h
score/scorer_combined.o: score/lennard_jones.h score/ribosome.h
score/rapdf.o: score/rapdf.h peptide/amino.h peptide/atom_id.h
score/rapdf.o: score/rapdf_impl.h peptide/peptide.h peptide/residue.h
score/rapdf.o: peptide/codon.h peptide/atom.h peptide/atom_type.h
score/rapdf.o: peptide/pdb_atom_rec.h main/point.h main/transform.h
score/rapdf.o: main/matrix.h peptide/conformation.h main/common.h
score/rapdf.o: score/scorer_combined.h score/scorer.h main/param_list.h
score/rapdf_impl.o: peptide/atom.h peptide/amino.h peptide/atom_id.h
score/rapdf_impl.o: peptide/atom_type.h peptide/pdb_atom_rec.h main/point.h
score/rapdf_impl.o: peptide/peptide.h peptide/residue.h peptide/codon.h
score/rapdf_impl.o: main/transform.h main/matrix.h peptide/conformation.h
score/rapdf_impl.o: main/common.h main/c_file.h score/scorer_combined.h
score/rapdf_impl.o: score/scorer.h main/param_list.h score/rapdf_impl.h
score/rapdf_impl.o: score/rapdf.h
score/solvation.o: score/solvation.h score/solvation_impl.h peptide/peptide.h
score/solvation.o: peptide/residue.h peptide/amino.h peptide/atom_id.h
score/solvation.o: peptide/codon.h peptide/atom.h peptide/atom_type.h
score/solvation.o: peptide/pdb_atom_rec.h main/point.h main/transform.h
score/solvation.o: main/matrix.h peptide/conformation.h main/common.h
score/solvation.o: score/scorer_combined.h score/scorer.h main/param_list.h
score/solvation_impl.o: peptide/peptide.h peptide/residue.h peptide/amino.h
score/solvation_impl.o: peptide/atom_id.h peptide/codon.h peptide/atom.h
score/solvation_impl.o: peptide/atom_type.h peptide/pdb_atom_rec.h
score/solvation_impl.o: main/point.h main/transform.h main/matrix.h
score/solvation_impl.o: peptide/conformation.h main/common.h main/c_file.h
score/solvation_impl.o: score/scorer_combined.h score/scorer.h
score/solvation_impl.o: main/param_list.h score/solvation.h
score/solvation_impl.o: score/solvation_impl.h
score/torsion.o: score/torsion.h score/torsion_impl.h peptide/amino.h
score/torsion.o: peptide/atom_id.h peptide/peptide.h peptide/residue.h
score/torsion.o: peptide/codon.h peptide/atom.h peptide/atom_type.h
score/torsion.o: peptide/pdb_atom_rec.h main/point.h main/transform.h
score/torsion.o: main/matrix.h peptide/conformation.h main/common.h
score/torsion.o: main/geom.h score/scorer_combined.h score/scorer.h
score/torsion.o: main/param_list.h
score/torsion_impl.o: peptide/peptide.h peptide/residue.h peptide/amino.h
score/torsion_impl.o: peptide/atom_id.h peptide/codon.h peptide/atom.h
score/torsion_impl.o: peptide/atom_type.h peptide/pdb_atom_rec.h main/point.h
score/torsion_impl.o: main/transform.h main/matrix.h peptide/conformation.h
score/torsion_impl.o: main/common.h main/geom.h main/c_file.h
score/torsion_impl.o: main/stream_printf.h score/scorer_combined.h
score/torsion_impl.o: score/scorer.h main/param_list.h score/torsion.h
score/torsion_impl.o: score/torsion_impl.h
score/hbond.o: peptide/peptide.h peptide/residue.h peptide/amino.h
score/hbond.o: peptide/atom_id.h peptide/codon.h peptide/atom.h
score/hbond.o: peptide/atom_type.h peptide/pdb_atom_rec.h main/point.h
score/hbond.o: main/transform.h main/matrix.h peptide/conformation.h
score/hbond.o: main/common.h score/scorer_combined.h score/scorer.h
score/hbond.o: main/param_list.h score/hbond.h
score/saulo.o: peptide/amino.h peptide/atom_id.h peptide/residue.h
score/saulo.o: peptide/codon.h peptide/atom.h peptide/atom_type.h
score/saulo.o: peptide/pdb_atom_rec.h main/point.h peptide/peptide.h
score/saulo.o: main/transform.h main/matrix.h peptide/conformation.h
score/saulo.o: main/common.h score/scorer_combined.h score/scorer.h
score/saulo.o: main/param_list.h score/saulo.h
score/eleanor.o: peptide/amino.h peptide/atom_id.h peptide/residue.h
score/eleanor.o: peptide/codon.h peptide/atom.h peptide/atom_type.h
score/eleanor.o: peptide/pdb_atom_rec.h main/point.h peptide/peptide.h
score/eleanor.o: main/transform.h main/matrix.h peptide/conformation.h
score/eleanor.o: main/common.h score/scorer_combined.h score/scorer.h
score/eleanor.o: main/param_list.h score/eleanor.h
score/predtor.o: peptide/amino.h peptide/atom_id.h peptide/residue.h
score/predtor.o: peptide/codon.h peptide/atom.h peptide/atom_type.h
score/predtor.o: peptide/pdb_atom_rec.h main/point.h peptide/peptide.h
score/predtor.o: main/transform.h main/matrix.h peptide/conformation.h
score/predtor.o: main/common.h score/scorer_combined.h score/scorer.h
score/predtor.o: main/param_list.h score/predtor.h
score/core.o: peptide/amino.h peptide/atom_id.h peptide/residue.h
score/core.o: peptide/codon.h peptide/atom.h peptide/atom_type.h
score/core.o: peptide/pdb_atom_rec.h main/point.h peptide/peptide.h
score/core.o: main/transform.h main/matrix.h peptide/conformation.h
score/core.o: main/common.h score/scorer_combined.h score/scorer.h
score/core.o: main/param_list.h score/core.h
score/core_impl.o: peptide/peptide.h peptide/residue.h peptide/amino.h
score/core_impl.o: peptide/atom_id.h peptide/codon.h peptide/atom.h
score/core_impl.o: peptide/atom_type.h peptide/pdb_atom_rec.h main/point.h
score/core_impl.o: main/transform.h main/matrix.h peptide/conformation.h
score/core_impl.o: main/common.h main/geom.h main/c_file.h
score/core_impl.o: main/stream_printf.h score/scorer_combined.h
score/core_impl.o: score/scorer.h main/param_list.h score/torsion.h
score/core_impl.o: score/torsion_impl.h
score/predss.o: peptide/amino.h peptide/atom_id.h peptide/residue.h
score/predss.o: peptide/codon.h peptide/atom.h peptide/atom_type.h
score/predss.o: peptide/pdb_atom_rec.h main/point.h peptide/peptide.h
score/predss.o: main/transform.h main/matrix.h peptide/conformation.h
score/predss.o: main/common.h score/scorer_combined.h score/scorer.h
score/predss.o: main/param_list.h score/predss.h
score/rgyr.o: peptide/peptide.h peptide/residue.h peptide/amino.h
score/rgyr.o: peptide/atom_id.h peptide/codon.h peptide/atom.h
score/rgyr.o: peptide/atom_type.h peptide/pdb_atom_rec.h main/point.h
score/rgyr.o: main/transform.h main/matrix.h peptide/conformation.h
score/rgyr.o: main/common.h score/scorer_combined.h score/scorer.h
score/rgyr.o: main/param_list.h score/rgyr.h
score/contact.o: peptide/peptide.h peptide/residue.h peptide/amino.h
score/contact.o: peptide/atom_id.h peptide/codon.h peptide/atom.h
score/contact.o: peptide/atom_type.h peptide/pdb_atom_rec.h main/point.h
score/contact.o: main/transform.h main/matrix.h peptide/conformation.h
score/contact.o: main/common.h score/scorer_combined.h score/scorer.h
score/contact.o: main/param_list.h score/contact.h
score/crowding.o: peptide/peptide.h peptide/residue.h peptide/amino.h
score/crowding.o: peptide/atom_id.h peptide/codon.h peptide/atom.h
score/crowding.o: peptide/atom_type.h peptide/pdb_atom_rec.h main/point.h
score/crowding.o: main/transform.h main/matrix.h peptide/conformation.h
score/crowding.o: main/common.h score/scorer_combined.h score/scorer.h
score/crowding.o: main/param_list.h score/crowding.h
score/randomscr.o: peptide/peptide.h peptide/residue.h peptide/amino.h
score/randomscr.o: peptide/atom_id.h peptide/codon.h peptide/atom.h
score/randomscr.o: peptide/atom_type.h peptide/pdb_atom_rec.h main/point.h
score/randomscr.o: main/transform.h main/matrix.h peptide/conformation.h
score/randomscr.o: main/common.h score/scorer_combined.h score/scorer.h
score/randomscr.o: main/param_list.h score/randomscr.h
score/orientation.o: score/orientation.h peptide/amino.h peptide/atom_id.h
score/orientation.o: score/orientation_impl.h peptide/peptide.h
score/orientation.o: peptide/residue.h peptide/codon.h peptide/atom.h
score/orientation.o: peptide/atom_type.h peptide/pdb_atom_rec.h main/point.h
score/orientation.o: main/transform.h main/matrix.h peptide/conformation.h
score/orientation.o: main/common.h main/geom.h score/scorer_combined.h
score/orientation.o: score/scorer.h main/param_list.h
score/orientation_impl.o: peptide/peptide.h peptide/residue.h peptide/amino.h
score/orientation_impl.o: peptide/atom_id.h peptide/codon.h peptide/atom.h
score/orientation_impl.o: peptide/atom_type.h peptide/pdb_atom_rec.h
score/orientation_impl.o: main/point.h main/transform.h main/matrix.h
score/orientation_impl.o: peptide/conformation.h main/common.h main/geom.h
score/orientation_impl.o: main/c_file.h score/scorer_combined.h
score/orientation_impl.o: score/scorer.h main/param_list.h
score/orientation_impl.o: score/orientation.h score/orientation_impl.h
score/lennard_jones.o: peptide/peptide.h peptide/residue.h peptide/amino.h
score/lennard_jones.o: peptide/atom_id.h peptide/codon.h peptide/atom.h
score/lennard_jones.o: peptide/atom_type.h peptide/pdb_atom_rec.h
score/lennard_jones.o: main/point.h main/transform.h main/matrix.h
score/lennard_jones.o: peptide/conformation.h main/common.h
score/lennard_jones.o: score/scorer_combined.h score/scorer.h
score/lennard_jones.o: main/param_list.h score/lennard_jones.h
score/ribosome.o: main/common.h peptide/peptide.h peptide/residue.h
score/ribosome.o: peptide/amino.h peptide/atom_id.h peptide/codon.h
score/ribosome.o: peptide/atom.h peptide/atom_type.h peptide/pdb_atom_rec.h
score/ribosome.o: main/point.h main/transform.h main/matrix.h
score/ribosome.o: peptide/conformation.h score/ribosome.h
peptide/peptide.o: /usr/include/sys/stat.h /usr/include/features.h
peptide/peptide.o: /usr/include/stdc-predef.h /usr/include/sys/cdefs.h
peptide/peptide.o: /usr/include/bits/wordsize.h /usr/include/gnu/stubs.h
peptide/peptide.o: /usr/include/bits/types.h /usr/include/bits/typesizes.h
peptide/peptide.o: /usr/include/time.h /usr/include/bits/stat.h
peptide/peptide.o: /usr/include/sys/types.h /usr/include/endian.h
peptide/peptide.o: /usr/include/bits/endian.h /usr/include/bits/byteswap.h
peptide/peptide.o: /usr/include/bits/byteswap-16.h /usr/include/sys/select.h
peptide/peptide.o: /usr/include/bits/select.h /usr/include/bits/sigset.h
peptide/peptide.o: /usr/include/bits/time.h /usr/include/sys/sysmacros.h
peptide/peptide.o: /usr/include/bits/pthreadtypes.h peptide/peptide.h
peptide/peptide.o: peptide/residue.h peptide/amino.h peptide/atom_id.h
peptide/peptide.o: peptide/codon.h peptide/atom.h peptide/atom_type.h
peptide/peptide.o: peptide/pdb_atom_rec.h main/point.h main/transform.h
peptide/peptide.o: main/matrix.h peptide/conformation.h main/common.h
peptide/peptide.o: peptide/sequence.h main/param_list.h main/c_file.h
peptide/peptide.o: main/config.h main/temp_file.h main/stream_printf.h
peptide/peptide.o: main/parse.h main/rmsd.h main/geom.h
peptide/residue.o: peptide/residue.h peptide/amino.h peptide/atom_id.h
peptide/residue.o: peptide/codon.h peptide/atom.h peptide/atom_type.h
peptide/residue.o: peptide/pdb_atom_rec.h main/point.h
peptide/atom.o: peptide/atom.h peptide/amino.h peptide/atom_id.h
peptide/atom.o: peptide/atom_type.h peptide/pdb_atom_rec.h main/point.h
peptide/sequence.o: main/common.h main/c_file.h peptide/peptide.h
peptide/sequence.o: peptide/residue.h peptide/amino.h peptide/atom_id.h
peptide/sequence.o: peptide/codon.h peptide/atom.h peptide/atom_type.h
peptide/sequence.o: peptide/pdb_atom_rec.h main/point.h main/transform.h
peptide/sequence.o: main/matrix.h peptide/conformation.h peptide/sequence.h
peptide/sequence.o: main/param_list.h main/config.h
peptide/amino.o: /usr/include/strings.h /usr/include/features.h
peptide/amino.o: /usr/include/stdc-predef.h /usr/include/sys/cdefs.h
peptide/amino.o: /usr/include/bits/wordsize.h /usr/include/gnu/stubs.h
peptide/amino.o: /usr/include/xlocale.h peptide/amino.h peptide/atom_id.h
peptide/amino.o: peptide/atom_type.h
peptide/codon.o: peptide/amino.h peptide/atom_id.h peptide/codon.h
peptide/atom_type.o: peptide/atom_type.h peptide/atom_id.h
peptide/pdb_atom_rec.o: /usr/include/stdio.h /usr/include/features.h
peptide/pdb_atom_rec.o: /usr/include/stdc-predef.h /usr/include/sys/cdefs.h
peptide/pdb_atom_rec.o: /usr/include/bits/wordsize.h /usr/include/gnu/stubs.h
peptide/pdb_atom_rec.o: /usr/include/bits/types.h
peptide/pdb_atom_rec.o: /usr/include/bits/typesizes.h /usr/include/libio.h
peptide/pdb_atom_rec.o: /usr/include/_G_config.h /usr/include/wchar.h
peptide/pdb_atom_rec.o: /usr/include/bits/stdio_lim.h
peptide/pdb_atom_rec.o: /usr/include/bits/sys_errlist.h
peptide/pdb_atom_rec.o: peptide/pdb_atom_rec.h main/point.h peptide/amino.h
peptide/pdb_atom_rec.o: peptide/atom_id.h peptide/atom_type.h
peptide/conformation.o: peptide/peptide.h peptide/residue.h peptide/amino.h
peptide/conformation.o: peptide/atom_id.h peptide/codon.h peptide/atom.h
peptide/conformation.o: peptide/atom_type.h peptide/pdb_atom_rec.h
peptide/conformation.o: main/point.h main/transform.h main/matrix.h
peptide/conformation.o: peptide/conformation.h main/common.h main/geom.h
extend/extender.o: extend/extender.h main/param_list.h
extend/extender.o: extend/extender_fixed.h extend/extender_codon.h
extend/extender.o: peptide/codon.h peptide/amino.h peptide/atom_id.h
extend/extender.o: main/common.h main/config.h main/runner.h
extend/extender.o: peptide/peptide.h peptide/residue.h peptide/atom.h
extend/extender.o: peptide/atom_type.h peptide/pdb_atom_rec.h main/point.h
extend/extender.o: main/transform.h main/matrix.h peptide/conformation.h
extend/extender_fixed.o: main/config.h main/param_list.h peptide/sequence.h
extend/extender_fixed.o: peptide/amino.h peptide/atom_id.h peptide/codon.h
extend/extender_fixed.o: main/c_file.h extend/extender.h
extend/extender_fixed.o: extend/extender_fixed.h
extend/extender_codon.o: main/config.h main/param_list.h main/common.h
extend/extender_codon.o: main/c_file.h main/stream_printf.h peptide/codon.h
extend/extender_codon.o: peptide/amino.h peptide/atom_id.h main/runner.h
extend/extender_codon.o: peptide/peptide.h peptide/residue.h peptide/atom.h
extend/extender_codon.o: peptide/atom_type.h peptide/pdb_atom_rec.h
extend/extender_codon.o: main/point.h main/transform.h main/matrix.h
extend/extender_codon.o: peptide/conformation.h peptide/sequence.h
extend/extender_codon.o: extend/extender.h extend/extender_codon.h
