CC = gcc -fPIC
CXX = g++
AR = ar
LD = g++
WINDRES = windres

INC = 
CFLAGS = -Wall -fexceptions -O3 -m64 -std=c++11 -w
RESINC = 
LIBDIR = 
LIB = 
LDFLAGS = 

INC_RELEASE = $(INC) -I.
CFLAGS_RELEASE = $(CFLAGS)
RESINC_RELEASE = $(RESINC)
RCFLAGS_RELEASE = $(RCFLAGS)
LIBDIR_RELEASE = $(LIBDIR)
LIB_RELEASE = $(LIB)
LDFLAGS_RELEASE = $(LDFLAGS) -s
OBJDIR_RELEASE = obj/Release
DEP_RELEASE = 
OUT_RELEASE = bin/Release/alifreefold

OBJ_RELEASE = $(OBJDIR_RELEASE)/ViennaRNA/model.o $(OBJDIR_RELEASE)/ViennaRNA/move_set.o $(OBJDIR_RELEASE)/ViennaRNA/multibranch_loops.o $(OBJDIR_RELEASE)/ViennaRNA/naview.o $(OBJDIR_RELEASE)/ViennaRNA/neighbor.o $(OBJDIR_RELEASE)/ViennaRNA/mm.o $(OBJDIR_RELEASE)/ViennaRNA/params.o $(OBJDIR_RELEASE)/ViennaRNA/part_func.o $(OBJDIR_RELEASE)/ViennaRNA/part_func_co.o $(OBJDIR_RELEASE)/ViennaRNA/part_func_up.o $(OBJDIR_RELEASE)/ViennaRNA/perturbation_fold.o $(OBJDIR_RELEASE)/ViennaRNA/hairpin_loops.o $(OBJDIR_RELEASE)/ViennaRNA/interior_loops.o $(OBJDIR_RELEASE)/ViennaRNA/inverse.o $(OBJDIR_RELEASE)/ViennaRNA/list.o $(OBJDIR_RELEASE)/ViennaRNA/mfe.o $(OBJDIR_RELEASE)/ViennaRNA/mfe_window.o $(OBJDIR_RELEASE)/ViennaRNA/ugly_bt.o $(OBJDIR_RELEASE)/ViennaRNA/stringdist.o $(OBJDIR_RELEASE)/ViennaRNA/structure_utils.o $(OBJDIR_RELEASE)/ViennaRNA/subopt.o $(OBJDIR_RELEASE)/ViennaRNA/treedist.o $(OBJDIR_RELEASE)/ViennaRNA/string_utils.o $(OBJDIR_RELEASE)/ViennaRNA/units.o $(OBJDIR_RELEASE)/ViennaRNA/unstructured_domains.o $(OBJDIR_RELEASE)/ViennaRNA/utils.o $(OBJDIR_RELEASE)/ViennaRNA/walk.o $(OBJDIR_RELEASE)/main.o $(OBJDIR_RELEASE)/ViennaRNA/read_epars.o $(OBJDIR_RELEASE)/ViennaRNA/plex.o $(OBJDIR_RELEASE)/ViennaRNA/plex_functions.o $(OBJDIR_RELEASE)/ViennaRNA/plot_aln.o $(OBJDIR_RELEASE)/ViennaRNA/plot_layouts.o $(OBJDIR_RELEASE)/ViennaRNA/plot_structure.o $(OBJDIR_RELEASE)/ViennaRNA/gquad.o $(OBJDIR_RELEASE)/ViennaRNA/ribo.o $(OBJDIR_RELEASE)/ViennaRNA/snofold.o $(OBJDIR_RELEASE)/ViennaRNA/snoop.o $(OBJDIR_RELEASE)/ViennaRNA/special_const.o $(OBJDIR_RELEASE)/ViennaRNA/RNAstruct.o $(OBJDIR_RELEASE)/ViennaRNA/ali_plex.o $(OBJDIR_RELEASE)/ViennaRNA/alifold.o $(OBJDIR_RELEASE)/ViennaRNA/alipfold.o $(OBJDIR_RELEASE)/ViennaRNA/aln_util.o $(OBJDIR_RELEASE)/ViennaRNA/alphabet.o $(OBJDIR_RELEASE)/ViennaRNA/ProfileDist.o $(OBJDIR_RELEASE)/ViennaRNA/boltzmann_sampling.o $(OBJDIR_RELEASE)/ViennaRNA/c_plex.o $(OBJDIR_RELEASE)/ViennaRNA/centroid.o $(OBJDIR_RELEASE)/ViennaRNA/cofold.o $(OBJDIR_RELEASE)/ViennaRNA/combinatorics.o $(OBJDIR_RELEASE)/RepWeightedMotif.o $(OBJDIR_RELEASE)/ViennaRNA/2Dfold.o $(OBJDIR_RELEASE)/ViennaRNA/2Dpfold.o $(OBJDIR_RELEASE)/ViennaRNA/LPfold.o $(OBJDIR_RELEASE)/ViennaRNA/commands.o $(OBJDIR_RELEASE)/ViennaRNA/Lfold.o $(OBJDIR_RELEASE)/ViennaRNA/MEA.o $(OBJDIR_RELEASE)/ViennaRNA/PS_dot.o $(OBJDIR_RELEASE)/ViennaRNA/ProfileAln.o $(OBJDIR_RELEASE)/ViennaRNA/file_formats.o $(OBJDIR_RELEASE)/ViennaRNA/energy_par.o $(OBJDIR_RELEASE)/ViennaRNA/equilibrium_probs.o $(OBJDIR_RELEASE)/ViennaRNA/eval.o $(OBJDIR_RELEASE)/ViennaRNA/exterior_loops.o $(OBJDIR_RELEASE)/ViennaRNA/file_formats_msa.o $(OBJDIR_RELEASE)/ViennaRNA/file_utils.o $(OBJDIR_RELEASE)/ViennaRNA/findpath.o $(OBJDIR_RELEASE)/ViennaRNA/fold.o $(OBJDIR_RELEASE)/ViennaRNA/constraints.o $(OBJDIR_RELEASE)/ViennaRNA/constraints_SHAPE.o $(OBJDIR_RELEASE)/ViennaRNA/constraints_hard.o $(OBJDIR_RELEASE)/ViennaRNA/constraints_ligand.o $(OBJDIR_RELEASE)/ViennaRNA/constraints_soft.o $(OBJDIR_RELEASE)/RepMotif.o $(OBJDIR_RELEASE)/ViennaRNA/convert_epars.o $(OBJDIR_RELEASE)/ViennaRNA/data_structures.o $(OBJDIR_RELEASE)/ViennaRNA/dist_vars.o $(OBJDIR_RELEASE)/ViennaRNA/dp_matrices.o $(OBJDIR_RELEASE)/ViennaRNA/duplex.o

all: release

clean: clean_release

before_release: 
	test -d bin/Release || mkdir -p bin/Release
	test -d $(OBJDIR_RELEASE)/ViennaRNA || mkdir -p $(OBJDIR_RELEASE)/ViennaRNA
	test -d $(OBJDIR_RELEASE) || mkdir -p $(OBJDIR_RELEASE)

after_release: 

release: before_release out_release after_release

out_release: before_release $(OBJ_RELEASE) $(DEP_RELEASE)
	$(LD) $(LIBDIR_RELEASE) -o $(OUT_RELEASE) $(OBJ_RELEASE)  $(LDFLAGS_RELEASE) $(LIB_RELEASE)

$(OBJDIR_RELEASE)/ViennaRNA/model.o: ViennaRNA/model.c
	$(CC) $(CFLAGS_RELEASE) $(INC_RELEASE) -c ViennaRNA/model.c -o $(OBJDIR_RELEASE)/ViennaRNA/model.o

$(OBJDIR_RELEASE)/ViennaRNA/move_set.o: ViennaRNA/move_set.c
	$(CC) $(CFLAGS_RELEASE) $(INC_RELEASE) -c ViennaRNA/move_set.c -o $(OBJDIR_RELEASE)/ViennaRNA/move_set.o

$(OBJDIR_RELEASE)/ViennaRNA/multibranch_loops.o: ViennaRNA/multibranch_loops.c
	$(CC) $(CFLAGS_RELEASE) $(INC_RELEASE) -c ViennaRNA/multibranch_loops.c -o $(OBJDIR_RELEASE)/ViennaRNA/multibranch_loops.o

$(OBJDIR_RELEASE)/ViennaRNA/naview.o: ViennaRNA/naview.c
	$(CC) $(CFLAGS_RELEASE) $(INC_RELEASE) -c ViennaRNA/naview.c -o $(OBJDIR_RELEASE)/ViennaRNA/naview.o

$(OBJDIR_RELEASE)/ViennaRNA/neighbor.o: ViennaRNA/neighbor.c
	$(CC) $(CFLAGS_RELEASE) $(INC_RELEASE) -c ViennaRNA/neighbor.c -o $(OBJDIR_RELEASE)/ViennaRNA/neighbor.o

$(OBJDIR_RELEASE)/ViennaRNA/mm.o: ViennaRNA/mm.c
	$(CC) $(CFLAGS_RELEASE) $(INC_RELEASE) -c ViennaRNA/mm.c -o $(OBJDIR_RELEASE)/ViennaRNA/mm.o

$(OBJDIR_RELEASE)/ViennaRNA/params.o: ViennaRNA/params.c
	$(CC) $(CFLAGS_RELEASE) $(INC_RELEASE) -c ViennaRNA/params.c -o $(OBJDIR_RELEASE)/ViennaRNA/params.o

$(OBJDIR_RELEASE)/ViennaRNA/part_func.o: ViennaRNA/part_func.c
	$(CC) $(CFLAGS_RELEASE) $(INC_RELEASE) -c ViennaRNA/part_func.c -o $(OBJDIR_RELEASE)/ViennaRNA/part_func.o

$(OBJDIR_RELEASE)/ViennaRNA/part_func_co.o: ViennaRNA/part_func_co.c
	$(CC) $(CFLAGS_RELEASE) $(INC_RELEASE) -c ViennaRNA/part_func_co.c -o $(OBJDIR_RELEASE)/ViennaRNA/part_func_co.o

$(OBJDIR_RELEASE)/ViennaRNA/part_func_up.o: ViennaRNA/part_func_up.c
	$(CC) $(CFLAGS_RELEASE) $(INC_RELEASE) -c ViennaRNA/part_func_up.c -o $(OBJDIR_RELEASE)/ViennaRNA/part_func_up.o

$(OBJDIR_RELEASE)/ViennaRNA/perturbation_fold.o: ViennaRNA/perturbation_fold.c
	$(CC) $(CFLAGS_RELEASE) $(INC_RELEASE) -c ViennaRNA/perturbation_fold.c -o $(OBJDIR_RELEASE)/ViennaRNA/perturbation_fold.o

$(OBJDIR_RELEASE)/ViennaRNA/hairpin_loops.o: ViennaRNA/hairpin_loops.c
	$(CC) $(CFLAGS_RELEASE) $(INC_RELEASE) -c ViennaRNA/hairpin_loops.c -o $(OBJDIR_RELEASE)/ViennaRNA/hairpin_loops.o

$(OBJDIR_RELEASE)/ViennaRNA/interior_loops.o: ViennaRNA/interior_loops.c
	$(CC) $(CFLAGS_RELEASE) $(INC_RELEASE) -c ViennaRNA/interior_loops.c -o $(OBJDIR_RELEASE)/ViennaRNA/interior_loops.o

$(OBJDIR_RELEASE)/ViennaRNA/inverse.o: ViennaRNA/inverse.c
	$(CC) $(CFLAGS_RELEASE) $(INC_RELEASE) -c ViennaRNA/inverse.c -o $(OBJDIR_RELEASE)/ViennaRNA/inverse.o

$(OBJDIR_RELEASE)/ViennaRNA/list.o: ViennaRNA/list.c
	$(CC) $(CFLAGS_RELEASE) $(INC_RELEASE) -c ViennaRNA/list.c -o $(OBJDIR_RELEASE)/ViennaRNA/list.o

$(OBJDIR_RELEASE)/ViennaRNA/mfe.o: ViennaRNA/mfe.c
	$(CC) $(CFLAGS_RELEASE) $(INC_RELEASE) -c ViennaRNA/mfe.c -o $(OBJDIR_RELEASE)/ViennaRNA/mfe.o

$(OBJDIR_RELEASE)/ViennaRNA/mfe_window.o: ViennaRNA/mfe_window.c
	$(CC) $(CFLAGS_RELEASE) $(INC_RELEASE) -c ViennaRNA/mfe_window.c -o $(OBJDIR_RELEASE)/ViennaRNA/mfe_window.o

$(OBJDIR_RELEASE)/ViennaRNA/ugly_bt.o: ViennaRNA/ugly_bt.c
	$(CC) $(CFLAGS_RELEASE) $(INC_RELEASE) -c ViennaRNA/ugly_bt.c -o $(OBJDIR_RELEASE)/ViennaRNA/ugly_bt.o

$(OBJDIR_RELEASE)/ViennaRNA/stringdist.o: ViennaRNA/stringdist.c
	$(CC) $(CFLAGS_RELEASE) $(INC_RELEASE) -c ViennaRNA/stringdist.c -o $(OBJDIR_RELEASE)/ViennaRNA/stringdist.o

$(OBJDIR_RELEASE)/ViennaRNA/structure_utils.o: ViennaRNA/structure_utils.c
	$(CC) $(CFLAGS_RELEASE) $(INC_RELEASE) -c ViennaRNA/structure_utils.c -o $(OBJDIR_RELEASE)/ViennaRNA/structure_utils.o

$(OBJDIR_RELEASE)/ViennaRNA/subopt.o: ViennaRNA/subopt.c
	$(CC) $(CFLAGS_RELEASE) $(INC_RELEASE) -c ViennaRNA/subopt.c -o $(OBJDIR_RELEASE)/ViennaRNA/subopt.o

$(OBJDIR_RELEASE)/ViennaRNA/treedist.o: ViennaRNA/treedist.c
	$(CC) $(CFLAGS_RELEASE) $(INC_RELEASE) -c ViennaRNA/treedist.c -o $(OBJDIR_RELEASE)/ViennaRNA/treedist.o

$(OBJDIR_RELEASE)/ViennaRNA/string_utils.o: ViennaRNA/string_utils.c
	$(CC) $(CFLAGS_RELEASE) $(INC_RELEASE) -c ViennaRNA/string_utils.c -o $(OBJDIR_RELEASE)/ViennaRNA/string_utils.o

$(OBJDIR_RELEASE)/ViennaRNA/units.o: ViennaRNA/units.c
	$(CC) $(CFLAGS_RELEASE) $(INC_RELEASE) -c ViennaRNA/units.c -o $(OBJDIR_RELEASE)/ViennaRNA/units.o

$(OBJDIR_RELEASE)/ViennaRNA/unstructured_domains.o: ViennaRNA/unstructured_domains.c
	$(CC) $(CFLAGS_RELEASE) $(INC_RELEASE) -c ViennaRNA/unstructured_domains.c -o $(OBJDIR_RELEASE)/ViennaRNA/unstructured_domains.o

$(OBJDIR_RELEASE)/ViennaRNA/utils.o: ViennaRNA/utils.c
	$(CC) $(CFLAGS_RELEASE) $(INC_RELEASE) -c ViennaRNA/utils.c -o $(OBJDIR_RELEASE)/ViennaRNA/utils.o

$(OBJDIR_RELEASE)/ViennaRNA/walk.o: ViennaRNA/walk.c
	$(CC) $(CFLAGS_RELEASE) $(INC_RELEASE) -c ViennaRNA/walk.c -o $(OBJDIR_RELEASE)/ViennaRNA/walk.o

$(OBJDIR_RELEASE)/main.o: main.cpp
	$(CXX) $(CFLAGS_RELEASE) $(INC_RELEASE) -c main.cpp -o $(OBJDIR_RELEASE)/main.o

$(OBJDIR_RELEASE)/ViennaRNA/read_epars.o: ViennaRNA/read_epars.c
	$(CC) $(CFLAGS_RELEASE) $(INC_RELEASE) -c ViennaRNA/read_epars.c -o $(OBJDIR_RELEASE)/ViennaRNA/read_epars.o

$(OBJDIR_RELEASE)/ViennaRNA/plex.o: ViennaRNA/plex.c
	$(CC) $(CFLAGS_RELEASE) $(INC_RELEASE) -c ViennaRNA/plex.c -o $(OBJDIR_RELEASE)/ViennaRNA/plex.o

$(OBJDIR_RELEASE)/ViennaRNA/plex_functions.o: ViennaRNA/plex_functions.c
	$(CC) $(CFLAGS_RELEASE) $(INC_RELEASE) -c ViennaRNA/plex_functions.c -o $(OBJDIR_RELEASE)/ViennaRNA/plex_functions.o

$(OBJDIR_RELEASE)/ViennaRNA/plot_aln.o: ViennaRNA/plot_aln.c
	$(CC) $(CFLAGS_RELEASE) $(INC_RELEASE) -c ViennaRNA/plot_aln.c -o $(OBJDIR_RELEASE)/ViennaRNA/plot_aln.o

$(OBJDIR_RELEASE)/ViennaRNA/plot_layouts.o: ViennaRNA/plot_layouts.c
	$(CC) $(CFLAGS_RELEASE) $(INC_RELEASE) -c ViennaRNA/plot_layouts.c -o $(OBJDIR_RELEASE)/ViennaRNA/plot_layouts.o

$(OBJDIR_RELEASE)/ViennaRNA/plot_structure.o: ViennaRNA/plot_structure.c
	$(CC) $(CFLAGS_RELEASE) $(INC_RELEASE) -c ViennaRNA/plot_structure.c -o $(OBJDIR_RELEASE)/ViennaRNA/plot_structure.o

$(OBJDIR_RELEASE)/ViennaRNA/gquad.o: ViennaRNA/gquad.c
	$(CC) $(CFLAGS_RELEASE) $(INC_RELEASE) -c ViennaRNA/gquad.c -o $(OBJDIR_RELEASE)/ViennaRNA/gquad.o

$(OBJDIR_RELEASE)/ViennaRNA/ribo.o: ViennaRNA/ribo.c
	$(CC) $(CFLAGS_RELEASE) $(INC_RELEASE) -c ViennaRNA/ribo.c -o $(OBJDIR_RELEASE)/ViennaRNA/ribo.o

$(OBJDIR_RELEASE)/ViennaRNA/snofold.o: ViennaRNA/snofold.c
	$(CC) $(CFLAGS_RELEASE) $(INC_RELEASE) -c ViennaRNA/snofold.c -o $(OBJDIR_RELEASE)/ViennaRNA/snofold.o

$(OBJDIR_RELEASE)/ViennaRNA/snoop.o: ViennaRNA/snoop.c
	$(CC) $(CFLAGS_RELEASE) $(INC_RELEASE) -c ViennaRNA/snoop.c -o $(OBJDIR_RELEASE)/ViennaRNA/snoop.o

$(OBJDIR_RELEASE)/ViennaRNA/special_const.o: ViennaRNA/special_const.c
	$(CC) $(CFLAGS_RELEASE) $(INC_RELEASE) -c ViennaRNA/special_const.c -o $(OBJDIR_RELEASE)/ViennaRNA/special_const.o

$(OBJDIR_RELEASE)/ViennaRNA/RNAstruct.o: ViennaRNA/RNAstruct.c
	$(CC) $(CFLAGS_RELEASE) $(INC_RELEASE) -c ViennaRNA/RNAstruct.c -o $(OBJDIR_RELEASE)/ViennaRNA/RNAstruct.o

$(OBJDIR_RELEASE)/ViennaRNA/ali_plex.o: ViennaRNA/ali_plex.c
	$(CC) $(CFLAGS_RELEASE) $(INC_RELEASE) -c ViennaRNA/ali_plex.c -o $(OBJDIR_RELEASE)/ViennaRNA/ali_plex.o

$(OBJDIR_RELEASE)/ViennaRNA/alifold.o: ViennaRNA/alifold.c
	$(CC) $(CFLAGS_RELEASE) $(INC_RELEASE) -c ViennaRNA/alifold.c -o $(OBJDIR_RELEASE)/ViennaRNA/alifold.o

$(OBJDIR_RELEASE)/ViennaRNA/alipfold.o: ViennaRNA/alipfold.c
	$(CC) $(CFLAGS_RELEASE) $(INC_RELEASE) -c ViennaRNA/alipfold.c -o $(OBJDIR_RELEASE)/ViennaRNA/alipfold.o

$(OBJDIR_RELEASE)/ViennaRNA/aln_util.o: ViennaRNA/aln_util.c
	$(CC) $(CFLAGS_RELEASE) $(INC_RELEASE) -c ViennaRNA/aln_util.c -o $(OBJDIR_RELEASE)/ViennaRNA/aln_util.o

$(OBJDIR_RELEASE)/ViennaRNA/alphabet.o: ViennaRNA/alphabet.c
	$(CC) $(CFLAGS_RELEASE) $(INC_RELEASE) -c ViennaRNA/alphabet.c -o $(OBJDIR_RELEASE)/ViennaRNA/alphabet.o

$(OBJDIR_RELEASE)/ViennaRNA/ProfileDist.o: ViennaRNA/ProfileDist.c
	$(CC) $(CFLAGS_RELEASE) $(INC_RELEASE) -c ViennaRNA/ProfileDist.c -o $(OBJDIR_RELEASE)/ViennaRNA/ProfileDist.o

$(OBJDIR_RELEASE)/ViennaRNA/boltzmann_sampling.o: ViennaRNA/boltzmann_sampling.c
	$(CC) $(CFLAGS_RELEASE) $(INC_RELEASE) -c ViennaRNA/boltzmann_sampling.c -o $(OBJDIR_RELEASE)/ViennaRNA/boltzmann_sampling.o

$(OBJDIR_RELEASE)/ViennaRNA/c_plex.o: ViennaRNA/c_plex.c
	$(CC) $(CFLAGS_RELEASE) $(INC_RELEASE) -c ViennaRNA/c_plex.c -o $(OBJDIR_RELEASE)/ViennaRNA/c_plex.o

$(OBJDIR_RELEASE)/ViennaRNA/centroid.o: ViennaRNA/centroid.c
	$(CC) $(CFLAGS_RELEASE) $(INC_RELEASE) -c ViennaRNA/centroid.c -o $(OBJDIR_RELEASE)/ViennaRNA/centroid.o

$(OBJDIR_RELEASE)/ViennaRNA/cofold.o: ViennaRNA/cofold.c
	$(CC) $(CFLAGS_RELEASE) $(INC_RELEASE) -c ViennaRNA/cofold.c -o $(OBJDIR_RELEASE)/ViennaRNA/cofold.o

$(OBJDIR_RELEASE)/ViennaRNA/combinatorics.o: ViennaRNA/combinatorics.c
	$(CC) $(CFLAGS_RELEASE) $(INC_RELEASE) -c ViennaRNA/combinatorics.c -o $(OBJDIR_RELEASE)/ViennaRNA/combinatorics.o

$(OBJDIR_RELEASE)/RepWeightedMotif.o: RepWeightedMotif.cpp
	$(CXX) $(CFLAGS_RELEASE) $(INC_RELEASE) -c RepWeightedMotif.cpp -o $(OBJDIR_RELEASE)/RepWeightedMotif.o

$(OBJDIR_RELEASE)/ViennaRNA/2Dfold.o: ViennaRNA/2Dfold.c
	$(CC) $(CFLAGS_RELEASE) $(INC_RELEASE) -c ViennaRNA/2Dfold.c -o $(OBJDIR_RELEASE)/ViennaRNA/2Dfold.o

$(OBJDIR_RELEASE)/ViennaRNA/2Dpfold.o: ViennaRNA/2Dpfold.c
	$(CC) $(CFLAGS_RELEASE) $(INC_RELEASE) -c ViennaRNA/2Dpfold.c -o $(OBJDIR_RELEASE)/ViennaRNA/2Dpfold.o

$(OBJDIR_RELEASE)/ViennaRNA/LPfold.o: ViennaRNA/LPfold.c
	$(CC) $(CFLAGS_RELEASE) $(INC_RELEASE) -c ViennaRNA/LPfold.c -o $(OBJDIR_RELEASE)/ViennaRNA/LPfold.o

$(OBJDIR_RELEASE)/ViennaRNA/commands.o: ViennaRNA/commands.c
	$(CC) $(CFLAGS_RELEASE) $(INC_RELEASE) -c ViennaRNA/commands.c -o $(OBJDIR_RELEASE)/ViennaRNA/commands.o

$(OBJDIR_RELEASE)/ViennaRNA/Lfold.o: ViennaRNA/Lfold.c
	$(CC) $(CFLAGS_RELEASE) $(INC_RELEASE) -c ViennaRNA/Lfold.c -o $(OBJDIR_RELEASE)/ViennaRNA/Lfold.o

$(OBJDIR_RELEASE)/ViennaRNA/MEA.o: ViennaRNA/MEA.c
	$(CC) $(CFLAGS_RELEASE) $(INC_RELEASE) -c ViennaRNA/MEA.c -o $(OBJDIR_RELEASE)/ViennaRNA/MEA.o

$(OBJDIR_RELEASE)/ViennaRNA/PS_dot.o: ViennaRNA/PS_dot.c
	$(CC) $(CFLAGS_RELEASE) $(INC_RELEASE) -c ViennaRNA/PS_dot.c -o $(OBJDIR_RELEASE)/ViennaRNA/PS_dot.o

$(OBJDIR_RELEASE)/ViennaRNA/ProfileAln.o: ViennaRNA/ProfileAln.c
	$(CC) $(CFLAGS_RELEASE) $(INC_RELEASE) -c ViennaRNA/ProfileAln.c -o $(OBJDIR_RELEASE)/ViennaRNA/ProfileAln.o

$(OBJDIR_RELEASE)/ViennaRNA/file_formats.o: ViennaRNA/file_formats.c
	$(CC) $(CFLAGS_RELEASE) $(INC_RELEASE) -c ViennaRNA/file_formats.c -o $(OBJDIR_RELEASE)/ViennaRNA/file_formats.o

$(OBJDIR_RELEASE)/ViennaRNA/energy_par.o: ViennaRNA/energy_par.c
	$(CC) $(CFLAGS_RELEASE) $(INC_RELEASE) -c ViennaRNA/energy_par.c -o $(OBJDIR_RELEASE)/ViennaRNA/energy_par.o

$(OBJDIR_RELEASE)/ViennaRNA/equilibrium_probs.o: ViennaRNA/equilibrium_probs.c
	$(CC) $(CFLAGS_RELEASE) $(INC_RELEASE) -c ViennaRNA/equilibrium_probs.c -o $(OBJDIR_RELEASE)/ViennaRNA/equilibrium_probs.o

$(OBJDIR_RELEASE)/ViennaRNA/eval.o: ViennaRNA/eval.c
	$(CC) $(CFLAGS_RELEASE) $(INC_RELEASE) -c ViennaRNA/eval.c -o $(OBJDIR_RELEASE)/ViennaRNA/eval.o

$(OBJDIR_RELEASE)/ViennaRNA/exterior_loops.o: ViennaRNA/exterior_loops.c
	$(CC) $(CFLAGS_RELEASE) $(INC_RELEASE) -c ViennaRNA/exterior_loops.c -o $(OBJDIR_RELEASE)/ViennaRNA/exterior_loops.o

$(OBJDIR_RELEASE)/ViennaRNA/file_formats_msa.o: ViennaRNA/file_formats_msa.c
	$(CC) $(CFLAGS_RELEASE) $(INC_RELEASE) -c ViennaRNA/file_formats_msa.c -o $(OBJDIR_RELEASE)/ViennaRNA/file_formats_msa.o

$(OBJDIR_RELEASE)/ViennaRNA/file_utils.o: ViennaRNA/file_utils.c
	$(CC) $(CFLAGS_RELEASE) $(INC_RELEASE) -c ViennaRNA/file_utils.c -o $(OBJDIR_RELEASE)/ViennaRNA/file_utils.o

$(OBJDIR_RELEASE)/ViennaRNA/findpath.o: ViennaRNA/findpath.c
	$(CC) $(CFLAGS_RELEASE) $(INC_RELEASE) -c ViennaRNA/findpath.c -o $(OBJDIR_RELEASE)/ViennaRNA/findpath.o

$(OBJDIR_RELEASE)/ViennaRNA/fold.o: ViennaRNA/fold.c
	$(CC) $(CFLAGS_RELEASE) $(INC_RELEASE) -c ViennaRNA/fold.c -o $(OBJDIR_RELEASE)/ViennaRNA/fold.o

$(OBJDIR_RELEASE)/ViennaRNA/constraints.o: ViennaRNA/constraints.c
	$(CC) $(CFLAGS_RELEASE) $(INC_RELEASE) -c ViennaRNA/constraints.c -o $(OBJDIR_RELEASE)/ViennaRNA/constraints.o

$(OBJDIR_RELEASE)/ViennaRNA/constraints_SHAPE.o: ViennaRNA/constraints_SHAPE.c
	$(CC) $(CFLAGS_RELEASE) $(INC_RELEASE) -c ViennaRNA/constraints_SHAPE.c -o $(OBJDIR_RELEASE)/ViennaRNA/constraints_SHAPE.o

$(OBJDIR_RELEASE)/ViennaRNA/constraints_hard.o: ViennaRNA/constraints_hard.c
	$(CC) $(CFLAGS_RELEASE) $(INC_RELEASE) -c ViennaRNA/constraints_hard.c -o $(OBJDIR_RELEASE)/ViennaRNA/constraints_hard.o

$(OBJDIR_RELEASE)/ViennaRNA/constraints_ligand.o: ViennaRNA/constraints_ligand.c
	$(CC) $(CFLAGS_RELEASE) $(INC_RELEASE) -c ViennaRNA/constraints_ligand.c -o $(OBJDIR_RELEASE)/ViennaRNA/constraints_ligand.o

$(OBJDIR_RELEASE)/ViennaRNA/constraints_soft.o: ViennaRNA/constraints_soft.c
	$(CC) $(CFLAGS_RELEASE) $(INC_RELEASE) -c ViennaRNA/constraints_soft.c -o $(OBJDIR_RELEASE)/ViennaRNA/constraints_soft.o

$(OBJDIR_RELEASE)/RepMotif.o: RepMotif.cpp
	$(CXX) $(CFLAGS_RELEASE) $(INC_RELEASE) -c RepMotif.cpp -o $(OBJDIR_RELEASE)/RepMotif.o

$(OBJDIR_RELEASE)/ViennaRNA/convert_epars.o: ViennaRNA/convert_epars.c
	$(CC) $(CFLAGS_RELEASE) $(INC_RELEASE) -c ViennaRNA/convert_epars.c -o $(OBJDIR_RELEASE)/ViennaRNA/convert_epars.o

$(OBJDIR_RELEASE)/ViennaRNA/data_structures.o: ViennaRNA/data_structures.c
	$(CC) $(CFLAGS_RELEASE) $(INC_RELEASE) -c ViennaRNA/data_structures.c -o $(OBJDIR_RELEASE)/ViennaRNA/data_structures.o

$(OBJDIR_RELEASE)/ViennaRNA/dist_vars.o: ViennaRNA/dist_vars.c
	$(CC) $(CFLAGS_RELEASE) $(INC_RELEASE) -c ViennaRNA/dist_vars.c -o $(OBJDIR_RELEASE)/ViennaRNA/dist_vars.o

$(OBJDIR_RELEASE)/ViennaRNA/dp_matrices.o: ViennaRNA/dp_matrices.c
	$(CC) $(CFLAGS_RELEASE) $(INC_RELEASE) -c ViennaRNA/dp_matrices.c -o $(OBJDIR_RELEASE)/ViennaRNA/dp_matrices.o

$(OBJDIR_RELEASE)/ViennaRNA/duplex.o: ViennaRNA/duplex.c
	$(CC) $(CFLAGS_RELEASE) $(INC_RELEASE) -c ViennaRNA/duplex.c -o $(OBJDIR_RELEASE)/ViennaRNA/duplex.o

clean_release: 
	rm -f $(OBJ_RELEASE) $(OUT_RELEASE)
	rm -rf bin/Release
	rm -rf $(OBJDIR_RELEASE)/ViennaRNA
	rm -rf $(OBJDIR_RELEASE)

.PHONY: before_debug after_debug clean_debug before_release after_release clean_release

