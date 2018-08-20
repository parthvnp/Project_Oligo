
################################################################################
# Contents:
#
#   Section 1: Object File Collections.
#              This section lists variables that contain groups of object files 
#              used in one or more RNAstructure programs.
#              Source file dependencies (*.h *.cpp etc) should NOT be listed here. They belong in Section 2.
#
#   Section 2: Source File Dependencies.
#              This section should list header-file dependencies for each object file.
#              At a later date, this section could be replaced with auto-generated dependencies.
#              Compilation recipes should NOT be listed here. They belong in Section 3.
#
#   Section 3: Special Compilation Rules.
#              Most object files can be compiled using the default compilation 
#              recipes (defined in compiler.h).
#              This section should only list compilation rules for object files 
#              that CANNOT be made using the default recipes.
################################################################################

# Note that the ROOTPATH variable must be defined in a Makefile before this file is included.
# Otherwise, the paths to the various dependencies may not resolve correctly.
R=${ROOTPATH}# The relative path from the client Makefile to the RNAstructure root directory.

################################################################################
#  Section 1: Object File Collections.
################################################################################
OUTDIR = $R/exe

# The text interface command line parser.
CMD_LINE_PARSER = \
	$R/src/ParseCommandLine.o \
	$R/src/common_utils.o

# The text interface configuration file parser.
CONF_PARSER = $R/src/configfile.o

# The dot plot handler.
PLOT_HANDLER = $R/src/DotPlotHandler.o

# The utility that creates a structure image comparer.
STRUCTURE_COMPARER = $R/src/StructureComparedImageHandler.o

# The utility that creates a structure image.
STRUCTURE_IMAGER = \
	$R/src/StructureImageHandler.o \
	$R/src/loop_utils.o

# The utility that scores a structure against another.
STRUCTURE_SCORER = $R/src/score.o

# Common files for the RNA library.
RNA_FILES = \
	$R/RNA_class/RNA.o \
	$R/RNA_class/thermodynamics.o \
	$R/RNA_class/design.o \
	$R/RNA_class/RsampleData.o \
	$R/src/algorithm.o \
	$R/src/alltrace.o \
	$R/src/DynProgArray.o \
	$R/src/dotarray.o \
	$R/src/draw.o \
	$R/src/extended_double.o \
	$R/src/forceclass.o \
	$R/src/MaxExpect.o \
	$R/src/MaxExpectStack.o \
	$R/src/outputconstraints.o \
	$R/src/pfunction.o \
	$R/src/probknot.o \
	$R/src/random.o \
	$R/src/rna_library.o \
	$R/src/stackclass.o \
	$R/src/stackstruct.o \
	$R/src/stochastic.o \
	$R/src/structure.o \
	$R/src/substructure.o \
	$R/src/TProgressDialog.o \
	$R/src/common_utils.o \
	${CUDA_RNA_FILES}

# Common files for the RNA library, for use with OPENMP/SMP.
# Simply filter-out the files that need to be replaced by SMP versions
RNA_FILES_SMP = $(filter-out \
	$R/src/algorithm.o \
	$R/src/pfunction.o \
	$R/src/stochastic.o \
	$R/src/MaxExpect.o \
	,$(RNA_FILES)) \
	$R/src/algorithm-smp.o \
	$R/src/MaxExpect-smp.o \
	$R/src/pfunction-smp.o \
	$R/src/stochastic-smp.o

# Files required for PARTS.
PARTS_FILES = \
	$R/PARTS/src/parts/alignment_priors.o \
	$R/PARTS/src/parts/array_file_manager.o \
	$R/PARTS/src/parts/array_mem_manager.o \
	$R/PARTS/src/parts/map_alignment.o \
	$R/PARTS/src/parts/map_mhr_info.o \
	$R/PARTS/src/parts/map_results.o \
	$R/PARTS/src/parts/map_structures.o \
	$R/PARTS/src/parts/pf_alignment.o \
	$R/PARTS/src/parts/pp_results.o \
	$R/PARTS/src/parts/ppf_cli.o \
	$R/PARTS/src/parts/ppf_loops.o \
	$R/PARTS/src/parts/ppf_operators.o \
	$R/PARTS/src/parts/ppf_progress_bar.o \
	$R/PARTS/src/parts/ppf_scale.o \
	$R/PARTS/src/parts/ppf_ss.o \
	$R/PARTS/src/parts/ppf_tb_stack.o \
	$R/PARTS/src/parts/ppf_timer.o \
	$R/PARTS/src/parts/ppf_v_mhe.o \
	$R/PARTS/src/parts/ppf_w_ext.o \
	$R/PARTS/src/parts/ppf_w.o \
	$R/PARTS/src/parts/ppf_w_l.o \
	$R/PARTS/src/parts/ppf_w_mb.o \
	$R/PARTS/src/parts/ppf_w_mbl.o \
	$R/PARTS/src/parts/ppf_w_mhi.o \
	$R/PARTS/src/parts/process_sequences.o \
	$R/PARTS/src/parts/single_pf_array.o \
	$R/PARTS/src/parts/template_pf_array.o \
	$R/PARTS/src/parts/stoch_tb/stoch_sampled_alignment.o \
	$R/PARTS/src/parts/stoch_tb/stoch_sampling_math.o \
	$R/PARTS/src/parts/stoch_tb/stoch_sampled_str_aln_sample_set.o \
	$R/PARTS/src/parts/stoch_tb/stoch_sampled_structures.o \
	$R/src/phmm/utils/rng/rng.o \
	$R/src/phmm/utils/rng/seed_manager.o

# Files for the PHMM calculation group used with both PARTS.
PHMM_FILES = \
	$R/src/phmm/aln_env_utils.o \
	$R/src/phmm/p_alignment.o \
	$R/src/phmm/phmm_aln.o \
	$R/src/phmm/phmm_array.o \
	$R/src/phmm/phmm.o \
	$R/src/phmm/phmm_ml_loops.o \
	$R/src/phmm/phmm_pp_loops.o \
	$R/src/phmm/structure/folding_constraints.o \
	$R/src/phmm/structure/structure_object.o \
	$R/src/phmm/utils/ansi_string/ansi_string.o \
	$R/src/phmm/utils/file/utils.o \
	$R/src/phmm/utils/xmath/linear/linear_math.o \
	$R/src/phmm/utils/xmath/log/xlog_math.o \
	$R/src/phmm/utils/xmath/matrix/matrix.o

# Files unique to ShapeKnots.
SHAPEKNOTS_FILES = \
	$R/ShapeKnots/ShapeKnots_Interface.o \
	$R/src/ShapeKnots.o \
	$R/src/pkHelix.o

####  TurboFold  ####
# Files common to both Turbofold and Turbofold-smp
TURBOFOLD_FILES = \
	$R/TurboFold/TurboFold_object.o \
	$R/TurboFold/Alignment.o \
	$R/TurboFold/MultiSequence.o \
	$R/TurboFold/Sequence.o \
	$R/TurboFold/ProbabilisticModel.o \
	$R/TurboFold/GuideTree.o \
	$R/TurboFold/SparseMatrix.o \
	${PHMM_FILES}

TURBOFOLD_FILES_SMP = \
	$(filter-out $R/TurboFold/TurboFold_object.o,${TURBOFOLD_FILES}) \
	$R/TurboFold/TurboFold_object-smp.o \
	$R/src/phmm/utils/ansi_thread/ansi_thread.o \
	$R/src/phmm/utils/mutex/ansi_mutex.o

####  CUDA  ####
CUDA_RNA_FILES=
# If the CUDA environment variable (or Make variable) is set, then
# include some CUDA files in with the list of RNA_FILES.
ifneq ($(CUDA),)
# Note: Fold-cuda now uses util.o from partition-smp 
# (but base.o is no longer compatible since partition-cuda was updated for the new datatables.)
  CUDA_RNA_FILES = \
	$R/fold-smp/frna.o \
	$R/fold-smp/fparam.o \
	$R/fold-smp/fbase.o \
	$R/partition-smp/prna.o \
	$R/partition-smp/base.o \
	$R/partition-smp/param.o \
	$R/partition-smp/util.o
endif

# Common files for ALL versions of the Dynalign library 
# (SERIAL and SMP versions of both Dynalign and Dynalign_II.)
# The filter-out section at the end includes files from RNA_FILES except those that 
# NOT needed for Dynalign_II (Dynalign does add back some of these.)
DYNALIGN_COMMON_FILES = \
	$R/RNA_class/TwoRNA.o \
	$R/src/dynalignarray.o \
	$R/src/dynalignheap.o \
	$R/src/observable.o \
	$R/src/observer.o \
	$R/src/varray.o \
	$R/src/wendarray.o \
	${PHMM_FILES} \
	${CONF_PARSER} \
	$(filter-out \
		$R/RNA_class/RNA.o \
		$R/RNA_class/design.o \
		$R/src/algorithm.o \
	, ${RNA_FILES})	

# Files for the Dynalign library when compiled serially.
DYNALIGN_SERIAL_FILES = \
	$R/RNA_class/Dynalign_object.o \
	$R/RNA_class/RNA.o \
	$R/src/algorithm.o \
	$R/src/dynalignstackclass.o \
	$R/src/dynalign.o \
	$R/src/bimol.o \
	${DYNALIGN_COMMON_FILES}

# Files for the SMP Dynalign library (NOT Dynalign_II)
# Note: the SMP version of Dynalign uses the SERIAL versions of RNA_FILES
DYNALIGN_SMP_FILES = \
	$R/RNA_class/Dynalign_object.o \
	$R/RNA_class/RNA.o \
	$R/src/algorithm.o \
	$R/src/dynalignstackclass.o \
	$R/src/dynalign-smp.o \
	$R/src/bimol-smp.o \
	${DYNALIGN_COMMON_FILES} \
	$R/src/rank.o \
	$R/src/rankconsumer.o \
	$R/src/rankmanager.o \
	$R/src/rankproducer.o \
	$R/src/observingtextprogressbar.o

# Files for the Dynalign_II library when compiled serially.
DYNALIGN_II_SERIAL_FILES = \
	$R/RNA_class/Dynalign_ii_object.o \
	$R/RNA_class/RNA_dynalign_ii.o \
	$R/src/algorithm_dynalign_ii.o \
	$R/src/dynalignstackclass_ii.o \
	$R/src/dynalign_ii.o \
	$R/src/bimol.o \
	${DYNALIGN_COMMON_FILES}

# Files for the Dynalign_II library when compiled for SMP.
DYNALIGN_II_SMP_FILES = \
	$R/RNA_class/Dynalign_ii_object.o \
	$R/RNA_class/RNA_dynalign_ii.o \
	$R/src/algorithm_dynalign_ii.o \
	$R/src/dynalignstackclass_ii.o \
	$R/src/dynalign_ii-smp.o \
	$R/src/bimol-smp.o \
	${DYNALIGN_COMMON_FILES} \
	$R/src/rank.o \
	$R/src/rankconsumer_ii.o \
	$R/src/rankmanager.o \
	$R/src/rankproducer.o \
	$R/src/observingtextprogressbar.o

# Common files for the SERIAL HybridRNA library.
HYBRID_FILES = \
	${RNA_FILES} \
	$R/RNA_class/HybridRNA.o \
	$R/RNA_class/TwoRNA.o \
	$R/src/bimol.o

# Common files for the HybridRNA library, for use with OPENMP SMP.
HYBRID_FILES_SMP = \
	${RNA_FILES_SMP} \
	$R/RNA_class/HybridRNA.o \
	$R/RNA_class/TwoRNA.o \
	$R/src/bimol-smp.o


# Common files for the Oligo library.
OLIGO_FILES = \
	${RNA_FILES} \
	$R/RNA_class/Oligowalk_object.o \
	$R/src/alltrace_intermolecular.o \
	$R/src/intermolecular.o \
	$R/src/OligoScreenCalc.o \
	$R/src/pclass.o \
	$R/src/siRNAfilter.o \
	$R/src/thermo.o

# Common files for the Oligo library.
OLIGO_FILES_SMP = \
	${RNA_FILES_SMP} \
	$R/RNA_class/Oligowalk_object.o \
	$R/src/alltrace_intermolecular.o \
	$R/src/intermolecular.o \
	$R/src/OligoScreenCalc-smp.o \
	$R/src/pclass.o \
	$R/src/siRNAfilter.o \
	$R/src/thermo.o

# Common files for the Pseudoknot library.
PSEUDOKNOT_FILES = \
	$R/src/basepair.o \
	$R/src/Pseudoknot.o \
	$R/src/PseudoParser.o

# Common files for ProbScan library
PROBSCAN_FILES = \
	$R/RNA_class/ProbScan.o

PROBSCAN_FILES_SMP = \
	${RNA_FILES_SMP} \
	$R/RNA_class/ProbScan.o

STRUCTURE_TOOLS_FILES = ### Future improvement, not yet implemented: $R/src/StructureTools.o $R/src/geometry2D.o

# Files for the Java RNAstructure interface (excluding the interface proxies).
# Updates to this list may be manually necessary if macros change or the Java drawing proxies change.
# Interface bindings used to make the native RNA library accessible from Java
# Updates to this list may be manually necessary if macros change or the Java drawing proxies change.
SWIG_OUT=$R/java_interface/SWIG
JAVA_SWIG_FILES = \
	$(SWIG_OUT)/drawing/DotPlotBackend.o \
	$(SWIG_OUT)/drawing/DotPlotBackend_wrap.o \
	$(SWIG_OUT)/drawing/StructureBackend.o \
	$(SWIG_OUT)/drawing/StructureBackend_wrap.o \
	$(SWIG_OUT)/RNAstructureBackendCalculator.o \
	$(SWIG_OUT)/RNAstructureBackendCalculator_wrap.o \
	$(SWIG_OUT)/RNA_wrap.o

 # Note that $(sort) removes duplicates.
JAVA_LIBRARY_FILES = $(sort \
	${RNA_FILES} \
	${PHMM_FILES} \
	${TURBOFOLD_FILES} \
	${DYNALIGN_SERIAL_FILES} \
	${HYBRID_FILES} \
	${OLIGO_FILES} \
	${PLOT_HANDLER} \
	${STRUCTURE_IMAGER} \
	${PROBSCAN_FILES} \
	$R/RNA_class/Multilign_object.o \
	${JAVA_SWIG_FILES} \
	${INCLUDE_GAMMA} \
	${STRUCTURE_TOOLS_FILES} \
)

# The SMP version of the Java library uses the SMP version of RNA_FILES,
# Which conflicts with DYNALIGN_SMP and TURBOFOLD_SMP because those use the SERIAL
# version of the RNA_FILES. 
# So the library uses the SERIAL versions of dynalign and turbofold.
JAVA_LIBRARY_FILES_SMP = $(sort \
	${RNA_FILES_SMP} \
	${PHMM_FILES} \
	${TURBOFOLD_FILES} \
	$(filter-out ${HYBRID_FILES},${DYNALIGN_SERIAL_FILES}) \
	${HYBRID_FILES_SMP} \
	${OLIGO_FILES_SMP} \
	${PLOT_HANDLER} \
	${STRUCTURE_IMAGER} \
	${PROBSCAN_FILES} \
	$R/RNA_class/Multilign_object.o \
	${JAVA_SWIG_FILES} \
	${INCLUDE_GAMMA} \
	${STRUCTURE_TOOLS_FILES} \
)

PYTHON_LIBRARY_FILES = $(filter-out $(JAVA_SWIG_FILES), $(JAVA_LIBRARY_FILES)) \
	$R/RNAstructure_python_interface/RNAstructure_wrap.o

CYCLEFOLD_FILES = \
	$R/CycleFold/NCM_parameters.o \
	$R/CycleFold/options.o \
	$R/CycleFold/arrays.o \
	$R/CycleFold/io.o \
	$R/CycleFold/sequence.o \
	$R/CycleFold/mainloop.o \
	$R/CycleFold/main.o \
	$R/CycleFold/logdouble.o \
	$R/CycleFold/constraints.o \
	$R/CycleFold/maxexpect.o \
	$R/CycleFold/extrinsic.o \
	$R/CycleFold/alignment.o \
	$R/CycleFold/turbo_calculation.o \
	${PHMM_FILES} \
	${CMD_LINE_PARSER}

################################################################################
#  Section 2: Source File Dependencies.
#
#  List additional header files that each object file depends on.
#
#  Not all files defined in dependency groups above need dependencies here.
#  Specifically if 'example.o' only depends on 'example.cpp' and 'example.h', 
#  then it does not require explicit dependencies (assuming the default 
#  compilation recipe is sufficient to compile it).
################################################################################

$R/TurboFold/TurboFold_object.o: \
	$R/TurboFold/Alignment.h \
	$R/TurboFold/MultiSequence.h \
	$R/TurboFold/Sequence.h \
	$R/TurboFold/SafeVector.h \
	$R/TurboFold/GuideTree.h \
	$R/TurboFold/ProbabilisticModel.h \
	$R/TurboFold/SparseMatrix.h

$R/TurboFold/Alignment.o: \
	$R/TurboFold/Alignment.h \
	$R/TurboFold/Alignment.cpp \

$R/TurboFold/TurboFold_Interface-smp.o: \
	$R/TurboFold/TurboFold_Interface.h \
	$R/TurboFold/TurboFold_Interface.cpp \
	$R/TurboFold/TurboFold_object.h \
	$R/src/configfile.h \
	$R/src/ErrorChecker.h \
	$R/src/ParseCommandLine.h \
	$R/src/TProgressDialog.h

$R/TurboFold/TurboFold_Interface.o: \
	$R/TurboFold/TurboFold_Interface.h \
	$R/TurboFold/TurboFold_Interface.cpp \
	$R/TurboFold/TurboFold_object.h \
	$R/src/configfile.h \
	$R/src/ErrorChecker.h \
	$R/src/ParseCommandLine.h \
	$R/src/TProgressDialog.h

$R/TurboFold/TurboFold_object-smp.o: \
	$R/TurboFold/TurboFold_object.h \
	$R/TurboFold/TurboFold_object.cpp

$R/AllSub/AllSub.o: \
	$R/AllSub/AllSub.cpp $R/AllSub/AllSub.h

$R/bifold/bifold.o: \
	$R/bifold/bifold.cpp $R/bifold/bifold.h

$R/CircleCompare/CircleCompare_Interface.o: \
	$R/CircleCompare/CircleCompare_Interface.cpp $R/CircleCompare/CircleCompare_Interface.h \
	$R/src/StructureComparedImageHandler.cpp $R/src/StructureComparedImageHandler.h \
	$R/src/StructureImageHandler.cpp $R/src/StructureImageHandler.h

$R/RNA_class/design.o: \
	$R/RNA_class/design.cpp $R/RNA_class/design.h \
	$R/RNA_class/RNA.h \
	$R/RNA_class/design.h \
	$R/RNA_class/thermodynamics.h \
	$R/src/ErrorChecker.h \
	$R/src/ParseCommandLine.h \
	$R/src/TProgressDialog.h \
	$R/src/algorithm.h \
	$R/src/DynProgArray.h \
	$R/src/defines.h \
	$R/src/dotarray.h \
	$R/src/draw.h \
	$R/src/forceclass.h \
	$R/src/pfunction.h \
	$R/src/random.h \
	$R/src/rna_library.h \
	$R/src/stackclass.h \
	$R/src/stackstruct.h \
	$R/src/structure.h \

$R/draw/DrawStructure.o: \
	$R/draw/DrawStructure.cpp $R/draw/DrawStructure.h \
	$R/src/StructureImageHandler.cpp $R/src/StructureImageHandler.h

$R/DuplexFold/DuplexFold.o: \
	$R/DuplexFold/DuplexFold.cpp $R/DuplexFold/DuplexFold.h

$R/DynalignDotPlot/DynalignDotPlot.o: \
	$R/DynalignDotPlot/DynalignDotPlot.cpp $R/DynalignDotPlot/DynalignDotPlot.h

$R/dynalign/dynaligninterface.o: $R/dynalign/dynaligninterface.cpp

$R/dynalign/dynaligninterface-smp.o: $R/dynalign/dynaligninterface.cpp
	@# Special Compilation Rule
	${COMPILE_SMP}  $<

$R/dynalign/dynaligninterface_ii.o: $R/dynalign/dynaligninterface.cpp
	@# Special Compilation Rule
	${COMPILE_DYNALIGN_II}  $<

$R/dynalign/dynaligninterface_ii-smp.o: $R/dynalign/dynaligninterface.cpp
	@# Special Compilation Rule
	${COMPILE_DYNALIGN_II_SMP}  $<

$R/efn2/efn2.o: \
	$R/efn2/efn2.cpp $R/efn2/efn2.h $R/src/PseudoParser.h

$R/efn2/efn2-smp.o: \
	$R/efn2/efn2.cpp $R/efn2/efn2.h $R/src/PseudoParser.h
	@# Special Compilation Rule
	${COMPILE_OMP}  $<

$R/EnergyPlot/EnergyPlot.o: \
	$R/EnergyPlot/EnergyPlot.cpp $R/EnergyPlot/EnergyPlot.h

$R/EnsembleEnergy/EnsembleEnergy_Interface.o: \
	$R/EnsembleEnergy/EnsembleEnergy_Interface.cpp $R/EnsembleEnergy/EnsembleEnergy_Interface.h

$R/fold/Fold.o: \
	$R/fold/Fold.cpp $R/fold/Fold.h $R/src/ErrorChecker.h

$R/src/loop_utils.o: \
	$R/src/loop_utils.cpp $R/src/loop_utils.h

$R/MaxExpect/MaxExpectInterface.o: \
	$R/MaxExpect/MaxExpectInterface.cpp $R/MaxExpect/MaxExpect.h

$R/multilign/Multilign_Interface.o: \
	$R/multilign/Multilign_Interface.cpp $R/multilign/Multilign_Interface.h

$R/multilign/Multilign_Interface-smp.o: \
	$R/multilign/Multilign_Interface.cpp $R/multilign/Multilign_Interface.h
	@# Special Compilation Rule
	${COMPILE_SMP}    $<

$R/Multifind/Multifind_Interface.o: \
	$R/Multifind/Multifind_Interface.cpp $R/Multifind/Multifind_Interface.h
	@# Special Compilation Rule
	${COMPILE_SVM}    $<

$R/Multifind/Multifind_Interface-smp.o: \
	$R/Multifind/Multifind_Interface.cpp $R/Multifind/Multifind_Interface.h
	@# Special Compilation Rule
	${COMPILE_SVM_SMP}    $<

$R/napss/NAPSS_Interface.o: \
    $R/napss/NAPSS_Interface.cpp $R/napss/NAPSS_Interface.h \
    $R/napss/napss.h

$R/napss/napss.o: \
	$R/napss/napss.cpp \
	$R/src/algorithm.h \
	$R/src/DynProgArray.h \
	$R/src/configfile.h \
	$R/src/defines.h \
	$R/src/dotarray.h \
	$R/src/forceclass.h \
	$R/src/platform.h \
	$R/src/rna_library.h \
	$R/src/stackclass.h \
	$R/src/stackstruct.h \
	$R/src/structure.h \
	$R/src/TProgressDialog.h

$R/oligoscreen/oligoscreen.o: \
	$R/oligoscreen/oligoscreen.cpp $R/oligoscreen/oligoscreen.h

$R/oligowalk/src/globals.o: \
	$R/oligowalk/src/globals.cpp $R/oligowalk/src/globals.h

$R/oligowalk/src/oligowalk.o: \
	$R/oligowalk/src/globals.h \
	$R/oligowalk/src/oligowalk.cpp \
	$R/src/algorithm.h \
	$R/src/alltrace.h \
	$R/src/defines.h \
	$R/src/intermolecular.h \
	$R/src/pclass.h \
	$R/src/platform.h \
	$R/src/rna_library.h \
	$R/src/siRNAfilter.h \
	$R/src/structure.h \
	$R/src/stochastic.h \
	$R/src/thermo.h \
	$R/src/TProgressDialog.h

$R/pfunction/partition.o: \
	$R/pfunction/partition.cpp $R/pfunction/partition.h

$R/ProbablePair/ProbablePair.o: \
	$R/ProbablePair/ProbablePair.cpp $R/ProbablePair/ProbablePair.h

$R/ProbabilityPlot/ProbabilityPlot.o: \
	$R/ProbabilityPlot/ProbabilityPlot.cpp $R/ProbabilityPlot/ProbabilityPlot.h

$R/ProbKnot/ProbKnot_Interface.o: \
	$R/ProbKnot/ProbKnot_Interface.cpp $R/ProbKnot/ProbKnot_Interface.h

$R/ProbKnot/ProbScan_Interface.o: \
	$R/ProbScan/ProbScan_Interface.cpp $R/ProbKnot/ProbScan_Interface.h

$R/refold/refold.o: \
	$R/refold/refold.cpp $R/refold/refold.h

$R/RemovePseudoknots/RemovePseudoknots.o: \
	$R/RemovePseudoknots/RemovePseudoknots.cpp $R/RemovePseudoknots/RemovePseudoknots.h

$R/RNA_class/Dynalign_class.o: RNA_class/Dynalign_class.cpp

$R/RNA_class/Dynalign_object.o: \
	$R/RNA_class/Dynalign_object.cpp $R/RNA_class/Dynalign_object.h \
	$R/src/platform.h

$R/RNA_class/Dynalign_ii_object.o: \
	$R/RNA_class/Dynalign_object.cpp $R/RNA_class/Dynalign_object.h \
	$R/src/platform.h
	@# Special Compilation Rule
	${COMPILE_DYNALIGN_II}    $<

$R/RNA_class/HybridRNA.o: \
	$R/RNA_class/HybridRNA.cpp $R/RNA_class/HybridRNA.h \
	$R/RNA_class/RNA.cpp $R/RNA_class/RNA.h \
	$R/RNA_class/thermodynamics.cpp $R/RNA_class/thermodynamics.h \
	$R/RNA_class/TwoRNA.cpp $R/RNA_class/TwoRNA.h \
	$R/src/algorithm.cpp $R/src/algorithm.h \
	$R/src/DynProgArray.h \
	$R/src/bimol.h \
	$R/src/defines.h \
	$R/src/dotarray.h \
	$R/src/forceclass.h \
	$R/src/pfunction.h \
	$R/src/platform.h \
	$R/src/rna_library.h \
	$R/src/stackclass.h \
	$R/src/stackstruct.h \
	$R/src/structure.h \
	$R/src/TProgressDialog.h

$R/RNA_class/HybridRNA_class.o: \
	$R/RNA_class/HybridRNA_class.cpp \
	$R/src/structure.h

$R/RNA_class/Multifind_object.o: \
	$R/RNA_class/Multifind_object.cpp $R/RNA_class/Multifind_object.h
	@# Special Compilation Rule
	${COMPILE_SVM_SMP}    $<

$R/RNA_class/Multilign_object.o: \
	$R/RNA_class/Dynalign_object.h \
	$R/RNA_class/Multilign_object.cpp $R/RNA_class/Multilign_object.h \
	$R/RNA_class/RNA.h \
	$R/RNA_class/thermodynamics.h \
	$R/RNA_class/TwoRNA.h \
	$R/src/algorithm.h \
	$R/src/DynProgArray.h \
	$R/src/defines.h \
	$R/src/dotarray.h \
	$R/src/draw.h \
	$R/src/dynalign.h \
	$R/src/dynalignarray.h \
	$R/src/dynalignheap.h \
	$R/src/forceclass.h \
	$R/src/pfunction.h \
	$R/src/platform.h \
	$R/src/random.h \
	$R/src/rna_library.h \
	$R/src/stackclass.h \
	$R/src/stackstruct.h \
	$R/src/structure.h \
	$R/src/TProgressDialog.h \
	$R/src/varray.h \
	$R/src/wendarray.h

#### addSuffixedObject ####
#$(call addSuffixedObject Multilign_object,-Multifind,-D MULTIFIND)

$R/RNA_class/Multilign_object-Multifind.o: \
	$R/RNA_class/Multilign_object.cpp $R/RNA_class/Multilign_object.h \
	$R/RNA_class/Dynalign_object.h \
	$R/RNA_class/RNA.h \
	$R/RNA_class/thermodynamics.h \
	$R/RNA_class/TwoRNA.h \
	$R/src/algorithm.h \
	$R/src/DynProgArray.h \
	$R/src/defines.h \
	$R/src/dotarray.h \
	$R/src/draw.h \
	$R/src/dynalign.h \
	$R/src/dynalignarray.h \
	$R/src/dynalignheap.h \
	$R/src/forceclass.h \
	$R/src/pfunction.h \
	$R/src/platform.h \
	$R/src/random.h \
	$R/src/rna_library.h \
	$R/src/stackclass.h \
	$R/src/stackstruct.h \
	$R/src/structure.h \
	$R/src/TProgressDialog.h \
	$R/src/varray.h \
	$R/src/wendarray.h
	@# Special Compilation Rule
	${COMPILE_MULTIFIND}    $<


$R/phmm/phmm_interface.o: $R/phmm/phmm_interface.cpp $R/phmm/phmm_interface.h

$R/RNA_class/OligoWalk_class.o: $R/RNA_class/OligoWalk_class.cpp

$R/RNA_class/Oligowalk_object.o: $R/RNA_class/Oligowalk_object.cpp $R/RNA_class/Oligowalk_object.h

$R/RNA_class/RNA.o: \
	$R/RNA_class/RNA.cpp $R/RNA_class/RNA.h \
	$R/RNA_class/thermodynamics.cpp $R/RNA_class/thermodynamics.h \
	$R/src/algorithm.cpp $R/src/algorithm.h \
	$R/src/alltrace.h \
	$R/src/DynProgArray.h \
	$R/src/bimol.h \
	$R/src/defines.h \
	$R/src/dotarray.h \
	$R/src/draw.h \
	$R/src/forceclass.h \
	$R/src/MaxExpect.h \
	$R/src/pfunction.h \
	$R/src/probknot.h \
	$R/src/platform.h \
	$R/src/random.h \
	$R/src/rna_library.h \
	$R/src/stackclass.h \
	$R/src/stackstruct.h \
	$R/src/stochastic.h \
	$R/src/structure.h \
	$R/src/TProgressDialog.h

$R/RNA_class/RNA_dynalign_ii.o: \
	$R/RNA_class/RNA.cpp $R/RNA_class/RNA.h \
	$R/RNA_class/thermodynamics.cpp $R/RNA_class/thermodynamics.h \
	$R/src/algorithm.cpp $R/src/algorithm.h \
	$R/src/alltrace.h \
	$R/src/DynProgArray.h \
	$R/src/bimol.h \
	$R/src/defines.h \
	$R/src/dotarray.h \
	$R/src/draw.h \
	$R/src/forceclass.h \
	$R/src/MaxExpect.h \
	$R/src/pfunction.h \
	$R/src/probknot.h \
	$R/src/platform.h \
	$R/src/random.h \
	$R/src/rna_library.h \
	$R/src/stackclass.h \
	$R/src/stackstruct.h \
	$R/src/stochastic.h \
	$R/src/structure.h \
	$R/src/TProgressDialog.h
	@# Special Compilation Rule
	${COMPILE_DYNALIGN_II}    $<

$R/RNA_class/RNA_class.o: \
	$R/RNA_class/RNA_class.cpp \
	$R/src/structure.h

$R/RNA_class/thermodynamics.o: \
	$R/RNA_class/thermodynamics.cpp $R/RNA_class/thermodynamics.h

$R/RNA_class/TwoRNA.o: \
	$R/RNA_class/RNA.cpp $R/RNA_class/RNA.h \
	$R/RNA_class/thermodynamics.cpp $R/RNA_class/thermodynamics.h \
	$R/RNA_class/TwoRNA.cpp $R/RNA_class/TwoRNA.h \
	$R/src/algorithm.cpp $R/src/algorithm.h \
	$R/src/DynProgArray.h \
	$R/src/bimol.h \
	$R/src/defines.h \
	$R/src/dotarray.h \
	$R/src/forceclass.h \
	$R/src/pfunction.h \
	$R/src/platform.h \
	$R/src/rna_library.h \
	$R/src/stackclass.h \
	$R/src/stackstruct.h \
	$R/src/structure.h \
	$R/src/TProgressDialog.h

$R/scorer/Scorer_Interface.o: \
	$R/scorer/Scorer_Interface.cpp $R/scorer/Scorer_Interface.h

$R/src/algorithm.o: \
	$R/src/algorithm.cpp $R/src/algorithm.h \
	$R/src/DynProgArray.h \
	$R/src/defines.h \
	$R/src/dotarray.h \
	$R/src/forceclass.h \
	$R/src/platform.h \
	$R/src/rna_library.h \
	$R/src/stackclass.h \
	$R/src/stackstruct.h \
	$R/src/structure.h \
	$R/src/TProgressDialog.h

$R/src/algorithm-smp.o: \
	$R/src/algorithm.cpp $R/src/algorithm.h \
	$R/src/DynProgArray.h \
	$R/src/defines.h \
	$R/src/dotarray.h \
	$R/src/forceclass.h \
	$R/src/platform.h \
	$R/src/rna_library.h \
	$R/src/stackclass.h \
	$R/src/stackstruct.h \
	$R/src/structure.h \
	$R/src/TProgressDialog.h
	@# Special Compilation Rule
	${COMPILE_OMP}    $<


$R/src/algorithm_dynalign_ii.o: \
	$R/src/algorithm.cpp $R/src/algorithm.h \
	$R/src/DynProgArray.h \
	$R/src/defines.h \
	$R/src/dotarray.h \
	$R/src/forceclass.h \
	$R/src/platform.h \
	$R/src/rna_library.h \
	$R/src/stackclass.h \
	$R/src/stackstruct.h \
	$R/src/structure.h \
	$R/src/TProgressDialog.h
	@# Special Compilation Rule
	${COMPILE_DYNALIGN_II}    $<

$R/src/algorithm_instrumented.o: \
	$R/src/algorithm.cpp $R/src/algorithm.h \
	$R/src/DynProgArray.h \
	$R/src/defines.h \
	$R/src/dotarray.h \
	$R/src/forceclass.h \
	$R/src/platform.h \
	$R/src/rna_library.h \
	$R/src/stackclass.h \
	$R/src/stackstruct.h \
	$R/src/structure.h \
	$R/src/TProgressDialog.h
	@# Special Compilation Rule
	${COMPILE_INSTRUMENTED}    $<

$R/src/alltrace.o: \
	$R/src/alltrace.cpp $R/src/alltrace.h \
	$R/src/defines.h \
	$R/src/structure.h

$R/src/alltrace_intermolecular.o: \
	$R/src/alltrace_intermolecular.cpp $R/src/alltrace_intermolecular.h

$R/src/DynProgArray.o: \
	$R/src/DynProgArray.cpp $R/src/DynProgArray.h \
	$R/src/defines.h

$R/src/bimol.o: \
	$R/src/bimol.cpp $R/src/bimol.h

$R/src/bimol-smp.o: \
	$R/src/bimol.cpp $R/src/bimol.h
	@# Special Compilation Rule
	${COMPILE_OMP}    $<

$R/src/configfile.o: \
	$R/src/configfile.cpp $R/src/configfile.h

$R/src/dotarray.o: \
	$R/src/defines.h \
	$R/src/dotarray.cpp $R/src/dotarray.h

$R/src/DotPlotHandler.o: \
	$R/src/DotPlotHandler.cpp $R/src/DotPlotHandler.h

$R/src/draw.o: \
	$R/src/draw.cpp $R/src/draw.h \
	$R/src/ErrorChecker.h \
	$R/src/substructure.h

$R/src/dynalign.o: \
	$R/src/algorithm.cpp

$R/src/dynalign.o: \
	$R/src/algorithm.h \
	$R/src/DynProgArray.h \
	$R/src/defines.h \
	$R/src/dynalign.cpp $R/src/dynalign.h \
	$R/src/dynalignarray.h \
	$R/src/dynalignheap.h \
	$R/src/dynalignstackclass.h \
	$R/src/forceclass.h \
	$R/src/platform.h \
	$R/src/rna_library.h \
	$R/src/structure.h \
	$R/src/varray.h \
	$R/src/wendarray.h \
	$R/src/TProgressDialog.h

$R/src/dynalign-smp.o: \
	$R/src/dynalign.cpp $R/src/dynalign.h \
	$R/src/algorithm.h \
	$R/src/DynProgArray.h \
	$R/src/defines.h \
	$R/src/dynalignarray.h \
	$R/src/dynalignheap.h \
	$R/src/dynalignstackclass.h \
	$R/src/forceclass.h \
	$R/src/observingtextprogressbar.h \
	$R/src/platform.h \
	$R/src/rankconsumer.h \
	$R/src/rankmanager.h \
	$R/src/rankproducer.h \
	$R/src/rna_library.h \
	$R/src/structure.h \
	$R/src/TProgressDialog.h \
	$R/src/varray.h \
	$R/src/wendarray.h
	@# Special Compilation Rule
	${COMPILE_SMP}    $<

$R/src/dynalign_ii.o: \
	$R/src/dynalign.cpp $R/src/dynalign.h \
	$R/src/algorithm.h \
	$R/src/DynProgArray.h \
	$R/src/defines.h \
	$R/src/dynalignarray.h \
	$R/src/dynalignheap.h \
	$R/src/dynalignstackclass.h \
	$R/src/forceclass.h \
	$R/src/platform.h \
	$R/src/rna_library.h \
	$R/src/structure.h \
	$R/src/varray.h \
	$R/src/wendarray.h \
	$R/src/TProgressDialog.h
	@# Special Compilation Rule
	${COMPILE_DYNALIGN_II}    $<

$R/src/dynalign_ii-smp.o: \
	$R/src/dynalign.cpp $R/src/dynalign.h \
	$R/src/algorithm.h \
	$R/src/DynProgArray.h \
	$R/src/defines.h \
	$R/src/dynalignarray.h \
	$R/src/dynalignheap.h \
	$R/src/dynalignstackclass.h \
	$R/src/forceclass.h \
	$R/src/platform.h \
	$R/src/rna_library.h \
	$R/src/structure.h \
	$R/src/varray.h \
	$R/src/wendarray.h \
	$R/src/TProgressDialog.h
	@# Special Compilation Rule
	${COMPILE_DYNALIGN_II_SMP}    $<

$R/src/dynalignarray.o: \
	$R/src/dynalignarray.cpp $R/src/dynalignarray.h \
	$R/src/defines.h \
	$R/src/dynalign.h

$R/src/dynalignheap.o: \
	$R/src/dynalignheap.cpp $R/src/dynalignheap.h \
	$R/src/defines.h

$R/src/dynalignstackclass.o: \
	$R/src/dynalignstackclass.cpp $R/src/dynalignstackclass.h \
	$R/src/defines.h

$R/src/dynalignstackclass_ii.o: \
	$R/src/dynalignstackclass.cpp $R/src/dynalignstackclass.h \
	$R/src/defines.h
	@# Special Compilation Rule
	${COMPILE_DYNALIGN_II}    $<

$R/src/extended_double.o: \
	$R/src/extended_double.cpp $R/src/extended_double.h

$R/src/forceclass.o: \
	$R/src/forceclass.cpp $R/src/forceclass.h

$R/src/intermolecular.o: \
	$R/src/intermolecular.cpp $R/src/intermolecular.h \
	$R/src/siRNAfilter.cpp $R/src/siRNAfilter.h

$R/src/MaxExpect.o: \
	$R/src/MaxExpect.cpp $R/src/MaxExpect.h \
	$R/src/defines.h

$R/src/MaxExpect-smp.o: \
	$R/src/MaxExpect.cpp $R/src/MaxExpect.h
	@# Special Compilation Rule
	${COMPILE_OMP}    $<

$R/src/MaxExpectStack.o: \
	$R/src/MaxExpectStack.cpp $R/src/MaxExpectStack.h

$R/src/observable.o: \
	$R/src/observable.cpp $R/src/observable.h \
	$R/src/observer.h

$R/src/observer.o: \
	$R/src/observer.cpp $R/src/observer.h

$R/src/observingtextprogressbar.o: \
	$R/src/observingtextprogressbar.cpp $R/src/observingtextprogressbar.h \
	$R/src/TProgressDialog.h

$R/src/OligoScreenCalc.o: \
	$R/src/OligoScreenCalc.cpp $R/src/OligoScreenCalc.h

$R/src/OligoScreenCalc-smp.o: \
	$R/src/OligoScreenCalc.cpp $R/src/OligoScreenCalc.h
	@# Special Compilation Rule
	${COMPILE_OMP}    $<


$R/src/outputconstraints.o: \
	$R/src/outputconstraints.cpp $R/src/outputconstraints.h

$R/src/ParseCommandLine.o: \
	$R/src/ParseCommandLine.cpp $R/src/ParseCommandLine.h \
	$R/src/ErrorChecker.h \
	$R/src/version.h

$R/src/pclass.o: \
	$R/src/pclass.cpp $R/src/pclass.h

$R/src/probknot.o: \
	$R/src/probknot.cpp $R/src/probknot.h

$R/src/pfunction.o: \
	$R/src/pfunction.cpp $R/src/pfunction.h $R/src/boltzmann.h \
	$R/src/algorithm.h $R/src/structure.h $R/src/DynProgArray.h

$R/src/pfunction-smp.o: \
	$R/src/pfunction.cpp $R/src/pfunction.h $R/src/boltzmann.h \
	$R/src/algorithm.h $R/src/structure.h  
	@# Special Compilation Rule
	${COMPILE_OMP}    $<

$R/src/phmm.o: \
	$R/src/phmm.cpp $R/src/phmm.h

$R/src/random.o: \
	$R/src/random.cpp $R/src/random.h

$R/src/rank.o: \
	$R/src/rank.cpp $R/src/rank.h \
	$R/src/workslice.h \
	$R/src/workunit.h

$R/src/rankconsumer.o: \
	$R/src/dynalign.h \
	$R/src/dynalignarray.h \
	$R/src/rankconsumer.cpp $R/src/rankconsumer.h \
	$R/src/rankmanager.h \
	$R/src/rna_library.h \
	$R/src/structure.h \
	$R/src/varray.h \
	$R/src/wendarray.h \
	$R/src/workslice.h

$R/src/rankconsumer_ii.o: \
	$R/src/rankconsumer.cpp $R/src/rankconsumer.h \
	$R/src/dynalign.h \
	$R/src/dynalignarray.h \
	$R/src/rankmanager.h \
	$R/src/rna_library.h \
	$R/src/structure.h \
	$R/src/varray.h \
	$R/src/wendarray.h \
	$R/src/workslice.h
	@# Special Compilation Rule
	${COMPILE_DYNALIGN_II}    $<

$R/src/rankmanager.o: \
	$R/src/rankmanager.cpp $R/src/rankmanager.h \
	$R/src/observable.h \
	$R/src/rank.h \
	$R/src/TProgressDialog.h \
	$R/src/workslice.h

$R/src/rankproducer.o: \
	$R/src/rank.h \
	$R/src/rankmanager.h \
	$R/src/rankproducer.cpp $R/src/rankproducer.h \
	$R/src/workunit.h

$R/src/rna_library.o: \
	$R/src/defines.h \
	$R/src/platform.h \
	$R/src/rna_library.cpp $R/src/rna_library.h \
	$R/src/structure.h

$R/src/pkHelix.o: \
	$R/src/pkHelix.cpp $R/src/pkHelix.h

$R/src/basepair.o: \
	$R/src/basepair.cpp $R/src/basepair.h

$R/src/Pseudoknot.o: \
	$R/src/Pseudoknot.cpp $R/src/Pseudoknot.h

$R/src/PseudoParser.o: \
	$R/src/PseudoParser.cpp $R/src/PseudoParser.h

$R/src/score.o: \
	$R/src/score.cpp $R/src/score.h

$R/src/siRNAfilter.o: \
	$R/src/siRNAfilter.cpp $R/src/siRNAfilter.h

$R/src/stackclass.o: \
	$R/src/defines.h \
	$R/src/stackclass.cpp $R/src/stackclass.h

$R/src/stackstruct.o: \
	$R/src/stackstruct.cpp $R/src/stackstruct.h

$R/src/stochastic.o: \
	$R/src/stochastic.cpp $R/src/stochastic.h

$R/src/stochastic-smp.o: \
	$R/src/stochastic.cpp $R/src/stochastic.h
	@# Special Compilation Rule
	${COMPILE_OMP}    $<

$R/src/structure.o: \
	$R/src/structure.cpp $R/src/structure.h \
	$R/src/defines.h \
	$R/src/platform.h

$R/src/StructureComparedImageHandler.o: \
	$R/src/StructureComparedImageHandler.cpp $R/src/StructureComparedImageHandler.h \
	$R/src/StructureImageHandler.cpp $R/src/StructureImageHandler.h

$R/src/StructureImageHandler.o: \
	$R/src/StructureImageHandler.cpp $R/src/StructureImageHandler.h

$R/src/thermo.o: \
	$R/src/thermo.cpp $R/src/thermo.h

$R/src/TProgressDialog.o: \
	$R/src/TProgressDialog.cpp $R/src/TProgressDialog.h

$R/src/varray.o: \
	$R/src/defines.h \
	$R/src/dynalign.h \
	$R/src/varray.cpp $R/src/varray.h

$R/src/wendarray.o: \
	$R/src/defines.h \
	$R/src/wendarray.cpp $R/src/wendarray.h

$R/stochastic/stochastic.o: \
	$R/stochastic/stochastic.cpp $R/stochastic/stochastic.h

$R/ShapeKnots/ShapeKnots_Interface.o: \
	$R/ShapeKnots/ShapeKnots_Interface.cpp $R/ShapeKnots/ShapeKnots_Interface.h \
	$R/src/ShapeKnots.h 

$R/src/ShapeKnots.o: \
	$R/src/ShapeKnots.cpp $R/src/ShapeKnots.h \
	$R/src/pkHelix.h \
	$R/src/PseudoParser.h \
	$R/src/rna_library.h \
	$R/src/structure.h \
	$R/src/algorithm.h \
	$R/src/ParseCommandLine.h \
	$R/RNA_class/RNA.h 

$R/TurboFold/TurboFold_Interface-smp.o:
	@# Special Compilation Rule
	${COMPILE_SMP} $R/TurboFold/TurboFold_Interface.cpp

$R/TurboFold/TurboFold_object-smp.o:
	@# Special Compilation Rule
	${COMPILE_SMP} $R/TurboFold/TurboFold_object.cpp

####  CUDA files ####
$R/fold-smp/frna.o: $R/fold-smp/frna.cu $R/fold-smp/frna.h \
	$R/fold-smp/int.h $R/fold-smp/fparam.h \
	$R/fold-smp/fbase.h $R/partition-smp/util.h
	@# Special Compilation Rule
	${COMPILE_CUDA} $R/fold-smp/frna.cu

$R/fold-smp/fparam.o: $R/fold-smp/fparam.c $R/fold-smp/fparam.h \
	$R/partition-smp/util.h
	@# Special Compilation Rule
	${COMPILE_CUDA} $R/fold-smp/fparam.c

$R/fold-smp/fbase.o: $R/fold-smp/fbase.c $R/fold-smp/fbase.h \
	$R/partition-smp/util.h
	@# Special Compilation Rule
	${COMPILE_CUDA} $R/fold-smp/fbase.c

# Note: Fold-cuda now uses util.o from partition-smp (but base.o is no longer compatible since partition-cuda was updated for the new datatables.)
$R/fold-smp/main.o: $R/fold-smp/main.c \
	$R/fold-smp/frna.h $R/fold-smp/fparam.h \
	$R/partition-smp/util.h $R/fold-smp/fbase.h
	@# Special Compilation Rule
	${COMPILE_CUDA} $R/fold-smp/main.c

$R/partition-smp/prna.o:  $R/partition-smp/prna.cu
	@# Special Compilation Rule
	${COMPILE_CUDA} $R/partition-smp/prna.cu

$R/partition-smp/param.o: $R/partition-smp/param.c
	@# Special Compilation Rule
	${COMPILE_CUDA} $R/partition-smp/param.c

$R/partition-smp/base.o: $R/partition-smp/base.c $R/partition-smp/base.h $R/partition-smp/util.h
	@# Special Compilation Rule
	${COMPILE_CUDA} $R/partition-smp/base.c

$R/partition-smp/util.o: $R/partition-smp/util.c $R/partition-smp/util.h
	@# Special Compilation Rule
	${COMPILE_CUDA} $R/partition-smp/util.c

$R/partition-smp/main.o: $R/partition-smp/main.c
	@# Special Compilation Rule
	${COMPILE_CUDA} $R/partition-smp/main.c

## CycleFold ##
$R/CycleFold/main.o: \
	$R/CycleFold/mainloop.h \
	$R/CycleFold/options.h \
	$R/CycleFold/arrays.h \
	$R/CycleFold/io.h \
	$R/CycleFold/NCM_parameters.h \
	$R/CycleFold/constraints.h \
	$R/CycleFold/maxexpect.h \
	$R/CycleFold/turbo_calculation.h

$R/CycleFold/main-test.o: \
	$R/CycleFold/main.cpp \
	$R/CycleFold/test.h \
	$(wildcard  $R/CycleFold/Tests/*.h )
	@# Special Compilation Rule
	${COMPILE_CPP}    $<

$R/CycleFold/arrays.o: $R/CycleFold/constants.h
$R/CycleFold/NCM_parameters.o: $R/CycleFold/constants.h
$R/CycleFold/io.o: \
	$R/CycleFold/io.h \
	$R/CycleFold/sequence.h

$R/CycleFold/mainloop.o: \
	$R/CycleFold/mainloop.h \
	$R/CycleFold/NCM_parameters.h \
	$R/CycleFold/extrinsic.h

PY_DIR=RNAstructure_python_interface
$R/$(PY_DIR)/RNAstructure_wrap.o: CXXFLAGS+=$(PYTHON_CXXFLAGS)
$R/$(PY_DIR)/RNAstructure_wrap.o: \
	$R/$(PY_DIR)/RNAstructure_wrap.cxx
	${COMPILE_CPP}   $<
