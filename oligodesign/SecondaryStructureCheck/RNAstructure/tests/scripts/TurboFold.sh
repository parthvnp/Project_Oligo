# This script is not intended to be executed directly.  It is sourced by test-runner.sh and most functions you see here are defined in test-tools.sh
# The variables EXENAME, EXEBASE, OKDIR, and EXEDIR are defined externally.

beginTestBlock  # Begin a group of tests.
runMake -r @EXE  # run `make` for the main exe and any other listed programs
EXT=.ct  # set extension for @OKFILE and @OUTFILE placeholders (defaults for runDiff)

SMP_SETTING=  # SMP_SETTING is blank unless this is the SMP variant of TurboFold in which case, we set the number of processors (minimum 2)
if [[ $SMP ]]; then 
	N_PROCS=$({ getconf _NPROCESSORS_ONLN || getconf NPROCESSORS_ONLN ; } 2>/dev/null) # get number of processors (getconf is POSIX, but _NPROCESSORS_ONLN is optional)
	[[ $N_PROCS -gt 1 ]] || N_PROCS=2 # If N_PROCS is empty or 0 or 1, then set it to 2.
	SMP_SETTING="Processors = $N_PROCS"
	isQuiet || warnmsg "Using '$SMP_SETTING' in SMP mode."
fi

SHAPEFILE1=testFiles/testFile_random_dummy.shape
SHAPEFILE2=testFiles/testFile_random_dummy_114.shape # remove bases over 114 due to shorter sequence
SHAPEFILE3=$SHAPEFILE1

# Run a TurboFold test. This includes:
#   1. Writing the turbofold configuration file
#   2. running TurboFold
#   3. Comparing the output (alignment file and CT files) with references ("OK" files)
# Usage: 
#   turboTest <TEST_NAME>  <CONF_TEXT> [--standard | <SEQ_NAMES...>]
# Where SEQ_NAMES is a list of sequence names for which the CT file output 
# should be compared to analogously-named reference ("OK") files.
function turboTest() { 
	# Initilize the test with the testing system
	initTest "$1" || return  # exit if the test is excluded etc.
	
	# Write the conf file
	# Note that whereas runTest and runDiff automatically replace test placeholders (such as @TEST and @OUT etc),
	# we have to do this explicitly with `replaceVarsInText` because we are writing to a file.
	replaceVarsInText "$2" > "$TEST.conf"  # Write output to $TEST.conf e.g. "TurboFold_general_without_options.conf"
	
	# Now run TurboFold with the given conf
	runTest @EXE "$TEST.conf"

	# Now compare the output alignment files and (optionally) CT files
	# We can compare with the "standard" output files (from the general_without_options test)
	# or with a custom set of reference files.
	local okfile sequences
	if [[ $3 == --standard ]]; then
		setTestRef 'general_without_options'  # use the  'general_without_options' OK files as references for this test. This affects @OKBASE, @OK, and @OKFILE
		set -- seq1 seq2 seq3  # set the argument list to these three names
	else
		shift 2   # sequences are listed as all remaining arguments (after the first two)
	fi

	# Run a diff on the alignment file 
	# diff output is written to ${TEST}_aln_diff.txt e.g.  "TurboFold_general_without_options_aln_diff.txt"
	runDiff @OUT.aln @OK.aln 'aln'  

	# Now for each of the listed input sequences, do a diff with the expected output CT
	local seqName 
	for seqName; do  # loop through each sequence name in the remaining arguments (i.e. "$@")
		runDiff "@TEST_${seqName}_output.ct"  "@OKBASE_${seqName}_OK.ct"  "${seqName}"  # diff output written to  ${TEST}_${seqName}_diff.txt
	done

	endTest # perform end-of-test duties
}

# Test TurboFold_without_options_singly.
# Note that the general test has a mode specified of MEA; MEA is the default mode.
turboTest 'general_without_options_singly' "
Seq1 = testFiles/testFile_ec5s.seq
Seq2 = testFiles/testFile_ca5s.seq
Seq3 = testFiles/testFile_5s_8.seq
CT1 = @TEST_seq1_output.ct
CT2 = @TEST_seq2_output.ct
CT3 = @TEST_seq3_output.ct
OutAln = @OUT.aln
SequenceNumber = 3
Mode = MEA
$SMP_SETTING"   --standard

# Test fasta input files.
turboTest 'general_without_options_fasta' "
Seq1 = testFiles/testFile_ec5s.fasta
Seq2 = testFiles/testFile_ca5s.fasta
Seq3 = testFiles/testFile_5s_8.fasta
CT1 = @TEST_seq1_output.ct
CT2 = @TEST_seq2_output.ct
CT3 = @TEST_seq3_output.ct
OutAln = @OUT.aln
SequenceNumber = 3
Mode = MEA
$SMP_SETTING"   --standard

# Test TurboFold_without_options_groups.
# Note that the general test has a mode specified of MEA; MEA is the default mode.
turboTest 'general_without_options_groups' "
InSeq = {testFiles/testFile_ec5s.seq;testFiles/testFile_ca5s.seq;testFiles/testFile_5s_8.seq;}
OutCT = {@TEST_seq1_output.ct;@TEST_seq2_output.ct;@TEST_seq3_output.ct;}
OutAln = @OUT.aln
Mode = MEA
$SMP_SETTING"   --standard

# Test an alternate TurboFold_general_without_options.
# Note that the general test has a mode specified of MEA; MEA is the default mode.
# Note also that these alternate test files are used to test ProbKnot mode only.
# This alternate test is needed to ensure that the ProbKnot_without_options test is different from a default test where ProbKnot mode isn't specified.
turboTest 'general_without_options_alternate' "
Seq1 = testFiles/testFile_met-vol.seq
Seq2 = testFiles/testFile_met-fer.seq
Seq3 = testFiles/testFile_met-the.seq
CT1 = @TEST_seq1_output.ct
CT2 = @TEST_seq2_output.ct
CT3 = @TEST_seq3_output.ct
OutAln = @OUT.aln
SequenceNumber = 3
$SMP_SETTING
Mode = MEA"  seq1 seq2 seq3

# Test TurboFold_general_gamma_option.
# Note that the general test has a mode specified of MEA; MEA is the default mode.
turboTest 'general_gamma_option' "
Seq1 = testFiles/testFile_ec5s.seq
Seq2 = testFiles/testFile_ca5s.seq
Seq3 = testFiles/testFile_5s_8.seq
CT1 = @TEST_seq1_output.ct
CT2 = @TEST_seq2_output.ct
CT3 = @TEST_seq3_output.ct
OutAln = @OUT.aln
SequenceNumber = 3
$SMP_SETTING
Mode = MEA
Gamma = 0.9"  seq1 seq2 seq3

# Test TurboFold_general_iterations_option.
# Note that the general test has a mode specified of MEA; MEA is the default mode.
turboTest 'general_iterations_option' "
Seq1 = testFiles/testFile_ec5s.seq
Seq2 = testFiles/testFile_ca5s.seq
Seq3 = testFiles/testFile_5s_8.seq
CT1 = @TEST_seq1_output.ct
CT2 = @TEST_seq2_output.ct
CT3 = @TEST_seq3_output.ct
OutAln = @OUT.aln
SequenceNumber = 3
$SMP_SETTING
Mode = MEA
Iterations = 1"  seq1 seq2 seq3

# Test TurboFold_general_shape_option.
# Note that the general test has a mode specified of MEA; MEA is the default mode.
turboTest 'general_shape_option' "
Seq1 = testFiles/testFile_ec5s.seq
Seq2 = testFiles/testFile_ca5s.seq
Seq3 = testFiles/testFile_5s_8.seq
CT1 = @TEST_seq1_output.ct
CT2 = @TEST_seq2_output.ct
CT3 = @TEST_seq3_output.ct
SequenceNumber = 3
$SMP_SETTING
Mode = MEA
OutAln = @OUT.aln
SHAPE1 = $SHAPEFILE1
SHAPE2 = $SHAPEFILE2
SHAPE3 = $SHAPEFILE3"  seq1 seq2 seq3


# Test TurboFold_general_shape_intercept_option.
# Note that the general test has a mode specified of MEA; MEA is the default mode.
turboTest 'general_shape_intercept_option' "
Seq1 = testFiles/testFile_ec5s.seq
Seq2 = testFiles/testFile_ca5s.seq
Seq3 = testFiles/testFile_5s_8.seq
CT1 = @TEST_seq1_output.ct
CT2 = @TEST_seq2_output.ct
CT3 = @TEST_seq3_output.ct
OutAln = @OUT.aln
SequenceNumber = 3
$SMP_SETTING
Mode = MEA
OutAln = @OUT.aln
SHAPE1 = $SHAPEFILE1
SHAPE2 = $SHAPEFILE2
SHAPE3 = $SHAPEFILE3
SHAPEintercept = 0.1"  seq1 seq2 seq3


# Test TurboFold_general_shape_slope_option.
# Note that the general test has a mode specified of MEA; MEA is the default mode.
turboTest 'general_shape_slope_option' "
Seq1 = testFiles/testFile_ec5s.seq
Seq2 = testFiles/testFile_ca5s.seq
Seq3 = testFiles/testFile_5s_8.seq
CT1 = @TEST_seq1_output.ct
CT2 = @TEST_seq2_output.ct
CT3 = @TEST_seq3_output.ct
SequenceNumber = 3
$SMP_SETTING
Mode = MEA
OutAln = @OUT.aln
SHAPE1 = $SHAPEFILE1
SHAPE2 = $SHAPEFILE2
SHAPE3 = $SHAPEFILE3
SHAPEslope = 0.1"   seq1 seq2 seq3


# Test TurboFold_general_temperature_option.
# Note that the general test has a mode specified of MEA; MEA is the default mode.
turboTest 'general_temperature_option' "
Seq1 = testFiles/testFile_ec5s.seq
Seq2 = testFiles/testFile_ca5s.seq
Seq3 = testFiles/testFile_5s_8.seq
CT1 = @TEST_seq1_output.ct
CT2 = @TEST_seq2_output.ct
CT3 = @TEST_seq3_output.ct
OutAln = @OUT.aln
SequenceNumber = 3
$SMP_SETTING
Mode = MEA
Temperature = 330"   seq1 seq2 seq3


[[ $SMP ]] &&  # If the executable is SMP, test the processors option.
turboTest 'processors_option' "
Seq1 = testFiles/testFile_ec5s.seq
Seq2 = testFiles/testFile_ca5s.seq
Seq3 = testFiles/testFile_5s_8.seq
CT1 = @TEST_seq1_output.ct
CT2 = @TEST_seq2_output.ct
CT3 = @TEST_seq3_output.ct
OutAln = @OUT.aln
SequenceNumber = 3
Processors = 3;
Mode = MEA"   --standard


# Test TurboFold_MEA-mode_without_options.
turboTest 'MEA-mode_without_options' "
Seq1 = testFiles/testFile_ec5s.seq
Seq2 = testFiles/testFile_ca5s.seq
Seq3 = testFiles/testFile_5s_8.seq
CT1 = @TEST_seq1_output.ct
CT2 = @TEST_seq2_output.ct
CT3 = @TEST_seq3_output.ct
OutAln = @OUT.aln
SequenceNumber = 3
$SMP_SETTING
Mode = MEA"   --standard

# Test an alternate TurboFold_MEA-mode_without_options.
# These test files are not used for all MEA mode tests, just max structures, percent difference, and window size.
turboTest 'MEA-mode_without_options_alternate' "
Seq1 = testFiles/testFile_RA7680.seq
Seq2 = testFiles/testFile_RD0260.seq
Seq3 = testFiles/testFile_RD0500.seq
CT1 = @TEST_seq1_output.ct
CT2 = @TEST_seq2_output.ct
CT3 = @TEST_seq3_output.ct
OutAln = @OUT.aln
SequenceNumber = 3
$SMP_SETTING
Mode = MEA"   seq1 seq2 seq3

# Test TurboFold_MEA-mode_gamma_option.
turboTest 'MEA-mode_gamma_option' "
Seq1 = testFiles/testFile_ec5s.seq
Seq2 = testFiles/testFile_ca5s.seq
Seq3 = testFiles/testFile_5s_8.seq
CT1 = @TEST_seq1_output.ct
CT2 = @TEST_seq2_output.ct
CT3 = @TEST_seq3_output.ct
OutAln = @OUT.aln
SequenceNumber = 3
$SMP_SETTING
Mode = MEA
MeaGamma = 0.2"   seq1 seq2 seq3

# Test TurboFold_MEA-mode_max_structures_option.
turboTest 'MEA-mode_max_structures_option' "
Seq1 = testFiles/testFile_RA7680.seq
Seq2 = testFiles/testFile_RD0260.seq
Seq3 = testFiles/testFile_RD0500.seq
CT1 = @TEST_seq1_output.ct
CT2 = @TEST_seq2_output.ct
CT3 = @TEST_seq3_output.ct
OutAln = @OUT.aln
SequenceNumber = 3
$SMP_SETTING
Mode = MEA
MaxStructures = 5
Window = 0"   seq1 seq2 seq3

# Test TurboFold_MEA-mode_percent_difference_option.
turboTest 'MEA-mode_percent_difference_option' "
Seq1 = testFiles/testFile_RA7680.seq
Seq2 = testFiles/testFile_RD0260.seq
Seq3 = testFiles/testFile_RD0500.seq
CT1 = @TEST_seq1_output.ct
CT2 = @TEST_seq2_output.ct
CT3 = @TEST_seq3_output.ct
OutAln = @OUT.aln
SequenceNumber = 3
$SMP_SETTING
Mode = MEA
MaxPercent = 4
Window = 0"   seq1 seq2 seq3

# Test TurboFold_MEA-mode_window_size_option.
turboTest 'MEA-mode_window_size_option' "
Seq1 = testFiles/testFile_RA7680.seq
Seq2 = testFiles/testFile_RD0260.seq
Seq3 = testFiles/testFile_RD0500.seq
CT1 = @TEST_seq1_output.ct
CT2 = @TEST_seq2_output.ct
CT3 = @TEST_seq3_output.ct
OutAln = @OUT.aln
SequenceNumber = 3
$SMP_SETTING
Mode = MEA
Window = 0"   seq1 seq2 seq3

# Test TurboFold_ProbKnot-mode_without_options.
turboTest 'ProbKnot-mode_without_options' "
Seq1 = testFiles/testFile_met-vol.seq
Seq2 = testFiles/testFile_met-fer.seq
Seq3 = testFiles/testFile_met-the.seq
CT1 = @TEST_seq1_output.ct
CT2 = @TEST_seq2_output.ct
CT3 = @TEST_seq3_output.ct
OutAln = @OUT.aln
SequenceNumber = 3
$SMP_SETTING
Mode = ProbKnot"   seq1 seq2 seq3

# Test TurboFold_ProbKnot-mode_iterations_option.
turboTest 'ProbKnot-mode_iterations_option' "
Seq1 = testFiles/testFile_met-vol.seq
Seq2 = testFiles/testFile_met-fer.seq
Seq3 = testFiles/testFile_met-the.seq
CT1 = @TEST_seq1_output.ct
CT2 = @TEST_seq2_output.ct
CT3 = @TEST_seq3_output.ct
OutAln = @OUT.aln
SequenceNumber = 3
$SMP_SETTING
Mode = ProbKnot
PkIterations = 10"   seq1 seq2 seq3

# Test TurboFold_ProbKnot-mode_min_helix_length_option.
turboTest 'ProbKnot-mode_min_helix_length_option' "
Seq1 = testFiles/testFile_met-vol.seq
Seq2 = testFiles/testFile_met-fer.seq
Seq3 = testFiles/testFile_met-the.seq
CT1 = @TEST_seq1_output.ct
CT2 = @TEST_seq2_output.ct
CT3 = @TEST_seq3_output.ct
OutAln = @OUT.aln
SequenceNumber = 3
$SMP_SETTING
Mode = ProbKnot
MinHelixLength = 1"   seq1 seq2 seq3

# Test TurboFold_Threshold-mode_without_options.
turboTest 'Threshold-mode_without_options' "
Seq1 = testFiles/testFile_ec5s.seq
Seq2 = testFiles/testFile_ca5s.seq
Seq3 = testFiles/testFile_5s_8.seq
CT1 = @TEST_seq1_output.ct
CT2 = @TEST_seq2_output.ct
CT3 = @TEST_seq3_output.ct
OutAln = @OUT.aln
SequenceNumber = 3
$SMP_SETTING
Mode = Threshold"   seq1 seq2 seq3

# Test TurboFold_Threshold-mode_threshold_option.
turboTest 'Threshold-mode_threshold_option' "
Seq1 = testFiles/testFile_ec5s.seq
Seq2 = testFiles/testFile_ca5s.seq
Seq3 = testFiles/testFile_5s_8.seq
CT1 = @TEST_seq1_output.ct
CT2 = @TEST_seq2_output.ct
CT3 = @TEST_seq3_output.ct
SequenceNumber = 3
$SMP_SETTING
Mode = Threshold
OutAln = @OUT.aln
Threshold = 0.92"   seq1 seq2 seq3

# Test TurboFold_Rsample_option.
turboTest 'Rsample_option' "
InSeq = {TurboFold/input/tRNA.seq;TurboFold/input/tRNA2.seq;TurboFold/input/RD0260.seq;TurboFold/input/RD0500.seq;TurboFold/input/RA7680.seq;}
OutCT = {@TEST_tRNA_output.ct;@TEST_tRNA2_output.ct;@TEST_RD0260_output.ct;@TEST_RD0500_output.ct;@TEST_RA7680_output.ct;}
Mode = MEA
UseRsample = 1
SHAPEFiles = {TurboFold/input/tRNA.shape;TurboFold/input/tRNA2.shape;;;;}
Seed = 1
OutAln = @OUT.aln
$SMP_SETTING"   tRNA tRNA2 RD0260 RD0500 RA7680

# Should give an error about a missing shape file.
initTest 'shape_file_missing' && {
replaceVarsInText "InSeq = {testFiles/testFile_ec5s.seq;testFiles/testFile_ca5s.seq;}
OutCT = {@TEST_1.ct;@TEST_2.ct}
Mode = MEA
SHAPEFiles = {missing_shape_file.shape;}
$SMP_SETTING" > "$TEST.conf"
# run the test, but expect an exit code of 1 and perform the diff on the output from stderr
runProgAndDiff "$TEST.conf"   ---exit=1 ---stderr 
endTest; }

# SHAPEFILE1 is too long for ca5s.seq -- should give an error about an invalid nucleotide position
initTest 'shape_file_bad' && {
replaceVarsInText "InSeq = {testFiles/testFile_ec5s.seq;testFiles/testFile_ca5s.seq;}
OutCT = {@TEST_1.ct;@TEST_2.ct}
OutAln = @OUT.aln
Mode = MEA
SHAPEFiles = {$SHAPEFILE1;$SHAPEFILE1}
$SMP_SETTING" > "$TEST.conf"
# run the test, but expect an exit code of 1 and perform the diff on the output from stderr
# set RNA_WARNINGS environment variable to ERR, which causes warnings to be written to STDERR instead of STDOUT.
RNA_WARNINGS=ERR \
runProgAndDiff "$TEST.conf"   ---stderr
endTest; }


endTestBlock # End a group of tests. Also cleans up orphaned files etc.
