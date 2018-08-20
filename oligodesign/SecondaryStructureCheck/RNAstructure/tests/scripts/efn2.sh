# This script is not intended to be executed directly.  It is sourced by test-runner.sh and most functions you see here are defined in test-tools.sh
# The variables EXENAME, EXEBASE, OKDIR, and EXEDIR are defined externally.

beginTestBlock  # Begin a group of tests.
runMake -r @EXE # run `make` for the main exe and any other listed programs

# Verify that required input files exist.
verifyInputFiles "$SINGLECT"

EXT=.out  # set extension for @OKFILE and @OUTFILE placeholders (defaults for runDiff)

# Test efn2_without_options.
runFullTest 'without_options'  $SINGLECT @OUT.out

# Test efn2_dna_option.
runFullTest 'dna_option'  $SINGLECT @OUT.out -d

# Test efn2_alphabet_option.
runFullTest 'alphabet_option'  $SINGLECT @OUT.out -a 'dna' ---ref='dna_option'

# Test efn2_print_option.
initTest 'print_option'
runTest @EXE $SINGLECT @OUT.out -p
runDiff @OUT.out @OKFILE  # run diff on @OUT.out
runDiff @STDO @OKBASE_screen_OK.txt 'screen'  # run diff on the stdout of the previous command 
endTest

# Test efn2_shape_option.
runFullTest 'shape_option'  $SINGLECT @OUT.out -sh testFiles/testFile_tRNA.shape

# Test efn2_shape_intercept_option.
runFullTest 'shape_intercept_option'  $SINGLECT @OUT.out -sh testFiles/testFile_tRNA.shape -si 0.2

# Test efn2_shape_slope_option.
runFullTest 'shape_slope_option'  $SINGLECT @OUT.out -sh testFiles/testFile_tRNA.shape -sm 1.2
# efn2/efn2_shape_slope_option_OK.out

# Test efn2_temperature_option.
runFullTest 'temperature_option'  $SINGLECT @OUT.out -t 150

# Test efn2_write_thermodynamic_file_option.
runFullTest 'write_thermodynamic_file_option'  $SINGLECT @OUT.out -w

# Test efn2_pseudoknot.
runFullTest 'knotted'  testFiles/testFile_knotted.ct @OUT.out

# Test efn2_pseudoknot with write option.
runFullTest 'knotted_write'  testFiles/testFile_knotted.ct @OUT.out -w

endTestBlock # End a group of tests. Also cleans up orphaned files etc.
