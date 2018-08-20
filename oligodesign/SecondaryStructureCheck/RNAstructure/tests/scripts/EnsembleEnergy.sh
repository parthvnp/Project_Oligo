# This script is not intended to be executed directly.  It is sourced by test-runner.sh and most functions you see here are defined in test-tools.sh
# The variables EXENAME, EXEBASE, OKDIR, and EXEDIR are defined externally.

beginTestBlock  # Begin a group of tests.
runMake -r @EXE # run `make` for the main exe and any other listed programs
EXT=.out  # set extension for @OKFILE and @OUTFILE placeholders (defaults for runDiff)

# Verify that required input files exist.
verifyInputFiles "$SINGLESEQ" "$SINGLEPFS"

# Test EnsembleEnergy_without_options.
runFullTest 'without_options' $SINGLEPFS ---stdout   # (Test output is redirected to @STDO)

# Test EnsembleEnergy_dna_option.
runFullTest 'dna_option'  $SINGLESEQ --sequence -d   ---stdout # (Test output is redirected to @STDO)

# Test EnsembleEnergy_sequence_option.
runFullTest 'sequence_option'  $SINGLESEQ --sequence    ---stdout # (Test output is redirected to @STDO)

# Test EnsembleEnergy_silent_option.
runFullTest 'silent_option'  $SINGLEPFS --silent    ---stdout # (Test output is redirected to @STDO)

endTestBlock # End a group of tests. Also cleans up orphaned files etc.
