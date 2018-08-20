# This script is not intended to be executed directly.  It is sourced by test-runner.sh and most functions you see here are defined in test-tools.sh
# The variables EXENAME, EXEBASE, OKDIR, and EXEDIR are defined externally.

beginTestBlock  # Begin a group of tests.
runMake -r Fold-cuda  # run `make` for the main exe and any other listed programs
EXT=.ct #ct file is final output from refold

# Attempts to use the MD5 hash to compare PFS files failed, because it seems that some
# binary content in the PFS file changes each time it is generated.
# function getMD5() {
#     if isExe md5; then
#         md5 -q "$1" > "$2"
#     elif isExe md5sum; then
#         local out=$(md5sum "$1")
#         echo ${out% *} > "$2"
#     else
#         errmsg "No MD5 utility found. Test cannot be performed."
#     fi
# }
#runSubTest 'md5' getMD5 @OUT.sav @OUT.md5
#runDiff @OUT.md5 @OK.md5 'md5'


# Test fold-cuda_without_options.
initTest 'without_options' && {
    runTest @EXE time/ivslsu.seq  @OUT.sav 
    runSubTest 'refold' refold @OUT.sav @OUT.ct
    runDiff #compare the CT files
endTest; }

#runFullTest 'dna_option' -d time/ivslsu.seq @OUT

EXT=.out
runFullTest 'array_option' -v time/ivslsu.seq @OUT.sav ---stdout  # File is quite large. 3.8 MB

endTestBlock # End a group of tests. Also cleans up orphaned files etc.
