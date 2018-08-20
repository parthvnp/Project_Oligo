/*
 * A program that converts a CT file to a dot bracket file.
 *
 * (c) 2010 Mathews Lab, University of Rochester Medical Center.
 * Written by Jessica S. Reuter
 * Rewritten by Richard M. Watson (2017) to use structure::writedotbracket
 */

#include "ct2dot.h"

///////////////////////////////////////////////////////////////////////////////
// Constructor.
///////////////////////////////////////////////////////////////////////////////
ct2dot_Interface::ct2dot_Interface() {

	// Initialize the calculation type description.
	calcType = "CT file conversion";
	format = DBN_FMT_MULTI_TITLE;
}

///////////////////////////////////////////////////////////////////////////////
// Parse the command line arguments.
///////////////////////////////////////////////////////////////////////////////
bool ct2dot_Interface::parse( int argc, char** argv ) {

	// Create the command line parser and build in its required parameters.
	ParseCommandLine* parser = new ParseCommandLine( "ct2dot" );
	parser->addParameterDescription( "ct file", "The name of a file containing the CT structure to convert." );
	parser->addParameterDescription( "structure number", "The number, one-indexed, of the structure to convert in the CT file (use -1 or \"ALL\" to convert all structures)." );
	parser->addParameterDescription( "bracket file", "The name of a dot bracket file to which output will be written." );

	vector<string> formatOption;
	formatOption.push_back("-f");
	formatOption.push_back("--format");
	parser->addOptionFlagsWithParameters( formatOption, "A number or name that indicates how subsequent sub-structures are formatted (relevant only when more than one structure is being written).\n"
		"Valid values are:\n"
		"\t(1) simple -- Susbequent structures (after the first one) are written with a Structure-Line  '(((....)))' only -- (no title or sequence)\n"
		"\t(2) side   -- A structure label is appended to the right side of each Structure-Line e.g. '(((....)))  ENERGY= -3.6  E.coli'.\n"
		"\t(3) multi  -- Susbequent structures are each written with a Title-Line '>TITLE' followed by a Structure-Line.\n"
		"\t(4) full   -- All structures written with a full header, including a '>TITLE' line followed by a Sequence-Line and then a Structure-Line.\n"
		"The default is 'multi'.");

	// Parse the command line into pieces.
	parser->parseLine( argc, argv );

	// Get required parameters from the parser.
	string numberString;
	if( !parser->isError() ) {
		ctFile = parser->getParameter( 1 );
		numberString = parser->getParameter( 2 );
		bracketFile = parser->getParameter( 3 );
	}

	// Convert the structure number parameter into an integer.
	// If that can't be done, set an error.
	toLower(numberString);
	if (numberString == "all")
		number = -1; // write all structures
	else {
		if (!parseInt(numberString, number))
			parser->setError( "structure number" );
	}

	// Convert the format specifier into a DotBracketFormat
	// If that can't be done, set an error.
	if ( !parser->isError() && parser->contains(formatOption) ) {
		string fmt = parser->getOptionString(formatOption, false);
		toLower(fmt);
		if (fmt=="simple"||fmt=="1") format = DBN_FMT_SINGLE_TITLE;
		else if (fmt=="side"||fmt=="2") format = DBN_FMT_SIDE_TITLES;
		else if (fmt=="multi"||fmt=="3") format = DBN_FMT_MULTI_TITLE;
		else if (fmt=="full"||fmt=="4") format = DBN_FMT_MULTI_TITLE_AND_SEQ;
		else 
			parser->setError( "format indicator" );
	}

	// Delete the parser and return whether the parser encountered an error.
	bool noError = ( parser->isError() == false );
	delete parser;
	return noError;
}

///////////////////////////////////////////////////////////////////////////////
// Run calculations.
///////////////////////////////////////////////////////////////////////////////
int ct2dot_Interface::run() {
	const char* const alphabet = DT_RNA;
	const bool allowUnkownBases = true;
	const bool skipThermo = true; // no need to load thermodynamic data tables.

	// Show a message saying that conversion has started.
	cout << "Converting CT file..." << flush;

	// Initialize and open the CT file to convert.
	RNA rna(ctFile.c_str(), FILE_CT, alphabet, allowUnkownBases, skipThermo);
	int error = rna.GetErrorCode();
	if (error!=0)
		cerr << rna.GetFullErrorMessage() << endl;
	else {
		error = rna.GetStructure()->writedotbracket(bracketFile.c_str(), number, format);
		if (error!=0)
			cerr << "Error writing file: " << rna.GetErrorMessage(error) << rna.GetStructure()->GetErrorDetails() << endl;
	}
	
	// Print confirmation of run finishing.
	if( error == 0 ) 
		cout << calcType << " complete." << endl;
	else 
		cerr << calcType << " encountered an error." << endl;

	return error == 0 ? 0 : 1;
}

///////////////////////////////////////////////////////////////////////////////
// Main method to run the program.
///////////////////////////////////////////////////////////////////////////////
int main( int argc, char* argv[] ) {
	ct2dot_Interface runner;
	if (!runner.parse( argc, argv )) return 1;
	return runner.run();
}
