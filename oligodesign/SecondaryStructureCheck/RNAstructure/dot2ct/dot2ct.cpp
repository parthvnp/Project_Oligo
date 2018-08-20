/*
 * A program that converts a dot bracket file to a CT file.
 *
 * (c) 2010 Mathews Lab, University of Rochester Medical Center.
 * Written by Jessica S. Reuter
 */

#include "dot2ct.h"


///////////////////////////////////////////////////////////////////////////////
// Constructor.
///////////////////////////////////////////////////////////////////////////////
dot2ct_Interface::dot2ct_Interface() {

	// Initialize the calculation type description.
	calcType = "Dot bracket file conversion";
}

///////////////////////////////////////////////////////////////////////////////
// Parse the command line arguments.
///////////////////////////////////////////////////////////////////////////////
bool dot2ct_Interface::parse( int argc, char** argv ) {

	// Create the command line parser and build in its required parameters.
	ParseCommandLine* parser = new ParseCommandLine( "dot2ct" );
	parser->addParameterDescription( "bracket file", "The name of a file containing the dot bracket structure to convert." );
	parser->addParameterDescription( "ct file", "The name of a CT file to which output will be written." );

	// Parse the command line into pieces.
	parser->parseLine( argc, argv );

	// Get required parameters from the parser.
	if( !parser->isError() ) {
		bracketFile = parser->getParameter( 1 );
		ctFile = parser->getParameter( 2 );
	}

	// Delete the parser and return whether the parser encountered an error.
	bool noError = ( parser->isError() == false );
	delete parser;
	return noError;
}

///////////////////////////////////////////////////////////////////////////////
// Run calculations.
///////////////////////////////////////////////////////////////////////////////
bool dot2ct_Interface::run() {
	const char* const alphabet = DT_RNA;
	const bool allowUnkownBases = true;
	const bool skipThermo = true; // no need to load thermodynamic data tables.

	// Show a message saying that conversion has started.
	cout << "Converting dot bracket file..." << flush;

	RNA rna(bracketFile.c_str(), FILE_DBN, alphabet, allowUnkownBases, skipThermo);
	int error = rna.GetErrorCode();
	if (error!=0)
		cerr << rna.GetFullErrorMessage() << endl;
	else {
		error = rna.WriteCt(ctFile.c_str(), false);
		if (error!=0)
			cerr << "Error writing CT file: " << rna.GetErrorMessage(error) << endl;
	}

	// Print confirmation of run finishing.
	if( error == 0 ) { cout << calcType << " complete." << endl; }
	else { cerr << calcType << " encountered an error." << endl; }

	return error == 0;
}

///////////////////////////////////////////////////////////////////////////////
// Main method to run the program.
///////////////////////////////////////////////////////////////////////////////
int main( int argc, char* argv[] ) {
	dot2ct_Interface runner;
	if(!runner.parse(argc, argv)) 
		return 1;
	return runner.run() ? 0 : 1;
}
