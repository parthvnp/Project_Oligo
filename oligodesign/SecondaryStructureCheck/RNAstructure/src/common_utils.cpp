#include "common_utils.h"
#include <iostream>
#include <string>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cstdarg>
#include <sstream>
#include <fstream>
#include <iterator>
#include <algorithm>
#include <vector>
#include <sys/stat.h>
#include <cerrno>
#include <climits>

using namespace std;

//Determine if the path <directory>/<file> exists and (is not itself a directory)
bool fileExists(const char* const directory, const char* const file) {
	if (directory==NULL||file==NULL) return false;
	string filePath(directory);
	return fileExists(filePath.append("/").append(file).c_str());
}
bool fileExists(const char* const fullPath, const bool verifyReadable) {
	if (verifyReadable) {
		// Verify that the file can be opened for reading.
		ifstream in_test(fullPath);
		return in_test.good();
	} else {
		// Verify only that the file exists.
		struct stat info;
		return fullPath != NULL && stat(fullPath, &info)==0 && (info.st_mode & S_IFDIR)==0;
	}
}
bool dirExists(const char* const fullDirPath) {
	struct stat info;
	return fullDirPath != NULL && stat(fullDirPath, &info)==0 && (info.st_mode&S_IFDIR)==S_IFDIR;
}

// Gets the portion of a path that represents the filename.
// I.e. removes the directory portion of a path. 
// If removeExtension is true, the file extension is also removed.
string getFileName(const char * const path, bool removeExtension) {
	string filename(path);
	size_t pos = filename.find_last_of("/\\");
	if (pos != string::npos) filename.erase(0, pos + 1);
	if (removeExtension) {
		pos = filename.rfind('.');
		if (pos != string::npos)
			filename.erase(pos);
	}
	return filename;
}

//! Gets the portion of a path that represents the file extension 
//! i.e. the part of the file name after the last dot (.)
string getFileExt(const string& filePath) {
	size_t slash = filePath.find_last_of("/\\"); // find last occurance of a slash
	if (slash == string::npos) slash = 0;
	size_t dot = filePath.rfind('.'); // find last dot
	if (dot == string::npos || dot < slash) return ""; // either there is no dot or the last dot is in the directory name. Either way, the file has no extension.
	return filePath.substr(dot+1);
}

// Gets the directory portion of a path. 
// It is assumed that the path represents a file.
string getDirName(const char * const path) {
	string dir(path);
	size_t pos = dir.find_last_of("/\\");
	if (pos==string::npos) return "."; // e.g. input is "hello.txt" the directory is "." -- the current directory
	dir.resize(pos); // set the length to just before the final slash
	return dir;
}

// Creates string and fills it by sprintf-formatting the arguments.
string sfmt(const char* const format, ...) {
    va_list args;
	int size = strlen(format) + 256;
	char* buf = new char[size];
    va_start(args, format);
    int req_size = vsnprintf(buf, size, format, args);
    va_end(args);
	if (req_size < 0)
		// negative return values indicate an error. 
		req_size = sprintf(buf, "Error formatting arguments: %d", req_size);
	else if (req_size >= size) {  // indicates required buffer size (not including terminating \0)
		// the formatted arguments could not all fit inside the buffer. So resize the buffer to the exact required amount and repeat.
		delete[] buf;
		size = req_size+1; // +1 for \0
		buf = new char[size];
		va_start(args, format);
		vsnprintf(buf, size, format, args);
		va_end(args);
	}
	string ret = buf; // copy the buffer and release it.
	delete[] buf;
	return ret;
}

// Find a char in a c-string
size_t findchr(const char* const subject, const char find) {
	const char* found = strchr(subject, find); // strchr returns a pointer to the char that was found or NULL (if not found)
	return found==NULL ? string::npos : found - subject;
}

// converts control characters in subject into escape sequences.
// (useful for debugging strings that might contain control characters.)
// e.g. "Hello\n" becomes "Hello\\n" 
// i.e. The linefeed is converted into a literal slash and the letter 'n'
string& escapeChars(string &subject) {
	string snew;
	snew.reserve((int)(subject.size() * 1.3));
	char numbuf[5];
	for(string::iterator it=subject.begin(); it != subject.end(); it++) {
		char &c = *it;
		if (c < 32 || c > 126)
			switch(c) {
				case '\n': snew+="\\n"; break; 
				case '\r': snew+="\\r"; break;
				case '\t': snew+="\\t"; break;
				case '\0': snew+="\\0"; break;
				default:
					snew += "\\x";
					sprintf(numbuf, "%02X", c);
					snew += numbuf;
					break;
			}
		else
			snew += c;
	}
	subject.swap(snew);
	return subject;
}

string& trimLeft(string &s) {
	string::iterator it;
	for (it = s.begin(); it!=s.end(); it++) if (!::isspace(*it)) break; // we are guaranteed to find a non-whitespace character (the final '\0' at s.end() if nothing else)
	s.erase(s.begin(), it);
	return s;
}
string& trimRight(string &s) {
	string::iterator it;
	for (it = s.end()-1; it>=s.begin(); it--) if (!::isspace(*it)) { it++; break; } // start with s.end()-1 to skip the final '\0' at s.end()
	if (it < s.begin()) it++; // increment so we don't delete the NON-whitespace character that was found.
	s.erase(it, s.end());
	return s;
}
string& trim(string &s) {
	trimLeft(s);
	if (!s.empty())
		trimRight(s);
	return s;
}
// Converts a string to lower-case (operates on subject in-place).
string& toLower(string &s) {
	if (!s.empty())
		transform(s.begin(), s.end(), s.begin(), ::tolower);
	return s;
}
// Converts a string to upper-case (operates on subject in-place).
string& toUpper(string &s) {
	if (!s.empty())
		transform(s.begin(), s.end(), s.begin(), ::toupper);
	return s;
}

// The following string functions operate on const strings and 
// return a copy of the string with the desired modifiations.
string trimLeft(const string &subject) { string copy(subject); return trimLeft(copy);  }
string trimRight(const string &subject) { string copy(subject); return trimRight(copy);  }
string trim(const string &subject) { string copy(subject); return trim(copy);  }
string toLower(const string &subject) { string copy(subject); return toLower(copy);  }
string toUpper(const string &subject) { string copy(subject); return toUpper(copy);  }

 char* copy_cstr(const char* src) {
    if (src == NULL) return NULL;
    char* ptr = new char[strlen(src)+1];
    strcpy(ptr, src);
    return ptr;
}

// Parse a c-string to get an integer.
// Returns true on success or false if the string could not be converted to an interger.
// Note that intResult is NOT modified if conversion fails.
bool parseInt(const char* const input, int& intResult, const bool readToEnd) {
	long lnum; char *endptr;
	errno = 0; // reset global error indicator
	// long int strtol (const char* input, char** endptrptr, int base);
	lnum = strtol(input, &endptr, 0); // 0 specifies parsing of standard radix options (including 0x0 077, etc)
	if ( endptr==input //if no characters were converted, these pointers are equal.
		  // explicitly check for overflows. note that long *might* be the same size as int on some systems.
		  || errno == ERANGE || (lnum > INT_MAX) || (lnum < INT_MIN)  // number out of range
		  || (readToEnd&&(*nextNonWhitespaceChar(endptr)!='\0')) // extra (non-whitepsace) text after number.
		) return false;
	 //Success. Convert the result to a plain int.
	intResult = (int)lnum;
	return true;
}

// Parse a c-string to get a double.
// Return true on success or false if the string could not be converted to a double.
// Note that dblResult is NOT modified if conversion fails.
bool parseDbl(const char* input, double& dblResult, const bool readToEnd) {
	double num; char *endptr;
	errno = 0; // reset global error indicator
	//double strtod(const char *input, char **endptr);
	num = strtod(input, &endptr);
	if ( endptr==input //if no characters were converted these pointers would be equal
		 || errno!=0 // check global error indicator for non-zero values (e.g. ERANGE)
		 || (readToEnd&&(*nextNonWhitespaceChar(endptr)!='\0'))
		 ) return false;
	//Success.
	dblResult = num;
	return true;
}

// split a string using a delimiter.
vector<string> split(const string& str, const string& delim, const bool includeEmptyValues) {
    vector<string> tokens;
    size_t prev = 0, pos = 0;
    while (pos < str.length()) {
        pos = str.find(delim, prev);
        if (pos == string::npos) pos = str.length();
        string token = str.substr(prev, pos-prev);
        //cout<<"\t"<<includeEmptyValues<<"\t"<<(pos>prev)<<endl;
        if (includeEmptyValues||pos>prev) tokens.push_back(token);
        prev = pos + delim.length();
    }
    return tokens;
}

NullStream NullStream::Default; // construct default member.

//###################### Operating System Compatibility ##################################
// Windows's printf function incorrectly outputs 3-digit exponents when only 2 digits 
// are required, violating the C99+ spec. e.g. Windows: "1e-002" vs Linux: "1e-02"
// This causes issues with regression tests due to different output on each OS.
//
// Calling this function (before relevant output) fixes the problem.
void FixWindowsExponentFormat() {
	#if defined( WIN32 )
		#if (_MSC_VER >= 1400) && (_MSC_VER < 1900)
			// The function _set_output_format was added in Visual Studio 2005 (_MSC_VER=1400)
			// and later removed in Visual Studio 2015 (_MSC_VER=1900) which defaults to conformant output.
			_set_output_format(_TWO_DIGIT_EXPONENT);
		#elif defined( __MINGW32__ )
			// _set_output_format works in MinGW, but the following is probably more reliable
			// (and at the very least should not cause problems in later versions)
			if( getenv( "PRINTF_EXPONENT_DIGITS" ) == NULL ) 
				_putenv( "PRINTF_EXPONENT_DIGITS=2" );
		#endif
	#endif
}
