/***********************************************************************
ParameterReader class is used to simplify the process of reading parameters from an input file and/or from the command line.
Version 1.01 (09-20-2011) Zhi Qiu
***********************************************************************/

#ifndef _ParameterReaderHeader
#define _ParameterReaderHeader

#include <vector>
#include <string>

using namespace std;

class ParameterReader
{
  private:
    // store all parameter names and values
    vector<string>* names; vector<double>* values; 

    // all substring after "symbol" in "str" will be removed
    string removeComments(string str, string commentSymbol); 

    // phrase an equation like "x=1", assume string has no comments
    void phraseEquationWithoutComments(string equation); 

    // give the index of parameter with "name", or -1 if it does not exist
    long find(string name); 

  public:
    ParameterReader();
    ~ParameterReader();

    // read and phrase one setting string like "x=1"
    void phraseOneLine(string str, string commentSymbol=(string)("#")); 

    // read in parameters from a file
    void readFromFile(string filename, string commentSymbol=(string)("#")); 

    // read in parameter from argument list. 
    // The process starts with index="start_from".
    void readFromArguments(long argc, char * argv[], 
                           string commentSymbol=(string)("#"), 
                           long start_from=1); 

    bool exist(string name); // check if parameter with "name" exists

    // set the parameter with "name" to value "value"
    void setVal(string name, double value); 

    double getVal(string name); // return the value for parameter with "name"

    void echo(); // print out all parameters to the screen
};

#endif

/***********************************************************************
Changelog:
09-20-2011: Ver1.01
 -- Bug fix: If the parameter file that is passed to the readFromFile function does not exist, the program stops instead of going into infinite loops.

***********************************************************************/
