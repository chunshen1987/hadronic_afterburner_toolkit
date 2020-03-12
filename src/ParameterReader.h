/***********************************************************************
ParameterReader class is used to simplify the process of reading parameters from an input file and/or from the command line.
Version 1.01 (09-20-2011) Zhi Qiu
***********************************************************************/

#ifndef _ParameterReaderHeader
#define _ParameterReaderHeader

#include <vector>
#include <string>

class ParameterReader {
  private:
    // store all parameter names and values
    std::vector<std::string> names;
    std::vector<double> values;

    // all substd::string after "symbol" in "str" will be removed
    std::string removeComments(std::string str, std::string commentSymbol); 

    // phrase an equation like "x=1", assume std::string has no comments
    void phraseEquationWithoutComments(std::string equation); 

    // give the index of parameter with "name", or -1 if it does not exist
    long find(std::string name) const; 

  public:
    ParameterReader();
    ~ParameterReader();

    // read and phrase one setting std::string like "x=1"
    void phraseOneLine(std::string str, std::string commentSymbol=(std::string)("#")); 

    // read in parameters from a file
    //void readFromFile(std::string filename, std::string commentSymbol=(std::string)("#")); 
    void readFromFile(std::string filename); 

    // read in parameter from argument list. 
    // The process starts with index="start_from".
    void readFromArguments(long argc, char * argv[], 
                           std::string commentSymbol=(std::string)("#"), 
                           long start_from=1); 

    bool exist(std::string name); // check if parameter with "name" exists

    // set the parameter with "name" to value "value"
    void setVal(std::string name, double value); 

    double getVal(std::string name) const; // return the value for parameter with "name"

    void echo(); // print out all parameters to the screen
};

#endif

/***********************************************************************
Changelog:
09-20-2011: Ver1.01
 -- Bug fix: If the parameter file that is passed to the readFromFile function does not exist, the program stops instead of going into infinite loops.

***********************************************************************/
