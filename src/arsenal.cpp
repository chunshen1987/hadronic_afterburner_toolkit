// Ver 1.6.2
// Zhi Qiu
/*===========================================================================
Change logs: see arsenal.h
============================================================================*/

#include <iostream>
#include <string>
#include <sstream>
#include <stdlib.h>
#include <iomanip>
#include <cstdarg>

#include "arsenal.h"

#define OUTPUT_PRECISION 10

using namespace std;


//**********************************************************************
vector<double> stringToDoubles(string str)
// Return a vector of doubles from the string "str". "str" should
// be a string containing a line of data.
{
  // add a blank at the end so the last data will be read
  stringstream sst(str+" "); 
  vector<double> valueList;
  double val;
  sst >> val;
  while (sst.eof()==false)
  {
    valueList.push_back(val);
    sst >> val;
  }
  return valueList;
}


//**********************************************************************
double stringToDouble(string str)
// Return the 1st doubles number read from the string "str". 
// "str" should be a string containing a line of data.
{
  // add a blank at the end so the last data will be read
  stringstream sst(str+" "); 
  double val;
  sst >> val;
  return val;
}



//**********************************************************************
vector< vector<double>* >* readBlockData(istream &stream_in)
// Return a nested vector of vector<double>* object. Each column of data
// is stored in a vector<double> array and the collection is the returned
// object. Data are read from the input stream "stream_in". Each line
// of data is processed by the stringToDoubles function. Note that the
// data block is dynamicall allocated and is not release within the
// function.
// Note that all "vectors" are "new" so don't forget to delete them.
// Warning that also check if the last line is read correctly. Some files
// are not endded properly and the last line is not read.
{
  vector< vector<double>* >* data;
  vector<double> valuesInEachLine;
  long lineSize;
  long i; // temp variable
  char buffer[99999]; // each line should be shorter than this

  // first line:
  stream_in.getline(buffer,99999);
  valuesInEachLine = stringToDoubles(buffer);
  // see if it is empty:
  lineSize = valuesInEachLine.size();
  if (lineSize==0)
  {
    // empty:
    cout << "readBlockData warning: input stream has empty first row; "
         << "no data read" << endl;
    return NULL;
  }
  else
  {
    // not empty; allocate memory:
    data = new vector< vector<double>* >(lineSize);
    for (i=0; i<lineSize; i++) (*data)[i] = new vector<double>;
  }

  // rest of the lines:
  while (stream_in.eof()==false)
  {
    // set values:
    for (i=0; i<lineSize; i++) (*(*data)[i]).push_back(valuesInEachLine[i]);
    // next line:
    stream_in.getline(buffer,99999);
    valuesInEachLine = stringToDoubles(buffer);
  }

  return data;
}

//**********************************************************************
string toLower(string str)
// Convert all character in string to lower case
{
  string tmp = str;
  for (string::iterator it=tmp.begin(); it<=tmp.end(); it++) 
      *it = tolower(*it);
  return tmp;
}

//**********************************************************************
string trim(string str)
// Convert all character in string to lower case
{
  string tmp = str;
  long number_of_char = 0;
  for (size_t ii=0; ii<str.size(); ii++)
    if (str[ii]!=' ' && str[ii]!='\t')
    {
      tmp[number_of_char]=str[ii];
      number_of_char++;
    }
  tmp.resize(number_of_char);
  return tmp;
}

