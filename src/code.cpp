#include <Rcpp.h>
#include <fstream>
#include <sstream>
#include <string>
using namespace Rcpp;

CharacterVector read_file_cpp(std::string path) { // we need to figure out how to read compressed files
  std::ifstream t(path.c_str());
  std::stringstream ss;
  ss << t.rdbuf();
  return ss.str();
}


/**
Given: a vector of files paths
Return: a list of subsetted data matrices
 **/
// [[Rcpp::export]]
CharacterVector subsetCGfiles(StringVector CGFiles)
{
  String cur = CGFiles(0);
  Rcout << cur.get_cstring();
  CharacterVector file_raw = read_file_cpp(cur);
  NumericMatrix output(2, 3);
  return(file_raw);

}

// NumericMatrix subsetGCfiles(CharacterVector GCFiles)
// {
//
// }
