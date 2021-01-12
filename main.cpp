#include <iostream>
#include <string>
#include "src.hpp"
#include "util.hpp"


int main(int argc, char **argv){

int wflag, bp;
std::string ftyp, atok;
setopts(argc, argv, &ftyp, &wflag, &bp, &atok);
runner(argv[argc-1], ftyp, bp, wflag, atok);
  return 0;
}
