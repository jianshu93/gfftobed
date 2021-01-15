#include <string>
#ifndef UTIL_H
#define UTIL_H

struct gfftobed_opts
{
  int interval, w_flag, gtf_flag;
  std::string feat, token;
};


gfftobed_opts setopts(int argc, char **argv);

void print_usage();

#endif
