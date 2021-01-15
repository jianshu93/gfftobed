#include <vector>
#include "util.hpp"
#ifndef SRC_H
#define SRC_H

struct gffs
{

  std::vector<std::string> chr, srcs, ftyps, atts, strnds, frames, scores, geneids;
  std::vector<int> starts, stops;
};



gffs gffread(const char* pfile, std::string att_tok);
gffs gtfread(const char* pfile, std::string att_tok);
void printgffs(gffs mygff);
void print_bed(gffs mygff, std::string feature);
void print_bed_window(gffs mygff, std::string feature, int wind);
void print_bed_upstream(gffs mygff, std::string feature, int wind);
void print_bed_downstream(gffs mygff, std::string feature, int wind);
void runner(const char* pfile, gfftobed_opts cmd_opts);

#endif
