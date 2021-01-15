#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <sstream>
#include "src.hpp"
#include "util.hpp"



gffs gffread(const char* pfile, std::string att_tok) {

  std::ifstream gff(pfile);
  std::string line;
  std::string hsh = "#";
  std::string seqid, src, ftyp, attribs, seqstrand, seqframe, seqscore;
  int seqstart, seqstop;
  std::string ttok = att_tok + "=";



  gffs gs;

if(gff.is_open()){
/* read the gff without hash lines */
  while(std::getline(gff, line)){
    if(line.find(hsh) == std::string::npos){
          std::stringstream myline(line);
          while(myline >> seqid >> src >> ftyp >> seqstart >> seqstop >> seqscore >> seqstrand >> seqframe && std::getline(myline,attribs)){

            gs.chr.push_back(seqid);
            gs.srcs.push_back(src);
            gs.ftyps.push_back(ftyp);
            gs.starts.push_back(seqstart);
            gs.stops.push_back(seqstop);
            gs.strnds.push_back(seqstrand);
            gs.frames.push_back(seqframe);
            gs.scores.push_back(seqscore);
            gs.atts.push_back(attribs);

            /* Search for token in attribs, if not found use whole attribs */
            if(attribs.find(ttok) == std::string::npos){
              gs.geneids.push_back(attribs);
            } else {
              int substart = attribs.find(ttok)+ttok.length();
              int subend = attribs.find(";", substart);
              gs.geneids.push_back(attribs.substr(substart, subend-substart));}

          }
        }
  }
}else {
  std::cout << "ERROR: File was not able to be opened." << std::endl;
  exit (EXIT_FAILURE);
}


return gs;
}//end gffread








gffs gtfread(const char* pfile, std::string att_tok) {

  std::ifstream gff(pfile);
  std::string line;
  std::string hsh = "#";
  std::string seqid, src, ftyp, attribs, seqstrand, seqframe, seqscore;
  int seqstart, seqstop;
  std::string ttok = att_tok + " \"";



  gffs gs;

if(gff.is_open()){
/* read the gff without hash lines */
  while(std::getline(gff, line)){
    if(line.find(hsh) == std::string::npos){
          std::stringstream myline(line);
          while(myline >> seqid >> src >> ftyp >> seqstart >> seqstop >> seqscore >> seqstrand >> seqframe && std::getline(myline,attribs)){

            gs.chr.push_back(seqid);
            gs.srcs.push_back(src);
            gs.ftyps.push_back(ftyp);
            gs.starts.push_back(seqstart);
            gs.stops.push_back(seqstop);
            gs.strnds.push_back(seqstrand);
            gs.frames.push_back(seqframe);
            gs.scores.push_back(seqscore);
            gs.atts.push_back(attribs);

            if(attribs.find(ttok) == std::string::npos){
              gs.geneids.push_back(attribs);
            } else {
              int substart = attribs.find(ttok)+ttok.length();
              int subend = attribs.find("\";", substart);
              gs.geneids.push_back(attribs.substr(substart, subend-substart));}

          }
        }
  }
}else {
  std::cout << "ERROR: File was not able to be opened." << std::endl;
  exit (EXIT_FAILURE);
}


return gs;
}//end gtfread










/* general print function */
 void printgffs(gffs mygff){

    for(size_t i=0; i<10; i++){
        std::cout << mygff.chr[i] << "\t"<< mygff.srcs[i]<< "\t"<< mygff.ftyps[i]<< "\t"<<mygff.starts[i]<< "\t"<<mygff.stops[i]<< "\t"<<mygff.scores[i]<< "\t"<<mygff.strnds[i]<< "\t"<<mygff.frames[i]<<"\t" << mygff.atts[i] << std::endl;
        }

}//end printgffs


/* simple print bed function */
void print_bed(gffs mygff, std::string feature){
if(feature == "ncRNA" || feature == "tRNA"){
    for(size_t i=0; i<mygff.chr.size(); i++){
      if(mygff.ftyps[i] == feature){
        std::cout << mygff.chr[i] <<"\t"<< mygff.starts[i]-1 <<"\t"<< mygff.stops[i] << "\t" << mygff.atts[i] << "\t" << mygff.scores[i] << "\t" << mygff.strnds[i] << std::endl;
      }
    }
  } else if(feature == "gene" || feature == "exon" || feature == "mRNA" || feature == "CDS"){
    for(size_t i=0; i<mygff.chr.size(); i++){
      if(mygff.ftyps[i] == feature){
        std::cout << mygff.chr[i] <<"\t"<< mygff.starts[i]-1 <<"\t"<< mygff.stops[i] << "\t" << mygff.geneids[i] << "\t" << mygff.scores[i] << "\t" << mygff.strnds[i] << std::endl;
      }
    }
  } else if(feature == "tss"){
    for(size_t i=0; i<mygff.chr.size(); i++){
      if(mygff.ftyps[i] == "gene"){
        if(mygff.strnds[i]=="+"){
        std::cout << mygff.chr[i] <<"\t"<< mygff.starts[i]-1 <<"\t"<< mygff.starts[i] << "\t" << mygff.geneids[i] << "\t" << mygff.scores[i] << "\t" << mygff.strnds[i] << std::endl;
      } else if(mygff.strnds[i]=="-"){
        std::cout << mygff.chr[i] <<"\t"<< mygff.stops[i]-1 <<"\t"<< mygff.stops[i] << "\t" << mygff.geneids[i] << "\t" << mygff.scores[i] << "\t" << mygff.strnds[i] << std::endl;
      }
      }
    }
  } else {
    for(size_t i=0; i<mygff.chr.size(); i++){
      if(mygff.ftyps[i] == feature){
        std::cout << mygff.chr[i] <<"\t"<< mygff.starts[i]-1 <<"\t"<< mygff.stops[i] << "\t" << mygff.geneids[i] << "\t" << mygff.scores[i] << "\t" << mygff.strnds[i] << std::endl;
      }
    }
  }
}//end print_bed




/* print bed window function */
void print_bed_window(gffs mygff, std::string feature, int wind){
if(feature == "ncRNA" || feature == "tRNA"){
    for(size_t i=0; i<mygff.chr.size(); i++){
      if(mygff.ftyps[i] == feature){
          std::cout << mygff.chr[i] <<"\t"<< (mygff.starts[i]-wind-1 >=0 ? mygff.starts[i]-wind-1 : 0) <<"\t"<< mygff.stops[i]+wind << "\t" << mygff.atts[i] << "\t" << mygff.scores[i] << "\t" << mygff.strnds[i] << std::endl;
      }
    }
  } else if(feature == "gene" || feature == "mRNA" ){
    for(size_t i=0; i<mygff.chr.size(); i++){
      if(mygff.ftyps[i] == feature){
        std::cout << mygff.chr[i] <<"\t"<< (mygff.starts[i]-wind-1 >=0 ? mygff.starts[i]-wind-1 : 0) <<"\t"<< mygff.stops[i]+wind << "\t" << mygff.geneids[i] << "\t" << mygff.scores[i] << "\t" << mygff.strnds[i] << std::endl;
      }
    }
  } else if(feature == "tss"){
    for(size_t i=0; i<mygff.chr.size(); i++){
      if(mygff.ftyps[i] == "gene"){
        if(mygff.strnds[i]=="+"){
        std::cout << mygff.chr[i] <<"\t"<< (mygff.starts[i]-wind-1 >=0 ? mygff.starts[i]-wind-1 : 0) <<"\t"<< mygff.starts[i]+wind << "\t" << mygff.geneids[i] << "\t" << mygff.scores[i] << "\t" << mygff.strnds[i] << std::endl;
      } else if(mygff.strnds[i]=="-"){
        std::cout << mygff.chr[i] <<"\t"<< (mygff.stops[i]-wind >=0 ? mygff.stops[i]-wind-1 : 0) <<"\t"<< mygff.stops[i]+wind << "\t" << mygff.geneids[i] << "\t" << mygff.scores[i] << "\t" << mygff.strnds[i] << std::endl;
      }
      }
    }
  } else if(feature == "exon" || feature == "CDS"){
    std::cout << "ERROR: Feature `" << feature << "' is not valid with window." << std::endl;
    exit (EXIT_FAILURE);
  } else {
    for(size_t i=0; i<mygff.chr.size(); i++){
      if(mygff.ftyps[i] == feature){
        std::cout << mygff.chr[i] <<"\t"<< (mygff.starts[i]-wind-1 >=0 ? mygff.starts[i]-wind-1 : 0) <<"\t"<< mygff.stops[i]+wind << "\t" << mygff.geneids[i] << "\t" << mygff.scores[i] << "\t" << mygff.strnds[i] << std::endl;
      }
    }
  }
}//end print_bed_window









/* print bed upstream function */
void print_bed_upstream(gffs mygff, std::string feature, int wind){
if(feature == "ncRNA" || feature == "tRNA"){
    for(size_t i=0; i<mygff.chr.size(); i++){
      if(mygff.ftyps[i] == feature){
        if(mygff.strnds[i] == "+"){
          std::cout << mygff.chr[i] <<"\t"<< (mygff.starts[i]-wind-1 >=0 ? mygff.starts[i]-wind-1 : 0) <<"\t"<< mygff.stops[i] << "\t" << mygff.atts[i] << "\t" << mygff.scores[i] << "\t" << mygff.strnds[i] << std::endl;
        } else if(mygff.strnds[i] == "-"){
          std::cout << mygff.chr[i] <<"\t"<< mygff.starts[i]-1 <<"\t"<< mygff.stops[i]+wind << "\t" << mygff.atts[i] << "\t" << mygff.scores[i] << "\t" << mygff.strnds[i] << std::endl;
        }
      }
    }
  } else if(feature == "gene" || feature == "mRNA"){
    for(size_t i=0; i<mygff.chr.size(); i++){
      if(mygff.ftyps[i] == feature){
        if(mygff.strnds[i] == "+"){
          std::cout << mygff.chr[i] <<"\t"<< (mygff.starts[i]-wind-1 >=0 ? mygff.starts[i]-wind-1 : 0) <<"\t"<< mygff.stops[i] << "\t" << mygff.geneids[i] << "\t" << mygff.scores[i] << "\t" << mygff.strnds[i] << std::endl;
        } else if(mygff.strnds[i] == "-"){
          std::cout << mygff.chr[i] <<"\t"<< mygff.starts[i]-1 <<"\t"<< mygff.stops[i]+wind << "\t" << mygff.geneids[i] << "\t" << mygff.scores[i] << "\t" << mygff.strnds[i] << std::endl;
        }
      }
    }
  } else if(feature == "tss"){
    for(size_t i=0; i<mygff.chr.size(); i++){
      if(mygff.ftyps[i] == "gene"){
        if(mygff.strnds[i]=="+"){
        std::cout << mygff.chr[i] <<"\t"<< (mygff.starts[i]-wind-1 >=0 ? mygff.starts[i]-wind-1 : 0) <<"\t"<< mygff.starts[i] << "\t" << mygff.geneids[i] << "\t" << mygff.scores[i] << "\t" << mygff.strnds[i] << std::endl;
      } else if(mygff.strnds[i]=="-"){
        std::cout << mygff.chr[i] <<"\t"<< mygff.stops[i]-1 <<"\t"<< mygff.stops[i]+wind << "\t" << mygff.geneids[i] << "\t" << mygff.scores[i] << "\t" << mygff.strnds[i] << std::endl;
      }
      }
    }
  } else if(feature == "exon" || feature == "CDS"){
    std::cout << "ERROR: Feature `" << feature << "' is not valid with window." << std::endl;
    exit (EXIT_FAILURE);
  } else {
    for(size_t i=0; i<mygff.chr.size(); i++){
      if(mygff.ftyps[i] == feature){
        if(mygff.strnds[i] == "+"){
          std::cout << mygff.chr[i] <<"\t"<< (mygff.starts[i]-wind-1 >=0 ? mygff.starts[i]-wind-1 : 0) <<"\t"<< mygff.stops[i] << "\t" << mygff.geneids[i] << "\t" << mygff.scores[i] << "\t" << mygff.strnds[i] << std::endl;
        } else if(mygff.strnds[i] == "-"){
          std::cout << mygff.chr[i] <<"\t"<< mygff.starts[i]-1 <<"\t"<< mygff.stops[i]+wind << "\t" << mygff.geneids[i] << "\t" << mygff.scores[i] << "\t" << mygff.strnds[i] << std::endl;
        }
      }
    }
  }
}//end print_bed_upstream








/* print bed downstream function */
void print_bed_downstream(gffs mygff, std::string feature, int wind){
if(feature == "ncRNA" || feature == "tRNA"){
    for(size_t i=0; i<mygff.chr.size(); i++){
      if(mygff.ftyps[i] == feature){
        if(mygff.strnds[i] == "+"){
          std::cout << mygff.chr[i] <<"\t"<< mygff.starts[i]-1 <<"\t"<< mygff.stops[i]+wind << "\t" << mygff.atts[i] << "\t" << mygff.scores[i] << "\t" << mygff.strnds[i] << std::endl;
        } else if(mygff.strnds[i] == "-"){
          std::cout << mygff.chr[i] <<"\t"<< (mygff.starts[i]-wind-1 >=0 ? mygff.starts[i]-wind-1 : 0) <<"\t"<< mygff.stops[i] << "\t" << mygff.atts[i] << "\t" << mygff.scores[i] << "\t" << mygff.strnds[i] << std::endl;
        }
      }
    }
  } else if(feature == "gene" || feature == "mRNA"){
    for(size_t i=0; i<mygff.chr.size(); i++){
      if(mygff.ftyps[i] == feature){
        if(mygff.strnds[i] == "+"){
          std::cout << mygff.chr[i] <<"\t"<< mygff.starts[i]-1 <<"\t"<< mygff.stops[i]+wind << "\t" << mygff.geneids[i] << "\t" << mygff.scores[i] << "\t" << mygff.strnds[i] << std::endl;
        } else if(mygff.strnds[i] == "-"){
          std::cout << mygff.chr[i] <<"\t"<< (mygff.starts[i]-wind-1 >=0 ? mygff.starts[i]-wind-1 : 0) <<"\t"<< mygff.stops[i] << "\t" << mygff.geneids[i] << "\t" << mygff.scores[i] << "\t" << mygff.strnds[i] << std::endl;
        }
    }
    }
  } else if(feature == "tss"){
    for(size_t i=0; i<mygff.chr.size(); i++){
      if(mygff.ftyps[i] == "gene"){
        if(mygff.strnds[i]=="+"){
        std::cout << mygff.chr[i] <<"\t"<< mygff.starts[i]-1 <<"\t"<< mygff.starts[i]+wind << "\t" << mygff.geneids[i] << "\t" << mygff.scores[i] << "\t" << mygff.strnds[i] << std::endl;
      } else if(mygff.strnds[i]=="-"){
        std::cout << mygff.chr[i] <<"\t"<< (mygff.stops[i]-wind-1 >=0 ? mygff.stops[i]-wind-1 : 0) <<"\t"<< mygff.stops[i] << "\t" << mygff.geneids[i] << "\t" << mygff.scores[i] << "\t" << mygff.strnds[i] << std::endl;
      }
      }
    }
  } else if(feature == "exon" || feature == "CDS"){
    std::cout << "ERROR: Feature `" << feature << "' is not valid with window." << std::endl;
    exit (EXIT_FAILURE);
  } else {
    for(size_t i=0; i<mygff.chr.size(); i++){
      if(mygff.ftyps[i] == feature){
        if(mygff.strnds[i] == "+"){
          std::cout << mygff.chr[i] <<"\t"<< mygff.starts[i]-1 <<"\t"<< mygff.stops[i]+wind << "\t" << mygff.geneids[i] << "\t" << mygff.scores[i] << "\t" << mygff.strnds[i] << std::endl;
        } else if(mygff.strnds[i] == "-"){
          std::cout << mygff.chr[i] <<"\t"<< (mygff.starts[i]-wind-1 >=0 ? mygff.starts[i]-wind-1 : 0) <<"\t"<< mygff.stops[i] << "\t" << mygff.geneids[i] << "\t" << mygff.scores[i] << "\t" << mygff.strnds[i] << std::endl;
        }
    }
    }
  }
}//end print_bed_downstream









void runner(const char* pfile, gfftobed_opts cmd_opts){

  if(cmd_opts.gtf_flag == 0){
            if(cmd_opts.w_flag == 0){
            print_bed(gffread(pfile, cmd_opts.token), cmd_opts.feat);
          } else if(cmd_opts.w_flag ==1){
            print_bed_window(gffread(pfile, cmd_opts.token),cmd_opts.feat,cmd_opts.interval);
          } else if(cmd_opts.w_flag ==2){
            print_bed_upstream(gffread(pfile, cmd_opts.token),cmd_opts.feat,cmd_opts.interval);
          } else if(cmd_opts.w_flag ==3){
            print_bed_downstream(gffread(pfile, cmd_opts.token),cmd_opts.feat,cmd_opts.interval);
          }
    } else if (cmd_opts.gtf_flag ==1){
            if(cmd_opts.w_flag == 0){
            print_bed(gtfread(pfile, cmd_opts.token), cmd_opts.feat);
          } else if(cmd_opts.w_flag ==1){
            print_bed_window(gtfread(pfile, cmd_opts.token),cmd_opts.feat,cmd_opts.interval);
          } else if(cmd_opts.w_flag ==2){
            print_bed_upstream(gtfread(pfile, cmd_opts.token),cmd_opts.feat,cmd_opts.interval);
          } else if(cmd_opts.w_flag ==3){
            print_bed_downstream(gtfread(pfile, cmd_opts.token),cmd_opts.feat,cmd_opts.interval);
          }
    }

}//end runner
