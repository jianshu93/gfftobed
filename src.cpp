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
}//end gffread7




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
        std::cout << mygff.chr[i] <<"\t"<< mygff.starts[i] <<"\t"<< mygff.stops[i] << "\t" << mygff.atts[i] << "\t" << mygff.scores[i] << "\t" << mygff.strnds[i] << std::endl;
      }
    }
  } else if(feature == "gene" || feature == "exon" || feature == "mRNA" || feature == "CDS"){
    for(size_t i=0; i<mygff.chr.size(); i++){
      if(mygff.ftyps[i] == feature){
        std::cout << mygff.chr[i] <<"\t"<< mygff.starts[i] <<"\t"<< mygff.stops[i] << "\t" << mygff.geneids[i] << "\t" << mygff.scores[i] << "\t" << mygff.strnds[i] << std::endl;
      }
    }
  } else if(feature == "tss"){
    for(size_t i=0; i<mygff.chr.size(); i++){
      if(mygff.ftyps[i] == "gene"){
        if(mygff.strnds[i]=="+"){
        std::cout << mygff.chr[i] <<"\t"<< mygff.starts[i] <<"\t"<< mygff.starts[i]+1 << "\t" << mygff.geneids[i] << "\t" << mygff.scores[i] << "\t" << mygff.strnds[i] << std::endl;
      } else if(mygff.strnds[i]=="-"){
        std::cout << mygff.chr[i] <<"\t"<< mygff.stops[i] <<"\t"<< mygff.stops[i]+1 << "\t" << mygff.geneids[i] << "\t" << mygff.scores[i] << "\t" << mygff.strnds[i] << std::endl;
      }
      }
    }
  } else {
    std::cout << "ERROR: `" << feature << "' is not a valid feature." << std::endl;
    exit (EXIT_FAILURE);
  }
}//end print_bed




/* print bed window function */
void print_bed_window(gffs mygff, std::string feature, int wind){
if(feature == "ncRNA" || feature == "tRNA"){
    for(size_t i=0; i<mygff.chr.size(); i++){
      if(mygff.ftyps[i] == feature){
          std::cout << mygff.chr[i] <<"\t"<< (mygff.starts[i]-wind >=0 ? mygff.starts[i]-wind : 0) <<"\t"<< mygff.stops[i]+wind << "\t" << mygff.atts[i] << "\t" << mygff.scores[i] << "\t" << mygff.strnds[i] << std::endl;
      }
    }
  } else if(feature == "gene" || feature == "mRNA" ){
    for(size_t i=0; i<mygff.chr.size(); i++){
      if(mygff.ftyps[i] == feature){
        std::cout << mygff.chr[i] <<"\t"<< (mygff.starts[i]-wind >=0 ? mygff.starts[i]-wind : 0) <<"\t"<< mygff.stops[i]+wind << "\t" << mygff.geneids[i] << "\t" << mygff.scores[i] << "\t" << mygff.strnds[i] << std::endl;
      }
    }
  } else if(feature == "tss"){
    for(size_t i=0; i<mygff.chr.size(); i++){
      if(mygff.ftyps[i] == "gene"){
        if(mygff.strnds[i]=="+"){
        std::cout << mygff.chr[i] <<"\t"<< (mygff.starts[i]-wind >=0 ? mygff.starts[i]-wind : 0) <<"\t"<< mygff.starts[i]+wind << "\t" << mygff.geneids[i] << "\t" << mygff.scores[i] << "\t" << mygff.strnds[i] << std::endl;
      } else if(mygff.strnds[i]=="-"){
        std::cout << mygff.chr[i] <<"\t"<< (mygff.stops[i]-wind >=0 ? mygff.stops[i]-wind : 0) <<"\t"<< mygff.stops[i]+wind << "\t" << mygff.geneids[i] << "\t" << mygff.scores[i] << "\t" << mygff.strnds[i] << std::endl;
      }
      }
    }
  } else if(feature == "exon" || feature == "CDS"){
    std::cout << "ERROR: Feature `" << feature << "' is not valid with window." << std::endl;
    exit (EXIT_FAILURE);
  } else {
    std::cout << "ERROR: `" << feature << "' is not a valid feature." << std::endl;
    exit (EXIT_FAILURE);
  }
}//end print_bed_window









/* print bed upstream function */
void print_bed_upstream(gffs mygff, std::string feature, int wind){
if(feature == "ncRNA" || feature == "tRNA"){
    for(size_t i=0; i<mygff.chr.size(); i++){
      if(mygff.ftyps[i] == feature){
        if(mygff.strnds[i] == "+"){
          std::cout << mygff.chr[i] <<"\t"<< (mygff.starts[i]-wind >=0 ? mygff.starts[i]-wind : 0) <<"\t"<< mygff.stops[i] << "\t" << mygff.atts[i] << "\t" << mygff.scores[i] << "\t" << mygff.strnds[i] << std::endl;
        } else if(mygff.strnds[i] == "-"){
          std::cout << mygff.chr[i] <<"\t"<< mygff.starts[i] <<"\t"<< mygff.stops[i]+wind << "\t" << mygff.atts[i] << "\t" << mygff.scores[i] << "\t" << mygff.strnds[i] << std::endl;
        }
      }
    }
  } else if(feature == "gene" || feature == "mRNA"){
    for(size_t i=0; i<mygff.chr.size(); i++){
      if(mygff.ftyps[i] == feature){
        if(mygff.strnds[i] == "+"){
          std::cout << mygff.chr[i] <<"\t"<< (mygff.starts[i]-wind >=0 ? mygff.starts[i]-wind : 0) <<"\t"<< mygff.stops[i] << "\t" << mygff.geneids[i] << "\t" << mygff.scores[i] << "\t" << mygff.strnds[i] << std::endl;
        } else if(mygff.strnds[i] == "-"){
          std::cout << mygff.chr[i] <<"\t"<< mygff.starts[i] <<"\t"<< mygff.stops[i]+wind << "\t" << mygff.geneids[i] << "\t" << mygff.scores[i] << "\t" << mygff.strnds[i] << std::endl;
        }
      }
    }
  } else if(feature == "tss"){
    for(size_t i=0; i<mygff.chr.size(); i++){
      if(mygff.ftyps[i] == "gene"){
        if(mygff.strnds[i]=="+"){
        std::cout << mygff.chr[i] <<"\t"<< (mygff.starts[i]-wind >=0 ? mygff.starts[i]-wind : 0) <<"\t"<< mygff.starts[i] << "\t" << mygff.geneids[i] << "\t" << mygff.scores[i] << "\t" << mygff.strnds[i] << std::endl;
      } else if(mygff.strnds[i]=="-"){
        std::cout << mygff.chr[i] <<"\t"<< mygff.stops[i] <<"\t"<< mygff.stops[i]+wind << "\t" << mygff.geneids[i] << "\t" << mygff.scores[i] << "\t" << mygff.strnds[i] << std::endl;
      }
      }
    }
  } else if(feature == "exon" || feature == "CDS"){
    std::cout << "ERROR: Feature `" << feature << "' is not valid with window." << std::endl;
    exit (EXIT_FAILURE);
  } else {
    std::cout << "ERROR: `" << feature << "' is not a valid feature." << std::endl;
    exit (EXIT_FAILURE);
  }
}//end print_bed_upstream








/* print bed downstream function */
void print_bed_downstream(gffs mygff, std::string feature, int wind){
if(feature == "ncRNA" || feature == "tRNA"){
    for(size_t i=0; i<mygff.chr.size(); i++){
      if(mygff.ftyps[i] == feature){
        if(mygff.strnds[i] == "+"){
          std::cout << mygff.chr[i] <<"\t"<< mygff.starts[i] <<"\t"<< mygff.stops[i]+wind << "\t" << mygff.atts[i] << "\t" << mygff.scores[i] << "\t" << mygff.strnds[i] << std::endl;
        } else if(mygff.strnds[i] == "-"){
          std::cout << mygff.chr[i] <<"\t"<< (mygff.starts[i]-wind >=0 ? mygff.starts[i]-wind : 0) <<"\t"<< mygff.stops[i] << "\t" << mygff.atts[i] << "\t" << mygff.scores[i] << "\t" << mygff.strnds[i] << std::endl;
        }
      }
    }
  } else if(feature == "gene" || feature == "mRNA"){
    for(size_t i=0; i<mygff.chr.size(); i++){
      if(mygff.ftyps[i] == feature){
        if(mygff.strnds[i] == "+"){
          std::cout << mygff.chr[i] <<"\t"<< mygff.starts[i] <<"\t"<< mygff.stops[i]+wind << "\t" << mygff.geneids[i] << "\t" << mygff.scores[i] << "\t" << mygff.strnds[i] << std::endl;
        } else if(mygff.strnds[i] == "-"){
          std::cout << mygff.chr[i] <<"\t"<< (mygff.starts[i]-wind >=0 ? mygff.starts[i]-wind : 0) <<"\t"<< mygff.stops[i] << "\t" << mygff.geneids[i] << "\t" << mygff.scores[i] << "\t" << mygff.strnds[i] << std::endl;
        }
    }
    }
  } else if(feature == "tss"){
    for(size_t i=0; i<mygff.chr.size(); i++){
      if(mygff.ftyps[i] == "gene"){
        if(mygff.strnds[i]=="+"){
        std::cout << mygff.chr[i] <<"\t"<< mygff.starts[i] <<"\t"<< mygff.starts[i]+wind << "\t" << mygff.geneids[i] << "\t" << mygff.scores[i] << "\t" << mygff.strnds[i] << std::endl;
      } else if(mygff.strnds[i]=="-"){
        std::cout << mygff.chr[i] <<"\t"<< (mygff.stops[i]-wind >=0 ? mygff.stops[i]-wind : 0) <<"\t"<< mygff.stops[i] << "\t" << mygff.geneids[i] << "\t" << mygff.scores[i] << "\t" << mygff.strnds[i] << std::endl;
      }
      }
    }
  } else if(feature == "exon" || feature == "CDS"){
    std::cout << "ERROR: Feature `" << feature << "' is not valid with window." << std::endl;
    exit (EXIT_FAILURE);
  } else {
    std::cout << "ERROR: `" << feature << "' is not a valid feature." << std::endl;
    exit (EXIT_FAILURE);
  }
}//end print_bed_downstream







/*runner function*/
void runner(const char* pfile, std::string ff, int wind, int wflag, std::string tok){

  if(wflag == 0){
  print_bed(gffread(pfile, tok), ff);
} else if(wflag ==1){
  print_bed_window(gffread(pfile, tok),ff,wind);
} else if(wflag ==2){
  print_bed_upstream(gffread(pfile, tok),ff,wind);
} else if(wflag ==3){
  print_bed_downstream(gffread(pfile, tok),ff,wind);
}

}//end runner
