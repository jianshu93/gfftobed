#include <iostream>
#include <unistd.h>
#include <getopt.h>
#include <string>
#include "util.hpp"

void print_usage(){

  printf("gfftobed by Jacob Bierstedt (jacob.bierstedt@und.edu)\n\n");
  printf("USAGE:\n\tgfftobed [options] <input_file_GFF3>\n\n");
  printf("Extracts genomic coordinates of features from GFF3\n");

  printf("-g/--gene\t\t\textract gene features in bed format\n");
  printf("-e/--exons\t\t\textract exon features in bed format\n");
  printf("-c/--cds\t\t\textract CDS features in bed format\n");
  printf("-m/--mrna\t\t\textract mRNA features in bed format\n");
  printf("-t/--tss\t\t\textract tss features in bed format\n");
  printf("-f/--feature <feat>\t\textract custom features (i.e. ncRNA)\n\n");

  printf("Extracts genomic coordinates of features with window around feature\n");
  printf("-w/--window <int>\t\tadd <int> basepairs upstream and downstream of feature\n");
  printf("-u/--upstream <int>\t\tadd <int> basepairs upstream/5' of feature\n");
  printf("-d/--downstream <int>\t\tadd <int> basepairs downstream/3' of feature\n\n");

  printf("-a/--attribute <string>\t\tSpecify attribute for name column ('note' by default)\n");
  printf("-G/--GTF\t\t\tInput file is in GTF format\n");


  printf("-h/--help\t\t\tPrint this help message\n\n");

  exit (EXIT_FAILURE);
}




gfftobed_opts setopts(int argc, char **argv){

  int c;
  gfftobed_opts opts;
  opts.w_flag =0;
  opts.gtf_flag=0;

    std::string feat_extract;
    std::string ttoken;
    int window,wflag=0,aflag=0, ftype_flag=0;

  while (1)
    {
      static struct option long_options[] =
        {
          {"help",       no_argument,             0,   'h'},
          {"gene",       no_argument,             0,   'g'},
          {"exon",       no_argument,             0,   'e'},
          {"cds",        no_argument,             0,   'c'},
          {"mrna",       no_argument,             0,   'm'},
          {"tss",        no_argument,             0,   't'},
          {"feature",    required_argument,             0,   'f'},
          {"window",        required_argument,             0,   'w'},
          {"upstream",        required_argument,             0,   'u'},
          {"downstream",        required_argument,             0,   'd'},
          {"attribute",        required_argument,             0,   'a'},
          {"GTF",        required_argument,             0,   'G'},
          {0, 0, 0, 0}
        };

      int option_index = 0;

      c = getopt_long (argc, argv, "hgecmtf:w:u:d:a:G",
                       long_options, &option_index);


      if (c == -1)
        break;

      switch (c)
        {
        case 'h':
          print_usage();
          break;
        case 'g':
          feat_extract = "gene";
          break;
        case 'e':
          feat_extract = "exon";
          break;
        case 'c':
          feat_extract = "CDS";
          break;
        case 'm':
          feat_extract = "mRNA";
          break;
        case 't':
          feat_extract = "tss";
          break;
        case 'f':
          feat_extract = optarg;
          break;
        case 'w':
          wflag =1;
          window = atoi(optarg);
          break;
        case 'u':
          wflag =2;
          window = atoi(optarg);
          break;
        case 'd':
          wflag = 3;
          window = atoi(optarg);
          break;
        case 'a':
          aflag =1;
          ttoken = optarg;
          break;
        case 'G':
          ftype_flag = 1;
          break;

        case '?':
          break;

        default:
          abort ();
        }
    }





    opts.feat = feat_extract;
    opts.w_flag = wflag;
    opts.gtf_flag = ftype_flag;

    if(wflag==1 || wflag ==2 || wflag ==3){
    opts.interval = window;
  } else {opts.interval = 0;}

    if(aflag==1){
      opts.token=ttoken;
    } else {opts.token = "note";}




    return opts;
} 
