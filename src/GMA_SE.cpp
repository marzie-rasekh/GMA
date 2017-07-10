#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <string.h>
#include <math.h>
#include <time.h> 
#include <errno.h>        /* errno */
#include <sstream>
#include <iomanip>

using namespace std;

/* 
 * in case of error explain how to run
 */
void explain () {
  // explain how to run
  fprintf(stderr, "simulates reads of fixes length for GMA and writes in fq files of fixed size, \ninput should be 1 chrom only.\n" );
  fprintf(stderr, "Usage: GMA_SE genome.fasta [options]\n\n" );
  fprintf(stderr, "Options : \n" );
  
  fprintf(stderr, "-SNP_RATE\tfloat : single nucleotide polymorphism error ratio [0, 1)\n" );
  fprintf(stderr, "-COVERAGE\tint : coverage of each base pair (default 100)\n" );
  fprintf(stderr, "-READ_LENGTH\tint : read length, (default 100)\n" );
  fprintf(stderr, "-FILE_BREAK\tdouble : maximum number of reads in each fastq file, (default 10E12)\n" );
}

/*
 * A simple code that reads a fasta file of one sequence only
 * and simulates reads with uniform coverage and read length
 */
int main (int argc, char *argv[]) {
  
  // set the parameter names 
  const char *SNP_RATE_PAR = "-SNP_RATE";
  const char *READ_LENGTH_PAR = "-READ_LENGTH";
  const char *COVERAGE_PAR = "-COVERAGE";
  const char *MAX_LINE_PAR = "-FILE_BREAK";
  
  // check the input parameters
  if (argc < 2) {
    explain();
    return 1;
  }
  
  const string FASTA_FILE = argv[1];
  const char *ALPHABET = "ACTG";
  
  double SNP_RATE = 0.00;
  int COVERAGE = 1;
  int READ_LENGTH = 100;
  double MAX_LINE = 10e12;
  
  bool error = 0;
  // set parameters
  for (int i = 2; i < argc; i++) {
    if (strcmp(SNP_RATE_PAR, argv[i]) == 0) {
      int r = 0;
      r = sscanf(argv[i+1], "%lf", &SNP_RATE);
      if (r <= 0) {
        // error handling.
        error = 1;
      }
    }
    else if (strcmp(READ_LENGTH_PAR, argv[i]) == 0) {
      int r = 0;
      r = sscanf(argv[i+1], "%df", &READ_LENGTH);
      if (r <= 0) {
        // error handling.
        error = 1;
      }
    }
    else if (strcmp(COVERAGE_PAR, argv[i]) == 0) {
      int r = 0;
      r = sscanf(argv[i+1], "%df", &COVERAGE);
      if (r <= 0) {
        // error handling.
        error = 1;
      }
    }
    else if (strcmp(MAX_LINE_PAR, argv[i]) == 0) {
      int r = 0;
      r = sscanf(argv[i+1], "%lf", &MAX_LINE);
      if (r <= 0) {
        // error handling.
        error = 1;
      }
    }
    else {
      error = 1;
    }
    if (error) {
      explain();
      return 1;
    }
    i++;
  }
  cout << "running: GMA " << FASTA_FILE << " ";
  cout << SNP_RATE_PAR << " " <<  SNP_RATE << " ";
  cout << READ_LENGTH_PAR << " " << READ_LENGTH << " ";
  cout << COVERAGE_PAR << " " << COVERAGE << " ";
  cout << MAX_LINE_PAR << " " << MAX_LINE << endl;

  // ok lets start
  char c, * read;
  int read_pivot, copy, file_count;
  double random_number, insertion_probability, deletion_probability, SNP_probability;
  unsigned long seq_pos;
  unsigned long long line_count;
  ifstream input;
  ofstream output;

  string seq_name, filename, output_base;
  
  // read the fasta file character by character
  read = new char[READ_LENGTH]; // simulated read 
  read_pivot = 0; // a pivot to circle around the read
  file_count = 0; // number of output files
  line_count = MAX_LINE + 1; // to make first output file

  input.open(FASTA_FILE.c_str());
  
  int start = (FASTA_FILE.find_last_of("/\\") > 0 ? (FASTA_FILE.find_last_of("/")+1) : 0);
  // get fasta basename
  output_base = FASTA_FILE.substr(start);
  output_base = output_base.substr(0, output_base.find(".fa"));

  if (input.is_open()) {
    while ( !input.eof() ) {
      input >> c; // read char by char
      if (c == '>') {
        getline(input, seq_name); // new chromosome/sequence
        seq_pos = 0; // reset the position
        while ( !input.eof() && read_pivot < READ_LENGTH - 1 ) {
          input >> c; // burn off READ_LENGTH characters
          read[read_pivot] = c;
          read_pivot ++; // will not overflow cause we said less than READ_LENGTH above
        }
      }
      else {
        // set the next bp
        read[read_pivot] = c;
        seq_pos ++;
        for (copy = 1; copy <= COVERAGE; copy++) {
          // check output
          if (line_count > MAX_LINE) {
            line_count = 0;
            // open a new output
            if (output.is_open()) {
              output.close();
            }
            std::ostringstream sstream;
            sstream << "fastq/" << output_base << "_"  << file_count << ".fq";
            filename = sstream.str();
            output.open(filename.c_str(), ios::out);
            if (!output.is_open()) {
              // error
              fprintf(stderr, "Could not open file to write %s \n", filename.c_str() );
              return 2;
            }
            cout << "Opening new fq file " << filename << endl;
            file_count++;
          }
          // write the seq name 
          output << "@se_" << seq_name << "_";
          output << setfill('0') << setw(10) << seq_pos;
          output << "(" << copy << ")" << endl; // read name is se_[chrom]_[pos]([copy])
          
          // write out the read
          for (int i = read_pivot + 1; i < READ_LENGTH; i++) {
            random_number = ((double)rand() / (double)RAND_MAX);
            if (random_number < SNP_RATE) {
              output << ALPHABET[rand() % 4];
            }
            else {
              output << read[i];
            }
          }
          for (int i = 0; i <= read_pivot; i++) {
            random_number = ((double)rand() / (double)RAND_MAX);
            if (random_number < SNP_RATE) {
              output << ALPHABET[rand() % 4];
            }
            else {
              output << read[i];
            }
          }
          output << endl;
          
          read_pivot = (read_pivot + 1) % READ_LENGTH; // to control overflow
          // print a dummy line for quality, all top qual
          output << "+se_" << seq_name << "_";
          output << setfill('0') << setw(10) << seq_pos;
          output << "(" << copy << ")" << endl;
          for (int i = 0; i < READ_LENGTH; i++) 
          {
            output << '~';
          }
          output << endl;
          line_count ++;
        }
      }
    }
    input.close();
  }
  else {
    fprintf(stderr, "Could not open file to read %s \n", FASTA_FILE.c_str() );
  }
  //cout << "read the chromosome: " << name << position << endl;
  

  delete read;
  read = NULL;

  return 0;
}

