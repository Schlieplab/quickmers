// FILE: kmc_npy.c
//
// Convert kmc output in text form containing k-mers and frequencies into a binary npy-file 
// with 2-bit encoded k-mers
//
// NOTE: do not compile with -O3. Breaks the output
//
// Author: alexander@schlieplab.org
// 
#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>

// Total header size 128 bytes
// - needs to be multiple of 64
// - needs to be large enough for dict describing array
// We only need two unsigned ints for the 2-bit encoded k-mer and its frequency
#define HEADER_LEN (128 - 6 - 1 - 1 - 2)
typedef struct npy_header_t {
  char magic[6];
  char major;
  char minor;
  uint16_t header_len;
  char dict[HEADER_LEN];
} npy_header_t;


// Pragma needed for proper binary format
// Needs to be adapted to different formats
#pragma pack(push,1)
typedef struct record_t {
  uint32_t kmer;
  uint16_t freq;
} record_t;
#pragma pack(pop)
#define BUFFER_NR_RECORDS 100000


void write_npy_header(int nr_kmer, FILE* outfile) {
  int i;
  npy_header_t npy_header;
  //
  // TODO: Make formats adaptable to k-mer length:
  // u2 suffices for 7-mers, u8 needed for 31-mers. Should be determined automatically or from command line
  // frequencies should be capped at maximal value and #bytes selected on command line
  char header_str[] = "{'descr': [('kmer', '<u4'), ('freq', '<u2')], 'fortran_order': False, 'shape': (%d,), }";
  int written;
    
  npy_header.magic[0] = '\x93';
  npy_header.magic[1] = 'N';
  npy_header.magic[2] = 'U';
  npy_header.magic[3] = 'M';
  npy_header.magic[4] = 'P';
  npy_header.magic[5] = 'Y';
  npy_header.major = '\x01';
  npy_header.minor = '\x00';
  npy_header.header_len = HEADER_LEN;
  written = sprintf(npy_header.dict, header_str, nr_kmer);
  for (i = written; i < HEADER_LEN; i++)
	npy_header.dict[i] = '\x20';
    
  fwrite((void*)&npy_header, sizeof(char), sizeof(npy_header_t),  outfile);
}


int kmer_len(char *line)
{
    int k = 0;
    while ((line[k] != ' ') && (line[k] != '\t')) {
        k++;
    }
    return k;
}

void peek_first_line(FILE *file, char *line)
{   
    long pos = ftell(file);
    fgets(line, 256, file);
    fseek(file, pos, SEEK_SET);
}

int count_nr_lines(FILE *file)
{
    long pos = ftell(file);
    int count = 0;
    char line[256];
    while (fgets(line, sizeof(line), file)) {
        count++;
    }
    fseek(file, pos, SEEK_SET);
    return count;
}

// From https://www.w3resource.com/c-programming-exercises/c-snippets/print-binary-format-in-c-using-printf-alternatives.php
void printBinaryWithPadding(uint64_t num) {
    for (int i = sizeof(uint64_t) * 8 - 1; i >= 0; i--) {
        printf("%d", (num >> i) & 1);
        if (i % 4 == 0) printf(" "); // Group by 4 bits for readability
    }
    printf("\n");
}


uint32_t kmer_to_uint32(char *s, int k)
{
    int i;
	uint64_t low_bits_mask = 0b0000001100000011000000110000001100000011000000110000001100000011; 
    uint32_t kmer = 0;
    uint64_t* u = (uint64_t*)s;

	printf("%.*s\n", k, s);

	printBinaryWithPadding(u[0]);
	printBinaryWithPadding(u[1]);
	
    for (i = 0; i < (k>>3); i++) {
        u[i] = (u[i] >> 1) & low_bits_mask;
        printBinaryWithPadding(u[i]);
        u[i] |= (u[i] >> 6);
		printBinaryWithPadding(u[i]);
        u[i] |= (u[i] >> 12);    
		printBinaryWithPadding(u[i]);
	}
	
    // Assuming k <= 16
    kmer = ((uint32_t)s[13]) << 24 | ((uint32_t)s[9]) << 16 | ((uint32_t)s[5]) << 8 | (uint32_t)s[0];
	printBinaryWithPadding(kmer);
    return kmer;
}

int main(int argc, char *argv[])
{
  int i = 0;
  int k;
  char line[256];
  int nr_lines;
  int total;
  record_t* buffer;
  
  if (argc != 3) {
	printf("Usage: %s <kmc-test-output.txt> <numpy-output.npy>\n", argv[0]);
	return 1;
  }

  FILE *infile = fopen(argv[1], "r");
  if (!infile) {
	printf("Error opening file '%s'\n", argv[1]);
	return 1;
  }

  peek_first_line(infile, line);
  k = kmer_len(line);
  nr_lines = count_nr_lines(infile); // Needed for header in output file
  printf("# Input file %s contains %d %d-mers\n", argv[1], nr_lines, k); 
  
  buffer = (record_t*)malloc(BUFFER_NR_RECORDS * sizeof(record_t));
  
  FILE *outfile = fopen(argv[2], "wb");
  if (!outfile) {
	printf("Error opening file '%s'\n", argv[2]);
	return 1;
  }

  write_npy_header(nr_lines, outfile);

  while (total < nr_lines) {
	for (i = 0; i < BUFFER_NR_RECORDS; i++) {
	  if (fgets(line, sizeof(line), infile)) {
		buffer[i].kmer = kmer_to_uint32(line, k);
		buffer[i].freq = (uint16_t)atoi(line + k + 1); 
	  } else // File ended before buffer was full 
		break;
	}
	// printf("# Writing %d additional records (total %d)\n", i, total);
	fwrite((void*)buffer, sizeof(char), i*sizeof(record_t),  outfile);
	total += i;
  }

  fclose(infile);
  fclose(outfile);
  return 0;
}
