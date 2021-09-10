#ifndef __CHR_H__
#define __CHR_H__

#include <math.h>
#include <glib.h>


typedef struct {
  int id;
  char *name;
  long len;
  char *assembly;

  // DAVID: partition coordinates 
  long p_start;
  long p_end;
} Chromosome;


void chr_array_free(Chromosome *chr, int n_chr);
Chromosome *chr_copy(const Chromosome *chr);
void chr_free(Chromosome *chr);
long long_min(long a, long b);

#endif
