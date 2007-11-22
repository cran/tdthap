#include <math.h>
#include <stdlib.h>

static int nloci;
static double *spacing = (double *) 0;
static int focus;
static double power;

/* S/R callable */

void set_tdt_similarity(long *nl, double *sp, long *foc, double *pow) {
  int i;
  if (spacing) free(spacing);
  nloci = *nl;
  focus = *foc;
  power = *pow;
  spacing = (double *) calloc(nloci+1, sizeof(double));
  for (i=0; i<=nloci; i++) spacing[i] = sp[i];
}

void get_tdt_similarity(long *nl, double *sp, long *foc, double *pow) {
  int i;
  *nl = nloci;
  *foc = focus;
  *pow = power;
  for (i=0; i<=nloci; i++) sp[i] = spacing[i];
}

double tdt_similarity(int *a, int *b) {
  double s;
  int i;
  /* Haplotypes must have focal locus in common */
  if (a[focus-1] != b[focus-1]) return(0.0); 
  /* Find beginning of common area */
  for (i=focus-2; i>=0 && a[i]==b[i]; i--) {}
  i++;
  s = spacing[i++]/2.0; /* off-end correction */
  /* Find end of common area */
  for (; i<nloci && a[i]==b[i]; i++) s += spacing[i];
  s += spacing[i]/2.0; /* off-end correction */
  if (power!=1.0) s = pow(s, power);
  return(s);
}
