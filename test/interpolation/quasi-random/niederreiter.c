#include <stdio.h>
#include <gsl/gsl_qrng.h>

int main (void)
{
  int i;
  gsl_qrng * q = gsl_qrng_alloc (gsl_qrng_niederreiter_2, 1);
  FILE * f = fopen("niederreiter_1d.txt", "wt");
  fprintf (f, "%.18f\n", 1.0);
  for (i = 0; i < 65536; i++) {
    double v[1];
    gsl_qrng_get (q, v);
    fprintf (f, "%.18f\n", v[0]);
  }
  fclose(f);
  gsl_qrng_free (q);

  q = gsl_qrng_alloc (gsl_qrng_niederreiter_2, 2);
  f = fopen("niederreiter_2d.txt", "wt");
  fprintf (f, "%.18f %.18f\n", 1.0, 1.0);
  fprintf (f, "%.18f %.18f\n", 0.0, 1.0);
  fprintf (f, "%.18f %.18f\n", 1.0, 0.0);
  for (i = 0; i < 65536; i++) {
    double v[2];
    gsl_qrng_get (q, v);
    fprintf (f, "%.18f %.18f\n", v[0], v[1]);
  }
  fclose(f);
  gsl_qrng_free (q);

  return 0;
}
