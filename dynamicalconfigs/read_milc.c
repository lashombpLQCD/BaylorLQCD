// read_milc.c
//
// Jim Hetrick
// 6 Mar 2017
//
// Read in a MILC v.7 lattice, report header, 
// and print the first row of U(x=0,0,0,0)[0]
//
///////////////////////////////////////////////////
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <errno.h>

//#define PRECISION 1
//#include "../include/precision.h"

//#include "../include/int32type.h"

#define MAXFILENAME  256   /* ASCII string length for all file names */

#define MILC_VERSION_NUMBER 20103  /* hex 00010203 */

/* This "should" be 80! but the checksums were added in v.7, 
   taking some away from the timestamp string  */
#define TIMESTAMP_STRING_LEN 68

void byterevn(int w[], int n);

int main(int argc, char *argv[]) {
   FILE *fp;
   char fname[FILENAME_MAX];
   int milc_magic_number[1];
   int dim[4];
   char startfile[MAXFILENAME];
   int total_sweeps[1];
   int byterevflag;
   float u[6];
   unsigned int check[2];


   /* get input filename */
   if (argc == 2) {
      strcpy(fname, argv[1]);
   }
   else {
      printf("Enter lattice filename       ");
      scanf("%s", fname);
   }


   /* open binary MILC lattice file and header information */
   if ( (fp = fopen(fname, "rb")) == NULL ) {
      fprintf(stderr, "Unable to open input file %s, ", fname);
      fprintf(stderr, "error %d.\n", errno);
      exit(1);
   }


   fread(milc_magic_number, sizeof(int), 1, fp);
   if(*milc_magic_number == MILC_VERSION_NUMBER) {
      byterevflag = 0;
   }
   else {
      byterevn(milc_magic_number, 1);
      if(*milc_magic_number == MILC_VERSION_NUMBER) {
         byterevflag = 1;
      } else {
         printf("Unrecognized version number: %d\nTerminating program.\n",
		*milc_magic_number);
         exit(1);
      }
   }
   fread(dim, sizeof(int), 4, fp);
   if (byterevflag == 1) {
      byterevn(dim, 4);
   }

   fread(startfile, sizeof(char), TIMESTAMP_STRING_LEN, fp);

   fread(check, sizeof(unsigned int), 2, fp);

   fread(u, sizeof(float), 6, fp);

   fclose(fp);

   /*
   if (byterevflag == 1)
   {
      byterevn((float *)u, 6);
   }
   */


   /* report lattice file description */
   printf("magic_number  %d\n", *milc_magic_number);
   printf("byterev flag is %d\n", byterevflag);
   printf("nx ny nz nt   %d %d %d %d\n", dim[0], dim[1], dim[2], dim[3]);
   printf("startfile     %s\n", startfile);
   printf("checksums: %lx %lx\n", check[0], check[1]);
   printf("First row of U(x=0,0,0,0)[mu=0]:\n");
   printf("%f\n", u[0]);
   printf("%f\n", u[1]);
   printf("%f\n", u[2]);
   printf("%f\n", u[3]);
   printf("%f\n", u[4]);
   printf("%f\n", u[5]);

   /* terminate program normally */
   return 0;
}

/*----------------------------------------------------------------------*/

/* For doing byte reversal on 32-bit words */
/* From io_lat4.c */

void byterevn(int w[], int n)
{
  register int old,new;
  int j;
  
  for(j=0; j<n; j++)
    {
      old = w[j];
      new = old >> 24 & 0x000000ff;
      new |= old >> 8 & 0x0000ff00;
      new |= old << 8 & 0x00ff0000;
      new |= old << 24 & 0xff000000;
      w[j] = new;
    }
} /* byterevn */
