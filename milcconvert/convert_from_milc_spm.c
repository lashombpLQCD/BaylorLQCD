//convert_from_milc_spm.c
//reads in the milc gauge field and creates files for input to main
//lewis program when gfconvflag = 1 (set in cfgspropsmain.f90)
//if gfconvflag = 1, then the main program will convert these files to 
//lewis format. Then rerun with gfconvflag = 0
//This is necessary due to the header info in the milc files
//Travis H. Whyte and James Hetrick
//9/15/18

#include <stdio.h>
#include <unistd.h>
#include <string.h>
#include <stdlib.h>
#include <errno.h>
#include "latparams.h"


//#define MAXFILENAME  424
//#define MILC_VERSION_NUMBER 20103
//#define TIMESTAMP_STRING_LEN 68

void byterevn(int w[], int n);

int main(int argc, char *argv[]) {

   FILE *in;
   FILE *out;
   char fnamein[FILENAME_MAX] = MILCIN;
   char fnameout[FILENAME_MAX];
   int milc_magic_number[1];
   int dim[4];
   int k = 0;
   int i = 0;
   int numprocs = NUMPROCS;
   int latdims[] = { LATDIMS };
   int sizelat = 18*4*latdims[0]*latdims[1]*latdims[2]*latdims[3];
   float *buffer;
   char startfile[MAXFILENAME];
   int total_sweeps[1];
   int byterevflag;
   unsigned int check[2];
   int gaction[] = {GACTION};
   float beta[] = {BETA};
   float alat[] = {ALAT};
   float tad[] = {TAD};
   int icfgsave[1] = {NCONFIG};
   int myidin[1];
   char errorString[512];

/* get header info */
   if ( (in = fopen(fnamein, "rb")) == NULL ) {
      fprintf(stderr, "Unable to open input file %s, ", fnamein);
      fprintf(stderr, "error %d.\n", errno);
      exit(1);
   }

   fread(milc_magic_number, sizeof(int), 1, in);
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
   fread(dim, sizeof(int), 4, in);
   if (byterevflag == 1) {
      byterevn(dim, 4);
   }

   fread(startfile, sizeof(char), TIMESTAMP_STRING_LEN, in);

   fread(check, sizeof(unsigned int), 2, in);

   buffer = (float *)malloc(sizeof(float)*sizelat);

// READ IN THE GAUGE FIELD FLOATS
   fread(buffer, sizeof(float), sizelat, in);

  while ( i < numprocs) { 
   myidin[0] = i;
   

   sprintf(fnameout, "l2464f211b600m0102m0509m635a.000%i%03d",nconfig, i);


   out = fopen(fnameout, "wb");
   if ((out = fopen(fnameout, "wb")) == NULL) {
     printf("Unable to open output file. error %d\n", errno);
     perror(errorString);
     printf("%s\n", errorString);
     exit(1);
   }


   fwrite(gaction, sizeof(int), 4, out);
   fwrite(beta, sizeof(float), 1, out);
   fwrite(alat, sizeof(float), 4, out);
   fwrite(tad, sizeof(float), 4, out);
   fwrite(icfgsave, sizeof(int), 1, out);
   fwrite(myidin, sizeof(int), 1, out);
   if ( i == 0 ) {
    // WRITE GAUGEFIELD FLOATS ONLY TO THE FIRST FILE
    fwrite(buffer, sizeof(float), sizelat, out);
   } 


    i++;
  } //end while

    free(buffer);
    fclose(out);
    fclose(in);

    return 0;

}

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
}
