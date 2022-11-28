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
   char fnamein[FILENAME_MAX] = "/data/whytet/qqcd/convert_from_milc_spm/l2464f211b600m0102m0509m635a.1000";
   char fnameout[FILENAME_MAX];
   int milc_magic_number[1];
   int dim[4];
   //int t, z, y, x, mu, icri, isite, ibl, ieo, err;
   //int tinc, zinc, yinc, xinc;
   //int tmod, zmod, ymod, xmod;
   //int tid, zid, yid, xid, sumid;
   int k = 0;
   int i = 0;
   int numprocs = NUMPROCS;
   int latdims[] = { LATDIMS }
   int sizelat = 18*4*latdims[0]*latdims[1]*latdims[2]*latdims[3]
   char startfile[MAXFILENAME];
   int total_sweeps[1];
   int byterevflag;
   //float u[442368];
   unsigned int check[2];
   //float myu[6];
   int gaction[4];
   float beta[1];
   float alat[4];
   float tad[4];
   //long int start;
   //long int startelem;
   int icfgsave = NCONFIG;
   int myidin[1];
   //int gactionout[4];
   //int icfgsaveout[1];
   //int myidout[1];
   //float tadout[4];
   //float betaout[1];
   //float alatout[4];
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

   fread(buffer, sizeof(float), sizelat, in);

   /*int nt = dim[3];
   int nz = dim[2];
   int ny = dim[1];
   int nx = dim[0];

   long int offset = (442368*sizeof(float));
   float buffer[442368];

   printf("before main for\n");

   //for (i = 0; i < 144; i++); {
    while (i < 144) {
    start = (offset * i);

    fseek(in, (start + i*sizeof(float)), SEEK_CUR); // THIS PART NEEDS WORK!!!! TW 3/12/18
    fread(buffer, sizeof(float), 442368, in);

    printf("i is %i, and my offset is %ld\n", i, (start + i*sizeof(float)));
    float gf[18][192][4][2][16];

    startelem = ((start + i*sizeof(float)) / sizeof(float));
    xinc = (startelem/72);
    yinc = (startelem/(72*nx));
    zinc = (startelem/(72*nx*ny));
    tinc = (startelem/(72*nx*ny*nz));
    tmod = (tinc % nt);
    zmod = (zinc % nz);
    ymod = (yinc % ny);
    xmod = (xinc % nx);

    tid = (int)tmod;
    zid = (int)zmod;
    yid = (int)ymod;
    xid = (int)xmod;
    sumid = (tid + zid + yid + xid);

         for (ibl = 0; ibl < 16; ibl++) {
            for (mu = 0; mu < 4; mu++) {
              for (isite = 0; isite < 192; isite++) {
               for (icri = 0; icri < 18; icri++) {
               if ((sumid % 2) == 0) {
                  gf[ibl][0][mu][isite][icri] = buffer[k];
                }
               if ((sumid % 2) != 0) {
                  gf[ibl][1][mu][isite][icri] = buffer[k];
                  }
                   k++;

                if ( ((startelem + k) % (72)) == 0 ) {
                      xid++;
                }
                if ( xid > nx) {
                      yid++;
                      xid = 0;
                }

                if ( yid > ny) {
                      zid++;
                      yid = 0;
                }

                if ( zid > nz) {
                      tid++;
                      zid = 0;
                }

                sumid = (tid + zid + yid + xid);
              }
            }
         }
     }

     k = 0;
    
     printf("after gf for\n");     

     for (icri = 0; icri < 18; icri++) {
      for (ieo = 0; ieo < 2; ieo++) {
       for (isite = 0; isite < 192; isite++) {
        for (mu = 0; mu < 4; mu++) {
          for (ibl = 0; ibl < 16; ibl++) {
             gf[icri][isite][mu][ieo][ibl] =  gf[ibl][ieo][mu][isite][icri];
         }
       }
     }
    }
  } 

   printf("after gf transpose\n"); */

   gaction[0]= 1;
   gaction[1]= 1;
   gaction[2]= 1;
   gaction[3]= 1;

   tad[0] = 0.86372;
   tad[1] = 0.86372;
   tad[2] = 0.86372;
   tad[3] = 0.86372;

   beta[0] = 6.000000;
   alat[0] = 1.000000;
   alat[1] = 1.000000;
   alat[2] = 1.000000;
   alat[3] = 1.000000;
   //icfgsave[0] = 1;
   
  while ( i < numprocs) { 
   myidin[0] = i;
   
   printf("before sprintf\n");

   sprintf(fnameout, "l2464f211b600m0102m0509m635a.0001%03d", i);

   printf("after sprintf\n");

   out = fopen(fnameout, "wb");
   if ((out = fopen(fnameout, "wb")) == NULL) {
     printf("Unable to open output file. error %d\n", errno);
     perror(errorString);
     printf("%s\n", errorString);
     exit(1);
   }

   printf("after output file open\n");

   fwrite(gaction, sizeof(int), 4, out);
   fwrite(beta, sizeof(float), 1, out);
   fwrite(alat, sizeof(float), 4, out);
   fwrite(tad, sizeof(float), 4, out);
   fwrite(icfgsave, sizeof(int), 1, out);
   fwrite(myidin, sizeof(int), 1, out);
   if ( i == 0 ) {
   fwrite(buffer, sizeof(float), sizelat, out);
   } 

   printf("after write\n");

    if ( (fwrite(buffer, sizeof(float), sizelat, out)) != sizelat ) {
      fprintf(stderr, "unable to write to file %s, ", out);
      fprintf(stderr, "error %d.\n", errno);
      exit(1);
    }
  
    //free(gf);

    //fclose(out);

    //printf("value of i at end = %i\n", i);
    i++;
    }
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
