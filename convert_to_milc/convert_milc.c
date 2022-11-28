// program to convert MILC binary gauge files to 
// a version compatible with the Fortran program
// Travis Whyte
// March 11th 2018


#include <stdio.h>
#include <unistd.h>
#include <string.h>
#include <stdlib.h>
#include <errno.h>
#include <mpi.h>

#define MAXFILENAME  424
#define MILC_VERSION_NUMBER 20103
#define TIMESTAMP_STRING_LEN 68


void byterevn(int w[], int n);

int main(int argc, char *argv[]) {
   MPI_File mpiin;
   long int myidoff;
   long int myidlong;
   long int start;
   long int startelem;
   long int header = 97;
   MPI_Status status;
   MPI_File mpiout;
   FILE *in;
   FILE *out;
   char fnamein[FILENAME_MAX] = "l2464f211b600m0102m0509m635a.1000";
   char fnameout[FILENAME_MAX];
   int milc_magic_number[1];
   int dim[4];
   int ierr, myid, numprocs, err;
   int i, t, z, y, x, mu, icri, isite, ibl, ieo;
   int tinc, zinc, yinc, xinc;
   int tmod, zmod, ymod, xmod;
   int tid, zid, yid, xid, sumid;
   int k = 0;
   char startfile[MAXFILENAME];
   int total_sweeps[1];
   int byterevflag;
   float u[442368];
   unsigned int check[2];
   float myu[6];
   int gaction[4];
   float beta[1];
   float alat[4];
   float tad[4];
   int icfgsave[1];
   int myidin[1];
   int gactionout[4];
   int icfgsaveout[1];
   int myidout[1];
   float tadout[4];
   float betaout[1];
   float alatout[4];

   //int gaction[1][1][1][1];
   //float beta[1];
   //float alat[1];
   //float tad[1][1][1][1];
   //int icfgsave[1];
   //int myidout[1];
   

   
   ierr = MPI_Init(&argc, &argv);
   ierr = MPI_Comm_rank(MPI_COMM_WORLD, &myid);
   ierr = MPI_Comm_size(MPI_COMM_WORLD, &numprocs);

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
 //  fclose(in);

  // fclose(in); 
   //end get header info

  // ierr = MPI_File_open(MPI_COMM_WORLD,"l2464f211b600m0102m0509m635a.1000", MPI_MODE_RDONLY, MPI_INFO_NULL, &mpiin);
   /* if ( (MPI_File_open(MPI_COMM_SELF,"l2464f211b600m0102m0509m635a.1000", MPI_MODE_RDONLY, MPI_INFO_NULL, &mpiin)) != MPI_SUCCESS) {
       printf("could not open input file\n");
       printf("error %d.\n", errno);
       MPI_Abort(MPI_COMM_SELF, err);
       printf("MPI error %d.\n", err);
    } */
   /*MPI_File_read(mpiin, milc_magic_number, 1, MPI_INT, &status);
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

   MPI_File_read(mpiin, dim, 4, MPI_INT, &status);
   if (byterevflag == 1) {
      byterevn(dim, 4);
   }

   MPI_File_read(mpiin, startfile, TIMESTAMP_STRING_LEN, MPI_INT, &status);
   MPI_File_read(mpiin, check, 2, MPI_UNSIGNED, &status);

   if (myid == 0) {
    printf("magic_number  %d\n", *milc_magic_number);
    printf("byterev flag is %d\n", byterevflag);
    printf("nx ny nz nt   %d %d %d %d\n", dim[0], dim[1], dim[2], dim[3]);
    printf("startfile     %s\n", startfile);
    printf("checksums: %lx %lx\n", check[0], check[1]);
   }

  // MPI_File_set_view(mpiin, 96, MPI_FLOAT, MPI_FLOAT, "native", MPI_INFO_NULL); */
  /* in = fopen(fnamein, "rb");
      if ( (in = fopen(fnamein, "rb")) == NULL ) {
      fprintf(stderr, "Unable to open input file %s, ", fnamein);
      fprintf(stderr, "error %d.\n", errno);
      exit(1);
      } */

// for (myid = 0; myid < 7; myid++) {
 /*  in = fopen(fnamein, "rb");
   if ((in = fopen(fnamein, "rb")) == NULL) {
     fprintf(stderr, "Unable to open input file %s, ", fnamein);
     fprintf(stderr, "error %d.\n", errno);
      exit(1);
   }*/
   int nt = dim[3];
   int nz = dim[2];
   int ny = dim[1];
   int nx = dim[0];
   int npt = 8;
   int npz = 3;
   int npy = 3;
   int npx = 2;
   int nps = (npt*npz*npy*npx);
   int ntzyx = nt*nz*ny*nx;
   int nvhalf = (ntzyx/32/nps);
   //int elmpproc = 442368; //((72*ntzyx)/144); //BIG COMMENT!!! NEED TO DECLARE NUMBER OF PROCESSORS HERE
   long int offset = (442368*sizeof(float));  //I WILL WORK ON USING MALLOC SO THIS IS NOT NEEDED -TW 3/15/18
   
   //float buffer[elmpproc];
   float buffer[442368];
   start = (offset * myid);
 //  printf("myid is %i and my value of start is %i\n", myid, start);
   //myidoff = (start + header);
 //  printf("myid is %i and my value of myidoffset is %i\n", myid, myidoff);
 // if (myid == 0) {
    // fseek(in, HEADER_SIZE, SEEK_SET); //HEADER_SIZE is not long int
 //    fread(buffer, sizeof(float), 442368, in);

   //ierr = MPI_File_read_at(mpiin, HEADER_SIZE, buffer, offset, MPI_FLOAT, &status);
  // } else {
     fseek(in, (start + myid*sizeof(float)), SEEK_CUR); // THIS PART NEEDS WORK!!!! TW 3/12/18
     fread(buffer, sizeof(float), 442368, in);
  // ierr = MPI_File_read_at(mpiin, (HEADER_SIZE + start), buffer, offset, MPI_FLOAT, &status);
 //  } 
     fclose(in);
// }
  // ierr = MPI_File_close(&mpiin);
   printf("myid is %i, and my offset is %ld\n", myid, (start + myid*sizeof(float)));  
   float gf[18][192][4][2][16];
  // double dgf[18][nvhalf][4][2][16];
  // float gf;

   startelem = ((start + myid*sizeof(float)) / sizeof(float));
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

   //printf("myid is %i and my initial indexes are %i, %i, %i, %i\n", myid, tid, zid, yid, xid);

  // while (k < 442368 ) {
   //for (t = 0; t < nt; t++) {  // THIS PART NEEDS WORK, EACH PROCESSOR SEES THE SAME VALUE FOR T,Z,Y,X
    // for (z = 0; z < nz; z++) {  // AND WILL NOT ASSIGN IEO CORRECTLY
     // for (y = 0; y < ny; y++) {
      // for (x = 0; x < nx; x++) {
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
               /* if ( ((startelem + k) % (72*nx*ny*nz)) == 0 ) {
                      tid++;
                }
                if ( tid > nt) {
                      tid = 0;
                }
                if ( ((startelem + k) % (72*nx*ny)) == 0 ) {
                      zid++;
                }
                if ( zid > nz) {
                      zid = 0;
                }
                if ( ((startelem + k) % (72*nx)) == 0 ) {
                      yid++;
                }
                if ( yid > ny) {
                      yid = 0;
                }
                if ( ((startelem + k) % (72)) == 0 ) {
                      xid++;
                }
                if ( xid > nx) {
                      xid = 0;
                } */
                sumid = (tid + zid + yid + xid);
                  }
                }
               // k++;
               }
              }
            // }
           // }
          // }
          //}
        // }
      // } 
    //   k++;
       //printf("on iteration %d",k);

   // }


  /* gaction[0][0][0][0] = 1;
   tad[0][0][0][0] = 0.86372;
   beta[1] = 6.00;
   alat[1] = 0.12;
   icfgsave[1] = 1;
   myidout[1] = myid; */

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

   gaction[0]= 1;
   gaction[1]= 1;
   gaction[2]= 1;
   gaction[3]= 1;

   tad[0] = 0.86372;
   tad[1] = 0.86372;
   tad[2] = 0.86372;
   tad[3] = 0.86372;

   beta[0] = 6.000000;
   alat[0] = 0.120000;
   alat[1] = 0.120000;
   alat[2] = 0.120000;
   alat[3] = 0.120000;
   icfgsave[0] = 1;
   myidin[0] = myid;


 if (myid < 10) {
   sprintf(fnameout, "l2464f211b600m0102m0509m635a.000100%i", myid);
 }
 if (myid >= 10) {
   if  (myid < 100)  {
     sprintf(fnameout, "l2464f211b600m0102m0509m635a.00010%i", myid); //creates a string stored in fnameout
   }  else {
   sprintf(fnameout, "l2464f211b600m0102m0509m635a.0001%i", myid);
   }
 }
   out = fopen(fnameout, "w+b");
   if ((out = fopen(fnameout, "w+b")) == NULL) {
     fprintf(stderr, "Unable to open input file %s, ", fnameout);
     fprintf(stderr, "error %d.\n", errno);
     exit(1);
   }
  /* MPI_File_open(MPI_COMM_SELF, "fnameout", MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &mpiout); // call to fnameout here has to be the problem
     if ( ( MPI_File_open(MPI_COMM_SELF, "fnameout", MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &mpiout)) != MPI_SUCCESS) {
         printf("could not open the file\n");
         MPI_Abort(MPI_COMM_SELF, err);
         printf("MPI error %d.\n", err);
     }
   MPI_File_write(mpiout, gf, 442368, MPI_FLOAT, &status);
 //  MPI_File_read(mpiout, u, 5, MPI_FLOAT, &status);*/
 //
   fwrite(gaction, sizeof(int), 4, out);
   fwrite(beta, sizeof(float), 1, out);
   fwrite(alat, sizeof(float), 4, out);
   fwrite(tad, sizeof(float), 4, out);
   fwrite(icfgsave, sizeof(int), 1, out);
   fwrite(myidin, sizeof(int), 1, out); 
   fwrite(gf, sizeof(float), 442368, out);
   
    if ( (fwrite(gf, sizeof(float), 442368, out)) != 442368 ) {
      fprintf(stderr, "unable to write to file %s, ", out);
      fprintf(stderr, "error %d.\n", errno);
      exit(1);
    }
   
   
  /* fwrite(gaction, sizeof(int), 4, out);
   fwrite(beta, sizeof(float), 1, out);
   fwrite(alat, sizeof(float), 1, out);
   fwrite(tad, sizeof(float), 4, out);
   fwrite(icfgsave, sizeof(int), 1, out);
   fwrite(myidout, sizeof(int), 1, out); */
   fclose(out);
  // MPI_File_close(&mpiout);

  if (myid == 0) {
  /* MPI_File_open(MPI_COMM_SELF, "l2464f211b600m0102m0509m635a.00010", MPI_MODE_RDONLY, MPI_INFO_NULL, &mpiout);
   MPI_File_read(mpiout, u, 5, MPI_FLOAT, &status);
   MPI_File_close(&mpiout); */
   fopen("l2464f211b600m0102m0509m635a.0001000", "rb");
   fread(gactionout, sizeof(int), 4, out);
   fread(betaout, sizeof(float), 1, out);
   fread(alatout, sizeof(float), 4, out);
   fread(tadout, sizeof(float), 4, out);
   fread(icfgsaveout, sizeof(int), 1, out);
   fread(myidout, sizeof(int), 1, out);
   fread(u, sizeof(float), 442368, out);
   fclose(out);

  // fclose(out);
  for (i = 0; i < 442368; i++) {
    printf("myid is %i and my floats are %f\n", myid, u[i]);
   }
   // printf("size of gauge field floats is %f\n", sizefl);
   //printf("gaction = %i %i %i %i\n", gactionout[0], gactionout[1], gactionout[2], gactionout[3]);
   //printf("beta = %f\n", betaout[0]);
  // printf("alat = %f %f %f %f\n", alatout[0], alatout[1], alatout[2], alatout[3]);
  // printf("tad = %f %f %f %f\n", tadout[0], tadout[1], tadout[2], tadout[3]);
  // printf("icfgsave = %i\n", icfgsaveout[0]);
   //printf("myidin = %i\n", myidout[0]);

   } 
   //printf("myid is %i and my first five floats are %f %f %f %f %f %f\n", myid, u[0], u[1], u[2], u[3], u[4], u[5]);
  
  
   /*if (myid == 1) {
     fopen("l2464f211b600m0102m0509m635a.0001001", "rb");
     fread(myu, sizeof(float), 6, out);
     fclose(out);
     printf("myid is %i and my first five floats are %f %f %f %f %f %f\n", myid, myu[0], myu[1], myu[2], myu[3], myu[4], myu[5]);
   } */
 
   ierr = MPI_Finalize();

  /* if (myid == 0) {
    printf("magic_number  %d\n", *milc_magic_number);
    printf("byterev flag is %d\n", byterevflag);
    printf("nx ny nz nt   %d %d %d %d\n", dim[0], dim[1], dim[2], dim[3]);
    printf("startfile     %s\n", startfile);
    printf("checksums: %lx %lx\n", check[0], check[1]);
    printf("%f\n", gf[0][0][0][0][0]);
    printf("%f\n", buffer[0]);
    printf("%f\n", gf[1][0][0][0][0]);
    printf("%f\n", buffer[1]);
    printf("%f\n", gf[2][0][0][0][0]);
    printf("%f\n", buffer[2]);
    printf("%f\n", gf[3][0][0][0][0]);
    printf("%f\n", buffer[3]);
    printf("%f\n", gf[4][0][0][0][0]);
    printf("%f\n", buffer[4]);
    printf("%f\n", gf[5][0][0][0][0]);
    printf("%f\n", buffer[5]);
    printf("%f\n", gf[0][0][0][0][1]);
   }

   if (myid == 1) {
    printf("%f\n", gf[0][0][0][0][0]);
    printf("%f\n", buffer[0]);
    printf("%f\n", gf[1][0][0][0][0]);
    printf("%f\n", buffer[1]);
    printf("%f\n", gf[2][0][0][0][0]);
    printf("%f\n", buffer[2]);
    printf("%f\n", gf[3][0][0][0][0]);
    printf("%f\n", buffer[3]);
    printf("%f\n", gf[4][0][0][0][0]);
    printf("%f\n", buffer[4]);
    printf("%f\n", gf[5][0][0][0][0]);
    printf("%f\n", buffer[5]);
  }  //end comment

   if (myid == 2) {
    printf("%f\n", gf[0][0][0][0][0]);
    printf("%f\n", buffer[0]);
    printf("%f\n", gf[1][0][0][0][0]);
    printf("%f\n", buffer[1]);
    printf("%f\n", gf[2][0][0][0][0]);
    printf("%f\n", buffer[2]);
    printf("%f\n", gf[3][0][0][0][0]);
    printf("%f\n", buffer[3]);
    printf("%f\n", gf[4][0][0][0][0]);
    printf("%f\n", buffer[4]);
    printf("%f\n", gf[5][0][0][0][0]);
    printf("%f\n", buffer[5]);
  } */

  //free(gf);

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
} // byterevn

















