// This program reads the gridded visibilities,and collapse the different channel and side bands

# include <stdio.h>
# include <math.h>
# include <stdlib.h>
# include <fitsio.h>
# include <nr.h>
# include <nrutil.h>
# include <fftw3.h>
# include <unistd.h>

double CRVAL[4],CDELT[4],CRPIX[4];
long naxes[4],naxesim[4];

void BREAD_HDR(char *);
void SWRITE_HDR(char *);
void SREAD_HDR(char *);

int main(int argc, char *argv[])
{ 
  char INFITS[128],OUTFITS[128];
  int ii,pp;
  
  if(argc!=3)
    {
      printf("Usage: %s <input FITS file1> <ouput FITS file>\n", argv[0]);
      return 1;
    }
  
  sscanf(argv[1],"%s",INFITS);
  sscanf(argv[2],"%s",OUTFITS);
  
  fitsfile *fptr;
  char OPT[1];
  long fpixel[4],lpixel[4],inc[4];
  int status=0,anynull;
  double *nulval;

  /*-----OPEN INPUT FITS FILE and READ HEADER-----*/

  if(access(INFITS, F_OK)!=0){
    printf("Input File %s does not exists\n",INFITS);
    exit(0);
  }

  BREAD_HDR(INFITS);
  printf("dU=%lf naxes1 vis1={%ld,%ld,%ld,%ld}\n",CDELT[0],naxes[0],naxes[1],naxes[2],naxes[3]);
 
  
  if(access(OUTFITS, F_OK)==0){
    printf("File %s exists\n",OUTFITS);
    return 0;
  }
  
  printf("Creating new file and writing header\n\n");
  
  naxesim[0]=naxes[0];naxesim[1]=naxes[1];naxesim[2]=naxes[2];naxesim[3]=1;
  //other values will be same as input fits file    
  SWRITE_HDR(OUTFITS);
       
  /*------------------------------------------------------------------------*/
  
  double *out,*datareim;
  long Ndata;
  Ndata=naxes[0]*naxes[1]*naxes[2];
  datareim= (double*)calloc(Ndata,sizeof(double*));
  out=(double*)calloc(Ndata,sizeof(double*));
  
  fits_open_file(&fptr,INFITS,READONLY,&status);//output fits file pointer already opened
  for(pp=0;pp<naxes[3];++pp)
    {
      fpixel[0]=1;fpixel[1]=1;fpixel[2]=1;fpixel[3]=pp+1;
      lpixel[0]=naxes[0];lpixel[1]=naxes[1];lpixel[2]=naxes[2];lpixel[3]=pp+1;
      inc[0]=1;inc[1]=1;inc[2]=1;inc[3]=1;
      
      printf("\nvis read  [%ld,%ld,%ld,%ld] [%ld,%ld,%ld,%ld]\n",fpixel[0],fpixel[1],fpixel[2],fpixel[3],lpixel[0],lpixel[1],lpixel[2],lpixel[3]);
      
      fits_read_subset(fptr,TDOUBLE,fpixel,lpixel,inc,nulval,datareim,&anynull,&status);

      for(ii=0;ii<Ndata;++ii)
	{
	  out[ii]+=datareim[ii];
	}

    }
  fits_close_file(fptr,&status);

  for(ii=0;ii<Ndata;++ii)
    {
      out[ii]=(out[ii]/(1.*naxes[3]));
    }
  
  fits_open_file(&fptr,OUTFITS,READWRITE,&status);
  fits_write_img(fptr, TDOUBLE, 1, Ndata, out, &status);
  fits_close_file(fptr, &status);

}
