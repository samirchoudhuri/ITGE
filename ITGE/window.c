//Calculate the window function in the image plane

# include <stdio.h>
# include <math.h>
# include <stdlib.h>
# include <fitsio.h>
# include <unistd.h>

void SREAD_HDR(char *);
void win(double * ,long ,long , char *);
void disk(double * ,long ,long , char *);
void inver_disk(double * ,long ,long , char *);
void exp_sinc(double * ,long ,long , char *);
void radial(double * ,long ,long , char *);
void butterworth(double * ,long ,long , char *);
void disk_gauss(double * ,long ,long , char *);
void butterworthannular(double * ,long ,long , char *);

double theta0,theta02;
double CRVAL[4],CDELT[4],CRPIX[4];
long naxes[4],naxesim[4];

int main(int argc, char *argv[])
{ 
  char INFITS[128],OUTFITS[128];
  int N,NL=0,AN[3]={0,0,0},ii,jj,pp,index,chan;
  FILE *fp;

  if(argc<5)
    {
      printf("Usage: %s <input FITS file> <ouput FITS file> <N> <theta0 arcmin> <theta02 arcmin> <pt catalogue>\n", argv[0]);
      return 1;
    }
  
  sscanf(argv[1],"%s",INFITS);
  sscanf(argv[2],"%s",OUTFITS);
  N=atoi(argv[3]);
  sscanf(argv[4],"%lf",&theta0);
  sscanf(argv[5],"%lf",&theta02);
  printf("N=%d theta0=%lf theta02=%lf\n",N,theta0,theta02);
 
  if(N%2==1)
    {
      AN[NL]=1;
      NL++;
    }
  if((N/2)%2==1)
    {
      AN[NL]=3;
      NL++;
    }
  if((N/4)%2==1)
    {
      AN[NL]=5;
      NL++;
    }
  for(ii=0;ii<3;++ii)
    printf("NL=%d AN[%d]=%d\n",NL,ii,AN[ii]);

  fitsfile *fptr,*fptrout;
  int naxis=4,nfound;
  char OPT[1];
  long dimleng[4],fpixel[naxis],lpixel[naxis],inc[naxis];
  int status=0,anynull;
  double *img,*nulval;
 
  /*-----OPEN INPUT FITS FILE and READ HEADER-----*/

  if(access(INFITS, F_OK)!=0){
    printf("Input File %s does not exists\n",INFITS);
    exit(0);
  }

  SREAD_HDR(INFITS);

  printf("pixel=%lf naxes image={%ld,%ld,%ld,%ld}\n",CDELT[0],naxesim[0],naxesim[1],naxesim[2],naxesim[3]);

  if(access(OUTFITS, F_OK)!=0){
    printf("Creating new file and writing header\n\n");
    
    fits_open_file(&fptr,INFITS,READONLY,&status);
    fits_create_file(&fptrout,OUTFITS,&status);
    fits_copy_header(fptr,fptrout,&status);
    fits_close_file(fptr,&status);
  }


  else{
    printf("\nFile %s already exists\n",OUTFITS);
    
    fits_open_file(&fptrout,OUTFITS,READWRITE,&status);
    fits_read_keys_lng(fptrout,"NAXIS",1,naxis,dimleng,&nfound,&status);

    if((dimleng[0]!=naxesim[0])||(dimleng[1]!=naxesim[1])||(dimleng[2]!=naxesim[2])||(dimleng[3]!=naxesim[3]))
      {
	printf("Dimension does not matches between existing image file\n");
	exit(0);
      }
  }

  img= (double*)calloc(naxesim[0]*naxesim[1],sizeof(double*));

  //win(img,naxesim[0],naxesim[1],argv[6]);
  //disk(img,naxesim[0],naxesim[1],argv[6]);
  //inver_disk(img,naxesim[0],naxesim[1],argv[6]);
  //exp_sinc(img,naxesim[0],naxesim[1],argv[6]);
  //radial(img,naxesim[0],naxesim[1],argv[6]);
  butterworth(img,naxesim[0],naxesim[1],argv[6]);
  //butterworthannular(img,naxesim[0],naxesim[1],argv[6]);
  //disk_gauss(img,naxesim[0],naxesim[1],argv[6]);
 
  for(chan=0;chan<naxesim[3];++chan)
    { 
      for(pp=0;pp<NL;++pp)
	{
	  fpixel[0]=1;fpixel[1]=1;fpixel[2]=(long)(1+(AN[pp]-1)/2);fpixel[3]=chan+1;
	  lpixel[0]=naxesim[0];lpixel[1]=naxesim[1];lpixel[2]=(long)(1+(AN[pp]-1)/2);lpixel[3]=chan+1;
	  inc[0]=1;inc[1]=1;inc[2]=1;inc[3]=1;
      
	  printf("\nimage write  [%ld,%ld,%ld,%ld] [%ld,%ld,%ld,%ld]\n\n",fpixel[0],fpixel[1],fpixel[2],fpixel[3],lpixel[0],lpixel[1],lpixel[2],lpixel[3]);
      
	  fits_write_subset(fptrout,TDOUBLE,fpixel,lpixel,img,&status);
	}
    }
  fits_close_file(fptrout,&status);
  free(img);
}
