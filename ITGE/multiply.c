// This program reads data from two images and multiply them, write in a new FITS file
//Image dim should be 4
# include <stdio.h>
# include <math.h>
# include <stdlib.h>
# include <fitsio.h>
# include <unistd.h>

int main(int argc, char *argv[])
{ 
  char INFITS[128],INFITS2[128],OUTFITS[128];
  int N,NL=0,AN[3]={0,0,0},ii,jj,pp,index,chan;
   
  if(argc!=5)
    {
      printf("Usage: %s <input FITS file> <input2 FITS FILE> <ouput FITS file> <N>\n", argv[0]);
      return 1;
    }
  
  sscanf(argv[1],"%s",INFITS);
  sscanf(argv[2],"%s",INFITS2);
  sscanf(argv[3],"%s",OUTFITS);
  sscanf(argv[4],"%d",&N);
  printf("N=%d\n",N);

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
  
  fitsfile *fptr,*fptr2,*fptrout;
  int naxis=4;
  char comment[FLEN_COMMENT],OPT[1];
  long naxes[naxis],naxes2[naxis],dimleng[naxis],fpixel[naxis],lpixel[naxis],inc[naxis];
  int status=0,anynull;
  double *img,*img2,*nulval;
 
  /*-----OPEN INPUT FITS FILE and READ HEADER-----*/

  if(access(INFITS, F_OK)!=0){
    printf("Input File %s does not exists\n",INFITS);
    exit(0);
  }

  if(access(INFITS2, F_OK)!=0){
    printf("Input File2 %s does not exists\n",INFITS2);
    exit(0);
  }

  fits_open_file(&fptr,INFITS,READONLY,&status);
  fits_read_keys_lng(fptr,"NAXIS",1,naxis,naxes,&naxis,&status);
  printf("naxes={%ld,%ld,%ld,%ld}\n",naxes[0],naxes[1],naxes[2],naxes[3]); 
  
  fits_open_file(&fptr2,INFITS2,READONLY,&status);
  fits_read_keys_lng(fptr2,"NAXIS",1,naxis,naxes2,&naxis,&status);
  printf("naxes2={%ld,%ld,%ld,%ld}\n",naxes2[0],naxes2[1],naxes2[2],naxes2[3]); 
  
  if((naxes[0]!=naxes2[0])||(naxes[1]!=naxes2[1])||(naxes[2]!=naxes2[2])||(naxes[3]!=naxes2[3]))
      {
	printf("Dimension does not matches between two input images\n");
	exit(0);
      }
 

  if(access(OUTFITS, F_OK)!=0){
    printf("Creating new file and writing header\n\n");
    
    fits_create_file(&fptrout,OUTFITS,&status);
    fits_copy_header(fptr,fptrout,&status);
  }

  /*----------------------IF EXIST CHECK DIMENSION WITH UVFITS------------------*/

  else{
    printf("\nFile %s already exists\n",OUTFITS);
    
    fits_open_file(&fptrout,OUTFITS,READWRITE,&status);
    fits_read_keys_lng(fptrout,"NAXIS",1,naxis,dimleng,&naxis,&status);
    if((dimleng[0]!=naxes[0])||(dimleng[1]!=naxes[1])||(dimleng[2]!=naxes[2])||(dimleng[3]!=naxes[3]))
      {
	printf("Dimension does not matche between existing image file\n");
	exit(0);
      }
  }
  /*------------------------------------------------------------------------*/

  img= (double*)calloc(naxes[0]*naxes[1],sizeof(double*));
  img2= (double*)calloc(naxes[0]*naxes[1],sizeof(double*));

  for(chan=0;chan<naxes[3];++chan)
    for(pp=0;pp<NL;++pp)
      {
	fpixel[0]=1;fpixel[1]=1;fpixel[2]=(long)(1+(AN[pp]-1)/2);fpixel[3]=chan+1;
	lpixel[0]=naxes[0];lpixel[1]=naxes[1];lpixel[2]=(long)(1+(AN[pp]-1)/2);lpixel[3]=chan+1;
	inc[0]=1;inc[1]=1;inc[2]=1;inc[3]=1;
	
	printf("\nimage read  [%ld,%ld,%ld,%ld] [%ld,%ld,%ld,%ld]\n\n",fpixel[0],fpixel[1],fpixel[2],fpixel[3],lpixel[0],lpixel[1],lpixel[2],lpixel[3]);
	
	fits_read_subset(fptr,TDOUBLE,fpixel,lpixel,inc,nulval,img,&anynull,&status);
	fits_read_subset(fptr2,TDOUBLE,fpixel,lpixel,inc,nulval,img2,&anynull,&status);

	for(jj=0;jj<naxes[1];++jj)
	  for(ii=0;ii<naxes[0];++ii)
	    {
	      index=jj*naxes[0]+ii;
	      img[index]*=img2[index];
	    }
	fits_write_subset(fptrout,TDOUBLE,fpixel,lpixel,img,&status);
      }
  
  fits_close_file(fptr,&status);
  fits_close_file(fptr2,&status);
  fits_close_file(fptrout,&status);
  free(img);
  free(img2);
}
