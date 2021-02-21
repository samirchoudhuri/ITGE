// This program reads the gridded visibilities,and make image write in FITS format 
# include <stdio.h>
# include <math.h>
# include <stdlib.h>
# include <fitsio.h>
# include <nr.h>
# include <nrutil.h>
# include <fftw3.h>
# include <unistd.h>

void SREAD_HDR(char *);
void BWRITE_HDR(char *);
void BREAD_HDR(char *);

double CRVAL[4],CDELT[4],CRPIX[4];
long naxes[4],naxesim[4];

int main(int argc, char *argv[])
{ 
  char INFITS[128],OUTFITS[128];
  int N,NL=0,AN[3]={0,0,0},ii,jj,pp,index,index1,chan;
   
  if(argc!=4){
    printf("Usage: %s <input FITS file> <ouput FITS file> N\n", argv[0]);
    return 1;
  }
  
  sscanf(argv[1],"%s",INFITS);
  sscanf(argv[2],"%s",OUTFITS);
  sscanf(argv[3],"%d",&N);
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
  
  fitsfile *fptr,*fptrim;
  int Nu,Nv;
  char OPT[1];
  long fpixel[4],lpixel[4],inc[4],fimpixel[4],limpixel[4];
  int status=0,anynull;
  double *nulval;
  double dU;

  /*-----OPEN INPUT FITS FILE and READ HEADER-----*/

  if(access(INFITS, F_OK)!=0){
    printf("Input File %s does not exists\n",INFITS);
    exit(0);
  }
  
  SREAD_HDR(INFITS);

  printf("pixel=%lf naxes image={%ld,%ld,%ld,%ld}\n",CDELT[0],naxesim[0],naxesim[1],naxesim[2],naxesim[3]);
  
  Nu=(int) naxesim[0];Nv=(int) (naxesim[1]/2)+1;
  printf("Nu=%d Nv=%d\n",Nu,Nv);

  /*------------------CHECK OUTPUT FILE-------------------------------------*/
  /*-------------IF NOT EXIST-----CREATE OUTPUT FITS FILE-------------------*/
  if(access(OUTFITS, F_OK)!=0){
    printf("Creating new file and writing header\n\n");
    
    dU=(double) ((180.*60.)/(M_PI*Nu*CDELT[0]));
    naxes[0]=(long) Nu;naxes[1]=(long) Nv;naxes[2]=(long) naxesim[2]*2;naxes[3]=naxesim[3];
    CRVAL[0]=0.;CRVAL[1]=0.;CRVAL[2]=1.;CRVAL[3]=CRVAL[3];
    CDELT[0]=dU;CDELT[1]=dU;CDELT[2]=1.;CDELT[3]=CDELT[3];
    CRPIX[0]=1.;CRPIX[1]=1.;CRPIX[2]=1.;CRPIX[3]=CRPIX[3];
    printf("naxes vis={%ld,%ld,%ld,%ld}\n",naxes[0],naxes[1],naxes[2],naxes[3]); 
  
    BWRITE_HDR(OUTFITS);
  } 

  /*----------------------IF EXIST CHECK DIMENSION WITH UVFITS------------------*/

  else{
    printf("\nFile %s already exists\n",OUTFITS);
    
    BREAD_HDR(OUTFITS);
    if((naxes[0]!=naxesim[0])||(naxes[1]!=(naxesim[1]/2)+1)||(naxes[2]!=naxesim[2]*2)||(naxes[3]!=naxesim[3]))
      {
	printf("Dimension does not matches with IMAGE file\n");
	exit(0);
      }
  }
  /*------------------------------------------------------------------------*/

  fits_open_file(&fptr,INFITS,READONLY,&status);
  fits_open_file(&fptrim,OUTFITS,READWRITE,&status);

  fftw_complex *in;
  double *out,*datareim,*img;
  fftw_plan p;

  datareim= (double*)calloc(Nu*Nu,sizeof(double*));
  out=(double*)calloc ((Nu*(Nu+2)), sizeof (double));
  in=(fftw_complex*)&out[0]; 
  img= (double*)calloc(Nu*Nv*2,sizeof(double*));
  
  p= fftw_plan_dft_r2c_2d ( Nu, Nu, out,in, FFTW_ESTIMATE );
  

  for(chan=0;chan<naxesim[3];++chan)
    {
      for(pp=0;pp<NL;++pp)
	{
	  fpixel[0]=1;fpixel[1]=1;fpixel[2]=(long) (1+(AN[pp]-1)/2);fpixel[3]=chan+1;
	  lpixel[0]=Nu;lpixel[1]=Nu;lpixel[2]=(long)(1+(AN[pp]-1)/2);lpixel[3]=chan+1;
	  inc[0]=1;inc[1]=1;inc[2]=1;inc[3]=1;

	  fimpixel[0]=1;fimpixel[1]=1;fimpixel[2]=AN[pp];fimpixel[3]=chan+1;
	  limpixel[0]=(long)Nu;limpixel[1]=(long)Nv;limpixel[2]=AN[pp]+1;limpixel[3]=chan+1;
	  
	  printf("\nimage read  [%ld,%ld,%ld,%ld] [%ld,%ld,%ld,%ld]\n",fpixel[0],fpixel[1],fpixel[2],fpixel[3],lpixel[0],lpixel[1],lpixel[2],lpixel[3]);
	  printf("vis write=[%ld %ld %ld %ld] [%ld %ld %ld %ld]\n",fimpixel[0],fimpixel[1],fimpixel[2],fimpixel[3],limpixel[0],limpixel[1],limpixel[2],limpixel[3]);
	  
	  fits_read_subset(fptr,TDOUBLE,fpixel,lpixel,inc,nulval,datareim,&anynull,&status);
	  
	  for(jj=0;jj<Nu;++jj)
	    for(ii=0;ii<Nu;++ii)
	      {
		index=jj*Nu+ii;
		index1=ii*(Nu+2)+jj;
		out[index1]=datareim[index];
	      }
	  
	  fftw_execute ( p);
	  
	  for(jj=0;jj<Nu;++jj)
	    for(ii=0;ii<Nv;++ii)
	      {
		index=ii*Nu+jj;
		index1=jj*Nv+ii; // transpose for making FITS from row major to column major format 
		img[0*Nu*Nv+index]= pow(-1.,ii+jj)*in[index1][0]/(Nu*Nu);
		img[1*Nu*Nv+index]= pow(-1.,ii+jj)*in[index1][1]/(Nu*Nu);
	      }
	  fits_write_subset(fptrim,TDOUBLE,fimpixel,limpixel,img,&status);
	}
    }
  
  fftw_destroy_plan ( p);
  fits_close_file(fptr,&status);
  fits_close_file(fptrim,&status);
  free(datareim);
  free(out);
  free(img);
}
