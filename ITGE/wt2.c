//Fill wt2
# include <stdio.h>
# include <math.h>
# include <stdlib.h>
# include <fitsio.h>
# include <unistd.h>

int main(int argc, char *argv[])
{ 
  char INFITS[128];
  int N,NL=0,AN[3]={0,0,0},ii,jj,pp,index,chan;
   
  if(argc!=3)
    {
      printf("Usage: %s <input FITS file> <N>\n", argv[0]);
      return 1;
    }
  
  sscanf(argv[1],"%s",INFITS);
  sscanf(argv[2],"%d",&N);
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

  fitsfile *fptr;
  int naxis=4,Nu,Nv;
  long naxes[naxis],fpixel[naxis],lpixel[naxis],inc[naxis];
  int status=0,anynull;
  double *datareim,*img,*nulval;
 
  if(access(INFITS, F_OK)!=0){
    printf("Input File %s does not exists\n",INFITS);
    exit(0);
  }
 

  fits_open_file(&fptr,INFITS,READWRITE,&status);
  fits_read_keys_lng(fptr,"NAXIS",1,naxis,naxes,&naxis,&status);
  printf("naxes={%ld,%ld,%ld,%ld}\n",naxes[0],naxes[1],naxes[2],naxes[3]);

  Nu=(int)naxes[0];Nv=(int) naxes[1];

  datareim= (double*)calloc(Nu*Nv*2,sizeof(double*));
  img= (double*)calloc(Nu*Nv*2,sizeof(double*));
  
  fpixel[0]=1;fpixel[1]=1;fpixel[2]=1;fpixel[3]=1;
  lpixel[0]=Nu;lpixel[1]=Nv;lpixel[2]=2;lpixel[3]=1;
  inc[0]=1;inc[1]=1;inc[2]=1;inc[3]=1;

  printf("\nvis read  [%ld,%ld,%ld,%ld] [%ld,%ld,%ld,%ld]\n\n",fpixel[0],fpixel[1],fpixel[2],fpixel[3],lpixel[0],lpixel[1],lpixel[2],lpixel[3]);

  fits_read_subset(fptr,TDOUBLE,fpixel,lpixel,inc,nulval,datareim,&anynull,&status);

  for(jj=0;jj<Nv;++jj)
    for(ii=0;ii<Nu;++ii)
      {
	index=jj*Nu+ii;
	img[0*Nu*Nv+index]=pow(datareim[0*Nu*Nv+index],2.)+pow(datareim[1*Nu*Nv+index],2.);
	img[1*Nu*Nv+index]=0.;
      }

  for(chan=0;chan<naxes[3];++chan)
    {
      for(pp=0;pp<NL;++pp)
	{
	  fpixel[0]=1;fpixel[1]=1;fpixel[2]=(long) AN[pp];fpixel[3]=chan+1;
	  lpixel[0]=(long) Nu;lpixel[1]=(long)Nv;lpixel[2]=(long)(AN[pp]+1);lpixel[3]=chan+1;
	  printf("\nvis^2 write  [%ld,%ld,%ld,%ld] [%ld,%ld,%ld,%ld]\n",fpixel[0],fpixel[1],fpixel[2],fpixel[3],lpixel[0],lpixel[1],lpixel[2],lpixel[3]);
	  
	  fits_write_subset(fptr,TDOUBLE,fpixel,lpixel,img,&status);
	}
    }
  fits_close_file(fptr,&status);
  free(img);
  free(datareim);
}
