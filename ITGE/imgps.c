// This program reads the gridded visibilities,and make image write in FITS format 
# include <stdio.h>
# include <math.h>
# include <stdlib.h>
# include <fitsio.h>
# include <unistd.h>

double CRVAL[4],CDELT[4],CRPIX[4];
long naxes[4],naxesim[4];

void BREAD_HDR(char *);
void BWRITE_HDR(char *outfile);


int main(int argc, char *argv[])
{ 
  char INFITS[128],OUTFITS[128];
  
  if(argc!=3){
    printf("Usage: %s <input FITS file> <output FITS file>\n", argv[0]);
    return 1;
  }
  
  sscanf(argv[1],"%s",INFITS);
  sscanf(argv[2],"%s",OUTFITS);
  
  fitsfile *fptr,*fptr1;
  int Nu,Nv,Nuv,n_ave,chan;
  long Nvis;
  int status=0,anynull;
  long fpixel[4],lpixel[4],inc[4]={1,1,1,1};
  double *nulval;  
  double *GV,*numerator,nullval;
  double Ag2,Bg;
  int ii,jj,index;
  
  if(access(INFITS, F_OK)!=0){
    printf("Input File %s does not exists\n",INFITS);
    exit(0);
  }
  
  BREAD_HDR(INFITS);
  
  printf("dU=%lf naxes vis={%ld,%ld,%ld,%ld}\n",CDELT[0],naxes[0],naxes[1],naxes[2],naxes[3]); 
  Nu=(int) naxes[0];Nv=(int) naxes[1]; n_ave=(int) naxes[3];
  Nvis=(long) Nu*Nv*6;Nuv= Nu*Nv;
  printf("Nu=%d Nv=%d Nvis=%ld Nuv=%d n_ave=%d\n",Nu,Nv,Nvis,Nuv,n_ave);

  if(access(OUTFITS, F_OK)==0){
    printf("Output File %s exists\n",OUTFITS);
    exit(0);
  }
  
  printf("Creating new file and writing header\n\n");
  
  naxes[2]=1;
  //other values will be same as input fits file    
  BWRITE_HDR(OUTFITS);
  
  GV = (double*) calloc(Nvis, sizeof(double));
  numerator = (double*) calloc(Nuv, sizeof(double));
  
  
  fits_open_file(&fptr,INFITS,READONLY,&status);
  fits_open_file(&fptr1,OUTFITS,READWRITE,&status);
  
  for(chan=0;chan<n_ave;++chan)
    { 
      printf("chan=%d\n",chan); 
      fpixel[0]=1;fpixel[1]=1;fpixel[2]=1;fpixel[3]=chan+1;
      lpixel[0]=(long)Nu;lpixel[1]=(long)Nv;lpixel[2]=6;lpixel[3]=chan+1;
      fits_read_subset(fptr,TDOUBLE,fpixel,lpixel,inc,nulval,GV,&anynull,&status);
      
      for(jj=0;jj<Nv;jj++)
      	for(ii=0;ii<Nu;ii++)
      	  {
	    index=jj*Nu+ii;
	    Ag2=GV[0*Nu*Nv+index]*GV[0*Nu*Nv+index]+GV[1*Nu*Nv+index]*GV[1*Nu*Nv+index];
	    Bg=GV[4*Nu*Nv+index];
		
	    numerator[index]=(Ag2-Bg);
	  }
      fpixel[0]=1;fpixel[1]=1;fpixel[2]=1;fpixel[3]=chan+1;
      lpixel[0]=(long)Nu;lpixel[1]=(long)Nv;lpixel[2]=1;lpixel[3]=chan+1;
      fits_write_subset(fptr1,TDOUBLE,fpixel,lpixel,numerator,&status);
    }
  
  fits_close_file(fptr,&status);
  fits_close_file(fptr1,&status);

}
