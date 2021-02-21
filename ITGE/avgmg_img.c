//This program reads all Mg.fits files and average them
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
  int ii,jj,pp,chan,index,Nu,Nv,n_ave,Nrea;
  
  if(argc!=4)
    {
      printf("Usage: %s <input file initial> <ouput FITS file> <Nrea>\n", argv[0]);
      return 1;
    }
  
  sprintf(INFITS, "%s%d.fits",argv[1],1);
  sscanf(argv[2],"%s",OUTFITS);
  Nrea=atoi(argv[3]);  

  fitsfile *fptr;
  int status=0,anynull;
  long fpixel[4],lpixel[4];
  double nullval;  


  if(access(INFITS, F_OK)!=0){
    printf("Input File %s does not exists\n",INFITS);
    exit(0);
  }
  if(access(OUTFITS, F_OK)==0){
    printf("Output File %s exists\n",OUTFITS);
    exit(0);
  }

  BREAD_HDR(INFITS);

  Nu=(int)naxes[0];Nv=(int)naxes[1];n_ave=(int)naxes[3];
  printf(" Nu=%d Nv=%d n_ave=%d\n",Nu,Nv,n_ave);

  //write output header
  printf("Creating new file and writing header\n\n");

  //all values will be same as input fits file    
  BWRITE_HDR(OUTFITS);

  /*------------------------------------------------------------------------*/
  
  double *out,*GV;
  long Ndata;
 
  Ndata=Nu*Nv*naxes[2]*n_ave;
  printf("Ndata=%ld\n",Ndata);
  GV= (double*)calloc(Ndata,sizeof(double*));
  out=(double*)calloc(Ndata,sizeof(double*));

  for(pp=0;pp<Nrea;++pp)
    { 
      sprintf(INFITS, "%s%d.fits",argv[1],pp+1);
      fits_open_file(&fptr,INFITS,READONLY,&status);
      
      fits_read_img(fptr,TDOUBLE,1,Ndata,&nullval,GV,&anynull,&status);
 
      for(ii=0;ii<Ndata;++ii)
	{
	  out[ii]+=GV[ii];
	}
      fits_close_file(fptr,&status);
    }
  
  for(ii=0;ii<Ndata;++ii)
    {
      out[ii]/=(1.*Nrea);
    }
  
 
  fits_open_file(&fptr,OUTFITS,READWRITE,&status);
  fits_write_img(fptr, TDOUBLE, 1, Ndata, out, &status);
  fits_close_file(fptr, &status);
}
