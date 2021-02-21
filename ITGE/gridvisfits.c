// This program reads a multichannel visibility data, grid them and write in a multichannel FITS image format 
# include <stdio.h>
# include <math.h>
# include <stdlib.h>
# include <fitsio.h>
# include <nr.h>
# include <nrutil.h>
# include <unistd.h>
# include "read_fits_func.h"

void BWRITE_HDR(char *);

long naxes[4],naxesim[4];
double CRVAL[4],CDELT[4],CRPIX[4];

extern float *uu,*vv,**RRre,**RRim,**LLre,**LLim,del_chan,nu_chan0,chan0;
extern int flagLL;

int main(int argc, char *argv[])
{ 
  char INFITS[128],OUTFITS[128],input[128],OPT[1];
  int chan,chan1,chan2,n_ave,Ng,Nu,Nv,ii,jj,ii1,jj1;
  int nbasln;
  long Nvis,index,index1;
  double Umax,Umin,dU,*corrfact,fac;
  double *GV;
  double uuc,vvc;
  double Re,Im;
  
  FILE *fp;
  if(argc!=4)
    {
      printf("Usage: %s <input FITS file> <ouput FITS file> <input>\n", argv[0]);
      return 1;
    }
  
  sscanf(argv[1],"%s",INFITS);
  sscanf(argv[2],"%s",OUTFITS);
  sscanf(argv[3],"%s",input);
 
  //reading input parameter
  fp=fopen(input,"r");
  fscanf(fp,"%d%d",&chan1,&chan2);
  fscanf(fp,"%lf%lf%lf",&Umax,&Umin,&dU);
  fscanf(fp,"%lf%d",&fac,&flagLL);
  fclose(fp);
  
  chan1=chan1-1; chan2=chan2-1;
  n_ave=chan2-chan1+1;
  
  printf("\n\nchan1-1=%d\tchan2-1=%d\tn_ave=%d\nUmax=%lf\tUmin=%lf\tdU=%lf\n\n",chan1,chan2,n_ave,Umax,Umin,dU);
  printf("flagLL=%d\n",flagLL);

  if(access(INFITS, F_OK)!=0){
    printf("Input File %s does not exists\n",INFITS);
    exit(0);
  }
  
  //It will fill RRre[][],RRim[][],LLre[][],LLim[][]
  nbasln=readfits(INFITS,Umax,Umin,chan1,chan2);
  printf("nbasln=%d\n",nbasln);
  
  Ng =(int) ceil(Umax/dU);
  Nu = 2*Ng +2;           
  Nv = Ng + 2;   
  Nvis =(long) Nu*Nv*6*n_ave;
  
  printf("\nNg=%d\tNu=%d\tNv=%d n_ave=%d Nvis=%ld\n",Ng,Nu,Nv,n_ave,Nvis);
  GV = (double*) calloc(Nvis, sizeof(double));
  
  corrfact = (double*) calloc(n_ave, sizeof(double));
  for(ii=0; ii<n_ave; ii++)
    {    
      corrfact[ii] = 1.+ fac*(del_chan/nu_chan0)*(chan1+ii+1.+0.5-chan0);//lamda[chan0]/lambda[ii]
    }
  
  /*Gridding the visibility data*/
  for(ii=0;ii<nbasln;ii++)
    {
      for(chan=0;chan<n_ave;++chan)
	{
	  uuc=corrfact[chan]*uu[ii];
	  vvc=corrfact[chan]*vv[ii];

	  if((uuc*uuc+vvc*vvc)<= Umax*Umax)
	    {	  
	      ii1 = (int)roundf(uuc/dU);
	      ii1=(ii1<0) ? Nu+ii1 : ii1 ; 
	      jj1 = (int)roundf(vvc/dU);
	     
	      index1=(long) chan*6*Nu*Nv+jj1*Nu+ii1;
	      if(RRre[ii][chan]>-1.e7)
		{
		  Re=RRre[ii][chan];
		  Im=-1.*RRim[ii][chan];//as vis calculated with +1
		 
		  GV[0*Nu*Nv+index1] += Re;
		  GV[1*Nu*Nv+index1] += Im; // complex conjugate for -1 in FFT
		  GV[2*Nu*Nv+index1] += 1.;	
		  GV[3*Nu*Nv+index1]=0.;
		  GV[4*Nu*Nv+index1] += (Re*Re+Im*Im);
		  GV[5*Nu*Nv+index1]=0.;	  
		}
	      
	      if(LLre[ii][chan]>-1.e7)
		{
		  Re=LLre[ii][chan];
		  Im=-1.*LLim[ii][chan];//as vis calculated with +1
		 
		  GV[0*Nu*Nv+index1] += Re;
		  GV[1*Nu*Nv+index1] += Im; // complex conjugate for -1 in FFT
		  GV[2*Nu*Nv+index1] += 1.;	
		  GV[3*Nu*Nv+index1]=0.;
		  GV[4*Nu*Nv+index1] += (Re*Re+Im*Im);
		  GV[5*Nu*Nv+index1]=0.;	  
		}
	    }
	}
    }
  

  for(chan=0;chan<n_ave;++chan)
    {
      for(jj1=0;jj1<Nv;jj1=jj1+Nv-1)//make complex conjugate for jj=0 and jj=Nv-1 axis
	for(ii1=1;ii1<=Ng; ++ii1)
	  {    
	    index = (long)chan*6*Nu*Nv+jj1*Nu +ii1;
	    index1 =(long)chan*6*Nu*Nv+jj1*Nu+(Nu-ii1);
	    GV[0*Nu*Nv+index] += GV[0*Nu*Nv+index1];
	    GV[1*Nu*Nv+index] -= GV[1*Nu*Nv+index1];
	    GV[2*Nu*Nv+index] += GV[2*Nu*Nv+index1];
	    GV[4*Nu*Nv+index] += GV[4*Nu*Nv+index1];
	    
	    GV[0*Nu*Nv+index] *=0.5;
	    GV[1*Nu*Nv+index] *=0.5;
	    GV[2*Nu*Nv+index] *=0.5;
	    GV[4*Nu*Nv+index] *=0.5;
	    
	    GV[0*Nu*Nv+index1] = GV[0*Nu*Nv+index];
	    GV[1*Nu*Nv+index1] = -1.*GV[1*Nu*Nv+index];
	    GV[2*Nu*Nv+index1] = GV[2*Nu*Nv+index];
	    GV[4*Nu*Nv+index1] = GV[4*Nu*Nv+index];
	  }
      
      
      for(ii1=0;ii1<Nv;ii1=ii1+Nv-1)//ii=0 and ii=Nv-1 (Ng+1)
	for(jj1=0;jj1<Nv;jj1=jj1+Nv-1)//jj=0 and jj=Nv-1 make 4 points real
	  {
	    index=(long)chan*6*Nu*Nv+jj1*Nu+ii1;
	    GV[1*Nu*Nv+index]=0.;
	  }
    }
  
  //Write FITS format

  fitsfile *fptr;
  int status=0;
  naxes[0]=(long)Nu;naxes[1]=(long)Nv;naxes[2]=6;naxes[3]=(long)n_ave;
  CRVAL[0]=0.;CRVAL[1]=0.;CRVAL[2]=1.;CRVAL[3]=nu_chan0;
  CDELT[0]=dU;CDELT[1]=dU;CDELT[2]=1.;CDELT[3]=del_chan;
  CRPIX[0]=1.;CRPIX[1]=1.;CRPIX[2]=1.;CRPIX[3]=chan0;
  long Ndata;

  if(access(OUTFITS, F_OK)==0){
    printf("Output File exist\n",OUTFITS);
    return 0;
  }
  
  BWRITE_HDR(OUTFITS);
  
  Ndata =naxes[0]*naxes[1]*naxes[2]*naxes[3];
  
  fits_open_file(&fptr,OUTFITS,READWRITE,&status);
  fits_write_img(fptr, TDOUBLE, 1, Ndata, GV, &status);
  fits_close_file(fptr, &status);

  free(GV);
}
