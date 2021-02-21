// This program reads the gridded visibilities,and make image write in FITS format 
# include <stdio.h>
# include <stdlib.h>
# include <math.h>
# include <string.h>
# include <fitsio.h>
# include <unistd.h>

double CRVAL[4],CDELT[4],CRPIX[4];

long naxes[4],naxesim[4];

void BREAD_HDR(char *);

int main(int argc, char *argv[])
{ 
  char INFITS[128];
  double Umax,Umin,kv,V1,theta0;
  double UvalGrid,*binval,*uval,*weight,Ag2,Bg,Cg2;
  int NUGrid,ii,jj,ii1,index;
  char param[128];
  int Nbin;
  FILE *fp;


  if(argc!=4){
      printf("Usage: %s <input FITS file> <input file> <outputfile>\n", argv[0]);
      return 1;
    }
  
  sscanf(argv[1],"%s",INFITS);
   
  fp=fopen(argv[2],"r");
  fscanf(fp,"%lf%lf%d",&Umax,&Umin,&Nbin);
  fclose(fp);
  
  printf("Umax=%.2f Umin=%.2f theta0=%.2fmin Nbin=%d\n",Umax,Umin,theta0,Nbin);
  kv=(1.*(Nbin+0))/log10(Umax/Umin);

  //theta0=(M_PI*theta0)/(180.*60.);//theta0 in rad
  //V1=M_PI*pow(theta0,2.)/2.;
  V1=1.;
  printf("V1=%e\n",V1);
  
  
  fitsfile *fptr;
  int Nu,Nv;
  long Nvis, *count, kk;
  int status=0,anynull;
  double *GV,nullval;
  
  if(access(INFITS, F_OK)!=0){
    printf("Input File %s does not exists\n",INFITS);
    return 1;
  }
  
  BREAD_HDR(INFITS);
  
  printf("dU=%lf naxes vis={%ld,%ld,%ld,%ld}\n",CDELT[0],naxes[0],naxes[1],naxes[2],naxes[3]); 
  Nu=(int) naxes[0];Nv=(int) naxes[1]; Nvis=(long) Nu*Nv*6;
  printf("Nu=%d Nv=%d Nvis=%ld\n",Nu,Nv,Nvis);
 
  GV = (double*) calloc(Nvis, sizeof(double));

  fits_open_file(&fptr,INFITS,READONLY,&status);
  fits_read_img(fptr,TDOUBLE,1,Nvis,&nullval,GV,&anynull,&status);
  fits_close_file(fptr,&status);
 
  binval=(double*)calloc(Nbin,sizeof(double));
  weight=(double*)calloc(Nbin,sizeof(double));
  uval=(double*)calloc(Nbin,sizeof(double));
  count = (long *)calloc(Nbin, sizeof(long)); 

  //Start binning
  for(jj=0;jj<Nv;jj++)
    for(ii=0;ii<Nu;ii++)
      {
	index=jj*Nu+ii;
	ii1=(ii>Nu/2) ? Nu-ii : ii; 
	UvalGrid=CDELT[0]*sqrt(ii1*ii1+jj*jj);
	if( (UvalGrid> Umin) && (UvalGrid< Umax) )
	  {
	    NUGrid=(int)floor(kv*log10(UvalGrid/Umin));
	    NUGrid = (NUGrid>0) ? NUGrid:0;
	    Ag2=GV[0*Nu*Nv+index]*GV[0*Nu*Nv+index]+GV[1*Nu*Nv+index]*GV[1*Nu*Nv+index];
	    Bg=GV[4*Nu*Nv+index];
	    
	    Cg2=GV[2*Nu*Nv+index]*GV[2*Nu*Nv+index]+GV[3*Nu*Nv+index]*GV[3*Nu*Nv+index];
	    binval[NUGrid] += (Ag2-Bg);
	    uval[NUGrid] += (UvalGrid*Cg2); 
	    weight[NUGrid] += (1.*Cg2);
	    count[NUGrid] ++;
	  }
      }
  fp = fopen(argv[3],"w");
  for(ii=0;ii<Nbin;++ii)
    {
      if(weight[ii]>0.)
	//fprintf(fp,"%e %e %e %ld\n",uval[ii]/weight[ii],binval[ii]/(weight[ii]*V1),weight[ii],count[ii]);
	fprintf(fp,"%e %e %e\n",uval[ii]/weight[ii],binval[ii]/(weight[ii]*V1),weight[ii]);
      else
	fprintf(fp,"%e %e %e\n",0.,-1.e10,0.);
    }
  fclose(fp);
  return EXIT_SUCCESS;
}
