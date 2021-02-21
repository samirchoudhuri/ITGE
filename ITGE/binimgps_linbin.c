// This program reads the gridded visibilities,and make image write in FITS format 
# include <stdio.h>
# include <math.h>
# include <stdlib.h>
# include <fitsio.h>
# include <unistd.h>

double CRVAL[4],CDELT[4],CRPIX[4];
long naxes[4],naxesim[4];

void BREAD_HDR(char *);

double P(double U)
{
  double pu,A=3.e5,betav=1.88;
  pu=A*pow((1000./(2.*M_PI*U)),betav);
  return pu;
}


int main(int argc, char *argv[])
{ 
  char INFITS[128],CONST[128];
  double Umax,Umin,kv,cutoff;
  int Nbin;
  FILE *fp;

  if(argc!=5){
    printf("Usage: %s <input FITS file> <input file> <outputfile> <const FITS file>\n", argv[0]);
      return 1;
  }
  
  sscanf(argv[1],"%s",INFITS);
  sscanf(argv[4],"%s",CONST);
  
  //reading input parameter
  fp=fopen(argv[2],"r");
  fscanf(fp,"%lf%lf%d%lf",&Umax,&Umin,&Nbin,&cutoff);
  fclose(fp);
  
  printf("Umax=%.2f Umin=%.2f Nbin=%d cutoff=%lf\n",Umax,Umin,Nbin,cutoff);
  kv=Nbin/(Umax-Umin);

 
  fitsfile *fptr,*fptr1;
  int Nu,Nv,Nuv,n_ave,chan;
  int status=0,anynull;
  long fpixel[4],lpixel[4],inc[4]={1,1,1,1};
  double *nulval;  
  double *GV,*GVconst,nullval;
  double UvalGrid,*binval,*binval_mod,*uval,*weight,*num,wg;
  int NUGrid,ii,jj,ii1,index;
  
  if(access(INFITS, F_OK)!=0){
    printf("Input File %s does not exists\n",INFITS);
    exit(0);
  }
  
  if(access(CONST, F_OK)!=0){
    printf("Input File %s does not exists\n",CONST);
    exit(0);
  }
  
  BREAD_HDR(INFITS);
  
  printf("dU=%lf naxes vis={%ld,%ld,%ld,%ld}\n",CDELT[0],naxes[0],naxes[1],naxes[2],naxes[3]); 
  Nu=(int) naxes[0];Nv=(int) naxes[1]; n_ave=(int) naxes[3];
  Nuv= Nu*Nv;
  printf("Nu=%d Nv=%d Nuv=%d n_ave=%d\n",Nu,Nv,Nuv,n_ave);
  
  BREAD_HDR(CONST);
  
  printf("dU=%lf naxes vis={%ld,%ld,%ld,%ld}\n",CDELT[0],naxes[0],naxes[1],naxes[2],naxes[3]); 
  Nu=(int) naxes[0];Nv=(int) naxes[1];
  printf("Nu=%d Nv=%d\n",Nu,Nv);

  GV = (double*) calloc(Nuv, sizeof(double));
  GVconst = (double*) calloc(Nuv, sizeof(double));
  binval=(double*)calloc(Nbin,sizeof(double));
  binval_mod=(double*)calloc(Nbin,sizeof(double));
  weight=(double*)calloc(Nbin,sizeof(double));
  uval=(double*)calloc(Nbin,sizeof(double));
  num=(double*)calloc(Nbin,sizeof(double));

  fits_open_file(&fptr,INFITS,READONLY,&status);
  fits_open_file(&fptr1,CONST,READONLY,&status);
  fp = fopen(argv[3],"w");  
  
  for(chan=0;chan<n_ave;++chan)
    { 
      printf("chan=%d\n",chan); 
      fpixel[0]=1;fpixel[1]=1;fpixel[2]=1;fpixel[3]=chan+1;
      lpixel[0]=(long)Nu;lpixel[1]=(long)Nv;lpixel[2]=1;lpixel[3]=chan+1;

      fits_read_subset(fptr,TDOUBLE,fpixel,lpixel,inc,&nulval,GV,&anynull,&status);
      fits_read_subset(fptr1,TDOUBLE,fpixel,lpixel,inc,&nulval,GVconst,&anynull,&status);

      for(ii=0;ii<Nbin;++ii)
      	{
      	  uval[ii]=0.;
      	  binval[ii]=0.;
      	  weight[ii]=0.;
      	}
      
      for(jj=0;jj<Nv;jj++)
      	for(ii=0;ii<Nu;ii++)
      	  {
      	    index=jj*Nu+ii;
      	    ii1=(ii>Nu/2) ? ii-Nu : ii;
      	    UvalGrid=CDELT[0]*sqrt(1.*(ii1*ii1+jj*jj));

      	    if(UvalGrid> Umin && UvalGrid<= Umax)
      	      {
		NUGrid=(int)floor(kv*(UvalGrid-Umin));
		NUGrid = (NUGrid>0) ? NUGrid:0;

		//wg=GVconst[index];
		wg=1.;
      		
		if(GVconst[index]>cutoff && NUGrid<Nbin)
		  {
		    binval[NUGrid]+=(wg*GV[index]/GVconst[index]);
		    binval_mod[NUGrid]+=P(UvalGrid);
		    uval[NUGrid]+=(UvalGrid*wg);
		    weight[NUGrid]+=wg;
		    num[NUGrid]+=1.;
		  }
      	      }
      	  }
      
      for(ii=0;ii<Nbin;++ii)
      	{
      	  if(weight[ii]>0.)
      	    fprintf(fp,"%e %e %e %e\n",uval[ii]/weight[ii],binval[ii]/weight[ii],num[ii],binval_mod[ii]/weight[ii]);
	  else
	    fprintf(fp,"%e %e %e %e\n",0.,-1.e10,0.,0.);
      	}
      fprintf(fp,"\n");
    }
  
  fits_close_file(fptr,&status);
  fclose(fp);
}
