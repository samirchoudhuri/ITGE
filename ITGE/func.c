//This programme for mask window
# include <stdio.h>
# include <fitsio.h>
#include <math.h>
# include <unistd.h>

extern double CDELT[],theta0,theta02;

void win(double *image,long xdim,long ydim,char *catlog)
{
  int ii,jj,ii1,jj1,ll,mm,nn=4,index,test=1;
  int xmin,xmax,ymin,ymax;
  double thetax,thetay,width;
  FILE *fp;

  for(jj=0;jj<ydim;++jj)
    for(ii=0;ii<xdim;++ii)
      {
	index=jj*xdim+ii;
	ll=ii-xdim/2;mm=jj-ydim/2;
	image[index]=exp(-1.*(ll*ll+mm*mm)*CDELT[0]*CDELT[0]/(theta0*theta0));
      }
  
  if(catlog==NULL)
    printf("No Input catalog File only tapering window applied\n");
  else
    {
      if(access(catlog, F_OK)!=0)
	{
	  printf("Catalog File %s does not exists\n",catlog);
	  printf("Remove created FITS Image File\n");
	  exit(0);
	}
      printf("Catalog File found, window modifying accordingly\n");

      fp=fopen(catlog,"r");
      while(test>0)
	{
	  test=fscanf(fp,"%lf%lf%*f%lf",&thetax,&thetay,&width);

	  if(test>0)
	    {
	      ii1=(int)roundf(thetax/CDELT[0]);
	      ii1+=xdim/2;
	      jj1=(int)roundf(thetay/CDELT[0]);
	      jj1+=ydim/2;
	      mm=(int)roundf(width/CDELT[0]);
	      xmin=ii1-nn*mm;xmax=ii1+nn*mm;
	      ymin=jj1-nn*mm;ymax=jj1+nn*mm;
	      printf("mask centre=[%d,%d] %dX%d window=[%d,%d %d,%d]\n",ii1,jj1,nn,mm,xmin,ymin,xmax,ymax);
	      for(jj=ymin;jj<=ymax;++jj)
		for(ii=xmin;ii<=xmax;++ii)
		  {
		    index=jj*xdim+ii;
		    image[index]*=(1.-exp(-1.*((ii-ii1)*(ii-ii1)+(jj-jj1)*(jj-jj1))*CDELT[0]*CDELT[0]/(width*width)));
		  }
	    }
	}
      fclose(fp);
    }
}

void disk(double *image,long xdim,long ydim,char *catlog)
{
  int ii,jj,ii1,jj1,ll,mm,nn=4,index,test=1;
  int xmin,xmax,ymin,ymax;
  double thetax,thetay,width,dist2;
  FILE *fp;

  for(jj=0;jj<ydim;++jj)
    for(ii=0;ii<xdim;++ii)
      {
	index=jj*xdim+ii;
	ll=ii-xdim/2;mm=jj-ydim/2;
	dist2=1.*(ll*ll+mm*mm)*CDELT[0]*CDELT[0];
	image[index]=(dist2<(theta0*theta0))?1.:0.;

	/* if(dist2<(theta0*theta0)) */
	/*   image[index]=1.; */
	/* else */
	/*   image[index]=0. ;*/
      }
  

  if(catlog==NULL)
    printf("No Input catalog File only tapering window applied\n");
  else
    {
      if(access(catlog, F_OK)!=0)
	{
	  printf("Catalog File %s does not exists\n",catlog);
	  printf("Remove created FITS Image File\n");
	  exit(0);
	}
      printf("Catalog File found, window modifying accordingly\n");

      fp=fopen(catlog,"r");
      while(test>0)
	{
	  test=fscanf(fp,"%lf%lf%*f%lf",&thetax,&thetay,&width);
	  
	  if(test>0)
	    {
	      ii1=(int)roundf(thetax/CDELT[0]);
	      ii1+=xdim/2;
	      jj1=(int)roundf(thetay/CDELT[0]);
	      jj1+=ydim/2;
	      mm=(int)roundf(width/CDELT[0]);
	      xmin=ii1-nn*mm;xmax=ii1+nn*mm;
	      ymin=jj1-nn*mm;ymax=jj1+nn*mm;
	      printf("mask centre=[%d,%d] %dX%d window=[%d,%d %d,%d]\n",ii1,jj1,nn,mm,xmin,ymin,xmax,ymax);
	      for(jj=ymin;jj<=ymax;++jj)
		for(ii=xmin;ii<=xmax;++ii)
		  {
		    index=jj*xdim+ii;
		    image[index]*=0.;
		  }
	    }
	}
      fclose(fp);
    }
  
  
}

void inver_disk(double *image,long xdim,long ydim,char *catlog)
{
  int ii,jj,ii1,jj1,ll,mm,nn=4,index,test=1;
  int xmin,xmax,ymin,ymax;
  double thetax,thetay,width,dist2;
  FILE *fp;

  for(jj=0;jj<ydim;++jj)
    for(ii=0;ii<xdim;++ii)
      {
	index=jj*xdim+ii;
	ll=ii-xdim/2;mm=jj-ydim/2;
	dist2=1.*(ll*ll+mm*mm)*CDELT[0]*CDELT[0];
	image[index]=(dist2<(theta0*theta0))?0.:1.;

	/* if(dist2<(theta0*theta0)) */
	/*   image[index]=1.; */
	/* else */
	/*   image[index]=0. ;*/
      }
  

  if(catlog==NULL)
    printf("No Input catalog File only tapering window applied\n");
  else
    {
      if(access(catlog, F_OK)!=0)
	{
	  printf("Catalog File %s does not exists\n",catlog);
	  printf("Remove created FITS Image File\n");
	  exit(0);
	}
      printf("Catalog File found, window modifying accordingly\n");

      fp=fopen(catlog,"r");
      while(test>0)
	{
	  test=fscanf(fp,"%lf%lf%*f%lf",&thetax,&thetay,&width);
	  
	  if(test>0)
	    {
	      ii1=(int)roundf(thetax/CDELT[0]);
	      ii1+=xdim/2;
	      jj1=(int)roundf(thetay/CDELT[0]);
	      jj1+=ydim/2;
	      mm=(int)roundf(width/CDELT[0]);
	      xmin=ii1-nn*mm;xmax=ii1+nn*mm;
	      ymin=jj1-nn*mm;ymax=jj1+nn*mm;
	      printf("mask centre=[%d,%d] %dX%d window=[%d,%d %d,%d]\n",ii1,jj1,nn,mm,xmin,ymin,xmax,ymax);
	      for(jj=ymin;jj<=ymax;++jj)
		for(ii=xmin;ii<=xmax;++ii)
		  {
		    index=jj*xdim+ii;
		    image[index]*=0.;
		  }
	    }
	}
      fclose(fp);
    }
  
  
}



void exp_sinc(double *image,long xdim,long ydim,char *catlog)
{
  int ii,jj,ii1,jj1,ll,mm,nn=4,index,test=1;
  int xmin,xmax,ymin,ymax;
  double thetax,thetay,width,dist2;
  double w1,w2,alpha=2.;
  FILE *fp;

  w1=theta0/CDELT[0];
  w2=theta0/(CDELT[0]*M_PI);

  for(jj=0;jj<ydim;++jj)
    for(ii=0;ii<xdim;++ii)
      {
	index=jj*xdim+ii;
	ll=ii-xdim/2;mm=jj-ydim/2;
	dist2=1.*(ll*ll+mm*mm)*CDELT[0]*CDELT[0];
	//image[index]=(dist2<(theta0*theta0))?1.:0.;

	if(dist2==0.)
	  image[index]=1.;

	if((0<dist2) && (dist2<(theta0*theta0)))
	  image[index]=exp(-1.*dist2/(w1*w1*CDELT[0]*CDELT[0]))*sin(sqrt(dist2)/(w2*CDELT[0]))/(sqrt(dist2)/(w2*CDELT[0]));
	else
	  image[index]=0. ;
      }
  

  if(catlog==NULL)
    printf("No Input catalog File only tapering window applied\n");
  else
    {
      if(access(catlog, F_OK)!=0)
	{
	  printf("Catalog File %s does not exists\n",catlog);
	  printf("Remove created FITS Image File\n");
	  exit(0);
	}
      printf("Catalog File found, window modifying accordingly\n");

      fp=fopen(catlog,"r");
      while(test>0)
	{
	  test=fscanf(fp,"%lf%lf%*f%lf",&thetax,&thetay,&width);
	  
	  if(test>0)
	    {
	      ii1=(int)roundf(thetax/CDELT[0]);
	      ii1+=xdim/2;
	      jj1=(int)roundf(thetay/CDELT[0]);
	      jj1+=ydim/2;
	      mm=(int)roundf(width/CDELT[0]);
	      xmin=ii1-nn*mm;xmax=ii1+nn*mm;
	      ymin=jj1-nn*mm;ymax=jj1+nn*mm;
	      printf("mask centre=[%d,%d] %dX%d window=[%d,%d %d,%d]\n",ii1,jj1,nn,mm,xmin,ymin,xmax,ymax);
	      for(jj=ymin;jj<=ymax;++jj)
		for(ii=xmin;ii<=xmax;++ii)
		  {
		    index=jj*xdim+ii;
		    image[index]*=0.;
		  }
	    }
	}
      fclose(fp);
    }
  
  
}


void radial(double *image,long xdim,long ydim,char *catlog)
{
  int ii,jj,ii1,jj1,ll,mm,nn=4,index,test=1;
  int xmin,xmax,ymin,ymax;
  double thetax,thetay,width,dist2,n=2.;
  FILE *fp;

  for(jj=0;jj<ydim;++jj)
    for(ii=0;ii<xdim;++ii)
      {
	index=jj*xdim+ii;
	ll=ii-xdim/2;mm=jj-ydim/2;
	dist2=1.*(ll*ll+mm*mm)*CDELT[0]*CDELT[0];
	image[index]=(dist2<(theta0*theta0))?(1.-pow((sqrt(dist2)/theta0),n)):0.;

	/* if(dist2<(theta0*theta0)) */
	/*   image[index]=1.; */
	/* else */
	/*   image[index]=0. ;*/
      }
  

  if(catlog==NULL)
    printf("No Input catalog File only tapering window applied\n");
  else
    {
      if(access(catlog, F_OK)!=0)
	{
	  printf("Catalog File %s does not exists\n",catlog);
	  printf("Remove created FITS Image File\n");
	  exit(0);
	}
      printf("Catalog File found, window modifying accordingly\n");

      fp=fopen(catlog,"r");
      while(test>0)
	{
	  test=fscanf(fp,"%lf%lf%*f%lf",&thetax,&thetay,&width);
	  
	  if(test>0)
	    {
	      ii1=(int)roundf(thetax/CDELT[0]);
	      ii1+=xdim/2;
	      jj1=(int)roundf(thetay/CDELT[0]);
	      jj1+=ydim/2;
	      mm=(int)roundf(width/CDELT[0]);
	      xmin=ii1-nn*mm;xmax=ii1+nn*mm;
	      ymin=jj1-nn*mm;ymax=jj1+nn*mm;
	      printf("mask centre=[%d,%d] %dX%d window=[%d,%d %d,%d]\n",ii1,jj1,nn,mm,xmin,ymin,xmax,ymax);
	      for(jj=ymin;jj<=ymax;++jj)
		for(ii=xmin;ii<=xmax;++ii)
		  {
		    index=jj*xdim+ii;
		    image[index]*=0.;
		  }
	    }
	}
      fclose(fp);
    }
  
  
}


void butterworth(double *image,long xdim,long ydim,char *catlog)
{
  int ii,jj,ii1,jj1,ll,mm,nn=4,index,test=1;
  int xmin,xmax,ymin,ymax;
  double thetax,thetay,width,dist2,n=4.;
  FILE *fp;

  for(jj=0;jj<ydim;++jj)
    for(ii=0;ii<xdim;++ii)
      {
	index=jj*xdim+ii;
	ll=ii-xdim/2;mm=jj-ydim/2;
	dist2=1.*(ll*ll+mm*mm)*CDELT[0]*CDELT[0];
	image[index]=1./(1.+pow((sqrt(dist2)/theta0),2.*n));

	/* if(dist2<(theta0*theta0)) */
	/*   image[index]=1.; */
	/* else */
	/*   image[index]=0. ;*/
      }
  

  if(catlog==NULL)
    printf("No Input catalog File only tapering window applied\n");
  else
    {
      if(access(catlog, F_OK)!=0)
	{
	  printf("Catalog File %s does not exists\n",catlog);
	  printf("Remove created FITS Image File\n");
	  exit(0);
	}
      printf("Catalog File found, window modifying accordingly\n");

      fp=fopen(catlog,"r");
      while(test>0)
	{
	  test=fscanf(fp,"%lf%lf%*f%lf",&thetax,&thetay,&width);
	  
	  if(test>0)
	    {
	      ii1=(int)roundf(thetax/CDELT[0]);
	      ii1+=xdim/2;
	      jj1=(int)roundf(thetay/CDELT[0]);
	      jj1+=ydim/2;
	      mm=(int)roundf(width/CDELT[0]);
	      xmin=ii1-nn*mm;xmax=ii1+nn*mm;
	      ymin=jj1-nn*mm;ymax=jj1+nn*mm;
	      printf("mask centre=[%d,%d] %dX%d window=[%d,%d %d,%d]\n",ii1,jj1,nn,mm,xmin,ymin,xmax,ymax);
	      for(jj=ymin;jj<=ymax;++jj)
		for(ii=xmin;ii<=xmax;++ii)
		  {
		    index=jj*xdim+ii;
		    image[index]*=(1.-(1./(1.+pow((sqrt(1.*((ii-ii1)*(ii-ii1)+(jj-jj1)*(jj-jj1))*CDELT[0]*CDELT[0])/width),2.*n))));
		    //image[index]*=0.;
		  }
	    }
	}
      fclose(fp);
    }
}


void butterworthannular(double *image,long xdim,long ydim,char *catlog)
{
  int ii,jj,ii1,jj1,ll,mm,nn=4,index,test=1;
  int xmin,xmax,ymin,ymax;
  double thetax,thetay,width,dist2,n=4.;
  FILE *fp;

  for(jj=0;jj<ydim;++jj)
    for(ii=0;ii<xdim;++ii)
      {
	index=jj*xdim+ii;
	ll=ii-xdim/2;mm=jj-ydim/2;
	dist2=1.*(ll*ll+mm*mm)*CDELT[0]*CDELT[0];
	image[index]=(1./(1.+pow((sqrt(dist2)/theta02),2.*n)))-(1./(1.+pow((sqrt(dist2)/theta0),2.*n)));

	/* if(dist2<(theta0*theta0)) */
	/*   image[index]=1.; */
	/* else */
	/*   image[index]=0. ;*/
      }
  

  if(catlog==NULL)
    printf("No Input catalog File only tapering window applied\n");
  else
    {
      if(access(catlog, F_OK)!=0)
	{
	  printf("Catalog File %s does not exists\n",catlog);
	  printf("Remove created FITS Image File\n");
	  exit(0);
	}
      printf("Catalog File found, window modifying accordingly\n");

      fp=fopen(catlog,"r");
      while(test>0)
	{
	  test=fscanf(fp,"%lf%lf%*f%lf",&thetax,&thetay,&width);
	  
	  if(test>0)
	    {
	      ii1=(int)roundf(thetax/CDELT[0]);
	      ii1+=xdim/2;
	      jj1=(int)roundf(thetay/CDELT[0]);
	      jj1+=ydim/2;
	      mm=(int)roundf(width/CDELT[0]);
	      xmin=ii1-nn*mm;xmax=ii1+nn*mm;
	      ymin=jj1-nn*mm;ymax=jj1+nn*mm;
	      printf("mask centre=[%d,%d] %dX%d window=[%d,%d %d,%d]\n",ii1,jj1,nn,mm,xmin,ymin,xmax,ymax);
	      for(jj=ymin;jj<=ymax;++jj)
		for(ii=xmin;ii<=xmax;++ii)
		  {
		    index=jj*xdim+ii;
		    image[index]*=0.;
		  }
	    }
	}
      fclose(fp);
    }
  
  
}




void disk_gauss(double *image,long xdim,long ydim,char *catlog)
{
  int ii,jj,ii1,jj1,ll,mm,nn=4,index,test=1;
  int xmin,xmax,ymin,ymax;
  double thetax,thetay,width,dist2,n=2.;
  FILE *fp;

  for(jj=0;jj<ydim;++jj)
    for(ii=0;ii<xdim;++ii)
      {
	index=jj*xdim+ii;
	ll=ii-xdim/2;mm=jj-ydim/2;
	dist2=1.*(ll*ll+mm*mm)*CDELT[0]*CDELT[0];
	//image[index]=(dist2<(theta0*theta0))?1:exp((-1.*(sqrt(dist2)-theta0)*(sqrt(dist2)-theta0))/(theta0*theta0));
	image[index]=(dist2>(theta0*theta0))?1:exp((-1.*(sqrt(dist2)-theta0)*(sqrt(dist2)-theta0))/(theta0*theta0));

	/* if(dist2<(theta0*theta0)) */
	/*   image[index]=1.; */
	/* else */
	/*   image[index]=0. ;*/
      }
  

  if(catlog==NULL)
    printf("No Input catalog File only tapering window applied\n");
  else
    {
      if(access(catlog, F_OK)!=0)
	{
	  printf("Catalog File %s does not exists\n",catlog);
	  printf("Remove created FITS Image File\n");
	  exit(0);
	}
      printf("Catalog File found, window modifying accordingly\n");

      fp=fopen(catlog,"r");
      while(test>0)
	{
	  test=fscanf(fp,"%lf%lf%*f%lf",&thetax,&thetay,&width);
	  
	  if(test>0)
	    {
	      ii1=(int)roundf(thetax/CDELT[0]);
	      ii1+=xdim/2;
	      jj1=(int)roundf(thetay/CDELT[0]);
	      jj1+=ydim/2;
	      mm=(int)roundf(width/CDELT[0]);
	      xmin=ii1-nn*mm;xmax=ii1+nn*mm;
	      ymin=jj1-nn*mm;ymax=jj1+nn*mm;
	      printf("mask centre=[%d,%d] %dX%d window=[%d,%d %d,%d]\n",ii1,jj1,nn,mm,xmin,ymin,xmax,ymax);
	      for(jj=ymin;jj<=ymax;++jj)
		for(ii=xmin;ii<=xmax;++ii)
		  {
		    index=jj*xdim+ii;
		    image[index]*=0.;
		  }
	    }
	}
      fclose(fp);
    }
  
  
}
