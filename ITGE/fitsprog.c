//This programme for FITS HEADER
# include <stdio.h>
# include <fitsio.h>
#include <math.h>

extern double CRVAL[],CDELT[],CRPIX[];
extern long naxes[],naxesim[];

int naxis=4,status=0;

void printerror(int status){
  if (status){
    fits_report_error(stderr, status);
    exit( status ); 
  }
}

void BWRITE_HDR(char *outfile)
{
  fitsfile *fitsp;
  int ii;
  char *CTYPE[4] = {"Nu", "Nv", "RE_IM3", "CHAN"};
  char keynam[6];

  fits_create_file(&fitsp, outfile, &status);
  fits_create_img(fitsp, DOUBLE_IMG, naxis, naxes, &status);
  
  for(ii=0; ii<naxis; ii++){
    
    if(fits_make_keyn("CTYPE", ii+1, keynam, &status))
      printerror(status);
    if(fits_write_key_str(fitsp, keynam, CTYPE[ii], " axis type", &status))
      printerror(status);
    
    if(fits_make_keyn("CRVAL", ii+1, keynam, &status))
      printerror(status);
    if(fits_write_key_dbl(fitsp, keynam, CRVAL[ii], 10, " ", &status))
      printerror(status);
    
    if(fits_make_keyn("CDELT", ii+1, keynam, &status))
      printerror(status);
    if(fits_write_key_dbl(fitsp, keynam, CDELT[ii],  9, " ", &status))
      printerror(status);
     
    if(fits_make_keyn("CRPIX", ii+1, keynam, &status))
      printerror(status);
    if(fits_write_key_dbl(fitsp, keynam, CRPIX[ii],  9, " ", &status))
      printerror(status);
  }
  fits_close_file(fitsp, &status);
}


void BWRITEtau_HDR(char *outfile)
{
  fitsfile *fitsp;
  int ii;
  char *CTYPE[4] = {"Nu", "Nv", "RE_IM3", "tau"};
  char keynam[6];

  fits_create_file(&fitsp, outfile, &status);
  fits_create_img(fitsp, DOUBLE_IMG, naxis, naxes, &status);
  
  for(ii=0; ii<naxis; ii++){
    
    if(fits_make_keyn("CTYPE", ii+1, keynam, &status))
      printerror(status);
    if(fits_write_key_str(fitsp, keynam, CTYPE[ii], " axis type", &status))
      printerror(status);
    
    if(fits_make_keyn("CRVAL", ii+1, keynam, &status))
      printerror(status);
    if(fits_write_key_dbl(fitsp, keynam, CRVAL[ii], 10, " ", &status))
      printerror(status);
    
    if(fits_make_keyn("CDELT", ii+1, keynam, &status))
      printerror(status);
    if(fits_write_key_dbl(fitsp, keynam, CDELT[ii],  9, " ", &status))
      printerror(status);
     
    if(fits_make_keyn("CRPIX", ii+1, keynam, &status))
      printerror(status);
    if(fits_write_key_dbl(fitsp, keynam, CRPIX[ii],  9, " ", &status))
      printerror(status);
  }
  fits_close_file(fitsp, &status);
}

void BREAD_HDR(char *infile)
{
  fitsfile *fitsp;
  int nfound;
  fits_open_file(&fitsp,infile,READONLY,&status);
  fits_read_keys_lng(fitsp,"NAXIS",1,naxis,naxes,&nfound,&status);
  fits_read_keys_dbl(fitsp,"CDELT",1,naxis,CDELT,&nfound,&status);
  fits_read_keys_dbl(fitsp,"CRVAL",1,naxis,CRVAL,&nfound,&status);
  fits_read_keys_dbl(fitsp,"CRPIX",1,naxis,CRPIX,&nfound,&status);
  fits_close_file(fitsp,&status);
}

void SWRITE_HDR(char *outfile)
{
  fitsfile *fitsp;
  int ii;
  char *CTYPE[4] = {"thetax", "thetay", "3IMAGE", "CHAN"};
  char keynam[6];

  fits_create_file(&fitsp, outfile, &status);
  fits_create_img(fitsp, DOUBLE_IMG, naxis, naxesim, &status);
  
  for(ii=0; ii<naxis; ii++){
    
    if(fits_make_keyn("CTYPE", ii+1, keynam, &status))
      printerror(status);
    if(fits_write_key_str(fitsp, keynam, CTYPE[ii], " axis type", &status))
      printerror(status);
    
    if(fits_make_keyn("CRVAL", ii+1, keynam, &status))
      printerror(status);
    if(fits_write_key_dbl(fitsp, keynam, CRVAL[ii], 10, " ", &status))
      printerror(status);
    
    if(fits_make_keyn("CDELT", ii+1, keynam, &status))
      printerror(status);
    if(fits_write_key_dbl(fitsp, keynam, CDELT[ii],  9, " ", &status))
      printerror(status);
     
    if(fits_make_keyn("CRPIX", ii+1, keynam, &status))
      printerror(status);
    if(fits_write_key_dbl(fitsp, keynam, CRPIX[ii],  9, " ", &status))
      printerror(status);
  }
  fits_close_file(fitsp, &status);
}

void SREAD_HDR(char *infile)
{
  fitsfile *fitsp;
  int nfound;
  fits_open_file(&fitsp,infile,READONLY,&status);
  fits_read_keys_lng(fitsp,"NAXIS",1,naxis,naxesim,&nfound,&status);
  fits_read_keys_dbl(fitsp,"CRVAL",1,naxis,CRVAL,&nfound,&status);
  fits_read_keys_dbl(fitsp,"CDELT",1,naxis,CDELT,&nfound,&status);
  fits_read_keys_dbl(fitsp,"CRPIX",1,naxis,CRPIX,&nfound,&status);
  fits_close_file(fitsp,&status);
}

