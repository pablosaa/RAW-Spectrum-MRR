/* This is a subroutine to read the MicroRainRadar (MRR) RAW data files by METEK
   This code is called by the main MRR4ADMI Fortran code.
   Part of the RAW-MRR-spectrum repository.
   
   Copyright 2012-2014 Pablo Saavedra Garfias (pablosaa@uni-bonn.de)
   See LICENSE.TXT
*/

#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>
#include<stdbool.h>

/* Following extern variables are defined at FORTRAN module: variables.mod */
extern int Nh, Nspec, Nt;
extern int PFIELD, PHEAD, NLINE;

/* INPUT VARIABLES:
   filename: (string) containing absolute path for the RAW data file,
   OUTPUT VARIABLES:
   Ntot: (int) total number of profiles,
   time: (double) time in hr corresponding to each profile,
   timeheight:
   Fnn: (integer) 3D array containing number of profiles, ranges, spectrum,
   TF: (double) Transfer function for the range profile,
   CC: (integer) instrument specific Calibration Constant
*/

void read_raw_cfile_(char *fname, int *Ndat, float *t_hr, int *height, int Fnn[Nt][Nspec][Nh], float *TF, int *CC, int dates[2][6], int *exitcode){

  FILE *fp;
  int i,j,k;
  int idat=0, status;
  int yy,mm,dd,hr,mi,se,dsn,bw, mdq[3];
  float dvs;
  char key[]="xx", spkey[]="xxx", buffer[NLINE+1];
  char head1[]="T:%2d%2d%2d%2d%2d%2d UTC DVS %4f DSN %6d CC %7d\n";
  char head2[]="MRR %2d%2d%2d%2d%2d%2d UTC DVS %4f DSN %9d BW %5d CC %7d MDQ %3d %2d %2d TYP RAW\n";
  bool OLDVER=false, NEWDATA=false;

  printf("Working at %s\n",fname);
  /* Open the RAW data file */
  fp = fopen(fname, "r");
  if(fp==NULL){
    printf("Input file cannot be open!\n");
    *exitcode = -1;
    return;
  }

  /* Finding out the time stamp line: */
  while(!feof(fp)){ /* !feof(fp)){ */
    fscanf(fp,"%2c",key);
    
    if(!strcmp(key,"T:") || !strcmp(key,"MR")) {fseek(fp,-2,SEEK_CUR);  NEWDATA=true;}
    else {fgets(buffer,NLINE+1,fp);
      /* printf("ERROR: Time stamp not found! %s\n",buffer); */
      continue;}

    /* Here starts the reading process to feed the variables */
    if(fgets(buffer,NLINE+1,fp)==NULL){
      printf("ERROR: First data file line not possible to read.\n");
      *exitcode = -1;
      return;
    }

    /* Finding out the Firmware version: */
    status = sscanf(buffer,head1,&yy,&mm,&dd,&hr,&mi,&se,&dvs,&dsn,CC);
    if(status==9 && dvs<4.5) OLDVER = true;
    else {
      status = sscanf(buffer,head2,&yy,&mm,&dd,&hr,&mi,&se,&dvs,&dsn,&bw,CC,&mdq[0],&mdq[1],&mdq[2]); 
      if(status==13 && dvs>=6) OLDVER=false;
      else {
	printf("ERROR: Firmware version not identified\n");
	*exitcode = -1;
	return;}
    }
    /* Feeding the hour in the time variable */
    if(idat==0) {
      dates[0][0]=dd;dates[0][1]=mm;dates[0][2]=yy;
      dates[0][3]=hr;dates[0][4]=mi;dates[0][5]=se;}
    t_hr[idat] = (float) hr + ((float) mi)/60 + ((float) se)/3600;

    /* 1.- Reading the RANGES for the radar profile */
    status = fscanf(fp,OLDVER?"M:%c =":"%c%*2c",key);
    key[1] = '\0';
    if(status!=1 && strcmp(key,OLDVER?"h":"H")) printf("ERROR: Range H header not found.\n");

    if(fgets(buffer,NLINE+1,fp)==NULL) printf("ERROR: bad buffer %s\n",buffer);

    for(i=0;i<Nh;++i){
      status = sscanf(buffer+i*PFIELD,"%d[0123456789]", &height[i]);
      if(status!=1)  printf("Status: %d %d\n",i, status);
    }

    /* 2.- Reading the TRANSFER FUNCTION profile*/
    buffer[0] = 0;
    status = fscanf(fp,OLDVER?"M:%02c =":"%02c ",key);
    if(status!=1 && strcmp(key,"TF")) printf("ERROR: not able to get M:TF header\n");
    if(fgets(buffer,NLINE+1,fp)==NULL)  printf("ERROR: bad buffer? %s\n",buffer);

    for(i=0;i<Nh;++i){
      status = sscanf(buffer+i*PFIELD,"%f[.0123456789]", &TF[i]);
      if(status!=1) printf("Status: %d %d\n",i, status);
    }

    /* 3.- Reading the 64 spectral profiles for every time step */
    j=0;
    while(j<Nspec){
      buffer[0]='\0';
      spkey[0]='*'; spkey[1]='\0';
      for(i=0;i<Nh;++i) Fnn[idat][j][i]=(int) NAN; /*Fnn[i+j*Nh+idat*(Nspec*Nh)] = (int) NAN; First filling profile with NAN*/
      status = fscanf(fp,OLDVER?"%3c":"%c",spkey);
      if(status!=1 || strcmp(spkey,OLDVER?"M:f":"F")) {
	printf("ERROR: (status=%d, key=%s) at %d %d %d. Spectral bin should:%d, but read:%d\n",status,spkey,hr,mi,se,j,k);
	/* Checking what kind of bad data has been found */
	spkey[OLDVER?2:1]='\0';
	if(!strcmp(spkey,OLDVER?"T:":"M")){
	  /* Meaning the spectrum is incomplete in the data file and new Time stamp has beed found */
	  fseek(fp,OLDVER?-3:-1,SEEK_CUR);
	  j++;
	  continue;
	}
	fgets(buffer,NLINE+1,fp);  /* retrieving whatever is left in the stream */ 	
	printf("Buffer within spectrum: %s\n", buffer);
	continue;
      }

      status = fscanf(fp,"%2d[0123456789]",&k);
      if(status!=1 || k<0 || k>=Nspec){
	fgets(buffer,NLINE+1,fp);  /* retrieving whatever is left in the stream */ 	
	printf("Buffer within spectrum: %s\n", buffer);
	continue;
      }
      
      if(fgets(buffer,NLINE+1,fp)==NULL)  printf("%s\n",buffer);
      /* Feeding the output variable Fnn with the power spectrum profile*/
      for(i=0;i<Nh;++i){
	status = sscanf(buffer+i*PFIELD,"%9d[0123456789]", &Fnn[idat][k][i]); /* &Fnn[i+k*Nh+idat*(Nspec*Nh)]); */
	/* if(status!=1) Fnn[i+k*Nh+idat*(Nspec*Nh)] = (int) NAN; */
      }
      ++j;
    } /* end over the spectrum feeding loop */
    ++idat;   /* getting ready for next data profile */

  }  /* end over the EOF while */

  fclose(fp);
  if(!NEWDATA) {
    printf("ERROR: Seems the file is empty or not an MRR RAW data file!.\n");
    *exitcode = -1;
    return;
  }
  if(idat!=0) {
    dates[1][0]=dd;dates[1][1]=mm;dates[1][2]=yy;
    dates[1][3]=hr;dates[1][4]=mi;dates[1][5]=se;}
  *Ndat = idat;
  *exitcode = 0;
  return;

  /* End of the MRR RAW data file procedure */
}


