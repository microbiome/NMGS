/*System includes*/
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <sys/stat.h>

/*GSL includes*/
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_sf.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_roots.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_statistics_double.h>
#include <gsl/gsl_fft_complex.h>
#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_multimin.h>
#include <gsl/gsl_deriv.h>
#include <gsl/gsl_randist.h>
#include <pthread.h>

/*User includes*/
#include "NMGS.h"

static char *usage[] = {"NMGS - Fits multisite HDP neutral model to a matrix of abundances\n",
                        "Required parameters:\n",
			"\t-in\tfilename\tcsv file\n",
			"\t-out\toutfilestub\n",
                        "Optional:\n",	
			"\t-b\tinteger\tnumber of burn iterations\n",
			"\t-e\tinteger\textrapolate size\n",
			"\t-l\tinteger\tseed\n",
			"\t-o\t\toutput samples\n",
			"\t-s\t\tbootstrap for neutral fit\n",
			"\t-t\tintger\tnumber of iterations\n",
			"\t-v\t\tverbose\n"};

static int nLines   = 13;

static int bVerbose = FALSE;

int main(int argc, char* argv[])
{
  int  ii = 0, i = 0, j = 0, jj = 0, k = 0, nN = 0, nS = 0;
  int nMaxX = 0, nIter = 0, nMaxIter = 0;
  t_Params tParams;
  t_Data   tData;
  gsl_rng            *ptGSLRNG     = NULL;
  const gsl_rng_type *ptGSLRNGType = NULL;
  double             *adLogNormalisation = NULL, *adLogProbV = NULL, *adLogProbAdjustV, *adProbV, *adCProbV; 
  double             **aadStirlingMatrix = NULL, *adStirlingVector = NULL;
  double             *adJ = NULL, *adM = NULL, **aadMStore = NULL, **aadIStore = NULL, *adThetaStore = NULL;
  int                **aanTStore = NULL;
  double             *adW = NULL, *adS = NULL;
  int                *anTC = NULL,*anTR = NULL, **aanT = NULL, **aanX = NULL, *anJ = NULL;
  double             dLogNormalisation = 0.0, dTheta = 0.0;
  char*              szSampleDir = (char *) malloc(MAX_LINE_LENGTH*sizeof(char));
  /*initialise GSL RNG*/
  gsl_rng_env_setup();

  gsl_set_error_handler_off();
  
  ptGSLRNGType = gsl_rng_default;
  ptGSLRNG     = gsl_rng_alloc(ptGSLRNGType);

  gsl_set_error_handler_off();
  
  /*get command line params*/
  getCommandLineParams(&tParams, argc, argv);
  nMaxIter = tParams.nMaxIter;
  /*set RNG seed*/
  gsl_rng_set (ptGSLRNG, tParams.nL);
  	
  /*read in abundance distribution*/
  readAbundanceData(tParams.szInputFile,&tData);
  
  nS = tData.nS;
  nN = tData.nN;

  nMaxX = maxX(&tData);

  if(tParams.nRarefy != -1){
    int nMinJ = minJ(&tData);
    t_Data tRData;
    char szRFile[MAX_LINE_LENGTH];
    
    if(bVerbose){
      printf("Rarefy to %d:\n",nMinJ);
      fflush(stderr);
    }

    rarefy(ptGSLRNG,tParams.nRarefy,&tRData,&tData);

    sprintf(szRFile,"%sR%s",tParams.szOutFileStub,CSV_FILE_SUFFIX);
    writeAbundanceData(szRFile, &tRData);
   
    free(tData.aanX);
    free(tData.aszSampleNames);
    free(tData.aszOTUNames);

    tData.nN = tRData.nN;
    tData.nS = tRData.nS;
    tData.nSize = tRData.nSize;
    tData.aanX = tRData.aanX;
    tData.aszSampleNames = tRData.aszSampleNames;
    tData.aszOTUNames = tRData.aszOTUNames;

    nS = tRData.nS;
  }
	
  aanT = (int **) malloc(nN*sizeof(int*));
  if(!aanT)
    goto memoryError;
	
  anJ = (int *) malloc(nN*sizeof(int));
  if(!anJ)
    goto memoryError;
	

  adJ = (double *) malloc(nN*sizeof(double));
  if(!adJ)
    goto memoryError;

  aanX = tData.aanX;
  for(i = 0; i < nN; i++){
    anJ[i] = 0;

    aanT[i] = (int *) malloc(nS*sizeof(int));
    if(!aanT[i])
      goto memoryError;

    for(j = 0; j < nS; j++){
      anJ[i] += aanX[i][j];
      /*initialise with one ancestor if species present*/
      if(aanX[i][j] > 0){
	aanT[i][j] = 1;
      }
      else{
	aanT[i][j] = 0;
      }
      
    }
    adJ[i] = (double) anJ[i];
    
  }
	
  //generateStirlingDMatrix(aadStirlingMatrix, nMaxX);
  Stirling(&aadStirlingMatrix,&adLogNormalisation,nMaxX);
	
  adStirlingVector = (double *) malloc((nMaxX + 1)*sizeof(double));
  if(!adStirlingVector)
    goto memoryError;

  adLogProbV = (double *) malloc((nMaxX + 1)*sizeof(double));
  if(!adLogProbV)
    goto memoryError;

  adLogProbAdjustV = (double *) malloc((nMaxX + 1)*sizeof(double));
  if(!adLogProbAdjustV)
    goto memoryError;

  adProbV = (double *) malloc((nMaxX + 1)*sizeof(double));
  if(!adProbV)
    goto memoryError;

  adCProbV = (double *) malloc((nMaxX + 1)*sizeof(double));
  if(!adCProbV)
    goto memoryError;

  adM = (double *) malloc((nS + 1)*sizeof(double));
  if(!adM)
    goto memoryError;

  aadMStore = (double **) malloc(nMaxIter*sizeof(double*));
  if(!aadMStore)
    goto memoryError;

  aadIStore = (double **) malloc(nMaxIter*sizeof(double*));
  if(!aadIStore)
    goto memoryError;

  adThetaStore = (double *) malloc(nMaxIter*sizeof(double));
  if(!adThetaStore)
    goto memoryError;

  aanTStore = (int **) malloc(nMaxIter*sizeof(int*));
  if(!aanTStore)
    goto memoryError;

  for(i = 0; i < nMaxIter; i++){
    aanTStore[i] = (int *) malloc(nS*sizeof(int));
    if(!aanTStore[i])
      goto memoryError;

    aadMStore[i] = (double *) malloc((nS + 1)*sizeof(double));
    if(!aadMStore[i])
      goto memoryError;

    aadIStore[i] = (double *) malloc(nN*sizeof(double));
    if(!aadIStore[i])
      goto memoryError;
  }
	

  anTC = (int *) malloc(nS*sizeof(int));
  if(!anTC)
    goto memoryError;

  anTR = (int *) malloc(nN*sizeof(int));
  if(!anTR)
    goto memoryError;

  adW = (double *) malloc(nN*sizeof(double));
  if(!adW)
    goto memoryError;

  adS = (double *) malloc(nN*sizeof(double));
  if(!adS)
    goto memoryError;

  while(nIter < nMaxIter){
    sumColumns(anTC, aanT, nN, nS);
    sumRows(anTR, aanT, nN, nS);

    /*set current theta*/
    if(nIter > 0){
      dTheta = adThetaStore[nIter - 1]; 
    }
    else{
      dTheta = THETA_INIT;
    }
	
    /*sample metacommunity*/    
    sampleMetacommunity(nIter, ptGSLRNG, adM, aadMStore, dTheta, nS, anTC);
    
    sampleImmigrationRates(nIter, ptGSLRNG, aadIStore, nN, adW, adS, adJ, anTR);
    
    adThetaStore[nIter] = sampleTheta(ptGSLRNG, dTheta, anTR, nN, nS);

    sampleT(nIter, ptGSLRNG, aanX, aanT, aadIStore, aadMStore, aadStirlingMatrix, nN, nS, adLogProbV, adProbV, adCProbV);

    for(i = 0; i < nS; i++){
      aanTStore[nIter][i] = anTC[i];
    }

    if(bVerbose == TRUE){
      printf("Iteration:%d\n",nIter + 1);
      printf("Theta:%f\n",adThetaStore[nIter]);
      printf("Immigration vector:");
      for(i = 0; i < nN; i++){
	printf("%f ",aadIStore[nIter][i]);
      }
      printf("\n");
    }

    nIter++;

  }

  if(1){
    FILE *ofp = NULL;
    char szOutFile[MAX_LINE_LENGTH];

    sprintf(szOutFile,"%s%s",tParams.szOutFileStub,CSV_FILE_SUFFIX);
  
    ofp = fopen(szOutFile, "w");

    if(ofp){
      for(i = 0; i < nMaxIter; i++){
	fprintf(ofp,"%d,%f,",i,adThetaStore[i]);
	
	for(j = 0; j < nN - 1; j++){
	  fprintf(ofp,"%f,",aadIStore[i][j]);
	}

	fprintf(ofp,"%f\n",aadIStore[i][nN - 1]);
      }
      
      fclose(ofp);
    }
    else{
      fprintf(stderr, "Failed to open %s for writing\n",szOutFile);
      fflush(stderr);
    }
  }

  if(tParams.bOutputSample){
    sprintf(szSampleDir,"%s%s",tParams.szOutFileStub,SAMPLE_DIR);

    mkdir(szSampleDir, S_IRWXU);
  }

  if(tParams.nExtrapolate > -1){
    //  extrapolateSamples(nIter, nMaxIter, ptGSLRNG, nN, nS, &tParams, &tData,
    //	       adThetaStore, anJ, aadIStore, aadMStore);
    extrapolateSamples2(szSampleDir,nIter, nMaxIter, ptGSLRNG, tParams.nExtrapolate, nN, nS, &tParams, &tData,
			 adThetaStore, anJ, aadIStore, aanTStore);
  }

  if(tParams.bSample){
    printf("Sampling fit...\n");

    outputSamples(szSampleDir,nIter, nMaxIter, ptGSLRNG, nN, nS, &tParams, &tData, adThetaStore, anJ, aadIStore, aadMStore);
  }
  /*free up memory*/
  free(adW);
  adW = NULL;

  free(adS);
  adS = NULL;
	
  free(anJ);
  anJ = NULL;

  free(adJ);
  adJ = NULL;

  for(i = 0; i < nMaxIter; i++){
    free(aadMStore[i]);
    aadMStore[i] = NULL;
    free(aadIStore[i]);
    aadIStore[i] = NULL;
    free(aanTStore[i]);
    aanTStore[i] = NULL;
  }
  free(aanTStore);
  aanTStore = NULL;

  free(aadMStore);
  aadMStore = NULL;
  
  free(aadIStore);
  aadIStore = NULL;

  free(adM);
  adM = NULL;
  
  free(anTC); free(anTR);
  anTC = NULL; anTR = NULL;
  
  for(i = 0; i < nMaxX + 1; i++){
    free(aadStirlingMatrix[i]);
    aadStirlingMatrix[i] = NULL;
  }
  free(aadStirlingMatrix);
  aadStirlingMatrix = NULL;

  free(adStirlingVector);
  adStirlingVector = NULL;

  free(adLogProbV);
  adLogProbV = NULL;

  free(adLogProbAdjustV);
  adLogProbAdjustV = NULL;

  free(adProbV);
  adProbV = NULL;

  free(adCProbV);
  adCProbV = NULL;

  free(adLogNormalisation);
  adLogNormalisation = NULL;

  free(adThetaStore);
  adThetaStore = NULL;

  free(szSampleDir);						   
  szSampleDir = NULL;
  
  for(i = 0; i < nN; i++){
    free(aanT[i]);
    aanT[i] = NULL;
  }
  free(aanT);
  aanT = NULL;

  for(i = 0; i < nN; i++){
    free(tData.aanX[i]);
  }
  free(tData.aanX);

  exit(EXIT_SUCCESS);

 memoryError:
  fprintf(stderr, "Failed allocating memory in main\n");
  fflush(stderr);
  exit(EXIT_FAILURE);
}

void writeUsage(FILE* ofp)
{
  int i = 0;
  char *line;

  for(i = 0; i < nLines; i++){
    line = usage[i];
    fputs(line,ofp);
  }
}

char *extractParameter(int argc, char **argv, char *param,int when)
{
  int i = 0;

  while((i < argc) && (strcmp(param,argv[i]))){
    i++;
  }

  if(i < argc - 1){
    return(argv[i + 1]);
  }

  if((i == argc - 1) && (when == OPTION)){
    return "";
  }

  if(when == ALWAYS){
    fprintf(stdout,"Can't find asked option %s\n",param);
  }

  return (char *) NULL;
}

void getCommandLineParams(t_Params *ptParams,int argc,char *argv[])
{
  char *szTemp = NULL;
  char *cError = NULL;

  /*get parameter file name*/
  ptParams->szInputFile  = extractParameter(argc,argv, INPUT_FILE,ALWAYS);  
  if(ptParams->szInputFile == NULL)
    goto error;
 
  /*get parameter file name*/
  ptParams->szOutFileStub  = extractParameter(argc,argv,OUT_FILE_STUB,ALWAYS);  
  if(ptParams->szOutFileStub == NULL)
    goto error;

  szTemp = extractParameter(argc,argv,VERBOSE,OPTION);
  if(szTemp != NULL){
    bVerbose=TRUE;
  }

  szTemp = extractParameter(argc,argv,N_BURN_ITER,OPTION);
  if(szTemp != NULL){
    ptParams->nBurnIter = strtol(szTemp,&cError,10);
    if(*cError != '\0'){
      goto error;
    }
  }
  else{
    ptParams->nBurnIter = DEF_BURN_ITER;
  }

  szTemp = extractParameter(argc,argv,RAREFY,OPTION);
  if(szTemp != NULL){
    ptParams->nRarefy = strtol(szTemp,&cError,10);
    if(*cError != '\0'){
      goto error;
    }
  }
  else{
    ptParams->nRarefy = -1;
  }

  szTemp = extractParameter(argc,argv,OUTPUT_SAMPLE,OPTION);
  if(szTemp != NULL){
    ptParams->bOutputSample=TRUE;
  }
  else{
    ptParams->bOutputSample=FALSE;
  }

  szTemp = extractParameter(argc,argv,SAMPLE,OPTION);
  if(szTemp != NULL){
    ptParams->bSample=TRUE;
  }
  else{
    ptParams->bSample=FALSE;
  }

  szTemp = extractParameter(argc,argv,SEED,OPTION);
  if(szTemp != NULL){
    ptParams->nL = strtol(szTemp,&cError,10);
    if(*cError != '\0'){
      goto error;
    }
  }
  else{
    ptParams->nL = 1;
  }

  szTemp = extractParameter(argc,argv,EXTRAPOLATE,OPTION);
  if(szTemp != NULL){
    ptParams->nExtrapolate = strtol(szTemp,&cError,10);
    if(*cError != '\0'){
      goto error;
    }
  }
  else{
    ptParams->nExtrapolate = -1;
  }


  szTemp = extractParameter(argc,argv,N_ITERATIONS,OPTION);
  if(szTemp != NULL){
    ptParams->nMaxIter = strtol(szTemp,&cError,10);
    if(*cError != '\0'){
      goto error;
    }
  }
  else{
    ptParams->nMaxIter = DEF_MAX_ITER;
  }

  return;

 error:
  writeUsage(stdout);
  exit(EXIT_FAILURE);
}

void writeAbundanceData(const char *szFile, t_Data *ptData)
{
  FILE* ofp = fopen(szFile,"w");
  int i = 0, j = 0;
  int nS = ptData->nS, nN = ptData->nN;

  if(ofp){
    fprintf(ofp,"OTUs,");
    
    for(j = 0; j < nN - 1; j++){
      fprintf(ofp,"%s,",ptData->aszSampleNames[j]);
    }
    fprintf(ofp,"%s\n",ptData->aszSampleNames[j]);

    for(i = 0; i < nS; i++){
      fprintf(ofp,"%s,",ptData->aszOTUNames[i]);

      for(j = 0; j < nN - 1; j++){
	fprintf(ofp,"%d,",ptData->aanX[j][i]);
      }
      fprintf(ofp,"%d\n",ptData->aanX[nN - 1][i]);
    }

    fclose(ofp);
  }
  else{
    fprintf(stderr,"Failed writing abundance data to file %s\n",szFile);
    fflush(stderr);
  }
}

void readAbundanceData(const char *szFile, t_Data *ptData)
{
  int  **aanX = NULL;
  int  i = 0, j = 0, nS = 0, nN = 0;
  char szLine[MAX_LINE_LENGTH];
  FILE* ifp = NULL;

  ifp = fopen(szFile, "r");

  if(ifp){
    char* szTok   = NULL;
    char* pcError = NULL;

    fgets(szLine, MAX_LINE_LENGTH, ifp);
    szTok = strtok(szLine, DELIM);

    while(strtok(NULL, DELIM) != NULL){
      
      nN++;
    }

    while(fgets(szLine, MAX_LINE_LENGTH, ifp) != NULL){
    	nS++;
    }
    fclose(ifp);

    ifp = fopen(szFile, "r");	
    fgets(szLine, MAX_LINE_LENGTH, ifp);

    ptData->aszSampleNames = (char **) malloc(nN*sizeof(char*));

    szTok = strtok(szLine, DELIM);

    for(i = 0; i < nN; i++){
      szTok = strtok(NULL, DELIM);
      ptData->aszSampleNames[i] = strdup(szTok);
    }
	
    aanX = (int **) malloc(nN*sizeof(int*));
    for(i = 0; i < nN; i++){
      aanX[i] = (int *) malloc(nS*sizeof(int));
    }
    ptData->aszOTUNames = (char **) malloc(nS*sizeof(char*));
    for(i = 0; i < nS; i++){
    
      fgets(szLine, MAX_LINE_LENGTH, ifp);
      szTok = strtok(szLine, DELIM);
      ptData->aszOTUNames[i] = strdup(szTok);
      for(j = 0; j < nN; j++){
	szTok = strtok(NULL, DELIM);

	aanX[j][i] = strtol(szTok,&pcError,10);

	if(*pcError != '\0'){
	  goto formatError;
	}
      }
    }
  }
  else{
    fprintf(stderr, "Failed to open abundance data file %s aborting\n", szFile);
    fflush(stderr);
    exit(EXIT_FAILURE);
  }

  ptData->nS = nS;
  ptData->nN = nN;
  ptData->aanX = aanX;
  return;

 formatError:
  fprintf(stderr, "Incorrectly formatted abundance data file\n");
  fflush(stderr);
  exit(EXIT_FAILURE);
}

int compare_doubles(const void* a, const void* b) 
{
  double* arg1 = (double *) a;
  double* arg2 = (double *) b;
  if( *arg1 < *arg2 ) return -1;
  else if( *arg1 == *arg2 ) return 0;
  else return 1;
}       

void Stirling(double ***paadStirlingMatrix,double** padNormMatrix,unsigned long n){
  double **aadStirlingMatrix = NULL, *adNormMatrix = NULL;
  unsigned long nN2 = n + 2, nN1 = n + 1;
  double *adS=(double *) malloc(nN2*sizeof(double));
  double *adL=(double *) malloc(nN2*sizeof(double));
  double dMaxVal = 0.0;
  unsigned long i = 0, k = 0, j = 0;

  aadStirlingMatrix = (double **) malloc(nN1*sizeof(double*));
  if(!aadStirlingMatrix)
    goto memoryError;

  for(i = 0; i < nN2; i++){
    adS[i] = 0.0;
    adL[i] = 0.0;
  }

  for(i = 0; i < nN1; i++){
    aadStirlingMatrix[i] = (double *) malloc((i+1)*sizeof(double));
    
    if(!aadStirlingMatrix[i])
      goto memoryError;
    
    for(j = 0; j < i + 1; j++){
        aadStirlingMatrix[i][j] = 0.0;
    }

  }

  adNormMatrix = (double *) malloc(nN1*sizeof(double));
  if(!adNormMatrix)
    goto memoryError;

  if(!adS)
    goto memoryError;
  
  if(!adL)
    goto memoryError;

  adS[1] = 1.0;
  adL[1] = 0.0;

  aadStirlingMatrix[0][0] = 1.0;
  adNormMatrix[0] = 0.0;

  for(k = 2; k < nN1; k++){
    adL[k]  = 0;
    dMaxVal = 0.0;
        	
    for(j = k; j > 0;--j)
      adS[j]=adS[j-1]+(k-1)*adS[j];
        	
    for(j = 1;j <= k;j++)
      if(adS[j] > dMaxVal)
	dMaxVal=adS[j];
       	 	
    for(j = 1; j <= k;j++)
      adS[j]/=dMaxVal;

    adL[k] = adL[k-1]+ log(dMaxVal);
        
    for(j=1;j<=k;j++){
      aadStirlingMatrix[k-1][j-1]=adS[j];
      adNormMatrix[k-1]=adL[k];
    }

  }

  (*paadStirlingMatrix) = aadStirlingMatrix;
  (*padNormMatrix) = adNormMatrix;

  free(adS);
  free(adL);

  return;

 memoryError:
  fprintf(stderr, "Failed allocating memory in main\n");
  fflush(stderr);
  exit(EXIT_FAILURE);
}

void generateStirlingMatrix(unsigned long **aalStirling, unsigned long lK)
{
  int i = 0, j = 0;

  aalStirling[0][0] = 1;
  aalStirling[0][1] = 0;
  aalStirling[1][0] = 0;
  aalStirling[1][1] = 1;

  for(i = 2; i <= lK; i++){

    aalStirling[i][0] = (i - 1)*aalStirling[i-1][0];

    for(j = 1; j <= i; j++){
      
      aalStirling[i][j] = (i - 1)*aalStirling[i-1][j] + aalStirling[i-1][j - 1];
    }
  }
}

void generateStirlingDMatrix(double **aadStirling, unsigned long lK)
{
  int i = 0, j = 0;

  aadStirling[0][0] = 1.0;
  aadStirling[0][1] = 0;
  aadStirling[1][0] = 0;
  aadStirling[1][1] = 1.0;

  for(i = 2; i <= lK; i++){

    aadStirling[i][0] = (i - 1)*aadStirling[i-1][0];

    for(j = 1; j <= i; j++){
      
      aadStirling[i][j] = (i - 1)*aadStirling[i-1][j] + aadStirling[i-1][j - 1];
    }
  }
}

int maxX(t_Data *ptData)
{
	int i = 0, j = 0, nMax = 0, nS = ptData->nS, nN = ptData->nN;

	for(i = 0; i < nN; i++){
		for(j = 0; j < nS; j++){
			if(ptData->aanX[i][j] > nMax){
				nMax = ptData->aanX[i][j];
			}
		}
	}

	return nMax;
}

int minJ(t_Data *ptData)
{
  int i = 0, j = 0, nMin = 1e9, nS = ptData->nS, nN = ptData->nN;

  for(i = 0; i < nN; i++){
    int nTotal = 0;

    for(j = 0; j < nS; j++){
      nTotal += ptData->aanX[i][j];
    }

    if(nTotal < nMin){
      nMin = nTotal;
    }
  }
  return nMin;
}

void sumColumns(int *anC, int **aanX, int nN, int nS)
{
	int i = 0, j = 0;

	for(j = 0; j < nS; j++){
		anC[j] = 0;
		for(i = 0; i < nN; i++){
			anC[j] += aanX[i][j];
		}
	}
}

void sumRows(int *anR, int **aanX, int nN, int nS)
{
	int i = 0, j = 0;

	for(i = 0; i < nN; i++){
		anR[i] = 0;
		for(j = 0; j < nS; j++){
			anR[i] += aanX[i][j];
		}
	}
}

int Sum(int *anT, int nN)
{
  int i = 0;
  int nRet = 0;

  for(i = 0; i < nN; i++){
    nRet += anT[i];
  }

  return nRet;
}

double safeexp(double x)
{
  double result = 0.0;

  if (x < log(1e-323))
    result = log(1e-323);
  else if(x > log(1e308))
    result = log(1e308);
  else
    result = x ;

  return exp(result);
}

int selectCat(gsl_rng* ptGSLRNG, int nN, double* adP)
{
  double *adCP = (double *) malloc(nN*sizeof(double));
  int l = 0;
  double dRand = gsl_rng_uniform(ptGSLRNG);

  if(!adCP)
    goto memoryError;

  for(l = 0; l < nN; l++){
    if(l > 0){
      adCP[l] = adCP[l - 1] + adP[l];
    }
    else{
      adCP[l] = adP[l];
    }
  }

  l = 0;
  while(dRand > adCP[l]){
    l++;
  }

  free(adCP);

  return l;

 memoryError:
  fprintf(stderr,"Failed allocating memory in selectCat\n");
  fflush(stderr);
  exit(EXIT_FAILURE);
}

int selectIntCat(gsl_rng* ptGSLRNG, int nN, int* anN)
{
  double *adP = (double *) malloc(nN*sizeof(double));
  double *adCP = (double *) malloc(nN*sizeof(double));
  int nT = 0, l = 0;
  double dRand = gsl_rng_uniform(ptGSLRNG);

  if(!adP)
    goto memoryError;

  if(!adCP)
    goto memoryError;

  for(l = 0; l < nN; l++){
    nT += anN[l];
  }

  for(l = 0; l < nN; l++){
    adP[l] = ((double) anN[l])/((double) nT);
    if(l > 0){
      adCP[l] = adCP[l - 1] + adP[l];
    }
    else{
      adCP[l] = adP[l];
    }
  }

  l = 0;
  while(dRand > adCP[l]){
    l++;
  }

  free(adP);
  free(adCP);

  return l;

 memoryError:
  fprintf(stderr,"Failed allocating memory in selectCat\n");
  fflush(stderr);
  exit(EXIT_FAILURE);
}

void addUnobserved(t_Data* ptDataR, t_Data *ptData)
{
  int    i = 0, j = 0;
  int    nN = ptData->nN, nS = ptData->nS;
  int    **aanX = (int **) malloc(nN*sizeof(int*));
  int    nSDash = nS + 1;

  if(!aanX)
    goto memoryError;

  for(i = 0; i < nN; i++){
    aanX[i] = (int *) malloc(nSDash*sizeof(int));

    if(!aanX[i])
      goto memoryError;

    for(j = 0; j < nS; j++){
      aanX[i][j] = ptData->aanX[i][j];
    }
    aanX[i][nS] = 0;
  }

  ptDataR->aanX = aanX;
  ptDataR->nS = nSDash;
  ptDataR->nN = nN;

  return;

 memoryError:
  fprintf(stderr,"Failed allocating memory in generateData\n");
  fflush(stderr);
  exit(EXIT_FAILURE);
}

void copyData(t_Data* ptDataR, t_Data *ptData)
{
  int    i = 0, j = 0;
  int    nN = ptData->nN, nS = ptData->nS;
  int    **aanX = (int **) malloc(nN*sizeof(int*));

  if(!aanX)
    goto memoryError;

  ptDataR->nSize = ptData->nS;

  for(i = 0; i < nN; i++){
    aanX[i] = (int *) malloc(ptDataR->nSize*sizeof(int));

    if(!aanX[i])
      goto memoryError;

    for(j = 0; j < ptDataR->nSize; j++){
      aanX[i][j] = 0;
    }

    for(j = 0; j < nS; j++){
      aanX[i][j] = ptData->aanX[i][j];
    }
  }

  ptDataR->aanX = aanX;
  ptDataR->nS = nS;
  ptDataR->nN = nN;

  return;

 memoryError:
  fprintf(stderr,"Failed allocating memory in generateData\n");
  fflush(stderr);
  exit(EXIT_FAILURE);
}

void generateData(gsl_rng* ptGSLRNG, t_Data *ptData, int nN, double dTheta, int *anJ, double* adI, int nS, double *adM)
{
  int    i = 0, j = 0, n = 0, nSize = nS*10;
  int    **aanX = (int **) malloc(nN*sizeof(int*));
  int    nSDash = nS + 1;

  if(!aanX)
    goto memoryError;

  for(i = 0; i < nN; i++){
    aanX[i] = (int *) malloc(nSDash*sizeof(int));

    if(!aanX[i])
      goto memoryError;

    for(j = 0; j < nSDash; j++){
      aanX[i][j] = 0;
    }

    for(n = 0; n < anJ[i]; n++){
      double dImmigrate = adI[i]/(adI[i] + (double) n);
      double dRand = gsl_rng_uniform(ptGSLRNG);
      int    l = -1;

      if(dRand < dImmigrate){
	l = selectCat(ptGSLRNG, nSDash, adM);
      }
      else{
	double adP[nSDash];

	for(j = 0; j < nSDash; j++){
	  adP[j] = ((double) aanX[i][j])/((double) n);
	}

	l = selectCat(ptGSLRNG, nSDash, adP);
      }

      aanX[i][l]++;
    }
  }

  ptData->aanX = aanX;
  ptData->nS = nSDash;
  ptData->nN = nN;

  return;

 memoryError:
  fprintf(stderr,"Failed allocating memory in generateData\n");
  fflush(stderr);
  exit(EXIT_FAILURE);
}

void reallocateData(t_Data *ptData)
{
  int i = 0, j = 0;
  int nN = ptData->nN, nS = ptData->nS, nSize = ptData->nSize;

  for(i = 0; i < nN; i++){
    ptData->aanX[i] = (int *) realloc(ptData->aanX[i], nSize*sizeof(int));
    if(!ptData->aanX[i])
      goto memoryError;
    for(j = nS; j < nSize; j++){
      ptData->aanX[i][j] = 0;
    }
  }
  return;
 memoryError:
  fprintf(stderr,"Failed allocating memory in reallocateData\n");
  fflush(stderr);
  exit(EXIT_FAILURE);
}

void extrapolateDataHDP(gsl_rng* ptGSLRNG, t_Data *ptData, int nN, int nE, double dTheta, double* adI, int* anTSample)
{
  int    i = 0, j = 0, k = 0, n = 0;
  int    *anT = NULL, nT = 0;
  int    *anJ = (int *) malloc(nN*sizeof(int));
  int    bAdd = TRUE;
  if(!anJ)
    goto memoryError;

  anT = (int *) malloc(ptData->nSize*sizeof(int));
  if(!anT)
    goto memoryError;

  for(i = 0; i < ptData->nSize; i++){
    anT[i] = 0;
  }

  for(i = 0; i < ptData->nS; i++){
    anT[i] = anTSample[i];
    nT += anT[i];
  }

  for(i = 0; i < nN; i++){
    anJ[i] = 0;
    for(j = 0; j < ptData->nS; j++){
      anJ[i] += ptData->aanX[i][j];
    }
  }

  while(bAdd == TRUE){
    bAdd = FALSE;
    for(i = 0; i < nN; i++){
      if(anJ[i] < nE){
	int n = anJ[i];
	double dImmigrate = adI[i]/(adI[i] + (double) n);
	double dRand = gsl_rng_uniform(ptGSLRNG);
	int    l = -1;

	bAdd = TRUE;

	if(dRand < dImmigrate){
	  double dSpeciate = dTheta/(dTheta + (double) nT);
	
	  dRand = gsl_rng_uniform(ptGSLRNG);
	  if(dRand < dSpeciate){
	  
	    if(ptData->nS >= ptData->nSize){
	      ptData->nSize = ptData->nSize*2;
	    
	      reallocateData(ptData);
	    
	      anT = (int *) realloc(anT, ptData->nSize*sizeof(int));
	      if(!anT)
		goto memoryError;

	      for(j = ptData->nS; j < ptData->nSize; j++){
		anT[j] = 0;
	      }
	    }
	  
	    anT[ptData->nS] = 1;

	    ptData->aanX[i][ptData->nS] = 1;

	    ptData->nS = ptData->nS + 1;
	  }
	  else{
	    l = selectIntCat(ptGSLRNG, ptData->nS, anT);
	    anT[l]++;
	    ptData->aanX[i][l]++;
	  }
	  nT++;
	}
	else{
	  l = selectIntCat(ptGSLRNG, ptData->nS, ptData->aanX[i]);
	  ptData->aanX[i][l]++;
	}
	anJ[i]++;
      }
    }
  }
  
  free(anT);
  free(anJ); 
  return;

 memoryError:
  fprintf(stderr,"Failed allocating memory in generateData\n");
  fflush(stderr);
  exit(EXIT_FAILURE);
}


void generateDataHDP(gsl_rng* ptGSLRNG, t_Data *ptData, int nN, double dTheta, int *anJ, double* adI,int** panT)
{
  int    i = 0, j = 0, k = 0, n = 0;
  int    *anT = NULL, nT = 0;

  ptData->aanX = (int **) malloc(nN*sizeof(int*));
  if(!ptData->aanX)
    goto memoryError;

  ptData->nSize = INIT_SIZE;
  ptData->nS = 0;
  ptData->nN = nN;

  anT = (int *) malloc(ptData->nSize*sizeof(int));
  if(!anT)
    goto memoryError;

  for(i = 0; i < ptData->nSize; i++){
    anT[i] = 0;
  }

  for(i = 0; i < nN; i++){
    ptData->aanX[i] = (int *) malloc(ptData->nSize*sizeof(int));

    if(!ptData->aanX[i])
      goto memoryError;

    for(j = 0; j < ptData->nSize; j++){
      ptData->aanX[i][j] = 0;
    }
  }

  for(i = 0; i < nN; i++){
    for(n = 0; n < anJ[i]; n++){
      double dImmigrate = adI[i]/(adI[i] + (double) n);
      double dRand = gsl_rng_uniform(ptGSLRNG);
      int    l = -1;

      if(dRand < dImmigrate){
	double dSpeciate = dTheta/(dTheta + (double) nT);
	
	dRand = gsl_rng_uniform(ptGSLRNG);
	if(dRand < dSpeciate){
	  
	  if(ptData->nS + 1 == ptData->nSize){
	    ptData->nSize = ptData->nSize*2;
	    
	    reallocateData(ptData);
	    
	    anT = (int *) realloc(anT, ptData->nSize*sizeof(int));
	    if(!anT)
	      goto memoryError;

	    for(j = ptData->nS; j < ptData->nSize; j++){
	      anT[j] = 0;
	    }
	  }
	  
	  anT[ptData->nS] = 1;

	  ptData->aanX[i][ptData->nS] = 1;

	  ptData->nS = ptData->nS + 1;
	}
	else{
	  l = selectIntCat(ptGSLRNG, ptData->nS, anT);
	  anT[l]++;
	  ptData->aanX[i][l]++;
	}
	nT++;
      }
      else{
	l = selectIntCat(ptGSLRNG, ptData->nS, ptData->aanX[i]);
	ptData->aanX[i][l]++;
      }
    }
  }
  *panT = anT;
  return;

 memoryError:
  fprintf(stderr,"Failed allocating memory in generateData\n");
  fflush(stderr);
  exit(EXIT_FAILURE);
}

void generateDataStick(gsl_rng* ptGSLRNG, t_Data *ptData, int nN, double dTheta, int *anJ, double* adI, int nS, double *adM, int *pnSDash, double **padMDash)
{
  int    i = 0, j = 0, k = 0, n = 0, nSize = nS*10;
  int    **aanX = (int **) malloc(nN*sizeof(int*));
  int    nSDash = nS;
  double *adMDash = (double *) malloc(nSize*sizeof(double));
  
  if(!adMDash)
    goto memoryError;

  for(i = 0; i < nS + 1; i++){
    adMDash[i] = adM[i];
  }

  if(!aanX)
    goto memoryError;

  for(i = 0; i < nN; i++){
    aanX[i] = (int *) malloc(nSize*sizeof(int));

    if(!aanX[i])
      goto memoryError;

    for(j = 0; j < nSize; j++){
      aanX[i][j] = 0;
    }
  }

  for(i = 0; i < nN; i++){
    for(n = 0; n < anJ[i]; n++){
      double dImmigrate = adI[i]/(adI[i] + (double) n);
      double dRand = gsl_rng_uniform(ptGSLRNG);
      int    l = -1;

      if(dRand < dImmigrate){
	l = selectCat(ptGSLRNG, nSDash + 1, adMDash);
	
	if(l == nSDash){
	  double dBeta = gsl_ran_beta(ptGSLRNG, 1.0, dTheta);
	  double dOld = adMDash[nSDash];

	  if(nSDash == nSize){
	    int nOldSize = nSize;

	    nSize = nSize*2;

	    adMDash = (double *) realloc(adM,nSize*sizeof(double));
	    
	    if(!adMDash)
	      goto memoryError;

	    for(j = 0; j < nN; j++){
	      aanX[j] = (int *) realloc(aanX[i],nSize*sizeof(int));

	      if(!aanX[i])
		goto memoryError;

	      for(k = nOldSize; k < nSize; k++){
		aanX[j][k] = 0;
	      }
	    }
	  }

	  adMDash[nSDash] = dBeta*dOld;
	  adMDash[nSDash + 1] = (1.0 - dBeta)*dOld;

	  aanX[i][l]++;
	  nSDash = nSDash + 1;
	}
	else{
	  aanX[i][l]++;
	}
      }
      else{
	double adP[nSDash];

	for(j = 0; j < nSDash; j++){
	  adP[j] = ((double) aanX[i][j])/((double) n);
	}

	l = selectCat(ptGSLRNG, nSDash, adP);

	aanX[i][l]++;
      }
    }
  }

  ptData->aanX = aanX;
  ptData->nS = nSDash;
  ptData->nN = nN;

  (*padMDash) = adMDash;
  (*pnSDash) = nSDash;

  return;

 memoryError:
  fprintf(stderr,"Failed allocating memory in generateData\n");
  fflush(stderr);
  exit(EXIT_FAILURE);
}

void generateDataFixed(gsl_rng* ptGSLRNG, t_Data *ptData, int nN, double dTheta, int *anJ, double* adI, int nS, double *adM)
{
  int    i = 0, j = 0, k = 0, n = 0, nSize = nS;
  int    **aanX = (int **) malloc(nN*sizeof(int*));
  double dMSum = 0.0, adMDash[nS];

  if(!aanX)
    goto memoryError;

  for(i = 0; i < nN; i++){
    aanX[i] = (int *) malloc(nSize*sizeof(int));

    if(!aanX[i])
      goto memoryError;

    for(j = 0; j < nSize; j++){
      aanX[i][j] = 0;
    }
  }

  for(i = 0; i < nS; i++){
    dMSum += adM[i];
  }

  for(i = 0; i < nS; i++){
    adMDash[i] = adM[i]/dMSum;
  }

  for(i = 0; i < nN; i++){
    for(n = 0; n < anJ[i]; n++){
      double dImmigrate = adI[i]/(adI[i] + (double) n);
      double dRand = gsl_rng_uniform(ptGSLRNG);
      int    l = -1;

      if(dRand < dImmigrate){
	l = selectCat(ptGSLRNG, nS, adMDash);
	  
	aanX[i][l]++;
	
      }
      else{
	double adP[nS];

	for(j = 0; j < nS; j++){
	  adP[j] = ((double) aanX[i][j])/((double) n);
	}

	l = selectCat(ptGSLRNG, nS, adP);

	aanX[i][l]++;
      }
    }
  }

  ptData->aanX = aanX;
  ptData->nS = nS;
  ptData->nN = nN;

  return;

 memoryError:
  fprintf(stderr,"Failed allocating memory in generateData\n");
  fflush(stderr);
  exit(EXIT_FAILURE);
}

/* calculates negative log likelihood of neutral community assembly from this metapopulation*/
double calcNLLikelihood(t_Data *ptData, double dTheta, double* adI, double *adM)
{
  int i = 0, j = 0;
  int nN = ptData->nN;
  int nS = ptData->nS;
  int **aanX = ptData->aanX;
  double dRet = 0.0;

  for(i = 0; i < nN; i++){
    double adAlpha[nS];
    double dSumAlpha = 0.0;
    double dJ = 0.0;
    
    for(j = 0; j < nS; j++){
      adAlpha[j] = adM[j]*adI[i];
      dSumAlpha += adAlpha[j];
      dJ += (double) aanX[i][j];
    }

    for(j = 0; j < nS; j++){
      
      dRet -= gsl_sf_lngamma((double) aanX[i][j] + 1);

      dRet += gsl_sf_lngamma(adAlpha[j] + (double) aanX[i][j]) - gsl_sf_lngamma(adAlpha[j]); 
    
      // printf("%d %d %f\n",i,j,dRet);
    }

    dRet +=  gsl_sf_lngamma(dSumAlpha) - gsl_sf_lngamma(dSumAlpha + dJ);

    dRet += gsl_sf_lngamma(dJ + 1);
  }

  return dRet;
}

int calcR(t_Data *ptData)
{
  int i = 0, j = 0;
  int nN = ptData->nN;
  int nS = ptData->nS;
  int nR = 0;
  int **aanX = ptData->aanX;
  
  for(j = 0; j < nS; j++){
    double dF = 0.0;

    for(i = 0; i < nN; i++){
      dF += (double) aanX[i][j];
    }

    if(dF > 0.0){
      nR++;
    }
  }  

  return nR;
}

double calcMeanSpeciesEntropy(t_Data *ptData)
{
  int i = 0, j = 0;
  int nN = ptData->nN;
  int nS = ptData->nS;
  int nR = 0;
  int **aanX = ptData->aanX;
  double dRet = 0.0;
  double dMean = 0;
  for(j = 0; j < nS; j++){
    double dF = 0.0;
    double dH = 0.0;

    for(i = 0; i < nN; i++){
      dF += (double) aanX[i][j];
    }
    if(dF > 0.0){
      nR ++;
      for(i = 0; i < nN; i++){
	if(aanX[i][j] > 0){
	  double dP = ((double) aanX[i][j])/dF;

	  dH += -dP*log(dP);
	}
      }
      dMean += exp(dH);
    }
  }  

  return dMean / (double) nR;
}

double sampleTheta(gsl_rng *ptGSLRNG, double dTheta, int* anTR, int nN, int nS)
{
  double dTT = (double) Sum(anTR,nN);

  double dPhi=gsl_ran_beta(ptGSLRNG,dTheta + 1.0, dTT);
      
  double dRho=gsl_ran_bernoulli(ptGSLRNG, dTT/(dTT + dTheta));
  
  double dG1 = ALPHA + (double) nS - dRho;
      
  double dG2 = BETA - log(dPhi);

  return gsl_ran_gamma(ptGSLRNG, dG1, 1.0/dG2);
    
}

void sampleMetacommunity(int nIter, gsl_rng *ptGSLRNG, double *adM, double **aadMStore, double dTheta, int nS, int *anTC)
{
  double dMSum = 0.0;
  int i = 0;

  for(i = 0; i < nS; i++){
    adM[i] = gsl_ran_gamma (ptGSLRNG, (double) anTC[i], 1.0);
    dMSum += adM[i];
  }

  adM[nS] = gsl_ran_gamma (ptGSLRNG, dTheta, 1.0);
  
  dMSum += adM[nS];

  for(i = 0; i < nS + 1; i++)
    aadMStore[nIter][i] = adM[i] / dMSum;
}

void sampleMetacommunityHDP(gsl_rng *ptGSLRNG, double *adM, double dTheta, int nS, int *anT)
{
  double dMSum = 0.0;
  int i = 0;

  for(i = 0; i < nS; i++){
    adM[i] = gsl_ran_gamma (ptGSLRNG, (double) anT[i], 1.0);
    dMSum += adM[i];
  }

  adM[nS] = gsl_ran_gamma (ptGSLRNG, dTheta, 1.0);
  
  dMSum += adM[nS];

  for(i = 0; i < nS + 1; i++)
    adM[i] = adM[i] / dMSum;
}

void sampleImmigrationRates(int nIter, gsl_rng *ptGSLRNG, double **aadIStore, int nN, double *adW, double *adS, double* adJ, int* anTR)
{
  int i = 0;

  /*calculate immigration rates*/
  if(nIter > 0){
    for(i = 0; i < nN; i++){
      double dT1 = 0.0, dT2 = 0.0;
      adW[i] = gsl_ran_beta (ptGSLRNG, aadIStore[nIter - 1][i] + 1.0, adJ[i]);

      adS[i] = gsl_ran_bernoulli (ptGSLRNG, adJ[i]/(adJ[i] +  aadIStore[nIter - 1][i]));
      dT1 = ETA + (double) anTR[i] - adS[i];
      dT2 =  1.0/(NU - log(adW[i]));
      aadIStore[nIter][i] = gsl_ran_gamma (ptGSLRNG, dT1, dT2);		
    }
  }
  else{
    for(i = 0; i < nN; i++){
      adW[i] = gsl_ran_beta (ptGSLRNG, I_INIT + 1.0, adJ[i]);
	
      adS[i] = gsl_ran_bernoulli (ptGSLRNG, adJ[i]/(adJ[i] +  I_INIT));

      aadIStore[nIter][i] = gsl_ran_gamma (ptGSLRNG, ETA + (double) anTR[i] - adS[i], 1.0/(NU - log(adW[i])));
    }
  }
}

void sampleT(int nIter, gsl_rng* ptGSLRNG, int** aanX, int **aanT, double** aadIStore, double** aadMStore, 
	     double** aadStirlingMatrix, int nN, int nS, double *adLogProbV, double *adProbV, double *adCProbV)
{
  int i = 0, ii = 0, jj = 0;

  for(ii = 0; ii < nN; ii++){

    for(jj = 0; jj < nS; jj++){
      int nXX = aanX[ii][jj];

      if(nXX == 0){
	aanT[ii][jj] = 0;
      }
      else if(nXX == 1){
	aanT[ii][jj] = 1;
      }
      else{
	double D = MIN_DBL, dSum = 0.0, dLogNormC = 0.0, u4 = 0.0;
				
	for(i = 0; i < nXX;i++){
	  double i1 = i + 1.0;
	  adLogProbV[i] = log(aadStirlingMatrix[nXX - 1][i]) + (i1*log(aadIStore[nIter][ii]*aadMStore[nIter][jj]));
	}
	
	dSum = 0.0;
	for(i = 0; i < nXX; i++){  
	  adProbV[i] = safeexp(adLogProbV[i]);
	  dSum += adProbV[i];
	}
			
	adProbV[0] = adProbV[0]/dSum;
	adCProbV[0] = adProbV[0];

	for(i = 1; i < nXX; i++){
	  adProbV[i] /= dSum;
	  adCProbV[i] = adProbV[i] + adCProbV[i - 1];
	}
	
	u4=gsl_rng_uniform(ptGSLRNG);

	i = 0;
	while(u4 > adCProbV[i]){
	  i++;
	}

	aanT[ii][jj] = i + 1;
				
      }
    }
  }
}

void outputSamples(char *szOutputDir, int nIter, int nMaxIter, gsl_rng* ptGSLRNG, int nN, int nS, t_Params *ptParams, t_Data *ptData,
		   double *adThetaStore, int *anJ, double **aadIStore, double** aadMStore)
{
  t_Data tTempData, tTempDataH;
  FILE *sfp = NULL, *mfp = NULL;
  char szOutFile[MAX_LINE_LENGTH];
  char szSampleFile[MAX_LINE_LENGTH];
  char szMetaFile[MAX_LINE_LENGTH];
  t_Data tDataR;
  int** aanX = ptData->aanX;
  int i = 0, j = 0, k = 0;

  sprintf(szSampleFile,"%s%s",ptParams->szOutFileStub,SAMPLE_FILE_SUFFIX);
  sprintf(szMetaFile,"%s%s",ptParams->szOutFileStub,META_FILE_SUFFIX);

  addUnobserved(&tDataR, ptData);

  sfp = fopen(szSampleFile, "w");
  mfp = fopen(szMetaFile,"w");

  if(sfp && mfp){
    int aanF[nS][nN + 1];
    int nTotal = 0;
    double adMMean[nS];

    for(i = 0; i < nS; i++){
      for(j = 0; j <= nN; j++){
	aanF[i][j] = 0;
      }
      adMMean[i] = 0.0;
    }

    for(i = ptParams->nBurnIter; i < nMaxIter; i++){
      if(i % N_SAMPLE == 0){
	int nSDash = 0, anF[nS];
	double dHS = 0.0, dHH = 0.0, dHR = 0.0, dHF = 0.0, dNLLS = 0.0, dNLLR = 0.0, dNLLH = 0.0, dNLLF = 0.0;
	double *adMDash = NULL;
	double *adMH = NULL;
	int    *anT = NULL;
	char   szStickSampleFile[MAX_LINE_LENGTH];
	
	generateDataStick(ptGSLRNG, &tTempData, nN, adThetaStore[i], anJ, aadIStore[i], nS, aadMStore[i], &nSDash, &adMDash);

	tTempData.aszOTUNames = (char **) malloc(tTempData.nS*sizeof(char*));
	for(j = 0; j < ptData->nS; j++){
	  tTempData.aszOTUNames[j] = strdup(ptData->aszOTUNames[j]);
	}
	for(j = ptData->nS; j < tTempData.nS;j++){
	  char*szTemp = malloc(MAX_LINE_LENGTH*sizeof(char));

	  sprintf(szTemp,"D%d",j);
	  tTempData.aszOTUNames[j] = szTemp; 
	}

	/*output sampled file*/
	if(ptParams->bOutputSample){
	  sprintf(szStickSampleFile,"%s/%s%d%s",szOutputDir,ptParams->szOutFileStub,i,CSV_FILE_SUFFIX);
	
	  tTempData.aszSampleNames = ptData->aszSampleNames;
	
	  writeAbundanceData(szStickSampleFile, &tTempData);
	}

	generateDataHDP(ptGSLRNG, &tTempDataH, nN, adThetaStore[j], anJ, aadIStore[i], &anT);

	adMH = (double *) malloc(sizeof(double)*(tTempDataH.nS + 1));

	if(bVerbose == TRUE){
	  printf("Sampling %d:\n",i);
	}
     
	sampleMetacommunityHDP(ptGSLRNG, adMH, adThetaStore[i], tTempDataH.nS, anT);

	dNLLS = calcNLLikelihood(&tTempData, adThetaStore[i], aadIStore[i], adMDash);
	dNLLR = calcNLLikelihood(&tDataR, adThetaStore[i], aadIStore[i], aadMStore[i]);
	dNLLH = calcNLLikelihood(&tTempDataH, adThetaStore[i], aadIStore[i],adMH);

	dHS = calcMeanSpeciesEntropy(&tTempData);
	dHR = calcMeanSpeciesEntropy(&tDataR);
	dHH = calcMeanSpeciesEntropy(&tTempDataH);

	fprintf(sfp,"%d,%f,%f,%f,%f,%f,%f,%d,%d,%d\n",i,dNLLH,dNLLS,dNLLR,dHH,dHS,dHR,calcR(&tTempDataH),calcR(&tTempData),calcR(&tDataR));

	fprintf(mfp,"%d,%d,%d\n",i,nS,tTempDataH.nS);
	for(j = 0; j < nS; j++){
	  fprintf(mfp,"%f,",aadMStore[i][j]);
	}
	fprintf(mfp,"%f\n",aadMStore[i][nS]);
	
	for(j = 0; j < tTempDataH.nS; j++){
	  fprintf(mfp,"%f,",adMH[j]);
	}
	fprintf(mfp,"%f\n",adMH[tTempDataH.nS]);

	for(j = 0; j < nS; j++){
	  adMMean[j] += aadMStore[i][j];

	  anF[j] = 0;
	  for(k = 0; k < nN; k++){
	    if(tTempData.aanX[k][j] > 0){
	      anF[j]++;
	    }
	  }

	  aanF[j][anF[j]]++;
	}
	nTotal++;
	for(j = 0; j < nN; j++){
	  free(tTempData.aanX[j]);
	  free(tTempDataH.aanX[j]);
	}
	for(j = 0; j < tTempData.nS; j++){
	  free(tTempData.aszOTUNames[j]);
	}
	free(tTempData.aszOTUNames);

	free(tTempData.aanX);
	free(tTempDataH.aanX);
	free(adMDash);
	free(adMH);
	free(anT);
     }
   }
   
   fclose(mfp);
   fclose(sfp);
 }
 else{
   fprintf(stderr, "Failed to open %s for writing\n",szOutFile);
   fflush(stderr);
 }

 for(i = 0; i < nN; i++){
   free(tDataR.aanX[i]);
 }
 free(tDataR.aanX);

 return;
}

void extrapolateSamples2(char *szOutputDir, int nIter, int nMaxIter, gsl_rng* ptGSLRNG, int nE, int nN, int nS, t_Params *ptParams, t_Data *ptData,
			double *adThetaStore, int *anJ, double **aadIStore, int** aanTStore)
{
 
 char szOutFile[MAX_LINE_LENGTH];
 int i = 0, j = 0, k = 0;
 char szExtraFile[MAX_LINE_LENGTH];
 FILE *ofp = NULL;
 sprintf(szOutFile,"%s%s",ptParams->szOutFileStub,EXTRAPOLATE_FILE_SUFFIX);
   
 for(i = ptParams->nBurnIter; i < nMaxIter; i++){
     if(i % N_SAMPLE == 0){
       int nSDash = 0, nO = 0, anF[nS];
       t_Data tDataR;
     
       copyData(&tDataR, ptData);

       extrapolateDataHDP(ptGSLRNG, &tDataR, nN, nE, adThetaStore[i], aadIStore[i], aanTStore[i]);
       
       ofp = fopen(szOutFile, "a");   
       if(ofp){
	        fprintf(ofp,"%d,",i);
	        for(j = 0; j < nN; j++){
	            int nSN = 0;
	            for(k = 0; k < tDataR.nS;k++){
	                if(tDataR.aanX[j][k] > 0){
	                    nSN++;
	                }
	            }
	            if(j < nN - 1){
	                fprintf(ofp,"%d,",nSN);
	            }
	            else{
	                fprintf(ofp,"%d\n",nSN);
	            }
	        }
	        fclose(ofp);
       }

       if(ptParams->bOutputSample){
	    tDataR.aszOTUNames = (char **) malloc(tDataR.nS*sizeof(char *));
	    for(k = 0; k < tDataR.nS; k++){
	        tDataR.aszOTUNames[k] = (char *) malloc(MAX_LINE_LENGTH*sizeof(char));
	         sprintf(tDataR.aszOTUNames[k],"D%d",k);
	    }

	    sprintf(szExtraFile,"%s/%s%de%s",szOutputDir,ptParams->szOutFileStub,i,CSV_FILE_SUFFIX);
		 
	    tDataR.aszSampleNames = ptData->aszSampleNames;
	
	    writeAbundanceData(szExtraFile, &tDataR);
	 
	    for(k = 0; k < tDataR.nS; k++){
	      free(tDataR.aszOTUNames[k]); 
	    }

	    free(tDataR.aszOTUNames);
       }

       for(j = 0; j < nN; j++){
	        free(tDataR.aanX[j]);
       }
       free(tDataR.aanX);
     }
   }

 return;
}

void extrapolateSamples(int nIter, int nMaxIter, gsl_rng* ptGSLRNG, int nN, int nS, t_Params *ptParams, t_Data *ptData,
			double *adThetaStore, int *anJ, double **aadIStore, double** aadMStore)
{
 t_Data tTempData;
 FILE *ofp = NULL;
 char szOutFile[MAX_LINE_LENGTH];
 t_Data tDataR;
 int** aanX = ptData->aanX;
 int i = 0, j = 0, k = 0;
 int anE[nN];

 for(i = 0; i < nN; i++){
   anE[i] = ptParams->nExtrapolate;
 }

 sprintf(szOutFile,"%s%s",ptParams->szOutFileStub,EXTRAPOLATE_FILE_SUFFIX);

 addUnobserved(&tDataR, ptData);
    
 ofp = fopen(szOutFile, "w");

 if(ofp){
   for(i = ptParams->nBurnIter; i < nMaxIter; i++){
     int nSDash = 0, nO = 0, anF[nS];
     double *adMDash = NULL;

     generateDataStick(ptGSLRNG, &tTempData, nN, adThetaStore[i], anE, aadIStore[i], nS, aadMStore[i], &nSDash, &adMDash);

     for(j = 0; j < nSDash; j++){
       for(k = 0; k < nN; k++){
	 if(tTempData.aanX[k][j] > 0){
	   nO++;
	   break;
	 }
       }
     }
     
     fprintf(ofp,"%d,%d,%d\n",i,nO,nSDash);
     for(j = 0; j < nN; j++){
       free(tTempData.aanX[j]);
     }
     free(tTempData.aanX);
     free(adMDash);
   }
   fclose(ofp);  
 }
 else{
   fprintf(stderr, "Failed to open %s for writing\n",szOutFile);
   fflush(stderr);
 }

 for(i = 0; i < nN; i++){
   free(tDataR.aanX[i]);
 }
 free(tDataR.aanX);

 return;
}

void rarefy(gsl_rng *ptGSLRNG,int nMaxJ, t_Data *ptDataR, t_Data *ptData)
{
  int n = 0, i = 0, j = 0, nS = ptData->nS, nN = ptData->nN;
  int **aanX = NULL;
  int nSDash = 0;
  int anC[nS];
  int nC = 0;

  aanX = (int **) malloc(nN*sizeof(int *));
  if(!aanX)
    goto memoryError;

  for(i = 0; i < nN; i++){
    aanX[i] = (int *) malloc(nS*sizeof(int));
    if(!aanX[i])
      goto memoryError;

    for(j = 0; j < nS; j++){
      aanX[i][j] = 0;
    }
  }

  for(i = 0; i < nN; i++){
    for(n = 0; n < nMaxJ; n++){
      int nI = selectIntCat(ptGSLRNG, nS, ptData->aanX[i]);
      aanX[i][nI]++;
    }
  }

  for(j = 0; j < nS; j++){
    int nJ = 0;
    for(i = 0; i < nN; i++){
      if(aanX[i][j] > 0){
	nJ = 1;
      }
    }
    if(nJ == 1){
      anC[j] = nC;

      nC++;
    }
    else{
      anC[j] = -1;
    }
    nSDash += nJ;
  }

  ptDataR->aanX = (int **) malloc(nN*sizeof(int *));
  if(!ptDataR->aanX)
    goto memoryError;

  ptDataR->nN = nN;
  ptDataR->nS = nSDash;
  ptDataR->nSize = nSDash;

  for(i = 0; i < nN; i++){
    ptDataR->aanX[i] = (int *) malloc(nSDash*sizeof(int));
    if(!ptDataR->aanX[i])
      goto memoryError;

    for(j = 0; j < nSDash; j++){
      ptDataR->aanX[i][j] = 0;
    }
  }

  for(i = 0; i < nN; i++){

    for(j = 0; j < nS; j++){
      if(anC[j] >= 0){
	ptDataR->aanX[i][anC[j]] = aanX[i][j];
      }
    }
  }
  
  ptDataR->aszSampleNames = (char **) malloc(nN*sizeof(char*));
  if(! ptDataR->aszSampleNames)
    goto memoryError;

  for(i = 0; i < nN; i++){
    ptDataR->aszSampleNames[i] = strdup(ptData->aszSampleNames[i]);
  }

  ptDataR->aszOTUNames = (char **) malloc(nSDash*sizeof(char*));
  if(! ptDataR->aszOTUNames)
    goto memoryError;

  for(i = 0; i < nSDash; i++){
    ptDataR->aszOTUNames[i] = (char *) malloc(MAX_LINE_LENGTH*sizeof(char));
    if(!ptDataR->aszOTUNames[i])
      goto memoryError;

    sprintf(ptDataR->aszOTUNames[i],"D%d",i);
  }

  /*free up memory*/
  for(i = 0; i < nN; i++){
    free(aanX[i]);
  }
  free(aanX);
  return;
 memoryError:
  fprintf(stderr,"Failed allocating memory in rarefy\n");
  fflush(stderr);
  exit(EXIT_FAILURE);
}
