#include <math.h>
#include <R.h>
#include <Rinternals.h>
#include <Rdefines.h>
#include <R_ext/Rdynload.h>
#include "UtilFunc.h"
/** Finding the peaks ina FIA file, using a first step of peak detection on mass traces
and then to rely line by lines all those peaks
**/

#define ALLOCATION_INCREMENT 2
#define MULTIPLIER_NUM_SCAN 1.5
#define DEBUGGING 0
#define MAX_BAND 8
#define NUM_EXCEPT 10


//V2 now only the peaks are supposed to be pushed in.
//Registering the routines.



struct scanStruct
{
    double  mz;
    double  intensity;
};

struct scanFull
{
    struct scanStruct * scan;
    unsigned int scanLength;
};


struct centroid
{
    int scan;
    double mz;
    double intensity;
};



struct band
{
    struct centroid * seqCentroid;
    int size;
    int maxSize;
    double meanMz;
    struct band * prevBand;
    struct band * nextBand;
};

struct bandList
{
    struct band * head;
    int size;
};

void debuggingBVizualisation(struct band * posBand)
{
    int i = 0;
    Rprintf("mz : %0.4f size : %d\nmz : ",posBand->meanMz,posBand->size);
    for(i=0; i<posBand->size; i++)
    {
        Rprintf("%0.2f",posBand->seqCentroid[i].mz);
    }
    Rprintf("\nScans : ");
    for(i=0; i<posBand->size; i++)
    {
        Rprintf("%d ",posBand->seqCentroid[i].scan);
    }
    Rprintf("\nIntensity : ");
    for(i=0; i<posBand->size; i++)
    {
        Rprintf("%0.2f ",posBand->seqCentroid[i].intensity);
    }
    Rprintf("\n\n");
}



void debuggingVizualisation(struct bandList * bL)
{
    struct band * posBand=bL->head;
    int countBand=1;
    while(posBand!=NULL)
    {
        Rprintf("Band : %d %0.5f scmin %d scmax %d size %d\n",countBand,posBand->meanMz,posBand->seqCentroid[0].scan,posBand->seqCentroid[posBand->size-1].scan,posBand->size);
        //debuggingBVizualisation(posBand);
        countBand++;
        posBand=posBand->nextBand;
    }

}


//PHD VISUS ONLY.
int size_bl(struct bandList *bL){
    struct band * posBand=bL->head;
    int countBand=0;
    while(posBand!=NULL)
    {
        countBand++;
        posBand=posBand->nextBand;
    }
    return countBand;
}


double ponderate_mean(double * weights,double * values,int left,int right)
{
    int i;
    double r=0;
    double accu_weights=0;
    for(i=left; i<right; i++)
    {
        r=r+weights[i]*values[i];
        accu_weights=accu_weights+weights[i];
    }
    return(r/accu_weights);
}

double trapzApprox(struct band * nBand, double * scantime)
{
    //debuggingBVizualisation(nBand);
    int size=nBand->size;
    double accumulator=0;
    int i=0;
    for(i=1; i<size; i++)
    {
        accumulator=accumulator+(nBand->seqCentroid[i-1].intensity+nBand->seqCentroid[i].intensity)*(scantime[nBand->seqCentroid[i].scan-1]-scantime[nBand->seqCentroid[i-1].scan-1]);
    }
    return(accumulator);
}


double distIntMz(struct centroid C1,struct centroid C2)
{
    double dist=fabs(C1.mz-C2.mz)*1000000/C1.mz+fabs(log10(C2.intensity)-log10(C1.intensity));
    return(dist);
}


struct band * fuseBand(struct band * B1,struct band * B2,int first_scan, int max_scan,int * fused,double * scantime,struct bandList * bL)
{
    int posB1=0;
    int posB2=0;
    int posN=0;
    int nexception=0;
    if(B1->size+B2->size>max_scan-first_scan+1)
    {
        *fused=0;
        return(B2);
    }
    struct band * nBand=malloc(sizeof(struct band));
    nBand->seqCentroid=malloc((B1->size+B2->size)*sizeof(struct centroid));
    nBand->size=B1->size+B2->size;
    *fused=1;
    while(posB1<B1->size&&posB2<B2->size)
    {
        if(B1->seqCentroid[posB1].scan==B2->seqCentroid[posB2].scan)
        {
            //Case fusion forced
            nexception++;
            if(nexception>=NUM_EXCEPT)
            {
                *fused=0;
                break;
            }
            if(posN==0)
            {
                if(B1->seqCentroid[posB1].intensity>B2->seqCentroid[posB2].intensity)
                {
                    nBand->seqCentroid[posN]=B1->seqCentroid[posB1];
                    posB1++;
                    posN++;
                }
                else
                {
                    //Rprintf("c");
                    nBand->seqCentroid[posN]=B2->seqCentroid[posB2];
                    posB2++;
                    posN++;
                }

            }
            else
            {
                //Rprintf("d");
                double D1=distIntMz(B1->seqCentroid[posB1],nBand->seqCentroid[posN-1]);
                double D2=distIntMz(B2->seqCentroid[posB2],nBand->seqCentroid[posN-1]);
                //Rprintf("dk");
                if(D2>D1)
                {
                    nBand->seqCentroid[posN]=B2->seqCentroid[posB2];
                    posB2++;
                    posN++;
                }
                else
                {

                    nBand->seqCentroid[posN]=B1->seqCentroid[posB1];
                    posB1++;
                    posN++;
                }
            }


        }
        else
        {
            if(B1->seqCentroid[posB1].scan>B2->seqCentroid[posB2].scan)
            {
                //Rprintf(" fC \n");
                nBand->seqCentroid[posN]=B2->seqCentroid[posB2];
                posB2++;
                posN++;
            }
            else
            {
                //Rprintf(" fD \n");
                nBand->seqCentroid[posN]=B1->seqCentroid[posB1];
                posB1++;
                posN++;
            }
        }
    }

    if(*fused)
    {
        while(posB1<B1->size)
        {
            //Rprintf(" B1 ");
            nBand->seqCentroid[posN]=B1->seqCentroid[posB1];
            posB1++;
            posN++;
        }
        //Rprintf(" EB1 ");
        while(posB2<B2->size)
        {
            //Rprintf(" B2 ");
            nBand->seqCentroid[posN]=B2->seqCentroid[posB2];
            posB2++;
            posN++;
        }
        //Rprintf(" fused ");
        //nBand->meanMz=(B1->meanMz*B1->intensity+B2->meanMz*B2->intensity)/(B1->intensity+B2->intensity);
        nBand->meanMz=(B1->meanMz*B1->size+B2->meanMz*B2->size)/(B1->size+B2->size);
        //nBand->intensity=trapzApprox(nBand,scantime);
        nBand->prevBand=B1->prevBand;
        if(nBand->prevBand!=NULL)
        {
            nBand->prevBand->nextBand=nBand;
        }
        else
        {
            nBand->prevBand=NULL;
            bL->head=nBand;
        }
        //Rprintf(" Prevband");
        nBand->nextBand=B2->nextBand;
        if(nBand->nextBand!=NULL)
        {
            nBand->nextBand->prevBand=nBand;
        }
        else
        {
            nBand->nextBand=NULL;
        }
        //Rprintf(" Pointers_Set ");
        //Rprintf("trpz OK !!\n");
        //Rprintf(" fB12 ");
        free(B1->seqCentroid);
        free(B2->seqCentroid);
        //Rprintf("free tPeaksMz OK !!\n");
        free(B1);
        free(B2);
        //Rprintf("free B1B2 and fusion OK !!\n");
    }
    else
    {
        //Rprintf(" fBN ");
        free(nBand->seqCentroid);
        free(nBand);
        //Rprintf("g\n\n");
        return(B2);
    }
    //Rprintf("h\n\n");
    return(nBand);
}


void extendBand(struct band * b,int n, int min_scan, int max_scan){
    if(n>=b->maxSize){
        int csize;
        if((max_scan-min_scan+1)<(b->maxSize*2)){
            csize = (max_scan-min_scan+1);
        }else{
            csize = (b->maxSize*2);
        }
        struct centroid * scentroid = calloc(csize,sizeof(struct centroid));

        int i=0;
        for(i =0;(i<b->size);i++){
            scentroid[i] = b->seqCentroid[i];
        }
        b->seqCentroid = scentroid;
        b->maxSize = csize;
    }
}


void insertBandAfter(struct bandList * bL, struct band ** i,struct centroid nCentroid, int first_scan, int max_scan)
{
    struct band *nBand=calloc(1,sizeof(struct band));
    nBand->seqCentroid=calloc(5,sizeof(struct centroid));
    nBand->seqCentroid[0]=nCentroid;
    nBand->size=1;
    nBand->maxSize=5;
    nBand->nextBand=(*i)->nextBand;
    nBand->prevBand=*i;
    nBand->meanMz=nCentroid.mz;
    if((*i)->nextBand!=NULL)
    {
        (*i)->nextBand->prevBand=nBand;
    }
    (*i)->nextBand=nBand;
    bL->size++;
}

void insertBandBefore(struct bandList * bL, struct band ** i,struct centroid nCentroid, int first_scan, int max_scan)
{
    struct band *nBand=calloc(1,sizeof(struct band));
    nBand->seqCentroid=calloc(5,sizeof(struct centroid));
    nBand->seqCentroid[0]=nCentroid;
    nBand->size=1;
    nBand->maxSize=5;
    nBand->nextBand=*i;
    nBand->prevBand=(*i)->prevBand;
    nBand->meanMz=nCentroid.mz;
    if((*i)->prevBand!=NULL)
    {
        (*i)->prevBand->nextBand=nBand;
    }
    else
    {
        bL->head=nBand;
    }
    (*i)->prevBand=nBand;
    (*i)=nBand;
    bL->size++;
}




void insertCentroidBandList(struct centroid nCentroid,struct band ** i, int max_scan, int first_scan, struct bandList * bL,double ppm, double rdmz)
{
    //Case where the peak is the first of the list
    if(bL->head==NULL)
    {
        //Rprintf("  FPeak\n");
        struct band * nBand=malloc(sizeof(struct band));
        nBand->nextBand=NULL;
        nBand->prevBand=NULL;
        nBand->meanMz=nCentroid.mz;
        nBand->seqCentroid=calloc(5,sizeof(struct centroid));
        nBand->size=1;
        //nBand->intensity=0;
        bL->head=nBand;
        bL->head->size=1;
        bL->head->maxSize=5;
        bL->head->seqCentroid[0]=nCentroid;
        bL->size=1;
        *i=bL->head;
        return;
    }
    int nfound=0;
    struct band * lfound[MAX_BAND];
    while(nfound<MAX_BAND)
    {
        double mass_diff = rdmz>ppm*0.000001*nCentroid.mz ? rdmz : (ppm*0.000001*nCentroid.mz);
        //Case where the peak is found.
        if(fabs((*i)->meanMz-nCentroid.mz)<=mass_diff)
        {
            lfound[nfound]=(*i);
            nfound=nfound+1;
        }
        //Case where it is the last of the bandList and not found
        if((*i)->nextBand==NULL)
        {
            break;
        }
        //Case where we passed throught the peak
        if((*i)->meanMz-nCentroid.mz>mass_diff)
        {
            break;
        }
        //If we still in the loop the peak is incremented
        *i=(*i)->nextBand;
    }
    //If the peak have not been found
    if(nfound==0)
    {
        if((*i)->nextBand==NULL) //End of the list, peak need to be added in last position.
        {
            insertBandAfter(bL,i,nCentroid,first_scan,max_scan);
        }
        else   //Case where we just went too far.
        {
            insertBandBefore(bL,i,nCentroid,first_scan,max_scan);
        }
        return; //Leaving the function.
    }

    //While the situation is not okay, bands can be modified.
    int j;
    struct centroid toinsert=nCentroid;
    int tab_trial[MAX_BAND];
    struct centroid memcentroid;
    for(j=0; j<nfound; j++)
    {
        tab_trial[j]=0;
    }
    int ntrial=0;
    int bp;
    while(ntrial!=nfound)
    {
        bp=0;

        //Initialised at the maximum limit.
        double minstate=ppm*toinsert.mz*0.000001;
        double norm=ppm*toinsert.mz*0.000001;
        double current_dmz;
        //Determining the best Pos which have not been tried to add the peak.
        for(j=0; j<nfound; j++)
        {
            if(tab_trial[j]==1) continue;
            //Distance is taken as the minimum between the centroid and the intensity
            current_dmz=fabs(lfound[j]->seqCentroid[lfound[j]->size-1].mz-toinsert.mz)/norm+fabs((log10(toinsert.intensity)-log10(lfound[j]->seqCentroid[lfound[j]->size-1].intensity)));
            if(current_dmz<=minstate)
            {
                minstate=current_dmz;
                bp=j;
            }
        }
        //Trying to add the point to his best position
        struct band * cband=lfound[bp];
        if((cband->size!=1)&(cband->seqCentroid[cband->size-1].scan==toinsert.scan)) //Case were there already is a point.
        {
            //Testing which point is the closest;
            double old_centroid=((cband->meanMz*cband->size)-cband->seqCentroid[cband->size-1].mz)/(cband->size-1);
            if(fabs(cband->seqCentroid[cband->size-1].mz-old_centroid)>fabs(toinsert.mz-old_centroid))
            {
                memcentroid = toinsert;
                cband->meanMz=(old_centroid*(cband->size-1)-cband->seqCentroid[cband->size-1].mz+memcentroid.mz)/(cband->size-1);
                toinsert = cband->seqCentroid[cband->size-1];
                cband->seqCentroid[cband->size-1] = memcentroid;

                //Reseting the trial tab
                for(j=0; j<nfound; j++)
                {
                    tab_trial[j]=0;
                }
                ntrial=0;

            }
            else
            {
                tab_trial[bp]=1;
                ntrial=ntrial+1;
            }

        }
        else   //Best case, the peak can be inserted in a line, in this case the function stop.
        {

            extendBand(cband,cband->size,first_scan,max_scan);
            cband->seqCentroid[cband->size]=toinsert;
            cband->meanMz=(cband->meanMz*(cband->size)+toinsert.mz)/(cband->size+1);
            cband->size=cband->size+1;
            (*i)=lfound[0];
            return;
        }
        //if(ntrial==nfound) Rprintf("   %d:%d Outcond   ",ntrial,nfound);
    }
    j=0;
    while(lfound[j]->meanMz<toinsert.mz)
    {
        j++;
        if(j==nfound) break;
    }

    if(j==nfound)
    {
        insertBandAfter(bL,&(lfound[j-1]),toinsert,first_scan,max_scan);
    }
    else
    {
        insertBandBefore(bL,&(lfound[j]),toinsert,first_scan,max_scan);
    }
//Placing the pointer on the first ok pos.
    (*i)=lfound[0];
}



//Return the scan given by the number I
void getScan(int scan_number, double * mz, double * intensity, int * scan_index, int last_scan, int num_mz, struct scanFull * this_scan)
{
    int b_1,b_2,i,N;
    if(scan_number>last_scan)
    {
        error("Invalid scan required !");
    }
    b_1=scan_index[scan_number-1];
    if(scan_number==last_scan)
    {
        b_2=num_mz-1;
    }
    else
    {
        b_2=scan_index[scan_number]-1;
    }
    N=b_2-b_1+1;
    if(this_scan->scan!=NULL)
    {
        free(this_scan->scan);
    }
    this_scan->scan=calloc(N,sizeof(struct scanStruct));
    this_scan->scanLength=N;
    for(i=b_1; i<=b_2; i++)
    {
        this_scan->scan[i-b_1].mz=mz[i];
        this_scan->scan[i-b_1].intensity=intensity[i];
    }
    this_scan->scanLength=N;
    //Rprintf("Scan number %d is %d numbers long b_1 %d b_2 %d !\n",scan_number,this_scan->scanLength,b_1,b_2);
}



void removeBand(struct bandList * bL,struct band * posBand)
{
    if(posBand->prevBand==NULL)
    {
        bL->head=posBand->nextBand;
        posBand->nextBand->prevBand=NULL;
    }
    else
    {
        if(posBand->nextBand==NULL)
        {
            posBand->prevBand->nextBand=NULL;
        }
        else
        {
            posBand->nextBand->prevBand=posBand->prevBand;
            posBand->prevBand->nextBand=posBand->nextBand;
        }
    }
    free(posBand->seqCentroid);
    free(posBand);
    posBand=NULL;
    bL->size--;
}




void fuseBandList(struct bandList * bL,double ppm,int first_scan, int max_scan, double * scantime, int *numFused, double dmz)
{
    *numFused=0;
    if(bL->size<=1)
    {
        Rprintf("No band detected impossible to process !\n");
        return;
    }
    struct band * posBand=bL->head->nextBand;
    struct band * previousBand=bL->head;
    int fused=0;
    int counter=0;
    //Rprintf("Beginning to fuse bands !\n");
    while(posBand!=NULL)
    {
        //Rprintf("\n%d ",counter);
        double mass_diff= dmz>previousBand->meanMz*ppm*1e-6 ? dmz : (previousBand->meanMz*ppm*1e-6);
        if(fabs(posBand->meanMz-previousBand->meanMz)<mass_diff)
        {
            //Rprintf("\ncounter %d A ",counter);
            //debuggingBVizualisation(posBand);
            //debuggingBVizualisation(previousBand);
            //Rprintf("fusing %0.5f and %0.5f",previousBand->meanMz,posBand->meanMz);
            posBand=fuseBand(previousBand,posBand,first_scan,max_scan,&fused,scantime,bL);
            //Rprintf("afterfusion %0.5f and %0.5f   ",previousBand->meanMz,posBand->meanMz);
            *numFused=*numFused+fused;
            //Rprintf(" C ");
            previousBand=posBand;
            //Rprintf(" D ");
            posBand=posBand->nextBand;
            //Rprintf(" Ee ");
            if(posBand==NULL||previousBand==NULL)
                break;
            //Rprintf("afterorder %0.5f and %0.5f\n",previousBand->meanMz,posBand->meanMz);
            continue;
        }
        previousBand=posBand;
        posBand=posBand->nextBand;
        counter++;
    }
    bL->size=bL->size-*numFused;
    //Rprintf(" OUT ");
}

/**void calculateIntensity(struct band * posBand,double * scantime){


}
**/

//Quick clean up to remove at least the lvl 1 of the peaks.
void cleanUpBandList(struct bandList * bL,int size_min, int injsc, int injend,double frac_min)
{
    int i=0;
    struct band * posBand =bL->head;
    struct band * to_remove;
    //Loop an all the detected Bands.
    /**TO DEBUG**/
    int counterBand=1;
    int injSize = injend-injsc+1;
    int limitInj = ceil(((double)injSize)*frac_min);
    //Rprintf("limitInj %d and size_min %d\n",limitInj,size_min);

    while(1)
    {
        //Rprintf("\nClean : %0.5f  %0.5f   %d  %d\n",posBand->tPeaksMz[0].low_limit,posBand->tPeaksMz[0].high_limit,posBand->tPeaksMz[0].scan,posBand->tPeaksMz[posBand->size-1].scan);
        if(posBand==NULL)
        {
            //Rprintf("c%db ",counterBand);
            break;
        }
        //Rprintf("\n The value posBand : %f %d      ",posBand->meanMz,posBand->size);
        if(posBand->size==1)
        {
            //Rprintf("c%dd ",counterBand);
            to_remove=posBand;
            if(posBand->nextBand!=NULL)
            {
                posBand=posBand->nextBand;
                removeBand(bL,to_remove);
            }
            else
            {
                removeBand(bL,to_remove);
                break;
            }
            counterBand++;
            continue;
        }
        int countSize=0;
        int old_scan=posBand->seqCentroid[0].scan;
        int current_scan=posBand->seqCentroid[0].scan;
        int maxSize = 0;
        int numInj = 0;
        //Loops on all the elements in the bands.
        for(i=0; i<posBand->size; i++)
        {
            //Rprintf("  %0.5f-%0.5f \n ",posBand->tPeaksMz[i].low_limit,posBand->tPeaksMz[i].high_limit);
            current_scan=posBand->seqCentroid[i].scan;
            if(current_scan-old_scan>1)
            {
                //Case where a signle scan is missing
                if(current_scan-old_scan==2)
                {
                    countSize++;

                }
                else
                {
                    countSize=0;

                }
            }
            else
            {
                countSize++;

            }
            if((current_scan>=injsc) & (current_scan<=injend))
            {
                numInj++;
            }

            old_scan=current_scan;
            if(countSize>maxSize)
            {
                maxSize=countSize;
            }
        }
        //Checking the faction of point present.
        if((maxSize<=size_min)&(numInj<limitInj))
        {
            //Rprintf(" %d ",counterBand);
            to_remove=posBand;
            posBand=posBand->nextBand;
            removeBand(bL,to_remove);
        }
        else
        {
            posBand=posBand->nextBand;
        }
        counterBand++;
    }
}


SEXP findBandsFIACentroids(SEXP mzVal, SEXP intVal, SEXP scanIndexVal,SEXP ScanTimeVal,SEXP firstScan, SEXP lastScan, SEXP maxScan, SEXP numMz, SEXP vPpm, SEXP nIso, SEXP injSc,SEXP endInj, SEXP vSizeMin, SEXP vDmz,SEXP rFracMin)
{
    double *intensity, *mz, mean_top, mean_bottom, ppm, *scantime, diffmz, frac_min;
    int *scanindex, lastscan, nummz, maxscan, firstscan,i,j,size_min, num_iso, injbeginning, injend;
    SEXP mat_to_return;
    //clock_t begin,end;
    mz=REAL(mzVal);
    intensity=REAL(intVal);
    scanindex=INTEGER(scanIndexVal);
    lastscan=INTEGER(lastScan)[0];
    num_iso=INTEGER(nIso)[0];
    maxscan=INTEGER(maxScan)[0];
    firstscan=INTEGER(firstScan)[0];
    injbeginning=INTEGER(injSc)[0];
    injend=INTEGER(endInj)[0];
    nummz=INTEGER(numMz)[0];
    ppm=REAL(vPpm)[0];
    scantime=REAL(ScanTimeVal);
    diffmz=REAL(vDmz)[0];
    frac_min=REAL(rFracMin)[0];
    //SEE IF THIS NEED TO BE PASSED IN PARAMETERS
    size_min=INTEGER(vSizeMin)[0];
    /** SECURITY CHECK **/
    if(lastscan>maxscan)
    {
        error("The maximum scan is superior to the scan limit.");
    }
    if(lastscan-firstscan<3)
    {
        error("Please pick more than 3 scans.");
    }
    struct scanFull * this_scan = calloc(1,sizeof(struct scanFull));
    this_scan->scan=NULL;
    struct bandList *bL=calloc(1,sizeof(struct bandList));
    //begin=clock();
    //Initializing the bandlist head.
    getScan(firstscan,mz,intensity,scanindex,maxscan,nummz,this_scan);
    struct centroid cmz;
    cmz.mz=this_scan->scan[0].mz;
    cmz.scan=1;
    cmz.intensity=this_scan->scan[0].intensity;
    struct band *nBand=calloc(1,sizeof(struct band));
    nBand->seqCentroid=calloc(5,sizeof(struct centroid));
    nBand->seqCentroid[0]=cmz;
    nBand->size=1;
    nBand->maxSize=5;
    nBand->nextBand=NULL;
    nBand->prevBand=NULL;
    nBand->meanMz=cmz.mz;
    bL->head=nBand;
    struct band* currentBand=nBand;
    bL->size=1;
    //Rprintf("BandList initialized : 0");
    //Creating the intial list of band with first scan
    for(j=1; j<this_scan->scanLength; j++)
    {
        cmz.mz=this_scan->scan[j].mz;
        cmz.scan=1;
        cmz.intensity=this_scan->scan[j].intensity;
        nBand=calloc(1,sizeof(struct band));
        nBand->seqCentroid=calloc(5,sizeof(struct centroid));
        nBand->seqCentroid[0]=cmz;
        nBand->size=1;
        nBand->maxSize=5;
        nBand->nextBand=NULL;
        nBand->prevBand=currentBand;
        nBand->meanMz=cmz.mz;
        currentBand->nextBand=nBand;
        currentBand=nBand;
        bL->size++;
    }
    for(i=firstscan+1; i<=lastscan; i++)
    {
        currentBand=bL->head;
        getScan(i,mz,intensity,scanindex,maxscan,nummz,this_scan);
        struct centroid cmz;
        for(j=0; j<this_scan->scanLength; j++)
        {

            cmz.mz=this_scan->scan[j].mz;
            cmz.scan=i;
            cmz.intensity=this_scan->scan[j].intensity;
            insertCentroidBandList(cmz,&currentBand,maxscan,firstscan,bL,ppm,diffmz);
        }
    }
    //TO DEBUG PRINTING ALL THE PEAKS BAND FOUND
    struct band *posBand=bL->head;
    free(this_scan->scan);
    free(this_scan);

    struct band * memBand;
//    while(posBand!=NULL)
//    {
//        //Rprintf("%0.5f %d  %d\n",posBand->meanMz,posBand->tPeaksMz[0].scan,posBand->tPeaksMz[posBand->size-1].scan);
//        posBand=posBand->nextBand;
//    }

//    Rprintf("IS:%d\n",size_bl(bL));

    cleanUpBandList(bL,3,injbeginning,injend,0);
//    Rprintf("AC:%d\n",size_bl(bL));
    int numFused=0;
    fuseBandList(bL,ppm,firstscan,lastscan,scantime,&numFused,diffmz);
//    Rprintf("AF:%d\n",size_bl(bL));
    cleanUpBandList(bL,size_min,injbeginning,injend,frac_min);
    //debuggingVizualisation(bL);
//    Rprintf("AC2:%d\n",size_bl(bL));
    double *mz_top=calloc(bL->size,sizeof(double));
    double *mz_bottom=calloc(bL->size,sizeof(double));
    int *scan_min=calloc(bL->size,sizeof(int));
    int *scan_max=calloc(bL->size,sizeof(int));
    int *num_points=calloc(bL->size,sizeof(int));
    double *intensity_v=calloc(bL->size,sizeof(double));
    int numFound = bL->size;
    posBand=bL->head;
    i=0;

    //Freeing the data structure.
    PROTECT(mat_to_return = allocMatrix(REALSXP,numFound,8));
    int solv;
    int countsol;
    int lastSc;
    int numS;
    double accu_mean;
    double maxInt;


    while((posBand!=NULL)&(i<numFound))
    {
        //Rprintf("%d ",i);
        solv=0;
        countsol=0;
        lastSc=0;
        numS=0;
        maxInt=0;
        accu_mean=0;
        mean_top=posBand->seqCentroid[posBand->size-1].mz;
        mean_bottom=posBand->seqCentroid[posBand->size-1].mz;
        //Rprintf("a%0.3f ",mean_bottom);
        for(j=0; j<posBand->size; j++)
        {
            //Rprintf("bb%0.3f ",posBand->seqCentroid[j].mz);
            if(maxInt<posBand->seqCentroid[j].intensity)
            {
                maxInt=posBand->seqCentroid[j].intensity;
            }
            //Determining if there is solvent.
            //Rprintf("%d:%d:%d:%d   ",injbeginning,posBand->seqCentroid[j].scan,countsol,j);
            if(posBand->seqCentroid[j].scan<injbeginning)
            {
                if(posBand->seqCentroid[j].scan-lastSc<=1)
                {
                    countsol++;
                    if(countsol>=num_iso)
                    {
                        solv=1;
                    }
                    numS++;
                }
                else
                {
                    countsol=0;
                    //Rprintf("No Sol");
                }
                accu_mean=accu_mean+posBand->seqCentroid[posBand->size-1].intensity;
                lastSc=posBand->seqCentroid[j].scan;
            }
            if(mean_top<posBand->seqCentroid[j].mz)
            {
                mean_top=posBand->seqCentroid[j].mz;
            }
            if(mean_bottom>posBand->seqCentroid[j].mz)
            {
                //Rprintf("bb%0.3f ",posBand->seqCentroid[j].mz);
                mean_bottom=posBand->seqCentroid[j].mz;
                //Rprintf("b%0.3f ",mean_bottom);
            }
        }

        if(solv)
        {
            REAL(mat_to_return)[i+numFound*6]=accu_mean/numS;
        }
        else
        {
            REAL(mat_to_return)[i+numFound*6]=0;
        }
        REAL(mat_to_return)[i+numFound*7]=maxInt;
        REAL(mat_to_return)[i+numFound*1]=mean_top;
        REAL(mat_to_return)[i]=mean_bottom;
        //Rprintf("c%0.3f ",mean_bottom);
        REAL(mat_to_return)[i+numFound*3]=posBand->seqCentroid[0].scan;
        REAL(mat_to_return)[i+numFound*4]=posBand->seqCentroid[posBand->size-1].scan;
        REAL(mat_to_return)[i+numFound*5]=posBand->size;
        REAL(mat_to_return)[i+numFound*2]=posBand->meanMz;

        //Rprintf("%0.5f %0.5f  %d  %d\n",mz_bottom[i],mz_top[i],scan_min[i],scan_max[i]);
        memBand=posBand;
        i++;
        posBand=posBand->nextBand;
        free(memBand->seqCentroid);
        free(memBand);
        //Rprintf("\n");
    }
    free(bL);
    free(mz_bottom);
    free(mz_top);
    free(scan_min);
    free(scan_max);
    free(num_points);
    free(intensity_v);
    UNPROTECT(1);
    return(mat_to_return);
}

double distPointLine2P(double xp,double yp,double x1,double y1,double x2, double y2)
{
    return(fabs(((y2-y1)*xp-(x2-x1)*yp+x2*y1-y2*x1))/sqrt((y2-y1)*(y2-y1)+(x2-x1)*(x2-x1)));
}

//The douglas peucker algoirthm.
int * segmentCurve(double * x,double *y, double epsilon,int np,int * sarray)
{

    if(np<=2)
    {
        *sarray=2;
        int * vres=malloc(sizeof(int)*2);
        vres[0]=0;
        vres[1]=np-1;
        return(vres);
    }

    double dist,dmax;
    int i,pmax;
    dmax=0;
    pmax=0;
    for(i=0; i<(np); i++)
    {
        dist=distPointLine2P(x[i],y[i],x[0],y[0],x[np-1],y[np-1]);
        //printf("dist:%d %0.0f  ",i,dist);
        if(dist>dmax)
        {
            dmax=dist;
            pmax=i;
        }
    }
    //printf("\n");
    if(dmax>epsilon)
    {


        int s1=0,s2=0;
        int * v1=segmentCurve(x, y, epsilon, pmax+1,&s1);
        int * v2=segmentCurve(&(x[pmax]),&(y[pmax]), epsilon, np-pmax,&s2);
        //printf("sum %d %d \n",s1,s2);
        int nsize=s1+s2-1;
        int * res = malloc(sizeof(int)*nsize);
        for(i=0; i<s1; i++)
        {
            res[i]=v1[i];
            //printf("%d ",res[i]);
        }
        for(i=1; i<s2; i++)
        {
            res[s1+i-1]=v2[i]+pmax;
            //printf("%d ",res[s1+i-1]);
        }
        *sarray=nsize;
        free(v1);
        free(v2);

        return(res);
    }
    else
    {
        *sarray=2;
        int * vres=malloc(sizeof(int)*2);
        vres[0]=0;
        vres[1]=np-1;
        return(vres);

    }
}

SEXP segmentCurveW(SEXP x, SEXP y, SEXP epsilon, SEXP n)
{
    SEXP vreturn;
    double *xx,*yy,*eps;
    int *seg,*num;
    num=INTEGER(n);
    xx=REAL(x);
    yy=REAL(y);
    eps=REAL(epsilon);
    int sizeArray=2;
    seg=segmentCurve(xx,yy,*eps,*num,&sizeArray);
    PROTECT(vreturn = allocVector(INTSXP,sizeArray));
    int i=0;
    for(i=0; i<sizeArray; i++)
    {
        INTEGER(vreturn)[i]=seg[i];
    }
    free(seg);
    UNPROTECT(1);
    return(vreturn);
}


static R_CallMethodDef callMethods[]  =
{
    {"findBandsFIACentroids", (DL_FUNC) &findBandsFIACentroids, 13},
    {"segmentCurveW", (DL_FUNC)&segmentCurveW, 4},
    {NULL, NULL, 0}
};

static R_NativePrimitiveArgType argfindLimDensity[] =
{
    REALSXP, INTSXP, INTSXP, INTSXP, INTSXP, INTSXP
};
static R_NativePrimitiveArgType argFindEqualGreaterM[] =
{
    REALSXP, INTSXP, REALSXP, INTSXP, INTSXP
};
static R_NativePrimitiveArgType arglineAboveNoise[] =
{
    INTSXP, INTSXP, INTSXP, REALSXP, INTSXP, INTSXP, REALSXP
};
static R_NativePrimitiveArgType argfindLim[] =
{
    REALSXP, INTSXP, INTSXP, INTSXP, INTSXP, INTSXP
};
static R_NativePrimitiveArgType argbinarySearch[] =
{
    REALSXP, REALSXP, INTSXP, INTSXP, INTSXP
};
static R_NativePrimitiveArgType argreplaceByZero[] =
{
    REALSXP, REALSXP, INTSXP
};
static R_NativePrimitiveArgType arglinearInterpolation[] =
{
    REALSXP, INTSXP, INTSXP, INTSXP
};
static R_NativePrimitiveArgType argcheckIso[] =
{
    REALSXP, INTSXP, INTSXP, INTSXP
};

static R_CMethodDef cMethods[] =
{
    {"findLimDensity", (DL_FUNC) &findLimDensity, 6, argfindLimDensity},
    {"FindEqualGreaterM", (DL_FUNC) &FindEqualGreaterM, 5, argFindEqualGreaterM},
    {"lineAboveNoise", (DL_FUNC) &lineAboveNoise, 7, arglineAboveNoise},
    {"findLim", (DL_FUNC) &findLim, 6, argfindLim},
    {"binarySearch", (DL_FUNC) &binarySearch, 5, argbinarySearch},
    {"replaceByZero", (DL_FUNC) &replaceByZero, 3, argreplaceByZero},
    {"linearInterpolation", (DL_FUNC) &linearInterpolation, 4, arglinearInterpolation},
    {"checkIso", (DL_FUNC) &checkIso, 4, argcheckIso},
    {NULL, NULL, 0}
};

void R_init_myLib(DllInfo *info)
{
    R_registerRoutines(info, cMethods, callMethods, NULL, NULL);
}
