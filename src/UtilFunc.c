#include <math.h>
#include <R.h>

int sign_d(double x)
{
    return((x>0)-(x<0));
}

void lineAboveNoise(int * beginning,int * posLeft,int * posRight,double * intval,int * sizeMz,int * size_min,double * noise)
{
    //Rprintf("line above noise %d %d %d %f\n",*sizeMz,*beginning,*size_min,*noise);
    int i=0;
    int currentSize=0;
    //Rprintf("\n%d and %d\n",*beginning,*sizeMz-1);
    if(*beginning==*sizeMz-1)
    {
        return;
    }
    for(i=*beginning; i<*sizeMz; i++)
    {
        //Rprintf("%d_%d ",i,currentSize);
        if(intval[i]>*noise)
        {
            currentSize++;
        }
        else
        {
            if(currentSize>*size_min)
            {
                (*beginning)=*posRight;
                (*posLeft)=i-currentSize;
                (*posRight)=i;
                break;
            }
            currentSize=0;
        }
    }
    if(currentSize>*size_min)
    {
        (*beginning)=*posRight;
        (*posLeft)=i-currentSize;
        (*posRight)=i;
    }
    //Rprintf("Zone found : %d - %d \n",*posLeft,*posRight);
}

void findingPeakLimit(int * posLeft, int * posRight,double * intval,int *sizeMz)
{
    //Rprintf("finding Peak Limit\n");
    while(intval[(*posLeft)]<intval[((*posLeft)+1)]&&(*posLeft)>0)
    {
        (*posLeft)--;
    }
    while(intval[(*posRight-1)]>intval[(*posRight)]&&(*posRight)<*sizeMz-1)
    {
        (*posRight)++;
    }
}



void findLimDensity(double * seq,int * size,int * istart,int * linflex, int * rinflex,int * state)
{
    //Rprintf("In the function \n");
    //Rprintf("here \n");
    int i,difference;
    *linflex=*istart;
    *rinflex=*istart;
    for(i=*istart; i<*size-1; i++)
    {
        difference=sign_d(seq[i+1]-seq[i])-sign_d(seq[i]-seq[i-1]);
        if(difference>=1)  //inflex point
        {
            if(*state==2)
            {
                *state=1;
                *rinflex=i;
                *istart=i;
                break;
            }
            else//First encounter with an inflex
            {
                *state=1;
                *linflex=i;
            }
        }
        if(difference==-2&&(*state)==1)
        {
            *state=2;
        }
    }
    //Rprintf("l %d and r %d and state %d\n",*linflex,*rinflex,*state);
}

void FindEqualGreaterM(const double *in, const int *size, const double *values,
                       const int *valsize, int *index)
{

    int i, idx = 0;

    for (i = 0; i < *valsize; i++)
    {
        while (idx < *size && in[idx] < values[i])
            idx++;
        index[i] = idx;
    }
}


void findLim(double *yseq,int *istartl,int *istartr,int *lengthmax, int *ileft, int *iright)
{
//Finding the left limit
//Rprintf("entering findLim %d %d %d\n ",*istartl,*istartr,*lengthmax);
    int i;
    for (i = *istartl; i > 1; i--)
    {
        //Rprintf("val i : %d \n",i);
        if (yseq[i-1] >= yseq[i])
            break;
    }
    *ileft = i;
    //Rprintf("left ok \n");
    for (i = *istartr; i < *lengthmax-2; i++)
        if (yseq[i+1] >= yseq[i])
            break;
    *iright = i;
    //Rprintf("right ok \n");
}


//Locate the closest local Max.
void findMax(double *yseq,int *istart,int *res)
{
    int i=*istart-1;
    double dLeft,dRight;
    dLeft=yseq[*istart-1]+yseq[*istart-2];
    dRight=yseq[*istart+1]+yseq[*istart+2];
    while((dLeft>=2*yseq[i])|(yseq[i+1]+yseq[i+2]>=2*yseq[i]))
    {
        dLeft=yseq[*istart-1]+yseq[*istart-2];
        dRight=yseq[*istart+1]+yseq[*istart+2];
        if(dLeft>dRight)
        {
            i--;
        }
        else
        {
            i++;
        }
    }
    *res=i;
}


//Finding the closest value to a number.
void binarySearch(double *xseq, double *val, int *imin, int *imax, int *imid)
{
    int vmin=*imin;
    int vmax=*imax;
    //Rprintf("entering binarySearch %f %d %d \n",*val,min,max);
    while (vmax-vmin>1)
    {
        // calculate the midpoint for roughly equal partition
        //Rprintf("div 2");
        *imid = (vmin+vmax) >> 1;
        //Rprintf("min : %d max : %d mid : %d \n",min,max,*imid);
        if (xseq[*imid] == *val)
        {
            *imin = *imid;
            *imax = *imid;
            return;
        }
        else if (xseq[*imid] < *val)
            vmin = *imid;
        else
            vmax = *imid;
    }
    //Then interval is of length 1
    if(fabs(*val-xseq[vmin])<fabs(xseq[vmax]-*val))
    {
        *imid=vmin;
    }
    else
    {
        *imid=vmax;
    }
}


void replaceByZero(double *rseq,double *yseq,int* sizex)
{
    int i;
    for(i=0; i<*sizex; i++)
    {
        if(yseq[i]<=0)
        {
            rseq[i]=0;
        }
    }
}

void linearInterpolation(double *intseq,int *imin, int *imax, int *counter)
{
    int i;
    *counter=0;
    for(i=*imin+1; i<*imax+1-1; i++)
    {
        if((intseq[i-1]>0)&(intseq[i+1]>0)&(intseq[i]==0))
        {
            intseq[i]=(intseq[i-1]+intseq[i+1])/2;
            (*counter)++;
        }
    }
}

//hotfix to prevent randoom crash.
void checkIso(double *intensity, int *num_iso, int *valid, int *smax)
{
    int i;
    int counter=0;
    *valid=1;
    if(*smax == 0){
        if(intensity[0] != 0){
            *valid = 0;
            return;
        }
    }
    for(i=0; i<*smax; i++)
    {
        if(intensity[i]>0)
        {
            counter++;
            if(counter>=*num_iso)
            {
                *valid=0;
                return;
            }
        }
        else
        {
            counter=0;
        }
    }
}
