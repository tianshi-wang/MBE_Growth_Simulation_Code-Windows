#import "solver.h"

#define kT (.026*temperature/300.0*q)
#define q 1.6e-19
#define hbar 1.055e-34
#define planck 6.02e-34
double freq, bondEnergy, hoppingEnergy, desorbEnergy, currTime, temperature;
double bondRate, hoppingRate, arrivalRate;
int steps;
int height[L][L];
int update[L][L];
double ratetotals[DEPTH][L][L];
double rates[L][L][6];
double bonds[5];
double desorbRates[5];
int arraySize;
double speedratio;
bool valueChanged = false;



unsigned int myrand();             /* returns a myrandom 32 bit integer */
void myrand_seed( unsigned int, unsigned int, unsigned int );

/* return a myrandom float >= 0 and < 1 */
#define myrand_float       ((double)myrand()/4294967296.0)

static unsigned int s1=390451501, s2=613566701, s3=858993401;  /* The seeds */
static unsigned mask1, mask2, mask3;
static int shft1, shft2, shft3, k1=31, k2=29, k3=28;

/* use either of the following two sets of parameters*/
static int q1=13, q2=2, q3=3, p1=12, p2=4, p3=17;
/* static int q1=3, q2=2, q3=13, p1=20, p2=16, p3=7; */
double myrandv;




unsigned int myrand()
{
    unsigned int b;

    b  = ((s1 << q1)^s1) >> shft1;
    s1 = ((s1 & mask1) << p1) ^ b;
    b  = ((s2 << q2) ^ s2) >> shft2;
    s2 = ((s2 & mask2) << p2) ^ b;
    b  = ((s3 << q3) ^ s3) >> shft3;
    s3 = ((s3 & mask3) << p3) ^ b;
    return (s1 ^ s2 ^ s3);
}


void myrand_seed( unsigned int a, unsigned int b, unsigned int c )
{
    static unsigned int x = 4294967295U;

    shft1 = k1-p1;
    shft2 = k2-p2;
    shft3 = k3-p3;
    mask1 = x << (32-k1);
    mask2 = x << (32-k2);
    mask3 = x << (32-k3);
    if (a > (1 << (32-k1))) s1 = a;
    if (b > (1 << (32-k2))) s2 = b;
    if (c > (1 << (32-k3))) s3 = c;
    myrand();
}

int calcBonds(int i, int j, int heightadj)
{
    int bonds = 0;
    j=(j+arraySize)%arraySize;
    i=(i+arraySize)%arraySize;
    int currheight = height[(i+arraySize)%arraySize][(j+arraySize)%arraySize] + heightadj;

    if (height[(i-1+arraySize)%arraySize][j] + (i==0?-steps:0) >= currheight) bonds++;
    if (height[(i+1+arraySize)%arraySize][j] + (i==arraySize-1?steps:0) >= currheight) bonds++;
    if (height[i][(j-1+arraySize)%arraySize] >= currheight) bonds++;
    if (height[i][(j+1+arraySize)%arraySize] >= currheight) bonds++;

//    printf( "currheight is %d  atombonds is %d + at (i,j) %d, %d\n", currheight, bonds, i, j);


    return bonds;
}

void calcLocalRate(int i, int j)
{
    i=(i+arraySize)%arraySize;
    j=(j+arraySize)%arraySize;
    int currbonds = calcBonds(i,j,0);

    int newbonds;
    int k;
//    printf( "currheight is %d   at (i,j) %d, %d\n", height[i][j], i, j);

    height[i][j]--;
//    printf( "currheight is %d   at (i,j) %d, %d\n", height[i][j], i, j);

    newbonds = calcBonds(i, j+1, 1);

    //printf("newbonds  is %d, at i= %d j=%d\n", newbonds, i,j);

    if (currbonds > newbonds)
        rates[i][j][0] = bonds[currbonds-newbonds];
    else
        rates[i][j][0] = bonds[0];

    newbonds = calcBonds(i, j-1, 1);
    if (currbonds > newbonds)
        rates[i][j][1] = bonds[currbonds-newbonds];

    else
        rates[i][j][1] = bonds[0];

    newbonds = calcBonds(i-1, j, 1);
    if (currbonds > newbonds)
        rates[i][j][2] = bonds[currbonds-newbonds];
    else
        rates[i][j][2] = bonds[0];

    newbonds = calcBonds(i+1, j, 1);
    if (currbonds > newbonds)
        rates[i][j][3] = bonds[currbonds-newbonds];
    else
        rates[i][j][3] = bonds[0];
    height[i][j]++;

    rates[i][j][4] = desorbRates[currbonds];

    rates[i][j][5] = arrivalRate;
    ratetotals[DEPTH-1][i][j]=0;
    for (k=0;k<6;k++)
        {ratetotals[DEPTH-1][i][j]+=rates[i][j][k];
        //if (k==5) printf("At k=%d, self.ratetotals is: %f:\n", k, ratetotals[DEPTH-1][i][j]);

    }

    update[i][j] = 1;
}

void updateTree(int i, int j)
{
    if (update[i][j]==0)
        return;
    int l,m,k;
    l=i/2;
    m=j/2;
    update[l*2][m*2]=0;
    update[l*2][m*2+1]=0;
    update[l*2+1][m*2]=0;
    update[l*2+1][m*2+1]=0;
    for (k=DEPTH-2; k>=0; k--) {
        ratetotals[k][l][m] = ratetotals[k+1][l*2][m*2] + ratetotals[k+1][l*2][m*2+1] +
                              ratetotals[k+1][l*2+1][m*2] + ratetotals[k+1][l*2+1][m*2+1];
        l/=2;
        m/=2;
    }
}

void initArrays()
{
    int i, j, k;
    for (i=0;i<arraySize;i++) {
        for (j=0; j<arraySize; j++) {
            height[i][j] = i*steps/(arraySize); //steps
            //	height[i][j]=0;
            //	if (i>arraySize/2)
            //		height[i][j]=arraySize/4-i/(arraySize/8);
            //	else
            //		height[i][j]=i/(arraySize/8);
//            printf("At i,j, (%d,%d),height[i][j] is %d \n", i, j, height[i][j]);

            update[i][j] = 0;
            for (k=0;k<DEPTH;k++)
                ratetotals[k][i][j]=0.0;
        }
   }

//    printf("initArrays done\n");

}
void updateRates()
{
    int i;
    int j;
    double desorbRate;
    freq = kT/planck;
    //freq = 1e13;
    bondRate = exp(-bondEnergy*q/kT);
    hoppingRate = freq*exp(-hoppingEnergy*q/kT);
    desorbRate = freq*exp(-desorbEnergy*q/kT);
//    printf("bondRate is %f \n",bondRate);
//    printf("freq is %f \n",freq);
//    printf("bondEnergy is %f \n", bondEnergy);

//	desorbRate = 0.0;
//	hoppingRate = 0.0;
//	bondRate = 0.0;
    desorbRates[0] = desorbRate;
    bonds[0]=hoppingRate;
    for (i=1;i<=4;i++) {
        bonds[i]=bonds[i-1]*bondRate;
        desorbRates[i] = desorbRates[i-1]*bondRate;
//        printf("Bond[%d] is %f", i,  bonds[i]);
    }

    for (i=0;i<arraySize;i++) {
        for (j=0; j<arraySize; j++) {
            calcLocalRate(i,j);
        }
    }
    for (i=0;i<arraySize;i++) {
        for (j=0; j<arraySize; j++) {
            updateTree(i,j);
        }
    }


}

void calcLocalRates(int i, int j)
{
    calcLocalRate((i+arraySize)%arraySize,(j+arraySize)%arraySize);

    if (!update[(i+arraySize+1)%arraySize][(j+arraySize+1)%arraySize])
        calcLocalRate((i+arraySize+1)%arraySize,(j+arraySize+1)%arraySize);
    if (!update[(i+arraySize-1)%arraySize][(j+arraySize+1)%arraySize])
        calcLocalRate((i+arraySize-1)%arraySize,(j+arraySize+1)%arraySize);
    if (!update[(i+arraySize+1)%arraySize][(j+arraySize-1)%arraySize])
        calcLocalRate((i+arraySize+1)%arraySize,(j+arraySize-1)%arraySize);
    if (!update[(i+arraySize-1)%arraySize][(j+arraySize-1)%arraySize])
        calcLocalRate((i+arraySize-1)%arraySize,(j+arraySize-1)%arraySize);
    if (!update[(i+arraySize)%arraySize][(j+arraySize+1)%arraySize])
        calcLocalRate((i+arraySize)%arraySize,(j+arraySize+1)%arraySize);
    if (!update[(i+arraySize+1)%arraySize][(j+arraySize)%arraySize])
        calcLocalRate((i+arraySize+1)%arraySize,(j+arraySize)%arraySize);
    if (!update[(i+arraySize-1)%arraySize][(j+arraySize)%arraySize])
        calcLocalRate((i+arraySize-1)%arraySize,(j+arraySize)%arraySize);
    if (!update[(i+arraySize)%arraySize][(j+arraySize-1)%arraySize])
        calcLocalRate((i+arraySize)%arraySize,(j+arraySize-1)%arraySize);
}
void setupSteps(int newSteps)
{
    if (newSteps == steps)
        return;

    int i, j, k;
    for (i=0;i<arraySize;i++) {
        for (j=0; j<arraySize; j++) {
            height[i][j] += (i*(newSteps))/arraySize-(i*steps)/arraySize;
            //	height[i][j]=0;
            //	if (i>arraySize/2)
            //		height[i][j]=arraySize/4-i/(arraySize/8);
            //	else
            //		height[i][j]=i/(arraySize/8);

            update[i][j] = 0;
            for (k=0;k<DEPTH;k++)
                ratetotals[k][i][j]=0.0;
        }
    }
    steps = newSteps;


}
void updateLocalTree(int i, int j)
{
    updateTree((i-1+arraySize)%arraySize,(j-1+arraySize)%arraySize);
    updateTree((i+1+arraySize)%arraySize,(j+1+arraySize)%arraySize);
    updateTree((i-1+arraySize)%arraySize,(j+1+arraySize)%arraySize);
    updateTree((i+1+arraySize)%arraySize,(j-1+arraySize)%arraySize);
}

double step(double val)
{
    double R, sumr;
    int i,l,m;
    R=0.0;
    l=0;
    m=0;

    R=ratetotals[0][0][0]+ratetotals[0][1][0]+ratetotals[0][0][1]+ratetotals[0][1][1];
    val=val*R;

    for (i=0;i<DEPTH-1; i++) {
        sumr = 0.0;
        if (val<0)
            printf("whoopps");
        if (val>ratetotals[i][l][m]+ratetotals[i][l+1][m]+ratetotals[i][l][m+1]+ratetotals[i][l+1][m+1])
            printf("got a problem here");
        if (val>sumr+ratetotals[i][l][m]) {
            sumr=ratetotals[i][l][m];
            if (val>sumr+ratetotals[i][l+1][m]) {
                sumr+=ratetotals[i][l+1][m];
                if (val>sumr+ratetotals[i][l+1][m+1]) {
                    sumr+=ratetotals[i][l+1][m+1];
                    l*=2;
                    m=(m+1)*2;
                    val -= sumr;
                }
                else {
                    l=(l+1)*2;
                    m=(m+1)*2;
                    val -= sumr;
                }
            }
            else {
                m*=2;
                l=(l+1)*2;
                val -= sumr;
            }
        }
        else {
            m*=2;
            l*=2;
        }

    }
    sumr = 0.0;
    if (val>sumr+ratetotals[i][l][m]) {
        sumr=ratetotals[i][l][m];
        if (val>sumr+ratetotals[i][l+1][m]) {
            sumr+=ratetotals[i][l+1][m];
            if (val>sumr+ratetotals[i][l+1][m+1]) {
                sumr+=ratetotals[i][l+1][m+1];

                m=m+1;
                val -= sumr;
            }
            else {
                l=l+1;
                m=m+1;
                val -= sumr;
            }
        }
        else {
            l=l+1;
            val -= sumr;
        }
    }
    i=0;
    sumr = rates[l][m][i];
    int move;
    int nl, nm;
    int dc, dn;
    dc=0;
    if (val > sumr) {
        i++;
        sumr += rates[l][m][i];
        if (val > sumr) {
            i++;
            sumr += rates[l][m][i];
            if (val > sumr) {
                i++;
                sumr += rates[l][m][i];
                if (val > sumr) {
                    i++;
                    sumr += rates[l][m][i];
                    if (val > sumr) {

                        //adsorb
//						printf("should be no adsorption\n");
                        move = 0;
                        dc = 1;
                    }
                    else {
                        //desorb
                        //printf("should be no desorption\n");
                        move = 0;
                        dc = -1;
                    }
                }
                else {
                    //right
                    move = 1;
                    dc = -1;
                    dn = 1;
                    nl = l+1;
                    nm = m;
                }
            }
            else {
                //left
                move = 1;
                dc = -1;
                dn = 1;
                nl = l-1;
                nm = m;
            }
        }
        else {
            //down
            move = 1;
            dc = -1;
            dn = 1;
            nl = l;
            nm = m-1;
        }
    }
    else {
        // up
        move = 1;
        dc = -1;
        dn = 1;
        nl = l;
        nm = m+1;
    }
    height[l][m]+=dc;


    if (move == 1) {
        height[(nl+arraySize)%arraySize][(nm+arraySize)%arraySize] += dn;
        calcLocalRates(nl,nm);
    }

    calcLocalRates(l,m);

    updateLocalTree(l,m);

    if (move==1)
        updateLocalTree(nl,nm);
//    printf("l,m and height[l][m]: %d, %d, %d\n", l, m, height[l][m]);

    return R;

}

void buildTree()
{
    int i;
    for (i=0;i<arraySize;i++) {
        for (int j=0; j<arraySize; j++) {
            calcLocalRate(i,j);
        }
    }
    for (i=0;i<arraySize;i++) {
        for (int j=0; j<arraySize; j++) {
            updateTree(i,j);
        }
    }
}





Solver::Solver(Results* myresult)
{
    results = myresult;
}


void Solver::run()
{
    QTime* runtimer = new QTime;
    runtimer->start();
    int i, j;
    currTime = results->time;
    arraySize = results->data->getSize();
    hoppingEnergy = results->hoppingEnergy;
    bondEnergy = results->bondEnergy;
    desorbEnergy = results->desorptionEnergy;
    speedratio = results->speedratio;


    arrivalRate = results->arrivalRate;
    temperature = results->temperature;
    steps = results->steps;

    myrand_seed(s1,s2,s3);
    initArrays();
    int iterations=0;
    int loopiterations=0;
//    printf("kT is %.12e\n:",kT);
    double runTimeLoopBegin;
    double growthTimeLoopBegin;
    bool resetRuntime = true;

//    double sleepingTime = 0.0; //sleepingTime is the total time the thread fall asleep. Use it to correct run_time.
    while (true) {
        // stop is true after clicking pause button
        if (results->stop){
            results->stop = false;
            results->resume = false;

            while(results->resume==false && results->reset==false){
                this->msleep(100); }
//                sleepingTime += 0.1; }
           }

        if (resetRuntime == true){        //this function to calculate timescale with
            resetRuntime=false;
            runtimer->restart();
            growthTimeLoopBegin = results->time;
        }
        arraySize = results->data->getSize();
        hoppingEnergy = results->hoppingEnergy;
        bondEnergy = results->bondEnergy;
        desorbEnergy = results->desorptionEnergy;
        arrivalRate = results->arrivalRate;
        temperature = results->temperature;
        speedratio = results->speedratio;
        steps=results->steps;
//        printf("desorbEnergy is :%f",desorbEnergy);

        //reset: activated by clicking the start_over button
        if (results->reset ) {
            results->reset = false;
            currTime = 0.0;
//            sleepingTime = 0.0;
            runtimer->restart();
            initArrays();
            results->initHeight();
            steps=results->steps;
            resetRuntime = true;     //

        }
        if (results->resetTimer) {
            currTime = 0.0;
            results->time = 0.0;
            results->resetTimer = false;
        }
        setupSteps(results->steps);


        updateRates();

        double  R;
        for (i=0;i<2000;i++) {
            currTime-=log(myrand_float)/step(myrand_float);
            iterations ++;
        }

        if (!results->reset) {
        for (i=0;i<arraySize;i++)
            for (j=0;j<arraySize; j++)
                results->data->setAt(i,j,height[i][j]);
            //results->data->setAt(0,0,10);
            //results->data->setAt(arraySize-1 ,0,10);

      }

        if (runtimer->elapsed()>500){               //update the "timescale" every 0.5 sec. Here we use runtime instead of iteration number for refreshing.
            results->timescale = (currTime-growthTimeLoopBegin)/(runtimer->elapsed()*0.001); //in second.
            resetRuntime = true;     //
        }
\
        results->time = currTime;   //update growth time in each loop.

        Solver::sleep(speedratio);  //Control run speed by setting a thread sleep time for each loop.
        loopiterations++;

     }



}
