#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include <stdbool.h>
#include <unistd.h>
#include "pfs.h"
//#include "motifproj.c"

#ifndef M_PI
    #define M_PI 3.14159265358979323846
#endif

double ESCALE,CISCALE;
int TOTCONST;

char CWD[256]; 

double epsRR,epsRS,epsSS;
double COMRP[3];
double beta,epsilon;
int iter,CHstep,totiter,maxiter,HB1,HB2;
double RM[3][3], EGS, RAD, EBT, EMIN, BRAD;
double (**ActualCoor)[3], (*ActualSolCoor)[3][3];

//double tau;

char *ProID;
int AA[91],lenAA[21], lenAAH[21],*ATR[21],NumSegs[21],*SegSize[21],**SegAtom[21];
int *lenRP,**ATRP,*NumSegsRP,**SegSizeRP,***SegAtomRP, **SegIntRP, *SegInt[21];
double **PQ,*PQAA[3][21],PQS[3],(**parRP)[3];
int ProLen,NumSolEx,TotAtoms,NumSol,*RT,NumCys,*CysIndex, NumAD, (*ADBK)[8][2];
int (*HBBK[2])[2], NumDN, NumAC;
// int (*ADkey)[3], (**ADCheck)[3];
short int BBMotifNum[11], *SCMotifNum[21], SSMotifNum, CTMotifNum[21];
double NTMotif[5][3], OXTMotif[4][3], SolMotif[3][3];

double (***SCStore)[3];
double *WeightsCT;
double **WeightsHB, (**HBStore)[4][3];
double **WeightsSS, (**SSStore)[4][3];
double *WeightsAD,(*ADStore)[8][3];

double (**CTMotifs[21])[3],(***SCMotifs[21])[3], (*SSMotifs)[4][3], ADMotif[4][3];

double (**HB)[4][3], (**HBA)[4][3], (**HBR)[4][3], (**HBB)[4][3];
double (**DH)[4][3], (**DHR)[4][3], (**DHA)[4][3], (**DHB)[4][3]; 
double (**SS)[4][3], (**SSR)[4][3], (**SSA)[4][3], (**SSB)[4][3]; 
double (**CT)[3], (**CTA)[3], (**CTR)[3], (**CTB)[3];
double (*AD)[8][3], (*ADA)[8][3], (*ADR)[8][3], (*ADB)[8][3];
double (***SC)[3], (***SCA)[3], (***SCR)[3], (***SCB)[3];
double (****VB)[2][3], (****VBA)[2][3], (****VBB)[2][3], (****VBR)[2][3];
double (***VS)[3][2][2][3], (***VSA)[3][2][2][3], (***VSB)[3][2][2][3], (***VSR)[3][2][2][3];
double (**VO)[3][3][2][3], (**VOA)[3][3][2][3], (**VOB)[3][3][2][3], (**VOR)[3][3][2][3];
double (***ES)[3][2], (***ESA)[3][2], (***ESB)[3][2], (***ESR)[3][2];
double (**EO)[3][3][2], (**EOA)[3][3][2], (**EOR)[3][3][2], (**EOB)[3][3][2];
double (****EB)[2], (****EBA)[2], (****EBB)[2], (****EBR)[2];
double (*OS)[3][3], (*OSA)[3][3], (*OSR)[3][3], (*OSB)[3][3];
double **ED, **EDA, **EDR, **EDB;
double (*BL)[2][3], (*BLA)[2][3], (*BLR)[2][3], (*BLB)[2][3];
double *BE, *BEA, *BER, *BEB;

double (****INTRP)[3];
double (**parRP)[3],parSol[3][3];
double (**atom)[3], (*solatom)[3][3];
int (*A2R)[2], **R2A, ****CONH;
double ****VXRR, (**VXRS)[3], VXSS[3][3];
int NumBL, (*BLAtoms)[2];
double (*BLpar)[2];
int NumDH, (*DHAtoms)[4], *NumDHpar;
double (**DHpar)[3];

double toterr;
double **etaHB, **HBerr, tHBerr;
double *etaAD, *ADerr, tADerr;
double **etaSS, **SSerr, tSSerr;
double **etaSC, **SCerr, tSCerr;
double *etaCT, *CTerr, tCTerr;
double *etaOS, *OSerr, tOSerr;
double *etaBL, *BLerr, tBLerr;
double **etaDH, **DHerr, tDHerr;
double (****etaEB)[2], (****EBerr)[2], tEBerr[2];
double (***etaES)[3][2], (***ESerr)[3][2], tESerr[2];
double (**etaEO)[3][3][2], (**EOerr)[3][3][2], tEOerr[2];

double (*LJP)[2];

char errfile[50],statsfile[50],solfile[50],ProFile[50],etafile[50],Enerfile[50],initfile[50];

int (*val_index)[9], max_size;
double *val;


int valinc(const void *a,const void *b)
        {
        double arg1 = val[*(const int*)a];
        double arg2 = val[*(const int*)b];

        return (arg1 > arg2) - (arg1 < arg2);
        }


double urand(double a, double b)
        {
        double num;
        num = (b-a)*(((double)rand())/RAND_MAX)+a;
        return num;
        }


double RELU(double a)
        {
        if (a < 0.)
                return 0.;
        else
                return a;
        }

int compare_error_asc_G( const void *pa, const void *pb)
        {
        const double *a = pa;
        const double *b = pb;

        if ( *a < *b ) return -1;
        if ( *a > *b ) return +1;
        return 0;
        }


int compare_error_asc( const void *pa, const void *pb )
        {
        const double (*a)[2] = pa;
        const double (*b)[2] = pb;
        if ( (*a)[1] < (*b)[1] ) return -1;
        if ( (*a)[1] > (*b)[1] ) return +1;
        return 0;
        }


int compare_error_asc2( const void *pa, const void *pb )
        {
        const double (*a)[3] = pa;
        const double (*b)[3] = pb;
        if ( (*a)[2] < (*b)[2] ) return -1;
        if ( (*a)[2] > (*b)[2] ) return +1;
        return 0;
        }


static inline double sq(double diff)
        {
        return diff*diff;
        }

static inline double norm(double u[3])
        {   
        return u[0]*u[0]+u[1]*u[1]+u[2]*u[2];
        }  

static inline double dot_prod(double u[3],double v[3])
        {   
        return u[0]*v[0]+u[1]*v[1]+u[2]*v[2];
        }

static inline double trip_prod(double n[3], double u[3],double v[3])
        {
        return n[0]*(u[1]*v[2]-u[2]*v[1]) + n[1]*(u[2]*v[0]-u[0]*v[2]) + n[2]*(u[0]*v[1]-u[1]*v[0]);
        }

void RotationMatrix(double ax[3], double theta)
	{
	RM[0][0] = cos(theta) + sq(ax[0])*(1-cos(theta));
	RM[0][1] = ax[0]*ax[1]*(1-cos(theta)) - ax[2]*sin(theta);
	RM[0][2] = ax[0]*ax[2]*(1-cos(theta)) + ax[1]*sin(theta);	
	RM[1][0] = ax[1]*ax[0]*(1-cos(theta)) + ax[2]*sin(theta);
	RM[1][1] = cos(theta) + sq(ax[1])*(1-cos(theta));
	RM[1][2] = ax[1]*ax[2]*(1-cos(theta)) - ax[0]*sin(theta);
	RM[2][0] = ax[2]*ax[0]*(1-cos(theta)) - ax[1]*sin(theta);
	RM[2][1] = ax[2]*ax[1]*(1-cos(theta)) + ax[0]*sin(theta);
	RM[2][2] = cos(theta) + sq(ax[2])*(1-cos(theta));
	}

double calc_dihed( double at[4][3])
	{
	int n;
	double vecA[3], vecD[3], normal[3], perpCBA[3], perpDCB[3];
	double dotp, trpp;

	for(n=0;n<3;++n)
		{
		vecA[n] = at[0][n] - at[1][n];
		vecD[n] = at[3][n] - at[2][n];
		normal[n] = at[1][n] - at[2][n];
		}

	cross(vecD,normal,perpDCB);
	cross(vecA,normal,perpCBA);

	normalize(normal);
	dotp = dot_prod(perpCBA,perpDCB); 
	trpp = trip_prod(normal,perpDCB,perpCBA);
	
	return atan2(trpp,dotp);
	}

void proj_bondNE( double at[2][3], double new_at[2][3], double parm[2])
	{
	double dis = 0.,proj_dis;
	int n;

	for(n=0;n<3;++n)
		dis += sq(at[0][n] - at[1][n]);

	dis = sqrt(dis);

	proj_dis = parm[1];

	for(n=0;n<3;++n)
		{
		new_at[0][n] = at[0][n] + ((proj_dis-dis)/(2*dis))*(at[0][n]-at[1][n]);
		new_at[1][n] = at[1][n] - ((proj_dis-dis)/(2*dis))*(at[0][n]-at[1][n]);
		}
	}

double proj_bond2( double at[2][3], double new_at[2][3], double Ein, double parm[2], double Emax, int count)
	{
	double dis = 0.,proj_dis,x,maxx,delx,dely,ms,xc,yc;
	int n;

	for(n=0;n<3;++n)
		dis += sq(at[0][n] - at[1][n]);

	dis = sqrt(dis);

	double xi = dis/sqrt(2);
	double yi = Ein;
	double bo = parm[1]/sqrt(2);
	double tempxy[2][2];

	xc = xi;

	for (int i = 0; i<count; ++i)
		{
		if (yi <  0)
			{
			xc = xi;
			yc = ESCALE*parm[0]*sq(xc - bo);
			}
		else
			{
			tempxy[0][0] = xi;
			tempxy[0][1] = ESCALE*parm[0]*sq(xc - bo);

			x = sqrt(yi/(ESCALE*parm[0]));

			if (xi < bo)
				tempxy[1][0] = bo - x;
			else
				tempxy[1][0] = bo + x;

			tempxy[1][1] = yi;

			if (fabs(tempxy[1][0] - xi) > fabs(tempxy[0][1] - yi))
				{
				xc = xi;
				yc = tempxy[0][1];
				}
			else
				{
				xc = tempxy[1][0];
				yc = yi;
				}
			}

		ms = ESCALE*2*parm[0]*(xc - bo);
		delx = (xi - xc + ms*(yi - yc))/(1+sq(ms));
		dely = ms*delx;

		xi = xc + delx;
		yi = yc + dely;
		}

	proj_dis = xc*sqrt(2);

	maxx = sqrt(2.*Emax/parm[0]);
	if (maxx > parm[1])
		maxx = parm[1];

	if (proj_dis < parm[1] - maxx) 
		proj_dis = parm[1] - maxx;

	if (proj_dis > parm[1] + maxx) 
		proj_dis = parm[1] + maxx;

	for(n=0;n<3;++n)
		{
		new_at[0][n] = at[0][n] + ((proj_dis-dis)/(2*dis))*(at[0][n]-at[1][n]);
		new_at[1][n] = at[1][n] - ((proj_dis-dis)/(2*dis))*(at[0][n]-at[1][n]);
		}

	double projE;
	projE = ESCALE*0.5*parm[0]*sq(parm[1] - proj_dis);

	return projE;
	}
/*
double proj_bond( double at[2][3], double new_at[2][3], double Ein, double parm[2], double Emax)
	{
	double dis = 0.,proj_dis,minweight,tempweight,x,dx,twp,twm,maxdx,new_dx;
	int n;
	double E1,E2,eng;

	for(n=0;n<3;++n)
		dis += sq(at[0][n] - at[1][n]);

	dis = sqrt(dis);

	minweight = sq(Ein) + 0.5*sq(dis - parm[1]);
	proj_dis = parm[1];
	E1 = 0.;

	dx = 0.01;

	maxdx = sqrt(2.*Emax/parm[0]);
	if (maxdx > parm[1])
		maxdx = parm[1];

	x = dx;

	while (x < maxdx)
		{
		tempweight = sq(Ein - 0.5*parm[0]*sq(x));

		//if (tempweight > minweight)
		//	break;
		
		twp = tempweight + 0.5*sq(dis + x - parm[1]);
		if (twp < minweight)
			{
			minweight = twp;
			proj_dis = parm[1] - x;
			}

		twm = tempweight + 0.5*sq(dis - x - parm[1]);
		if (twm < minweight)
			{
			minweight = twm;
			proj_dis = parm[1] + x;
			}

		x += dx;
		}

	double maxx, minx;

	minx = proj_dis - dx;
	if (minx < parm[1] - maxdx)
		minx = parm[1] - maxdx;

	maxx = proj_dis + dx;
	if (maxx > parm[1] + maxdx)
		maxx = parm[1] + maxdx;

	x = minx;

	while (x < maxx)
		{
		tempweight = sq(Ein - 0.5*parm[0]*sq(parm[1] - x)) + 0.5*sq(dis - x);

		if (tempweight < minweight)
			{
			minweight = tempweight;
			proj_dis = x;
			}

		x += dx/10.; 
		}

	new_dx = dx/10.;

	minx = proj_dis - new_dx;
	if (minx < parm[1] - maxdx)
		minx = parm[1] - maxdx;

	maxx = proj_dis + new_dx;
	if (maxx > parm[1] + maxdx)
		maxx = parm[1] + maxdx;

	x = minx;

	while (x < maxx)
		{
		tempweight = sq(Ein - 0.5*parm[0]*sq(parm[1] - x)) + 0.5*sq(dis - x);

		if (tempweight < minweight)
			{
			minweight = tempweight;
			proj_dis = x;
			}

		x += new_dx/10.; 
		}

	new_dx /= 10.;

	minx = proj_dis - new_dx;
	if (minx < parm[1] - maxdx)
		minx = parm[1] - maxdx;

	maxx = proj_dis + new_dx;
	if (maxx > parm[1] + maxdx)
		maxx = parm[1] + maxdx;

	x = minx;

	while (x < maxx)
		{
		tempweight = sq(Ein - 0.5*parm[0]*sq(parm[1] - x)) + 0.5*sq(dis - x);

		if (tempweight < minweight)
			{
			minweight = tempweight;
			proj_dis = x;
			}

		x += new_dx/10.; 
		}

	//proj_dis = parm[1];

	for(n=0;n<3;++n)
		{
		new_at[0][n] = at[0][n] + ((proj_dis-dis)/(2*dis))*(at[0][n]-at[1][n]);
		new_at[1][n] = at[1][n] - ((proj_dis-dis)/(2*dis))*(at[0][n]-at[1][n]);
		}

	double projE;
	projE = 0.5*parm[0]*sq(parm[1] - proj_dis);
	//projE = 0.;

	return projE;
	}
*/

double proj_dihed( double at[4][3], double ENRG, double new_at[4][3], double parm[3])
	{
	int m, n, c;
	double new_ENRG;
	double vecA[3], vecD[3], normal[3], perpCBA[3], perpDCB[3];
	double dotp, trpp, DHD, curweight, tempweight, arm1, arm2, armr;

	if (fabs(parm[0]) < 0.0001)
		{
		for(m=0;m<4;++m)
			for(n=0;n<3;++n)
				new_at[m][n] = at[m][n];

		return 0.;
		}
		

	for(n=0;n<3;++n)
		{
		vecA[n] = at[0][n] - at[1][n];
		vecD[n] = at[3][n] - at[2][n];
		normal[n] = at[1][n] - at[2][n];
		}

	normalize(normal);

	cross(vecD,normal,perpDCB);
	cross(vecA,normal,perpCBA);

	arm1 = norm(perpCBA);
	arm2 = norm(perpDCB);
	//printf("%lf %lf\n",arm1,arm2);

	if (arm1 < 0.01 && arm2 < 0.01) 
		{
		for(m=0;m<4;++m)
			for(n=0;n<3;++n)
				new_at[m][n] = at[m][n];

		return 0.;
		}

	dotp = dot_prod(perpCBA,perpDCB); 
	trpp = trip_prod(normal,perpDCB,perpCBA);
	
	DHD = atan2(trpp,dotp);

	//printf("%lf %lf %lf %lf %lf\n",DHD,parm[0],parm[1],parm[2],ENRG);
	

	c = (arm1 > arm2) ? 0 : 1;	
	armr = (c == 0) ? arm2/arm1 : arm1/arm2;

	//printf("%d %lf %lf %lf\n",c,arm1,arm2,armr);

	curweight = FLT_MAX;
	new_ENRG = ENRG;

	double dPhi, dPhi1, dPhi2;
	double Ep,Em;
	double pdp1,pdp2;

	dPhi1 = 0.;
	dPhi2 = 0.;
	pdp1 = 0.;
	pdp2 = 0.;
	dPhi = 0.;

	while (dPhi < M_PI)
		{
		if (c == 0)	
			tempweight = 2*sq(arm2*2.*sin(dPhi1/2.));
		else
			tempweight = 2*sq(arm1*2.*sin(dPhi1/2.));

		//printf("%lf %lf\n",dPhi1,tempweight);
		if (tempweight > curweight)
			break;

		Ep = ESCALE*(parm[0]*(1 + cos(parm[1]*(DHD+dPhi) - parm[2])));
		Em = ESCALE*(parm[0]*(1 + cos(parm[1]*(DHD-dPhi) - parm[2])));

		if (fabs(ENRG - Ep) < fabs(ENRG - Em))
			{
			tempweight += sq(ENRG - Ep); 
			//tempweight += 10.*sq(ENRG - Ep); 

			if (tempweight < curweight)
				{
				pdp1 = dPhi1;
				pdp2 = dPhi2;
				new_ENRG = Ep;
				curweight = tempweight;

				//printf("%lf %lf %lf %lf\n",new_ENRG,pdp1,pdp2,curweight);
				}
			}
	
		else
			{
			tempweight += sq(ENRG - Em); 
			//tempweight += 10.*sq(ENRG - Em); 

			if (tempweight < curweight)
				{
				pdp1 = -dPhi1;
				pdp2 = -dPhi2;
				new_ENRG = Em;
				curweight = tempweight;

				//printf("%lf %lf %lf %lf\n",new_ENRG,pdp1,pdp2,curweight);
				}

			}
	
		dPhi1 += M_PI/180.;
		dPhi2 = asin(armr*sin(dPhi1/2.));

		dPhi = dPhi1 + dPhi2;

		}

	if (c == 0)
		{
		dPhi2 = pdp1;
		dPhi1 = pdp2;
		}
	else
		{
		dPhi1 = pdp1;
		dPhi2 = pdp2;
		}
		
	
	RotationMatrix(normal, dPhi1);
	
	for(n=0;n<3;++n)
		new_at[0][n] = RM[n][0]*vecA[0] + RM[n][1]*vecA[1] + RM[n][2]*vecA[2] + at[1][n];

	RotationMatrix(normal, -dPhi2);
	
	for(n=0;n<3;++n)
		new_at[3][n] = RM[n][0]*vecD[0] + RM[n][1]*vecD[1] + RM[n][2]*vecD[2] + at[2][n];

	for(n=0;n<3;++n)
		{
		new_at[1][n] = at[1][n];
		new_at[2][n] = at[2][n];
		}

	double new_dhd;
	new_dhd = calc_dihed( new_at);
	new_ENRG = ESCALE*(parm[0]*(1 + cos(parm[1]*(new_dhd) - parm[2])));
	return new_ENRG;

/*
	double new_dhd;
	new_dhd = calc_dihed( new_at);

	printf("%lf %lf\n",new_dhd, DHD);
	Ep = parm[0]*(1 + cos(parm[1]*(DHD) - parm[2]));
	Em = parm[0]*(1 + cos(parm[1]*(new_dhd) - parm[2]));
	printf("%lf %lf %lf\n",parm[0],parm[1],parm[2]);
	printf("%lf %lf\n",Ep, Em);
	tempweight = RMSD(4,at,new_at);	
	printf("%lf\n",tempweight);

	for(int m=0;m<4;++m)
		{
		for(n=0;n<3;++n)
			printf("%lf ",at[m][n]);
		printf("\n");
		}

	printf("\n");

	for(int m=0;m<4;++m)
		{
		for(n=0;n<3;++n)
			printf("%lf ",new_at[m][n]);
		printf("\n");
		}
*/
	}

void ShortProj(int num_atoms, double (*out)[3], double (*in)[3], double (*EP)[3], double shortrad)
        {
        int m,n;
        double shortproj,AAdis;

	AAdis = 0;

        for(m=0;m<num_atoms;++m)
                for(n=0;n<3;++n)
                        AAdis += sq(EP[m][n] - in[m][n]);

	AAdis = sqrt(AAdis/num_atoms);

	shortproj = shortrad/AAdis;

	if (shortproj < 1.)
		{
		for(m=0;m<num_atoms;++m)
			for(n=0;n<3;++n)
				out[m][n] = (1.-shortproj)*EP[m][n] + shortproj*in[m][n];
		}
	else
		{
		for(m=0;m<num_atoms;++m)
			for(n=0;n<3;++n)
				out[m][n] = in[m][n];
		}
        }


void changeVar(double (****VBo)[2][3], double (***VSo)[3][2][2][3], double (**VOo)[3][3][2][3], double (*OSo)[3][3], double (***SCo)[3], double (**CTo)[3], double (**SSo)[4][3], double (*ADo)[8][3], double (**DHo)[4][3], double (*BLo)[2][3], double (**HBo)[4][3]) 
	{
	int i,j,m,n,k,c;

	for(j=0;j<NumBL;++j)
		for(i=0;i<2;++i)
			for(k=0; k<3; ++k)
				BLo[j][i][k] = atom[A2R[BLAtoms[j][i]][0]][A2R[BLAtoms[j][i]][1]][k];

	for(j=0;j<NumDH;++j)
		for(m=0;m<NumDHpar[j];++m)
			for(i=0;i<4;++i)
				for(k=0; k<3; ++k)
					DHo[j][m][i][k] = atom[A2R[DHAtoms[j][i]][0]][A2R[DHAtoms[j][i]][1]][k];

        for (i=0;i<NumCys;++i)
                for (j=0;j<NumCys;++j)
			for (n=0;n<3;++n)
				{
				SSo[i][j][0][n] = atom[CysIndex[i]][4][n];
				SSo[i][j][1][n] = atom[CysIndex[i]][5][n];
				SSo[i][j][2][n] = atom[CysIndex[j]][5][n];
				SSo[i][j][3][n] = atom[CysIndex[j]][4][n];
				}

        for (i=0;i<NumAD;++i)
		for (m=0;m<8;++m)
			for (n=0;n<3;++n)
				ADo[i][m][n] = atom[ADBK[i][m][0]][ADBK[i][m][1]][n];
		
        for (i=0;i<NumDN;++i)
                for (j=0;j<NumAC;++j)
                        for (n=0;n<3;++n)
                                {
                                HBo[i][j][0][n] = atom[HBBK[0][i][0]][0][n];
                                HBo[i][j][1][n] = atom[HBBK[0][i][0]][HBBK[0][i][1]][n];

                                HBo[i][j][2][n] = atom[HBBK[1][j][0]][HBBK[1][j][1]][n];
                                HBo[i][j][3][n] = atom[HBBK[1][j][0]][2][n];
                                }

        for (i=0;i<ProLen-1;++i)
		for (n=0;n<3;++n)
			{
			CTo[i][0][n] = atom[i][1][n];
			CTo[i][1][n] = atom[i][2][n];
			CTo[i][2][n] = atom[i][3][n];

			CTo[i][3][n] = atom[i+1][0][n];
			CTo[i][4][n] = atom[i+1][1][n];

			if (RT[i+1] != 11)
				CTo[i][5][n] = atom[i+1][lenAAH[RT[i+1]]-1][n];
			}


        for (i=0;i<ProLen;++i)
		for (j=0;j<NumSegsRP[i];++j)
			for (k=0;k<SegSizeRP[i][j];++k)
				for (n=0;n<3;++n)
					SCo[i][j][k][n] = atom[i][SegAtomRP[i][j][k]][n];

        for(i=0;i<ProLen;++i)
                for(j=i;j<ProLen;++j)
			for(m=0;m<lenRP[i];++m)
				for(n=0;n<lenRP[j];++n)
					for (c=0;c<2;++c)
						for (k=0;k<3;++k)
							{
							VBo[i][j][m][n][c][k] = atom[i][m][k];
							VBo[j][i][n][m][c][k] = atom[j][n][k];
							}

	for(i=0;i<ProLen;++i)
		for(m=0;m<lenRP[i];++m)
			for(j=0;j<NumSol;++j)
				for(n=0;n<3;++n)
					for(c=0;c<2;++c)
						for(k=0;k<3;++k)
							{
							VSo[i][m][j][n][c][0][k] = atom[i][m][k];
							VSo[i][m][j][n][c][1][k] = solatom[j][n][k];
							}

	for(i=0;i<NumSol;++i)
		for(j=0;j<NumSol;++j)
			for(m=0;m<3;++m)
				for(n=0;n<3;++n)
					for(c=0;c<2;++c)
						for(k=0;k<3;++k)
							VOo[i][j][m][n][c][k] = solatom[i][m][k];

	for(i=0;i<NumSol;++i)
		for(m=0;m<3;++m)
			for(k=0;k<3;++k)
				OSo[i][m][k] = solatom[i][m][k];
	} 

/*
double VEX(double disin, double sigma, double eps, double Emax)
	{
	if (disin > sigma)
		return disin;

	int p;

	for(p=0;p<501;p++)
		if (eps*LJP[p][1] < Emax)
			break;

	if (LJP[p][0]*sigma < disin)
		return disin;
	else
		return LJP[p][0]*sigma;
	}
*/

double VEX(double sigma, double eps, double Emax)
	{
	int p;

	for(p=0;p<501;p++)
		if (eps*LJP[p][1] < Emax)
			break;

	return LJP[p][0]*sigma;
	}

/*
double LJProj2(double Ein, double disin, double sigma, double eps, int projtype, int acc, double Emax)
	{
	int min_index,p,pin,pfi,p_min;
	double minweight,tempweight,dis,dis2;
	double r1,r2,r,dr,tempE,mind;

	dis = disin/sigma;

	minweight = FLT_MAX;

	if (Emax > 0.)
		{
		mind = (-2. + sqrt(4.*(1. + Emax/eps)))/(2.*Emax/eps);
		mind = pow(mind,1./6.);
		}
	else
		{
		printf("Error: Emax for lennard jones potential can't be negative");
		return disin/sigma;
		}

	if (mind > 0.5)
		p_min = (int) ((mind-0.5)*1000);
	else 
		p_min = 0;

	for(p=p_min;p<2501;p+=100)
		{
		tempweight = sq((Ein-ESCALE*eps*LJP[p][1]));

		if (projtype)
			tempweight += sq(sigma*(dis-LJP[p][0]));
		else
			tempweight += 0.5*sq(sigma*(dis-LJP[p][0]));

		if (tempweight < minweight)
			{
			minweight = tempweight;
			min_index = p;
			}
		}

	pin = (min_index == p_min) ? p_min : min_index-100;
	pfi = (min_index == 2500) ? 2500 : min_index+100;

	for(p=pin;p<pfi;p+=10)
		{
		tempweight = sq((Ein-ESCALE*eps*LJP[p][1]));

		if (projtype)
			tempweight += sq(sigma*(dis-LJP[p][0]));
		else
			tempweight += 0.5*sq(sigma*(dis-LJP[p][0]));

		if (tempweight < minweight)
			{
			minweight = tempweight;
			min_index = p;
			}
		}
	
	pin = (min_index == p_min) ? p_min : min_index-10;
	pfi = (min_index == 2500) ? 2500 : min_index+10;

	for(p=pin;p<pfi;p++)
		{
		tempweight = sq((Ein-ESCALE*eps*LJP[p][1]));

		if (projtype)
			tempweight += sq(sigma*(dis-LJP[p][0]));
		else
			tempweight += 0.5*sq(sigma*(dis-LJP[p][0]));

		if (tempweight < minweight)
			{
			minweight = tempweight;
			min_index = p;
			}
		}

	dis2 = LJP[min_index][0];

	if (dis > mind)
		{
		dis2 = dis;
		tempE = ESCALE*eps*(pow(dis2,-12.)-2.*pow(dis2,-6.)); 
		tempweight = sq(Ein - tempE);
		if (tempweight < minweight)
			{
			minweight = tempweight;
			dis2 = dis;
			}
		}

	else
		{
		tempweight = (projtype) ? sq((dis - mind)*sigma): 0.5*sq((dis - mind)*sigma);
		tempweight += sq(Ein - Emax);
		if (tempweight < minweight)
			{
			minweight = tempweight;
			dis2 = mind;
			}
		}

	if (Ein < Emax && Ein > 0)
		{
		if (Ein/eps < 0.00000001)
			r = pow(0.5,1./6.); 
		else	
			{
			r = (-2. + sqrt(4.*(1. + Ein/eps)))/(2.*Ein/eps);
			r = pow(r,1./6.);
			}

		tempweight = (projtype) ? sq((dis - r)*sigma): 0.5*sq((dis - r)*sigma);

		if (tempweight < minweight)
			{
			minweight = tempweight;
			dis2 = r;
			}
		} 

	else if (Ein > -eps)
		{
		if (Ein/eps > -0.00000001)
			{
			r1 = pow(0.5,1./6.); 
			r2 = 10.*pow(200,1./6.);; 
			}
		else	
			{
			r1 = (-2. + sqrt(4.*(1. + Ein/eps)))/(2.*Ein/eps);
			r1 = pow(r1,1./6.);

			r2 = (-2. - sqrt(4.*(1. + Ein/eps)))/(2.*Ein/eps);
			r2 = pow(r2,1./6.);
			}

		tempweight = (projtype) ? sq((dis - r1)*sigma): 0.5*sq((dis - r1)*sigma);

		if (tempweight < minweight)
			{
			minweight = tempweight;
			dis2 = r1;
			}

		tempweight = (projtype) ? sq((dis - r2)*sigma): 0.5*sq((dis - r2)*sigma);

		if (tempweight < minweight)
			{
			minweight = tempweight;
			dis2 = r2;
			}
		}

	return dis2;
	}

double LJProj(double Ein, double disin, double sigma, double eps, int projtype, int acc, double Emax)
	{
	int min_index,p,pin,pfi,p_min;
	double minweight,tempweight,dis,dis2;
	double r1,r2,r,dr,mind;

	dis = disin/sigma;

	if (dis > 3.0)
		{
		minweight = sq(Ein);
		min_index = -1;
		}

	else
		minweight = FLT_MAX;

	if (Emax > 0.)
		{
		mind = (-2. + sqrt(4.*(1. + Emax/eps)))/(2.*Emax/eps);
		mind = pow(mind,1./6.);
		}
	else
		{
		printf("Error: Emax for lennard jones potential can't be negative");
		return disin/sigma;
		}

	if (mind > 0.5)
		p_min = (int) ((mind-0.5)*1000);
	else 
		p_min = 0;

	for(p=p_min;p<2501;p+=100)
		{
		tempweight = sq((Ein-ESCALE*eps*LJP[p][1]));

		if (projtype)
			tempweight += sq(sigma*(dis-LJP[p][0]));
		else
			tempweight += 0.5*sq(sigma*(dis-LJP[p][0]));

		if (tempweight < minweight)
			{
			minweight = tempweight;
			min_index = p;
			}
		}

	if (min_index == -1)
		return dis;
	else
		{
		pin = (min_index == p_min) ? p_min : min_index-100;
		pfi = (min_index == 2500) ? 2500 : min_index+100;

		for(p=pin;p<pfi;p+=10)
			{
			tempweight = sq((Ein-ESCALE*eps*LJP[p][1]));

			if (projtype)
				tempweight += sq(sigma*(dis-LJP[p][0]));
			else
				tempweight += 0.5*sq(sigma*(dis-LJP[p][0]));

			if (tempweight < minweight)
				{
				minweight = tempweight;
				min_index = p;
				}
			}
		
		pin = (min_index == p_min) ? p_min : min_index-10;
		pfi = (min_index == 2500) ? 2500 : min_index+10;

		for(p=pin;p<pfi;p++)
			{
			tempweight = sq((Ein-ESCALE*eps*LJP[p][1]));

			if (projtype)
				tempweight += sq(sigma*(dis-LJP[p][0]));
			else
				tempweight += 0.5*sq(sigma*(dis-LJP[p][0]));

			if (tempweight < minweight)
				{
				minweight = tempweight;
				min_index = p;
				}
			}

		dis2 = LJP[min_index][0];

		r1 = (min_index == p_min) ? LJP[min_index][0] : LJP[min_index-1][0];
		r2 = (min_index == 2500) ? LJP[2500][0] : LJP[min_index+1][0];

		dr = (r2-r1)/20.;

		r = r1;
		while (r < r2)
			{
			tempweight = sq((Ein-ESCALE*eps*(pow(r,-12.)-2.*pow(r,-6.))));

			if (projtype)
				tempweight += sq(sigma*(dis-r));
			else
				tempweight += 0.5*sq(sigma*(dis-r));

			if (tempweight < minweight)
				{
				minweight = tempweight;
				dis2 = r;
				}
			
			r += dr;
			}

		r1 = (dis2 - dr < LJP[p_min][0]) ? LJP[p_min][0] : dis2 - dr;
		r2 = (dis2 + dr > LJP[2500][0]) ? LJP[2500][0] : dis2 + dr;

		dr = (r2-r1)/20.;

		r = r1;
		while (r < r2)
			{
			tempweight = sq((Ein-ESCALE*eps*(pow(r,-12.)-2.*pow(r,-6.))));

			if (projtype)
				tempweight += sq(sigma*(dis-r));
			else
				tempweight += 0.5*sq(sigma*(dis-r));

			if (tempweight < minweight)
				{
				minweight = tempweight;
				dis2 = r;
				}
			
			r += dr;
			}

		return dis2;

		}
	}
*/

double LJProjNew(double Ein, double disin, double sigma, double eps, int projtype, double Emax, int count)
	{
	double xc,yc,mind;
	double cutoff, tx, ms,delx, dely;
	double xi = disin/sigma;
	double yi = Ein/(eps*ESCALE);

	cutoff = pow(0.5,1./6.);

	for (int i = 0; i<count; ++i)
		{
		if (yi < 0.000000000001)
			{
			xc = xi;
			yc = (xi > cutoff) ? 0 : (pow(xi,-12.)-2.*pow(xi,-6.));
			}
		else
			{
			xc = (-1. + sqrt(1. + yi))/yi; 
			xc = pow(xc,1./6.);
			yc = (xi > cutoff) ? 0 : (pow(xi,-12.)-2.*pow(xi,-6.));

			if ((sigma*fabs(xc - xi)/sqrt(2)) > (eps*ESCALE*fabs(yc - yi)))
				xc = xi;
			else 
				yc = yi;

			if (xc > cutoff)
				{
				xi = xc;
				yi = 0; 
				break;
				}
			else
				{
				ms = -12.*(pow(xc,-13.) - pow(xc,-7.));
				tx = (yi - yc - ms*(xi - xc));
				delx = (ms*tx)/(0.5*sq(sigma/(eps*ESCALE))+sq(ms));
				dely = tx/(1. + 2.*sq(ms*eps*ESCALE/sigma));
				}

			xi += delx;
			yi -= dely;
			}
		}

	if (Emax > 0.)
		{
		mind = (-1. + sqrt(1. + Emax/eps))/(Emax/eps);
		mind = pow(mind,1./6.);

		if (xi < mind)
			xi = mind;
		}

	return xi;
	}

double CoulombProjTrue(double Ein, double disin, double rmin, double qmul, int projtype, int count, double epsDE)
	{
	double xc,yc,delx,dely,qeff;

	qeff = qmul/epsDE;

	if (fabs(qeff) < 0.00000001) 
		return disin;

	double xi = disin/sqrt(2);
	double yi = Ein;

	for (int i = 0; i<count; ++i)
		{
		xc = ESCALE*331.21*qeff/(yi*sqrt(2));
		yc = ESCALE*331.21*qeff/(xi*sqrt(2));
		
		if (xc < 0)
			xc = xi;
		else
			{
			if (fabs(xc - xi) > fabs(yc - yi))
				xc = xi;
			else
				yc = yi;
			}

		delx = (xc/(sq(xc) + sq(yc)))*((xi - xc)*xc - (yi - yc)*yc); 
		dely = (-yc/xc)*delx;

		xi = xc + delx;
		yi = yc + dely;
		}

	xi *= sqrt(2);

	if (xi < rmin)
		return rmin;

	return xi;
	}


void printeta(double avgerror, int itercount)
	{
	int i,j,m,n,c,z,p,Count = 0;

        for(i=0; i<ProLen; ++i)
		for(j=0; j<ProLen; ++j)
			for(m=0; m<lenRP[i]; ++m)
				for(n=0; n<lenRP[j]; ++n)
					for (c=0;c<2;++c)
						if (CONH[i][j][m][n] != 0)
							{
							val_index[Count][0] = Count;
							val_index[Count][1] = c;
							val_index[Count][2] = i;
							val_index[Count][3] = j;
							val_index[Count][4] = m;
							val_index[Count][5] = n;

							val[Count] = etaEB[i][j][m][n][c];
							Count++;
							}
/*
        for(i=0;i<ProLen;++i)
                for(m=0;m<lenRP[i];++m)
                        for(j=0;j<NumSol;++j)
                                for(n=0;n<3;++n)
                                        for(c=0;c<2;++c)
                                                etaES[i][m][j][n][c] += epsilon*((ESerr[i][m][j][n][c]/avgerr)-etaES[i][m][j][n][c]);

        for(i=0;i<NumSol;++i)
                for(j=i+1;j<NumSol;++j)
                        for(m=0;m<3;++m)
                                for(n=0;n<3;++n)
                                        for(c=0;c<2;++c)
						etaEO[i][j][m][n][c] += epsilon*((EOerr[i][j][m][n][c]/avgerr)-etaEO[i][j][m][n][c]);

        for (i=0;i<NumCys;++i)
                for (j=0;j<NumCys;++j)
			etaSS[i][j] += epsilon*((SSerr[i][j]/avgerr) - etaSS[i][j]);	

*/

        for(i=0; i<NumSol; ++i)
		{
		val_index[Count][0] = Count;
		val_index[Count][1] = 7;
		val_index[Count][2] = i;
		val[Count] = etaOS[i];
		Count++;
		}

        for (i=0;i<ProLen-1;++i)
		{
		val_index[Count][0] = Count;
		val_index[Count][1] = 8;
		val_index[Count][2] = i;
		val[Count] = etaCT[i];
		Count++;
		}

        for (i=0;i<ProLen;++i)
		for(j=0;j<NumSegsRP[i];++j)
			{
			val_index[Count][0] = Count;
			val_index[Count][1] = 9;
			val_index[Count][2] = i;
			val_index[Count][3] = j;
			val[Count] = etaSC[i][j];
			Count++;
			}

        //for (i=0;i<NumDN;++i)
        //        for (j=0;j<NumAC;++j)
	//		{
	//		val_index[Count][0] = Count;
	//		val_index[Count][1] = 10;
	//		val_index[Count][2] = i;
	//		val_index[Count][3] = j;
	//		val[Count] = etaAD[i][j];
	//		Count++;
	//		}

        for(i=0; i<NumDH; ++i)
		for(m=0;m<NumDHpar[i];++m)
			{
			val_index[Count][0] = Count;
			val_index[Count][1] = 11;
			val_index[Count][2] = i;
			val_index[Count][3] = m;
			val[Count] = etaDH[i][m];
			Count++;
			}

        for(i=0; i<NumBL; ++i)
		{
		val_index[Count][0] = Count;
		val_index[Count][1] = 12;
		val_index[Count][2] = i;
		val[Count] = etaBL[i];
		Count++;
		}

	qsort(val_index, Count, 9*sizeof(int), valinc);

        FILE *fp;

        fp=fopen(etafile,"a");

	fprintf(fp,"%lf %d \n",avgerror,itercount);
        for(p=0; p<20; ++p)
		{
		z = Count-1-p;
		if (val_index[z][1] == 0)
			{
			i = val_index[z][2];
			j = val_index[z][3];
			m = val_index[z][4];
			n = val_index[z][5];
			c = 0;
			fprintf(fp,"LJ \t %lf \t %lf %d %d %d %d\n",etaEB[i][j][m][n][c],EBerr[i][j][m][n][c],i,m,j,n);
			}

		else if (val_index[z][1] == 1)
			{
			i = val_index[z][2];
			j = val_index[z][3];
			m = val_index[z][4];
			n = val_index[z][5];
			c = 1;
			fprintf(fp,"CP \t %lf \t %lf %d %d %d %d\n",etaEB[i][j][m][n][c],EBerr[i][j][m][n][c],i,m,j,n);
			}

		else if (val_index[z][1] == 8)
			{
			i = val_index[z][2];
			fprintf(fp,"CT \t %lf \t %lf %d\n",etaCT[i],CTerr[i],i);
			}

		else if (val_index[z][1] == 9)
			{
			i = val_index[z][2];
			j = val_index[z][3];
			fprintf(fp,"SC \t %lf \t %lf %d %d\n",etaSC[i][j],SCerr[i][j],i,j);
			}

	//	else if (val_index[z][1] == 10)
	//		{
	//		i = val_index[z][2];
	//		j = val_index[z][3];
	//		fprintf(fp,"AD \t %lf \t %lf %d %d\n",etaAD[i][j],ADerr[i][j],i,j);
	//		}

		else if (val_index[z][1] == 11)
			{
			i = val_index[z][2];
			j = val_index[z][3];
			fprintf(fp,"DH \t %lf \t %lf %d %d \n",etaDH[i][j],DHerr[i][j],i,j);
			}

		else if (val_index[z][1] == 12)
			{
			i = val_index[z][2];
			fprintf(fp,"BL \t %lf \t %lf %d %d %d %d %d\n",etaBL[i],BLerr[i],i,A2R[BLAtoms[i][0]][0],A2R[BLAtoms[i][0]][1],A2R[BLAtoms[i][1]][0],A2R[BLAtoms[i][1]][1]);
			}
		}

        fprintf(fp,"\n");

        fclose(fp);
	}

/*
DEBUGGING
	for(j=0;j<NumBL;++j)
		{
		BEA[j] = BEo[j];
		
		for(p=0;p<2;++p)
			for(n=0;n<3;++n)
				BLA[j][p][n] = BLo[j][p][n];
		}

	for(j=0;j<NumDH;++j)
		for(m=0;m<NumDHpar[j];++m)
			{
			EDA[j][m] = EDo[j][m];

			for(p=0;p<4;++p)
				for(n=0;n<3;++n)
					DHA[j][m][p][n] = DHo[j][m][p][n];
			}

        for(i=0; i<ProLen; ++i)
		for(j=0; j<ProLen; ++j)
			for(m=0; m<lenRP[i]; ++m)
				for(n=0; n<lenRP[j]; ++n)
					{
					if (CONH[i][j][m][n] == 0)
						continue;

					for(k=0;k<3;++k)
						{
						VBA[i][j][m][n][1][k] = VBo[i][j][m][n][1][k]; 
						VBA[j][i][n][m][1][k] = VBo[j][i][n][m][1][k];
						}

					EBA[i][j][m][n][1] = EBo[i][j][m][n][1];
					}

	if (dis1 < VXRR[i][j][m][n])
		{						
		dis2 = VXRR[i][j][m][n];
		for(k=0;k<3;++k)
			{
		VBA[i][j][m][n][0][k] += ((dis2-dis1)/(2*dis1))*(VBo[i][j][m][n][0][k]-VBo[j][i][n][m][0][k]);
		VBA[j][i][n][m][0][k] -= ((dis2-dis1)/(2*dis1))*(VBo[i][j][m][n][0][k]-VBo[j][i][n][m][0][k]);
			}
		}

	EBA[i][j][m][n][0] = EBo[i][j][m][n][0];

        for(i=0; i<ProLen; ++i)
                for(m=0; m<lenRP[i]; ++m)
                        for(j=0; j<NumSol; ++j)
                                for(n=0; n<3; ++n)
                                        {
                                        for(k=0;k<3;++k)
                                                {
						VSA[i][m][j][n][1][0][k] = VSo[i][m][j][n][1][0][k];
						VSA[i][m][j][n][1][1][k] = VSo[i][m][j][n][1][1][k];
						}

					ESA[i][m][j][n][1] = ESo[i][m][j][n][1];
					}

        for(i=0; i<NumSol; ++i)
		for(j=i+1; j<NumSol; ++j)
			for(m=0; m<3; ++m)
                                for(n=0; n<3; ++n)
                                        {
                                        for(k=0;k<3;++k)
                                                {
						VOA[i][j][m][n][1][k] = VOo[i][j][m][n][1][k];
						VOA[j][i][n][m][1][k] = VOo[j][i][n][m][1][k];
						}

					EOA[i][j][m][n][1] = EOo[i][j][m][n][1];
					}
*/

void projA(double (****VBo)[2][3], double (***VSo)[3][2][2][3], double (**VOo)[3][3][2][3], double (****EBo)[2],  double (***ESo)[3][2], double (**EOo)[3][3][2], double (*OSo)[3][3], double (***SCo)[3], double (**CTo)[3], double (**SSo)[4][3], double (*ADo)[8][3], double (**DHo)[4][3], double **EDo, double (*BLo)[2][3], double *BEo, double (**HBo)[4][3]) 
	{
	int i,j,m,n,k,p,q,r,z,sizeCT;
	double dis1,dis2,qmul,rmin,eps, srad,tempweight,curweight;
	//int min_index;

	// Bond Length Constraint

	for(j=0;j<NumBL;++j)
		{
		proj_bondNE( BLo[j], BLA[j], BLpar[j]);
		BEA[j] = BEo[j];
		}

        // DiHed Constraint

	for(j=0;j<NumDH;++j)
		for(m=0;m<NumDHpar[j];++m)
			{
			EDA[j][m] = EDo[j][m];

			for(p=0;p<4;++p)
				for(n=0;n<3;++n)
					DHA[j][m][p][n] = DHo[j][m][p][n];
			}

	//for(j=0;j<NumDH;++j)
	//	for(m=0;m<NumDHpar[j];++m)
	//		EDA[j][m] = proj_dihed(DHo[j][m], EDo[j][m], DHA[j][m], DHpar[j][m]);

        // H2O Constraint
 
	double tempOS[3][3];

        for(i=0;i<NumSol;++i)
                {
                totiter = 20;
                motifproj(3,SolMotif,OSo[i],tempOS);
		tempweight = RMSD(3,OSo[i],tempOS);
                ShortProj(3, OSA[i], OSo[i], tempOS, 0.05);
                }

	//Disulphide Bond Constraint

//	double tempSS[4][3];
//	int CountSS = 0;
//
//        for (i=0;i<NumCys;++i)
//                for (j=0;j<NumCys;++j)
//                        {
//			if (i == j)
//				continue;
//
//			WeightsSS[i][j] = FLT_MAX;
//
//                        for(m=0;m<4;++m)
//                                for(n=0;n<3;++n)
//                                        SSA[i][j][m][n] = SSo[i][j][m][n];
//
//			for(z=0;z<SSMotifNum;++z)
//				{
//				totiter = 10;
//				motifproj(4,SSMotifs[z],SSo[i][j],tempSS);
//				tempweight = RMSD(4,SSo[i][j],tempSS);
//
//				if (tempweight < WeightsSS[i][j])
//					{
//					WeightsSS[i][j] = tempweight;
//
//					for(m=0;m<4;++m)
//						for(n=0;n<3;++n)
//							SSStore[i][j][m][n] = tempSS[m][n];
//					}
//				}
//
//			val_index[CountSS][0] = CountSS;
//			val_index[CountSS][1] = i;
//			val_index[CountSS][2] = j;
//
//			val[CountSS] = etaSS[i][j]*WeightsSS[i][j];
//			CountSS++;
//			}
//
//	qsort(val_index, CountSS, 9*sizeof(int), valinc);
//
//        CountSS = 0;
//
//        r = 0;
//
//        while (CountSS < 6 && r < NumCys*(NumCys-1))
//                {
//		i = val_index[r][0];
//                p = val_index[r][1];
//                q = val_index[r][2];
//
//		for(m=0;m<4;++m)
//			for(n=0;n<3;++n)
//				SSA[p][q][m][n] = SSStore[p][q][m][n];
//
//		CountSS++;
//                r++;
//                }

        // Hydrogen Bonding Constraint

        int CountHB = 0, resdiff;
        double tempHB[4][3];

        for (i=0;i<NumDN;++i)
                for (j=0;j<NumAC;++j)
                        for(m=0;m<4;++m)
                                for(n=0;n<3;++n)
                                        HBA[i][j][m][n] = HBo[i][j][m][n];

        // step 1: get weights

	int MINRSD = 3;
	//double BRAD1 = 0.

        for (i=0;i<NumDN;++i)
                for (j=0;j<NumAC;++j)
                        {
                        WeightsHB[i][j] = FLT_MAX;

                        resdiff = HBBK[0][i][0]-HBBK[1][j][0];

                        if (abs(resdiff) < MINRSD)
                                continue;

                        totiter = 10;
                        motifproj(4,ADMotif,HBo[i][j],tempHB);
                        tempweight = RMSD(4,HBo[i][j],tempHB);

                        if (tempweight < BRAD)
                                tempweight = 0.;
                        else
                                tempweight -= BRAD;

                        if (tempweight < WeightsHB[i][j])
                                {
                                WeightsHB[i][j] = tempweight;

                                for(m=0;m<4;++m)
                                        for(n=0;n<3;++n)
                                                HBStore[i][j][m][n] = tempHB[m][n];
                                }

                        val_index[CountHB][0] = CountHB;
                        val_index[CountHB][1] = i;
                        val_index[CountHB][2] = j;

                        val[CountHB] = etaHB[i][j]*WeightsHB[i][j];
                        CountHB++;
                        }

        qsort(val_index, CountHB, 9*sizeof(int), valinc);

        CountHB = 0;

        r = 0;

        int DNYN[NumDN],ACYN[NumAC];

        for (i=0;i<NumDN;i++)
                DNYN[i] = 0;
        for (i=0;i<NumAC;i++)
                ACYN[i] = 0;

        while (CountHB < HB1 && r < NumDN*NumAC)
                {
                i = val_index[r][0];
                p = val_index[r][1];
                q = val_index[r][2];

                //if ((DNYN[p] == 0) && (ACYN[q] == 0))
                //    {
                      ShortProj(4, HBA[p][q], HBo[p][q], HBStore[p][q], BRAD);
                      CountHB++;
                //    DNYN[p] = 1;
                //    ACYN[q] = 1;
                //    }
                r++;
                }







	int CountAD = 0;
	double tempAD[2][4][3], tempADo[2][4][3], totcost, shortAD[4][3];

        for (i=0;i<NumAD;++i)
		for(m=0;m<8;++m)
			for(n=0;n<3;++n)
				ADA[i][m][n] = ADo[i][m][n];

        // step 1: get weights

	totiter = 10;

        for (i=0;i<NumAD;++i)
		{
		WeightsAD[i] = FLT_MAX;

		for(m=0;m<4;++m)
			for(n=0;n<3;++n)
				{
				tempADo[0][m][n] = ADo[i][m][n];
				tempADo[1][m][n] = ADo[i][m+4][n];
				}

		totcost = 0.;

		motifproj(4,ADMotif,tempADo[0],tempAD[0]);
		motifproj(4,ADMotif,tempADo[1],tempAD[1]);

		tempweight = RMSD(4,tempADo[0],tempAD[0]);

		if (tempweight > BRAD)
			totcost += tempweight - BRAD;

		tempweight = RMSD(4,tempADo[1],tempAD[1]);

		if (tempweight > BRAD)
			totcost += tempweight - BRAD;

		if (totcost < WeightsAD[i])
			{
			WeightsAD[i] = totcost;

			ShortProj(4, shortAD, tempADo[0], tempAD[0], BRAD);

			for(m=0;m<4;++m)
				for(n=0;n<3;++n)
					ADStore[i][m][n] = shortAD[m][n];

			ShortProj(4, shortAD, tempADo[1], tempAD[1], BRAD);

			for(m=0;m<4;++m)
				for(n=0;n<3;++n)
					ADStore[i][m+4][n] = shortAD[m][n];
			}

		val_index[i][0] = i;
		val[i] = etaAD[i]*WeightsAD[i];
		}

	qsort(val_index, NumAD, 9*sizeof(int), valinc);

	int DApairs[ProLen][ProLen], p1,p2,q1,q2;

	for(i=0;i<ProLen;++i)
		for(j=0;j<ProLen;++j)
			DApairs[i][j] = 0;

	int DNs[ProLen], ACs[ProLen];

	for(i=0;i<ProLen;++i)
		{
		DNs[i] = 0;
		ACs[i] = 0;
		}

	CountAD = 0;
        r = 0;

	//HBB = 4; 

        while (CountAD < HB2 && r < NumAD)
                {
		i = val_index[r][0];

		p1 = ADBK[i][0][0];
		p2 = ADBK[i][4][0];
		q1 = ADBK[i][2][0];
		q2 = ADBK[i][6][0];

		//printf("%d %d %d %d\n",p1,p2,q1,q2);

		if (((DNs[p1] == 0 && ACs[q1] == 0) || DApairs[p1][q1] == 1) && ((DNs[p2] == 0 && ACs[q2] == 0) || DApairs[p2][q2] == 1) && (DApairs[p1][q1] != 1 || DApairs[p2][q2] != 1)) 
			{
			for(m=0;m<8;++m)
				for(n=0;n<3;++n)
					ADA[i][m][n] = ADStore[i][m][n];

			DApairs[p1][q1] = 1;
			DApairs[p2][q2] = 1;
			DNs[p1] = 1;
			DNs[p2] = 1;
			ACs[q1] = 1;
			ACs[q2] = 1;

			//printf("%d %d %d %d %lf\n",p1,p2,q1,q2,val[i]);

			CountAD++;
			}

                r++;
                }


	// CT Constraint

	double tempCT[6][3];

	for(i=0;i<ProLen-1;++i)
		{
		WeightsCT[i] = FLT_MAX;
		totiter = 10;
		sizeCT = (RT[i+1] != 11) ? 6 : 5;

		for(z=0;z<CTMotifNum[RT[i+1]];++z)
			{
			motifproj(sizeCT,CTMotifs[RT[i+1]][z],CTo[i],tempCT);
			tempweight = RMSD(sizeCT,CTo[i],tempCT);
			//printf("%lf ",tempweight);

			if (tempweight < WeightsCT[i])
				{
				WeightsCT[i] = tempweight; 
				for(j=0;j<sizeCT;++j)
					for(n=0;n<3;++n)
						CTA[i][j][n] = tempCT[j][n];
				}
			}

		for(j=0;j<sizeCT;++j)
			for(n=0;n<3;++n)
				tempCT[j][n] = CTA[i][j][n];

		ShortProj(sizeCT, CTA[i], CTo[i], tempCT, 0.05);

		//printf("%d %lf \n",i,WeightsCT[i]);

		//srad = 0.01;
		//shortproj = (srad < WeightsCT[i]) ? srad/WeightsCT[i] : 1.;

		
		//for(j=0;j<sizeCT;++j)
		//	for(n=0;n<3;++n)
		//		CTA[i][j][n] = (1.-shortproj)*CTA[i][j][n] + shortproj*CTo[i][j][n];
		}

	// SC Constraint
	
	double tempSC[16][3];

	srad = 0.05;
	
	for(i=0;i<ProLen;++i)
		for(j=0;j<NumSegs[RT[i]];++j)
			{
			curweight = FLT_MAX;

			//for(z=0;z<SCMotifNum[RT[i]][j];++z)
			for(z=0;z<1;++z)
				{
				totiter = 10;

				motifproj(SegSize[RT[i]][j],SCMotifs[RT[i]][j][z],SCo[i][j],tempSC);

				tempweight = RMSD(SegSize[RT[i]][j],SCo[i][j],tempSC);

				if (tempweight < curweight)
					{
					curweight  = tempweight;

					for(m=0;m<SegSize[RT[i]][j];++m)
						for(n=0;n<3;++n)
							SCStore[i][j][m][n] = tempSC[m][n];
							//SCA[i][j][m][n] = tempSC[m][n];
							//SCA[i][j][m][n] = SCo[i][j][m][n];
					}
				}

			ShortProj(SegSize[RT[i]][j], SCA[i][j], SCo[i][j], SCStore[i][j], srad);
			}

	motifproj(5,NTMotif,SCo[0][NumSegs[RT[0]]],SCStore[0][NumSegs[RT[0]]]);
	ShortProj(5, SCA[0][NumSegs[RT[0]]], SCo[0][NumSegs[RT[0]]], SCStore[0][NumSegs[RT[0]]], srad);

	motifproj(4,OXTMotif,SCo[ProLen-1][NumSegs[RT[ProLen-1]]],SCStore[ProLen-1][NumSegs[RT[ProLen-1]]]);
	ShortProj(4,SCA[ProLen-1][NumSegs[RT[ProLen-1]]],SCo[ProLen-1][NumSegs[RT[ProLen-1]]],SCStore[ProLen-1][NumSegs[RT[ProLen-1]]],srad);

        //Coulomb's Potential for Protein-Protein interaction

	int HBiter = 5;

        for(i=0; i<ProLen; ++i)
		for(j=0; j<ProLen; ++j)
			for(m=0; m<lenRP[i]; ++m)
				for(n=0; n<lenRP[j]; ++n)
					{
					if (CONH[i][j][m][n] == 0)
						continue;

					dis1 = 0.;
					for(k=0;k<3;++k)
						{
						VBA[i][j][m][n][1][k] = VBo[i][j][m][n][1][k];
						VBA[j][i][n][m][1][k] = VBo[j][i][n][m][1][k];
						dis1 += sq(VBo[i][j][m][n][1][k] - VBo[j][i][n][m][1][k]);
						}

					//EBA[i][j][m][n][1] = EBo[i][j][m][n][1];

					dis1 = sqrt(dis1);
					
					qmul = INTRP[i][j][m][n][0];
					rmin = 0.5*INTRP[i][j][m][n][2];

					dis2 =  CoulombProjTrue( EBo[i][j][m][n][1], dis1, rmin, qmul, 0, HBiter,epsRR);

					for(k=0;k<3;++k)
						{
						VBA[i][j][m][n][1][k] += ((dis2-dis1)/(2*dis1))*(VBo[i][j][m][n][1][k]-VBo[j][i][n][m][1][k]);
						VBA[j][i][n][m][1][k] -= ((dis2-dis1)/(2*dis1))*(VBo[i][j][m][n][1][k]-VBo[j][i][n][m][1][k]);
						}

					EBA[i][j][m][n][1] = ESCALE*331.21*qmul/(epsRR*dis2); 
					}

        //LJ Potential for Protein-Protein interaction

        for(i=0; i<ProLen; ++i)
		for(j=0; j<ProLen; ++j)
			for(m=0; m<lenRP[i]; ++m)
				for(n=0; n<lenRP[j]; ++n)
					{
					if (CONH[i][j][m][n] == 0)
						continue;

					dis1 = 0.;
					for(k=0;k<3;++k)
						{
						VBA[i][j][m][n][0][k] = VBo[i][j][m][n][0][k];
						VBA[j][i][n][m][0][k] = VBo[j][i][n][m][0][k];
						dis1 += sq(VBo[i][j][m][n][0][k] - VBo[j][i][n][m][0][k]);
						}

					dis1 = sqrt(dis1);

					eps = INTRP[i][j][m][n][1];
					rmin = INTRP[i][j][m][n][2];

					//dis2 = LJProj(EBo[i][j][m][n][0], dis1, rmin, eps, 0, 1, 4.0);
					dis2 = LJProjNew(EBo[i][j][m][n][0], dis1, rmin, eps, 1, 4.0,4);
					EBA[i][j][m][n][0] = ESCALE*eps*(pow(dis2,-12.)-2.*pow(dis2,-6.));

					if (EBA[i][j][m][n][0] < 0.)
						EBA[i][j][m][n][0] = 0.;

					dis2 *= rmin;

					for(k=0;k<3;++k)
						{
						VBA[i][j][m][n][0][k] += ((dis2-dis1)/(2*dis1))*(VBo[i][j][m][n][0][k]-VBo[j][i][n][m][0][k]);
						VBA[j][i][n][m][0][k] -= ((dis2-dis1)/(2*dis1))*(VBo[i][j][m][n][0][k]-VBo[j][i][n][m][0][k]);
						}
					}

//        for(i=0; i<ProLen; ++i)
//		for(j=0; j<ProLen; ++j)
//			for(m=0; m<lenRP[i]; ++m)
//				for(n=0; n<lenRP[j]; ++n)
//					{
//					if (CONH[i][j][m][n] == 0)
//						continue;
//
//					dis1 = 0.;
//					for(k=0;k<3;++k)
//						{
//						VBA[i][j][m][n][0][k] = VBo[i][j][m][n][0][k];
//						VBA[j][i][n][m][0][k] = VBo[j][i][n][m][0][k];
//						dis1 += sq(VBo[i][j][m][n][0][k] - VBo[j][i][n][m][0][k]);
//						}
//
//					dis1 = sqrt(dis1);
//
//					eps = INTRP[i][j][m][n][1];
//					rmin = INTRP[i][j][m][n][2];
//
//					dis1 /= rmin;
//					if (dis1 < VXRR[i][j][m][n])
//						{						
//						dis2 = VXRR[i][j][m][n];
//						for(k=0;k<3;++k)
//							{
//							VBA[i][j][m][n][0][k] += ((dis2-dis1)/(2*dis1))*(VBo[i][j][m][n][0][k]-VBo[j][i][n][m][0][k]);
//							VBA[j][i][n][m][0][k] -= ((dis2-dis1)/(2*dis1))*(VBo[i][j][m][n][0][k]-VBo[j][i][n][m][0][k]);
//							}
//						}
//
//					EBA[i][j][m][n][0] = EBo[i][j][m][n][0];
//					}

        //Coulomb's Potential for Solvent-Protein interaction

        for(i=0; i<ProLen; ++i)
                for(m=0; m<lenRP[i]; ++m)
                        for(j=0; j<NumSol; ++j)
                                for(n=0; n<3; ++n)
                                        {
                                        dis1 = 0.;
                                        for(k=0;k<3;++k)
                                                {
						VSA[i][m][j][n][1][0][k] = VSo[i][m][j][n][1][0][k];
						VSA[i][m][j][n][1][1][k] = VSo[i][m][j][n][1][1][k];
						dis1 += sq(VSo[i][m][j][n][1][0][k] - VSo[i][m][j][n][1][1][k]);
						}
			
					dis1 = sqrt(dis1);

					qmul = parRP[i][m][0]*PQS[n];
					rmin = 0.5*(parRP[i][m][2] + parSol[n][2]);

					HBiter = (dis1 > 10.) ? 2 : 3;
					dis2 =  CoulombProjTrue( ESo[i][m][j][n][1], dis1, rmin, qmul, 0, HBiter,epsRS);

					for(k=0;k<3;++k)
						{
						VSA[i][m][j][n][1][0][k] += ((dis2-dis1)/(2*dis1))*(VSo[i][m][j][n][1][0][k] - VSo[i][m][j][n][1][1][k]);
						VSA[i][m][j][n][1][1][k] -= ((dis2-dis1)/(2*dis1))*(VSo[i][m][j][n][1][0][k] - VSo[i][m][j][n][1][1][k]);
						}

					ESA[i][m][j][n][1] = ESCALE*331.21*qmul/(epsRS*dis2); 
					}

        //Lennard Jones Potential for Solvent-Protein interaction


        for(i=0; i<ProLen; ++i)
                for(m=0; m<lenRP[i]; ++m)
                        for(j=0; j<NumSol; ++j)
                                for(n=0; n<3; ++n)
                                        {
                                        dis1 = 0.;
                                        for(k=0;k<3;++k)
                                                {
						VSA[i][m][j][n][0][0][k] = VSo[i][m][j][n][0][0][k];
						VSA[i][m][j][n][0][1][k] = VSo[i][m][j][n][0][1][k];
						dis1 += sq(VSo[i][m][j][n][0][0][k] - VSo[i][m][j][n][0][1][k]);
						}
			
					dis1 = sqrt(dis1);

					ESA[i][m][j][n][0] = ESo[i][m][j][n][0];

					//if (dis1 < VXRS[i][m][n] && n != 0)
					if (dis1 < VXRS[i][m][n])
						{						
						dis2 = VXRS[i][m][n];
						for(k=0;k<3;++k)
							{
							VSA[i][m][j][n][0][0][k] += ((dis2-dis1)/(2*dis1))*(VSo[i][m][j][n][0][0][k] - VSo[i][m][j][n][0][1][k]);
							VSA[i][m][j][n][0][1][k] -= ((dis2-dis1)/(2*dis1))*(VSo[i][m][j][n][0][0][k] - VSo[i][m][j][n][0][1][k]);
							}
						}


//					if (n == 0)
//						{
//						rmin = (parRP[i][m][2] + parSol[n][2]);
//						eps = sqrt(parRP[i][m][1]*parSol[n][1]);
//
//						//dis2 = LJProj(ESo[i][m][j][n][0], dis1, rmin, eps, 0, 1, 4.);
//						dis2 = LJProjNew(ESo[i][m][j][n][0], dis1, rmin, eps, 1, 4.0, 4);
//
//						ESA[i][m][j][n][0] = ESCALE*eps*(pow(dis2,-12.)-2.*pow(dis2,-6.));
//						dis2 *= rmin;
//
//						for(k=0;k<3;++k)
//							{
//							VSA[i][m][j][n][0][0][k] += ((dis2-dis1)/(2*dis1))*(VSo[i][m][j][n][0][0][k] - VSo[i][m][j][n][0][1][k]);
//							VSA[i][m][j][n][0][1][k] -= ((dis2-dis1)/(2*dis1))*(VSo[i][m][j][n][0][0][k] - VSo[i][m][j][n][0][1][k]);
//							}
//						}
					}

        //Coulomb's Potential for Solvent-Solvent interaction

        for(i=0; i<NumSol; ++i)
		for(j=i+1; j<NumSol; ++j)
			for(m=0; m<3; ++m)
                                for(n=0; n<3; ++n)
                                        {
                                        for(k=0;k<3;++k)
                                                {
						VOA[i][j][m][n][1][k] = VOo[i][j][m][n][1][k];
						VOA[j][i][n][m][1][k] = VOo[j][i][n][m][1][k];
						}

					EOA[i][j][m][n][1] = EOo[i][j][m][n][1];
					}

//        for(i=0; i<NumSol; ++i)
//		for(j=i+1; j<NumSol; ++j)
//			for(m=0; m<3; ++m)
//                                for(n=0; n<3; ++n)
//                                        {
//                                        dis1 = 0.;
//                                        for(k=0;k<3;++k)
//                                                {
//						VOA[i][j][m][n][1][k] = VOo[i][j][m][n][1][k];
//						VOA[j][i][n][m][1][k] = VOo[j][i][n][m][1][k];
//						dis1 += sq(VOo[i][j][m][n][1][k] - VOo[j][i][n][m][1][k]);
//						}
//			
//					dis1 = sqrt(dis1);
//
//					qmul = PQS[m]*PQS[n];
//					rmin = 0.5*(parSol[m][2] + parSol[n][2]);
//
//					HBiter = (dis1 > 10.) ? 4 : 5;
//					dis2 = CoulombProjTrue(EOo[i][j][m][n][1], dis1, rmin, qmul, 0, HBiter, epsSS);
//
//					for(k=0;k<3;++k)
//						{
//						VOA[i][j][m][n][1][k] += ((dis2-dis1)/(2.*dis1))*(VOo[i][j][m][n][1][k] - VOo[j][i][n][m][1][k]);
//						VOA[j][i][n][m][1][k] -= ((dis2-dis1)/(2.*dis1))*(VOo[i][j][m][n][1][k] - VOo[j][i][n][m][1][k]);
//						}
//
//					EOA[i][j][m][n][1] = ESCALE*331.21*qmul/(dis2*epsSS); 
//					}
//
        //LJ Potential for Solvent-Solvent interaction

        for(i=0; i<NumSol; ++i)
		for(j=i+1; j<NumSol; ++j)
			for(m=0; m<3; ++m)
                                for(n=0; n<3; ++n)
                                        {
                                        dis1 = 0.;
                                        for(k=0;k<3;++k)
                                                {
						VOA[i][j][m][n][0][k] = VOo[i][j][m][n][0][k];
						VOA[j][i][n][m][0][k] = VOo[j][i][n][m][0][k];
						dis1 += sq(VOo[i][j][m][n][0][k] - VOo[j][i][n][m][0][k]);
						}
			
					dis1 = sqrt(dis1);

					EOA[i][j][m][n][0] = EOo[i][j][m][n][0];

					//if (dis1 < VXSS[m][n] && (n != 0 || m != 0))
					if (dis1 < VXSS[m][n])
						{						
						dis2 = VXSS[m][n];
						for(k=0;k<3;++k)
							{
							VOA[i][j][m][n][0][k] += ((dis2-dis1)/(2.*dis1))*(VOo[i][j][m][n][0][k] - VOo[j][i][n][m][0][k]);
							VOA[j][i][n][m][0][k] -= ((dis2-dis1)/(2.*dis1))*(VOo[i][j][m][n][0][k] - VOo[j][i][n][m][0][k]);
							}
						}

//					if (n == 0 && m == 0)
//						{
//						rmin = (parSol[m][2] + parSol[n][2]);
//						eps = sqrt(parSol[m][1]*parSol[n][1]);
//
//						min_index = LJProj(EOo[i][j][m][n][0], dis1, rmin, eps, 0, 0, 5.);
//
//						if (min_index == -1)
//							{
//							EOA[i][j][m][n][0] = 0.;
//							dis2 = dis1;
//							}
//						else
//							{
//							EOA[i][j][m][n][0] = ESCALE*eps*LJP[min_index][1];
//							dis2 = rmin*LJP[min_index][0];
//							}
//
//						for(k=0;k<3;++k)
//							{
//							VOA[i][j][m][n][0][k] += ((dis2-dis1)/(2.*dis1))*(VOo[i][j][m][n][0][k] - VOo[j][i][n][m][0][k]);
//							VOA[j][i][n][m][0][k] -= ((dis2-dis1)/(2.*dis1))*(VOo[i][j][m][n][0][k] - VOo[j][i][n][m][0][k]);
//							}
//						}
					}
	}


void reflect(double (****VBo)[2][3], double (***VSo)[3][2][2][3], double (**VOo)[3][3][2][3], double (****EBo)[2],  double (***ESo)[3][2], double (**EOo)[3][3][2], double (*OSo)[3][3], double (***SCo)[3], double (**CTo)[3], double (**SSo)[4][3], double (*ADo)[8][3], double (**DHo)[4][3], double **EDo, double (*BLo)[2][3], double *BEo, double (**HBo)[4][3]) 
	{
	int i,j,m,n,c,k,p,sizeCT;

        for(i=0;i<ProLen;++i)
                for(j=0;j<ProLen;++j)
			for(m=0;m<lenRP[i];++m)
				for(n=0;n<lenRP[j];++n)
					for(c=0;c<2;++c)
						{
						if (CONH[i][j][m][n] == 0)
							continue;

						EBR[i][j][m][n][c] = 2.*EBo[i][j][m][n][c] - EB[i][j][m][n][c];

						for(k=0;k<3;++k)
							{
							VBR[i][j][m][n][c][k] = 2.*VBo[i][j][m][n][c][k] - VB[i][j][m][n][c][k];
							VBR[j][i][n][m][c][k] = 2.*VBo[j][i][n][m][c][k] - VB[j][i][n][m][c][k];
							}
						}

        for(i=0;i<ProLen;++i)
                for(m=0;m<lenRP[i];++m)
                        for(j=0;j<NumSol;++j)
                                for(n=0;n<3;++n)
                                        for(c=0;c<2;++c)
						{
						ESR[i][m][j][n][c] = 2.*ESo[i][m][j][n][c] - ES[i][m][j][n][c];

                                                for(k=0;k<3;++k)
							{
                                                        VSR[i][m][j][n][c][0][k] = 2.*VSo[i][m][j][n][c][0][k] - VS[i][m][j][n][c][0][k];
                                                        VSR[i][m][j][n][c][1][k] = 2.*VSo[i][m][j][n][c][1][k] - VS[i][m][j][n][c][1][k];
							}
						}

        for(i=0;i<NumSol;++i)
                for(j=i+1;j<NumSol;++j)
                        for(m=0;m<3;++m)
                                for(n=0;n<3;++n)
                                        for(c=0;c<2;++c)
						{
						EOR[i][j][m][n][c] = 2.*EOo[i][j][m][n][c] - EO[i][j][m][n][c];

						for (k=0;k<3;++k)
							{
							VOR[i][j][m][n][c][k] = 2.*VOo[i][j][m][n][c][k] - VO[i][j][m][n][c][k];
							VOR[j][i][n][m][c][k] = 2.*VOo[j][i][n][m][c][k] - VO[j][i][n][m][c][k];
							}
						}

        for (i=0;i<NumSol;++i)
                for(p=0;p<3;++p)
                        for(n=0;n<3;++n)
                                OSR[i][p][n] = 2.*OSo[i][p][n] - OS[i][p][n];

/*
        for (i=0;i<NumCys;++i)
                for (j=0;j<NumCys;++j)
                        for(p=0;p<4;++p)
                                for(n=0;n<3;++n)
					SSR[i][j][p][n] = 2.*SSo[i][j][p][n] - SS[i][j][p][n];
*/

        for (i=0;i<ProLen-1;++i)
		{
		sizeCT = (RT[i+1] != 11) ? 6 : 5;
		for(p=0;p<sizeCT;++p)
			for(n=0;n<3;++n)
				CTR[i][p][n] = 2.*CTo[i][p][n] - CT[i][p][n];
		}

	for(i=0;i<ProLen;++i)
		for(j=0;j<NumSegsRP[i];++j)
			for(m=0;m<SegSizeRP[i][j];++m)
				for(n=0;n<3;++n)
					SCR[i][j][m][n] = 2.*SCo[i][j][m][n] - SC[i][j][m][n]; 

        for (i=0;i<NumAD;++i)
		for(m=0;m<8;++m)
			for(n=0;n<3;++n)
				ADR[i][m][n] = 2.*ADo[i][m][n] - AD[i][m][n];

        for (i=0;i<NumDN;++i)
                for (j=0;j<NumAC;++j)
                        for(p=0;p<4;++p)
                                for(n=0;n<3;++n)
                                        HBR[i][j][p][n] = 2.*HBo[i][j][p][n] - HB[i][j][p][n];

	for(j=0;j<NumDH;++j)
		for(m=0;m<NumDHpar[j];++m)
			{
			EDR[j][m] = 2.*EDo[j][m] - ED[j][m];

			for(p=0;p<4;++p)
				for(n=0;n<3;++n)
					DHR[j][m][p][n] = 2.*DHo[j][m][p][n] - DH[j][m][p][n];
			}

	for(j=0;j<NumBL;++j)
		{
		BER[j] = 2.*BEo[j] - BE[j];

		for(p=0;p<2;++p)
			for(n=0;n<3;++n)
				BLR[j][p][n] = 2.*BLo[j][p][n] - BL[j][p][n];
		}
	}


void projB(double (****VBo)[2][3], double (***VSo)[3][2][2][3], double (**VOo)[3][3][2][3], double (****EBo)[2],  double (***ESo)[3][2], double (**EOo)[3][3][2], double (*OSo)[3][3], double (***SCo)[3], double (**CTo)[3], double (**SSo)[4][3], double (*ADo)[8][3], double (**DHo)[4][3], double **EDo, double (*BLo)[2][3], double *BEo, double (**HBo)[4][3]) 
	{
	int i,j,k,m,n,c,atom_type;
	double NumVar[ProLen][26],NumVarSol[NumSol][3],Etot,lambda,COM[3];

        for (i=0;i<ProLen;++i)
                for (j=0;j<lenRP[i];++j)
			{
			NumVar[i][j] = 0.;

                        for (n=0;n<3;++n)
                                atom[i][j][n] = 0.;
			}

        for (i=0;i<NumSol;++i)
                for (j=0;j<3;++j)
                        {
                        NumVarSol[i][j] = 0.;

                        for (n=0;n<3;++n)
                                solatom[i][j][n] = 0.;
                        }

	Etot = 0.;
	lambda = 0.;

	for(j=0;j<NumBL;++j)
		{
		//Etot += BEo[j];
		//lambda += 1./sq(etaBL[j]);

		BEB[j] = BEo[j];

		for(i=0;i<2;++i)
			{
			NumVar[A2R[BLAtoms[j][i]][0]][A2R[BLAtoms[j][i]][1]] += sq(etaBL[j]);

			for(k=0; k<3; ++k)
				atom[A2R[BLAtoms[j][i]][0]][A2R[BLAtoms[j][i]][1]][k] += sq(etaBL[j])*BLo[j][i][k];
			}
		}

	for(j=0;j<NumDH;++j)
		for(m=0;m<NumDHpar[j];++m)
			{
			//Etot += EDo[j][m];
			//lambda += 1./sq(etaDH[j][m]);

			EDB[j][m] = EDo[j][m];

			for(i=0;i<4;++i)
				{
				NumVar[A2R[DHAtoms[j][i]][0]][A2R[DHAtoms[j][i]][1]] += sq(etaDH[j][m]);

				for(k=0; k<3; ++k)
					atom[A2R[DHAtoms[j][i]][0]][A2R[DHAtoms[j][i]][1]][k] += sq(etaDH[j][m])*DHo[j][m][i][k];
				}
			}

        for (i=0;i<NumAD;++i)
		for (m=0;m<8;++m)
			{
			NumVar[ADBK[i][m][0]][ADBK[i][m][1]] += sq(etaAD[i]);

			for (n=0;n<3;++n)
				atom[ADBK[i][m][0]][ADBK[i][m][1]][n] += sq(etaAD[i])*ADo[i][m][n];
			}

        for (i=0;i<NumDN;++i)
                for (j=0;j<NumAC;++j)
                        {
                        NumVar[HBBK[0][i][0]][0] += sq(etaHB[i][j]);
                        NumVar[HBBK[0][i][0]][HBBK[0][i][1]] += sq(etaHB[i][j]);
                        NumVar[HBBK[1][j][0]][HBBK[1][j][1]] += sq(etaHB[i][j]);
                        NumVar[HBBK[1][j][0]][2] += sq(etaHB[i][j]);

                        for (n=0;n<3;++n)
                                {
                                atom[HBBK[0][i][0]][0][n] += sq(etaHB[i][j])*HBo[i][j][0][n];
                                atom[HBBK[0][i][0]][HBBK[0][i][1]][n] += sq(etaHB[i][j])*HBo[i][j][1][n];
                                atom[HBBK[1][j][0]][HBBK[1][j][1]][n] += sq(etaHB[i][j])*HBo[i][j][2][n];
                                atom[HBBK[1][j][0]][2][n] += sq(etaHB[i][j])*HBo[i][j][3][n];
                                }
                        }

        for(i=0;i<ProLen;++i)
                for(j=0;j<ProLen;++j)
			for(m=0;m<lenRP[i];++m)
				for(n=0;n<lenRP[j];++n)
					for(c=0;c<2;++c)
						{
						if (CONH[i][j][m][n] == 0)
							continue;	

						//if (c == 1)
						//	{
							Etot += EBo[i][j][m][n][c];
							lambda += 1./sq(etaEB[i][j][m][n][c]);
						//	}

						EBB[i][j][m][n][c] = EBo[i][j][m][n][c];

						NumVar[i][m] += sq(etaEB[i][j][m][n][c]);
						NumVar[j][n] += sq(etaEB[i][j][m][n][c]);

						for (k=0;k<3;++k)
							{
							atom[i][m][k] += sq(etaEB[i][j][m][n][c])*VBo[i][j][m][n][c][k];
							atom[j][n][k] += sq(etaEB[i][j][m][n][c])*VBo[j][i][n][m][c][k];
							}
						}

        for(i=0;i<ProLen;++i)
                for(m=0;m<lenRP[i];++m)
                        for(j=0;j<NumSol;++j)
                                for(n=0;n<3;++n)
                                        for(c=0;c<2;++c)
						{
						ESB[i][m][j][n][c] = ESo[i][m][j][n][c];

						NumVarSol[j][n] += sq(etaES[i][m][j][n][c]);
						NumVar[i][m] += sq(etaES[i][m][j][n][c]);

						for (k=0;k<3;++k)
							{
							atom[i][m][k] += sq(etaES[i][m][j][n][c])*VSo[i][m][j][n][c][0][k];
							solatom[j][n][k] += sq(etaES[i][m][j][n][c])*VSo[i][m][j][n][c][1][k];
							}

						//if (c == 0 && n != 0)
						if (c == 0)
							continue;
						else
							{
							Etot += ESo[i][m][j][n][c];
							lambda += 1./sq(etaES[i][m][j][n][c]);
							}
						}

        for(i=0;i<NumSol;++i)
                for(j=i+1;j<NumSol;++j)
                        for(m=0;m<3;++m)
                                for(n=0;n<3;++n)
                                        for(c=0;c<2;++c)
						{
						EOB[i][j][m][n][c] = EOo[i][j][m][n][c];

						NumVarSol[i][m] += sq(etaEO[i][j][m][n][c]);
						NumVarSol[j][n] += sq(etaEO[i][j][m][n][c]);

						for (k=0;k<3;++k)
							{
							solatom[i][m][k] += sq(etaEO[i][j][m][n][c])*VOo[i][j][m][n][c][k];
							solatom[j][n][k] += sq(etaEO[i][j][m][n][c])*VOo[j][i][n][m][c][k];
							}
/*

						if (c == 0 && ( n!= 0 || m != 0))
							continue;
						else
							{
							Etot += EOo[i][j][m][n][c];
							lambda += 1./sq(etaEO[i][j][m][n][c]);
							}
*/
						}

	EBT = Etot;

        for(i=0;i<NumSol;++i)
                for(m=0;m<3;++m)
                        {
                        NumVarSol[i][m] += sq(etaOS[i]);

                        for (k=0;k<3;++k)
                                solatom[i][m][k] += sq(etaOS[i])*OSo[i][m][k];
                        }
/*
        for (i=0;i<NumCys;++i)
                for (j=0;j<NumCys;++j)
			{
			if (i==j)
				continue;

			NumVar[CysIndex[i]][4] += sq(etaSS[i][j]);
			NumVar[CysIndex[i]][5] += sq(etaSS[i][j]);
			NumVar[CysIndex[j]][5] += sq(etaSS[i][j]);
			NumVar[CysIndex[j]][4] += sq(etaSS[i][j]);

			for (n=0;n<3;++n)
				{
				atom[CysIndex[i]][4][n] += sq(etaSS[i][j])*SSo[i][j][0][n];
				atom[CysIndex[i]][5][n] += sq(etaSS[i][j])*SSo[i][j][1][n];
				atom[CysIndex[j]][5][n] += sq(etaSS[i][j])*SSo[i][j][2][n];
				atom[CysIndex[j]][4][n] += sq(etaSS[i][j])*SSo[i][j][3][n];
				}
			}
*/

        for (i=0;i<ProLen-1;++i)
		{
		NumVar[i][1] += sq(etaCT[i]);
		NumVar[i][2] += sq(etaCT[i]);
		NumVar[i][3] += sq(etaCT[i]);
		NumVar[i+1][0] += sq(etaCT[i]);
		NumVar[i+1][1] += sq(etaCT[i]);

		for (n=0;n<3;++n)
			{
			atom[i][1][n] += sq(etaCT[i])*CTo[i][0][n];
			atom[i][2][n] += sq(etaCT[i])*CTo[i][1][n];
			atom[i][3][n] += sq(etaCT[i])*CTo[i][2][n];
			atom[i+1][0][n] += sq(etaCT[i])*CTo[i][3][n];
			atom[i+1][1][n] += sq(etaCT[i])*CTo[i][4][n];
			}

		if (RT[i+1] != 11)
			{
			NumVar[i+1][lenAAH[RT[i+1]]-1] += sq(etaCT[i]);

			for (n=0;n<3;++n)
				atom[i+1][lenAAH[RT[i+1]]-1][n] += sq(etaCT[i])*CTo[i][5][n];
			}
		}

	for(i=0;i<ProLen;++i)
		for(j=0;j<NumSegsRP[i];++j)
			for(m=0;m<SegSizeRP[i][j];++m)
				{
				atom_type = SegAtomRP[i][j][m];
				NumVar[i][atom_type] += sq(etaSC[i][j]); 

				for(n=0;n<3;++n)
					atom[i][atom_type][n] += sq(etaSC[i][j])*SCo[i][j][m][n]; 
				}

	int tot_atoms = 0;

	for (n=0;n<3;++n)
		COM[n] = 0.;

        for (i=0;i<ProLen;++i)
		for(m=0; m<lenRP[i]; ++m)
			{
			tot_atoms++;
                        for (n=0;n<3;++n)
				{
                                atom[i][m][n] /= NumVar[i][m];
				COM[n] += atom[i][m][n];
				}
			}

        for (i=0;i<NumSol;++i)
                for(m=0; m<3; ++m)
                        for (n=0;n<3;++n)
				{
                                solatom[i][m][n] /= NumVarSol[i][m];
				COM[n] += solatom[i][m][n];
				}

	tot_atoms += NumSol*3;

	for (n=0;n<3;++n)
		COM[n] /= tot_atoms;

	// Total Energy Constraint

	if (Etot > EGS*ESCALE)
		{
		for(i=0;i<ProLen;++i)
			for(j=0;j<ProLen;++j)
				for(m=0;m<lenRP[i];++m)
					for(n=0;n<lenRP[j];++n)
						for(c=0;c<2;++c) /////// 
							{
							if (CONH[i][j][m][n] == 0)
								continue;	

							EBB[i][j][m][n][c] = EBo[i][j][m][n][c] + ((EGS*ESCALE-Etot)/(lambda*sq(etaEB[i][j][m][n][c])));
							}

		for(i=0;i<ProLen;++i)
			for(m=0;m<lenRP[i];++m)
				for(j=0;j<NumSol;++j)
					for(n=0;n<3;++n)
						for(c=0;c<2;++c) /////// 
							{
							//if (c == 0 && n != 0)
							if (c == 0)
								continue;
							else
								ESB[i][m][j][n][c] = ESo[i][m][j][n][c] + ((EGS*ESCALE-Etot)/(lambda*sq(etaES[i][m][j][n][c])));
							}
		} 

/*
		for(j=0;j<NumBL;++j)
			BEB[j] = BEo[j] + ((EGS*ESCALE-Etot)/(lambda*sq(etaBL[j])));

		for(j=0;j<NumDH;++j)
			for(m=0;m<NumDHpar[j];++m)
				EDB[j][m] = EDo[j][m] + ((EGS*ESCALE-Etot)/(lambda*sq(etaDH[j][m])));

		for(i=0;i<NumSol;++i)
			for(j=i+1;j<NumSol;++j)
				for(m=0;m<3;++m)
					for(n=0;n<3;++n)
						for(c=0;c<2;++c) /////// 
							{
							if (c == 0 && ( n!= 0 || m != 0))
								continue;
							else
								EOB[i][j][m][n][c] = EOo[i][j][m][n][c] + ((EGS*ESCALE-Etot)/(lambda*sq(etaEO[i][j][m][n][c])));
							}
*/

	//if (EGS < -300 && iter > 5000 && (iter%100 == 0))
	//	EGS += 2.;

	//inertia(atom, solatom, temp);

	double dis;

	for(j=0;j<NumSol;++j)
		for(n=0;n<3;++n)
			{
			dis = 0.;
                        for (k=0;k<3;++k)
				dis += sq(solatom[j][n][k] - COM[k]);

			dis = sqrt(dis);
			
			if (dis > RAD)
				for (k=0;k<3;++k)
					solatom[j][n][k] = (RAD/dis)*(solatom[j][n][k] - COM[k]) + COM[k];
			}

	for(i=0;i<ProLen;++i)
		for(m=0;m<lenRP[i];++m)
			{
			dis = 0.;
                        for (k=0;k<3;++k)
				dis += sq(atom[i][m][k] - COM[k]);

			dis = sqrt(dis);
			
			if (dis > RAD)
				for (k=0;k<3;++k)
					atom[i][m][k] = (RAD/dis)*(atom[i][m][k] - COM[k]) + COM[k];
			}

	changeVar(VBB, VSB, VOB, OSB, SCB, CTB, SSB, ADB, DHB, BLB, HBB);
	}


void RRR()
	{
	int i,j,k,m,n,c,p;
	double diff,avgerr,totCons,totVar;

	projA(VB, VS, VO, EB, ES, EO, OS, SC, CT, SS, AD, DH, ED, BL, BE, HB);
	reflect(VBA, VSA, VOA, EBA, ESA, EOA, OSA, SCA, CTA, SSA, ADA, DHA, EDA, BLA, BEA, HBA);
	projB(VBR, VSR, VOR, EBR, ESR, EOR, OSR, SCR, CTR, SSR, ADR, DHR, EDR, BLR, BER, HBR);

	tEBerr[0] = 0.;
	tEBerr[1] = 0.;
	tESerr[0] = 0.;
	tESerr[1] = 0.;
	tEOerr[0] = 0.;
	tEOerr[1] = 0.;
	tOSerr = 0.;
	tSCerr = 0.;
	tCTerr = 0.;
	tSSerr = 0.;
	tADerr = 0.;
	tDHerr = 0.;
	tBLerr = 0.;
	tHBerr = 0.;

        totVar = 0.;
        totCons = 0.;
        toterr = 0.;
        avgerr = 0.;

	for(j=0;j<NumBL;++j)
		{
		BLerr[j] = 0.; 

		diff = BEB[j] - BEA[j];
		BE[j] += beta*diff;
		BLerr[j] += sq(diff);

		for(p=0;p<2;++p)
			for(n=0;n<3;++n)
				{
				diff = BLB[j][p][n] - BLA[j][p][n];
				BL[j][p][n] += beta*diff;
				BLerr[j] += sq(diff);
				}

		tBLerr += BLerr[j];
		BLerr[j] /= 7.;
		avgerr += BLerr[j];
		}

        toterr += tBLerr;
        tBLerr /= NumBL*7.;
	totVar += NumBL*7; 
	totCons += NumBL;

	int totDHPar = 0;

	for(j=0;j<NumDH;++j)
		for(m=0;m<NumDHpar[j];++m)
			{
			totDHPar++;

			DHerr[j][m] = 0.;

			diff = EDB[j][m] - EDA[j][m];
			ED[j][m] += beta*diff;
			DHerr[j][m] += sq(diff);

			for(p=0;p<4;++p)
				for(n=0;n<3;++n)
					{
					diff = DHB[j][m][p][n] - DHA[j][m][p][n];
					DH[j][m][p][n] += beta*diff;
					DHerr[j][m] += sq(diff);
					}

			tDHerr += DHerr[j][m];
			DHerr[j][m] /= 13.;
			avgerr += DHerr[j][m];
			}

        toterr += tDHerr;
        tDHerr /= totDHPar*13.;
	totVar += totDHPar*13; 
	totCons += totDHPar;

	// Hydrogen Bonding Constraint

        for (i=0;i<NumDN;++i)
                for (j=0;j<NumAC;++j)
                        {
                        HBerr[i][j] = 0.;

                        for(p=0;p<4;++p)
                                for(n=0;n<3;++n)
                                        {
                                        diff = HBB[i][j][p][n] - HBA[i][j][p][n];
                                        HB[i][j][p][n] += beta*diff;
                                        HBerr[i][j] += sq(diff);
                                        }
                        tHBerr += HBerr[i][j];
                        HBerr[i][j] /= 12.;
                        avgerr += HBerr[i][j];
                        }

        toterr += tHBerr;
        tHBerr /= (NumDN*NumAC*12.);
        totVar += (NumDN*NumAC*12.);
        totCons += (NumDN*NumAC);

	// Secondary Structure Constraint
        for (i=0;i<NumAD;++i)
		{
		ADerr[i] = 0.;

		for(p=0;p<8;++p)
			for(n=0;n<3;++n)
				{
				diff = ADB[i][p][n] - ADA[i][p][n];
				AD[i][p][n] += beta*diff;
				ADerr[i] += sq(diff);
				}

		tADerr += ADerr[i];
		ADerr[i] /= 24.;
		avgerr += ADerr[i];
		}

        toterr += tADerr;
        tADerr /= (NumAD*24.);
	totVar += (NumAD*24.); 
	totCons += NumAD;


	int NumVBVar = 0;

	for(i=0;i<ProLen;++i)
		for(j=0;j<ProLen;++j)
			for(m=0;m<lenRP[i];++m)
				for(n=0;n<lenRP[j];++n)
					{
					if (CONH[i][j][m][n] == 0)
						continue;	
					
					NumVBVar++;

					for(c=0;c<2;++c)
						{
						EBerr[i][j][m][n][c] = 0.;

						diff = EBB[i][j][m][n][c] - EBA[i][j][m][n][c];
						EB[i][j][m][n][c] += beta*diff;
						EBerr[i][j][m][n][c] += sq(diff);

						for (k=0;k<3;++k)
							{
							diff = VBB[i][j][m][n][c][k] - VBA[i][j][m][n][c][k];
							VB[i][j][m][n][c][k] += beta*diff;
							EBerr[i][j][m][n][c] += sq(diff);

							diff = VBB[j][i][n][m][c][k] - VBA[j][i][n][m][c][k];
							VB[j][i][n][m][c][k] += beta*diff;
							EBerr[i][j][m][n][c] += sq(diff);
							}

						tEBerr[c] += EBerr[i][j][m][n][c];

						EBerr[i][j][m][n][c] /= 7.;
						avgerr += EBerr[i][j][m][n][c];
						}
					}

        toterr += tEBerr[0];
        tEBerr[0] /= (NumVBVar*7.);

        toterr += tEBerr[1];
        tEBerr[1] /= (NumVBVar*7.);

	totCons += NumVBVar*2.;
	totVar += NumVBVar*2.*7.;


        int NumVSVar = 0;

        for(i=0;i<ProLen;++i)
                for(m=0;m<lenRP[i];++m)
                        for(j=0;j<NumSol;++j)
                                for(n=0;n<3;++n)
					{
					NumVSVar += 1;

                                        for(c=0;c<2;++c)
                                                {
                                                ESerr[i][m][j][n][c] = 0.;

						diff = ESB[i][m][j][n][c] - ESA[i][m][j][n][c];
						ES[i][m][j][n][c] += beta*diff;
						ESerr[i][m][j][n][c] += sq(diff);
	
                                                for (k=0;k<3;++k)
                                                        {
                                                        diff = VSB[i][m][j][n][c][0][k] - VSA[i][m][j][n][c][0][k];
                                                        VS[i][m][j][n][c][0][k] += beta*diff;
							ESerr[i][m][j][n][c] += sq(diff);

                                                        diff = VSB[i][m][j][n][c][1][k] - VSA[i][m][j][n][c][1][k];
                                                        VS[i][m][j][n][c][1][k] += beta*diff;
							ESerr[i][m][j][n][c] += sq(diff);
							}

                                                tESerr[c] += ESerr[i][m][j][n][c];
                                                ESerr[i][m][j][n][c] /= 7.;
                                                avgerr += ESerr[i][m][j][n][c];
                                                }
					}

	if (NumSol >= 1)
		{
		toterr += tESerr[0];
		tESerr[0] /= (NumVSVar*7.);

		toterr += tESerr[1];
		tESerr[1] /= (NumVSVar*7.);

		totCons += NumVSVar*2.;
		totVar += NumVSVar*2.*7.;
		}
	else
		{
		tESerr[0] = 0.;
		tESerr[1] = 0.;
		}

        for(i=0;i<NumSol;++i)
                for(j=i+1;j<NumSol;++j)
                        for(m=0;m<3;++m)
                                for(n=0;n<3;++n)
                                        for(c=0;c<2;++c)
						{
						EOerr[i][j][m][n][c] = 0.;

						diff = EOB[i][j][m][n][c] - EOA[i][j][m][n][c];
						EO[i][j][m][n][c] += beta*diff;
						EOerr[i][j][m][n][c] += sq(diff);

						for (k=0;k<3;++k)
							{
							diff = VOB[i][j][m][n][c][k] - VOA[i][j][m][n][c][k];
							VO[i][j][m][n][c][k] += beta*diff;
							EOerr[i][j][m][n][c] += sq(diff);

							diff = VOB[j][i][n][m][c][k] - VOA[j][i][n][m][c][k];
							VO[j][i][n][m][c][k] += beta*diff;
							EOerr[i][j][m][n][c] += sq(diff);
							}

						tEOerr[c] += EOerr[i][j][m][n][c];
						EOerr[i][j][m][n][c] /= 7.;
						avgerr += EOerr[i][j][m][n][c];
						}

	if (NumSol >= 2)
		{
		toterr += tEOerr[0];
		tEOerr[0] /= (0.5*NumSol*(NumSol-1)*9*7);

		toterr += tEOerr[1];
		tEOerr[1] /= (0.5*NumSol*(NumSol-1)*9*7);

		totCons += NumSol*(NumSol-1)*9.;
		totVar += NumSol*(NumSol-1)*9*7;
		}
	else
		{
		tEOerr[0] = 0.;
		tEOerr[1] = 0.;
		}


        for (i=0;i<NumSol;++i)
                {
                OSerr[i] = 0.; 
                for (j=0;j<3;++j)
                        for (k=0;k<3;++k)
                                {
                                diff = OSB[i][j][k] - OSA[i][j][k];
                                OS[i][j][k] += beta*diff;
                                OSerr[i] += sq(diff);
                                }
                tOSerr += OSerr[i];
                OSerr[i] /= 9.;
                avgerr += OSerr[i];
                }

	if (NumSol >= 1)
		{
		toterr += tOSerr;
		tOSerr /= (NumSol*9.);

		totCons += NumSol;
		totVar += NumSol*9;
		}
	else
		tOSerr = 0.;
/*
        for (i=0;i<NumCys;++i)
                for (j=0;j<NumCys;++j)
			{
			if (i==j)
				continue;

			SSerr[i][j] = 0.;

			for(p=0;p<4;++p)
				for(n=0;n<3;++n)
					{
					diff = SSB[i][j][p][n] - SSA[i][j][p][n];
					SS[i][j][p][n] += beta*diff;
					SSerr[i][j] += sq(diff);
					}

			tSSerr += SSerr[i][j];
			SSerr[i][j] /= 12.;
			avgerr += SSerr[i][j];
			}

        toterr += tSSerr;
        tSSerr /= (NumCys*(NumCys-1)*12.);
	totVar += (NumCys*(NumCys-1)*12); 
	totCons += (NumCys*(NumCys-1));
*/
	int NumCTVar = 0,sizeCT;

        for (i=0;i<ProLen-1;++i)
		{
		CTerr[i] = 0.; 

		sizeCT = (RT[i+1] != 11) ? 6 : 5;
		NumCTVar += sizeCT;

		for (j=0;j<sizeCT;++j)
			for (k=0;k<3;++k)
				{
				diff = CTB[i][j][k] - CTA[i][j][k];
				CT[i][j][k] += beta*diff;
				CTerr[i] += sq(diff);
				}

		tCTerr += CTerr[i];
		CTerr[i] /= sizeCT*3.;
		avgerr += CTerr[i];
		}

        toterr += tCTerr;
        tCTerr /= (NumCTVar*3.);

	totVar += (NumCTVar*3);
	totCons += (ProLen-1);

	int NumSCVar = 0;

        for (i=0;i<ProLen;++i)
		{
		totCons += NumSegsRP[i];
		for(j=0;j<NumSegsRP[i];++j)
			{
			SCerr[i][j] = 0.;
			for(m=0;m<SegSizeRP[i][j];++m)
				for(n=0;n<3;++n)
					{
					diff = SCB[i][j][m][n] - SCA[i][j][m][n];
					SC[i][j][m][n] += beta*diff;
					SCerr[i][j] += sq(diff);
					}

			tSCerr += SCerr[i][j];
			SCerr[i][j] /= (SegSizeRP[i][j]*3.);
			NumSCVar += SegSizeRP[i][j];
			avgerr += SCerr[i][j];
			}
		}

        toterr += tSCerr;
        tSCerr /= (NumSCVar*3.);

	totVar += (NumSCVar*3);

	////////
        toterr /= totVar;
        avgerr /= totCons;
	////////

	//if ((iter%100==0) ) 
	//	printeta(avgerr,iter);

        for(i=0; i<ProLen; ++i)
		for(j=0; j<ProLen; ++j)
			for(m=0; m<lenRP[i]; ++m)
				for(n=0; n<lenRP[j]; ++n)
					for (c=0;c<2;++c)
						if (CONH[i][j][m][n] != 0)
							etaEB[i][j][m][n][c] += epsilon*((EBerr[i][j][m][n][c]/avgerr)-etaEB[i][j][m][n][c]);

        for(i=0;i<ProLen;++i)
                for(m=0;m<lenRP[i];++m)
                        for(j=0;j<NumSol;++j)
                                for(n=0;n<3;++n)
                                        for(c=0;c<2;++c)
                                                etaES[i][m][j][n][c] += epsilon*((ESerr[i][m][j][n][c]/avgerr)-etaES[i][m][j][n][c]);

        for(i=0;i<NumSol;++i)
                for(j=i+1;j<NumSol;++j)
                        for(m=0;m<3;++m)
                                for(n=0;n<3;++n)
                                        for(c=0;c<2;++c)
						etaEO[i][j][m][n][c] += epsilon*((EOerr[i][j][m][n][c]/avgerr)-etaEO[i][j][m][n][c]);

        for(i=0; i<NumSol; ++i)
                etaOS[i] += epsilon*((OSerr[i]/avgerr)-etaOS[i]);

        for (i=0;i<NumCys;++i)
                for (j=0;j<NumCys;++j)
			etaSS[i][j] += epsilon*((SSerr[i][j]/avgerr) - etaSS[i][j]);	

        for (i=0;i<ProLen-1;++i)
		etaCT[i] += epsilon*((CTerr[i]/avgerr)-etaCT[i]);

        for (i=0;i<ProLen;++i)
		for(j=0;j<NumSegsRP[i];++j)
			etaSC[i][j] += epsilon*((SCerr[i][j]/avgerr)-etaSC[i][j]);

        for (i=0;i<NumAD;++i)
		etaAD[i] += epsilon*((ADerr[i]/avgerr) - etaAD[i]);

        for (i=0;i<NumDN;++i)
                for (j=0;j<NumAC;++j)
                        etaHB[i][j] += epsilon*((HBerr[i][j]/avgerr) - etaHB[i][j]);

        for(i=0; i<NumDH; ++i)
		for(m=0;m<NumDHpar[i];++m)
			etaDH[i][m] += epsilon*((DHerr[i][m]/avgerr)-etaDH[i][m]);

        for(i=0; i<NumBL; ++i)
		etaBL[i] += epsilon*((BLerr[i]/avgerr)-etaBL[i]);
	

        tEBerr[0] = sqrt(tEBerr[0]);
        tEBerr[1] = sqrt(tEBerr[1]);
        tESerr[0] = sqrt(tESerr[0]);
        tESerr[1] = sqrt(tESerr[1]);
        tEOerr[0] = sqrt(tEOerr[0]);
        tEOerr[1] = sqrt(tEOerr[1]);
        tOSerr = sqrt(tOSerr);
        tSSerr = sqrt(tSSerr);
        tSCerr = sqrt(tSCerr);
        tCTerr = sqrt(tCTerr);
        tADerr = sqrt(tADerr);
        tDHerr = sqrt(tDHerr);
        tBLerr = sqrt(tBLerr);
        tHBerr = sqrt(tHBerr);

        toterr = sqrt(toterr);
	}


void printsol(char *solfile)
        {
        FILE *fp;
        int i,j,k;

        fp=fopen(solfile,"a");

        for(i=0;i<5;++i)
                for(j=0;j<ProLen;++j)
                        {    
                        for(k=0;k<3;++k)
                                fprintf(fp,"%lf ",atom[j][i][k]);
                        fprintf(fp,"\n");
                        }   

	for(j=0;j<3;++j)
		for(i=0;i<NumSol;++i)
			{    
			for(k=0;k<3;++k)
				fprintf(fp,"%lf ",solatom[i][j][k]);
			fprintf(fp,"\n");
			}    

        fprintf(fp,"\n");

        fclose(fp);
        }

void printinit(char *initfile, double Eout)
        {
        FILE *fp;
        int i,j,k,n;

        fp=fopen(initfile,"w");

	fprintf(fp,"%d %d\n",ProLen,NumSol);
	fprintf(fp,"%lf\n\n",Eout);

	for(i=0;i<ProLen;++i)
		{
		for (j=0;j<lenRP[i];++j)
			{
                        for(k=0;k<3;++k)
                                fprintf(fp,"%lf ",atom[i][j][k]);
                        fprintf(fp,"\n");
			}
		fprintf(fp,"\n");
		}

	for(i=0;i<NumSol;++i)
		{
		for(n=0;n<3;++n)
			{
			for(k=0;k<3;++k)
                                fprintf(fp,"%lf ",solatom[i][n][k]);
                        fprintf(fp,"\n");
			}
		fprintf(fp,"\n");
		}

        fclose(fp);
        }

int solve(int maxiter,int iterstride,double stoperr)
        {
        FILE *fperr,*fpen;
        //FILE *fpsol;
        int i,j,m,n,k;
	double dis,rmin,eps,qmul,LJRS,CPRS,LJSS,CPSS,LJRR,CPRR,TE,tmpE, ETX;  
	double tmp[4][3],dhd,Edhd, Ebond;

	EMIN = FLT_MAX;
        for(iter=0;iter<=maxiter;++iter)
                {
                RRR();

		//if(iter%500==0)
		//	printf("%.2e\t%.2e %.2e %.2e %.2e %.2e %.2e %.2e %.2e %.2e %.2e %.2e %.2e %.2e\n",toterr,tBBerr,tOSerr,tCTerr,tSCerr,tSSerr,tEBerr[0],tEBerr[1],tESerr[0],tESerr[1],tEOerr[0],tEOerr[1],tDHerr,tBLerr);

                if((iter%iterstride==0) || (toterr < stoperr))
                        {
			LJRR = 0.; 
			CPRR = 0.; 
			LJRS = 0.; 
			CPRS = 0.; 
			LJSS = 0.;
			CPSS = 0.;
			ETX = 0.;

			for(i=0; i<ProLen; ++i)
				for(j=i; j<ProLen; ++j)
					for(m=0; m<lenRP[i]; ++m)
						for(n=0; n<lenRP[j]; ++n)
							if (CONH[i][j][m][n] != 0)
								{
								dis = 0.;
								for(k=0; k<3; ++k)
									dis += sq(atom[j][n][k] - atom[i][m][k]);
								dis = sqrt(dis);

								qmul = INTRP[i][j][m][n][0];
								eps = INTRP[i][j][m][n][1];
								rmin = INTRP[i][j][m][n][2];

								tmpE = eps*(pow((rmin/dis), 12) - 2.*pow((rmin/dis), 6));
								if (tmpE < 0.)
									tmpE = 0.;

								LJRR += tmpE;
								ETX += tmpE;

								tmpE = 331.21*qmul/(epsRR*dis);
								CPRR += tmpE;
								ETX += tmpE;
								}

			for (i=0;i<ProLen;++i)
				for(m=0;m<lenRP[i];++m)
					for (j=0;j<NumSol;++j)
						for(n=0; n<3; ++n)
							{
							dis = 0.;
							for(k=0; k<3; ++k)
								dis += sq(solatom[j][n][k] - atom[i][m][k]);
							dis = sqrt(dis);

							rmin = (parRP[i][m][2] + parSol[n][2]);
							eps = sqrt(parRP[i][m][1]*parSol[n][1]);

							tmpE = eps*(pow((rmin/dis), 12) - 2.*pow((rmin/dis), 6));
							//if (tmpE < 0.)
							//	tmpE = 0.;

							LJRS += tmpE;

							//if (n == 0)
							//	ETX += tmpE;

							tmpE = 331.21*parRP[i][m][0]*parSol[n][0]/(epsRS*dis);
							CPRS += tmpE;
							ETX += tmpE;
							}

			for(i=0;i<NumSol;++i)
				for(j=i+1;j<NumSol;++j)
					for(m=0;m<3;++m)
						for(n=0;n<3;++n)
							{
							dis = 0.;
							for(k=0; k<3; ++k)
								dis += sq(solatom[j][n][k] - solatom[i][m][k]);
							dis = sqrt(dis);
							rmin = (parSol[m][2] + parSol[n][2]);
							eps = sqrt(parSol[m][1]*parSol[n][1]);

							tmpE =  eps*(pow((rmin/dis), 12) - 2.*pow((rmin/dis), 6));
							LJSS += tmpE;
							
							//if (m == 0 && n == 0)
							//	ETX += tmpE;

							tmpE = 331.21*parSol[m][0]*parSol[n][0]/(epsSS*dis);
							CPSS += tmpE;
							//ETX += tmpE;
							}

			TE = CPRR+LJRR+CPRS+LJRS+CPSS+LJSS;

			Edhd = 0.;
			for(j=0;j<NumDH;++j)
				{
				for(i=0;i<4;++i)
					for(k=0; k<3; ++k)
						tmp[i][k] = atom[A2R[DHAtoms[j][i]][0]][A2R[DHAtoms[j][i]][1]][k];

				dhd = calc_dihed(tmp);
			
				for(m=0;m<NumDHpar[j];++m)
					Edhd += DHpar[j][m][0]*(1 + cos(DHpar[j][m][1]*dhd - DHpar[j][m][2]));
				}

			TE += Edhd;

			Ebond = 0.;

			for(j=0;j<NumBL;++j)
				{
				dis = 0.;

                                for(k=0;k<3;++k)
					dis += sq(atom[A2R[BLAtoms[j][0]][0]][A2R[BLAtoms[j][0]][1]][k] - atom[A2R[BLAtoms[j][1]][0]][A2R[BLAtoms[j][1]][1]][k]);
				dis = sqrt(dis);

				Ebond += 0.5*BLpar[j][0]*sq(dis - BLpar[j][1]);
				}

			TE += Ebond;
			//ETX += Ebond;
			//ETX += Edhd;

			//ETX -= LJRR;

			if (TE < EMIN && toterr < 0.005)
				{
				EMIN = TE;
				printinit(initfile,ETX);
				}


			fperr=fopen(errfile,"a");

                        fprintf(fperr,"%.2e\t%.2e %.2e %.2e %.2e %.2e %.2e %.2e %.2e %.2e %.2e %.2e %.2e %.2e %.2e\n",toterr,tADerr,tHBerr,tOSerr,tCTerr,tSCerr,tSSerr,tEBerr[0],tEBerr[1],tESerr[0],tESerr[1],tEOerr[0],tEOerr[1],tDHerr,tBLerr);

                        fclose(fperr);

			fpen = fopen(Enerfile,"a");

                        fprintf(fpen,"%.1f %.1f %.1f %.1f %.1f \t  %.1f %.1f %.1f %.1f %.1f %.1f\n",ETX, EBT/ESCALE, TE, Edhd, Ebond, CPRR, CPRS, CPSS, LJRR, LJRS, LJSS);

                        fclose(fpen);
                        printsol(solfile);

			
			if(toterr<stoperr)
				return iter;
                        }
                }

        return 0;
        }


void makevars()
        {
        int i,j,k,m,sizeCT;

	AD = malloc((NumAD)*sizeof(double*[8][3]));
	ADA = malloc((NumAD)*sizeof(double*[8][3]));
	ADR = malloc((NumAD)*sizeof(double*[8][3]));
	ADB = malloc((NumAD)*sizeof(double*[8][3]));

        HB = malloc((NumDN)*sizeof(double*));
        HBA = malloc((NumDN)*sizeof(double*));
        HBR = malloc((NumDN)*sizeof(double*));
        HBB = malloc((NumDN)*sizeof(double*));

        for(i=0;i<NumDN;i++)
                {
                HB[i] = malloc((NumAC)*sizeof(double*[4][3]));
                HBA[i] = malloc((NumAC)*sizeof(double*[4][3]));
                HBR[i] = malloc((NumAC)*sizeof(double*[4][3]));
                HBB[i] = malloc((NumAC)*sizeof(double*[4][3]));
                }

        BL = malloc(NumBL*sizeof(double*[2][3]));
        BLA = malloc(NumBL*sizeof(double*[2][3]));
        BLR = malloc(NumBL*sizeof(double*[2][3]));
        BLB = malloc(NumBL*sizeof(double*[2][3]));

        BE = malloc(NumBL*sizeof(double));
        BEA = malloc(NumBL*sizeof(double));
        BER = malloc(NumBL*sizeof(double));
        BEB = malloc(NumBL*sizeof(double));

        DH = malloc(NumDH*sizeof(double*));
        DHA = malloc(NumDH*sizeof(double*));
        DHR = malloc(NumDH*sizeof(double*));
        DHB = malloc(NumDH*sizeof(double*));

	for(i=0;i<NumDH;i++)
		{
		DH[i] = malloc(NumDHpar[i]*sizeof(double*[4][3]));
		DHA[i] = malloc(NumDHpar[i]*sizeof(double*[4][3]));
		DHR[i] = malloc(NumDHpar[i]*sizeof(double*[4][3]));
		DHB[i] = malloc(NumDHpar[i]*sizeof(double*[4][3]));
		}

        ED = malloc(NumDH*sizeof(double*));
        EDA = malloc(NumDH*sizeof(double*));
        EDR = malloc(NumDH*sizeof(double*));
        EDB = malloc(NumDH*sizeof(double*));

	for(i=0;i<NumDH;i++)
		{
		ED[i] = malloc(NumDHpar[i]*sizeof(double));
		EDA[i] = malloc(NumDHpar[i]*sizeof(double));
		EDR[i] = malloc(NumDHpar[i]*sizeof(double));
		EDB[i] = malloc(NumDHpar[i]*sizeof(double));
		}

        OS = malloc(NumSol*sizeof(double*[3][3]));
        OSA = malloc(NumSol*sizeof(double*[3][3]));
        OSB = malloc(NumSol*sizeof(double*[3][3]));
        OSR = malloc(NumSol*sizeof(double*[3][3]));

	SS = malloc(NumCys*sizeof(double*));
	SSA = malloc(NumCys*sizeof(double*));
	SSR = malloc(NumCys*sizeof(double*));
	SSB = malloc(NumCys*sizeof(double*));

	for(i=0;i<NumCys;i++)
		{
		SS[i] = malloc(NumCys*sizeof(double*[4][3]));
		SSA[i] = malloc(NumCys*sizeof(double*[4][3]));
		SSR[i] = malloc(NumCys*sizeof(double*[4][3]));
		SSB[i] = malloc(NumCys*sizeof(double*[4][3]));
		}

        CT = malloc((ProLen-1)*sizeof(double*));
        CTA = malloc((ProLen-1)*sizeof(double*));
        CTR = malloc((ProLen-1)*sizeof(double*));
        CTB = malloc((ProLen-1)*sizeof(double*));

	for(i=0;i<ProLen-1;i++)
		{
		sizeCT = (RT[i+1] != 11) ? 6 : 5;
		CT[i] = malloc(sizeCT*sizeof(double*[3]));
		CTA[i] = malloc(sizeCT*sizeof(double*[3]));
		CTR[i] = malloc(sizeCT*sizeof(double*[3]));
		CTB[i] = malloc(sizeCT*sizeof(double*[3]));
		}

	SC = malloc(ProLen*sizeof(double**));
	SCA = malloc(ProLen*sizeof(double**));
	SCR = malloc(ProLen*sizeof(double**));
	SCB = malloc(ProLen*sizeof(double**));

	for(i=0;i<ProLen;i++)
		{
		SC[i] = malloc(NumSegsRP[i]*sizeof(double*));
		SCA[i] = malloc(NumSegsRP[i]*sizeof(double*));
		SCR[i] = malloc(NumSegsRP[i]*sizeof(double*));
		SCB[i] = malloc(NumSegsRP[i]*sizeof(double*));

		for(j=0;j<NumSegsRP[i];j++)
			{
			SC[i][j] = malloc(SegSizeRP[i][j]*sizeof(double*[3]));
			SCA[i][j] = malloc(SegSizeRP[i][j]*sizeof(double*[3]));
			SCR[i][j] = malloc(SegSizeRP[i][j]*sizeof(double*[3]));
			SCB[i][j] = malloc(SegSizeRP[i][j]*sizeof(double*[3]));
			}
		}


        VB = malloc(ProLen*sizeof(double***));
        VBA = malloc(ProLen*sizeof(double***));
        VBR = malloc(ProLen*sizeof(double***));
        VBB = malloc(ProLen*sizeof(double***));

        for(i=0; i<ProLen; ++i)
                {
                VB[i] = malloc(ProLen*sizeof(double**));
                VBA[i] = malloc(ProLen*sizeof(double**));
                VBR[i] = malloc(ProLen*sizeof(double**));
                VBB[i] = malloc(ProLen*sizeof(double**));
		
		for(j=0; j<ProLen; ++j)
			{
			VB[i][j] = malloc(lenRP[i]*sizeof(double*));
			VBA[i][j] = malloc(lenRP[i]*sizeof(double*));
			VBR[i][j] = malloc(lenRP[i]*sizeof(double*));
			VBB[i][j] = malloc(lenRP[i]*sizeof(double*));

			for(k=0; k<lenRP[i]; ++k)
				{
				VB[i][j][k] = malloc(lenRP[j]*sizeof(double*[2][3]));
				VBA[i][j][k] = malloc(lenRP[j]*sizeof(double*[2][3]));
				VBR[i][j][k] = malloc(lenRP[j]*sizeof(double*[2][3]));
				VBB[i][j][k] = malloc(lenRP[j]*sizeof(double*[2][3]));
				}
			}
		}

        VS = malloc(ProLen*sizeof(double**));
        VSA = malloc(ProLen*sizeof(double**));
        VSR = malloc(ProLen*sizeof(double**));
        VSB = malloc(ProLen*sizeof(double**));

        for(i=0; i<ProLen; ++i)
                {
                VS[i] = malloc(lenRP[i]*sizeof(double*));
                VSA[i] = malloc(lenRP[i]*sizeof(double*));
                VSR[i] = malloc(lenRP[i]*sizeof(double*));
                VSB[i] = malloc(lenRP[i]*sizeof(double*));

                for(j=0; j<lenRP[i]; ++j)
                        {
                        VS[i][j] = malloc(NumSol*sizeof(double*[3][2][2][3]));
                        VSA[i][j] = malloc(NumSol*sizeof(double*[3][2][2][3]));
                        VSR[i][j] = malloc(NumSol*sizeof(double*[3][2][2][3]));
                        VSB[i][j] = malloc(NumSol*sizeof(double*[3][2][2][3]));
                        }
                }

        VO = malloc(NumSol*sizeof(double*));
        VOA = malloc(NumSol*sizeof(double*));
        VOR = malloc(NumSol*sizeof(double*));
        VOB = malloc(NumSol*sizeof(double*));

        for(i=0; i<NumSol; ++i)
                {
                VO[i] = malloc(NumSol*sizeof(double*[3][3][2][3]));
                VOA[i] = malloc(NumSol*sizeof(double*[3][3][2][3]));
                VOR[i] = malloc(NumSol*sizeof(double*[3][3][2][3]));
                VOB[i] = malloc(NumSol*sizeof(double*[3][3][2][3]));
                }

        EB = malloc(ProLen*sizeof(double***));
        EBA = malloc(ProLen*sizeof(double***));
        EBR = malloc(ProLen*sizeof(double***));
        EBB = malloc(ProLen*sizeof(double***));

        for(i=0; i<ProLen; ++i)
                {
                EB[i] = malloc(ProLen*sizeof(double**));
                EBA[i] = malloc(ProLen*sizeof(double**));
                EBR[i] = malloc(ProLen*sizeof(double**));
                EBB[i] = malloc(ProLen*sizeof(double**));

		for(j=0; j<ProLen; ++j)
			{
			EB[i][j] = malloc(lenRP[i]*sizeof(double*));
			EBA[i][j] = malloc(lenRP[i]*sizeof(double*));
			EBR[i][j] = malloc(lenRP[i]*sizeof(double*));
			EBB[i][j] = malloc(lenRP[i]*sizeof(double*));
			
			for(m=0;m<lenRP[i];++m)
				{
				EB[i][j][m] = malloc(lenRP[j]*sizeof(double*[2]));
				EBA[i][j][m] = malloc(lenRP[j]*sizeof(double*[2]));
				EBR[i][j][m] = malloc(lenRP[j]*sizeof(double*[2]));
				EBB[i][j][m] = malloc(lenRP[j]*sizeof(double*[2]));
				}
			}
		}

        ES = malloc(ProLen*sizeof(double**));
        ESA = malloc(ProLen*sizeof(double**));
        ESR = malloc(ProLen*sizeof(double**));
        ESB = malloc(ProLen*sizeof(double**));

        for(i=0; i<ProLen; ++i)
                {
                ES[i] = malloc(lenRP[i]*sizeof(double*));
                ESA[i] = malloc(lenRP[i]*sizeof(double*));
                ESR[i] = malloc(lenRP[i]*sizeof(double*));
                ESB[i] = malloc(lenRP[i]*sizeof(double*));

                for(m=0;m<lenRP[i];++m)
                        {
                        ES[i][m] = malloc(NumSol*sizeof(double*[3][2]));
                        ESA[i][m] = malloc(NumSol*sizeof(double*[3][2]));
                        ESR[i][m] = malloc(NumSol*sizeof(double*[3][2]));
                        ESB[i][m] = malloc(NumSol*sizeof(double*[3][2]));
                        }
                }

        EO = malloc(NumSol*sizeof(double*));
        EOA = malloc(NumSol*sizeof(double*));
        EOR = malloc(NumSol*sizeof(double*));
        EOB = malloc(NumSol*sizeof(double*));

        for(i=0; i<NumSol; ++i)
                {
                EO[i] = malloc(NumSol*sizeof(double*[3][3][2]));
                EOA[i] = malloc(NumSol*sizeof(double*[3][3][2]));
                EOR[i] = malloc(NumSol*sizeof(double*[3][3][2]));
                EOB[i] = malloc(NumSol*sizeof(double*[3][3][2]));
                }

        etaEB = malloc(ProLen*sizeof(double***));
        EBerr = malloc(ProLen*sizeof(double***));

        for (i=0;i<ProLen;++i)
                {
                etaEB[i] = malloc(ProLen*sizeof(double**));
                EBerr[i] = malloc(ProLen*sizeof(double**));

		for(j=0;j<ProLen;++j)
			{
			etaEB[i][j] =  malloc(lenRP[i]*sizeof(double*));
			EBerr[i][j] =  malloc(lenRP[i]*sizeof(double*));

			for(m=0;m<lenRP[i];++m)
				{
				etaEB[i][j][m] =  malloc(lenRP[j]*sizeof(double*[2]));
				EBerr[i][j][m] =  malloc(lenRP[j]*sizeof(double*[2]));
				}
			}
                }

        etaES = malloc(ProLen*sizeof(double**));
        ESerr = malloc(ProLen*sizeof(double**));

        for (i=0;i<ProLen;++i)
                {
                etaES[i] = malloc(lenRP[i]*sizeof(double*));
                ESerr[i] = malloc(lenRP[i]*sizeof(double*));

                for(m=0;m<lenRP[i];++m)
                        {
                        etaES[i][m] = malloc(NumSol*sizeof(double*[3][2]));
                        ESerr[i][m] = malloc(NumSol*sizeof(double*[3][2]));
                        }
                }


        etaEO = malloc(NumSol*sizeof(double*));
        EOerr = malloc(NumSol*sizeof(double*));

        for (i=0;i<NumSol;++i)
                {
                etaEO[i] = malloc(NumSol*sizeof(double*[3][3][2]));
                EOerr[i] = malloc(NumSol*sizeof(double*[3][3][2]));
                }

        etaHB = malloc(NumDN*sizeof(double*));
        HBerr = malloc(NumDN*sizeof(double*));

        for(i=0;i<NumDN;i++)
                {
                etaHB[i] = malloc(NumAC*sizeof(double));
                HBerr[i] = malloc(NumAC*sizeof(double));
                }

        etaAD = malloc(NumAD*sizeof(double));
        ADerr = malloc(NumAD*sizeof(double));

        etaOS = malloc(NumSol*sizeof(double));
        OSerr = malloc(NumSol*sizeof(double));

        etaSS = malloc(NumCys*sizeof(double*));
        SSerr = malloc(NumCys*sizeof(double*));

	for(i=0;i<NumCys;i++)
		{
		etaSS[i] = malloc(NumCys*sizeof(double));
		SSerr[i] = malloc(NumCys*sizeof(double));
		}

        etaBL = malloc(NumBL*sizeof(double));
        BLerr = malloc(NumBL*sizeof(double));

        etaDH = malloc(NumDH*sizeof(double*));
        DHerr = malloc(NumDH*sizeof(double*));

	for(i=0;i<NumDH;i++)
		{
		etaDH[i] = malloc(NumDHpar[i]*sizeof(double));
		DHerr[i] = malloc(NumDHpar[i]*sizeof(double));
		}

        etaCT = malloc((ProLen-1)*sizeof(double));
        CTerr = malloc((ProLen-1)*sizeof(double));

        etaSC = malloc(ProLen*sizeof(double*));
        SCerr = malloc(ProLen*sizeof(double*));

	for(i=0;i<ProLen;i++)
		{
		etaSC[i] = malloc(NumSegsRP[i]*sizeof(double));
		SCerr[i] = malloc(NumSegsRP[i]*sizeof(double));
		}

        HBStore = malloc(NumDN*sizeof(double*));
        for (i=0;i<NumDN;++i)
                HBStore[i] =  malloc(NumAC*sizeof(double*[4][3]));

	ADStore = malloc(NumAD*sizeof(double*[8][3]));

	SSStore = malloc(NumCys*sizeof(double*));
	for (i=0;i<NumCys;++i)
		SSStore[i] =  malloc(NumCys*sizeof(double*[4][3]));

	SCStore = malloc(ProLen*sizeof(double**));

	for (i=0;i<ProLen;++i)
		{
		SCStore[i] = malloc(NumSegsRP[i]*sizeof(double*));
		
		for (j=0;j<NumSegsRP[i];++j)	
			SCStore[i][j] = malloc(SegSizeRP[i][j]*sizeof(double*[3]));
		}

	}


int getRigidMotifs()
        {
        int i,j,k,r1,z,m,n,sizeCT;
        FILE *fp;
        char buf[1024];
	unsigned char a1;

        //char file[50];

	chdir("../Files4RRR");

        fp = fopen("CTMotifs.txt","r");
        if(!fp)
                {
                printf("CT motif file not found\n");
                return 0;
                }

        while(fgets(buf, sizeof buf, fp))
                {
                if (strlen(buf)<10 && buf[0] != '\n')
                        {
                        sscanf(buf, "%c%*[ ]%d",&a1,&z);
                        r1 = AA[a1];
			CTMotifNum[r1] = z;

			sizeCT = (r1 != 11) ? 6 : 5;
			CTMotifs[r1] = malloc(CTMotifNum[r1]*sizeof(double*));
			for (z=0;z<CTMotifNum[r1];++z)
				CTMotifs[r1][z] = malloc(sizeCT*sizeof(double*[3]));

			for (z=0;z<CTMotifNum[r1];++z)
				for(m=0;m<sizeCT;++m)
					for(n=0;n<3;++n)
						fscanf(fp,"%lf%*[ \t ]%*[\n]",&CTMotifs[r1][z][m][n]);
			}
		}

        fclose(fp);

	chdir("../GenFiles");

        fp=fopen("OtherMotifs.txt","r");
        if(!fp)
                {    
                printf("Other Motifs file not found\n");
                return 0;
                }    

        for(j=0;j<5;++j)
                for(k=0;k<3;++k)
                        fscanf(fp,"%lf%*[ \t ]%*[\n]",&NTMotif[j][k]);

        for(j=0;j<4;++j)
                for(k=0;k<3;++k)
                        fscanf(fp,"%lf%*[ \t ]%*[\n]",&OXTMotif[j][k]);

        for(j=0;j<3;++j)
                for(k=0;k<3;++k)
                        fscanf(fp,"%lf%*[ \t ]%*[\n]",&SolMotif[j][k]);

        for(j=0;j<4;++j)
                for(k=0;k<3;++k)
                        fscanf(fp,"%lf%*[ \t ]%*[\n]",&ADMotif[j][k]);

        fclose(fp);


	//Disulphide bonds
        fp = fopen("SSBonds.txt","r");
        if(!fp)
                {
                printf("SSBonds.txt not found\n");
                return 0;
                }

	fscanf(fp,"%hd%*[\n]",&SSMotifNum);

	SSMotifs = malloc(SSMotifNum*sizeof(double*[4][3]));

	for(int z=0;z<SSMotifNum;++z)
		for(i=0;i<4;++i)
			for(j=0;j<3;++j)
				fscanf(fp,"%lf%*[ \t ]%*[\n]",&SSMotifs[z][i][j]);


        fclose(fp);

	chdir("../SSHB");

        return 1;
	}


int getMotifsSC()
        {
        int i,j,m,n,z,r1,tot_motifs;
        FILE *fp;
        char buf[1024];
	unsigned char a1;

	chdir("../Files4RRR");

        fp=fopen("SCMotifs.txt","r");
        if(!fp)
                {
                printf("SC motif_file not found\n");
                return 0;
                }

	for (i=0;i<21;++i)
		{
		SCMotifNum[i] = malloc((NumSegs[i]+1)*sizeof(int));
		SCMotifs[i] = malloc((NumSegs[i]+1)*sizeof(double**));
		}

	tot_motifs = 0;

        while(fgets(buf, sizeof buf, fp))
                {
                if (strlen(buf)<10 && buf[0] != '\n')
                        {
                        sscanf(buf, "%c%*[ ]%d%*[ ]%d",&a1,&j,&z);
                        r1 = AA[a1];
			SCMotifNum[r1][j] = z;
			tot_motifs += z;

			SCMotifs[r1][j] = malloc(SCMotifNum[r1][j]*sizeof(double*));

			for (z=0;z<SCMotifNum[r1][j];++z)
				{
				SCMotifs[r1][j][z] = malloc(SegSize[r1][j]*sizeof(double*[3]));

				for(m=0;m<SegSize[r1][j];++m)
					for(n=0;n<3;++n)
						fscanf(fp,"%lf%*[ \t ]%*[\n]",&SCMotifs[r1][j][z][m][n]);
				}
			}
		}

	printf("\ntotal SC motifs:%d\n",tot_motifs);
	printf("avg SC motifs per type: %lf\n",(double) tot_motifs/60.);

        fclose(fp);

	chdir("../SSHB");

        return 1;
        }

int getprotein(char *protein_file)
        {
        int i,j,k,m,n,CountC,p,q;
        FILE *fp;
        char buf[1024];
        char file[50];
        //unsigned char aa;

	chdir("../Files4RRR");

        fp=fopen(protein_file,"r");
        if(!fp)
                {
                printf("protein_file not found\n");
                return 0;
                }

        fscanf(fp,"%d%*[ ]%d%*[\n]",&ProLen,&NumSolEx);

        RT=malloc(ProLen*sizeof(int));

        fgets(buf, sizeof buf, fp);

        fclose(fp);

	NumCys = 0;

        for(i=0;i<ProLen;++i)
                {
                RT[i] = AA[(unsigned char)buf[i]];
                if (RT[i] == 8 )
                        NumCys++;
                }

        CysIndex = malloc(NumCys*sizeof(int));

        CountC = 0;
        for(i=0;i<ProLen;++i)
                {
                if (RT[i] == 8)
                        {
                        CysIndex[CountC] = i;
                        CountC++;
                        }
                }

        lenRP = malloc(ProLen*sizeof(double));
        NumSegsRP = malloc(ProLen*sizeof(int));
        SegSizeRP = malloc(ProLen*sizeof(int*));
        PQ = malloc(ProLen*sizeof(double*));
        ATRP = malloc(ProLen*sizeof(int*));

        for(i=0;i<ProLen;++i)
                {
                if (i == 0)
                        {
                        lenRP[i] = lenAAH[RT[i]]+2;
                        NumSegsRP[i] = NumSegs[RT[i]]+1;
                        }

                else if (i == ProLen-1)
                        {
                        lenRP[i] = lenAAH[RT[i]]+1;
                        NumSegsRP[i] = NumSegs[RT[i]]+1;
                        }

                else
                        {
                        lenRP[i] = lenAAH[RT[i]];
                        NumSegsRP[i] = NumSegs[RT[i]];
                        }
		}

        for(i=0;i<ProLen;++i)
                {
                ATRP[i] = malloc((lenRP[i])*sizeof(int));
                for(j=0;j<lenAAH[RT[i]];++j)
                        ATRP[i][j] = ATR[RT[i]][j];


                SegSizeRP[i] = malloc((NumSegsRP[i])*sizeof(int));

                for(j=0;j<NumSegs[RT[i]];++j)
                        SegSizeRP[i][j] = SegSize[RT[i]][j];
		}

	SegSizeRP[0][NumSegs[RT[0]]] = 5;
	SegSizeRP[ProLen-1][NumSegs[RT[ProLen-1]]] = 4;

        SegAtomRP = malloc(ProLen*sizeof(int**));
        for(i=0;i<ProLen;++i)
		{
                SegAtomRP[i] = malloc((NumSegsRP[i])*sizeof(int*));

                for(j=0;j<NumSegsRP[i];++j)
                        SegAtomRP[i][j] = malloc((SegSizeRP[i][j])*sizeof(int));
		}

        for(i=0;i<ProLen;++i)
                for(j=0;j<NumSegs[RT[i]];++j)
                        for(k=0;k<SegSize[RT[i]][j];++k)
                                SegAtomRP[i][j][k] = SegAtom[RT[i]][j][k];

	memcpy(SegAtomRP[0][NumSegs[RT[0]]], (const int []) {0,1,lenRP[0]-3,lenRP[0]-2,lenRP[0]-1}, 5 * sizeof(int));
	memcpy(SegAtomRP[ProLen-1][NumSegs[RT[ProLen-1]]], (const int []) {1,2,3,lenRP[ProLen-1]-1}, 4 * sizeof(int));

        ATRP[0][lenRP[0]-1] = 4;
        ATRP[0][lenRP[0]-2] = 4;
        ATRP[ProLen-1][lenRP[ProLen-1]-1] = 2;

        PQS[1] = 0.417;
        PQS[2] = 0.417;
        PQS[0] = -2.0*PQS[1];

	parSol[0][0] = -0.834;
	parSol[1][0] = 0.417;
	parSol[2][0] = 0.417;

	parSol[0][1] = -0.1521;
	parSol[1][1] = -0.046;
	parSol[2][1] = -0.046;

	parSol[0][2] = 1.7682;
	parSol[1][2] = 0.2245;
	parSol[2][2] = 0.2245;

        SegIntRP = malloc(ProLen*sizeof(int*));

        for(i=0;i<ProLen;++i)
                SegIntRP[i] = malloc(NumSegsRP[i]*sizeof(int));

        for(i=0;i<ProLen;++i)
                for(j=0;j<NumSegs[RT[i]];++j)
                        SegIntRP[i][j] = SegInt[RT[i]][j];

        SegIntRP[0][NumSegs[RT[0]]] = 1;
        SegIntRP[ProLen-1][NumSegs[RT[ProLen-1]]] = 1;

        atom = malloc(ProLen*sizeof(double*));

        for (i=0;i<ProLen;++i)
                atom[i] = malloc(lenRP[i]*sizeof(double*[3]));

        solatom = malloc(NumSol*sizeof(double*[3][3]));

	// Parameters and Protein information

        parRP = malloc(ProLen*sizeof(double*));

        for (i=0;i<ProLen;++i)
                parRP[i] = malloc(lenRP[i]*sizeof(double*[3]));

        snprintf(file, 50, "%s.par",ProID);

        fp = fopen(file,"r");
        if(!fp)
                {
                printf("Parameter file not found\n");
                return 0;
                }

	double temp[3];

        while(fgets(buf, sizeof buf, fp))
                {
		sscanf(buf, "%d%*[ ]%d%*[ ]%lf%*[ ]%lf%*[ ]%lf",&i,&j,&temp[0],&temp[1],&temp[2]);
		parRP[i][j][0] = temp[0];
		parRP[i][j][1] = temp[1];
		parRP[i][j][2] = temp[2];
		}

        fclose(fp);

        CONH = malloc(ProLen*sizeof(int***));
        INTRP = malloc(ProLen*sizeof(double***));
        VXRR = malloc(ProLen*sizeof(double***));

        for(i=0; i<ProLen; ++i)
		{
                CONH[i] = malloc(ProLen*sizeof(int**));
		INTRP[i] = malloc(ProLen*sizeof(double**));
		VXRR[i] = malloc(ProLen*sizeof(double**));

		for(j=0; j<ProLen; ++j)
			{
			CONH[i][j] = malloc(lenRP[i]*sizeof(int*));
			INTRP[i][j] = malloc(lenRP[i]*sizeof(double*));
			VXRR[i][j] = malloc(lenRP[i]*sizeof(double*));

			for(k=0; k<lenRP[i]; ++k)
				{
				CONH[i][j][k] = malloc(lenRP[j]*sizeof(int));
				INTRP[i][j][k] = malloc(lenRP[j]*sizeof(double*[3]));
				VXRR[i][j][k] = malloc(lenRP[j]*sizeof(double));

				for(m=0; m<lenRP[j]; ++m)
					{
					CONH[i][j][k][m] = 0;
					INTRP[i][j][k][m][0] = 0.;
					INTRP[i][j][k][m][1] = 0.;
					INTRP[i][j][k][m][2] = 0.;

					}	
				}
			}
		}

        snprintf(file, 50, "%s.a2r",ProID);

        fp = fopen(file,"r");
        if(!fp)
                {
                printf("atom mapping file not found\n");
                return 0;
                }

        fscanf(fp,"%d%*[\n]",&TotAtoms);

        A2R = malloc(TotAtoms*sizeof(int*[2]));

        R2A = malloc(ProLen*sizeof(int*));

        for(i=0; i<ProLen; ++i)
		R2A[i] = malloc(lenRP[i]*sizeof(int));

        while(fgets(buf, sizeof buf, fp))
                {
		sscanf(buf, "%d%*[ ]%d%*[ ]%d",&j,&i,&m);
		A2R[j][0] = i;
		A2R[j][1] = m;
		R2A[i][m] = j;
		}

        fclose(fp);

        snprintf(file, 50, "%s.nonb",ProID);

        fp = fopen(file,"r");
        if(!fp)
                {
                printf("interaction file not found\n");
                return 0;
                }

        while(fgets(buf, sizeof buf, fp))
                {
		sscanf(buf, "%d%*[ ]%d%*[ ]%lf%*[ ]%lf%*[ ]%lf",&p,&q,&temp[0],&temp[1],&temp[2]);
		i = A2R[p][0];
		m = A2R[p][1];
		j = A2R[q][0];
		n = A2R[q][1];

		CONH[i][j][m][n] = 1;
		INTRP[i][j][m][n][0] = temp[0];
		INTRP[i][j][m][n][1] = temp[1];
		INTRP[i][j][m][n][2] = temp[2];
		}

        fclose(fp);

	chdir("../SSHB");

        LJP = malloc(2501*sizeof(double*[2]));

        for(i=0;i<2501;++i)
                {
                LJP[i][0] = 0.5+(i*0.001);
                LJP[i][1] = pow(1/LJP[i][0], 12) - 2.*pow(1/LJP[i][0], 6);
                }

        WeightsSS = malloc(NumCys*sizeof(double*));
	for(i=0;i<NumCys;++i)
		WeightsSS[i] = malloc(NumCys*sizeof(double));

        //WeightsRS = malloc((ProLen-2)*sizeof(double));

        WeightsCT = malloc((ProLen-1)*sizeof(double));

	// Volume exclusion based on maximum LJ potential energy for each atom pair 

	double eps,sigma;
	double Emax = 4.;

	for(i=0; i<ProLen; ++i)
		for(j=i; j<ProLen; ++j)
			for(m=0; m<lenRP[i]; ++m)
				for(n=0; n<lenRP[j]; ++n)
					if (CONH[i][j][m][n] != 0)
						{
						sigma = (parRP[i][m][2] + parRP[j][n][2]);
						eps = sqrt(parRP[i][m][1]*parRP[j][n][1]);
						//VXRR[i][j][m][n] = VEX(sigma, eps, Emax);
						VXRR[i][j][m][n] = (-2. + sqrt(4.*(1. + Emax/eps)))/(2.*Emax/eps);
						VXRR[i][j][m][n] = pow(VXRR[i][j][m][n],1./6.);
						}


        VXRS = malloc(ProLen*sizeof(double*));
        for(i=0; i<ProLen; ++i)
		VXRS[i] = malloc(lenRP[i]*sizeof(double*[3]));

	Emax = 2.0;

        for (i=0;i<ProLen;++i)
                for(m=0;m<lenRP[i];++m)
			for(n=0; n<3; ++n)
				{
				sigma = (parRP[i][m][2] + parSol[n][2]);
				eps = sqrt(parRP[i][m][1]*parSol[n][1]);
				VXRS[i][m][n] = VEX(sigma, eps, Emax);
				}

	Emax = 5.0;

	for(m=0; m<3; ++m)
		for(n=0; n<3; ++n)
			{
			sigma = (parSol[m][2] + parSol[n][2]);
			eps = sqrt(parSol[m][1]*parSol[n][1]);
			VXSS[m][n] = VEX(sigma, eps, Emax);
			//printf("%d %d %lf\n",m,n,VXSS[m][n]);
			}

	VXSS[0][0] = 3.1; 

	// Hydrogen Bonding Book-keeping

        NumDN = 2;
        NumAC = ProLen+1;
 
        for(i=0;i<ProLen;++i)
                if (RT[i] != 11)
                        NumDN++;
                        
        HBBK[0] = malloc(NumDN*sizeof(int*[2]));
        HBBK[1] = malloc(NumAC*sizeof(int*[2]));

        int count;
        memcpy(HBBK[0][0], (const int []) {0,lenRP[0]-3},  2*sizeof(int));
        memcpy(HBBK[0][1], (const int []) {0,lenRP[0]-2},  2*sizeof(int));
        memcpy(HBBK[0][2], (const int []) {0,lenRP[0]-1},  2*sizeof(int));
        
        count = 3;
        for(i=1;i<ProLen;++i)
                if (RT[i] != 11)
                        {
                        memcpy(HBBK[0][count], (const int []) {i,lenAAH[RT[i]]-1},  2*sizeof(int));
                        count++;
                        }

        for(i=0;i<ProLen;++i)
                memcpy(HBBK[1][i], (const int []) {i,3},  2*sizeof(int));

        memcpy(HBBK[1][ProLen], (const int []) {ProLen-1,lenRP[ProLen-1]-1},  2*sizeof(int));

        WeightsHB = malloc(NumDN*sizeof(double*));
        for(i=0;i<NumDN;++i)
                WeightsHB[i] = malloc(NumAC*sizeof(double));

	// Secondary Structure Book-keeping

	int CountAD = 0;

	// Anti-parallel beta sheets 

        for(i=0;i<ProLen;++i)
		for(j=0;j<ProLen;++j)
			if ((i+2 < ProLen) && (j-2 >= 0) && (RT[j] != 11) && (RT[j-2] != 11) && (abs(i-j) > 2) && (abs(i-j+4) > 2))
				CountAD++;
	// Alpha-helix

        for(i=0;i<ProLen-5;++i)
		if (RT[i+5] != 11 && RT[i+4] != 11)
			CountAD++;

	NumAD = CountAD;
	ADBK = malloc(NumAD*sizeof(int*[8][2]));

	CountAD = 0;

        for(i=0;i<ProLen;++i)
		for(j=0;j<ProLen;++j)
			if ((i+2 < ProLen) && (j-2 >= 0) && (RT[j] != 11) && (RT[j-2] != 11) && (abs(i-j) > 2) && (abs(i-j+4) > 2))
				{
				memcpy(ADBK[CountAD][0], (const int []) {j,0},  2*sizeof(int));
				memcpy(ADBK[CountAD][1], (const int []) {j,lenAAH[RT[j]]-1},  2*sizeof(int));
				memcpy(ADBK[CountAD][2], (const int []) {i,3},  2*sizeof(int));
				memcpy(ADBK[CountAD][3], (const int []) {i,2},  2*sizeof(int));
	 
				memcpy(ADBK[CountAD][4], (const int []) {j-2,0},  2*sizeof(int));
				memcpy(ADBK[CountAD][5], (const int []) {j-2,lenAAH[RT[j-2]]-1},  2*sizeof(int));
				memcpy(ADBK[CountAD][6], (const int []) {i+2,3},  2*sizeof(int));
				memcpy(ADBK[CountAD][7], (const int []) {i+2,2},  2*sizeof(int));

				CountAD++;
				}
	
        for(i=0;i<ProLen-5;++i)
		if (RT[i+5] != 11 && RT[i+4] != 11)
			{
			memcpy(ADBK[CountAD][0], (const int []) {i+4,0},  2*sizeof(int));
			memcpy(ADBK[CountAD][1], (const int []) {i+4,lenAAH[RT[i+4]]-1},  2*sizeof(int));
			memcpy(ADBK[CountAD][2], (const int []) {i,3},  2*sizeof(int));
			memcpy(ADBK[CountAD][3], (const int []) {i,2},  2*sizeof(int));
 
			memcpy(ADBK[CountAD][4], (const int []) {i+5,0},  2*sizeof(int));
			memcpy(ADBK[CountAD][5], (const int []) {i+5,lenAAH[RT[i+5]]-1},  2*sizeof(int));
			memcpy(ADBK[CountAD][6], (const int []) {i+1,3},  2*sizeof(int));
			memcpy(ADBK[CountAD][7], (const int []) {i+1,2},  2*sizeof(int));

			CountAD++;
			}





/*

	int CountAD = 0;

	// Anti-parallel beta sheets 

        for(i=0;i<ProLen;++i)
		for(j=i+4;j<ProLen;++j)
			if (RT[i] != 11 && RT[j] != 11)
				CountAD++;

	// Alpha-helix

        for(i=0;i<ProLen-5;++i)
		if (RT[i+5] != 11 && RT[i+4] != 11)
			CountAD++;

	NumAD = CountAD;
	ADBK = malloc(NumAD*sizeof(int*[8][2]));

	CountAD = 0;

        for(i=0;i<ProLen;++i)
		for(j=i+4;j<ProLen;++j)
			if (RT[i] != 11 && RT[j] != 11)
				{
				memcpy(ADBK[CountAD][0], (const int []) {i,0},  2*sizeof(int));
				memcpy(ADBK[CountAD][1], (const int []) {i,lenAAH[RT[i]]-1},  2*sizeof(int));
				memcpy(ADBK[CountAD][2], (const int []) {j,3},  2*sizeof(int));
				memcpy(ADBK[CountAD][3], (const int []) {j,2},  2*sizeof(int));
	 
				memcpy(ADBK[CountAD][4], (const int []) {j,0},  2*sizeof(int));
				memcpy(ADBK[CountAD][5], (const int []) {j,lenAAH[RT[j]]-1},  2*sizeof(int));
				memcpy(ADBK[CountAD][6], (const int []) {i,3},  2*sizeof(int));
				memcpy(ADBK[CountAD][7], (const int []) {i,2},  2*sizeof(int));

				CountAD++;
				}
	
        for(i=0;i<ProLen-5;++i)
		if (RT[i+5] != 11 && RT[i+4] != 11)
			{
			memcpy(ADBK[CountAD][0], (const int []) {i+4,0},  2*sizeof(int));
			memcpy(ADBK[CountAD][1], (const int []) {i+4,lenAAH[RT[i+4]]-1},  2*sizeof(int));
			memcpy(ADBK[CountAD][2], (const int []) {i,3},  2*sizeof(int));
			memcpy(ADBK[CountAD][3], (const int []) {i,2},  2*sizeof(int));
 
			memcpy(ADBK[CountAD][4], (const int []) {i+5,0},  2*sizeof(int));
			memcpy(ADBK[CountAD][5], (const int []) {i+5,lenAAH[RT[i+5]]-1},  2*sizeof(int));
			memcpy(ADBK[CountAD][6], (const int []) {i+1,3},  2*sizeof(int));
			memcpy(ADBK[CountAD][7], (const int []) {i+1,2},  2*sizeof(int));

			CountAD++;
			}





	int CountAD = 0,DN[3][2],AC[3][2];
	int maxRSD = 3, m1,m2;

        for(i=0;i<ProLen-1;++i)
		for(j=i+1;j<ProLen-1;++j)
			{
			DN[0][0] = i+1;   
			AC[0][0] = j;

			DN[0][1] = j+1;
			AC[0][1] = i;

			DN[1][0] = i+2;   
			AC[1][0] = j+1;

			DN[1][1] = j+2;
			AC[1][1] = i+1;

			DN[2][0] = i+2;   
			AC[2][0] = j-1;

			DN[2][1] = j;
			AC[2][1] = i+1;

			for (m1=0;m1<2;++m1)
				for (m2=0;m2<2;++m2)
					{
					if (RT[DN[0][m1]] == 11)
						continue;

					if ((abs(DN[0][m1] - AC[0][m1]) > maxRSD) && (abs(DN[1][m2] - AC[1][m2]) > maxRSD))
						if ((i+1 < ProLen-1) && (j+1 < ProLen-1) && (RT[DN[1][m2]] != 11)) 
							CountAD++;

					if ((abs(DN[0][m1] - AC[0][m1]) > maxRSD) && (abs(DN[2][m2] - AC[2][m2]) > maxRSD))
						if ((i+1 != ProLen-1) && (RT[DN[2][m2]] != 11))
							CountAD++;
					}
			}

	NumAD = CountAD;

	ADBK = malloc(NumAD*sizeof(int*[8][2]));
	ADkey = malloc(NumAD*sizeof(int*[3]));

	ADCheck = malloc((ProLen-1)*sizeof(int*));

        for(i=0;i<ProLen-1;++i)
		ADCheck[i] = malloc((ProLen-1)*sizeof(int*[2]));


	CountAD = 0;
        for(i=0;i<ProLen-1;++i)
		for(j=i+1;j<ProLen-1;++j)
			{
			DN[0][0] = i+1;   
			AC[0][0] = j;

			DN[0][1] = j+1;
			AC[0][1] = i;

			DN[1][0] = i+2;   
			AC[1][0] = j+1;

			DN[1][1] = j+2;
			AC[1][1] = i+1;

			DN[2][0] = i+2;   
			AC[2][0] = j-1;

			DN[2][1] = j;
			AC[2][1] = i+1;

			for (m1=0;m1<2;++m1)
				for (m2=0;m2<2;++m2)
					{
					if (RT[DN[0][m1]] == 11)
						continue;

					if ((abs(DN[0][m1] - AC[0][m1]) > maxRSD) && (abs(DN[1][m2] - AC[1][m2]) > maxRSD))
						if ((i+1 < ProLen-1) && (j+1 < ProLen-1) && (RT[DN[1][m2]] != 11)) 
							{
							memcpy(ADBK[CountAD][0], (const int []) {DN[0][m1],0},  2*sizeof(int));
							memcpy(ADBK[CountAD][1], (const int []) {DN[0][m1],lenAAH[RT[DN[0][m1]]]-1},  2*sizeof(int));
							memcpy(ADBK[CountAD][2], (const int []) {AC[0][m1],3},  2*sizeof(int));
							memcpy(ADBK[CountAD][3], (const int []) {AC[0][m1],2},  2*sizeof(int));
				 
							memcpy(ADBK[CountAD][4], (const int []) {DN[1][m2],0},  2*sizeof(int));
							memcpy(ADBK[CountAD][5], (const int []) {DN[1][m2],lenAAH[RT[DN[1][m2]]]-1},  2*sizeof(int));
							memcpy(ADBK[CountAD][6], (const int []) {AC[1][m2],3},  2*sizeof(int));
							memcpy(ADBK[CountAD][7], (const int []) {AC[1][m2],2},  2*sizeof(int));

							ADkey[CountAD][0] = i;
							ADkey[CountAD][1] = j;
							ADkey[CountAD][2] = 0;

							CountAD++;
							}
	
					if ((abs(DN[0][m1] - AC[0][m1]) > maxRSD) && (abs(DN[2][m2] - AC[2][m2]) > maxRSD))
						if ((i+1 != ProLen-1) && (RT[DN[2][m2]] != 11))
							{
							memcpy(ADBK[CountAD][0], (const int []) {DN[0][m1],0},  2*sizeof(int));
							memcpy(ADBK[CountAD][1], (const int []) {DN[0][m1],lenAAH[RT[DN[0][m1]]]-1},  2*sizeof(int));
							memcpy(ADBK[CountAD][2], (const int []) {AC[0][m1],3},  2*sizeof(int));
							memcpy(ADBK[CountAD][3], (const int []) {AC[0][m1],2},  2*sizeof(int));
				 
							memcpy(ADBK[CountAD][4], (const int []) {DN[2][m2],0},  2*sizeof(int));
							memcpy(ADBK[CountAD][5], (const int []) {DN[2][m2],lenAAH[RT[DN[2][m2]]]-1},  2*sizeof(int));
							memcpy(ADBK[CountAD][6], (const int []) {AC[2][m2],3},  2*sizeof(int));
							memcpy(ADBK[CountAD][7], (const int []) {AC[2][m2],2},  2*sizeof(int));

							ADkey[CountAD][0] = i;
							ADkey[CountAD][1] = j;
							ADkey[CountAD][2] = 1;

							CountAD++;
							}
					}
			}

	int CountAD = 0., countDN[2], countAC[2];
	int maxRSD = 4;
	int RD[4] = {-2,-1,1,2}; 

	NumAD = 0;

        for(i=0;i<ProLen;++i)
		for(j=0;j<ProLen;++j)
			{
			if (RT[i] == 11)
				continue;

			countDN[0] = (i == 0) ? 3 : 1;
			countAC[0] = (j == ProLen-1) ? 2 : 1;

			for (k=0;k<4;++k)
				{
				p = i + abs(RD[k]); 
				q = j + RD[k]; 

				if (p >= ProLen)
					continue;
				if ((q >= ProLen) || (q < 0))
					continue;

				if (RT[p] == 11)
					continue;

				countDN[1] = 1;
				countAC[1] = (q == ProLen-1) ? 2 : 1;

	
				if ((abs(i-j) >= maxRSD) && (abs(p-q) >= maxRSD))
					NumAD += countDN[0]*countAC[0]*countDN[1]*countAC[1]; 
				}
			}

	ADBK = malloc(NumAD*sizeof(int*[8][2]));

	int m1,m2,m3,typeOH[3][3];

        for(i=0;i<ProLen;++i)
		for(j=0;j<ProLen;++j)
			{
			if (RT[i] == 11)
				continue;

			countDN[0] = (i == 0) ? 3 : 1;
			countAC[0] = (j == ProLen-1) ? 2 : 1;

			for (k=0;k<4;++k)
				{
				p = i + abs(RD[k]); 
				q = j + RD[k]; 

				if (p >= ProLen)
					continue;
				if ((q >= ProLen) || (q < 0))
					continue;

				if (RT[p] == 11)
					continue;

				countDN[1] = 1;
				countAC[1] = (q == ProLen-1) ? 2 : 1;

				if (i == 0)
					{
					typeOH[0][0] = lenRP[0]-3;
					typeOH[0][1] = lenRP[0]-2;
					typeOH[0][2] = lenRP[0]-1;
					}
				else
					typeOH[0][0] = lenAAH[RT[i]]-1;

				if (j == ProLen-1)
					{
					typeOH[1][0] = 3;
					typeOH[1][1] = lenRP[ProLen-1]-1;
					}
				else
					typeOH[1][0] = 3;

				if (q == ProLen-1)
					{
					typeOH[2][0] = 3;
					typeOH[2][1] = lenRP[ProLen-1]-1;
					}
				else
					typeOH[2][0] = 3;
	
				if ((abs(i-j) >= maxRSD) && (abs(p-q) >= maxRSD))
					{
					for (m1=0;m1<countDN[0];++m1)
					for (m2=0;m2<countAC[0];++m2)
					for (m3=0;m3<countAC[1];++m3)
						{
						memcpy(ADBK[CountAD][0], (const int []) {i,0},  2*sizeof(int));
						memcpy(ADBK[CountAD][1], (const int []) {i,typeOH[0][m1]},  2*sizeof(int));
						memcpy(ADBK[CountAD][2], (const int []) {j,typeOH[1][m2]},  2*sizeof(int));
						memcpy(ADBK[CountAD][3], (const int []) {j,2},  2*sizeof(int));
						
						memcpy(ADBK[CountAD][4], (const int []) {p,0},  2*sizeof(int));
						memcpy(ADBK[CountAD][5], (const int []) {p,lenAAH[RT[p]]-1},  2*sizeof(int));
						memcpy(ADBK[CountAD][6], (const int []) {q,typeOH[1][m3]},  2*sizeof(int));
						memcpy(ADBK[CountAD][7], (const int []) {q,2},  2*sizeof(int));

						CountAD++;
						}
					}
				}
			}

	printf("%d\n",NumAD);

	for (i=0;i<NumAD;++i)
		{
		for (m=0;m<8;++m)
			printf("%d %d \t",ADBK[i][m][0],ADBK[i][m][1]);
		printf("\n");
		}
	
	NumDN = 0;
        for(i=1;i<ProLen;++i)
                if (RT[i] != 11)
			NumDN++;

	NumAC = ProLen-1;

	ADBK[0] = malloc(NumDN*sizeof(int*[2]));
	ADBK[1] = malloc(NumAC*sizeof(int*[2]));

	int count = 0;

        for(i=1;i<ProLen;++i)
                if (RT[i] != 11)
			{
			memcpy(ADBK[0][count], (const int []) {i,lenAAH[RT[i]]-1},  2*sizeof(int));
			count++;
			}

        for(i=0;i<ProLen-1;++i)
		memcpy(ADBK[1][i], (const int []) {i,3},  2*sizeof(int));

	NumDN = 2;
	NumAC = ProLen+1;
 
        for(i=0;i<ProLen;++i)
                if (RT[i] != 11)
			NumDN++;
			
	ADBK[0] = malloc(NumDN*sizeof(int*[2]));
	ADBK[1] = malloc(NumAC*sizeof(int*[2]));

	int count;
	memcpy(ADBK[0][0], (const int []) {0,lenRP[0]-3},  2*sizeof(int));
	memcpy(ADBK[0][1], (const int []) {0,lenRP[0]-2},  2*sizeof(int));
	memcpy(ADBK[0][2], (const int []) {0,lenRP[0]-1},  2*sizeof(int));
	
	count = 3;
        for(i=1;i<ProLen;++i)
                if (RT[i] != 11)
			{
			memcpy(ADBK[0][count], (const int []) {i,lenAAH[RT[i]]-1},  2*sizeof(int));
			count++;
			}

        for(i=0;i<ProLen;++i)
		memcpy(ADBK[1][i], (const int []) {i,3},  2*sizeof(int));

	memcpy(ADBK[1][ProLen], (const int []) {ProLen-1,lenRP[ProLen-1]-1},  2*sizeof(int));
*/
        WeightsAD = malloc(NumAD*sizeof(double));

	ActualCoor = malloc(ProLen*sizeof(double*));

        for(i=0;i<ProLen;++i)
		ActualCoor[i] = malloc(lenRP[i]*sizeof(double*[3]));

        snprintf(file, 50, "%s.coor",ProID);

	chdir("../Files4RRR");

        fp = fopen(file,"r");
        if(!fp)
                printf("coor_file not found\n");

	for(i=0;i<ProLen;++i)
		for (j=0;j<lenRP[i];++j)
                        for(k=0;k<3;++k)
                                fscanf(fp,"%lf%*[ ]%*[\n]",&ActualCoor[i][j][k]);

        fclose(fp);

	chdir("../SSHB");

	int totatoms;

        totatoms = 0;
        for(i=0;i<ProLen;++i)
		for(j=0;j<lenRP[i];++j)
                        {
                        ++totatoms;
                        for(n=0;n<3;++n)
                                COMRP[n] += ActualCoor[i][j][n];
                        }

        for(n=0;n<3;++n)
                COMRP[n] /= totatoms;

        for(i=0;i<ProLen;++i)
		for(j=0;j<lenRP[i];++j)
                        for(n=0;n<3;++n)
				ActualCoor[i][j][n] -= COMRP[n];

        return 1;
        }


int getDHinfo()
        {
        int i,m;
	//int j,n;
        FILE *fp;
        char buf[1024];

        char file[50];
        snprintf(file, 50, "%s.dhd",ProID);

	chdir("../Files4RRR");

        fp=fopen(file,"r");
        if(!fp)
                {
                printf("dihedral file not found\n");
                return 0;
                }

        fscanf(fp,"%d%*[\n]",&NumBL);

	BLAtoms = malloc(NumBL*sizeof(int*[2]));
	BLpar = malloc(NumBL*sizeof(double*[2]));

	for(i=0;i<NumBL;++i)
		fscanf(fp,"%d%*[ ]%d%*[ ]%lf%*[ ]%lf%*[\n]",&BLAtoms[i][0],&BLAtoms[i][1],&BLpar[i][0],&BLpar[i][1]);

        fscanf(fp,"%d%*[\n]",&NumDH);

	DHAtoms = malloc(NumDH*sizeof(int*[4]));
	NumDHpar = malloc(NumDH*sizeof(int));
	DHpar = malloc(NumDH*sizeof(double*));

	i = 0;

        while(fgets(buf, sizeof buf, fp))
                {
		sscanf(buf, "%d%*[ \t ]%d%*[ ]%d%*[ ]%d%*[ ]%d",&NumDHpar[i],&DHAtoms[i][0],&DHAtoms[i][1],&DHAtoms[i][2],&DHAtoms[i][3]);

		DHpar[i] = malloc(NumDHpar[i]*sizeof(double*[3]));

		for(m=0;m<NumDHpar[i];m++)
			{
			fgets(buf, sizeof buf, fp);
			sscanf(buf,"%lf%*[ ]%lf%*[ ]%lf",&DHpar[i][m][0],&DHpar[i][m][1],&DHpar[i][m][2]);
			}

		++i;
		}

        fclose(fp);

	chdir("../SSHB");
/*
	for(i=0;i<NumBL;++i)
		printf("%d %d %d %lf %lf\n",i,BLAtoms[i][0],BLAtoms[i][1],BLpar[i][0],BLpar[i][1]);
 
	for(i=0;i<NumDH;++i)
		{
		printf("\n");
		printf("%d %d %d %d %d\n",i,DHAtoms[i][0],DHAtoms[i][1],DHAtoms[i][2],DHAtoms[i][3]);
		for(j=0;j<NumDHpar[i];++j)
			printf("%lf %lf %lf\n",DHpar[i][j][0],DHpar[i][j][1],DHpar[i][j][2]);
		}	
*/
        return 1;
        }


void AAmapping()
        {
        int i,j;
	//char a1;

        for(i=65;i<=90;++i)
                AA[i]=0;

        AA['R']=1;
        AA['N']=2;
        AA['D']=3;
        AA['Q']=4;
        AA['E']=5;
        AA['H']=6;
        AA['K']=7;

        AA['C']=8;
        AA['M']=9;
        AA['G']=10;
        AA['P']=11;
        AA['S']=12;
        AA['T']=13;
        AA['W']=14;
        AA['Y']=15;

        AA['I']=16;
        AA['L']=17;
        AA['V']=18;
        AA['F']=19;
        AA['A']=20;

	// number of heavy atoms in all the amino acids
	
	lenAA[1] = 11; 
	lenAA[2] = 8; 
	lenAA[3] = 8; 
	lenAA[4] = 9; 
	lenAA[5] = 9; 
	lenAA[6] = 10; 
	lenAA[7] = 9; 

	lenAA[8] = 6; 
	lenAA[9] = 8; 
	lenAA[10] = 4; 
	lenAA[11] = 7; 
	lenAA[12] = 6; 
	lenAA[13] = 7; 
	lenAA[14] = 14; 
	lenAA[15] = 12; 

	lenAA[16] = 8; 
	lenAA[17] = 8; 
	lenAA[18] = 7; 
	lenAA[19] = 11; 
	lenAA[20] = 5; 

	// number of total atoms in all the amino acids after peptide bonds on both sides
	
	lenAAH[1] = 24; 
	lenAAH[2] = 14; 
	lenAAH[3] = 12; 
	lenAAH[4] = 17; 
	lenAAH[5] = 15; 
	lenAAH[6] = 18; 
	lenAAH[7] = 22; 

	lenAAH[8] = 10; 
	lenAAH[9] = 17; 
	lenAAH[10] = 7; 
	lenAAH[11] = 14; 
	lenAAH[12] = 11; 
	lenAAH[13] = 14; 
	lenAAH[14] = 24; 
	lenAAH[15] = 21; 

	lenAAH[16] = 19; 
	lenAAH[17] = 19; 
	lenAAH[18] = 16; 
	lenAAH[19] = 20; 
	lenAAH[20] = 10; 

	// atom types in all the amino acids. Atoms "CNOS" are indexed by 0,1,2,3 respectively 

	for(i=1;i<21;++i)
		{
		ATR[i] = malloc(lenAAH[i]*sizeof(int));

		for(j=0;j<lenAA[i];++j)
			ATR[i][j] = 0;

		for(j=lenAA[i];j<lenAAH[i];++j)
			ATR[i][j] = 4;

		ATR[i][0] = 1;
		ATR[i][3] = 2;
		}
	
	ATR[AA['C']][5] = 3;
	ATR[AA['W']][8] = 1;
	ATR[AA['D']][6] = 2;
	ATR[AA['D']][7] = 2;
	ATR[AA['N']][6] = 2;
	ATR[AA['N']][7] = 1;
	ATR[AA['Y']][11] = 2;
	ATR[AA['S']][5] = 2;
	ATR[AA['Q']][7] = 2;
	ATR[AA['Q']][8] = 1;
	ATR[AA['E']][7] = 2;
	ATR[AA['E']][8] = 2;
	ATR[AA['K']][8] = 1;
	ATR[AA['R']][7] = 1;
	ATR[AA['R']][9] = 1;
	ATR[AA['R']][10] = 1;
	ATR[AA['T']][5] = 2;
	ATR[AA['M']][6] = 3;
	ATR[AA['H']][6] = 1;
	ATR[AA['H']][9] = 1;
        }


void SCAtomMapping()
	{
        int i,j;

	NumSegs[0] = 0;
	NumSegs[1] = 5;
	NumSegs[2] = 3;
	NumSegs[3] = 3;
	NumSegs[4] = 4;
	NumSegs[5] = 4;
	NumSegs[6] = 3;
	NumSegs[7] = 6;
	NumSegs[8] = 2;
	NumSegs[9] = 5;
	NumSegs[10] = 1;
	NumSegs[11] = 2;
	NumSegs[12] = 3;
	NumSegs[13] = 4;
	NumSegs[14] = 3;
	NumSegs[15] = 3;
	NumSegs[16] = 5;
	NumSegs[17] = 5;
	NumSegs[18] = 4;
	NumSegs[19] = 3;
	NumSegs[20] = 2;

	for (i=0;i<21;++i)
		{
		SegSize[i] = malloc((NumSegs[i])*sizeof(int));
		SegAtom[i] = malloc((NumSegs[i])*sizeof(int*));

		SegInt[i] = malloc((NumSegs[i])*sizeof(int));
		}
	
	for (i=0;i<21;++i)
		for (j=0;j<NumSegs[i];++j)
			SegInt[i][j] = 0;

	SegInt[1][4] = 1;	
	SegInt[2][2] = 1;	
	SegInt[3][2] = 1;	
	SegInt[4][3] = 1;	
	SegInt[5][3] = 1;	
	SegInt[6][2] = 1;	
	SegInt[7][5] = 1;	
	SegInt[12][2] = 1;	
	SegInt[13][2] = 1;	
	SegInt[14][2] = 1;	
	SegInt[15][2] = 1;	

	memcpy(SegSize[1], (const int []) {5,5,5,5,10}, (NumSegs[1]) * sizeof(int));
	memcpy(SegSize[2], (const int []) {5,5,6}, (NumSegs[2]) * sizeof(int));
	memcpy(SegSize[3], (const int []) {5,5,4}, (NumSegs[3]) * sizeof(int));
	memcpy(SegSize[4], (const int []) {5,5,5,6}, (NumSegs[4]) * sizeof(int));
	memcpy(SegSize[5], (const int []) {5,5,5,4}, (NumSegs[5]) * sizeof(int));
	memcpy(SegSize[6], (const int []) {5,5,10}, (NumSegs[6]) * sizeof(int));
	memcpy(SegSize[7], (const int []) {5,5,5,5,5,5}, (NumSegs[7]) * sizeof(int));
	memcpy(SegSize[8], (const int []) {5,5}, (NumSegs[8]) * sizeof(int));
	memcpy(SegSize[9], (const int []) {5,5,5,3,5}, (NumSegs[9]) * sizeof(int));
	memcpy(SegSize[10], (const int []) {5}, (NumSegs[10]) * sizeof(int));
	memcpy(SegSize[11], (const int []) {5,11}, (NumSegs[11]) * sizeof(int));
	memcpy(SegSize[12], (const int []) {5,5,3}, (NumSegs[12]) * sizeof(int));
	memcpy(SegSize[13], (const int []) {5,5,3,5}, (NumSegs[13]) * sizeof(int));
	memcpy(SegSize[14], (const int []) {5,5,16}, (NumSegs[14]) * sizeof(int));
	memcpy(SegSize[15], (const int []) {5,5,13}, (NumSegs[15]) * sizeof(int));
	memcpy(SegSize[16], (const int []) {5,5,5,5,5}, (NumSegs[16]) * sizeof(int));
	memcpy(SegSize[17], (const int []) {5,5,5,5,5}, (NumSegs[17]) * sizeof(int));
	memcpy(SegSize[18], (const int []) {5,5,5,5}, (NumSegs[18]) * sizeof(int));
	memcpy(SegSize[19], (const int []) {5,5,12}, (NumSegs[19]) * sizeof(int));
	memcpy(SegSize[20], (const int []) {5,5}, (NumSegs[20]) * sizeof(int));

	for (i=0;i<21;++i)
		for (j=0;j<NumSegs[i];++j)
			SegAtom[i][j] = malloc(SegSize[i][j]*sizeof(int));

	memcpy(SegAtom[1][0], (const int []) {0,1,2,4,11}, SegSize[1][0] * sizeof(int));
	memcpy(SegAtom[1][1], (const int []) {1,4,5,12,13}, SegSize[1][1] * sizeof(int));
	memcpy(SegAtom[1][2], (const int []) {4,5,6,14,15}, SegSize[1][2] * sizeof(int));
	memcpy(SegAtom[1][3], (const int []) {5,6,7,16,17}, SegSize[1][3] * sizeof(int));
	memcpy(SegAtom[1][4], (const int []) {6,7,8,9,10,18,19,20,21,22}, SegSize[1][4] * sizeof(int));

	memcpy(SegAtom[2][0], (const int []) {0,1,2,4,8}, SegSize[2][0] * sizeof(int));
	memcpy(SegAtom[2][1], (const int []) {1,4,5,9,10}, SegSize[2][1] * sizeof(int));
	memcpy(SegAtom[2][2], (const int []) {4,5,6,7,11,12}, SegSize[2][2] * sizeof(int));

	memcpy(SegAtom[3][0], (const int []) {0,1,2,4,8}, SegSize[3][0] * sizeof(int));
	memcpy(SegAtom[3][1], (const int []) {1,4,5,9,10}, SegSize[3][1] * sizeof(int));
	memcpy(SegAtom[3][2], (const int []) {4,5,6,7}, SegSize[3][2] * sizeof(int));

	memcpy(SegAtom[4][0], (const int []) {0,1,2,4,9}, SegSize[4][0] * sizeof(int));
	memcpy(SegAtom[4][1], (const int []) {1,4,5,10,11}, SegSize[4][1] * sizeof(int));
	memcpy(SegAtom[4][2], (const int []) {4,5,6,12,13}, SegSize[4][2] * sizeof(int));
	memcpy(SegAtom[4][3], (const int []) {5,6,7,8,14,15}, SegSize[4][3] * sizeof(int));

	memcpy(SegAtom[5][0], (const int []) {0,1,2,4,9}, SegSize[5][0] * sizeof(int));
	memcpy(SegAtom[5][1], (const int []) {1,4,5,10,11}, SegSize[5][1] * sizeof(int));
	memcpy(SegAtom[5][2], (const int []) {4,5,6,12,13}, SegSize[5][2] * sizeof(int));
	memcpy(SegAtom[5][3], (const int []) {5,6,7,8}, SegSize[5][3] * sizeof(int));

	memcpy(SegAtom[6][0], (const int []) {0,1,2,4,10}, SegSize[6][0] * sizeof(int));
	memcpy(SegAtom[6][1], (const int []) {1,4,5,11,12}, SegSize[6][1] * sizeof(int));
	memcpy(SegAtom[6][2], (const int []) {4,5,6,7,8,9,13,14,15,16}, SegSize[6][2] * sizeof(int));

	memcpy(SegAtom[7][0], (const int []) {0,1,2,4,9}, SegSize[7][0] * sizeof(int));
	memcpy(SegAtom[7][1], (const int []) {1,4,5,10,11}, SegSize[7][1] * sizeof(int));
	memcpy(SegAtom[7][2], (const int []) {4,5,6,12,13}, SegSize[7][2] * sizeof(int));
	memcpy(SegAtom[7][3], (const int []) {5,6,7,14,15}, SegSize[7][3] * sizeof(int));
	memcpy(SegAtom[7][4], (const int []) {6,7,8,16,17}, SegSize[7][4] * sizeof(int));
	memcpy(SegAtom[7][5], (const int []) {7,8,18,19,20}, SegSize[7][5] * sizeof(int));

	memcpy(SegAtom[8][0], (const int []) {0,1,2,4,6}, SegSize[8][0] * sizeof(int));
	memcpy(SegAtom[8][1], (const int []) {1,4,5,7,8}, SegSize[8][1] * sizeof(int));

	memcpy(SegAtom[9][0], (const int []) {0,1,2,4,8}, SegSize[9][0] * sizeof(int));
	memcpy(SegAtom[9][1], (const int []) {1,4,5,9,10}, SegSize[9][1] * sizeof(int));
	memcpy(SegAtom[9][2], (const int []) {4,5,6,11,12}, SegSize[9][2] * sizeof(int));
	memcpy(SegAtom[9][3], (const int []) {5,6,7}, SegSize[9][3] * sizeof(int));
	memcpy(SegAtom[9][4], (const int []) {6,7,13,14,15}, SegSize[9][4] * sizeof(int));

	memcpy(SegAtom[10][0], (const int []) {0,1,2,4,5}, SegSize[10][0] * sizeof(int));

	memcpy(SegAtom[11][0], (const int []) {0,1,2,4,7}, SegSize[11][0] * sizeof(int));
	memcpy(SegAtom[11][1], (const int []) {0,1,4,5,6,8,9,10,11,12,13}, SegSize[11][1] * sizeof(int));

	memcpy(SegAtom[12][0], (const int []) {0,1,2,4,6}, SegSize[12][0] * sizeof(int));
	memcpy(SegAtom[12][1], (const int []) {1,4,5,7,8}, SegSize[12][1] * sizeof(int));
	memcpy(SegAtom[12][2], (const int []) {4,5,9}, SegSize[12][2] * sizeof(int));

	memcpy(SegAtom[13][0], (const int []) {0,1,2,4,7}, SegSize[13][0] * sizeof(int));
	memcpy(SegAtom[13][1], (const int []) {1,4,5,6,8}, SegSize[13][1] * sizeof(int));
	memcpy(SegAtom[13][2], (const int []) {4,5,9}, SegSize[13][2] * sizeof(int));
	memcpy(SegAtom[13][3], (const int []) {4,6,10,11,12}, SegSize[13][3] * sizeof(int));

	memcpy(SegAtom[14][0], (const int []) {0,1,2,4,14}, SegSize[14][0] * sizeof(int));
	memcpy(SegAtom[14][1], (const int []) {1,4,5,15,16}, SegSize[14][1] * sizeof(int));
	memcpy(SegAtom[14][2], (const int []) {4,5,6,7,8,9,10,11,12,13,17,18,19,20,21,22}, SegSize[14][2] * sizeof(int));

	memcpy(SegAtom[15][0], (const int []) {0,1,2,4,12}, SegSize[15][0] * sizeof(int));
	memcpy(SegAtom[15][1], (const int []) {1,4,5,13,14}, SegSize[15][1] * sizeof(int));
	memcpy(SegAtom[15][2], (const int []) {4,5,6,7,8,9,10,11,15,16,17,18,19}, SegSize[15][2] * sizeof(int));

	memcpy(SegAtom[16][0], (const int []) {0,1,2,4,8}, SegSize[16][0] * sizeof(int));
	memcpy(SegAtom[16][1], (const int []) {1,4,5,6,9}, SegSize[16][1] * sizeof(int));
	memcpy(SegAtom[16][2], (const int []) {4,5,7,10,11}, SegSize[16][2] * sizeof(int));
	memcpy(SegAtom[16][3], (const int []) {4,6,12,13,14}, SegSize[16][3] * sizeof(int));
	memcpy(SegAtom[16][4], (const int []) {5,7,15,16,17}, SegSize[16][4] * sizeof(int));

	memcpy(SegAtom[17][0], (const int []) {0,1,2,4,8}, SegSize[17][0] * sizeof(int));
	memcpy(SegAtom[17][1], (const int []) {1,4,5,9,10}, SegSize[17][1] * sizeof(int));
	memcpy(SegAtom[17][2], (const int []) {4,5,6,7,11}, SegSize[17][2] * sizeof(int));
	memcpy(SegAtom[17][3], (const int []) {5,6,12,13,14}, SegSize[17][3] * sizeof(int));
	memcpy(SegAtom[17][4], (const int []) {5,7,15,16,17}, SegSize[17][4] * sizeof(int));

	memcpy(SegAtom[18][0], (const int []) {0,1,2,4,7}, SegSize[18][0] * sizeof(int));
	memcpy(SegAtom[18][1], (const int []) {1,4,5,6,8}, SegSize[18][1] * sizeof(int));
	memcpy(SegAtom[18][2], (const int []) {4,5,9,10,11}, SegSize[18][2] * sizeof(int));
	memcpy(SegAtom[18][3], (const int []) {4,6,12,13,14}, SegSize[18][3] * sizeof(int));

	memcpy(SegAtom[19][0], (const int []) {0,1,2,4,11}, SegSize[19][0] * sizeof(int));
	memcpy(SegAtom[19][1], (const int []) {1,4,5,12,13}, SegSize[19][1] * sizeof(int));
	memcpy(SegAtom[19][2], (const int []) {4,5,6,7,8,9,10,14,15,16,17,18}, SegSize[19][2] * sizeof(int));

	memcpy(SegAtom[20][0], (const int []) {0,1,2,4,5}, SegSize[20][0] * sizeof(int));
	memcpy(SegAtom[20][1], (const int []) {1,4,6,7,8}, SegSize[20][1] * sizeof(int));
	
	}

/*
void init2()
	{
        int i,j,k,l,m,n,c;
        FILE *fp;
        double dis,rmin,eps,qmul,tot_energy = 0.,COM[3],theta,phi,a;
        char buf[1024];

        fp=fopen("t0.init","r");
        if(!fp)
                {
                printf("init file not found\n");
                }

	fgets(buf, sizeof buf, fp);
	fgets(buf, sizeof buf, fp);

	for(i=0;i<ProLen;++i)
		for (j=0;j<lenRP[i];++j)
                        for(k=0;k<3;++k)
                                fscanf(fp,"%lf%*[ ]%*[\n]",&atom[i][j][k]);

        fclose(fp);

        for(i=0;i<NumSol;++i)
		{
		theta = 0.5*urand(-M_PI,M_PI);
		phi = urand(0,2*M_PI);
		a = urand(RAD/5.,RAD);

		solatom[i][0][0] = a*cos(theta)*cos(phi);
		solatom[i][0][1] = a*cos(theta)*sin(phi);
		solatom[i][0][2] = a*sin(theta);

                for(n=1;n<3;++n)
                        for(k=0;k<3;++k)
				solatom[i][n][k] = solatom[i][0][k] + urand(-1,1);
		}

	for(i=0; i<ProLen; ++i)
		for(j=i; j<ProLen; ++j)
			for(m=0; m<lenRP[i]; ++m)
				for(n=0; n<lenRP[j]; ++n)
					if (CONH[i][j][m][n] != 0)
						{
						for(c=0; c<2; ++c)
							etaEB[i][j][m][n][c] = 1.;

						dis = 0.;
						for(k=0; k<3; ++k)
							dis += sq(atom[j][n][k] - atom[i][m][k]);
						dis = sqrt(dis);

						qmul = INTRP[i][j][m][n][0];
						eps = INTRP[i][j][m][n][1];
						rmin = INTRP[i][j][m][n][2];

						EB[i][j][m][n][0] = eps*(pow((rmin/dis), 12) - 2.*pow((rmin/dis), 6));
						EB[i][j][m][n][1] = 331.21*qmul/(dis*epsRR);
						
						tot_energy += EB[i][j][m][n][0];
						tot_energy += EB[i][j][m][n][1];

						//EB[i][j][m][n][0] += lambda*urand(-1,1);
						//EB[i][j][m][n][1] += lambda*urand(-1,1);

						if (EB[i][j][m][n][0] > 11.)
							EB[i][j][m][n][0] = 11.;
						}
						
	printf("Energy RR: %lf\n",tot_energy);

        for (i=0;i<ProLen;++i)
                for(m=0;m<lenRP[i];++m)
                        for (j=0;j<NumSol;++j)
                                for(n=0; n<3; ++n)
					for(c=0; c<2; ++c)
						{
						etaES[i][m][j][n][c] = 1.;
						ES[i][m][j][n][c] = urand(-1.,1.);
						}

        for(i=0;i<NumSol;++i)
                for(j=i+1;j<NumSol;++j)
                        for(m=0;m<3;++m)
                                for(n=0;n<3;++n)
                                        for(c=0;c<2;++c)
						{
						etaEO[i][j][m][n][c] = 1.;
						EO[i][j][m][n][c] = urand(-1.,1.);
						}

        for(i=0;i<NumSol;++i)
		etaOS[i] = 1.;

        for (i=0;i<ProLen-1;++i)
		etaCT[i] = 1.;

        for (i=0;i<ProLen;++i)
		for (j=0;j<NumSegsRP[i];++j)
			etaSC[i][j] = 1.;

	for (i=0;i<NumCys;++i)
		for (j=0;j<NumCys;++j)
			etaSS[i][j] = 1.;

	for (i=0;i<ProLen-2;++i)
		for (j=0;j<ProLen-2;++j)
			etaBB[i][j] = 1.;

        for(i=0;i<NumDH;++i)
		for(m=0;m<NumDHpar[i];++m)
			etaDH[i][m] = 1.;

	for(j=0;j<NumBL;++j)
		{
		dis = 0.;
  
		for(k=0;k<3;++k)
			dis += sq(atom[A2R[BLAtoms[j][0]][0]][A2R[BLAtoms[j][0]][1]][k] - atom[A2R[BLAtoms[j][1]][0]][A2R[BLAtoms[j][1]][1]][k]);
		dis = sqrt(dis);

		BE[j] = 0.5*BLpar[j][0]*sq(dis - BLpar[j][1]);
		tot_energy += BE[j];

		//printf("%d %d %d %lf %lf %lf %lf\n",j,BLAtoms[j][0],BLAtoms[j][1],BLpar[j][0],BLpar[j][1],dis,BE[j]);
		etaBL[j] = 1.;
		}

	double tmp[4][3],dhd;

	for(j=0;j<NumDH;++j)
		{
		for(i=0;i<4;++i)
			for(k=0; k<3; ++k)
				tmp[i][k] = atom[A2R[DHAtoms[j][i]][0]][A2R[DHAtoms[j][i]][1]][k];

		dhd = calc_dihed(tmp);
	
		for(m=0;m<NumDHpar[j];++m)
			{
			ED[j][m] = DHpar[j][m][0]*(1 + cos(DHpar[j][m][1]*dhd - DHpar[j][m][2]));
			tot_energy += ED[j][m];
			}
			//ED[j] += DHpar[j][0]*lambda*urand(-1,1);
		}

	printf("total initial energy: %lf\n",tot_energy);

	changeVar(VB, VS, VO, OS, SC, CT, SS, BB, DH, BL);
	}
*/

void init()
        {
        int i,j,k,m,n,c;
        double a,theta,phi;

	TOTCONST = 0;

        for(i=0;i<NumSol;++i)
		{
		theta = 0.5*urand(-M_PI,M_PI);
		phi = urand(0,2*M_PI);
		a = urand(RAD/5.,RAD);

		solatom[i][0][0] = a*cos(theta)*cos(phi);
		solatom[i][0][1] = a*cos(theta)*sin(phi);
		solatom[i][0][2] = a*sin(theta);

                for(n=1;n<3;++n)
                        for(k=0;k<3;++k)
				solatom[i][n][k] = solatom[i][0][k] + urand(-1,1);
		}

	for(i=0;i<ProLen;++i)
		for (j=0;j<lenRP[i];++j)
			{
			a = RAD*urand(0.,1.);
			theta = 0.5*urand(-M_PI,M_PI);
			phi = urand(0,2*M_PI);
						
			atom[i][j][0] = a*cos(theta)*cos(phi);
			atom[i][j][1] = a*cos(theta)*sin(phi);
			atom[i][j][2] = a*sin(theta);
			}

	for(i=0; i<ProLen; ++i)
		for(j=0; j<ProLen; ++j)
			for(m=0; m<lenRP[i]; ++m)
				for(n=0; n<lenRP[j]; ++n)
					if (CONH[i][j][m][n] != 0)
						for(c=0; c<2; ++c)
							{
							TOTCONST++;

							etaEB[i][j][m][n][c] = 1.;
							if (c == 1)
								{
								if (INTRP[i][j][m][n][0] > 0.)
									EB[i][j][m][n][c] = urand(0.,10.);
								else
									EB[i][j][m][n][c] = -1.*urand(0.,10.);
								}

							else
								EB[i][j][m][n][c] = urand(INTRP[i][j][m][n][1],1);
							}

        for (i=0;i<ProLen;++i)
                for(m=0;m<lenRP[i];++m)
                        for (j=0;j<NumSol;++j)
                                for(n=0; n<3; ++n)
					for(c=0; c<2; ++c)
						{
						TOTCONST++;
						etaES[i][m][j][n][c] = 1.;
						ES[i][m][j][n][c] = urand(-1.,1.);
						}

        for(i=0;i<NumSol;++i)
                for(j=i+1;j<NumSol;++j)
                        for(m=0;m<3;++m)
                                for(n=0;n<3;++n)
                                        for(c=0;c<2;++c)
						{
						TOTCONST++;
						etaEO[i][j][m][n][c] = 1.;
						EO[i][j][m][n][c] = urand(-1.,1.);
						}

        for(i=0;i<NumSol;++i)
		etaOS[i] = 1.;

	TOTCONST += NumSol;

        for (i=0;i<ProLen-1;++i)
		etaCT[i] = 1.;

	TOTCONST += ProLen-1;

        for (i=0;i<ProLen;++i)
		for (j=0;j<NumSegsRP[i];++j)
			{
			TOTCONST++;
			etaSC[i][j] = 1.;
			}

	for (i=0;i<NumCys;++i)
		for (j=0;j<NumCys;++j)
			etaSS[i][j] = 1.;

	TOTCONST += (NumCys-1)*NumCys;

        for (i=0;i<NumAD;++i)
		etaAD[i] = 1.;

	TOTCONST += NumAD;

        for (i=0;i<NumDN;++i)
                for (j=0;j<NumAC;++j)
                        etaHB[i][j] = 1.;

        TOTCONST += NumDN*NumAC;

        for(i=0;i<NumDH;++i)
		for(m=0;m<NumDHpar[i];++m)
			{
			TOTCONST++;
			etaDH[i][m] = 1.;
			ED[i][m] = urand(0,2.*ESCALE*DHpar[i][m][0]);;
			}

	for(j=0;j<NumBL;++j)
		{
		etaBL[j] = 1.;
		BE[j] = urand(0,10.*ESCALE);
		}

	TOTCONST += NumBL;
	
	changeVar(VB, VS, VO, OS, SC, CT, SS, AD, DH, BL, HB);

	//max_size = 36;
	//max_size = (ProLen-2)*(ProLen-2);
	val = malloc(TOTCONST*sizeof(double));
	val_index = malloc(TOTCONST*sizeof(int*[9]));
	}

void initSol(double lambda)
        {
        int i,j,k,m,n,c;
        //FILE *fp;
        double dis,rmin,eps,qmul,COM[3],theta,phi,LJRR,CPRR;
	TOTCONST = 0;

        //char file[50];

	double LJRS,CPRS,LJSS,CPSS;

	//double maxweight = 0.,tempweight;

	for(i=0;i<NumSolEx;++i)
		for(n=0;n<3;++n)
			for(k=0;k<3;++k)
				solatom[i][n][k] = ActualSolCoor[i][n][k]+lambda*urand(-1.,1.);


	for(i=0;i<ProLen;++i)
		for (j=0;j<lenRP[i];++j)
                        for(k=0;k<3;++k)
				{
                                atom[i][j][k] = ActualCoor[i][j][k]+lambda*urand(-1.,1.);
				COM[k] += atom[i][j][k];
				}
/*
	double max_dis = 0.;
	for(i=0;i<ProLen;++i)
		for (j=0;j<lenRP[i];++j)
			{
			dis = 0.;
                        for(k=0;k<3;++k)
				dis += sq(atom[i][j][k] - COM[k]);
			dis = sqrt(dis);
		
			if (dis > max_dis)
				max_dis = dis;
			}

	printf("%lf\n",max_dis);
*/

	if (NumSol > NumSolEx)
		for(i=NumSolEx;i<NumSol;++i)
			{
			theta = 0.5*urand(-M_PI,M_PI);
			phi = urand(0,2*M_PI);

			solatom[i][0][0] = COM[0]+RAD*cos(theta)*cos(phi);
			solatom[i][0][1] = COM[1]+RAD*cos(theta)*sin(phi);
			solatom[i][0][2] = COM[2]+RAD*sin(theta);

			for(n=1;n<3;++n)
				for(k=0;k<3;++k)
					solatom[i][n][k] = solatom[i][0][k] + urand(-1,1);
			}

	double LJmax = 4;

	LJRR = 0.;
	CPRR = 0.;

	for(i=0; i<ProLen; ++i)
		for(j=i; j<ProLen; ++j)
			for(m=0; m<lenRP[i]; ++m)
				for(n=0; n<lenRP[j]; ++n)
					if (CONH[i][j][m][n] != 0)
						{
						TOTCONST += 2;
						for(c=0; c<2; ++c)
							etaEB[i][j][m][n][c] = 1.;

						dis = 0.;
						for(k=0; k<3; ++k)
							dis += sq(atom[j][n][k] - atom[i][m][k]);
						dis = sqrt(dis);

						qmul = INTRP[i][j][m][n][0];
						eps = INTRP[i][j][m][n][1];
						rmin = INTRP[i][j][m][n][2];

						EB[i][j][m][n][0] = ESCALE*eps*(pow((rmin/dis), 12) - 2.*pow((rmin/dis), 6));
						if (EB[i][j][m][n][0] < 0)
							EB[i][j][m][n][0] = 0.;

						EB[i][j][m][n][1] = ESCALE*331.21*qmul/dis;
						
						LJRR += EB[i][j][m][n][0];
						CPRR += EB[i][j][m][n][1];

						EB[i][j][m][n][0] += lambda*urand(-1,1);
						EB[i][j][m][n][1] += lambda*urand(-1,1);

						if (EB[i][j][m][n][0] > LJmax)
							{
							//printf("%d %d %d %d %lf\n",i,j,m,n,EB[i][j][m][n][0]);
							EB[i][j][m][n][0] = LJmax;
							}
						}

	LJRS = 0.;
	CPRS = 0.;

        for (i=0;i<ProLen;++i)
                for(m=0;m<lenRP[i];++m)
                        for (j=0;j<NumSol;++j)
                                for(n=0; n<3; ++n)
					{
					TOTCONST += 2;

					for(c=0; c<2; ++c)
						etaES[i][m][j][n][c] = 1.;

					dis = 0.;
					for(k=0; k<3; ++k)
						dis += sq(solatom[j][n][k] - atom[i][m][k]);
					dis = sqrt(dis);

					rmin = (parRP[i][m][2] + parSol[n][2]);
					eps = sqrt(parRP[i][m][1]*parSol[n][1]);

					ES[i][m][j][n][0] = eps*(pow((rmin/dis), 12) - 2.*pow((rmin/dis), 6));
					ES[i][m][j][n][1] = 331.21*parRP[i][m][0]*parSol[n][0]/dis;

					if (j < NumSolEx)
						{
						LJRS += ES[i][m][j][n][0];
						CPRS += ES[i][m][j][n][1];
						}

					ES[i][m][j][n][0] += lambda*urand(-1,1);
					ES[i][m][j][n][1] += lambda*urand(-1,1);

					if (ES[i][m][j][n][0] > 0.)
						ES[i][m][j][n][0] = 0.;
					}

	LJSS = 0.;
	CPSS = 0.;

        for(i=0;i<NumSol;++i)
                for(j=i+1;j<NumSol;++j)
                        for(m=0;m<3;++m)
                                for(n=0;n<3;++n)
					{
					TOTCONST += 2;

                                        for(c=0;c<2;++c)
						etaEO[i][j][m][n][c] = 1.;

					dis = 0.;
					for(k=0; k<3; ++k)
						dis += sq(solatom[j][n][k] - solatom[i][m][k]);
					dis = sqrt(dis);

					rmin = (parSol[m][2] + parSol[n][2]);
					eps = sqrt(parSol[m][1]*parSol[n][1]);

					EO[i][j][m][n][0] = eps*(pow((rmin/dis), 12) - 2.*pow((rmin/dis), 6));
					EO[i][j][m][n][1] = 331.21*parSol[m][0]*parSol[n][0]/dis;

					if ((i < NumSolEx) && (j < NumSolEx))
						{
						LJSS += EO[i][j][m][n][0];
						CPSS += EO[i][j][m][n][1];
						}

					EO[i][j][m][n][0] += lambda*urand(-1,1);
					EO[i][j][m][n][1] += lambda*urand(-1,1);

					if (EO[i][j][m][n][0] > 0.)
						EO[i][j][m][n][0] = 0.;
					}

	double tmp[4][3],dhd,Edhd = 0.;

	for(j=0;j<NumDH;++j)
		{
		for(i=0;i<4;++i)
			for(k=0; k<3; ++k)
				tmp[i][k] = atom[A2R[DHAtoms[j][i]][0]][A2R[DHAtoms[j][i]][1]][k];

		dhd = calc_dihed(tmp);
	
		for(m=0;m<NumDHpar[j];++m)
			{
			ED[j][m] = ESCALE*DHpar[j][m][0]*(1 + cos(DHpar[j][m][1]*dhd - DHpar[j][m][2]));
			Edhd += ED[j][m];
			//printf("%d %d %lf %lf\n",j,m,dhd,ED[j][m]);
			etaDH[j][m] = 1.;
			TOTCONST++;
			}
	
		}

	double Ebond = 0.;

	for(j=0;j<NumBL;++j)
		{
		dis = 0.;

		for(k=0;k<3;++k)
			dis += sq(atom[A2R[BLAtoms[j][0]][0]][A2R[BLAtoms[j][0]][1]][k] - atom[A2R[BLAtoms[j][1]][0]][A2R[BLAtoms[j][1]][1]][k]);
		dis = sqrt(dis);

		BE[j] = ESCALE*0.5*BLpar[j][0]*sq(dis - BLpar[j][1]);

		Ebond += BE[j];

		etaBL[j] = 1.;
		TOTCONST++;
		}

        for(i=0;i<NumSol;++i)
		{
		etaOS[i] = 1.;
		TOTCONST++;
		}

        for (i=0;i<ProLen-1;++i)
		{
		etaCT[i] = 1.;
		TOTCONST++;
		}

        for (i=0;i<ProLen;++i)
		for (j=0;j<NumSegsRP[i];++j)
			{
			etaSC[i][j] = 1.;
			TOTCONST++;
			}

	for (i=0;i<NumCys;++i)
		for (j=0;j<NumCys;++j)
			{
			etaSS[i][j] = 1.;
			TOTCONST++;
			}

        for (i=0;i<NumAD;++i)
		etaAD[i] = 1.;

	TOTCONST += NumAD;

        for (i=0;i<NumDN;++i)
                for (j=0;j<NumAC;++j)
                        etaHB[i][j] = 1.;

        TOTCONST += NumDN*NumAC;

	changeVar(VB, VS, VO, OS, SC, CT, SS, AD, DH, BL, HB);

	printf("\nTotal Energy: %lf %lf %lf %lf %lf\n",CPRR,LJRR,Edhd,Ebond,CPRR+LJRR+Edhd+Ebond);
	printf("\nTotal Energy: %lf %lf %lf %lf\n",CPRS,LJRS,CPSS,LJSS);

	val = malloc(TOTCONST*sizeof(double));
	val_index = malloc(TOTCONST*sizeof(int*[9]));
	}

int main(int argc,char* argv[])
        {
        //char *id,MotifFile[20],MotifFileBB[20],MotifFileSol[20],MotifFileRot[20],buf[1024];
        //int c,iterstride,i,j,k,m,n,z,p,q,r,seed;
	int c,iterstride,seed;
	double stoperr;
	char *id;

        FILE *fp;

        if(argc==20)
                {
                ProID = argv[1];
                NumSol = atoi(argv[2]);
                EGS = atof(argv[3]);
                CHstep = atoi(argv[4]);
                id = argv[5];
                beta = atof(argv[6]);
                maxiter = atoi(argv[7]);
                iterstride = atoi(argv[8]);
                stoperr = atof(argv[9]);
                epsilon = atof(argv[10]);
                seed = atoi(argv[11]);
                RAD = atof(argv[12]);
                HB1 = atoi(argv[13]);
                HB2 = atoi(argv[14]);
                epsRR = atof(argv[15]);
                epsRS = atof(argv[16]);
                epsSS = atof(argv[17]);
                ESCALE = atof(argv[18]);
                BRAD = atof(argv[19]);
                }
        else
                {
                printf("expected more/less arguments \n");
                return 1;
                }

        snprintf(errfile, 50, "%s.err",id);
        snprintf(statsfile, 50, "%s.stats",id);
        snprintf(solfile, 50, "%s.sol",id);
        snprintf(etafile, 50, "%s.eta",id);
        snprintf(Enerfile, 50, "%s.enrg",id);
        snprintf(initfile, 50, "%s.init",id);

        snprintf(ProFile, 50, "%s.seq",ProID);

        AAmapping();
        SCAtomMapping();

        //if(!getContactPairs())
        //        return 1;

        if(!getRigidMotifs())
                return 1;

        if(!getprotein(ProFile))
                return 1;

        if(!getDHinfo())
                return 1;

        if(!getMotifsSC())
                return 1;

        //if(!getMotifsBB())
        //        return 1;

        fp=fopen(errfile,"w");
        fclose(fp);

        fp=fopen(solfile,"w");
        fclose(fp);

        fp=fopen(Enerfile,"w");
        fclose(fp);

        fp=fopen(etafile,"w");
        fclose(fp);

        fp=fopen(statsfile,"w");
        for(c=0;c<argc;++c)
                fprintf(fp,"%s ",argv[c]);
        fprintf(fp,"\n\n");

        fprintf(fp,"Protein Size: %d\n\n",ProLen);
        fprintf(fp,"Residue Types:\n");

        for(c=0;c<ProLen;++c)
                fprintf(fp,"%d ",RT[c]);
        fprintf(fp,"\n\n");

        fclose(fp);

        makevars();
        srand(seed);


        //initSol(2.);
        init();
        //init2();


	//projA(VB, VS, VO, EB, ES, EO, OS, SC, CT, SS, AD, DH, ED, BL, BE, HB);
	//reflect(VBA, VSA, VOA, EBA, ESA, EOA, OSA, SCA, CTA, SSA, ADA, DHA, EDA, BLA, BEA, HBA);
	//projB(VBR, VSR, VOR, EBR, ESR, EOR, OSR, SCR, CTR, SSR, ADR, DHR, EDR, BLR, BER, HBR);
        //RRR();
        //iter=solve(maxiter,iterstride,stoperr);

        double deltat, iterpersec;
        clock_t start;
        start=clock();
        iter=solve(maxiter,iterstride,stoperr);

        deltat=((double)(clock()-start))/CLOCKS_PER_SEC;

        fp=fopen(statsfile,"a");
        if(iter)
                {
                iterpersec=iter/deltat;
                fprintf(fp,"\nTime elapsed: %10.2lf \nNumber of iterations: %d \n",deltat,iter);
                }
        else
                {
                iterpersec=maxiter/deltat;
                fprintf(fp,"\nTime elapsed: %10.2lf \n",deltat);
                }

        fprintf(fp,"iterations/sec:%10.2lf\n\n",iterpersec);

        fprintf(fp,"\n\nWeightsAD:\n\n");

	int i,m;
        for(i=0;i<NumAD;++i)
		{
		if (WeightsAD[i] < 0.5)
			{
			for(m=0;m<4;++m)
				fprintf(fp,"%d %d     ", ADBK[i][m][0],ADBK[i][m][1]);
			fprintf(fp,"\n");

			for(m=4;m<8;++m)
				fprintf(fp,"%d %d     ", ADBK[i][m][0],ADBK[i][m][1]);
			fprintf(fp,"\n");

			double tempweight;
			tempweight = RMSD(8,ADA[i],ADB[i]);
			fprintf(fp,"%lf \t %lf \t %lf  \t %lf\n\n",WeightsAD[i],etaAD[i],ADerr[i], tempweight);
			
			}
		}

        fclose(fp);

/*
	double tempAD[4][3],tempweight;

        for (int i=0;i<NumDN;++i)
                for (int j=0;j<NumAC;++j)
                        {
                        c = ADBK[0][i][0]-ADBK[1][j][0];

                        if (abs(c) < 3)
				continue;

			totiter = 10;
			motifproj(4,ADMotif,ADB[i][j],tempAD);
			tempweight = RMSD(4,ADB[i][j],tempAD);

			if (tempweight < 0.5)
				fprintf(fp,"%d %d %d %d %10.2lf\n",ADBK[0][i][0],ADBK[0][i][1],ADBK[1][j][0],ADBK[1][j][1],tempweight);
			}


        for (int i=0;i<NumDN;++i)
                for (int j=0;j<NumAC;++j)
                        {
			if (WeightsAD[i][j] < 0.2)
				fprintf(fp,"%d %d %d %d %10.2lf\n",ADBK[0][i][0],ADBK[0][i][1],ADBK[1][j][0],ADBK[1][j][1],WeightsAD[i][j]);
			}







	ActualCoor = malloc(ProLen*sizeof(double*));

        for(i=0;i<ProLen;++i)
		ActualCoor[i] = malloc(lenRP[i]*sizeof(double*[3]));

	char file[50];
        snprintf(file,"%s.coor",ProID);

	chdir("../Files4RRR");

        fp = fopen(file,"r");
        if(!fp)
                printf("coor_file not found\n");

	for(i=0;i<ProLen;++i)
		for (j=0;j<lenRP[i];++j)
                        for(k=0;k<3;++k)
                                fscanf(fp,"%lf%*[ ]%*[\n]",&ActualCoor[i][j][k]);

        fclose(fp);
	chdir("../SSHB");

	double allatom1[4*ProLen][3], allatom2[4*ProLen][3],rmsd,tempproj[4*ProLen][3];	
	double COMAC[3];

	int countatom = 0;

	for(n=0;n<3;++n)
		COMAC[n] = 0.;

        for (i=0;i<ProLen;++i)
                for(m=0;m<4;++m)
			{
			for(n=0;n<3;++n)
				{
				allatom1[countatom][n] = atom[i][m][n];
				allatom2[countatom][n] = ActualCoor[i][m][n];
				COMAC[n] += ActualCoor[i][m][n];
				}
			countatom++;
			}

	for(n=0;n<3;++n)
		COMAC[n] /= countatom;

        for (i=0;i<countatom;++i)
		for(n=0;n<3;++n)
			allatom2[i][n] -= COMAC[n]; 

	totiter = 200;

	motifproj(4*ProLen,allatom2,allatom1,tempproj);
			
	rmsd = RMSD(4*ProLen,tempproj,allatom1);

        fp=fopen(statsfile,"a");
        fprintf(fp,"\nRMSD:%10.2lf\n",rmsd);
*/

        return 0;
        }




