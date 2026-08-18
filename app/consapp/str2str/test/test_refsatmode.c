/*------------------------------------------------------------------------------
* test_refsatmode.c : unit tests of reference satellite selection (rtkpos.c)
*
* Covers: the opt.refsatmode selection rules used by ddres()/selrefsat()
*   0 elevation (legacy), 1 rover-SNR, 2 2nd elevation, 3 random,
*   4 az-el mask, 5 pinned satellite,
* the cycle-slip skip/fallback rule, SBAS exclusion, the rtk->refsat[][]
* logging used for the $REFSAT status record, QZSS/GPS DD group merging
* (opt.qzsmerge), and the config option wiring in options.c.
*
* note  : ddres() has internal linkage, so this file includes rtkpos.c
*         directly instead of linking against rtkpos.o (see makefile).
*
* usage : ./test_refsatmode
*         exit code 0 on all pass, 1 otherwise.
*-----------------------------------------------------------------------------*/
#include <stdio.h>
#include "rtkpos.c"

static int n_pass=0,n_fail=0;

#define CHECK(cond) do { \
    if (cond) n_pass++; \
    else {n_fail++; printf("FAIL %s:%d: %s\n",__FILE__,__LINE__,#cond);} \
} while (0)

/* base rtk_t/prcopt_t setup shared by ddres()-level tests ---------------------*/
static void init_test_opt(prcopt_t *opt)
{
    *opt=prcopt_default;
    opt->mode=PMODE_DGPS;      /* code-only DD residuals: skip amb/iono/tropo states */
    opt->nf=1;
    opt->navsys=SYS_GPS;
    opt->ionoopt=IONOOPT_OFF;
    opt->tropopt=TROPOPT_OFF;
    opt->glomodear=GLO_ARMODE_OFF;
    opt->baseline[0]=0.0;
    opt->maxinno[0]=opt->maxinno[1]=1000.0;
}

static void init_test_rtk(rtk_t *rtk, const prcopt_t *opt)
{
    rtkinit(rtk,opt);
    /* arbitrary but realistic ECEF position, rover==base (bl=0) */
    rtk->x[0]=-3961905.0; rtk->x[1]=3348994.0; rtk->x[2]=3697205.0;
    rtk->rb[0]=rtk->x[0]; rtk->rb[1]=rtk->x[1]; rtk->rb[2]=rtk->x[2];
}

/* find the reference satellite implied by vflg[] (upper byte of every entry) --*/
static int refsat_of(const int *vflg, int nv)
{
    int i,refsat;
    if (nv<=0) return -1;
    refsat=(vflg[0]>>16)&0xFF;
    for (i=1;i<nv;i++) {
        if (((vflg[i]>>16)&0xFF)!=refsat) return -1; /* inconsistent, shouldn't happen */
    }
    return refsat;
}

/* three-satellite scenario: highest-elevation sat != highest-SNR sat ----------*/
static void test_refsat_elevation_vs_snr(void)
{
    prcopt_t opt;
    rtk_t rtk;
    obsd_t obs[6]={{0}};
    double y[12]={0},e[18]={0},azel[12]={0},freq[6]={0};
    double v[8],R[64];
    int sat[3]={1,2,3},iu[3]={0,1,2},ir[3]={3,4,5},vflg[8];
    int j,nv,ref;

    printf("test_refsat_elevation_vs_snr\n");

    init_test_opt(&opt);
    /* sat A(1): highest elevation, low SNR
       sat B(2): highest SNR, lower elevation
       sat C(3): baseline (neither extreme) */
    for (j=0;j<3;j++) {
        obs[j].sat=obs[j+3].sat=(uint8_t)sat[j];
        freq[iu[j]]=FREQL1;
        azel[1+iu[j]*2]=30.0*D2R;   /* default elevation, overridden below */
        y[1+iu[j]*2]=100.0+j;       /* rover code residual */
        y[1+ir[j]*2]=98.0+j;        /* base code residual (DD stays small/constant) */
    }
    azel[1+iu[0]*2]=80.0*D2R; /* sat A: highest elevation */
    azel[1+iu[1]*2]=20.0*D2R; /* sat B: lowest elevation */
    azel[1+iu[2]*2]=50.0*D2R; /* sat C */

    /* opt.refsatmode==0: elevation-based -> sat A(1) must win */
    opt.refsatmode=0;
    init_test_rtk(&rtk,&opt);
    for (j=0;j<3;j++) rtk.ssat[sat[j]-1].sys=SYS_GPS;
    rtk.ssat[0].snr_rover[0]=20.0f; /* sat A: low SNR */
    rtk.ssat[1].snr_rover[0]=45.0f; /* sat B: high SNR */
    rtk.ssat[2].snr_rover[0]=35.0f; /* sat C */
    nv=ddres(&rtk,obs,0.0,rtk.x,rtk.P,sat,y,e,azel,freq,iu,ir,3,v,NULL,R,vflg);
    ref=refsat_of(vflg,nv);
    CHECK(nv==2);
    CHECK(ref==1);
    rtkfree(&rtk);

    /* opt.refsatmode==1: rover-SNR-based -> sat B(2) must win */
    opt.refsatmode=1;
    init_test_rtk(&rtk,&opt);
    for (j=0;j<3;j++) rtk.ssat[sat[j]-1].sys=SYS_GPS;
    rtk.ssat[0].snr_rover[0]=20.0f;
    rtk.ssat[1].snr_rover[0]=45.0f;
    rtk.ssat[2].snr_rover[0]=35.0f;
    nv=ddres(&rtk,obs,0.0,rtk.x,rtk.P,sat,y,e,azel,freq,iu,ir,3,v,NULL,R,vflg);
    ref=refsat_of(vflg,nv);
    CHECK(nv==2);
    CHECK(ref==2);
    rtkfree(&rtk);
}

/* default refsatmode (from rtkinit/prcopt_default) matches elevation mode -----*/
static void test_refsat_default_is_elevation(void)
{
    prcopt_t opt;
    rtk_t rtk;

    printf("test_refsat_default_is_elevation\n");

    init_test_opt(&opt);
    /* do NOT touch opt.refsatmode: verify prcopt_default's own value */
    CHECK(prcopt_default.refsatmode==0);
    init_test_rtk(&rtk,&opt);
    CHECK(rtk.opt.refsatmode==0);
    rtkfree(&rtk);
}

/* cycle-slip skip: a slipped, better-metric sat loses to a valid non-slip sat -*/
static void test_refsat_slip_skip(void)
{
    prcopt_t opt;
    rtk_t rtk;
    obsd_t obs[4]={{0}};
    double y[8]={0},e[12]={0},azel[8]={0},freq[4]={0};
    double v[8],R[64];
    int sat[2]={1,2},iu[2]={0,1},ir[2]={2,3},vflg[8];
    int j,nv,ref;

    printf("test_refsat_slip_skip\n");

    init_test_opt(&opt);
    opt.refsatmode=0; /* elevation-based; slip rule is mode-independent */
    init_test_rtk(&rtk,&opt);

    for (j=0;j<2;j++) {
        obs[j].sat=obs[j+2].sat=(uint8_t)sat[j];
        freq[iu[j]]=FREQL1;
        y[1+iu[j]*2]=100.0+j;
        y[1+ir[j]*2]=98.0+j;
        rtk.ssat[sat[j]-1].sys=SYS_GPS;
    }
    /* sat 1: lower elevation, no slip -> listed first in iteration order */
    azel[1+iu[0]*2]=20.0*D2R;
    /* sat 2: higher elevation, but slipped -> must be skipped in favor of sat 1 */
    azel[1+iu[1]*2]=80.0*D2R;
    rtk.ssat[1].slip[0]=LLI_SLIP;

    nv=ddres(&rtk,obs,0.0,rtk.x,rtk.P,sat,y,e,azel,freq,iu,ir,2,v,NULL,R,vflg);
    ref=refsat_of(vflg,nv);
    CHECK(nv==1);
    CHECK(ref==1); /* sat 1 wins despite worse elevation, because sat 2 is slipped */

    rtkfree(&rtk);
}

/* cycle-slip fallback: the only valid sat is used even though it is slipped ---*/
static void test_refsat_slip_fallback_only_sat(void)
{
    prcopt_t opt;
    rtk_t rtk;
    obsd_t obs[2]={{0}};
    double y[4]={0},e[6]={0},azel[4]={0},freq[2]={0};
    double v[4],R[16];
    int sat[1]={1},iu[1]={0},ir[1]={1},vflg[4];
    int nv;

    printf("test_refsat_slip_fallback_only_sat\n");

    init_test_opt(&opt);
    init_test_rtk(&rtk,&opt);

    obs[0].sat=obs[1].sat=1;
    freq[0]=FREQL1;
    y[1]=100.0; y[3]=98.0;
    azel[1]=45.0*D2R;
    rtk.ssat[0].sys=SYS_GPS;
    rtk.ssat[0].slip[0]=LLI_SLIP;

    /* single sat, nothing to double-difference against -> no residuals, but
       must not crash, and the reference-satellite search must accept the
       lone (slipped) candidate rather than bailing out early */
    nv=ddres(&rtk,obs,0.0,rtk.x,rtk.P,sat,y,e,azel,freq,iu,ir,1,v,NULL,R,vflg);
    CHECK(nv==0); /* no second sat to pair with -> zero DD residuals */

    rtkfree(&rtk);
}

/* four-satellite scenario shared by the mode 2/3/4/5 tests --------------------
 * elevations: A(1)=80, B(2)=20, C(3)=50, D(4)=65 deg
 * azimuths  : A due north (inside the mask), the rest due south (outside).
 * sys[] selects the constellation of each satellite.                          */
static void setup_four_sats(rtk_t *rtk, obsd_t *obs, double *y, double *azel,
                            double *freq, const int *sat, const int *iu,
                            const int *ir, const int *sys)
{
    static const double el[4]={80.0,20.0,50.0,65.0};
    int j;

    for (j=0;j<4;j++) {
        obs[j].sat=obs[j+4].sat=(uint8_t)sat[j];
        freq[iu[j]]=FREQL1;
        azel[iu[j]*2]=j==0?0.0:180.0*D2R;  /* sat A due north, others due south */
        azel[1+iu[j]*2]=el[j]*D2R;
        y[1+iu[j]*2]=100.0+j;
        y[1+ir[j]*2]=98.0+j;
        rtk->ssat[sat[j]-1].sys=sys[j];
        rtk->ssat[sat[j]-1].snr_rover[0]=40.0f;
    }
}

/* mode 2 (2nd highest elevation) and mode 4 (az-el mask) ----------------------*/
static void test_refsat_el2nd_and_mask(void)
{
    prcopt_t opt;
    rtk_t rtk;
    obsd_t obs[8]={{0}};
    double y[16]={0},e[24]={0},azel[16]={0},freq[8]={0};
    double v[16],R[256];
    int sat[4]={1,2,3,4},iu[4]={0,1,2,3},ir[4]={4,5,6,7},vflg[16];
    int sys[4]={SYS_GPS,SYS_GPS,SYS_GPS,SYS_GPS};
    int nv;

    printf("test_refsat_el2nd_and_mask\n");

    /* mode 0 picks the highest, sat A(1) at 80 deg, as the baseline */
    init_test_opt(&opt);
    opt.refsatmode=REFSAT_MAXEL;
    init_test_rtk(&rtk,&opt);
    setup_four_sats(&rtk,obs,y,azel,freq,sat,iu,ir,sys);
    nv=ddres(&rtk,obs,0.0,rtk.x,rtk.P,sat,y,e,azel,freq,iu,ir,4,v,NULL,R,vflg);
    CHECK(nv==3);
    CHECK(refsat_of(vflg,nv)==1);
    rtkfree(&rtk);

    /* mode 2: the 2nd highest is sat D(4) at 65 deg */
    opt.refsatmode=REFSAT_EL2ND;
    init_test_rtk(&rtk,&opt);
    setup_four_sats(&rtk,obs,y,azel,freq,sat,iu,ir,sys);
    nv=ddres(&rtk,obs,0.0,rtk.x,rtk.P,sat,y,e,azel,freq,iu,ir,4,v,NULL,R,vflg);
    CHECK(nv==3);
    CHECK(refsat_of(vflg,nv)==4);
    rtkfree(&rtk);

    /* mode 4: sat A(1) is due north above 35 deg, so it sits inside the mask
       and the highest remaining satellite, D(4), is used instead */
    CHECK(inrefsatmask(0.0,80.0*D2R)==1);        /* due north, high: blocked */
    CHECK(inrefsatmask(180.0*D2R,80.0*D2R)==0);  /* due south, high: clear */
    CHECK(inrefsatmask(0.0,20.0*D2R)==0);        /* due north, low: clear */
    opt.refsatmode=REFSAT_MASK;
    init_test_rtk(&rtk,&opt);
    setup_four_sats(&rtk,obs,y,azel,freq,sat,iu,ir,sys);
    nv=ddres(&rtk,obs,0.0,rtk.x,rtk.P,sat,y,e,azel,freq,iu,ir,4,v,NULL,R,vflg);
    CHECK(nv==3);
    CHECK(refsat_of(vflg,nv)==4);
    rtkfree(&rtk);

    /* mode 4 fallback: if every candidate is masked, the highest wins anyway
       so the sys group is never dropped just because of the selection rule */
    opt.refsatmode=REFSAT_MASK;
    init_test_rtk(&rtk,&opt);
    setup_four_sats(&rtk,obs,y,azel,freq,sat,iu,ir,sys);
    azel[iu[1]*2]=azel[iu[2]*2]=azel[iu[3]*2]=0.0;  /* all due north */
    azel[1+iu[1]*2]=40.0*D2R;                       /* and all above 35 deg */
    nv=ddres(&rtk,obs,0.0,rtk.x,rtk.P,sat,y,e,azel,freq,iu,ir,4,v,NULL,R,vflg);
    CHECK(nv==3);
    CHECK(refsat_of(vflg,nv)==1);                   /* fell back to max elevation */
    rtkfree(&rtk);
}

/* mode 3 (random): stays above the elevation threshold and is reproducible ----*/
static void test_refsat_random(void)
{
    prcopt_t opt;
    rtk_t rtk;
    obsd_t obs[8]={{0}};
    double y[16]={0},e[24]={0},azel[16]={0},freq[8]={0};
    double v[16],R[256];
    int sat[4]={1,2,3,4},iu[4]={0,1,2,3},ir[4]={4,5,6,7},vflg[16];
    int sys[4]={SYS_GPS,SYS_GPS,SYS_GPS,SYS_GPS};
    int i,nv,ref,first=-1,seen2=0,stable=1;

    printf("test_refsat_random\n");

    init_test_opt(&opt);
    opt.refsatmode=REFSAT_RANDOM;
    opt.refsatelmin=30.0*D2R;  /* excludes sat B(2) at 20 deg */

    /* same epoch time -> same seed -> same choice on every call.  ddres() runs
       several times per epoch, so this property is what keeps H and v
       consistent between the float update and the fixed-solution check. */
    for (i=0;i<8;i++) {
        init_test_rtk(&rtk,&opt);
        rtk.sol.time=gpst2time(2300,86400.0);
        setup_four_sats(&rtk,obs,y,azel,freq,sat,iu,ir,sys);
        nv=ddres(&rtk,obs,0.0,rtk.x,rtk.P,sat,y,e,azel,freq,iu,ir,4,v,NULL,R,vflg);
        ref=refsat_of(vflg,nv);
        CHECK(nv==3);
        CHECK(ref!=2);              /* sat B is below refsatelmin */
        if (first<0) first=ref; else if (ref!=first) stable=0;
        rtkfree(&rtk);
    }
    CHECK(stable);

    /* different epochs reseed, so the choice does move around */
    for (i=0;i<64&&!seen2;i++) {
        init_test_rtk(&rtk,&opt);
        rtk.sol.time=gpst2time(2300,86400.0+i*30.0);
        setup_four_sats(&rtk,obs,y,azel,freq,sat,iu,ir,sys);
        nv=ddres(&rtk,obs,0.0,rtk.x,rtk.P,sat,y,e,azel,freq,iu,ir,4,v,NULL,R,vflg);
        if (refsat_of(vflg,nv)!=first) seen2=1;
        rtkfree(&rtk);
    }
    CHECK(seen2);
}

/* mode 5 (pinned): the named satellite wins, with fallback when it is absent --*/
static void test_refsat_pinned(void)
{
    prcopt_t opt;
    rtk_t rtk;
    obsd_t obs[8]={{0}};
    double y[16]={0},e[24]={0},azel[16]={0},freq[8]={0};
    double v[16],R[256];
    int sat[4],iu[4]={0,1,2,3},ir[4]={4,5,6,7},vflg[16];
    int sys[4]={SYS_GPS,SYS_GPS,SYS_GPS,SYS_QZS};
    int nv,qzs=satid2no("J03");

    printf("test_refsat_pinned\n");
    CHECK(qzs>0);

    /* sat D is a QZSS satellite at 65 deg; the merged group puts it in m=0 */
    sat[0]=1; sat[1]=2; sat[2]=3; sat[3]=qzs;

    init_test_opt(&opt);
    opt.navsys=SYS_GPS|SYS_QZS;
    opt.qzsmerge=1;
    opt.refsatmode=REFSAT_PINNED;
    opt.refsatprn=qzs;
    opt.refsatelmin=10.0*D2R;
    init_test_rtk(&rtk,&opt);
    qzs_group=QZSGRP_GPS;   /* normally set by relpos(); ddres() is called directly here */
    setup_four_sats(&rtk,obs,y,azel,freq,sat,iu,ir,sys);
    nv=ddres(&rtk,obs,0.0,rtk.x,rtk.P,sat,y,e,azel,freq,iu,ir,4,v,NULL,R,vflg);
    CHECK(nv==3);
    CHECK(refsat_of(vflg,nv)==qzs); /* pinned sat wins over the higher sat A(1) */
    rtkfree(&rtk);

    /* the pinned satellite is not tracked: fall back to the highest QZSS */
    opt.refsatprn=satid2no("J07");
    init_test_rtk(&rtk,&opt);
    qzs_group=QZSGRP_GPS;
    setup_four_sats(&rtk,obs,y,azel,freq,sat,iu,ir,sys);
    nv=ddres(&rtk,obs,0.0,rtk.x,rtk.P,sat,y,e,azel,freq,iu,ir,4,v,NULL,R,vflg);
    CHECK(refsat_of(vflg,nv)==qzs);
    rtkfree(&rtk);

    /* no QZSS at all: fall back to the highest elevation, sat A(1) */
    sys[3]=SYS_GPS; sat[3]=4;
    init_test_rtk(&rtk,&opt);
    qzs_group=QZSGRP_GPS;
    setup_four_sats(&rtk,obs,y,azel,freq,sat,iu,ir,sys);
    nv=ddres(&rtk,obs,0.0,rtk.x,rtk.P,sat,y,e,azel,freq,iu,ir,4,v,NULL,R,vflg);
    CHECK(refsat_of(vflg,nv)==1);
    rtkfree(&rtk);

    qzs_group=QZSGRP_OWN;   /* restore for the remaining tests */
}

/* SBAS is never used as reference satellite, in any mode --------------------*/
static void test_refsat_sbas_excluded(void)
{
    prcopt_t opt;
    rtk_t rtk;
    obsd_t obs[4]={{0}};
    double y[8]={0},e[12]={0},azel[8]={0},freq[4]={0};
    double v[8],R[64];
    int sat[2],iu[2]={0,1},ir[2]={2,3},vflg[8];
    int j,nv,modes[3]={REFSAT_MAXEL,REFSAT_MAXSNR,REFSAT_MASK},sbs=satid2no("S20");

    printf("test_refsat_sbas_excluded\n");
    CHECK(sbs>0);
    sat[0]=1; sat[1]=sbs;

    for (j=0;j<3;j++) {
        init_test_opt(&opt);
        opt.navsys=SYS_GPS|SYS_SBS;
        opt.refsatmode=modes[j];
        init_test_rtk(&rtk,&opt);

        obs[0].sat=obs[2].sat=(uint8_t)sat[0];
        obs[1].sat=obs[3].sat=(uint8_t)sat[1];
        freq[0]=freq[1]=FREQL1;
        y[1]=100.0; y[3]=101.0; y[5]=98.0; y[7]=99.0;
        azel[1]=20.0*D2R;                  /* GPS sat: low */
        azel[2]=180.0*D2R; azel[3]=80.0*D2R; /* SBAS sat: high and unmasked */
        rtk.ssat[sat[0]-1].sys=SYS_GPS;  rtk.ssat[sat[0]-1].snr_rover[0]=30.0f;
        rtk.ssat[sat[1]-1].sys=SYS_SBS;  rtk.ssat[sat[1]-1].snr_rover[0]=50.0f;

        /* SBAS wins on both elevation and SNR, yet must not become reference */
        nv=ddres(&rtk,obs,0.0,rtk.x,rtk.P,sat,y,e,azel,freq,iu,ir,2,v,NULL,R,vflg);
        CHECK(nv==1);
        CHECK(refsat_of(vflg,nv)==sat[0]);
        rtkfree(&rtk);
    }
}

/* ddres() records the reference satellite in rtk->refsat[][] for $REFSAT ------*/
static void test_refsat_logged(void)
{
    prcopt_t opt;
    rtk_t rtk;
    obsd_t obs[8]={{0}};
    double y[16]={0},e[24]={0},azel[16]={0},freq[8]={0};
    double v[16],R[256];
    int sat[4]={1,2,3,4},iu[4]={0,1,2,3},ir[4]={4,5,6,7},vflg[16];
    int sys[4]={SYS_GPS,SYS_GPS,SYS_GPS,SYS_GPS};
    int m,f,nv,other=0;

    printf("test_refsat_logged\n");

    init_test_opt(&opt);
    opt.refsatmode=REFSAT_MAXEL;
    init_test_rtk(&rtk,&opt);
    setup_four_sats(&rtk,obs,y,azel,freq,sat,iu,ir,sys);
    /* the shared fixture uses a constant rover-base offset, which makes every
       DD exactly zero; offset sat 2 so its residual is distinguishable */
    y[1+iu[1]*2]+=0.5;
    nv=ddres(&rtk,obs,0.0,rtk.x,rtk.P,sat,y,e,azel,freq,iu,ir,4,v,NULL,R,vflg);
    CHECK(nv==3);

    /* PMODE_DGPS -> code only, so f runs over [nf,2*nf) = [1,2) with nf=1 */
    CHECK(rtk.refsat[0][1]==1);
    CHECK(rtk.refsat[0][1]==refsat_of(vflg,nv));
    /* every other sys group / frequency slot stays cleared */
    for (m=0;m<6;m++) for (f=0;f<NFREQ*2;f++) {
        if (m==0&&f==1) continue;
        if (rtk.refsat[m][f]!=0) other++;
    }
    CHECK(other==0);

    /* the reference satellite keeps resc/resp at zero, which is what makes the
       "$SAT with vsat==1 and resc==0" cross-check in the plan work */
    CHECK(rtk.ssat[0].resp[0]==0.0);
    CHECK(rtk.ssat[0].resc[0]==0.0);
    CHECK(rtk.ssat[1].resp[0]!=0.0);
    rtkfree(&rtk);
}

/* opt.qzsmerge moves QZSS between its own DD group and the GPS group ---------*/
static void test_qzsmerge_group(void)
{
    printf("test_qzsmerge_group\n");

    qzs_group=QZSGRP_OWN;   /* default: QZSS has its own group m=4 */
    CHECK(test_sys(SYS_QZS,4)==1);
    CHECK(test_sys(SYS_QZS,0)==0);
    CHECK(test_sys(SYS_GPS,0)==1);

    qzs_group=QZSGRP_GPS;   /* merged: QZSS joins the GPS group m=0 */
    CHECK(test_sys(SYS_QZS,0)==1);
    CHECK(test_sys(SYS_QZS,4)==0);
    /* the other constellations are unaffected either way */
    CHECK(test_sys(SYS_GPS,0)==1);
    CHECK(test_sys(SYS_SBS,0)==1);
    CHECK(test_sys(SYS_GLO,1)==1);
    CHECK(test_sys(SYS_GAL,2)==1);
    CHECK(test_sys(SYS_CMP,3)==1);
    CHECK(test_sys(SYS_IRN,5)==1);

    qzs_group=QZSGRP_OWN;
}

/* "pos1-refsatmode" config option round-trips through options.c ---------------*/
static void test_refsatmode_option_parsing(void)
{
    static const char *file="test_refsatmode.conf";
    prcopt_t popt;
    solopt_t sopt;
    filopt_t fopt;
    FILE *fp;

    printf("test_refsatmode_option_parsing\n");

    fp=fopen(file,"w");
    fputs("pos1-refsatmode  =1\n",fp);
    fclose(fp);

    resetsysopts();
    CHECK(loadopts(file,sysopts)!=0);
    getsysopts(&popt,&sopt,&fopt);
    CHECK(popt.refsatmode==1);

    /* the new modes and their parameters round-trip too */
    fp=fopen(file,"w");
    fputs("pos1-refsatmode  =5\n"
          "pos1-refsatelmin =12.5\n"
          "pos1-refsatprn   =J03\n"
          "pos1-qzsmerge    =1\n",fp);
    fclose(fp);

    resetsysopts();
    CHECK(loadopts(file,sysopts)!=0);
    getsysopts(&popt,&sopt,&fopt);
    CHECK(popt.refsatmode==REFSAT_PINNED);
    CHECK(fabs(popt.refsatelmin-12.5*D2R)<1E-9);
    CHECK(popt.refsatprn==satid2no("J03"));
    CHECK(popt.qzsmerge==1);

    /* an unset pinned satellite means "auto", not satellite 0 */
    fp=fopen(file,"w");
    fputs("pos1-refsatmode  =4\n",fp);
    fclose(fp);

    resetsysopts();
    CHECK(loadopts(file,sysopts)!=0);
    getsysopts(&popt,&sopt,&fopt);
    CHECK(popt.refsatprn==0);
    CHECK(popt.qzsmerge==0);

    remove(file);
}

int main(void)
{
    test_refsat_elevation_vs_snr();
    test_refsat_default_is_elevation();
    test_refsat_slip_skip();
    test_refsat_slip_fallback_only_sat();
    test_refsat_el2nd_and_mask();
    test_refsat_random();
    test_refsat_pinned();
    test_refsat_sbas_excluded();
    test_refsat_logged();
    test_qzsmerge_group();
    test_refsatmode_option_parsing();

    printf("%d passed, %d failed\n",n_pass,n_fail);
    return n_fail?1:0;
}
