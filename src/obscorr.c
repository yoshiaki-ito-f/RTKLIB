/*------------------------------------------------------------------------------
* obscorr.c : observation data correction functions
*
* Apply carrier-phase and pseudorange corrections defined on satellite/
* frequency-band/elevation-azimuth grid cells to decoded observation data
* before re-encoding (used by the stream server conversion in str2str).
*
* correction file (CSV):
*     # comment ('#' or '%' at column 1), header lines are skipped
*     prn,band,el_lo,el_hi,az_lo,az_hi,ph_corr_cyc,pr_corr_m,count
*     G01,L1,10,15,0,5,0.0123,0.45,1523
*
*   prn  : RINEX satellite id (G01,R07,E12,J01,C23,...)
*   band : L1|L2 (or 1|2), mapped to frequency index 0/1 of obsd_t
*   el/az: grid cell bounds (deg), corrections applied when
*          lo <= angle < hi (upper bound inclusive for the last cell)
*   ph   : carrier-phase correction (cycle), subtracted from L
*   pr   : pseudorange correction (m), subtracted from P
*   count: number of data in the grid cell (kept for quality check)
*
* correction file (binary, auto-detected by magic):
*     header : "OCB1"(4) version(u16) endian-marker 0x1234(u16) nrec(u32)
*              reserved(16)                                     = 28 bytes
*     record : sys(u8 'G'|'R'|'E'|'J'|'C'...) prn(u8) band(u8 1|2) pad(u8)
*              el_lo,el_hi,az_lo,az_hi,ph,pr (f32 x 6) count(i32) = 32 bytes
*
* The correction file is reloaded when its modification time changes
* (checked every reload_int seconds). On any error the previous table is
* kept and observations are passed through uncorrected: this module never
* stops the stream.
*-----------------------------------------------------------------------------*/
#define _POSIX_C_SOURCE 200112L
#include <sys/stat.h>
#include "rtklib.h"

#define OBSCORR_MAXWARN 20     /* max warnings per type per minute */
#define OBSCORR_SETTLE  2      /* delay reload if file modified within (s) */

/* warning types for rate limiting */
#define WARN_POS  0            /* receiver position not available */
#define WARN_EPH  1            /* no ephemeris / satpos error */
#define WARN_CELL 2            /* no matching grid cell */

/* correction record (intermediate for file loading) --------------------------*/
typedef struct {
    int sat,band,seq;
    float el_lo,el_hi,az_lo,az_hi,ph,pr;
    int32_t cnt;
} corrrec_t;

/* free correction grid ------------------------------------------------------*/
static void free_grid(corrgrid_t *grid)
{
    if (!grid) return;
    free(grid->el_lo); free(grid->el_hi);
    free(grid->az_lo); free(grid->az_hi);
    free(grid->ph); free(grid->pr); free(grid->cnt);
    free(grid);
}
/* free correction table -----------------------------------------------------*/
static void free_tab(corrtab_t *tab)
{
    int i,j;
    if (!tab) return;
    for (i=0;i<MAXSAT;i++) for (j=0;j<OBSCORR_NBAND;j++) {
        free_grid(tab->grid[i][j]);
    }
    free(tab);
}
/* parse band field ("L1"|"L2"|"1"|"2") --------------------------------------*/
static int parse_band(const char *str)
{
    if (!strcmp(str,"L1")||!strcmp(str,"l1")||!strcmp(str,"1")) return 0;
    if (!strcmp(str,"L2")||!strcmp(str,"l2")||!strcmp(str,"2")) return 1;
    return -1;
}
/* add correction record -----------------------------------------------------*/
static int add_rec(corrrec_t **recs, int *n, int *nmax, const corrrec_t *rec)
{
    corrrec_t *p;

    if (*n>=*nmax) {
        *nmax=*nmax<=0?4096:*nmax*2;
        if (!(p=(corrrec_t *)realloc(*recs,sizeof(corrrec_t)*(*nmax)))) {
            free(*recs); *recs=NULL; *n=*nmax=0;
            return 0;
        }
        *recs=p;
    }
    (*recs)[*n]=*rec;
    (*recs)[*n].seq=*n;
    (*n)++;
    return 1;
}
/* load correction records from CSV ------------------------------------------*/
static int load_csv(FILE *fp, const char *file, corrrec_t **recs, int *n,
                    int *nmax)
{
    corrrec_t rec={0};
    char buff[256],prn[16],band[16];
    int line=0,cnt;

    while (fgets(buff,sizeof(buff),fp)) {
        line++;
        if (buff[0]=='#'||buff[0]=='%'||buff[0]=='\r'||buff[0]=='\n') continue;
        if (sscanf(buff," %15[^, \t\r\n],%15[^, \t\r\n],%f,%f,%f,%f,%f,%f,%d",
                   prn,band,&rec.el_lo,&rec.el_hi,&rec.az_lo,&rec.az_hi,
                   &rec.ph,&rec.pr,&cnt)<9) {
            /* skip header lines silently */
            if (strstr(buff,"prn")||strstr(buff,"PRN")) continue;
            trace(2,"obscorr: invalid record %s:%d\n",file,line);
            continue;
        }
        if (!(rec.sat=satid2no(prn))) {
            trace(2,"obscorr: invalid satellite id %s:%d prn=%s\n",file,line,
                  prn);
            continue;
        }
        if ((rec.band=parse_band(band))<0) {
            trace(2,"obscorr: invalid band %s:%d band=%s\n",file,line,band);
            continue;
        }
        rec.cnt=(int32_t)cnt;
        if (!add_rec(recs,n,nmax,&rec)) return 0;
    }
    return 1;
}
/* load correction records from binary file ----------------------------------*/
static int load_bin(FILE *fp, const char *file, corrrec_t **recs, int *n,
                    int *nmax)
{
    corrrec_t rec={0};
    uint8_t head[28],buff[32];
    uint32_t i,nrec;
    char id[8];

    if (fread(head,1,28,fp)<28) {
        trace(2,"obscorr: binary header read error %s\n",file);
        return 0;
    }
    if ((uint16_t)(head[6]|head[7]<<8)!=0x1234) {
        trace(2,"obscorr: binary endian mismatch %s\n",file);
        return 0;
    }
    nrec=head[8]|head[9]<<8|head[10]<<16|(uint32_t)head[11]<<24;

    for (i=0;i<nrec;i++) {
        if (fread(buff,1,32,fp)<32) {
            trace(2,"obscorr: binary record read error %s n=%u/%u\n",file,i,
                  nrec);
            break;
        }
        sprintf(id,"%c%02d",(char)buff[0],buff[1]);
        if (!(rec.sat=satid2no(id))) {
            trace(2,"obscorr: invalid satellite %s rec=%u\n",file,i);
            continue;
        }
        rec.band=buff[2]-1;
        if (rec.band<0||rec.band>=OBSCORR_NBAND) {
            trace(2,"obscorr: invalid band %s rec=%u\n",file,i);
            continue;
        }
        memcpy(&rec.el_lo,buff+ 4,4); memcpy(&rec.el_hi,buff+ 8,4);
        memcpy(&rec.az_lo,buff+12,4); memcpy(&rec.az_hi,buff+16,4);
        memcpy(&rec.ph   ,buff+20,4); memcpy(&rec.pr   ,buff+24,4);
        memcpy(&rec.cnt  ,buff+28,4);
        if (!add_rec(recs,n,nmax,&rec)) return 0;
    }
    return 1;
}
/* compare records by (sat,band,seq) ------------------------------------------*/
static int cmp_rec(const void *p1, const void *p2)
{
    const corrrec_t *r1=(const corrrec_t *)p1,*r2=(const corrrec_t *)p2;
    if (r1->sat !=r2->sat ) return r1->sat -r2->sat;
    if (r1->band!=r2->band) return r1->band-r2->band;
    return r1->seq-r2->seq;
}
/* collect sorted unique lower bounds and their upper bounds ------------------*/
static int uniq_bounds(const corrrec_t *recs, int n, int az, float **lo,
                       float **hi)
{
    int i,j,k,m=0;
    float lo_i,hi_i;

    if (!(*lo=(float *)malloc(sizeof(float)*n))||
        !(*hi=(float *)malloc(sizeof(float)*n))) {
        free(*lo); *lo=NULL;
        return 0;
    }
    for (i=0;i<n;i++) {
        lo_i=az?recs[i].az_lo:recs[i].el_lo;
        hi_i=az?recs[i].az_hi:recs[i].el_hi;
        for (j=0;j<m;j++) if ((*lo)[j]>=lo_i) break;
        if (j<m&&(*lo)[j]==lo_i) {
            (*hi)[j]=hi_i;                    /* last record wins */
            continue;
        }
        for (k=m;k>j;k--) {
            (*lo)[k]=(*lo)[k-1]; (*hi)[k]=(*hi)[k-1];
        }
        (*lo)[j]=lo_i; (*hi)[j]=hi_i; m++;
    }
    return m;
}
/* find index of sorted lower bounds not exceeding value ----------------------*/
static int find_bound(const float *lo, int n, double v)
{
    int i=0,j=n-1,k;
    if (n<=0||v<lo[0]) return -1;
    while (i<j) {
        k=(i+j+1)/2;
        if (lo[k]<=v) i=k; else j=k-1;
    }
    return i;
}
/* build correction grid for one satellite/band -------------------------------*/
static corrgrid_t *build_grid(const corrrec_t *recs, int n)
{
    corrgrid_t *grid;
    int i,ie,ia,k,ncell;

    if (n<=0) return NULL;

    if (!(grid=(corrgrid_t *)calloc(1,sizeof(corrgrid_t)))) return NULL;

    if ((grid->nel=uniq_bounds(recs,n,0,&grid->el_lo,&grid->el_hi))<=0||
        (grid->naz=uniq_bounds(recs,n,1,&grid->az_lo,&grid->az_hi))<=0) {
        free_grid(grid);
        return NULL;
    }
    ncell=grid->nel*grid->naz;
    if (!(grid->ph =(float *)malloc(sizeof(float)*ncell))||
        !(grid->pr =(float *)malloc(sizeof(float)*ncell))||
        !(grid->cnt=(int32_t *)malloc(sizeof(int32_t)*ncell))) {
        free_grid(grid);
        return NULL;
    }
    for (k=0;k<ncell;k++) {
        grid->ph[k]=grid->pr[k]=(float)NAN;
        grid->cnt[k]=0;
    }
    for (i=0;i<n;i++) {
        if ((ie=find_bound(grid->el_lo,grid->nel,recs[i].el_lo))<0) continue;
        if ((ia=find_bound(grid->az_lo,grid->naz,recs[i].az_lo))<0) continue;
        k=ie*grid->naz+ia;
        grid->ph [k]=recs[i].ph;
        grid->pr [k]=recs[i].pr;
        grid->cnt[k]=recs[i].cnt;
    }
    /* test uniform grid for O(1) lookup */
    grid->el0=grid->el_lo[0];
    grid->del=grid->el_hi[0]-grid->el_lo[0];
    grid->az0=grid->az_lo[0];
    grid->daz=grid->az_hi[0]-grid->az_lo[0];
    grid->uniform=grid->del>0.0f&&grid->daz>0.0f;
    for (i=0;i<grid->nel&&grid->uniform;i++) {
        if (fabsf(grid->el_lo[i]-(grid->el0+i*grid->del))>1e-3f||
            fabsf(grid->el_hi[i]-grid->el_lo[i]-grid->del)>1e-3f) {
            grid->uniform=0;
        }
    }
    for (i=0;i<grid->naz&&grid->uniform;i++) {
        if (fabsf(grid->az_lo[i]-(grid->az0+i*grid->daz))>1e-3f||
            fabsf(grid->az_hi[i]-grid->az_lo[i]-grid->daz)>1e-3f) {
            grid->uniform=0;
        }
    }
    return grid;
}
/* load correction table from file --------------------------------------------*/
static corrtab_t *load_tab(const char *file)
{
    FILE *fp;
    corrtab_t *tab;
    corrrec_t *recs=NULL;
    uint8_t magic[4];
    int i,j,n=0,nmax=0,stat,ngrid=0;

    if (!(fp=fopen(file,"rb"))) {
        trace(2,"obscorr: correction file open error %s\n",file);
        return NULL;
    }
    if (fread(magic,1,4,fp)==4&&!memcmp(magic,"OCB1",4)) {
        rewind(fp);
        stat=load_bin(fp,file,&recs,&n,&nmax);
    }
    else {
        rewind(fp);
        stat=load_csv(fp,file,&recs,&n,&nmax);
    }
    fclose(fp);

    if (!stat||n<=0) {
        trace(2,"obscorr: no correction data %s\n",file);
        free(recs);
        return NULL;
    }
    if (!(tab=(corrtab_t *)calloc(1,sizeof(corrtab_t)))) {
        free(recs);
        return NULL;
    }
    /* sort by (sat,band), keeping file order within a cell group */
    qsort(recs,n,sizeof(corrrec_t),cmp_rec);

    for (i=0;i<n;i=j) {
        for (j=i+1;j<n;j++) {
            if (recs[j].sat!=recs[i].sat||recs[j].band!=recs[i].band) break;
        }
        tab->grid[recs[i].sat-1][recs[i].band]=build_grid(recs+i,j-i);
        if (tab->grid[recs[i].sat-1][recs[i].band]) ngrid++;
    }
    free(recs);

    if (ngrid<=0) {
        trace(2,"obscorr: no valid correction grid %s\n",file);
        free_tab(tab);
        return NULL;
    }
    tab->loadtime=timeget();
    trace(3,"obscorr: correction file loaded %s nrec=%d ngrid=%d\n",file,n,
          ngrid);
    return tab;
}
/* load rinex navigation file --------------------------------------------------*/
static int load_navf(obscorr_t *oc)
{
    nav_t nav={0};

    if (readrnx(oc->navfile,0,"",NULL,&nav,NULL)<1||(nav.n<=0&&nav.ng<=0)) {
        trace(2,"obscorr: nav file read error %s\n",oc->navfile);
        freenav(&nav,0xFF);
        return 0;
    }
    uniqnav(&nav);
    freenav(&oc->navf,0xFF);
    oc->navf=nav;
    trace(3,"obscorr: nav file loaded %s n=%d ng=%d\n",oc->navfile,nav.n,
          nav.ng);
    return 1;
}
/* test file update and settle time --------------------------------------------*/
static int file_updated(const char *file, time_t *mtime)
{
    struct stat st;

    if (stat(file,&st)) return 0;
    if (st.st_mtime==*mtime) return 0;

    /* delay reload if file may be being written */
    if ((double)(time(NULL)-st.st_mtime)<OBSCORR_SETTLE) return 0;

    *mtime=st.st_mtime;
    return 1;
}
/* rate limited warning test ---------------------------------------------------*/
static int warn_ok(obscorr_t *oc, int type)
{
    uint32_t tick=tickget();

    if ((int)(tick-oc->tick_warn[type])>=60000) {
        if (oc->nwarn[type]>OBSCORR_MAXWARN) {
            trace(3,"obscorr: %u warnings suppressed (type=%d)\n",
                  oc->nwarn[type]-OBSCORR_MAXWARN,type);
        }
        oc->tick_warn[type]=tick;
        oc->nwarn[type]=0;
    }
    return ++oc->nwarn[type]<=OBSCORR_MAXWARN;
}
/* initialize observation data correction ---------------------------------------
* initialize observation data correction control
* args   : obscorr_t *oc    IO  observation correction control
*          char   *file     I   correction file path
*          char   *navfile  I   rinex navigation file path ("":no file)
*          int    reload_int I  reload check interval (s) (<=0: default 30)
* return : status (1:ok,0:correction file not loaded (pass-through until
*          the file becomes readable))
* notes  : never aborts: on load error the stream runs uncorrected and the
*          file is retried by obscorr_reload()
*-----------------------------------------------------------------------------*/
extern int obscorr_init(obscorr_t *oc, const char *file, const char *navfile,
                        int reload_int)
{
    struct stat st;
    int i;

    trace(3,"obscorr_init: file=%s navfile=%s reload_int=%d\n",file,navfile,
          reload_int);

    memset(oc,0,sizeof(obscorr_t));
    sprintf(oc->file,"%.*s",(int)sizeof(oc->file)-1,file);
    sprintf(oc->navfile,"%.*s",(int)sizeof(oc->navfile)-1,navfile);
    oc->reload_int=reload_int<=0?30:reload_int;
    oc->tick_reload=tickget();
    for (i=0;i<4;i++) oc->tick_warn[i]=oc->tick_reload;

    if (!stat(oc->file,&st)) oc->mtime=st.st_mtime;
    if ((oc->tab=load_tab(oc->file))) oc->stat=1; else oc->stat=-1;

    if (*oc->navfile) {
        if (!stat(oc->navfile,&st)) oc->navmtime=st.st_mtime;
        load_navf(oc);
    }
    return oc->stat==1;
}
/* free observation data correction ---------------------------------------------
* free observation correction control
* args   : obscorr_t *oc    IO  observation correction control
* return : none
*-----------------------------------------------------------------------------*/
extern void obscorr_free(obscorr_t *oc)
{
    free_tab(oc->tab);
    oc->tab=NULL;
    freenav(&oc->navf,0xFF);
}
/* reload correction data --------------------------------------------------------
* test update of correction/navigation files and reload them. to be called
* cyclically (interval test is done internally)
* args   : obscorr_t *oc    IO  observation correction control
* return : status (1:reloaded,0:not reloaded)
* notes  : on reload error the previous correction table is kept
*-----------------------------------------------------------------------------*/
extern int obscorr_reload(obscorr_t *oc)
{
    corrtab_t *tab;
    uint32_t tick=tickget();
    int stat=0;

    if ((int)(tick-oc->tick_reload)<oc->reload_int*1000) return 0;
    oc->tick_reload=tick;

    if (file_updated(oc->file,&oc->mtime)) {
        if ((tab=load_tab(oc->file))) {
            free_tab(oc->tab);
            oc->tab=tab;
            oc->stat=1;
            stat=1;
        }
        else if (oc->tab) {
            trace(2,"obscorr: reload error, previous table kept %s\n",
                  oc->file);
        }
    }
    if (*oc->navfile&&file_updated(oc->navfile,&oc->navmtime)) {
        stat|=load_navf(oc);
    }
    return stat;
}
/* lookup correction values --------------------------------------------------------
* lookup correction values of a satellite/band/el-az grid cell
* args   : corrtab_t *tab   I   correction table
*          int    sat       I   satellite number
*          int    band      I   frequency band (0:L1,1:L2)
*          double el,az     I   elevation/azimuth (deg)
*          double *ph       O   carrier-phase correction (cycle)
*          double *pr       O   pseudorange correction (m)
*          int    *cnt      O   number of data in the grid cell
* return : status (1:found,0:no correction)
*-----------------------------------------------------------------------------*/
extern int obscorr_lookup(const corrtab_t *tab, int sat, int band, double el,
                          double az, double *ph, double *pr, int *cnt)
{
    const corrgrid_t *g;
    int ie,ia,k;

    if (!tab||sat<1||sat>MAXSAT||band<0||band>=OBSCORR_NBAND) return 0;
    if (!(g=tab->grid[sat-1][band])) return 0;

    az=fmod(az,360.0);
    if (az<0.0) az+=360.0;

    if (g->uniform) {
        ie=(int)floor((el-g->el0)/g->del);
        ia=(int)floor((az-g->az0)/g->daz);
        if (ie==g->nel&&el<=g->el_hi[g->nel-1]) ie=g->nel-1;
        if (ia==g->naz&&az<=g->az_hi[g->naz-1]) ia=g->naz-1;
        if (ie<0||ie>=g->nel||ia<0||ia>=g->naz) return 0;
    }
    else {
        if ((ie=find_bound(g->el_lo,g->nel,el))<0) return 0;
        if ((ia=find_bound(g->az_lo,g->naz,az))<0) return 0;
        if (el>=g->el_hi[ie]&&!(ie==g->nel-1&&el==g->el_hi[ie])) return 0;
        if (az>=g->az_hi[ia]&&!(ia==g->naz-1&&az==g->az_hi[ia])) return 0;
    }
    k=ie*g->naz+ia;
    if (isnan(g->ph[k])) return 0;

    *ph=g->ph[k];
    *pr=g->pr[k];
    *cnt=g->cnt[k];
    return 1;
}
/* apply corrections to observation data -------------------------------------------
* compute satellite elevation/azimuth of an epoch and subtract grid
* corrections from carrier-phase and pseudorange of decoded observation
* data in the stream converter output buffer
* args   : obscorr_t *oc    IO  observation correction control
*          strconv_t *conv  IO  stream converter
* return : none
* notes  : on any error observations are kept unchanged (pass-through)
*-----------------------------------------------------------------------------*/
extern void obscorr_apply(obscorr_t *oc, strconv_t *conv)
{
    obs_t *obs=&conv->out.obs;
    obsd_t *data;
    nav_t *nav;
    gtime_t time;
    double rr[3],pos[3],e[3],azel[2],el,az,ph,pr,prng;
    double rs[MAXOBS*6],dts[MAXOBS*2],var[MAXOBS];
    char id[8],tstr[40];
    int i,f,n=obs->n,svh[MAXOBS],cnt;

    if (!oc->tab||n<=0||n>MAXOBS) return;

    nav=conv->itype==STRFMT_RTCM2||conv->itype==STRFMT_RTCM3?
        &conv->rtcm.nav:&conv->raw.nav;

    /* receiver position (local station position or rtcm 1005/1006) */
    matcpy(rr,conv->out.sta.pos,3,1);
    if (norm(rr,3)<=0.0) matcpy(rr,conv->rtcm.sta.pos,3,1);
    if (norm(rr,3)<=0.0) {
        if (warn_ok(oc,WARN_POS)) {
            trace(3,"obscorr: receiver position not available, epoch skipped\n");
        }
        return;
    }
    ecef2pos(rr,pos);

    time=obs->data[0].time;
    satposs(time,obs->data,n,nav,EPHOPT_BRDC,rs,dts,var,svh);

    for (i=0;i<n;i++) {
        data=obs->data+i;
        satno2id(data->sat,id);
        time2str(data->time,tstr,0);

        /* retry with navigation file if no ephemeris in stream */
        if (svh[i]<0&&(oc->navf.n>0||oc->navf.ng>0)) {
            for (f=0,prng=0.0;f<NFREQ&&(prng=data->P[f])==0.0;f++) ;
            if (prng>0.0) {
                satpos(timeadd(data->time,-prng/CLIGHT),data->time,data->sat,
                       EPHOPT_BRDC,&oc->navf,rs+i*6,dts+i*2,var+i,svh+i);
            }
        }
        if (svh[i]<0||geodist(rs+i*6,rr,e)<=0.0) {
            if (warn_ok(oc,WARN_EPH)) {
                trace(3,"obscorr: no ephemeris sat=%s time=%s\n",id,tstr);
            }
            continue;
        }
        satazel(pos,e,azel);
        az=azel[0]*R2D;
        el=azel[1]*R2D;

        for (f=0;f<OBSCORR_NBAND&&f<NFREQ;f++) {
            if (data->code[f]==CODE_NONE) continue;
            if (data->L[f]==0.0&&data->P[f]==0.0) continue;

            if (!obscorr_lookup(oc->tab,data->sat,f,el,az,&ph,&pr,&cnt)) {
                if (warn_ok(oc,WARN_CELL)) {
                    trace(3,"obscorr: no correction sat=%s freq=L%d time=%s "
                          "el=%.1f az=%.1f\n",id,f+1,tstr,el,az);
                }
                continue;
            }
            if (data->L[f]!=0.0) data->L[f]-=ph;
            if (data->P[f]!=0.0) data->P[f]-=pr;

            trace(4,"obscorr: corrected sat=%s freq=L%d el=%.1f az=%.1f "
                  "ph=%.4f pr=%.3f cnt=%d\n",id,f+1,el,az,ph,pr,cnt);
        }
    }
}
