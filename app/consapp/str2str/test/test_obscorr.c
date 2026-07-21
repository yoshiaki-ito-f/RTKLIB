/*------------------------------------------------------------------------------
* test_obscorr.c : unit tests of observation data correction (obscorr.c)
*
* Covers: CSV/binary loading, grid construction (uniform/non-uniform),
* cell lookup boundaries, azimuth normalization, duplicate/missing cells,
* invalid records, and file reload behavior.
*
* usage : ./test_obscorr
*         exit code 0 on all pass, 1 otherwise. temporary files are created
*         in the current directory (test_obscorr_*.csv/.bin) and removed.
*
* note  : the binary-format test assumes a little-endian host.
*-----------------------------------------------------------------------------*/
#include <stdio.h>
#include <utime.h>
#include "rtklib.h"

static int n_pass=0,n_fail=0;

#define CHECK(cond) do { \
    if (cond) n_pass++; \
    else {n_fail++; printf("FAIL %s:%d: %s\n",__FILE__,__LINE__,#cond);} \
} while (0)

#define CHECK_NEAR(a,b) CHECK(fabs((a)-(b))<1e-6)

static const char *CSVFILE="test_obscorr_grid.csv";
static const char *BINFILE="test_obscorr_grid.bin";

/* write text file -------------------------------------------------------------*/
static void write_file(const char *file, const char *text)
{
    FILE *fp=fopen(file,"w");
    fputs(text,fp);
    fclose(fp);
}
/* set file mtime into the past so reload settle/mtime checks pass -------------*/
static void backdate(const char *file, int sec)
{
    struct utimbuf ut;
    ut.actime=ut.modtime=time(NULL)-sec;
    utime(file,&ut);
}
/* test basic csv load and uniform-grid lookup ---------------------------------*/
static void test_csv_uniform(void)
{
    obscorr_t oc;
    double ph,pr,el,az;
    int ie,ia,cnt;
    char buff[4096],*p=buff;

    printf("test_csv_uniform\n");

    /* G01 L1: el cells [10,15),[15,20), az cells [0,5),[5,10),[10,15)
       ph=ie*10+ia+0.5, pr=ie*100+ia, cnt=ie*10+ia */
    p+=sprintf(p,"# comment line\n");
    p+=sprintf(p,"prn,band,el_lo,el_hi,az_lo,az_hi,ph,pr,count\n");
    for (ie=0;ie<2;ie++) for (ia=0;ia<3;ia++) {
        p+=sprintf(p,"G01,L1,%d,%d,%d,%d,%.1f,%.1f,%d\n",10+ie*5,15+ie*5,
                   ia*5,ia*5+5,ie*10+ia+0.5,(double)(ie*100+ia),ie*10+ia);
    }
    p+=sprintf(p,"%% another comment\n");
    p+=sprintf(p,"G01,L2,10,20,0,15,7.5,75.0,42\n");   /* single cell, band 2 */
    p+=sprintf(p,"R07,2,10,20,0,15,8.5,85.0,43\n");    /* numeric band */
    p+=sprintf(p,"J01,L1,10,20,0,15,1.5,15.0,44\n");   /* QZSS */
    p+=sprintf(p,"C05,L1,10,20,0,15,2.5,25.0,45\n");   /* BeiDou */
    p+=sprintf(p,"E11,L1,10,20,0,15,3.5,35.0,46\n");   /* Galileo */
    write_file(CSVFILE,buff);

    CHECK(obscorr_init(&oc,CSVFILE,"",0)==1);
    CHECK(oc.stat==1);
    CHECK(oc.tab!=NULL);
    CHECK(oc.reload_int==30);                 /* default applied */

    CHECK(oc.tab->grid[satid2no("G01")-1][0]!=NULL);
    CHECK(oc.tab->grid[satid2no("G01")-1][0]->uniform==1);
    CHECK(oc.tab->grid[satid2no("G01")-1][0]->nel==2);
    CHECK(oc.tab->grid[satid2no("G01")-1][0]->naz==3);

    /* interior lookup: all cells */
    for (ie=0;ie<2;ie++) for (ia=0;ia<3;ia++) {
        el=12.5+ie*5.0; az=2.5+ia*5.0;
        CHECK(obscorr_lookup(oc.tab,satid2no("G01"),0,el,az,&ph,&pr,&cnt)==1);
        CHECK_NEAR(ph,ie*10+ia+0.5);
        CHECK_NEAR(pr,ie*100+ia);
        CHECK(cnt==ie*10+ia);
    }
    /* boundaries: lower bound inclusive */
    CHECK(obscorr_lookup(oc.tab,satid2no("G01"),0,10.0,0.0,&ph,&pr,&cnt)==1);
    CHECK_NEAR(ph,0.5);
    /* inner upper bound belongs to the next cell */
    CHECK(obscorr_lookup(oc.tab,satid2no("G01"),0,15.0,5.0,&ph,&pr,&cnt)==1);
    CHECK_NEAR(ph,11.5);
    /* outermost upper bound inclusive (el=20, az=15) */
    CHECK(obscorr_lookup(oc.tab,satid2no("G01"),0,20.0,15.0,&ph,&pr,&cnt)==1);
    CHECK_NEAR(ph,12.5);
    /* out of range */
    CHECK(obscorr_lookup(oc.tab,satid2no("G01"),0, 9.99,2.5,&ph,&pr,&cnt)==0);
    CHECK(obscorr_lookup(oc.tab,satid2no("G01"),0,20.01,2.5,&ph,&pr,&cnt)==0);
    CHECK(obscorr_lookup(oc.tab,satid2no("G01"),0,12.5,15.01,&ph,&pr,&cnt)==0);

    /* azimuth normalization: -355 -> 5, 367.5 -> 7.5, -5 -> 355 (outside) */
    CHECK(obscorr_lookup(oc.tab,satid2no("G01"),0,12.5,-355.0,&ph,&pr,&cnt)==1);
    CHECK_NEAR(ph,1.5);
    CHECK(obscorr_lookup(oc.tab,satid2no("G01"),0,12.5,367.5,&ph,&pr,&cnt)==1);
    CHECK_NEAR(ph,1.5);
    CHECK(obscorr_lookup(oc.tab,satid2no("G01"),0,12.5,-5.0,&ph,&pr,&cnt)==0);

    /* band 2 and other systems */
    CHECK(obscorr_lookup(oc.tab,satid2no("G01"),1,15.0,7.0,&ph,&pr,&cnt)==1);
    CHECK_NEAR(ph,7.5);
    CHECK(obscorr_lookup(oc.tab,satid2no("R07"),1,15.0,7.0,&ph,&pr,&cnt)==1);
    CHECK_NEAR(ph,8.5);
    CHECK(obscorr_lookup(oc.tab,satid2no("J01"),0,15.0,7.0,&ph,&pr,&cnt)==1);
    CHECK_NEAR(ph,1.5);
    CHECK(obscorr_lookup(oc.tab,satid2no("C05"),0,15.0,7.0,&ph,&pr,&cnt)==1);
    CHECK_NEAR(ph,2.5);
    CHECK(obscorr_lookup(oc.tab,satid2no("E11"),0,15.0,7.0,&ph,&pr,&cnt)==1);
    CHECK_NEAR(ph,3.5);

    /* no grid for that satellite/band */
    CHECK(obscorr_lookup(oc.tab,satid2no("R07"),0,15.0,7.0,&ph,&pr,&cnt)==0);
    CHECK(obscorr_lookup(oc.tab,satid2no("G02"),0,15.0,7.0,&ph,&pr,&cnt)==0);
    /* invalid arguments */
    CHECK(obscorr_lookup(oc.tab,0,0,15.0,7.0,&ph,&pr,&cnt)==0);
    CHECK(obscorr_lookup(oc.tab,satid2no("G01"),2,15.0,7.0,&ph,&pr,&cnt)==0);
    CHECK(obscorr_lookup(NULL,satid2no("G01"),0,15.0,7.0,&ph,&pr,&cnt)==0);

    obscorr_free(&oc);
    CHECK(oc.tab==NULL);
    remove(CSVFILE);
}
/* test duplicate rows (last wins) and missing cells ---------------------------*/
static void test_duplicate_missing(void)
{
    obscorr_t oc;
    double ph,pr;
    int cnt;

    printf("test_duplicate_missing\n");

    /* 2x3 bounds derived from partial cells; cell (0,1) is missing.
       cell (0,0) duplicated: the later row must win */
    write_file(CSVFILE,
        "G01,L1,10,15,0,5,1.0,10.0,1\n"
        "G01,L1,10,15,10,15,3.0,30.0,3\n"
        "G01,L1,15,20,5,10,5.0,50.0,5\n"
        "G01,L1,10,15,0,5,9.0,90.0,9\n");   /* duplicate of the first row */

    CHECK(obscorr_init(&oc,CSVFILE,"",0)==1);

    CHECK(obscorr_lookup(oc.tab,satid2no("G01"),0,12.0,2.0,&ph,&pr,&cnt)==1);
    CHECK_NEAR(ph,9.0);                     /* last record wins */
    CHECK_NEAR(pr,90.0);
    CHECK(cnt==9);
    CHECK(obscorr_lookup(oc.tab,satid2no("G01"),0,12.0,12.0,&ph,&pr,&cnt)==1);
    CHECK_NEAR(ph,3.0);
    /* missing cell (el 10-15, az 5-10) */
    CHECK(obscorr_lookup(oc.tab,satid2no("G01"),0,12.0,7.0,&ph,&pr,&cnt)==0);

    obscorr_free(&oc);
    remove(CSVFILE);
}
/* test non-uniform grid and gap between cells ---------------------------------*/
static void test_nonuniform(void)
{
    obscorr_t oc;
    double ph,pr;
    int cnt;

    printf("test_nonuniform\n");

    /* el cells [0,10),[10,15),[15,90]; az cells [0,10),[20,30) (gap 10-20) */
    write_file(CSVFILE,
        "G01,L1,0,10,0,10,1.0,1.0,1\n"
        "G01,L1,0,10,20,30,2.0,2.0,2\n"
        "G01,L1,10,15,0,10,3.0,3.0,3\n"
        "G01,L1,10,15,20,30,4.0,4.0,4\n"
        "G01,L1,15,90,0,10,5.0,5.0,5\n"
        "G01,L1,15,90,20,30,6.0,6.0,6\n");

    CHECK(obscorr_init(&oc,CSVFILE,"",0)==1);
    CHECK(oc.tab->grid[satid2no("G01")-1][0]->uniform==0);
    CHECK(oc.tab->grid[satid2no("G01")-1][0]->nel==3);
    CHECK(oc.tab->grid[satid2no("G01")-1][0]->naz==2);

    CHECK(obscorr_lookup(oc.tab,satid2no("G01"),0, 5.0, 5.0,&ph,&pr,&cnt)==1);
    CHECK_NEAR(ph,1.0);
    CHECK(obscorr_lookup(oc.tab,satid2no("G01"),0,50.0,25.0,&ph,&pr,&cnt)==1);
    CHECK_NEAR(ph,6.0);
    /* boundaries of non-uniform search */
    CHECK(obscorr_lookup(oc.tab,satid2no("G01"),0,10.0, 5.0,&ph,&pr,&cnt)==1);
    CHECK_NEAR(ph,3.0);
    CHECK(obscorr_lookup(oc.tab,satid2no("G01"),0,90.0,25.0,&ph,&pr,&cnt)==1);
    CHECK_NEAR(ph,6.0);                     /* last-cell upper bound incl. */
    CHECK(obscorr_lookup(oc.tab,satid2no("G01"),0,90.01,25.0,&ph,&pr,&cnt)==0);
    /* azimuth gap [10,20) */
    CHECK(obscorr_lookup(oc.tab,satid2no("G01"),0,50.0,15.0,&ph,&pr,&cnt)==0);
    /* azimuth upper bound of last az cell inclusive */
    CHECK(obscorr_lookup(oc.tab,satid2no("G01"),0,50.0,30.0,&ph,&pr,&cnt)==1);
    CHECK_NEAR(ph,6.0);
    CHECK(obscorr_lookup(oc.tab,satid2no("G01"),0,50.0,30.01,&ph,&pr,&cnt)==0);

    obscorr_free(&oc);
    remove(CSVFILE);
}
/* test invalid records are skipped while valid ones load ----------------------*/
static void test_invalid_records(void)
{
    obscorr_t oc;
    double ph,pr;
    int cnt;

    printf("test_invalid_records\n");

    write_file(CSVFILE,
        "prn,band,el_lo,el_hi,az_lo,az_hi,ph,pr,count\n"
        "X99,L1,10,20,0,10,1.0,1.0,1\n"     /* invalid satellite */
        "G01,L5,10,20,0,10,1.0,1.0,1\n"     /* invalid band */
        "G01,L1,10,20\n"                    /* truncated */
        "not a record at all\n"
        "G02,L1,10,20,0,10,2.0,2.0,2\n");   /* the only valid record */

    CHECK(obscorr_init(&oc,CSVFILE,"",0)==1);
    CHECK(obscorr_lookup(oc.tab,satid2no("G02"),0,15.0,5.0,&ph,&pr,&cnt)==1);
    CHECK_NEAR(ph,2.0);
    CHECK(obscorr_lookup(oc.tab,satid2no("G01"),0,15.0,5.0,&ph,&pr,&cnt)==0);
    obscorr_free(&oc);

    /* file with no valid record at all */
    write_file(CSVFILE,"# empty\nX99,L1,10,20,0,10,1.0,1.0,1\n");
    CHECK(obscorr_init(&oc,CSVFILE,"",0)==0);
    CHECK(oc.stat==-1);
    CHECK(oc.tab==NULL);
    CHECK(obscorr_lookup(oc.tab,satid2no("G01"),0,15.0,5.0,&ph,&pr,&cnt)==0);
    obscorr_free(&oc);

    /* missing file */
    CHECK(obscorr_init(&oc,"test_obscorr_no_such_file.csv","",0)==0);
    CHECK(oc.stat==-1);
    obscorr_free(&oc);

    remove(CSVFILE);
}
/* test binary format (little-endian host assumed) -----------------------------*/
static void test_binary(void)
{
    obscorr_t oc;
    FILE *fp;
    uint8_t head[28]={0},rec[32]={0};
    uint16_t version=1,endian=0x1234;
    uint32_t nrec=1;
    float v[6]={30.0f,35.0f,100.0f,105.0f,0.25f,1.5f};
    int32_t cnt_in=77;
    double ph,pr;
    int cnt;

    printf("test_binary\n");

    memcpy(head,"OCB1",4);
    memcpy(head+4,&version,2);
    memcpy(head+6,&endian,2);
    memcpy(head+8,&nrec,4);
    rec[0]='G'; rec[1]=5; rec[2]=1;         /* G05 band L1 */
    memcpy(rec+4,v,24);
    memcpy(rec+28,&cnt_in,4);

    fp=fopen(BINFILE,"wb");
    fwrite(head,1,28,fp);
    fwrite(rec,1,32,fp);
    fclose(fp);

    CHECK(obscorr_init(&oc,BINFILE,"",0)==1);
    CHECK(obscorr_lookup(oc.tab,satid2no("G05"),0,32.0,102.0,&ph,&pr,&cnt)==1);
    CHECK_NEAR(ph,0.25);
    CHECK_NEAR(pr,1.5);
    CHECK(cnt==77);
    CHECK(obscorr_lookup(oc.tab,satid2no("G05"),1,32.0,102.0,&ph,&pr,&cnt)==0);
    obscorr_free(&oc);

    /* endian marker mismatch is rejected */
    head[6]=0x12; head[7]=0x34;             /* big-endian marker */
    fp=fopen(BINFILE,"wb");
    fwrite(head,1,28,fp);
    fwrite(rec,1,32,fp);
    fclose(fp);
    CHECK(obscorr_init(&oc,BINFILE,"",0)==0);
    CHECK(oc.tab==NULL);
    obscorr_free(&oc);

    remove(BINFILE);
}
/* test reload on file update and keep-old-table on broken update --------------*/
static void test_reload(void)
{
    obscorr_t oc;
    double ph,pr;
    int cnt;

    printf("test_reload (takes a few seconds)\n");

    write_file(CSVFILE,"G01,L1,10,20,0,10,1.0,10.0,1\n");
    backdate(CSVFILE,30);
    CHECK(obscorr_init(&oc,CSVFILE,"",1)==1);
    CHECK(oc.reload_int==1);

    /* not reloaded before the check interval elapses */
    CHECK(obscorr_reload(&oc)==0);

    /* update the file: reload picks up new values */
    write_file(CSVFILE,"G01,L1,10,20,0,10,2.0,20.0,2\n");
    backdate(CSVFILE,10);
    sleepms(1100);
    CHECK(obscorr_reload(&oc)==1);
    CHECK(obscorr_lookup(oc.tab,satid2no("G01"),0,15.0,5.0,&ph,&pr,&cnt)==1);
    CHECK_NEAR(ph,2.0);

    /* unchanged mtime: no reload */
    sleepms(1100);
    CHECK(obscorr_reload(&oc)==0);

    /* broken update: previous table is kept */
    write_file(CSVFILE,"X99,broken\n");
    backdate(CSVFILE,5);
    sleepms(1100);
    CHECK(obscorr_reload(&oc)==0);
    CHECK(oc.tab!=NULL);
    CHECK(obscorr_lookup(oc.tab,satid2no("G01"),0,15.0,5.0,&ph,&pr,&cnt)==1);
    CHECK_NEAR(ph,2.0);

    /* fixed again: reload recovers */
    write_file(CSVFILE,"G01,L1,10,20,0,10,3.0,30.0,3\n");
    backdate(CSVFILE,5);
    sleepms(1100);
    CHECK(obscorr_reload(&oc)==1);
    CHECK(obscorr_lookup(oc.tab,satid2no("G01"),0,15.0,5.0,&ph,&pr,&cnt)==1);
    CHECK_NEAR(ph,3.0);

    /* file modified "just now": reload deferred by settle time */
    write_file(CSVFILE,"G01,L1,10,20,0,10,4.0,40.0,4\n");
    sleepms(1100);
    CHECK(obscorr_reload(&oc)==0);
    CHECK(obscorr_lookup(oc.tab,satid2no("G01"),0,15.0,5.0,&ph,&pr,&cnt)==1);
    CHECK_NEAR(ph,3.0);
    /* after settle time it is loaded */
    sleepms(2200);
    CHECK(obscorr_reload(&oc)==1);
    CHECK(obscorr_lookup(oc.tab,satid2no("G01"),0,15.0,5.0,&ph,&pr,&cnt)==1);
    CHECK_NEAR(ph,4.0);

    obscorr_free(&oc);
    remove(CSVFILE);
}
/* test startup without file then recovery by reload ---------------------------*/
static void test_late_file(void)
{
    obscorr_t oc;
    double ph,pr;
    int cnt;

    printf("test_late_file\n");

    remove(CSVFILE);
    CHECK(obscorr_init(&oc,CSVFILE,"",1)==0);
    CHECK(oc.tab==NULL);

    write_file(CSVFILE,"G01,L1,10,20,0,10,1.5,15.0,1\n");
    backdate(CSVFILE,10);
    sleepms(1100);
    CHECK(obscorr_reload(&oc)==1);
    CHECK(oc.stat==1);
    CHECK(obscorr_lookup(oc.tab,satid2no("G01"),0,15.0,5.0,&ph,&pr,&cnt)==1);
    CHECK_NEAR(ph,1.5);

    obscorr_free(&oc);
    remove(CSVFILE);
}
/* main ------------------------------------------------------------------------*/
int main(void)
{
    test_csv_uniform();
    test_duplicate_missing();
    test_nonuniform();
    test_invalid_records();
    test_binary();
    test_reload();
    test_late_file();

    printf("---\n%s: pass=%d fail=%d\n",n_fail?"NG":"OK",n_pass,n_fail);
    return n_fail?1:0;
}
