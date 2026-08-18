/*------------------------------------------------------------------------------
* snrmp.c : output SNR and multipath linear combination from RINEX files
*
* usage: snrmp [-el elmask] [-ts y/m/d h:m:s] [-te y/m/d h:m:s]
*              obsfile navfile outfile
*
*   reads a RINEX observation file and a RINEX navigation file and writes
*   a text file with, for each observation record:
*     GPS time (YYYY-mm-DD HH:MM:SS), satellite id, azimuth (deg), elevation (deg),
*     L1 SNR (dBHz), MP1 (m), L2 SNR (dBHz), MP2 (m)
*
*   ported from rtkplot_qt Plot::updateObservation()/updateMp()/saveSnrMp()
*   (app/qtapp/rtkplot_qt/plotdata.cpp)
*-----------------------------------------------------------------------------*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "rtklib.h"

#define SQR(x)  ((x)*(x))
#define TTOL    DTTOL   /* time-difference tolerance (s) */

/* compute azimuth/elevation for all observations ----------------------------
 * (port of Plot::updateObservation, receiver position = RINEX header position
 *  with fallback to single point positioning) */
static void update_azel(const obs_t *obs, const nav_t *nav, const sta_t *sta,
                        double *az, double *el)
{
    prcopt_t opt = prcopt_default;
    sol_t sol = {{0}};
    double rr[3], pos[3], azel[2*MAXOBS], e[3], rs[6], dts[2], var;
    char msg[128] = "";
    int i, j, k, sat, svh, posvalid = 0;

    if (norm(sta->pos, 3) > 0.0) { /* RINEX header position */
        matcpy(rr, sta->pos, 3, 1);
        ecef2pos(rr, pos);
        posvalid = 1;
    }
    for (i = 0; i < obs->n; i = j) {
        gtime_t time = obs->data[i].time;

        for (j = i; j < obs->n; j++)
            if (timediff(obs->data[j].time, time) > TTOL) break;

        if (!posvalid) { /* fallback: single point position */
            opt.err[0] = 900.0;
            if (pntpos(obs->data + i, j - i, nav, &opt, &sol, azel, NULL, msg) &&
                norm(sol.rr, 3) > 0.0) {
                matcpy(rr, sol.rr, 3, 1);
                ecef2pos(rr, pos);
                posvalid = 1;
                fprintf(stderr, "receiver position by single point positioning: "
                        "%.4f %.4f %.4f\n", pos[0]*R2D, pos[1]*R2D, pos[2]);
            } else {
                continue;
            }
        }
        for (k = 0; k < j - i; k++) {
            sat = obs->data[i+k].sat;

            if (!satpos(time, time, sat, EPHOPT_BRDC, nav, rs, dts, &var, &svh))
                continue;
            if (geodist(rs, rr, e) > 0.0) {
                satazel(pos, e, azel);
                if (azel[0] < 0.0) azel[0] += 2.0*PI;
            } else {
                azel[0] = azel[1] = 0.0;
            }
            az[i+k] = azel[0];
            el[i+k] = azel[1];
        }
    }
}
/* compute multipath linear combination ---------------------------------------
 * (port of Plot::updateMp)
 * mp[freq][obs-index], 0.0 means no value */
static void update_mp(const obs_t *obs, const nav_t *nav,
                      double *mp[NFREQ+NEXOBS])
{
    int i, j, k, m, n, sat;

    for (i = 0; i < obs->n; i++) {
        obsd_t *data = obs->data + i;
        /* choose two frequencies to calculate reference I */
        double freq1 = 0.0, freq2 = 0.0, I = 0.0;
        for (j = 0; j < NFREQ + NEXOBS; j++) {
            freq1 = sat2freq(data->sat, data->code[j], nav);
            if (data->L[j] == 0.0 || freq1 == 0.0) continue;
            for (k = j + 1; k < NFREQ + NEXOBS; k++) {
                freq2 = sat2freq(data->sat, data->code[k], nav);
                if (data->L[k] == 0.0 || freq2 == 0.0 || freq1 == freq2) continue;
                I = -CLIGHT * (data->L[j] / freq1 - data->L[k] / freq2) /
                    (1.0 - SQR(freq1 / freq2));
                break;
            }
            break;
        }
        if (freq1 == 0.0 || freq2 == 0.0) continue;

        for (j = 0; j < NFREQ + NEXOBS; j++) {
            double freq = sat2freq(data->sat, data->code[j], nav);
            if (data->P[j] == 0.0 || data->L[j] == 0.0 || freq == 0.0) continue;
            mp[j][i] = data->P[j] - CLIGHT * data->L[j] / freq -
                       2.0 * SQR(freq1 / freq) * I;
        }
    }
    /* remove per-arc bias (cycle-slip or 5 m jump starts a new arc) */
    for (sat = 1; sat <= MAXSAT; sat++) {
        for (j = 0; j < NFREQ + NEXOBS; j++) {
            double B = 0.0;
            m = 0;
            for (i = 0, n = 0; i < obs->n; i++) {
                obsd_t *data = obs->data + i;
                if (data->sat != sat) continue;
                if ((data->LLI[j] & 1) || (data->LLI[0] & 1) || (data->LLI[1] & 1) ||
                    fabs(mp[j][i] - B) > 5.0) {
                    for (k = m; k < i; k++) {
                        if (obs->data[k].sat == sat && mp[j][k] != 0.0) mp[j][k] -= B;
                    }
                    n = 0; m = i; B = 0.0;
                }
                if (mp[j][i] != 0.0) B += (mp[j][i] - B) / ++n;
            }
            for (k = m; k < obs->n; k++) {
                if (obs->data[k].sat == sat && mp[j][k] != 0.0) mp[j][k] -= B;
            }
        }
    }
}
/* write output file (format based on Plot::saveSnrMp) -----------------------*/
static int write_snrmp(const char *file, const obs_t *obs,
                       double *mp[NFREQ+NEXOBS], const double *az,
                       const double *el, double elmask)
{
    FILE *fp;
    char id[8];
    double ep[6];
    int i;

    if (!(fp = fopen(file, "w"))) {
        fprintf(stderr, "output file open error: %s\n", file);
        return 0;
    }
    fprintf(fp, "%% %-19s %6s %8s %8s %9s %10s %9s %10s\n",
            "TIME(GPST)", "SAT", "AZ(deg)", "EL(deg)",
            "SNR1(dBHz)", "MP1(m)", "SNR2(dBHz)", "MP2(m)");

    for (i = 0; i < obs->n; i++) {
        const obsd_t *data = obs->data + i;

        if (data->SNR[0] == 0.0 && data->SNR[1] == 0.0 &&
            mp[0][i] == 0.0 && mp[1][i] == 0.0) continue;
        if (elmask > 0.0 && el[i] * R2D < elmask) continue;

        time2epoch(data->time, ep);
        satno2id(data->sat, id);
        fprintf(fp, "%04.0f-%02.0f-%02.0f %02.0f:%02.0f:%02.0f %6s %8.1f %8.1f "
                "%9.2f %10.4f %9.2f %10.4f\n",
                ep[0], ep[1], ep[2], ep[3], ep[4], floor(ep[5] + 0.5),
                id, az[i] * R2D, el[i] * R2D,
                data->SNR[0], mp[0][i], data->SNR[1], mp[1][i]);
    }
    fclose(fp);
    return 1;
}
/* print usage ---------------------------------------------------------------*/
static void print_usage(void)
{
    fprintf(stderr,
        "usage: snrmp [-el elmask] [-ts y/m/d h:m:s] [-te y/m/d h:m:s] "
        "obsfile navfile outfile\n"
        "  obsfile : RINEX observation file\n"
        "  navfile : RINEX navigation file\n"
        "  outfile : output text file\n"
        "  -el     : elevation mask (deg) (default: none)\n"
        "  -ts/-te : time start/end (GPST)\n");
}
/* main ----------------------------------------------------------------------*/
int main(int argc, char **argv)
{
    gtime_t ts = {0}, te = {0};
    obs_t obs = {0};
    nav_t nav = {0};
    sta_t sta = {{0}};
    double *mp[NFREQ+NEXOBS] = {NULL}, *az = NULL, *el = NULL, ep[6];
    double elmask = 0.0;
    const char *infile[2] = {NULL, NULL}, *outfile = NULL;
    int i, j, n = 0, stat = 0;

    for (i = 1; i < argc; i++) {
        if (!strcmp(argv[i], "-el") && i + 1 < argc) {
            elmask = atof(argv[++i]);
        } else if (!strcmp(argv[i], "-ts") && i + 2 < argc) {
            sscanf(argv[++i], "%lf/%lf/%lf", ep, ep+1, ep+2);
            sscanf(argv[++i], "%lf:%lf:%lf", ep+3, ep+4, ep+5);
            ts = epoch2time(ep);
        } else if (!strcmp(argv[i], "-te") && i + 2 < argc) {
            sscanf(argv[++i], "%lf/%lf/%lf", ep, ep+1, ep+2);
            sscanf(argv[++i], "%lf:%lf:%lf", ep+3, ep+4, ep+5);
            te = epoch2time(ep);
        } else if (argv[i][0] == '-' && argv[i][1]) {
            print_usage();
            return 1;
        } else if (n < 2) {
            infile[n++] = argv[i];
        } else if (!outfile) {
            outfile = argv[i];
        } else {
            print_usage();
            return 1;
        }
    }
    if (n < 2 || !outfile) {
        print_usage();
        return 1;
    }
    /* read RINEX observation and navigation files */
    if (readrnxt(infile[0], 1, ts, te, 0.0, "", &obs, &nav, &sta) < 0) {
        fprintf(stderr, "observation file read error: %s\n", infile[0]);
        return 1;
    }
    if (readrnxt(infile[1], 1, ts, te, 0.0, "", NULL, &nav, NULL) < 0) {
        fprintf(stderr, "navigation file read error: %s\n", infile[1]);
        return 1;
    }
    if (obs.n <= 0) {
        fprintf(stderr, "no observation data: %s\n", infile[0]);
        return 1;
    }
    if (nav.n <= 0 && nav.ng <= 0 && nav.ns <= 0) {
        fprintf(stderr, "no navigation data: %s\n", infile[1]);
        return 1;
    }
    uniqnav(&nav);
    sortobs(&obs);

    /* azimuth/elevation */
    az = (double *)calloc(obs.n, sizeof(double));
    el = (double *)calloc(obs.n, sizeof(double));
    for (j = 0; j < NFREQ + NEXOBS; j++)
        mp[j] = (double *)calloc(obs.n, sizeof(double));
    if (!az || !el || !mp[NFREQ+NEXOBS-1]) {
        fprintf(stderr, "memory allocation error\n");
        return 1;
    }
    update_azel(&obs, &nav, &sta, az, el);

    /* multipath linear combination */
    update_mp(&obs, &nav, mp);

    /* write output */
    stat = write_snrmp(outfile, &obs, mp, az, el, elmask);

    free(az); free(el);
    for (j = 0; j < NFREQ + NEXOBS; j++) free(mp[j]);
    freeobs(&obs);
    freenav(&nav, 0xFF);

    return stat ? 0 : 1;
}
