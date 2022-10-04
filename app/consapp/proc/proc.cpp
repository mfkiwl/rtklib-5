#include <iostream>
#include <cmath>
#include <cstdio>
#include <vector>
#include <string>
#include <algorithm>
#include <ctime>
#include <functional> 

#include "rtklib.h"

#include "obs.h"
#include "engine.h"

#ifdef _PLOT_
#include "..\matplotlib-cpp\matplotlibcpp.h"

namespace plt = matplotlibcpp;

#endif

#ifndef PI
#define	PI 3.14159265358979
#endif

#ifndef D2R
#define D2R (PI/180.0)
#endif

#ifndef R2D
#define R2D (180.0/PI)
#endif

#ifndef TWOPI
#define	TWOPI (PI+PI)
#endif

#ifndef ae_WGS84
#define ae_WGS84 6378137.0
#endif

#ifndef finv_WGS84
#define finv_WGS84 298.257223563
#endif

#ifndef grav_WGS84
#define	grav_WGS84 9.7803267714e0
#endif

#define MAXFIELD 100

static int parse_fields(char* const buffer, char** val)
{
	char* p, *q;
	int n = 0;

	/* parse fields */
	for (p = buffer; *p && n < MAXFIELD; p = q + 1) {
		if (p == NULL) break;
		if ((q = strchr(p, ',')) || (q = strchr(p, '*')) || (q = strchr(p, '\n')) || (q = strchr(p, '\r'))) {
			val[n++] = p; *q = '\0';
		}
		else break;
	}
	return n;
}

static double lat2local(double lat, double* lat2north)
{
	double f_WGS84 = (1.0 / finv_WGS84);
	double e2WGS84 = (2.0 * f_WGS84 - f_WGS84 * f_WGS84);
	double slat = sin(lat);
	double clat = cos(lat);
	double one_e2_slat2 = 1.0 - e2WGS84 * slat * slat;
	double Rn = ae_WGS84 / sqrt(one_e2_slat2);
	double Rm = Rn * (1.0 - e2WGS84) / (one_e2_slat2);
	*lat2north = Rm;
	return Rn * clat;
}

static void set_output_file_name(const char* fname, const char* key, char *outfilename)
{
	char filename[255] = { 0 };
	strcpy(filename, fname);
	char* temp = strrchr(filename, '.');
	if (temp) temp[0] = '\0';
	sprintf(outfilename, "%s-%s", filename, key);
}

static FILE* set_output_file(const char* fname, const char* key)
{
	char filename[255] = { 0 }, outfilename[255] = { 0 };
	strcpy(filename, fname);
	char* temp = strrchr(filename, '.');
	if (temp) temp[0] = '\0';
	sprintf(outfilename, "%s-%s", filename, key);
	return fopen(outfilename, "w");
}

static void deg2dms(double deg, double* dms)
{
	double sign = deg < 0.0 ? (-1.0) : (1.0), a = fabs(deg);
	dms[0] = floor(a); a = (a - dms[0]) * 60.0;
	dms[1] = floor(a); a = (a - dms[1]) * 60.0;
	dms[2] = a; dms[0] *= sign;
}

extern int outnmea_gga(unsigned char* buff, float time, int type, double* blh, int ns, float dop, float age)
{
	double h, ep[6], dms1[3], dms2[3];
	char* p = (char*)buff, * q, sum;

	if (type != 1 && type != 4 && type != 5) {
		p += sprintf(p, "$GPGGA,,,,,,,,,,,,,,");
		for (q = (char*)buff + 1, sum = 0; *q; q++) sum ^= *q;
		p += sprintf(p, "*%02X%c%c", sum, 0x0D, 0x0A);
		return (int)(p - (char*)buff);
	}
	time -= 18.0;
	ep[2] = floor(time / (24 * 3600));
	time -= (float)(ep[2] * 24 * 3600.0);
	ep[3] = floor(time / 3600);
	time -= (float)(ep[3] * 3600);
	ep[4] = floor(time / 60);
	time -= (float)(ep[4] * 60);
	ep[5] = time;
	h = 0.0;
	deg2dms(fabs(blh[0]) * 180 / PI, dms1);
	deg2dms(fabs(blh[1]) * 180 / PI, dms2);
	p += sprintf(p, "$GPGGA,%02.0f%02.0f%06.3f,%02.0f%010.7f,%s,%03.0f%010.7f,%s,%d,%02d,%.1f,%.3f,M,%.3f,M,%.1f,",
		ep[3], ep[4], ep[5], dms1[0], dms1[1] + dms1[2] / 60.0, blh[0] >= 0 ? "N" : "S",
		dms2[0], dms2[1] + dms2[2] / 60.0, blh[1] >= 0 ? "E" : "W", type,
		ns, dop, blh[2] - h, h, age);
	for (q = (char*)buff + 1, sum = 0; *q; q++) sum ^= *q; /* check-sum */
	p += sprintf(p, "*%02X%c%c", sum, 0x0D, 0x0A);
	return (int)(p - (char*)buff);
}

/* read the data until return 1, the ephemeris will also save to nav_t */
static int read_rtcm_data(rtcm_t* rtcm, nav_t* nav, FILE *fRTCM)
{
	int data = 0;
	int ret = 0;
	int sys = 0;
	int prn = 0;
	while (fRTCM && !feof(fRTCM) && (data = fgetc(fRTCM)) != EOF)
	{
		ret = input_rtcm3(rtcm, (uint8_t)data);
		if (ret == 2)
		{
			sys = satsys(rtcm->ephsat, &prn);
			if (sys == SYS_GLO)
			{
				nav->geph[prn - 1] = rtcm->nav.geph[prn - 1];
			}
			else 
			{
				if (rtcm->ephset)
					prn = prn;
				nav->eph[rtcm->ephsat - 1 + MAXSAT * rtcm->ephset] = rtcm->nav.eph[rtcm->ephsat - 1 + MAXSAT * rtcm->ephset];
			}
		}
		else if (ret == 10)
		{
			/* SSR */
		}
		else if (ret == 5)
		{
			/* station coordinate */
		}
		if (ret == 1) break;
	}
	return ret;
}

// trim from start
static inline std::string& ltrim(std::string& s) {
	s.erase(s.begin(), std::find_if(s.begin(), s.end(),
		std::not1(std::ptr_fun<int, int>(std::isspace))));
	return s;
}

// trim from end
static inline std::string& rtrim(std::string& s) {
	s.erase(std::find_if(s.rbegin(), s.rend(),
		std::not1(std::ptr_fun<int, int>(std::isspace))).base(), s.end());
	return s;
}

// trim from both ends
static inline std::string& trim(std::string& s) {
	return ltrim(rtrim(s));
}

/* satellite system G,R,E,C,J,I,S */
static char sat2char(int sat, int* prn)
{
	int sys = satsys(sat, prn);
	if (sys == SYS_QZS) *prn -= MINPRNQZS;
	if (sys == SYS_GPS) return 'G';
	else if (sys == SYS_GLO) return 'R';
	else if (sys == SYS_GAL) return 'E';
	else if (sys == SYS_CMP) return 'C';
	else if (sys == SYS_SBS) return 'S';
	else if (sys == SYS_QZS) return 'J';
	else if (sys == SYS_IRN) return 'I';
	else return ' ';
}

static int convert_epoch_data(obs_t* obs, double* pos, nav_t* nav, epoch_t* epoch)
{
	double* rs = 0, * dts = 0, * var = 0;
	int i, j, f, svh[MAXOBS * 2] = { 0 }, n = obs->n;
	gtime_t time = { 0 };
	double ws = 0.0;
	int wk = 0;
	const obsd_t* obsd = obs->data + 0;
	int nf = 0;
	int nfloc[NFREQ + NEXOBS] = { 0 };
	double freq[NFREQ + NEXOBS] = { 0 };
	int prn = 0;
	char sys = 0;

	memset(epoch, 0, sizeof(epoch_t));

	if (n < 1) return 0;

	time = obs->data[0].time;
	rs = mat(6, n); dts = mat(2, n); var = mat(1, n);

	/* satellite positions/clocks */
	satposs(time, obs->data, n, nav, EPHOPT_BRDC, rs, dts, var, svh);

	/* output data */
	ws = time2gpst(time, &wk);

	epoch->time = ws;
	epoch->pos[0] = pos[0];
	epoch->pos[1] = pos[1];
	epoch->pos[2] = pos[2];

	for (i = 0, obsd = obs->data + i; i < n; ++i, ++obsd)
	{
		ws = time2gpst(obsd->time, &wk);
		sys = sat2char(obsd->sat, &prn);
		nf = 0;
		memset(freq, 0, sizeof(freq));
		memset(nfloc, 0, sizeof(nfloc));
		for (f = 0; f < (NFREQ + NEXOBS); ++f)
		{
			if (obsd->code[f] > 0)
			{
				freq[f] = sat2freq(obsd->sat, obsd->code[f], nav);
				if (freq[f] > 0.01)
				{
					nfloc[nf++] = f;
				}
			}
		}
		if (nf == 0) continue;
		sat_obs_t* sat_obs = epoch->obs + epoch->n;
		sat_obs->sys = sys;
		sat_obs->prn = prn;
		sat_obs->sat = obsd->sat;
		for (j = 0; j < nf; ++j)
		{
			f = nfloc[j];
			sat_obs->code[j] = obsd->code[f];
			sat_obs->P[j] = obsd->P[f];
			sat_obs->L[j] = obsd->L[f];
			sat_obs->D[j] = obsd->D[f];
			sat_obs->S[j] = obsd->SNR[f] * SNR_UNIT;
			sat_obs->wave[j] = freq[f] > 0.1 ? CLIGHT / freq[f] : 0;
		}
		sat_vec_t* sat_vec = epoch->vec + epoch->n;
		if (fabs(rs[i * 6 + 0]) < 0.1 || fabs(rs[i * 6 + 1]) < 0.1 || fabs(rs[i * 6 + 2]) < 0.1)
		{

		}
		else
		{
			sat_vec->rs[0] = rs[i * 6 + 0];
			sat_vec->rs[1] = rs[i * 6 + 1];
			sat_vec->rs[2] = rs[i * 6 + 2];
			sat_vec->rs[3] = rs[i * 6 + 3];
			sat_vec->rs[4] = rs[i * 6 + 4];
			sat_vec->rs[5] = rs[i * 6 + 5];
			sat_vec->dts[0] = dts[i * 2 + 0] * CLIGHT;
			sat_vec->dts[1] = dts[i * 2 + 1] * CLIGHT;
		}
		epoch->n++;
	}
	free(rs); free(dts); free(var);
	return epoch->n;
}

static int print_obs_data(int staid, epoch_t* epoch, FILE* fOUT)
{
	if (fOUT && epoch->n > 0)
	{
		char buff[1024] = { 0 };
		char* p = (char*)buff, * q, sum;
		/* output position $GNSSPOS, staid, nsat, wk, ws, X,Y,Z */
		if (fabs(epoch->pos[0]) < 0.1 || fabs(epoch->pos[1]) < 0.1 || fabs(epoch->pos[2]) < 0.1)
		{

		}
		else
		{
			p += sprintf(p, "$GNSSPOS,%04i,%03i,%4i,%10.4f,%14.4f,%14.4f,%14.4f", staid, epoch->n, epoch->wk, epoch->time, epoch->pos[0], epoch->pos[1], epoch->pos[2]);
			for (q = (char*)buff + 1, sum = 0; *q; q++) sum ^= *q; /* check-sum */
			p += sprintf(p, "*%02X%c%c", sum, 0x0D, 0x0A);
			fprintf(fOUT, "%s", buff);
		}
		sat_obs_t* sat_obs = epoch->obs + 0;
		sat_vec_t* sat_vec = epoch->vec + 0;
		int i = 0, j = 0, nf = 0, loc[MAX_FRQ] = { 0 }, f = 0;
		for (; i < epoch->n; ++i, ++sat_obs, ++sat_vec)
		{
			nf = 0;
			for (j = 0; j < MAX_FRQ; ++j)
			{
				if (sat_obs->code[j])
					loc[nf++] = j;
			}
			if (nf == 0) continue;
			memset(buff, 0, sizeof(buff));
			p = (char*)buff;
			p += sprintf(p, "$GNSSOBS,%04i,%c%02i,%03i,%10.4f", staid, sat_obs->sys, sat_obs->prn, sat_obs->sat, epoch->time);

			if (fabs(sat_vec->rs[0]) < 0.1 || fabs(sat_vec->rs[1]) < 0.1 || fabs(sat_vec->rs[2]) < 0.1)
			{
				p += sprintf(p, ",0");
			}
			else
			{
				p += sprintf(p, ",1,%14.4f,%14.4f,%14.4f,%10.4f,%10.4f,%10.4f,%14.4f,%10.4f", sat_vec->rs[0], sat_vec->rs[1], sat_vec->rs[2], sat_vec->rs[3], sat_vec->rs[4], sat_vec->rs[5], sat_vec->dts[0], sat_vec->dts[1]);
			}

			p += sprintf(p, ",%i", nf);
			for (j = 0; j < nf; ++j)
			{
				f = loc[j];
				p += sprintf(p, ",%3i,%14.4f,%14.4f,%10.4f,%3i,%.0f", sat_obs->code[f], sat_obs->P[f], sat_obs->L[f], sat_obs->D[f], sat_obs->S[f], sat_obs->wave[f] > 0.001 ? CLIGHT / sat_obs->wave[f] : 0);
			}
			for (q = (char*)buff + 1, sum = 0; *q; q++) sum ^= *q; /* check-sum */
			p += sprintf(p, "*%02X%c%c", sum, 0x0D, 0x0A);
			fprintf(fOUT, "%s", buff);
			fflush(fOUT);
		}
	}
	return 1;
}

static int read_log_file(const char* fname, std::vector<std::string>& logfnames, int *year, int *mon, int *day)
{
	FILE* fINI = fopen(fname, "r");
	char buffer[512] = { 0 };
	char rtcmfname[512] = { 0 };
	char* val[MAXFIELD];
	int type = 0;
	int num = 0;
	char* temp = 0;
	std::time_t t = std::time(0);   // get time now
	std::tm* now = std::localtime(&t);
	*year = (now->tm_year + 1900);
	*mon = (now->tm_mon + 1);
	*day = now->tm_mday;
	while (fINI && !feof(fINI) && fgets(buffer, sizeof(buffer), fINI) != NULL)
	{
		temp = strchr(buffer, '=');
		if (temp) temp[0] = ',';
		num = parse_fields(buffer, val);
		if (num < 2) continue;
		if (strstr(val[0], "date"))
		{
			if (num > 3)
			{
				*year = atoi(val[1]);
				*mon = atoi(val[2]);
				*day = atoi(val[3]);
			}
		}
		else if (strstr(val[0], "rtcm"))
		{
			std::string rtcmfname(val[1]);
			logfnames.push_back(trim(rtcmfname));
		}
	}
	if (fINI) fclose(fINI);
	return (int)logfnames.size();
}

static int free_nav(nav_t* nav)
{
	/* free memory for observation and ephemeris buffer */
	free(nav->eph); nav->eph = NULL; nav->n = 0;
	free(nav->geph); nav->geph = NULL; nav->ng = 0;
	return 1;
}

static int init_nav(nav_t* nav)
{
	eph_t  eph0 = { 0,-1,-1 };
	geph_t geph0 = { 0,-1 };
	ssr_t ssr0 = { {{0}} };
	int i=0;

	trace(3, "init_rtcm:\n");

	for (i = 0; i < MAXSAT; i++) {
		nav->ssr[i] = ssr0;
	}
	nav->eph = NULL;
	nav->geph = NULL;

	/* reallocate memory for observation and ephemeris buffer */
	if (!(nav->eph = (eph_t*)malloc(sizeof(eph_t) * MAXSAT * 2)) ||
		!(nav->geph = (geph_t*)malloc(sizeof(geph_t) * MAXPRNGLO))) {
		free_nav(nav);
		return 0;
	}
	nav->n = MAXSAT * 2;
	nav->ng = MAXPRNGLO;
	for (i = 0; i < MAXSAT * 2; i++) nav->eph[i] = eph0;
	for (i = 0; i < MAXPRNGLO; i++) nav->geph[i] = geph0;
	return 1;
}

/* process the data from rtcm log files */
static int proc(const char* fname)
{
	std::vector<std::string> logfnames;
	/* load the ini file to get the log file lists */
	int year = 0, mon = 0, day = 0, numrcv = read_log_file(fname, logfnames, &year, &mon, &day), i = 0, j = 0, ret = 0, wk = 0;
	/* read through file and process data */
	std::vector<FILE*> rtcmfs(numrcv);
	std::vector<rtcm_t> rtcms(numrcv);
	gtime_t cur_time = { 0 }, nex_time = { 0 };
	unsigned long numofepoch = 0;
	double dt = 0, ws = 0, ep[6] = { year, mon, day, 0, 0, 0 };
	nav_t nav = { 0 }; /* global */
	epoch_t epoch = { 0 };
	/* init ephemeris data struct */
	init_nav(&nav);
	/* output raw data file for output */
	char buffer[255] = { 0 };
	sprintf(buffer, "%04i-%02i-%02i-obs.log", year, mon, day);
	FILE *fOBS = fopen(buffer, "w");
	FILE* fTMP = fopen("rtcm_log.csv", "w");
	/* open file, init and read the first epoch */
	j = 0;
	for (i = 0; i < numrcv; ++i)
	{
		rtcmfs[i] = fopen(logfnames[i].c_str(), "rb");
		init_rtcm(&rtcms[i]);
		rtcms[i].time = rtcms[i].time_s = epoch2time(ep);
		/* read start time */
		if (rtcmfs[i] && !feof(rtcmfs[i]))
		{
			ret = read_rtcm_data(&rtcms[i], &nav, rtcmfs[i]);
			if (ret == 1)
			{
				if (j == 0)
					cur_time = rtcms[i].time;
				else
				{
					dt = timediff(cur_time, rtcms[i].time);
					if (dt > 0.0)
						cur_time = rtcms[i].time;
				}
				++j;
			}
		}
		/**/
	}
	/* data processing and read more data */
	while (1)
	{
		/* process current epoch data */
		j = 0;
		for (i = 0; i < numrcv; ++i)
		{
			ws = time2gpst(rtcms[i].time, &wk);
			if (fabs(dt = timediff(cur_time, rtcms[i].time)) < 0.001)
			{
				++j;
				/* convert data and add to engine */
				printf("%10.4f,%3i,%3i,*\n", ws, rtcms[i].obs.n, i);
				if (convert_epoch_data(&rtcms[i].obs, rtcms[i].sta.pos, &nav, &epoch))
				{
					print_obs_data(i + 1, &epoch, fOBS);
					if (fTMP) fprintf(fTMP, "%10.4f,%3i,%3i,%6i\n", ws, rtcms[i].obs.n, i + 1, numofepoch);
					add_rcv_data(i + 1, &epoch); /* use the station index as the station ID */
				}
			}
			else
			{
				printf("%10.4f,%3i,%3i,-\n", ws, rtcms[i].obs.n, i);
			}
		}
		/* data process */
		if (j == 0) break; /* no more data */
		ws = time2gpst(cur_time, &wk);
		printf("%10.4f,%6i,%i/%i\n", ws, numofepoch, j, numrcv);

		/* processing */
		engine_proc();

		++numofepoch;
		/* read more data */
		j = 0;
		for (i = 0; i < numrcv; ++i)
		{
			dt = timediff(cur_time, rtcms[i].time);
			if (dt < 0.0)
			{
				if (j == 0)
					nex_time = rtcms[i].time;
				else
				{
					if (timediff(nex_time, rtcms[i].time) > 0.0)
						nex_time = rtcms[i].time;
				}
				++j;
			}
			else if (rtcmfs[i] && !feof(rtcmfs[i]))
			{
				ret = read_rtcm_data(&rtcms[i], &nav, rtcmfs[i]);
				if (ret == 1)
				{
					if (j == 0)
						nex_time = rtcms[i].time;
					else
					{
						if (timediff(nex_time, rtcms[i].time) > 0.0)
							nex_time = rtcms[i].time;
					}
					++j;
				}
			}
		}
		if (j == 0) break; /* no more new data */
		cur_time = nex_time;
	}
	/* free rtcm data struct and close rtcm files */
	for (i = 0; i < numrcv; ++i)
	{
		free_rtcm(&rtcms[i]);
		if (rtcmfs[i]) fclose(rtcmfs[i]);
	}
	/* free ephemeris data struct */
	free_nav(&nav);
	/* close the raw data file */
	if (fOBS) fclose(fOBS); fOBS = NULL;
	if (fTMP) fclose(fTMP); fTMP = NULL;
	return 0;
}


/* process the data from log csv file */
static int proc_csv(const char* fname)
{
	FILE* fLOG = fopen(fname, "r");
	char buffer[512] = { 0 };
	char* val[MAXFIELD];
	int ret = 0, num = 0, staid =-1;
	epoch_t epoch = { 0 };
	FILE* fTMP = fopen("log.csv", "w");
	unsigned long numofepoch = 0;
	double cur_time = 0;
	while (fLOG && !feof(fLOG) && (fgets(buffer, sizeof(buffer), fLOG)) != NULL)
	{
		num = parse_fields(buffer, val);
		if (num < 2) continue;
		if (strstr(val[0], "GNSSPOS") != NULL)
		{
			/* $GNSSPOS,0556,047,   0,25217.0000, -2348583.1770,  5355378.8290,  2538133.6170*44 */
			if (num < 8) continue;
			int cur_staid = atoi(val[1]);
			int nsat = atoi(val[2]);
			int wk = atoi(val[3]);
			double ws = atof(val[4]);
			double pos[3] = { atof(val[5]), atof(val[6]), atof(val[7]) };
			if ((cur_staid != staid && staid >= 0) || fabs(ws - epoch.time) > 0.01)
			{
				/* new data */
				if (staid >= 0 && epoch.n > 0)
				{
					if (fabs(cur_time - epoch.time) > 0.01)
					{
						cur_time = epoch.time;
						++numofepoch;
					}
					printf("%10.4f,%3i,%3i,%6i\n", epoch.time, epoch.n, staid, numofepoch);
					if (fTMP) fprintf(fTMP, "%10.4f,%3i,%3i,%6i\n", epoch.time, epoch.n, staid, numofepoch);
					add_rcv_data(staid, &epoch);
				}
				memset(&epoch, 0, sizeof(epoch_t));
			}
			epoch.pos[0] = pos[0];
			epoch.pos[1] = pos[1];
			epoch.pos[2] = pos[2];
			epoch.time = ws;
			staid = cur_staid;
		}
		else if (strstr(val[0], "GNSSOBS") != NULL)
		{
			/*
$GNSSOBS,0233,E09,068,25217.0000,0,4, 12, 25998903.6908,136625249.2254,    0.0000, 45,1575420000, 29, 25998919.2548,104686996.5750,    0.0000, 47,1207140000, 26, 25998922.5963,102025482.6126,    0.0000, 46,1176450000, 33, 25998914.5195,110897205.2252,    0.0000, 46,1278750000*34
$GNSSOBS,0233,E10,069,25217.0000,1,  -835034.9258, 28663447.6591, -7316489.1719,   62.1004, -749.9842,-2949.8190,  -144144.1597,   -0.0008,4, 12, 25431865.8957,133645220.4514,    0.0000, 47,1575420000, 29, 25431883.6218,102403474.7888,    0.0000, 47,1207140000, 26, 25431887.2492, 99800004.7328,    0.0000, 47,1176450000, 33, 25431879.5477,108478260.4761,    0.0000, 47,1278750000*3A
			*/
			if (num < 5) continue;
			int cur_staid = atoi(val[1]);
			if (strlen(val[2]) < 3) continue;
			uint8_t sys = val[2][0];
			uint8_t prn = atoi(val[2] + 1);
			uint8_t sat = atoi(val[3]);
			double ws = atof(val[4]);
			if ((cur_staid != staid && staid >= 0) || fabs(ws - epoch.time) > 0.01)
			{
				/* new data */
				if (staid >= 0 && epoch.n > 0)
				{
					if (fabs(cur_time - epoch.time) > 0.01)
					{
						cur_time = epoch.time;
						++numofepoch;
					}
					printf("%10.4f,%3i,%3i,%6i\n", epoch.time, epoch.n, staid, numofepoch);
					if (fTMP) fprintf(fTMP, "%10.4f,%3i,%3i,%6i\n", epoch.time, epoch.n, staid, numofepoch);
					add_rcv_data(staid, &epoch);
				}
				memset(&epoch, 0, sizeof(epoch_t));
			}
			epoch.time = ws;
			staid = cur_staid;
			int sat_pv_flag = atoi(val[5]);
			sat_vec_t* sat_vec = epoch.vec + epoch.n;
			sat_obs_t* sat_obs = epoch.obs + epoch.n;
			memset(sat_vec, 0, sizeof(sat_vec_t));
			memset(sat_obs, 0, sizeof(sat_obs_t));
			sat_obs->sys = sys;
			sat_obs->prn = prn;
			sat_obs->sat = sat;
			int loc = 14;
			if (sat_pv_flag)
			{
				/* wth satellite pos, vel, clock */
				if (num < 14) continue;
				sat_vec->rs[0] = atof(val[6]);
				sat_vec->rs[1] = atof(val[7]);
				sat_vec->rs[2] = atof(val[8]);
				sat_vec->rs[3] = atof(val[9]);
				sat_vec->rs[4] = atof(val[10]);
				sat_vec->rs[5] = atof(val[11]);
				sat_vec->dts[0] = atof(val[12]);
				sat_vec->dts[1] = atof(val[13]);
			}
			else
			{
				if (num < 6) continue;
				loc = 6;
			}
			int nf = atof(val[loc]), f = 0;
			int tnum = (loc + 1 + nf * 6);
			if (num < tnum) continue;
			if (nf > 0)
			{
				for (f = 0; f < nf; ++f)
				{
					int cur_loc = loc + 1 + f * 6;
					sat_obs->code[f] = atoi(val[cur_loc + 0]);
					sat_obs->P[f] = atof(val[cur_loc + 1]);
					sat_obs->L[f] = atof(val[cur_loc + 2]);
					sat_obs->D[f] = atof(val[cur_loc + 3]);
					sat_obs->S[f] = atoi(val[cur_loc + 4]);
					double cur_frq = atof(val[cur_loc + 5]);
					if (cur_frq > 0.1)
						sat_obs->wave[f] = CLIGHT / cur_frq;
				}
			}
			if (epoch.n < MAX_SAT)
				++epoch.n;
			else
				epoch.n = epoch.n;
		}
	}
	/* last data */
	if (staid >= 0 && epoch.n > 0)
	{
		if (fabs(cur_time - epoch.time) > 0.01)
		{
			cur_time = epoch.time;
			++numofepoch;
		}
		printf("%10.4f,%3i,%3i,%6i\n", epoch.time, epoch.n, staid, numofepoch);
		if (fTMP) fprintf(fTMP, "%10.4f,%3i,%3i,%6i\n", epoch.time, epoch.n, staid, numofepoch);
		add_rcv_data(staid, &epoch);
	}
	if (fLOG) fclose(fLOG);
	if (fTMP) fclose(fTMP);
	return 0;
}

/* need to install matplotlib-cpp at upper directory https://github.com/yydgis/matplotlib-cpp.git */
/* need to install python and add the path => normally go to cmd, and use path command to locate the python path, for example C:\Users\xxx\AppData\Local\Programs\Python\Python310 */
/* need to install numpy and add the path for example, C:\Users\yudan\AppData\Local\Programs\Python\Python310\Lib\site-packages\numpy\core\include */
/* comment out some lines in matplotlibcpp.h about NPY_INT64, NPY_UINT64*/
/* add python310.lib in the path and dependencies, note this only works for release build, need to find the python310_d.lib for debug version */

/* input file will be rtcm file, please use the convbin to convert other format (rinex, ublox, septentrio, hemisphere, swiftnav, etc.) to rtcm */

/* example ini file
date = 2022,10,2
rtcm = D:\gnss\multi_caster_logger\2022-10-2\2022-10-2-0-0-0-23km.log
rtcm = D:\gnss\multi_caster_logger\2022-10-2\2022-10-2-0-0-0-5000.log
rtcm = D:\gnss\multi_caster_logger\2022-10-2\2022-10-2-0-0-0-brdc.log
rtcm = D:\gnss\multi_caster_logger\2022-10-2\2022-10-2-0-0-0-chqg.log
rtcm = D:\gnss\multi_caster_logger\2022-10-2\2022-10-2-0-0-0-cong.log
rtcm = D:\gnss\multi_caster_logger\2022-10-2\2022-10-2-0-0-0-dlr0.log
rtcm = D:\gnss\multi_caster_logger\2022-10-2\2022-10-2-0-0-0-fsnh.log
rtcm = D:\gnss\multi_caster_logger\2022-10-2\2022-10-2-0-0-0-gylc.log
rtcm = D:\gnss\multi_caster_logger\2022-10-2\2022-10-2-0-0-0-gzzc.log
rtcm = D:\gnss\multi_caster_logger\2022-10-2\2022-10-2-0-0-0-pync.log
rtcm = D:\gnss\multi_caster_logger\2022-10-2\2022-10-2-0-0-0-qyqc.log
*/

int main(int argc, char** argv)
{
	if (argc < 2) /* 1 */
	{
		/* */
		printf("%s inifname\n", argv[0]);
		proc("network.ini");
		//proc_csv("D:\\gnss\\rtklib\\app\\consapp\\proc\\2022-10-02-obs.log");
	}
	else 
	{
		proc(argv[1]);
	}
	return 0;
}