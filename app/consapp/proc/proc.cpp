#include <iostream>
#include <cmath>
#include <cstdio>
#include <vector>
#include <string>
#include <algorithm>
#include <ctime>

#include "rtklib.h"

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

static int read_log_file(const char* fname, std::vector<std::string>& logfnames, int *year, int *mon, int *day)
{
	FILE* fINI = fopen(fname, "r");
	char buffer[512] = { 0 };
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
			logfnames.push_back(val[1]);
		}
	}
	if (fINI) fclose(fINI);
	return (int)logfnames.size();
}

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
	double dt = 0, ws = 0;
	nav_t nav; /* global */
	/* open file, init and read the first epoch */
	j = 0;
	for (i = 0; i < numrcv; ++i)
	{
		rtcmfs[i] = fopen(logfnames[i].c_str(), "rb");
		init_rtcm(&rtcms[i]);
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
		ws = time2gpst(cur_time, &wk);
		printf("%10.4f,%6i\n", ws, numofepoch);
		for (i = 0; i < numrcv; ++i)
		{
			ws = time2gpst(rtcms[i].time, &wk);
			if (fabs(dt = timediff(cur_time, rtcms[i].time)) < 0.001)
			{
				++j;
				/* convert data and add to engine */
				printf("%10.4f,%3i,%3i,*\n", ws, rtcms[i].obs.n, i);
			}
			else
			{
				printf("%10.4f,%3i,%3i,-\n", ws, rtcms[i].obs.n, i);
			}
		}
		/* data process */
		if (j == 0) break; /* no more data */
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
	/* close files */
	for (i = 0; i < numrcv; ++i)
	{
		free_rtcm(&rtcms[i]);
		if (rtcmfs[i]) fclose(rtcmfs[i]);
	}
	return 0;
}

/* need to install matplotlib-cpp at upper directory https://github.com/yydgis/matplotlib-cpp.git */
/* need to install python and add the path => normally go to cmd, and use path command to locate the python path, for example C:\Users\xxx\AppData\Local\Programs\Python\Python310 */
/* need to install numpy and add the path for example, C:\Users\yudan\AppData\Local\Programs\Python\Python310\Lib\site-packages\numpy\core\include */
/* comment out some lines in matplotlibcpp.h about NPY_INT64, NPY_UINT64*/
/* add python310.lib in the path and dependencies, note this only works for release build, need to find the python310_d.lib for debug version */

/* input file will be rtcm file, please use the convbin to convert other format (rinex, ublox, septentrio, hemisphere, swiftnav, etc.) to rtcm */

int main(int argc, char** argv)
{
	if (argc < 2) /* 1 */
	{
		/* */
		printf("%s inifname\n", argv[0]);
	}
	else 
	{
		proc(argv[1]);
	}
	return 0;
}