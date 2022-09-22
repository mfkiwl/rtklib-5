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

#ifndef MAX_BUF_LEN
#define MAX_BUF_LEN 4096
#endif

typedef struct
{
	uint16_t type;
	uint32_t count;
}type_t;

typedef struct
{
	gtime_t time;
	uint16_t staid;
	double pos[3];
}coord_t;

typedef struct
{
	uint8_t data[MAX_BUF_LEN];
	uint16_t nbyte;
	uint16_t nlen;
	uint16_t type;
	uint32_t numofcrc;
	uint32_t numofmsg;
	uint32_t numofepo;
	uint32_t numofpos;
	uint32_t numofeph[6];/* GPS,GLO,GAL,BDS,QZS,unknown*/
}buff_t;

static rtcm_t rtcm_rove;

static int qc(const char* rovefname)
{
	FILE* fROV = fopen(rovefname, "rb");
	FILE* fOUT = NULL;
	int ret_rove = 0, data = 0, sys = 0, prn = 0, i = 0, j = 0;
	buff_t status = { 0 };
	std::vector<type_t> rtcm_type;
	std::vector< coord_t> coords;
	std::vector<obsd_t> obs;
	while (fROV && !feof(fROV) && (data = fgetc(fROV)) != EOF)
	{
		ret_rove = input_rtcm3(&rtcm_rove, data);
		if (rtcm_rove.type > 0)
		{
			int wk;
			double ws = time2gpst(rtcm_rove.time, &wk);
			printf("%10.3f,%4i,%4i,%i,%i\n", ws, rtcm_rove.type, rtcm_rove.staid, rtcm_rove.sync, rtcm_rove.crc);
		}
		if (ret_rove == 5)
		{
			/* pos 1005/1006 */
			++status.numofpos;
		}
		else if (ret_rove == 2)
		{
			/* eph 1019,1020,1045/1046,1042,1041 */
			sys = satsys(rtcm_rove.ephsat, &prn);
			if (sys==SYS_GPS)
				++status.numofeph[0];
			else if (sys==SYS_GLO)
				++status.numofeph[1];
			else if (sys == SYS_GAL)
				++status.numofeph[2];
			else if (sys == SYS_CMP)
				++status.numofeph[3];
			else if (sys == SYS_QZS)
				++status.numofeph[4];
			else
				++status.numofeph[5];
		}
		else if (ret_rove == 1)
		{
			/* obs MSM4,MSM7 */
			for (i = 0; i < rtcm_rove.obs.n; ++i)
			{
				obs.push_back(rtcm_rove.obs.data[i]);
			}
		}
	}
	if (fROV) fclose(fROV);
	if (fOUT) fclose(fOUT);
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
	int ret_rtcm_rove = init_rtcm(&rtcm_rove);

	std::time_t t = std::time(0);   // get time now
	std::tm* now = std::localtime(&t);
	int year = (now->tm_year + 1900);
	int month = (now->tm_mon + 1);
	int day = now->tm_mday;
	double ep[6] = { year, month, day, 0, 0,0 };

	if (argc < 2) /* 1 */
	{
		/* */
		printf("%s rtcmfname yyy-mm-day\n", argv[0]);
	}
	else if (argc < 3) /* 2 */
	{
		/* with input file name only */
		rtcm_rove.time = rtcm_rove.time_s = epoch2time(ep);
		qc(argv[1]);
	}
	else
	{
		int num = sscanf(argv[2], "%i-%i-%i", &year, &month, &day);
		if (num>=3)
		{
			ep[0] = year;
			ep[1] = month;
			ep[2] = day;
		}
		rtcm_rove.time = rtcm_rove.time_s = epoch2time(ep);
		qc(argv[1]);
	}

	free_rtcm(&rtcm_rove);
	return 0;
}