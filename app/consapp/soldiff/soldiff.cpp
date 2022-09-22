#include <iostream>
#include <cmath>
#include <cstdio>
#include <vector>
#include <string>
#include <algorithm>

#ifdef _PLOT_
#include "..\matplotlib-cpp\matplotlibcpp.h"

namespace plt = matplotlibcpp;

static int g_plot = 0;
static double g_lat0 = 0;
static double g_lon0 = 0;
static double g_l2n = 0;
static double g_l2e = 0;
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

typedef struct
{
	double time;
	double lat;
	double lon;
	double ht;
	double vn;
	double ve;
	double vd;
	double roll;
	double pitch;
	double yaw;
	uint8_t nsat;
	uint8_t type;
}pvt_t;

#ifndef MAX_BUF_LEN
#define MAX_BUF_LEN 4096
#endif

typedef struct
{
	uint8_t dat[MAX_BUF_LEN];
	uint32_t nbyte, len;
}nmea_buff_t;

static int add_buff(nmea_buff_t* buff, uint8_t data)
{
	int ret = 0;
	if (buff->nbyte >= MAX_BUF_LEN) buff->nbyte = 0;
	if (buff->nbyte == 1 && data != 'G') buff->nbyte = 0;
	if (buff->nbyte == 0)
	{
		memset(buff, 0, sizeof(nmea_buff_t));
		if (data == '$')
		{
			buff->dat[buff->nbyte++] = data;
		}
	}
	else
	{
		buff->dat[buff->nbyte++] = data;
	}
	if (buff->nbyte > 2 && buff->dat[buff->nbyte - 1] == '\n' && buff->dat[buff->nbyte - 2] == '\r')
	{
		ret = 1;
		buff->len = buff->nbyte;
		buff->nbyte = 0;
	}
	else if (buff->nbyte > 4 && buff->dat[buff->nbyte - 1] == '\n' && buff->dat[buff->nbyte - 4] == '*')
	{
		ret = 1;
		buff->len = buff->nbyte;
		buff->nbyte = 0;
	}
	else if (buff->nbyte > 4 && buff->dat[buff->nbyte - 1] == '\r' && buff->dat[buff->nbyte - 4] == '*')
	{
		ret = 1;
		buff->len = buff->nbyte;
		buff->nbyte = 0;
	}
	return ret;
}

double lat2local(double lat, double* lat2north)
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

static int read_nmea(const char* fname, std::vector<pvt_t>& vpvt)
{
	FILE* fLOG = fopen(fname, "rb");
	if (!fLOG)
	{
		printf("cannot open %s\n", fname);
		return 0;
	}
	nmea_buff_t nmea_buff = { };
	char buffer[512] = { 0 };
	FILE* fCSV = NULL;
	unsigned long numofline = 0;
	unsigned long numofspp = 0;
	unsigned long numofrtk = 0;
	unsigned long numoffix = 0;
	unsigned long numofdgps = 0;
	char* val[MAXFIELD];
	int data = 0;
	while (fLOG != NULL && !feof(fLOG) && (data=fgetc(fLOG))!=EOF)
	{
		if (!add_buff(&nmea_buff, data)) continue;
		++numofline;
		int num = parse_fields((char*)nmea_buff.dat, val);
		if (num < 14 || strstr(val[0], "GGA") == NULL || strlen(val[3]) == 0 || strlen(val[5]) == 0) continue;
		pvt_t pvt = { 0 };
		pvt.time = atof(val[1]);
		int hh = pvt.time / 10000;
		pvt.time = pvt.time - hh * 10000;
		int mm = pvt.time / 100;
		pvt.time = pvt.time - mm * 100;
		pvt.time = hh * 3600 + mm * 60 + pvt.time + 18.0;
		pvt.type = atoi(val[6]);
		pvt.nsat = atoi(val[7]);
		pvt.lat = atof(val[2]);
		int deg = pvt.lat / 100.0;
		pvt.lat = (pvt.lat - deg * 100);
		pvt.lat = (deg + pvt.lat / 60.0) * D2R;
		if (val[3][0] == 'S' || val[3][0] == 's')
			pvt.lat = -pvt.lat;
		pvt.lon = atof(val[4]);
		deg = pvt.lon / 100.0;
		pvt.lon = (pvt.lon - deg * 100);
		pvt.lon = (deg + pvt.lon / 60.0) * D2R;
		if (val[5][0] == 'W' || val[5][0] == 'w')
			pvt.lon = -pvt.lon;
		pvt.ht = atof(val[9]) + atof(val[11]);
		if (pvt.type == 1)
			++numofspp;
		else if (pvt.type == 2)
			++numofdgps;
		else if (pvt.type == 4|| pvt.type == 5)
		{
			if (pvt.type == 4)
			{
				++numoffix;
			}
			++numofrtk;
		}
		vpvt.push_back(pvt);
		if (!fCSV) fCSV = set_output_file(fname, "-gps.csv");
		if (fCSV)
		{
			fprintf(fCSV, "%10.3f,%14.9f,%14.9f,%10.4f,%i,%i\n", pvt.time, pvt.lat * R2D, pvt.lon * R2D, pvt.ht, pvt.nsat, pvt.type);
		}
	}
	if (fLOG) fclose(fLOG);
	if (fCSV) fclose(fCSV);
	printf("%6i,%6i,%6i,%6i,%6i,%7.3f,%s\n", numofline, numofspp, numofdgps, numofrtk, numoffix, numofrtk > 0 ? 100.0 * numoffix / numofrtk : 0, fname);
	return vpvt.size();
}

static int sol_diff(const char* fname, std::vector<pvt_t>& pvt1, std::vector<pvt_t>& pvt2)
{
	int numofmatch = 0;
	int i = 0;
	int j = 0;
	int wd1 = 0;
	int wd2 = 0;
	double time1 = 0;
	double time2 = 0;
	FILE* fDIF = NULL;
	std::vector<pvt_t>::iterator p1 = pvt1.begin();
	std::vector<pvt_t>::iterator p2 = pvt2.begin();
	for (; p1 != pvt1.end(); ++p1)
	{
		wd1 = p1->time / (24 * 3600);
		time1 = p1->time - wd1 * 24 * 3600;
		while (p2 != pvt2.end())
		{
			wd2 = p2->time / (24 * 3600);
			time2 = p2->time - wd2 * 24 * 3600;
			if (time1 > (time2 + 0.005))
			{
				++p2;
				continue;
			}
			break;
		}
		if (p2 == pvt2.end()) break;
		if (time1 < (time2 - 0.005))
			continue;
		double l2n = 0, l2e = lat2local(p1->lat, &l2n);
		double dn = (p1->lat - p2->lat) * l2n;
		double de = (p1->lon - p2->lon) * l2e;
		double dh = (p1->ht - p2->ht);
		if (!fDIF) fDIF = set_output_file(fname, "-diff.csv");
		if (fDIF)
		{
			fprintf(fDIF, "%10.3f,%10.3f,%10.3f,%10.3f,%7.3f,%i,%i\n", p1->time, dn, de, dh, time1 - time2, p1->type, p2->type);
		}
	}
	if (fDIF) fclose(fDIF);
	return numofmatch;
}

static void plot_sol(std::vector<pvt_t>& pvt)
{
#ifdef _PLOT_
	size_t i = 0, numoffix = 0, numofsol = 0;
	double lat0_fix = 0, lon0_fix = 0, ht0_fix = 0;
	double lat0 = 0, lon0 = 0, ht0 = 0;
	for (; i < pvt.size(); ++i)
	{
		if (pvt[i].type == 0) continue;
		if (pvt[i].type == 4)
		{
			++numoffix;
			lat0_fix += pvt[i].lat;
			lon0_fix += pvt[i].lon;
			ht0_fix += pvt[i].ht;
		}
		++numofsol;
		lat0 += pvt[i].lat;
		lon0 += pvt[i].lon;
		ht0 += pvt[i].ht;
	}
	if (numoffix > 0)
	{
		lat0_fix /= numoffix;
		lon0_fix /= numoffix;
		ht0_fix /= numoffix;
		if (!g_plot)
		{
			g_l2e = lat2local(lat0_fix, &g_l2n);
			g_lat0 = lat0_fix;
			g_lon0 = lon0_fix;
			g_plot = 1;
		}
	}
	if (numofsol > 0)
	{
		lat0 /= numofsol;
		lon0 /= numofsol;
		ht0 /= numofsol;
		if (!g_plot)
		{
			g_l2e = lat2local(lat0, &g_l2n);
			g_lat0 = lat0;
			g_lon0 = lon0;
			g_plot = 1;
		}
	}
	if (g_plot)
	{
		std::vector<double> x, y, x_fix, y_fix;
		double sn = 0, se = 0, sn_fix = 0, se_fix = 0;
		numoffix = 0;
		for (i = 0; i < pvt.size(); ++i)
		{
			double de = (pvt[i].lon - g_lon0) * g_l2e;
			double dn = (pvt[i].lat - g_lat0) * g_l2n;
			sn += dn * dn;
			se += de * de;
			x.push_back(de);
			y.push_back(dn);
			if (pvt[i].type == 4)
			{
				x_fix.push_back(de);
				y_fix.push_back(dn);
				sn_fix += dn * dn;
				se_fix += de * de;
				++numoffix;
			}
		}
		if (numoffix > 0)
		{
			sn_fix = sqrt(sn_fix / numoffix);
			se_fix = sqrt(se_fix / numoffix);
		}
		if (pvt.size() > 0)
		{
			sn = sqrt(sn / pvt.size());
			se = sqrt(se / pvt.size());
		}

		char fix_title[255];
		char flt_title[255];
		sprintf(fix_title, "std=%7.3f,%7.3f FIX", sn_fix, se_fix);
		sprintf(flt_title, "std=%7.3f,%7.3f RTK", sn, se);
		plt::named_plot(flt_title, x, y, ".");
		plt::named_plot(fix_title, x_fix, y_fix, ".");
		plt::grid(true);
		plt::xlabel("East /Longitude [m]");
		plt::ylabel("North/Laitude   [m]");
		plt::title("Groud Track");
		plt::legend();
		plt::show();
	}
#endif
}

/* need to install matplotlib-cpp at upper directory https://github.com/yydgis/matplotlib-cpp.git */
/* need to install python and add the path => normally go to cmd, and use path command to locate the python path, for example C:\Users\xxx\AppData\Local\Programs\Python\Python310 */
/* need to install numpy and add the path for example, C:\Users\yudan\AppData\Local\Programs\Python\Python310\Lib\site-packages\numpy\core\include */
/* comment out some lines in matplotlibcpp.h about NPY_INT64, NPY_UINT64*/
/* add python310.lib in the path and dependencies, note this only works for release build, need to find the python310_d.lib for debug version */
int main(int argc, char** argv)
{
#ifdef _PLOT_
	g_plot = 0;
#endif
	//test_plot();
	std::vector<pvt_t> vpvt1;
	std::vector<pvt_t> vpvt2;
	if (argc < 3)
	{
		read_nmea(argv[1], vpvt1);
		plot_sol(vpvt1);
	}
	else
	{
		read_nmea(argv[2], vpvt2);
		plot_sol(vpvt2);
		sol_diff(argv[1], vpvt1, vpvt2);
		//plot_sol_diff(vpvt1,vpvt2);
	}
	return 0;
}