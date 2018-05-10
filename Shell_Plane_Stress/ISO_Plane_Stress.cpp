#ifndef UMATCPP
#define UMATCPP
#include <iostream>
#include <cmath>
//#include <aba_for_c.h>
using namespace std;

void test1(double *stress, double *statev, double *ddsdde, double *sse, double *spd, double *scd,
           double *rpl, double *ddsddt, double *drplde, double *drpldt,
           double *stran, double *dstran, double *time, double *dtime, double *temp, double *dtemp, double *predef, double *dpred, char *cmname,
           int *ndi, int *nshr, int *ntens, int *nstatv, double *props, int *nprops, double *coords, double *drot, double *pnewdt,
           double *celent, double *dfgrd0, double *dfgrd1, int *noel, int *npt, int *layer, int *kspt, int *kstep, int *kinc);

extern "C" void umat_(double *stress, double *statev, double *ddsdde, double *sse, double *spd, double *scd,
           double *rpl, double *ddsddt, double *drplde, double *drpldt,
           double *stran, double *dstran, double *time, double *dtime, double *temp, double *dtemp, double *predef, double *dpred, char *cmname,
           int *ndi, int *nshr, int *ntens, int *nstatv, double *props, int *nprops, double *coords, double *drot, double *pnewdt,
           double *celent, double *dfgrd0, double *dfgrd1, int *noel, int *npt, int *layer, int *kspt, int *kstep, int *kinc)
{
    test1(stress, statev, ddsdde, sse, spd, scd,
          rpl, ddsddt, drplde, drpldt,
          stran, dstran, time, dtime, temp, dtemp, predef, dpred, cmname,
           ndi, nshr, ntens, nstatv, props, nprops, coords, drot, pnewdt,
           celent, dfgrd0, dfgrd1, noel, npt, layer, kspt, kstep, kinc);
}

void test1(double *stress, double *statev, double *ddsdde, double *sse, double *spd, double *scd,
           double *rpl, double *ddsddt, double *drplde, double *drpldt,
           double *stran, double *dstran, double *time, double *dtime, double *temp, double *dtemp, double *predef, double *dpred, char *cmname,
           int *ndi, int *nshr, int *ntens, int *nstatv, double *props, int *nprops, double *coords, double *drot, double *pnewdt,
           double *celent, double *dfgrd0, double *dfgrd1, int *noel, int *npt, int *layer, int *kspt, int *kstep, int *kinc)
{
    //std::cout<<<<std::endl;
    double E, anu, para1;
    E = props[0];
    anu = props[1];
    para1 = E / (1.0 - pow(anu,2));

    for (int i=0; i<3; i++) {
        for (int j=0; j<3; j++) {
            ddsdde[i*3+j] = 0.0;
        }
    }
    
    ddsdde[0*3 + 0] = para1;
    ddsdde[1*3 + 1] = para1;
    ddsdde[2*3 + 2] = para1 * (1.0 - anu) / 2.0;
    ddsdde[0*3 + 1] = para1 * anu;
    ddsdde[1*3 + 0] = para1 * anu;



    for (int i=0; i<3; i++) {
        for (int j=0; j<3; j++) {
            stress[i] = stress[i] + ddsdde[i*3+j] * dstran[j];
        }
    }
}
#endif
