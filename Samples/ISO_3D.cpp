#ifndef UMATCPP
#define UMATCPP
#include <iostream>
#include <aba_for_c.h>

void test1(double *stress, double *statev, double *ddsdde, double *sse, double *spd, double *scd,
           double *rpl, double *ddsddt, double *drplde, double *drpldt,
           double *stran, double *dstran, double *time, double *dtime, double *temp, double *dtemp, double *predef, double *dpred, char *cmname,
           int *ndi, int *nshr, int *ntens, int *nstatv, double *props, int *nprops, double *coords, double *drot, double *pnewdt,
           double *celent, double *dfgrd0, double *dfgrd1, int *noel, int *npt, int *layer, int *kspt, int *kstep, int *kinc);

extern "C" void FOR_NAME(umat)(double *stress, double *statev, double *ddsdde, double *sse, double *spd, double *scd,
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
    double E, anu, para1;
    E = props[0];
    anu = props[1];
    para1 = E / (1.0 + anu) /(1.0 - 2.0 * anu);

    for (int i=0; i<*ntens; i++) {
        for (int j=0; j<*ntens; j++) {
            ddsdde[i*(*ntens)+j] = 0.0;
        }
    }
    
    for (int i=0; i<3; i++) {
        ddsdde[i*(*ntens)+i] = para1 * (1.0 - anu);
    }

    for (int i=3; i<6; i++) {
        ddsdde[i*(*ntens)+i] = para1 * (1.0 -2.0 * anu);
    }

    ddsdde[0*(*ntens)+1] = para1 * anu;
    ddsdde[0*(*ntens)+2] = para1 * anu;
    ddsdde[1*(*ntens)+2] = para1 * anu;
    ddsdde[1*(*ntens)+0] = para1 * anu;
    ddsdde[2*(*ntens)+0] = para1 * anu;
    ddsdde[2*(*ntens)+1] = para1 * anu;
    //std::cout<<sizeof(ddsdde[0])<<std::endl;


    for (int i=0; i<*ntens; i++) {
        for (int j=0; j<*ntens; j++) {
            stress[i] = stress[i] + ddsdde[i*(*ntens)+j] * dstran[j];
        }
    }
}
#endif
