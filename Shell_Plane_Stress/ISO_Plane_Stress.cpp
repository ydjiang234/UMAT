#ifndef UMATCPP
#define UMATCPP
#include <iostream>
#include <cmath>
#include <Dense>
//#include <aba_for_c.h>
using namespace std;
using namespace Eigen;

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
    MatrixXd Stiff;
    Vector3d Stress, dStrain;
    Stress = Map<Vector3d>(stress);
    dStrain = Map<Vector3d>(dstran);
    //cout<<Stress<<endl;
    

    Stiff = MatrixXd::Zero(3,3);
    E = props[0];
    anu = props[1];
    para1 = E / (1.0 - pow(anu,2));

    Stiff(0,0) = para1;
    Stiff(1,1) = para1;
    Stiff(2,2) = para1 * (1.0 - anu) / 2.0;
    Stiff(0,1) = para1 * anu;
    Stiff(1,0) = para1 * anu;

    Stress += Stiff * dStrain;
    for (int i=0;i<*ntens;i++) {
        for (int j=0;j<*ntens;j++) {
            ddsdde[i*3+j] = Stiff(i,j);
        }
        stress[i] = Stress(i);
        dstran[i] = dStrain(i);
    }
}
#endif
