#include<vector>
#include<cmath>
#include<iostream>
#include<complex>
using namespace std;

class sl_parameters{
public:
    double gamma;
    double omega;
    double lambda;
    double delay_no;
    double G;

    sl_parameters(double SLomega, double SLgamma, int SLlambda, double SLdelay_no, double SLG){
        omega = SLomega;
        gamma = SLgamma;
        lambda = SLlambda;
        delay_no = SLdelay_no;
        G = SLG;
    }
};

class fb_parameters{
public:
    complex<double> ephi;
    double K;

    fb_parameters(double FBephi, double FBK){
        complex<double> I(0,1);
        ephi = exp(I*FBephi);
        K = FBK;
    }
};

void deriv(complex<double> & z, vector<complex<double> >  & zdelay, double & J, sl_parameters *SL, vector<fb_parameters*> &FP, complex<double> & dznew)
{
    complex<double> I(0.0, 1.0);
    dznew = (SL->lambda+SL->G*J+I*SL->omega+SL->gamma*(z.real()*z.real()+z.imag()*z.imag()))*z;
    for( int i =0; i<SL->delay_no; i++){
        dznew+=FP[i]->K*FP[i]->ephi*zdelay[i];
    }
    //return dznew;
    //dznew = (SL->lambda+I*SL->omega-(1.0+I*SL->gamma)*(z.real()*z.real()+z.imag()*z.imag()))*z+FP->K1*zdelay1+FP->K2*zdelay2;
}

void mid_step(complex<double> & older_state, complex<double> & older_deriv, complex<double> & newer_state, complex<double> & newer_deriv, complex<double> & ms_new)
{
    ms_new= (1 + 2. * (0.5 - 0.0)/(1.0 - 0.0) ) * ((1.0 - 0.5)/(1.0 - 0.0))*((1.0 - 0.5)/(1.0 - 0.0)) * older_state
            +(0.5 - 0.0) * ((1.0 - 0.5)/(1.0 - 0.0))*((1.0 - 0.5)/(1.0 - 0.0)) * older_deriv
            +(1 + 2. * (0.5 - 1.0)/(0.0 - 1.0) ) * ((0.0 - 0.5)/(0.0 - 1.0)) * ((0.0 - 0.5)/(0.0 - 1.0)) * newer_state
            + (0.5 - 1.0) * ((0.0 - 0.5)/(0.0 - 1.0)) * ((0.0 - 0.5)/(0.0 - 1.0)) * newer_deriv;
}

complex<double> rk_step(complex<double> & z,vector<complex<double> >  & zdelay, vector<complex<double> >  & zdelayplus1, complex<double> & dz, vector<complex<double> > & dzdelay, vector<complex<double> > & dzdelayplus1, double & J, sl_parameters *SL, vector<fb_parameters*> & FP, double & h)
{
    complex<double> dznew;

    vector<complex<double> > ms_new(SL->delay_no,0);
    for( int i =0; i<SL->delay_no; i++){
        mid_step(zdelay[i], dzdelay[i], zdelayplus1[i], dzdelayplus1[i], ms_new[i]);
    }
    //cout << "mid_step_1_" << ms1_new << "mid_step_2_" << ms2_new << '\n';
    deriv(z, zdelay, J, SL, FP, dznew);
    dz=dznew;
    complex<double> k1= h*dznew;

    complex<double> dummy=z+k1/2.0;
    deriv(dummy, ms_new, J, SL, FP, dznew);
    complex<double> k2= h*dznew;

    dummy=z+k2/2.0;
    deriv(dummy, ms_new, J, SL, FP, dznew);
    complex<double> k3= h*dznew;

    dummy=z+k3;
    deriv(dummy, zdelayplus1, J, SL, FP, dznew);
    complex<double> k4= h*dznew;

    complex<double> znew=z+1.0/6.0*(k1+2.0*k2+2.0*k3+k4);

    return znew;


}



/*int main () {
    sl_parameters TestSL(1,-0.1,0.6,2176,0.01);
    sl_parameters *p = &TestSL;
    fb_parameters TestFB(3.14159,0.1);
    fb_parameters *p1 = &TestFB;
    vector<fb_parameters*> Vec_TestFB(2176,p1);

    vector<complex<double>> zdelay(2176, 0.1);
    vector<complex<double>> zdelayplus1(2176, 0.1);
    vector<complex<double>> dzdelay(2176, 0.1);
    vector<complex<double>> dzdelayplus1(2176, 0.1);
    complex<double> z(0.1,0);
    complex<double> dz(0,0);
    double h = 0.01;
    double J = 1;

    for (int i = 0; i < 10000; i++) {
        cout << rk_step(z, zdelay, zdelayplus1, dz, dzdelay, dzdelayplus1,  J, p, Vec_TestFB, h) << "\n";
        //z = rk_step(z, zdelay, zdelayplus1, dz, dzdelay, dzdelayplus1,  J, p, Vec_TestFB, h);
    }*/
//}
