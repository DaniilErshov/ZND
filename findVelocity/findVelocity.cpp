// findVelocity.cpp : Этот файл содержит функцию "main". Здесь начинается и заканчивается выполнение программы.
//
#include <consts.cpp>
#include <string>
#include <cmath>
#include <fstream>
#include <iostream>
#include <cvode/cvode.h>               /* prototypes for CVODE fcts., consts.  */
#include <nvector/nvector_serial.h>    /* access to serial N_Vector            */
#include <sunmatrix/sunmatrix_dense.h> /* access to dense SUNMatrix            */
#include <sunlinsol/sunlinsol_dense.h> /* access to dense SUNLinearSolver      */
using namespace std;

static int num_gas_species = 9;
static int num_react = 22;

#define Ith(v,i)    NV_Ith_S(v,i-1)         /* i-th vector component i=1..NEQ */
#define IJth(A,i,j) SM_ELEMENT_D(A,i-1,j-1) /* (i,j)-th matrix component i,j=1..NEQ */
#define NEQ   10               /* number of equations  */
#define RTOL  RCONST(1.0e-9)
#define RTOL2  RCONST(1.0e-7)/* scalar relative tolerance            */
#define ATOL1 RCONST(1.0e-11)   /* vector absolute tolerance components */
#define ATOL2 RCONST(1.0e-8)   /* vector absolute tolerance components */
#define T0    RCONST(0.0)      /* initial time           */
#define T1    RCONST(pow(10, -6))     /* first output time      */
#define TMULT RCONST(1.006)     /* output time factor     */
#define NOUT  1800
#define ZERO  RCONST(0.0)

ofstream fout;
double sigma;
double eta;
double M2;
double Tcurr;
double gamma;
double koeff_topl = 0.0;


// Number of chemical reactions
const int N_X = 100;
const int Nk = 20;

static int check_retval(void* returnvalue, const char* funcname, int opt)
{
    int* retval;

    /* Check if SUNDIALS function returned NULL pointer - no memory allocated */
    if (opt == 0 && returnvalue == NULL) {
        fprintf(stderr, "\nSUNDIALS_ERROR: %s() failed - returned NULL pointer\n\n",
            funcname);
        return(1);
    }

    /* Check if retval < 0 */
    else if (opt == 1) {
        retval = (int*)returnvalue;
        if (*retval < 0) {
            fprintf(stderr, "\nSUNDIALS_ERROR: %s() failed with retval = %d\n\n",
                funcname, *retval);
            return(1);
        }
    }

    /* Check if function returned NULL pointer - no memory allocated */
    else if (opt == 2 && returnvalue == NULL) {
        fprintf(stderr, "\nMEMORY_ERROR: %s() failed - returned NULL pointer\n\n",
            funcname);
        return(1);
    }

    return(0);
}

static int init_right_part_H2_multiple_reaction_22(realtype t, N_Vector y, N_Vector ydot, void* user_data) {
    // chemistry using a 22 reversible reactions kinetics for H2 burning from
    // Alan Keromnes et al. An experimental and detailed chemical kinetic modeling study of hydrogen and syngas mixture oxidation at elevated pressures // Combustion and Flame 160 (2013) 995–1011

    double k_0_f[3], k_inf_f[3], k_0_r[3], k_inf_r[3];
    double c[3], m[3], d = 0.14;
    double Pr_f[3], Pr_r[3];
    int k = 0, l = 0;
    double logF_f, logF_core_f, logF_r, logF_core_r;
    double* ytmp = new double[num_gas_species];
    double sum1, sum2;


    double* forward = new double[num_react];
    double* reverse = new double[num_react];
    double* equilib = new double[num_react];

    for (int i = 0; i < num_react; i++) {
        if (i != 8 && i != 15 && i != 16) {
            forward[i] = chec.kPrex_f[i] * pow(Ith(y, 10), chec.kPow_f[i])
                * exp(-chec.kE_f[i] / Ith(y, 10) / phyc.kRc);
            reverse[i] = chec.kPrex_r[i] * pow(Ith(y, 10), chec.kPow_r[i])
                * exp(-chec.kE_r[i] / Ith(y, 10) / phyc.kRc);
        }
        else {
            if (i == 8) k = 0;
            if (i == 15) k = 1;
            if (i == 16) k = 2;

            k_inf_f[k] = chec.kPrex_f[i] * pow(Ith(y, 10), chec.kPow_f[i])
                * exp(-chec.kE_f[i] / Ith(y, 10) / phyc.kRc);
            k_inf_r[k] = chec.kPrex_r[i] * pow(Ith(y, 10), chec.kPow_r[i])
                * exp(-chec.kE_r[i] / Ith(y, 10) / phyc.kRc);
            k_0_f[k] = chec.kPrex_f_lp[i] * pow(Ith(y, 10), chec.kPow_f_lp[i])
                * exp(-chec.kE_f_lp[i] / Ith(y, 10) / phyc.kRc);
            k_0_r[k] = chec.kPrex_r_lp[i] * pow(Ith(y, 10), chec.kPow_r_lp[i])
                * exp(-chec.kE_r_lp[i] / Ith(y, 10) / phyc.kRc);

            c[k] = -0.4 - 0.67 * log10(chec.Fcent[i]);
            m[k] = 0.75 - 1.27 * log10(chec.Fcent[i]);
        }
    }

    Pr_f[0] = (k_0_f[0] * (1.3 * Ith(y, 1) + 10 * Ith(y, 7) + Ith(y, 2) + Ith(y, 3) + Ith(y, 4) + Ith(y, 5) + Ith(y, 6) + Ith(y, 8) + Ith(y, 9))) / k_inf_f[0];
    Pr_r[0] = (k_0_r[0] * (1.3 * Ith(y, 1) + 10 * Ith(y, 7) + Ith(y, 2) + Ith(y, 3) + Ith(y, 4) + Ith(y, 5) + Ith(y, 6) + Ith(y, 8) + Ith(y, 9))) / k_inf_r[0];
    Pr_f[1] = (k_0_f[1] * (3.7 * Ith(y, 1) + 1.5 * Ith(y, 9) + 1.2 * Ith(y, 3) + 7.7 * Ith(y, 8) + Ith(y, 2) + Ith(y, 4) + Ith(y, 5) + Ith(y, 6) + Ith(y, 7))) / k_inf_f[1];
    Pr_r[1] = (k_0_r[1] * (3.7 * Ith(y, 1) + 1.5 * Ith(y, 9) + 1.2 * Ith(y, 3) + 7.7 * Ith(y, 8) + Ith(y, 2) + Ith(y, 4) + Ith(y, 5) + Ith(y, 6) + Ith(y, 7))) / k_inf_r[1];
    Pr_f[2] = (k_0_f[2] * Ith(y, 7)) / k_inf_f[2];
    Pr_r[2] = (k_0_r[2] * Ith(y, 7)) / k_inf_r[2];

    for (k = 0; k < 3; k++) {
        if (k == 0) l = 8;
        if (k == 1) l = 15;
        if (k == 2) l = 16;

        if (Pr_f[k] == 0) forward[l] = k_inf_f[k];
        else {
            logF_core_f = pow((log10(Pr_f[k]) + c[k]) / (m[k] - d * (log10(Pr_f[k]) + c[k])), 2);
            logF_f = pow(1.0 + logF_core_f, -1) * log10(chec.Fcent[l]);
            chec.F_f[l] = pow(10, logF_f);
            forward[l] = k_inf_f[k] * (Pr_f[k] / (1 + Pr_f[k])) * chec.F_f[l];
        }

        if (Pr_r[k] == 0) reverse[l] = k_inf_r[k];
        else {
            logF_core_r = pow((log10(Pr_r[k]) + c[k]) / (m[k] - d * (log10(Pr_r[k]) + c[k])), 2);
            logF_r = pow(1.0 + logF_core_r, -1) * log10(chec.Fcent[l]);
            chec.F_r[l] = pow(10, logF_r);
            reverse[l] = k_inf_r[k] * (Pr_r[k] / (1 + Pr_r[k])) * chec.F_r[l];
        }
    }

    equilib[0] = forward[0] * Ith(y, 2) * Ith(y, 3) - reverse[0] * Ith(y, 4) * Ith(y, 5);
    equilib[1] = forward[1] * Ith(y, 1) * Ith(y, 4) - reverse[1] * Ith(y, 2) * Ith(y, 5);
    equilib[2] = forward[2] * Ith(y, 1) * Ith(y, 5) - reverse[2] * Ith(y, 2) * Ith(y, 7);
    equilib[3] = forward[3] * Ith(y, 7) * Ith(y, 4) - reverse[3] * Ith(y, 5) * Ith(y, 5);
    equilib[4] = (forward[4] * Ith(y, 1) - reverse[4] * Ith(y, 2) * Ith(y, 2)) *
        (2.5 * Ith(y, 1) + 12 * Ith(y, 7) + Ith(y, 2) + Ith(y, 3) + Ith(y, 4) + Ith(y, 5) + Ith(y, 6) + Ith(y, 8) + Ith(y, 9));
    equilib[5] = (forward[5] * Ith(y, 4) * Ith(y, 4) - reverse[5] * Ith(y, 3)) *
        (2.5 * Ith(y, 1) + 12 * Ith(y, 7) + Ith(y, 2) + Ith(y, 3) + Ith(y, 4) + Ith(y, 5) + Ith(y, 6) + Ith(y, 8) + Ith(y, 9));
    equilib[6] = (forward[6] * Ith(y, 4) * Ith(y, 2) - reverse[6] * Ith(y, 5)) *
        (2.5 * Ith(y, 1) + 12 * Ith(y, 7) + Ith(y, 2) + Ith(y, 3) + Ith(y, 4) + Ith(y, 5) + Ith(y, 6) + Ith(y, 8) + Ith(y, 9));
    equilib[7] = (forward[7] * Ith(y, 2) * Ith(y, 5) - reverse[7] * Ith(y, 7)) *
        (0.73 * Ith(y, 1) + 3.65 * Ith(y, 7) + Ith(y, 2) + Ith(y, 3) + Ith(y, 4) + Ith(y, 5) + Ith(y, 6) + Ith(y, 8) + Ith(y, 9));
    equilib[8] = forward[8] * Ith(y, 2) * Ith(y, 3) - reverse[8] * Ith(y, 6);
    equilib[9] = forward[9] * Ith(y, 1) * Ith(y, 3) - reverse[9] * Ith(y, 2) * Ith(y, 6);
    equilib[10] = forward[10] * Ith(y, 6) * Ith(y, 2) - reverse[10] * Ith(y, 5) * Ith(y, 5);
    equilib[11] = forward[11] * Ith(y, 6) * Ith(y, 4) - reverse[11] * Ith(y, 5) * Ith(y, 3);
    equilib[12] = forward[12] * Ith(y, 6) * Ith(y, 5) - reverse[12] * Ith(y, 7) * Ith(y, 3);
    equilib[13] = forward[13] * Ith(y, 6) * Ith(y, 6) - reverse[13] * Ith(y, 8) * Ith(y, 3);
    equilib[14] = forward[14] * Ith(y, 6) * Ith(y, 6) - reverse[14] * Ith(y, 8) * Ith(y, 3);
    equilib[15] = forward[15] * Ith(y, 8) - reverse[15] * Ith(y, 5) * Ith(y, 5);
    equilib[16] = forward[16] * Ith(y, 8) - reverse[16] * Ith(y, 5) * Ith(y, 5);
    equilib[17] = forward[17] * Ith(y, 8) * Ith(y, 2) - reverse[17] * Ith(y, 7) * Ith(y, 5);
    equilib[18] = forward[18] * Ith(y, 8) * Ith(y, 2) - reverse[18] * Ith(y, 6) * Ith(y, 1);
    equilib[19] = forward[19] * Ith(y, 8) * Ith(y, 4) - reverse[19] * Ith(y, 5) * Ith(y, 6);
    equilib[20] = forward[20] * Ith(y, 8) * Ith(y, 5) - reverse[20] * Ith(y, 6) * Ith(y, 7);
    equilib[21] = forward[21] * Ith(y, 8) * Ith(y, 5) - reverse[21] * Ith(y, 6) * Ith(y, 7);


    Ith(ydot, 1) = -equilib[1] - equilib[2] - equilib[4] - equilib[9] + equilib[18];
    Ith(ydot, 2) = -equilib[0] + equilib[4] - equilib[6] - equilib[7] - equilib[8] -
        equilib[10] - equilib[17] - Ith(ydot, 1);
    Ith(ydot, 3) = -equilib[0] + equilib[5] - equilib[8] - equilib[9] + equilib[11] +
        equilib[12] + equilib[13] + equilib[14];
    Ith(ydot, 4) = equilib[0] - equilib[1] - equilib[3] - 2. * equilib[5] - equilib[6] -
        equilib[11] - equilib[19];
    Ith(ydot, 5) = equilib[0] + equilib[1] - equilib[2] + 2. * equilib[3] + equilib[6] -
        equilib[7] + 2. * equilib[10] + equilib[11] - equilib[12] + 2. * equilib[15] +
        2. * equilib[16] + equilib[17] + equilib[19] - equilib[20] - equilib[21];
    Ith(ydot, 6) = equilib[8] + equilib[9] - equilib[10] - equilib[11] - equilib[12] -
        2. * equilib[13] - 2. * equilib[14] + equilib[18] + equilib[19] + equilib[20] +
        equilib[21];
    Ith(ydot, 7) = equilib[2] - equilib[3] + equilib[7] + equilib[12] + equilib[17] +
      equilib[20] + equilib[21];
    Ith(ydot, 8) = equilib[13] + equilib[14] - equilib[15] - equilib[16] - equilib[17] -
        equilib[18] - equilib[19] - equilib[20] - equilib[21];
    //Ith(ydot, 7) = (- phyc.mol_weight[0] * Ith(ydot, 1) - phyc.mol_weight[1] * Ith(ydot, 2) - phyc.mol_weight[2] * Ith(ydot, 3) - phyc.mol_weight[3] * Ith(ydot, 4) - phyc.mol_weight[4] * Ith(ydot, 5) - phyc.mol_weight[5] * Ith(ydot, 6) - phyc.mol_weight[7] * Ith(ydot, 8))/ phyc.mol_weight[6];
    Ith(ydot, 9) = 0;
    Ith(ydot, 10) = 0;



    delete[] reverse;
    delete[] equilib;
    delete[] ytmp;
    return 0;
}

void set_initial_concentration(N_Vector y, double* Temperature, double* rho2, double* Y0)
{
    double* ytmp = new double[num_gas_species];
    double konc;
    double checkrho = 0;
    double Msmes = 0.0;
    double rho = 0.0;
    //cout << "konc = " << konc << endl;
    //Mol fractions
    Ith(y, 1) = Y0[0]; //H2
    Ith(y, 2) = Y0[1];
    Ith(y, 3) = Y0[2];  //O2
    Ith(y, 4) = Y0[3];
    Ith(y, 5) = Y0[4]; //OH
    Ith(y, 6) = Y0[5];
    Ith(y, 7) = Y0[6];
    Ith(y, 8) = Y0[7];
    Ith(y, 9) = Y0[8];    //N2
    // Initial temperature
    Ith(y, 10) = *Temperature;

    for (int i = 0; i < num_gas_species; i++)
    {
        Msmes += Ith(y, i + 1) * phyc.mol_weight[i];
    }
    //Mass fractions
    for (int i = 0; i < num_gas_species; i++)
    {
        Ith(y, i + 1) = Ith(y, i + 1) * phyc.mol_weight[i] / Msmes;
        ytmp[i] = Ith(y, i + 1);
    }
    //MOL. CONCENTR.
    for (int i = 0; i < num_gas_species; i++)
    {
        Ith(y, i + 1) = Ith(y, i + 1) * (*rho2) / phyc.mol_weight[i];
        if (fabs(Ith(y, i + 1)) < ATOL1 && Ith(y, i + 1) < 0)
        {
            Ith(y, i + 1) = 0;
        }
    }
    delete[] ytmp;
    return;
}



void integrate(double* Temperature, double* Y2, double* rho2, double* Y0, N_Vector y, void *cvode_mem)
{
    double t, tout;
    int i;
    double checkrho = 0;
    int retval, iout;
    t = 0.0;

    iout = 0;  tout = T1;
    while (1) {
        retval = CVode(cvode_mem, tout, y, &t, CV_NORMAL);
        if (check_retval(&retval, "CVode", 1)) break;
        if (retval == CV_SUCCESS) {
            iout++;
            tout *= TMULT;
        }
        if (iout == NOUT) break;
    }
    for (i = 0; i < num_gas_species; i++)
    {
        checkrho += Ith(y, i + 1) * phyc.mol_weight[i];
    }
    //mass fraction
    for (i = 0; i < num_gas_species; i++)
    {
        Y2[i] = Ith(y, i + 1) * phyc.mol_weight[i]/(checkrho);
    }


}

double H_u_T(double u1, double h1, double Temperature, double X, double *Y2)
{
    double h2 = get_enthalpy(num_gas_species, Y2, Temperature);
    return (h2 - h1) + ((pow(u1/X, 2)/2) - (pow(u1, 2) / 2));
}

double P_u_T(double u1, double Temperature, double X, double* Y2, double p1, double rho1, double rho2)
{
    double R2 = get_gas_constant(num_gas_species, Y2);
    double p2 = X * rho1 * R2 * Temperature;
    return p2 - p1 + X * rho1 * pow(u1 / X, 2) - rho1 * pow(u1, 2);
}

double dHdu(double u1, double du1, double h1, double Temperature, double X, double* Y2)
{
    return (H_u_T(u1 + du1, h1, Temperature, X, Y2) - H_u_T(u1, h1, Temperature, X, Y2)) / du1;
}

double dHdT(double u1, double dT, double h1, double Temperature, double X, double* Y2, double* Y2_dT)
{
    return (H_u_T(u1, h1, Temperature + dT, X, Y2_dT) - H_u_T(u1, h1, Temperature, X, Y2)) / dT;
}

double dPdu(double u1, double du1, double Temperature, double X, double* Y2, double p1, double rho1, double rho2)
{
    return (P_u_T(u1+du1,Temperature, X, Y2, p1, rho1, rho2) - P_u_T(u1, Temperature, X, Y2, p1, rho1, rho2)) / du1;
}

double dPdT(double u1, double dT, double Temperature, double X, double* Y2, double p1, double rho1, double rho2, double* Y2_dT)
{
    return (P_u_T(u1, Temperature + dT, X, Y2_dT, p1, rho1, rho2) - P_u_T(u1, Temperature, X, Y2, p1, rho1, rho2)) / dT;
}

double detA(double u1, double du1, double T2, double dT, double X, double* Y2, double p1, double rho1, double rho2, double h1,double* Y2_dT)
{
    return dHdu(u1, du1, h1, T2, X, Y2) * dPdT(u1, dT, T2, X, Y2, p1, rho1, rho2, Y2_dT) - dHdT(u1, dT, h1, T2, X, Y2, Y2_dT) * dPdu(u1, du1, T2, X, Y2, p1, rho1, rho2);
}

void MatrixMul(double A[3][Nk], double B[Nk][3], double C[3][3])
{
    for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++)
        {
            C[i][j] = 0;
            for (int k = 0; k < Nk; k++)
                C[i][j] += A[i][k] * B[k][j];
        }
}

void Inverse_matrix(double A[3][3], double result[3][3])
{
    double  detA = A[0][0] * A[1][1] * A[2][2] + A[0][1] * A[1][2] * A[2][0]
        + A[1][0] * A[2][1] * A[0][2] - A[2][0] * A[1][1] * A[0][2] -
        A[2][1] * A[1][2] * A[0][0] - A[1][0] * A[0][1] * A[2][2];
    cout << "det = " << detA << endl;
    double invdet = 1. / detA;
    result[0][0] = (A[1][1] * A[2][2] - A[2][1] * A[1][2]) * invdet;
    result[0][1] = -(A[0][1] * A[2][2] - A[0][2] * A[2][1]) * invdet;
    result[0][2] = (A[0][1] * A[1][2] - A[0][2] * A[1][1]) * invdet;
    result[1][0] = -(A[1][0] * A[2][2] - A[1][2] * A[2][0]) * invdet;
    result[1][1] = (A[0][0] * A[2][2] - A[0][2] * A[2][0]) * invdet;
    result[1][2] = -(A[0][0] * A[1][2] - A[1][0] * A[0][2]) * invdet;
    result[2][0] = (A[1][0] * A[2][1] - A[2][0] * A[1][1]) * invdet;
    result[2][1] = -(A[0][0] * A[2][1] - A[2][0] * A[0][1]) * invdet;
    result[2][2] = (A[0][0] * A[1][1] - A[1][0] * A[0][1]) * invdet;

}


void Find_velocity(double*Y, double p1, double Tin, double &u, double &Xmin)
{

    double X_u1[N_X][2];
    //rho initial
    double rho1 = 0;
    double Msmes = 0;
    double Y0[9];
    for (int i = 0; i < num_gas_species; i++)
    {
        Y0[i] = Y[i];
    }
    for (int i = 0; i < num_gas_species; i++)
    {
        Msmes += Y0[i] * phyc.mol_weight[i];
    }
    //Mass fractions
    for (int i = 0; i < num_gas_species; i++)
    {
        Y0[i] = Y0[i] * phyc.mol_weight[i] / Msmes;
    }

    //initial data

    double R1 = get_gas_constant(num_gas_species, Y0);
    rho1 = p1 / R1 / Tin;
    double h1 = get_enthalpy(num_gas_species, Y0, Tin);
    int check_ret;


    //zero approximation
    double X = 1.6;
    double u1 = 0;
    double T2 = 0;
    double dx = (2.0 - 1.6) / (double)N_X;


    double* Y2 = new double[num_gas_species];
    double* Y_dT = new double[num_gas_species];
    int j;
    double rho2 = X * rho1;
    double H;
    double P;
    double u1n = 1.993;
    double T2n = 2650.0;
    double T2_dT = 2650.0;
    double du1 = 0.00001;
    double dT = 0.0001;
    double err = pow(10, -6);
    double umin = 100;
    int ind_umin;


    SUNContext sunctx;
    N_Vector y;
    N_Vector abstol;
    SUNMatrix A;
    SUNLinearSolver LS;
    int i;
    double checkrho = 0;
    void* cvode_mem;
    int retval, iout;
    y = NULL;
    abstol = NULL;
    A = NULL;
    LS = NULL;
    cvode_mem = NULL;
    retval = SUNContext_Create(NULL, &sunctx);
    if (check_retval(&retval, "SUNContext_Create", 1)) return;

    /* Initial conditions */
    y = N_VNew_Serial(NEQ, sunctx);
    if (check_retval((void*)y, "N_VNew_Serial", 0)) return;

    /* Initialize y */
    set_initial_concentration(y, &T2n, &rho2, Y0);
    /* Set the vector absolute tolerance */
    abstol = N_VNew_Serial(NEQ, sunctx);

    if (check_retval((void*)abstol, "N_VNew_Serial", 0)) return;
    Ith(abstol, 1) = ATOL1;
    Ith(abstol, 2) = ATOL1;
    Ith(abstol, 3) = ATOL1;
    Ith(abstol, 4) = ATOL1;
    Ith(abstol, 5) = ATOL1;
    Ith(abstol, 6) = ATOL1;
    Ith(abstol, 7) = ATOL1;
    Ith(abstol, 8) = ATOL1;
    Ith(abstol, 9) = ATOL1;
    Ith(abstol, 10) = ATOL1;

    /* Call CVodeCreate to create the solver memory and specify the
     * Backward Differentiation Formula */
    cvode_mem = CVodeCreate(CV_BDF, sunctx);
    if (check_retval((void*)cvode_mem, "CVodeCreate", 0)) return;

    /* Call CVodeInit to initialize the integrator memory and specify the
     * user's right hand side function in y'=f(t,y), the initial time T0, and
     * the initial dependent variable vector y. */
    retval = CVodeInit(cvode_mem, init_right_part_H2_multiple_reaction_22, T0, y);
    if (check_retval(&retval, "CVodeInit", 1)) return;

    /* Call CVodeSVtolerances to specify the scalar relative tolerance
     * and vector absolute tolerances */
    retval = CVodeSVtolerances(cvode_mem, RTOL, abstol);
    if (check_retval(&retval, "CVodeSVtolerances", 1)) return;

    /* Call CVodeRootInit to specify the root function g with 2 components */

    /* Create dense SUNMatrix for use in linear solves */
    A = SUNDenseMatrix(NEQ, NEQ, sunctx);
    if (check_retval((void*)A, "SUNDenseMatrix", 0)) return;

    /* Create dense SUNLinearSolver object for use by CVode */
    LS = SUNLinSol_Dense(y, A, sunctx);
    if (check_retval((void*)LS, "SUNLinSol_Dense", 0)) return;

    /* Attach the matrix and linear solver */
    retval = CVodeSetLinearSolver(cvode_mem, LS, A);
    if (check_retval(&retval, "CVodeSetLinearSolver", 1)) return;
    cout << " X   u1" << endl;



    for (int i = 0; i < N_X; i++)
    {
        X += dx;
        rho2 = X * rho1;
        u1n = 1.993;
        T2n = 2600.0;
        while (abs(u1n - u1) > err)
        {
            checkrho = 0;
            u1 = u1n;
            T2 = T2n;
            T2_dT = T2n + dT;
            integrate(&T2, Y2, &rho2, Y, y, cvode_mem);
            /*for (int i = 0; i < num_gas_species; i++)
            {
                cout << "Y2[i] = " << Y2[i] << endl;
            }*/

  

            H = H_u_T(u1, h1, T2, X, Y2);
            P = P_u_T(u1, T2, X, Y2, p1, rho1, rho2);

            set_initial_concentration(y, &T2_dT, &rho2, Y);
            retval = CVodeReInit(cvode_mem, T0, y);
            integrate(&T2_dT, Y_dT, &rho2, Y, y, cvode_mem);


            u1n = u1 - (dPdT(u1, dT, T2, X, Y_dT, p1, rho1, rho2, Y_dT) * H - dHdT(u1, dT, h1, T2, X, Y2, Y_dT) * P) / detA(u1, du1, T2, dT, X, Y2, p1, rho1, rho2, h1, Y_dT);
            T2n = T2 - (-dPdu(u1, du1, T2, X, Y2, p1, rho1, rho2) * H + dHdu(u1, du1, h1, T2, X, Y2) * P) / detA(u1, du1, T2, dT, X, Y2, p1, rho1, rho2, h1, Y_dT);
       
            set_initial_concentration(y, &T2n, &rho2, Y);
            retval = CVodeReInit(cvode_mem, T0, y);
        }
        X_u1[i][0] = X;
        X_u1[i][1] = u1n;
        if (umin > u1n)
        {
            umin = u1n;
            ind_umin = i;
        }
        cout << X_u1[i][0] << "  " << X_u1[i][1] << endl;
    }
    double Psi[3][Nk];
    double PsiT[Nk][3];
    double PsiU[3];
    cout << "ind_min = " << X_u1[ind_umin][0] << endl;
    cout << endl;
    cout << "PSi" << endl;
    for (int i = 0; i < 3; i++)
    {
        for (int j = 0; j < Nk; j++)
        {
            if (i == 0) Psi[i][j] = X_u1[ind_umin - Nk / 2 + j][0] * X_u1[ind_umin - Nk / 2 + j][0];
            if (i == 1) Psi[i][j] = X_u1[ind_umin - Nk / 2 + j][0];
            if (i == 2) Psi[i][j] = 1;
            cout << Psi[i][j] << " ";
        }
        cout << endl;

    }
    cout << "PSiT" << endl;
    for (int i = 0; i < Nk; i++)
    {
        for (int j = 0; j < 3; j++)
        {
            if (j == 0) PsiT[i][j] = X_u1[ind_umin - Nk / 2 + i][0] * X_u1[ind_umin - Nk / 2 + i][0];
            if (j == 1) PsiT[i][j] = X_u1[ind_umin - Nk / 2 + i][0];
            if (j == 2) PsiT[i][j] = 1;
            cout << PsiT[i][j] << " ";
        }
        cout << endl;
    }

    for (int i = 0; i < 3; i++)
    {
        PsiU[i] = 0;
        for (int j = 0; j < Nk; j++)
        {
            PsiU[i] += Psi[i][j] * X_u1[ind_umin - Nk / 2 + j][1];
        }
    }
    double PsiPsiT[3][3];
    double PsiPsiTinv[3][3];
    MatrixMul(Psi, PsiT, PsiPsiT);
    cout << endl;
    cout << "PsiPsiT" << endl;
    for (int i = 0; i < 3; i++)
    {
        for (int j = 0; j < 3; j++)
        {
            cout << PsiPsiT[i][j] << " ";
        }
        cout << endl;
    }

    cout << "inv" << endl;
    Inverse_matrix(PsiPsiT, PsiPsiTinv);
    for (int i = 0; i < 3; i++)
    {
        for (int j = 0; j < 3; j++)
        {
            cout << PsiPsiTinv[i][j] << " ";
        }
        cout << endl;
    }
    double a = PsiPsiTinv[0][0] * PsiU[0] + PsiPsiTinv[0][1] * PsiU[1] + PsiPsiTinv[0][2] * PsiU[2];
    double b = PsiPsiTinv[1][0] * PsiU[0] + PsiPsiTinv[1][1] * PsiU[1] + PsiPsiTinv[1][2] * PsiU[2];
    double c = PsiPsiTinv[2][0] * PsiU[0] + PsiPsiTinv[2][1] * PsiU[1] + PsiPsiTinv[2][2] * PsiU[2];
    u = (-b * b + 4. * a * c) / (4. * a);
    Xmin = -b / (2. * a);
    /* Free memory */
    N_VDestroy(y);                            /* Free y vector */
    N_VDestroy(abstol);                       /* Free abstol vector */
    CVodeFree(&cvode_mem);                    /* Free CVODE memory */
    SUNLinSolFree(LS);                        /* Free the linear solver memory */
    SUNMatDestroy(A);                         /* Free the matrix memory */
    SUNContext_Free(&sunctx);                 /* Free the SUNDIALS context */

}

double Func(double u, double u1, double P1, double rho1, double T, double* Y0)
{
    double I1 = rho1 * u1;
    double I2 = rho1 * u1 * u1 + P1;
    double I3 = get_enthalpy(num_gas_species, Y0, T) + u1 * u1 / 2.;
    double Tu = (u * I2 - u * u * I1) / (get_gas_constant(num_gas_species, Y0) * I1);
    cout << "Tu = " << Tu << endl;
    return get_enthalpy(num_gas_species, Y0, Tu) + u * u / 2. - I3;
}

double Func_prime(double u, double u1, double P1, double rho1, double T, double* Y0)
{
    double I1 = rho1 * u1;
    double I2 = rho1 * u1 * u1 + P1;
    double I3 = get_enthalpy(num_gas_species, Y0, T) + u1 * u1 / 2.;
    double Tu = (u * I2 - u * u * I1) / (get_gas_constant(num_gas_species, Y0) * I1);
    return (I2 - 2. * u * I1) / (get_gas_constant(num_gas_species, Y0) * I1) * get_Cp(num_gas_species, Y0, Tu)  + u;
}



static int init_right_part_ZND(realtype t, N_Vector y, N_Vector ydot, void* user_data) {
    // chemistry using a 22 reversible reactions kinetics for H2 burning from
    // Alan Keromnes et al. An experimental and detailed chemical kinetic modeling study of hydrogen and syngas mixture oxidation at elevated pressures // Combustion and Flame 160 (2013) 995–1011

    double k_0_f[3], k_inf_f[3], k_0_r[3], k_inf_r[3];
    double c[3], m[3], d = 0.14;
    double Pr_f[3], Pr_r[3];
    int k = 0, l = 0;
    double logF_f, logF_core_f, logF_r, logF_core_r;
    double sum1, sum2;
    sum1 = 0;

    for (int i = 0; i < num_gas_species; i++)
        sum1 += Ith(y, i + 1);
    //double* Ymassfr = new double[num_gas_species];

   // for (int i = 0; i < num_gas_species; i++)
     //   Ymassfr[i] = Ith(y, i + 1) * phyc.mol_weight[i] / Ith(y, num_gas_species + 2);

    Tcurr = Ith(y, num_gas_species + 1) / phyc.kR / sum1;

    //cout << "Tcurr = " << Tcurr << endl;
    double* forward = new double[num_react];
    double* reverse = new double[num_react];
    double* equilib = new double[num_react];

    for (int i = 0; i < num_react; i++) {
        if (i != 8 && i != 15 && i != 16) {
            forward[i] = chec.kPrex_f[i] * pow(Tcurr, chec.kPow_f[i])
                * exp(-chec.kE_f[i] / Tcurr / phyc.kRc);
            reverse[i] = chec.kPrex_r[i] * pow(Tcurr, chec.kPow_r[i])
                * exp(-chec.kE_r[i] / Tcurr / phyc.kRc);
        }
        else {
            if (i == 8) k = 0;
            if (i == 15) k = 1;
            if (i == 16) k = 2;

            k_inf_f[k] = chec.kPrex_f[i] * pow(Tcurr, chec.kPow_f[i])
                * exp(-chec.kE_f[i] / Tcurr / phyc.kRc);
            k_inf_r[k] = chec.kPrex_r[i] * pow(Tcurr, chec.kPow_r[i])
                * exp(-chec.kE_r[i] / Tcurr / phyc.kRc);
            k_0_f[k] = chec.kPrex_f_lp[i] * pow(Tcurr, chec.kPow_f_lp[i])
                * exp(-chec.kE_f_lp[i] / Tcurr / phyc.kRc);
            k_0_r[k] = chec.kPrex_r_lp[i] * pow(Tcurr, chec.kPow_r_lp[i])
                * exp(-chec.kE_r_lp[i] / Tcurr / phyc.kRc);

            c[k] = -0.4 - 0.67 * log10(chec.Fcent[i]);
            m[k] = 0.75 - 1.27 * log10(chec.Fcent[i]);
        }
    }

    Pr_f[0] = (k_0_f[0] * (1.3 * Ith(y, 1) + 10 * Ith(y, 7) + Ith(y, 2) + Ith(y, 3) + Ith(y, 4) + Ith(y, 5) + Ith(y, 6) + Ith(y, 8) + Ith(y, 9))) / k_inf_f[0];
    Pr_r[0] = (k_0_r[0] * (1.3 * Ith(y, 1) + 10 * Ith(y, 7) + Ith(y, 2) + Ith(y, 3) + Ith(y, 4) + Ith(y, 5) + Ith(y, 6) + Ith(y, 8) + Ith(y, 9))) / k_inf_r[0];
    Pr_f[1] = (k_0_f[1] * (3.7 * Ith(y, 1) + 1.5 * Ith(y, 9) + 1.2 * Ith(y, 3) + 7.7 * Ith(y, 8) + Ith(y, 2) + Ith(y, 4) + Ith(y, 5) + Ith(y, 6) + Ith(y, 7))) / k_inf_f[1];
    Pr_r[1] = (k_0_r[1] * (3.7 * Ith(y, 1) + 1.5 * Ith(y, 9) + 1.2 * Ith(y, 3) + 7.7 * Ith(y, 8) + Ith(y, 2) + Ith(y, 4) + Ith(y, 5) + Ith(y, 6) + Ith(y, 7))) / k_inf_r[1];
    Pr_f[2] = (k_0_f[2] * Ith(y, 7)) / k_inf_f[2];
    Pr_r[2] = (k_0_r[2] * Ith(y, 7)) / k_inf_r[2];

    for (k = 0; k < 3; k++) {
        if (k == 0) l = 8;
        if (k == 1) l = 15;
        if (k == 2) l = 16;

        if (Pr_f[k] == 0) forward[l] = k_inf_f[k];
        else {
            logF_core_f = pow((log10(Pr_f[k]) + c[k]) / (m[k] - d * (log10(Pr_f[k]) + c[k])), 2);
            logF_f = pow(1.0 + logF_core_f, -1) * log10(chec.Fcent[l]);
            chec.F_f[l] = pow(10, logF_f);
            forward[l] = k_inf_f[k] * (Pr_f[k] / (1 + Pr_f[k])) * chec.F_f[l];
        }

        if (Pr_r[k] == 0) reverse[l] = k_inf_r[k];
        else {
            logF_core_r = pow((log10(Pr_r[k]) + c[k]) / (m[k] - d * (log10(Pr_r[k]) + c[k])), 2);
            logF_r = pow(1.0 + logF_core_r, -1) * log10(chec.Fcent[l]);
            chec.F_r[l] = pow(10, logF_r);
            reverse[l] = k_inf_r[k] * (Pr_r[k] / (1 + Pr_r[k])) * chec.F_r[l];
        }
    }

    equilib[0] = forward[0] * Ith(y, 2) * Ith(y, 3) - reverse[0] * Ith(y, 4) * Ith(y, 5);
    equilib[1] = forward[1] * Ith(y, 1) * Ith(y, 4) - reverse[1] * Ith(y, 2) * Ith(y, 5);
    equilib[2] = forward[2] * Ith(y, 1) * Ith(y, 5) - reverse[2] * Ith(y, 2) * Ith(y, 7);
    equilib[3] = forward[3] * Ith(y, 7) * Ith(y, 4) - reverse[3] * Ith(y, 5) * Ith(y, 5);
    equilib[4] = (forward[4] * Ith(y, 1) - reverse[4] * Ith(y, 2) * Ith(y, 2)) *
        (2.5 * Ith(y, 1) + 12 * Ith(y, 7) + Ith(y, 2) + Ith(y, 3) + Ith(y, 4) + Ith(y, 5) + Ith(y, 6) + Ith(y, 8) + Ith(y, 9));
    equilib[5] = (forward[5] * Ith(y, 4) * Ith(y, 4) - reverse[5] * Ith(y, 3)) *
        (2.5 * Ith(y, 1) + 12 * Ith(y, 7) + Ith(y, 2) + Ith(y, 3) + Ith(y, 4) + Ith(y, 5) + Ith(y, 6) + Ith(y, 8) + Ith(y, 9));
    equilib[6] = (forward[6] * Ith(y, 4) * Ith(y, 2) - reverse[6] * Ith(y, 5)) *
        (2.5 * Ith(y, 1) + 12 * Ith(y, 7) + Ith(y, 2) + Ith(y, 3) + Ith(y, 4) + Ith(y, 5) + Ith(y, 6) + Ith(y, 8) + Ith(y, 9));
    equilib[7] = (forward[7] * Ith(y, 2) * Ith(y, 5) - reverse[7] * Ith(y, 7)) *
        (0.73 * Ith(y, 1) + 3.65 * Ith(y, 7) + Ith(y, 2) + Ith(y, 3) + Ith(y, 4) + Ith(y, 5) + Ith(y, 6) + Ith(y, 8) + Ith(y, 9));
    equilib[8] = forward[8] * Ith(y, 2) * Ith(y, 3) - reverse[8] * Ith(y, 6);
    equilib[9] = forward[9] * Ith(y, 1) * Ith(y, 3) - reverse[9] * Ith(y, 2) * Ith(y, 6);
    equilib[10] = forward[10] * Ith(y, 6) * Ith(y, 2) - reverse[10] * Ith(y, 5) * Ith(y, 5);
    equilib[11] = forward[11] * Ith(y, 6) * Ith(y, 4) - reverse[11] * Ith(y, 5) * Ith(y, 3);
    equilib[12] = forward[12] * Ith(y, 6) * Ith(y, 5) - reverse[12] * Ith(y, 7) * Ith(y, 3);
    equilib[13] = forward[13] * Ith(y, 6) * Ith(y, 6) - reverse[13] * Ith(y, 8) * Ith(y, 3);
    equilib[14] = forward[14] * Ith(y, 6) * Ith(y, 6) - reverse[14] * Ith(y, 8) * Ith(y, 3);
    equilib[15] = forward[15] * Ith(y, 8) - reverse[15] * Ith(y, 5) * Ith(y, 5);
    equilib[16] = forward[16] * Ith(y, 8) - reverse[16] * Ith(y, 5) * Ith(y, 5);
    equilib[17] = forward[17] * Ith(y, 8) * Ith(y, 2) - reverse[17] * Ith(y, 7) * Ith(y, 5);
    equilib[18] = forward[18] * Ith(y, 8) * Ith(y, 2) - reverse[18] * Ith(y, 6) * Ith(y, 1);
    equilib[19] = forward[19] * Ith(y, 8) * Ith(y, 4) - reverse[19] * Ith(y, 5) * Ith(y, 6);
    equilib[20] = forward[20] * Ith(y, 8) * Ith(y, 5) - reverse[20] * Ith(y, 6) * Ith(y, 7);
    equilib[21] = forward[21] * Ith(y, 8) * Ith(y, 5) - reverse[21] * Ith(y, 6) * Ith(y, 7);


    Ith(ydot, 1) = -equilib[1] - equilib[2] - equilib[4] - equilib[9] + equilib[18];
    Ith(ydot, 2) = -equilib[0] + equilib[4] - equilib[6] - equilib[7] - equilib[8] -
        equilib[10] - equilib[17] - Ith(ydot, 1);
    Ith(ydot, 3) = -equilib[0] + equilib[5] - equilib[8] - equilib[9] + equilib[11] +
        equilib[12] + equilib[13] + equilib[14];
    Ith(ydot, 4) = equilib[0] - equilib[1] - equilib[3] - 2. * equilib[5] - equilib[6] -
        equilib[11] - equilib[19];
    Ith(ydot, 5) = equilib[0] + equilib[1] - equilib[2] + 2. * equilib[3] + equilib[6] -
        equilib[7] + 2. * equilib[10] + equilib[11] - equilib[12] + 2. * equilib[15] +
        2. * equilib[16] + equilib[17] + equilib[19] - equilib[20] - equilib[21];
    Ith(ydot, 6) = equilib[8] + equilib[9] - equilib[10] - equilib[11] - equilib[12] -
        2. * equilib[13] - 2. * equilib[14] + equilib[18] + equilib[19] + equilib[20] +
        equilib[21];
    Ith(ydot, 7) = equilib[2] - equilib[3] + equilib[7] + equilib[12] + equilib[17] +
        equilib[20] + equilib[21];
    Ith(ydot, 8) = equilib[13] + equilib[14] - equilib[15] - equilib[16] - equilib[17] -
        equilib[18] - equilib[19] - equilib[20] - equilib[21];
    Ith(ydot, 9) = 0;


    double sigma_i;
    // Molar weight of mixture
    double W;
    // Cp
    double Cp;
    // Cv
    double Cv;
    // eta = 1 - M^2, M = w / af2 -- Mach number
    // Mass fraction
    double Yi, Yiprime;

    sum1 = 0.0;
    sum2 = 0.0;
    //cout << endl;
    for (int i = 0; i < num_gas_species; i++)
    {
        Yi = Ith(y, i + 1) * phyc.mol_weight[i] / Ith(y, num_gas_species + 2);
        //cout << "Yi = " << Yi << endl;
        sum1 += Ith(y, i + 1);
        sum2 += Yi * get_Cpi(i, Tcurr);
    }
    //cout << endl;
    W = Ith(y, num_gas_species + 2) / sum1;
    //cout << "W = " << W << endl;
    Cp = sum2;
    Cv = Cp - phyc.kR / W;
    gamma = Cp / Cv;
    double af2;
    //cout << "gamma = " << gamma << endl;
    af2 = gamma * Ith(y, num_gas_species + 1) / Ith(y, num_gas_species + 2);
    //cout << "af = " << pow(af2, 0.5) << endl;
    M2 = pow(Ith(y,num_gas_species + 3), 2) / af2;
    eta = 1.0 - M2;

    sum1 = 0.0;
    for (int i = 0; i < num_gas_species; i++)
    {
        sigma_i = (W / phyc.mol_weight[i]) - (get_Hi(i, Tcurr) / Cp /Tcurr);
        Yiprime = phyc.mol_weight[i] * Ith(ydot, i + 1) / Ith(y, num_gas_species + 2);
        sum1 += sigma_i * Yiprime;
    }
    sigma = sum1;
    for (int i = 0; i < num_gas_species; i++) Ith(ydot, i+1) = Ith(ydot, i+1) - Ith(y, i+1) * sigma / eta;
    //cout << "sigma = " << sigma << endl;
    // dP/dt
    Ith(ydot, 10) = -1.0 * Ith(y, 11) * pow(Ith(y, 12), 2) * sigma / eta;
    // drho/dt
    Ith(ydot, 11) = -1.0 * Ith(y, 11) * sigma / eta;
    // dw/dt
    Ith(ydot, 12) = Ith(y, 12) * sigma / eta;
    // dx/dt
    Ith(ydot, 13) = Ith(y, 12);

    delete[] forward;
    delete[] reverse;
    delete[] equilib;
    //delete[] Ymassfr;
    return 0;
}


void set_initial_ZND(N_Vector y, double rho, double P, double v, double x, double* Y0)
{
    double* ytmp = new double[num_gas_species];
    double konc;
    double checkrho = 0;
    double Msmes = 0.0;
    //cout << "konc = " << konc << endl;
    //Mol fractions
    Ith(y, 1) = Y0[0]; //H2
    Ith(y, 2) = Y0[1];
    Ith(y, 3) = Y0[2];  //O2
    Ith(y, 4) = Y0[3];
    Ith(y, 5) = Y0[4]; //OH
    Ith(y, 6) = Y0[5];
    Ith(y, 7) = Y0[6];
    Ith(y, 8) = Y0[7];
    Ith(y, 9) = Y0[8];    //N2
    Ith(y, 10) = P;
    Ith(y, 11) = rho;
    Ith(y, 12) = v;
    Ith(y, 13) = x;

    for (int i = 0; i < num_gas_species; i++)
    {
        Msmes += Ith(y, i + 1) * phyc.mol_weight[i];
    }
    //Mass fractions
    for (int i = 0; i < num_gas_species; i++)
    {
        Ith(y, i + 1) = Ith(y, i + 1) * phyc.mol_weight[i] / Msmes;
        ytmp[i] = Ith(y, i + 1);
    }

    //MOL. CONCENTR.
    for (int i = 0; i < num_gas_species; i++)
    {
        Ith(y, i + 1) = Ith(y, i + 1) * (rho) / phyc.mol_weight[i];
        if (fabs(Ith(y, i + 1)) < ATOL1 && Ith(y, i + 1) < 0)
        {
            Ith(y, i + 1) = 0;
        }
    }
    delete[] ytmp;
    return;
}


void integrate_ZND(double P, double rho, double v, double x, double* Y0)
{
    SUNContext sunctx;
    double t, tout;
    N_Vector y;
    N_Vector abstol;
    SUNMatrix A;
    SUNLinearSolver LS;
    void* cvode_mem;
    int retval, iout;
    t = 0.0;
    y = NULL;
    abstol = NULL;
    A = NULL;
    LS = NULL;
    cvode_mem = NULL;
    retval = SUNContext_Create(NULL, &sunctx);
    if (check_retval(&retval, "SUNContext_Create", 1)) return ;

    /* Initial conditions */
    y = N_VNew_Serial(num_gas_species + 4, sunctx);
    if (check_retval((void*)y, "N_VNew_Serial", 0)) return ;

    /* Initialize y */
    set_initial_ZND(y, rho,P, v, x, Y0);
    /* Set the vector absolute tolerance */
    abstol = N_VNew_Serial(num_gas_species + 4, sunctx);

    if (check_retval((void*)abstol, "N_VNew_Serial", 0)) return ;
    Ith(abstol, 1) = ATOL1;
    Ith(abstol, 2) = ATOL1;
    Ith(abstol, 3) = ATOL1;
    Ith(abstol, 4) = ATOL1;
    Ith(abstol, 5) = ATOL1;
    Ith(abstol, 6) = ATOL1;
    Ith(abstol, 7) = ATOL1;
    Ith(abstol, 8) = ATOL1;
    Ith(abstol, 9) = ATOL1;
    Ith(abstol, 10) = ATOL2;
    Ith(abstol, 11) = ATOL2;
    Ith(abstol, 12) = ATOL2;
    Ith(abstol, 13) = ATOL2;
    /* Call CVodeCreate to create the solver memory and specify the
     * Backward Differentiation Formula */
    cvode_mem = CVodeCreate(CV_BDF, sunctx);
    if (check_retval((void*)cvode_mem, "CVodeCreate", 0)) return;

    /* Call CVodeInit to initialize the integrator memory and specify the
     * user's right hand side function in y'=f(t,y), the initial time T0, and
     * the initial dependent variable vector y. */
    retval = CVodeInit(cvode_mem, init_right_part_ZND, T0, y);
    if (check_retval(&retval, "CVodeInit", 1)) return;

    /* Call CVodeSVtolerances to specify the scalar relative tolerance
     * and vector absolute tolerances */
    retval = CVodeSVtolerances(cvode_mem, RTOL, abstol);
    if (check_retval(&retval, "CVodeSVtolerances", 1)) return;

    /* Call CVodeRootInit to specify the root function g with 2 components */

    /* Create dense SUNMatrix for use in linear solves */
    A = SUNDenseMatrix(num_gas_species + 4, num_gas_species + 4, sunctx);
    if (check_retval((void*)A, "SUNDenseMatrix", 0)) return;

    /* Create dense SUNLinearSolver object for use by CVode */
    LS = SUNLinSol_Dense(y, A, sunctx);
    if (check_retval((void*)LS, "SUNLinSol_Dense", 0)) return;

    ofstream cons;
    auto str = std::to_string(koeff_topl);
    cons.open(str + "cons.dat");
    cons << "TITLE=\"" << "Graphics" << "\"" << endl;
    cons << R"(VARIABLES= "x", "H2", "H", "O2", "O", "OH", "HO2", "H2O", "H2O2", "N")" << endl;
    /* Attach the matrix and linear solver */
    retval = CVodeSetLinearSolver(cvode_mem, LS, A);
    if (check_retval(&retval, "CVodeSetLinearSolver", 1)) return;
    fout.open(str + "ZND.dat");
    fout << "TITLE=\"" << "Graphics" << "\"" << endl;
    fout << R"(VARIABLES= "x", "P", "rho", "v", "M", "sigma", "Tcurr", "gamma", "sound speed")" << endl;
    double dt = pow(10, -6);
    tout = pow(10, -6);
    iout = 0;
    double* Ymassfrac = new double[num_gas_species];
    double Tcheck;

   /*cout << "before integ" << endl;
   cout << Ith(y, 1) << endl;
   cout << Ith(y, 2) << endl;
   cout << Ith(y, 3) << endl;
   cout << Ith(y, 4) << endl;
   cout << Ith(y, 5) << endl;
   cout << Ith(y, 6) << endl;
   cout << Ith(y, 7) << endl;
   cout << Ith(y, 8) << endl;
   cout << Ith(y, 9) << endl;
   cout << Ith(y, 10) << endl; 
   cout << Ith(y, 11) << endl;
   cout << Ith(y, 12) << endl;
   cout << Ith(y, 13) << endl;*/
   double Yi;
   double max_sigma_x;
   cout << "M = " << pow(M2, 0.5) << endl;
    while (pow(M2, 0.5) < 0.998) {
        retval = CVode(cvode_mem, tout, y, &t, CV_NORMAL);
        if (check_retval(&retval, "CVode", 1)) break;
        if (retval == CV_SUCCESS && pow(M2, 0.5) < 0.998) {
            //cout << endl << "YmassFr" << endl;
            for (int i = 0; i < num_gas_species; i++)
            {
                Ymassfrac[i] = Ith(y, i + 1) * phyc.mol_weight[i] / Ith(y, num_gas_species + 2);
                //cout << "Yi = " << Yi << endl;
            }

            cons << Ith(y, 13) << " " << Ymassfrac[0] << " " << Ymassfrac[1] << " " << Ymassfrac[2] << " " << Ymassfrac[3] << " " << Ymassfrac[4] << " " << Ymassfrac[5] << " " << Ymassfrac[6] << " " << Ymassfrac[7] << " " << Ymassfrac[8] << endl;
            cout << "M = " << pow(M2, 0.5)  << endl;
            fout << Ith(y, num_gas_species + 4) << " " << Ith(y, num_gas_species + 1) << " " << Ith(y, num_gas_species + 2) << " " 
                << Ith(y, num_gas_species + 3) <<  " " << pow(M2, 0.5) << " " << sigma << " " << Tcurr << " " << gamma << " " 
                << pow(gamma * Ith(y, num_gas_species + 1) / Ith(y, num_gas_species + 2), 0.5) << endl;

            tout += dt;
            iout++;
        }
    }
    cout << "M = " << pow(M2, 0.5) << endl;
    //printf("\nFinal Statistics:\n");
    //cout << "t = " << OHt << endl;
    //retval = CVodePrintAllStats(cvode_mem, stdout, SUN_OUTPUTFORMAT_TABLE);

    /* check the solution error */
    //retval = check_ans(y, t, RTOL, abstol);

    /* Free memory */


    N_VDestroy(y);                            /* Free y vector */
    N_VDestroy(abstol);                       /* Free abstol vector */
    CVodeFree(&cvode_mem);                    /* Free CVODE memory */
    SUNLinSolFree(LS);                        /* Free the linear solver memory */
    SUNMatDestroy(A);                         /* Free the matrix memory */
    SUNContext_Free(&sunctx);                 /* Free the SUNDIALS context */
    fout.close();
    cons.close();

}
int main()
{
    koeff_topl = 2.0;
    double X;
    double uCj;
    double Msmes = 0.;
    double Tin = 295;
    double Y0[9] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
    Y0[0] = koeff_topl * 2.0 / (koeff_topl * 2. + 1. + 3.76);
    Y0[2] = 1.0 / (koeff_topl * 2. + 1. + 3.76);
    Y0[8] = 3.76 / (koeff_topl * 2. + 1. + 3.76);
    cout << "koef topl = " << koeff_topl << " " << Y0[0] << " " << Y0[2] << " " << Y0[8] << endl;
    double Y0mf[9];
    init_consts(num_gas_species, num_react);
    for (int i = 0; i < num_gas_species; i++) Msmes += Y0[i] * phyc.mol_weight[i];
    //Mass fractions
    for (int i = 0; i < num_gas_species; i++) Y0mf[i] = Y0[i] * phyc.mol_weight[i] / Msmes;
    double p1 = 1. * 0.101325;


    Find_velocity(Y0, p1, Tin, uCj, X);
    //uCj = 1.97422;
    cout << "uCj =  " << uCj << endl;

    double rho1 = p1 / get_gas_constant(num_gas_species, Y0mf)/Tin;
    double uVn2 = 0.5;
    double uVn1 = 0.0;
    while (abs(uVn2 - uVn1) > pow(10, -8))
    {
        uVn1 = uVn2;
        uVn2 = uVn1 - Func(uVn1, uCj, p1, rho1, Tin, Y0mf)/Func_prime(uVn1, uCj, p1, rho1, Tin, Y0mf);
    }
    cout << "uVnL =  " << uCj - uVn2 << endl;
    double rhoVn = rho1 * (uCj / uVn2);
    double pVn = ((rho1 * (uCj * uCj)) - (rhoVn * (uVn2 * uVn2))) + p1;

    cout << "integrate" << endl;
    integrate_ZND(pVn, rhoVn, uVn2, 0, Y0);
    cout << "uCj =  " << uCj << endl;
    cout << "uVnL =  " << uCj - uVn2 << endl;
    cout << "rhoVn =  " << rhoVn << endl;
    cout << "PVn =  " << pVn << endl;
    return 0;
}

// Запуск программы: CTRL+F5 или меню "Отладка" > "Запуск без отладки"
// Отладка программы: F5 или меню "Отладка" > "Запустить отладку"

// Советы по началу работы 
//   1. В окне обозревателя решений можно добавлять файлы и управлять ими.
//   2. В окне Team Explorer можно подключиться к системе управления версиями.
//   3. В окне "Выходные данные" можно просматривать выходные данные сборки и другие сообщения.
//   4. В окне "Список ошибок" можно просматривать ошибки.
//   5. Последовательно выберите пункты меню "Проект" > "Добавить новый элемент", чтобы создать файлы кода, или "Проект" > "Добавить существующий элемент", чтобы добавить в проект существующие файлы кода.
//   6. Чтобы снова открыть этот проект позже, выберите пункты меню "Файл" > "Открыть" > "Проект" и выберите SLN-файл.
