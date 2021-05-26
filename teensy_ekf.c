/*
 * TeensyEKF: Extended Kalman Filter for embedded processors
 *
 * Copyright (C) 2021 Morten Veng
 *
 * MIT License
 */

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>

/* Fills matrix with zeros */
static void zeros(double *a, int m, int n)
{
    int i;
    for (i = 0; i < m * n; ++i)
        a[i] = 0;
}

/* Fetches the k'th row from a matrix with width n */
static void row(double *a, double *b, int k, int n)
{
    int i;
    for (i = 0; i < n; ++i)  // Looping cols
        b[i] = a[n * k + i]; // Idx = ncols * row + col
}

/* Matrix transpose */
static void transpose(double *a, double *at, int m, int n)
{
    int i, j;

    for (i = 0; i < m; ++i)
        for (j = 0; j < n; ++j)
            at[j * m + i] = a[i * n + j];
}

/* C <- A * B */
static void mulmat(double *a, double *b, double *c, int arows, int acols, int bcols)
{
    int i, j, l;
    for (i = 0; i < arows; ++i)     // Looping rows
        for (j = 0; j < bcols; ++j) // Looping cols
        {
            c[i * bcols + j] = 0;
            for (l = 0; l < acols; ++l)
                c[i * bcols + j] += a[i * acols + l] * b[l * bcols + j];
        }
}

/* C <- A_mxn * B_nx1 */
static void mulvec(double *a, double *b, double *c, int m, int n)
{
    int i, j;
    for (i = 0; i < m; ++i) // Looping rows
    {
        c[i] = 0;
        for (j = 0; j < n; ++j)          // Looping cols
            c[i] += a[n * i + j] * b[j]; // Idx = ncols * row + col
    }
}

/* C <- B_1xn * A_mxn^T  */
static void mulvecflipped(double *a, double *b, double *c, int m, int n)
{
    int i, j;
    for (i = 0; i < m; ++i) // Looping rows
    {
        c[i] = 0;
        for (j = 0; j < n; ++j)          // Looping cols
            c[i] += a[m * j + i] * b[j]; // Idx = ncols * row + col
    }
}

/* sum <- A_1xn * B_nx1 */
static double dotprod(double *a, double *b, int n)
{
    int i;
    double sum = 0;
    for (i = 0; i < n; ++i)
        sum += a[i] * b[i];

    return sum;
}

/* C_nxn <- A_mx1 * B_1xn */
static void matprod(double *a, double *b, double *c, int m, int n)
{
    int i, j;
    for (i = 0; i < m; ++i)             // Looping rows
        for (j = 0; j < n; ++j)         // Looping cols
            c[n * i + j] = a[i] * b[j]; // Idx = ncols * row + col
}

/* B <- A_nx1 * factor | Note: works for matrices and vectors */
static void mulfac(double *a, double factor, double *b, int n)
{
    int i;
    for (i = 0; i < n; ++i)
        b[i] = a[i] * factor;
}

/* B <- A */
static void copy(double *a, double *b, int n)
{
    int i;
    for (i = 0; i < n; ++i)
        b[i] = a[i];
}

/* A <- A + B */
static void accum(double *a, double *b, int m, int n)
{
    int i;

    for (i = 0; i < m * n; ++i)
        a[i] += b[i];
}

/* A <- -A */
static void negate(double *a, int m, int n)
{
    int i;
    for (i = 0; i < m * n; ++i)
        a[i] = -a[i];
}

/* ------------------------------------------------
 * TeensyEKF Core 
 * ------------------------------------------------*/

#include "teensy_ekf.h"

typedef struct
{
    double *x; // State vector
    double *P; // Prediction error covariance
    double *Q; // Process noise covariance
    double *R; // Measurement error covariance

    double *fx; // Output of user defined f() state-transition function
    double *hx; // Output of user defined h() measurement function
    double *F;  // Jacobian of process model
    double *H;  // Jacobian of measurement model

    double *Hrow; // Single row of measurement Jacobian
    double *Krow; // Single row of Kalman gain
    double *P_xy; // Part of Kalman gain

    // Temporary storage
    double *tmp0;
    double *tmp1;

} ekf_t;

static void unpack(void *v, ekf_t *ekf, int n, int m)
{
    // Skip 'n' and 'm' in data structure
    char *cptr = (char *)v;
    cptr += 2 * sizeof(int);

    // Load double structures
    double *dptr = (double *)cptr;
    ekf->x = dptr;
    dptr += n;
    ekf->P = dptr;
    dptr += n * n;
    ekf->Q = dptr;
    dptr += n * n;
    ekf->R = dptr;
    dptr += m * m;

    ekf->fx = dptr;
    dptr += n;
    ekf->hx = dptr;
    dptr += m;
    ekf->F = dptr;
    dptr += n * n;
    ekf->H = dptr;
    dptr += m * n;

    ekf->Hrow = dptr;
    dptr += n;
    ekf->Krow = dptr;
    dptr += n;
    ekf->P_xy = dptr;
    dptr += n;

    ekf->tmp0 = dptr;
    dptr += n * n;
    ekf->tmp1 = dptr;
}

void ekf_init(void *v, int n, int m)
{
    // Retrieve 'n' and 'm'
    int *ptr = (int *)v;
    *ptr = n;
    ptr++;
    *ptr = m;

    // Unpack rest of incoming structure for initialization
    ekf_t ekf;
    unpack(v, &ekf, n, m);

    // Zero-out structures as default
    zeros(ekf.x, n, 1);
    zeros(ekf.P, n, n);
    zeros(ekf.Q, n, n);
    zeros(ekf.R, m, m);
    zeros(ekf.F, n, n);
    zeros(ekf.H, m, n);
}

bool ekf_step(void *v, double *z, bool *z_update)
{
    // Unpack incoming structure
    int *ptr = (int *)v;
    int n = *ptr;
    ptr++;
    int m = *ptr;

    ekf_t ekf;
    unpack(v, &ekf, n, m);

    // P = F * P * F^T + Q                                          | Covariance prediction
    mulmat(ekf.F, ekf.P, ekf.tmp0, n, n, n);
    transpose(ekf.F, ekf.tmp1, n, n);
    mulmat(ekf.tmp0, ekf.tmp1, ekf.P, n, n, n);
    accum(ekf.P, ekf.Q, n, n);

    double P_yy, y; // Scalars in iterative EKF correction
    for (int i = 0; i < m; ++i)
    {
        if (z_update[i] == false)
            continue;

        // Copy i'th row of H
        row(ekf.H, ekf.Hrow, i, n);
        /* TO DO:
         * We can rewrite the MeasurementModel with an additional input; the row number to be evaluated.
         * Currently, this implementation just computes H at once - this is equivalent of a traditional EKF.
         */

        // P_xy = P * H_i'                                          | Part of Kalman gain (the nominator)
        mulvec(ekf.P, ekf.Hrow, ekf.P_xy, n, n);

        // P_yy = (H_i * P * H_i' + R[i][i]);                       | Part of Kalman gain (the denominator)
        P_yy = dotprod(ekf.Hrow, ekf.P_xy, n) + ekf.R[i * m + i];

        if (P_yy > 0.0) // Ensure a nonzero denominator
        {
            // K_i = P * H_i' (H_i * P * H_i' + R[i][i])^-1         | Kalman Gain
            mulfac(ekf.P_xy, 1.0 / P_yy, ekf.Krow, n);

            // y_i = z_i - zhat_i                                   | Innovation
            y = z[i] - ekf.hx[i];

            // xhat_i = xhat_i + K_i * y_i                          | Updated state estimate
            mulfac(ekf.Krow, y, ekf.tmp1, n);
            accum(ekf.fx, ekf.tmp1, n, 1);

            // P = P - K_i * H_i * P                                | Updated covariance estimate
            // Note that 'tmp1' is used as a vector here
            mulvecflipped(ekf.P, ekf.Hrow, ekf.tmp1, n, n);
            matprod(ekf.Krow, ekf.tmp1, ekf.tmp0, n, n);
            negate(ekf.tmp0, n, n);
            accum(ekf.P, ekf.tmp0, n, n);

            /* Notes:
             * The EKF can be subject to numerical round-off errors in the covariance update step. 
             * If your filter diverges consider: 
             * 1) Scaling your states to increase their respective covariances Q and R
             * 2) Implement Joseph-form update of the covariance estimate
             */
        }
    }

    // Write corrected estimate to the state vector
    copy(ekf.fx, ekf.x, n);

    return true; // success
}