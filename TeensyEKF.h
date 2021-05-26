/*
 * TeensyEKF: Extended Kalman Filter intended for Arduino Teensy devices.
 *
 * Copyright (C) 2021 Morten Veng
 *
 * MIT License
 */

#include <stdio.h>
#include <stdlib.h>

typedef struct
{
    int n; // Number of state values
    int m; // Number of observables

    double x[Nsta];       // State vector
    double P[Nsta][Nsta]; // Prediction error covariance
    double Q[Nsta][Nsta]; // Process noise covariance
    double R[Mobs][Mobs]; // Measurement error covariance

    double fx[Nsta];      // Output of user defined f() state-transition function
    double hx[Mobs];      // Output of user defined h() measurement function
    double F[Nsta][Nsta]; // Jacobian of process model
    double H[Mobs][Nsta]; // Jacobian of measurement model

    double Hrow[Nsta]; // Single row of measurement Jacobian
    double Krow[Nsta]; // Single row of Kalman gain
    double P_xy[Nsta]; // Part of Kalman gain

    // Temporary storage
    double tmp0[Nsta][Nsta];
    double tmp1[Nsta][Nsta];

} ekf_t;

// Support both Arduino and command-line versions
#ifndef MAIN
extern "C"
{
#endif
    void ekf_init(void *, int, int);
    bool ekf_step(void *, double *, bool *);
#ifndef MAIN
}
#endif

/*
 * A header-only class for the Extended Kalman Filter.  Your implementing class should #define the constants Nsta and Mobs
 * and then #include <TeensyEKF.h>  You will also need to implement a model() method for your application.
 */
class TeensyEKF
{

private:
    ekf_t ekf;

protected:
    /**
     * The current state.
     */
    double *x;

    /**
     * Initializes a TeensyEKF object.
     */
    TeensyEKF()
    {
        ekf_init(&this->ekf, Nsta, Mobs);
        this->x = this->ekf.x;
    }

    /**
     * Deallocates memory for a TeensyEKF object.
     */
    ~TeensyEKF() {}

    /**
     * Implement these functions for your EKF model.
     * @param fx gets output of state-transition function <i>f(x<sub>0 .. n-1</sub>)</i>
     * @param F gets <i>n &times; n</i> Jacobian of <i>f(x)</i>
     * @param hx gets output of observation function <i>h(x<sub>0 .. n-1</sub>)</i>
     * @param H gets <i>m &times; n</i> Jacobian of <i>h(x)</i>
     */
    virtual void StateTransitionModel(double fx[Nsta], double F[Nsta][Nsta]) = 0;
    virtual void MeasurementModel(double hx[Mobs], double H[Mobs][Nsta]) = 0;

    /**
     * Sets the specified value of the prediction error covariance. <i>P<sub>i,j</sub> = value</i>
     * @param i row index
     * @param j column index
     * @param value value to set
     */
    void setP(int i, int j, double value)
    {
        this->ekf.P[i][j] = value;
    }

    /**
     * Sets the specified value of the process noise covariance. <i>Q<sub>i,j</sub> = value</i>
     * @param i row index
     * @param j column index
     * @param value value to set
     */
    void setQ(int i, int j, double value)
    {
        this->ekf.Q[i][j] = value;
    }

    /**
     * Sets the specified value of the observation noise covariance. <i>R<sub>i,j</sub> = value</i>
     * @param i row index
     * @param j column index
     * @param value value to set
     */
    void setR(int i, int j, double value)
    {
        this->ekf.R[i][j] = value;
    }

    /**
     * Utility function to set Q triplets
     * @param start_index first triplet index 
     * @param values triplet values
     * @param scale scaling parameter for convenience in tuning (scales all values equally)
     */
    void setQ3(int start_index, double values[3], double scale = 1.0)
    {
        setQ(start_index + 0, start_index + 0, scale * values[0]);
        setQ(start_index + 1, start_index + 1, scale * values[1]);
        setQ(start_index + 2, start_index + 2, scale * values[2]);
    };

    /**
     * Utility function to set R triplets
     * @param start_index first triplet index 
     * @param values triplet values
     * @param scale scaling parameter for convenience in tuning (scales all values equally)
     */
    void setR3(int start_index, double values[3], double scale = 1.0)
    {
        setR(start_index + 0, start_index + 0, scale * values[0]);
        setR(start_index + 1, start_index + 1, scale * values[1]);
        setR(start_index + 2, start_index + 2, scale * values[2]);
    };

    // Prints an 'm x n' matrix to serial. The input is the address to the first element i.e. &this->ekf.H[0][0]
    void mdump(double *a, int m, int n)
    {
        Serial.println("(");
        for (int i = 0; i < m; ++i)
        {
            for (int j = 0; j < n; ++j)
            {
                Serial.print(a[i * n + j], 6);
                Serial.print(" ");
            }
            Serial.println();
        }
        Serial.println(")");
    };

    // Prints an 'n' element vector to serial (forced row format). The input is the address to the first element i.e. &this->ekf.Hrow[0]
    void vdump(double *a, int n)
    {
        Serial.print("( ");
        for (int i = 0; i < n; ++i)
        {
            Serial.print(a[i], 6);
            Serial.print(" ");
        }
        Serial.println(")");
    }

    // Prints update vector
    void print_update(bool *a, int n)
    {
        Serial.print("( ");
        for (int i = 0; i < n; ++i)
        {
            Serial.print(a[i]);
            Serial.print(" ");
        }
        Serial.println(")");
    }

public:
    /**
     * Returns a reference to the state vector at a given index.
     * @param i the index (at least 0 and less than <i>n</i>
     * @return reference to state vector at starting index
     */
    double *getXref(int i)
    {
        return &this->ekf.x[i];
    }

    /**
     * Returns the state element at a given index.
     * @param i the index (at least 0 and less than <i>n</i>
     * @return state value at index
     */
    double getX(int i)
    {
        return this->ekf.x[i];
    }

    /**
     * Sets the state element at a given index.
     * @param i the index (at least 0 and less than <i>n</i>
     * @param value value to set
     */
    void setX(int i, double value)
    {
        this->ekf.x[i] = value;
    }

    /**
     * Performs one step of the prediction and update.
     * @param z observation vector, length <i>m</i>
     * @return true on success, false on failure caused by non-positive-definite matrix.
     */
    bool step(double *z, bool *z_update)
    {
        this->StateTransitionModel(this->ekf.fx, this->ekf.F);
        this->MeasurementModel(this->ekf.hx, this->ekf.H);
        bool success = ekf_step(&this->ekf, z, z_update);

        /* Example usage of dump functions: */
        //mdump(&this->ekf.P[0][0],Nsta,Nsta);
        //vdump(&this->ekf.Hrow[0], Nsta);
        //print_update(z_update, Mobs);

        return success;
    }

    /**
     * Utility function that checks an update vector for available updates
     * @param z_update logical vector that is true if new sensor data has arrived, length <i>m</i>
     * @return true if any new information is available
     */
    bool dataAvailable(bool *z_update)
    {
        int i;
        for (i = 0; i < Mobs; i++)
            if (z_update[i])
                return true;

        return false;
    }
};
