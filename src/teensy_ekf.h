/*
 * TeensyEKF: Extended Kalman Filter for embedded processors.
 *
 * Copyright (C) 2021 Morten Veng
 *
 * MIT License
 */

/**
  * Initializes an EKF structure.
  * @param ekf pointer to EKF structure to initialize
  * @param n number of state variables
  * @param m number of observables
  */
void ekf_init(void *ekf, int n, int m);

/**
  * Runs one step of EKF prediction and update. Your code should first build a model, setting
  * the contents of <tt>ekf.fx</tt>, <tt>ekf.F</tt>, <tt>ekf.hx</tt>, and <tt>ekf.H</tt> to appropriate values.
  * @param ekf pointer to structure EKF 
  * @param z array of measurement (observation) values
  * @param z_update array containing measurement update logic
  * @return true on success, false on failure caused by non-positive-definite matrix.
  */
bool ekf_step(void *ekf, double *z, bool *z_update);
