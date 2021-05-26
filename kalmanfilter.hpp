/*
   Sensor Fusion Code
*/

#define Mobs (12) // Measurement dimensions (used within TeensyEKF)
#define Nsta (21) // State dimensions (used within TeensyEKF)
#define EKF_RATE (400) // Target rate the EKF is to be executed at [Hz] (reduced effects of sudden sample rate fluctuations)

#include "src/TeensyEKF.h"



class KalmanFilter : public TeensyEKF {
  private:
    // Temporary variables used in the EKF model
    double trig_temp0, trig_temp1, trig_temp2, trig_temp3, trig_temp4, trig_temp5 = 0.0;

    // Sensor variances
    double R_pose[3] = {1e-6, 1e-6, 1e-6};
    double R_accel[3] = {0.05, 0.05, 0.05};
    double R_eta[3] = {7.6e-7, 7.6e-7, 7.6e-7};
    double R_omega[3] = {3.05e-6, 3.05e-6, 3.05e-6};

    // Time variables to run the EKF at a (constant) real-time rate
    const unsigned ekf_duration_micros = 1000000 / EKF_RATE;  // Sample time in microseconds
    elapsedMicros ekf_timer_micros;                           // Time since last EKF call in microseconds
    double Ts = 1.0 / EKF_RATE;                               // Discrete sample time value

    // Timer variable
    unsigned elapsed;

    // Function that sets the Q values based on an integration model
    void setQ_empirical(double translation_variance, double rotation_variance, double v_thrust, double v_tau)  {
      setQ(0, 0, ((Ts * Ts * Ts * Ts)*translation_variance) / 4.0);
      setQ(0, 3, ((Ts * Ts * Ts)*translation_variance) / 2.0);
      setQ(0, 6, ((Ts * Ts)*translation_variance) / 2.0);
      setQ(1, 1, ((Ts * Ts * Ts * Ts)*translation_variance) / 4.0);
      setQ(1, 4, ((Ts * Ts * Ts)*translation_variance) / 2.0);
      setQ(1, 7, ((Ts * Ts)*translation_variance) / 2.0);
      setQ(2, 2, ((Ts * Ts * Ts * Ts)*translation_variance) / 4.0);
      setQ(2, 5, ((Ts * Ts * Ts)*translation_variance) / 2.0);
      setQ(2, 8, ((Ts * Ts)*translation_variance) / 2.0);
      setQ(3, 0, ((Ts * Ts * Ts)*translation_variance) / 2.0);
      setQ(3, 3, (Ts * Ts)*translation_variance);
      setQ(3, 6, Ts * translation_variance);
      setQ(4, 1, ((Ts * Ts * Ts)*translation_variance) / 2.0);
      setQ(4, 4, (Ts * Ts)*translation_variance);
      setQ(4, 7, Ts * translation_variance);
      setQ(5, 2, ((Ts * Ts * Ts)*translation_variance) / 2.0);
      setQ(5, 5, (Ts * Ts)*translation_variance);
      setQ(5, 8, Ts * translation_variance);
      setQ(6, 0, ((Ts * Ts)*translation_variance) / 2.0);
      setQ(6, 3, Ts * translation_variance);
      setQ(6, 6, translation_variance);
      setQ(7, 1, ((Ts * Ts)*translation_variance) / 2.0);
      setQ(7, 4, Ts * translation_variance);
      setQ(7, 7, translation_variance);
      setQ(8, 2, ((Ts * Ts)*translation_variance) / 2.0);
      setQ(8, 5, Ts * translation_variance);
      setQ(8, 8, translation_variance);
      setQ(9, 9, ((Ts * Ts * Ts * Ts)*rotation_variance) / 4.0);
      setQ(9, 12, ((Ts * Ts * Ts)*rotation_variance) / 2.0);
      setQ(10, 10, ((Ts * Ts * Ts * Ts)*rotation_variance) / 4.0);
      setQ(10, 13, ((Ts * Ts * Ts)*rotation_variance) / 2.0);
      setQ(11, 11, ((Ts * Ts * Ts * Ts)*rotation_variance) / 4.0);
      setQ(11, 14, ((Ts * Ts * Ts)*rotation_variance) / 2.0);
      setQ(12, 9, ((Ts * Ts * Ts)*rotation_variance) / 2.0);
      setQ(12, 12, (Ts * Ts)*rotation_variance);
      setQ(13, 10, ((Ts * Ts * Ts)*rotation_variance) / 2.0);
      setQ(13, 13, (Ts * Ts)*rotation_variance);
      setQ(14, 11, ((Ts * Ts * Ts)*rotation_variance) / 2.0);
      setQ(14, 14, (Ts * Ts)*rotation_variance);
      setQ(15, 15, v_thrust);
      setQ(16, 16, v_thrust);
      setQ(17, 17, v_thrust);
      setQ(18, 18, v_tau);
      setQ(19, 19, v_tau);
      setQ(20, 20, v_tau);
    };

  public:
    bool Z_update[Mobs];
    double Z[Mobs];

    // Constructor - Sets initial values for x, Q, R and P
    KalmanFilter() {

      // Set sensor covariance
      setR3(0, R_pose);
      setR3(3, R_accel);
      setR3(6, R_eta);
      setR3(9, R_omega);

      // Set process noise based on time sample (translational variance, rotational variance, fext variance, text variance)
      setQ_empirical(1.0e-2, 1.0, 1.0, 1.0);

      // Set estimate covariance diagonal to 'large' value (the remaining cells are initialized at zero within TeensyEKF)
      for (byte i = 0; i < 15; i++)
        setP(i, i, 1.0);
      for (byte i = 15; i < Nsta; i++)
        setP(i, i, 1e-4);

    };

    unsigned update() {

      if (!dataAvailable(Z_update) || ekf_timer_micros < ekf_duration_micros)
        return 0; // Debounce timer to constrain the EKF update rate to the specified interval

      // Compensate time step by using interval since last EKF update (The compensated time step is within [Ts, 2*Ts])
      Ts = constrain((unsigned)ekf_timer_micros, ekf_duration_micros, 2 * ekf_duration_micros);
      Ts /= 1000000.0; // microseconds -> seconds
      ekf_timer_micros = 0; // Reset debounce timer


      elapsed = micros();

      // Update the EKF with the current measurements
      step(Z, Z_update);

      return micros() - elapsed;
    };

  
    // Prints 'n' states of the state vector starting from the given start index
    void printvector(int idx_start, int n) {
        if (idx_start + n > Nsta)
          return;

        // Fetch states by reference
        double* x_vector = getXref(idx_start);
        for (int i = 0; i < n; i++) {
          Serial.print(x_vector[i]);
          Serial.print(",");
        }   
    };


  protected:

    void StateTransitionModel(double fx[Nsta], double F[Nsta][Nsta]) {
      trig_temp0 = sin(this->x[9]), trig_temp1 = sin(this->x[10]), trig_temp2 = sin(this->x[11]),
      trig_temp3 = cos(this->x[9]), trig_temp4 = cos(this->x[10]), trig_temp5 = cos(this->x[11]);

      /* Prediction */
      fx[0] = this->x[0] + Ts * this->x[3] + ((Ts * Ts) * this->x[6]) / 2.0;
      fx[1] = this->x[1] + Ts * this->x[4] + ((Ts * Ts) * this->x[7]) / 2.0;
      fx[2] = this->x[2] + Ts * this->x[5] + ((Ts * Ts) * this->x[8]) / 2.0;
      fx[3] = this->x[3] + Ts * this->x[6];
      fx[4] = this->x[4] + Ts * this->x[7];
      fx[5] = this->x[5] + Ts * this->x[8];
      fx[6] = this->x[15] / (MASS * 1.0E+3) + (Thrust * (trig_temp0 * trig_temp2 + trig_temp3 * trig_temp5 * trig_temp1)) / MASS;
      fx[7] = this->x[16] / (MASS * 1.0E+3) - (Thrust * (trig_temp5 * trig_temp0 - trig_temp3 * trig_temp2 * trig_temp1)) / MASS;
      fx[8] = -GRAVITY + this->x[17] / (MASS * 1.0E+3) + (Thrust * trig_temp3 * trig_temp4) / MASS;
      fx[9] = this->x[9] + Ts * (this->x[12] + (this->x[14] * trig_temp3 * trig_temp1) / (trig_temp3 * trig_temp3 * trig_temp4 + trig_temp4 * trig_temp0 * trig_temp0) + (this->x[13] * trig_temp0 * trig_temp1) / (trig_temp3 * trig_temp3 * trig_temp4 + trig_temp4 * trig_temp0 * trig_temp0));
      fx[10] = this->x[10] + Ts * ((this->x[13] * trig_temp3) / (trig_temp3 * trig_temp3 + trig_temp0 * trig_temp0) - (this->x[14] * trig_temp0) / (trig_temp3 * trig_temp3 + trig_temp0 * trig_temp0));
      fx[11] = this->x[11] + Ts * ((this->x[14] * trig_temp3) / (trig_temp3 * trig_temp3 * trig_temp4 + trig_temp4 * trig_temp0 * trig_temp0) + (this->x[13] * trig_temp0) / (trig_temp3 * trig_temp3 * trig_temp4 + trig_temp4 * trig_temp0 * trig_temp0));
      fx[12] = this->x[12] + (Ts * (Tau[0] + this->x[18] / 1.0E+3 + I_YY * this->x[13] * this->x[14] - I_ZZ * this->x[13] * this->x[14])) / I_XX;
      fx[13] = this->x[13] + (Ts * (Tau[1] + this->x[19] / 1.0E+3 - I_XX * this->x[12] * this->x[14] + I_ZZ * this->x[12] * this->x[14])) / I_YY;
      fx[14] = this->x[14] + (Ts * (Tau[2] + this->x[20] / 1.0E+3 + I_XX * this->x[12] * this->x[13] - I_YY * this->x[12] * this->x[13])) / I_ZZ;
      fx[15] = this->x[15];
      fx[16] = this->x[16];
      fx[17] = this->x[17];
      fx[18] = this->x[18];
      fx[19] = this->x[19];
      fx[20] = this->x[20];

      /* State Transition Jacobian */
      F[0][0] = 1.0;
      F[0][3] = Ts;
      F[0][6] = (Ts * Ts) / 2.0;
      F[1][1] = 1.0;
      F[1][4] = Ts;
      F[1][7] = (Ts * Ts) / 2.0;
      F[2][2] = 1.0;
      F[2][5] = Ts;
      F[2][8] = (Ts * Ts) / 2.0;
      F[3][3] = 1.0;
      F[3][6] = Ts;
      F[4][4] = 1.0;
      F[4][7] = Ts;
      F[5][5] = 1.0;
      F[5][8] = Ts;
      F[6][9] = (Thrust * (trig_temp3 * trig_temp2 - trig_temp5 * trig_temp0 * trig_temp1)) / MASS;
      F[6][10] = (Thrust * trig_temp3 * trig_temp5 * trig_temp4) / MASS;
      F[6][11] = (Thrust * (trig_temp5 * trig_temp0 - trig_temp3 * trig_temp2 * trig_temp1)) / MASS;
      F[6][15] = 1.0 / (MASS * 1.0E+3);
      F[7][9] = -(Thrust * (trig_temp3 * trig_temp5 + trig_temp0 * trig_temp2 * trig_temp1)) / MASS;
      F[7][10] = (Thrust * trig_temp3 * trig_temp4 * trig_temp2) / MASS;
      F[7][11] = (Thrust * (trig_temp0 * trig_temp2 + trig_temp3 * trig_temp5 * trig_temp1)) / MASS;
      F[7][16] = 1.0 / (MASS * 1.0E+3);
      F[8][9] = -(Thrust * trig_temp4 * trig_temp0) / MASS;
      F[8][10] = -(Thrust * trig_temp3 * trig_temp1) / MASS;
      F[8][17] = 1.0 / (MASS * 1.0E+3);
      F[9][9] = (trig_temp4 + Ts * this->x[13] * trig_temp3 * trig_temp1 - Ts * this->x[14] * trig_temp0 * trig_temp1) / trig_temp4;
      F[9][10] = Ts * 1.0 / trig_temp4 * trig_temp4 * (this->x[14] * trig_temp3 + this->x[13] * trig_temp0);
      F[9][12] = Ts;
      F[9][13] = (Ts * trig_temp0 * trig_temp1) / trig_temp4;
      F[9][14] = (Ts * trig_temp3 * trig_temp1) / trig_temp4;
      F[10][9] = -Ts * (this->x[14] * trig_temp3 + this->x[13] * trig_temp0);
      F[10][10] = 1.0;
      F[10][13] = Ts * trig_temp3;
      F[10][14] = -Ts * trig_temp0;
      F[11][9] = (Ts * (this->x[13] * trig_temp3 - this->x[14] * trig_temp0)) / trig_temp4;
      F[11][10] = Ts * 1.0 / trig_temp4 * trig_temp4 * trig_temp1 * (this->x[14] * trig_temp3 + this->x[13] * trig_temp0);
      F[11][11] = 1.0;
      F[11][13] = (Ts * trig_temp0) / trig_temp4;
      F[11][14] = (Ts * trig_temp3) / trig_temp4;
      F[12][12] = 1.0;
      F[12][13] = (Ts * this->x[14] * (I_YY - I_ZZ)) / I_XX;
      F[12][14] = (Ts * this->x[13] * (I_YY - I_ZZ)) / I_XX;
      F[12][18] = Ts / (I_XX * 1.0E+3);
      F[13][12] = -(Ts * this->x[14] * (I_XX - I_ZZ)) / I_YY;
      F[13][13] = 1.0;
      F[13][14] = -(Ts * this->x[12] * (I_XX - I_ZZ)) / I_YY;
      F[13][19] = Ts / (I_YY * 1.0E+3);
      F[14][12] = (Ts * this->x[13] * (I_XX - I_YY)) / I_ZZ;
      F[14][13] = (Ts * this->x[12] * (I_XX - I_YY)) / I_ZZ;
      F[14][14] = 1.0;
      F[14][20] = Ts / (I_ZZ * 1.0E+3);
      F[15][15] = 1.0;
      F[16][16] = 1.0;
      F[17][17] = 1.0;
      F[18][18] = 1.0;
      F[19][19] = 1.0;
      F[20][20] = 1.0;

    };

    void MeasurementModel(double hx[Mobs], double H[Mobs][Nsta]) {
      trig_temp0 = sin(this->x[9]), trig_temp1 = sin(this->x[10]), trig_temp2 = sin(this->x[11]),
      trig_temp3 = cos(this->x[9]), trig_temp4 = cos(this->x[10]), trig_temp5 = cos(this->x[11]);

      /* Measurement Function */
      hx[0] = this->x[0];
      hx[1] = this->x[1];
      hx[2] = this->x[2];
      hx[3] = -this->x[8] * trig_temp1 + this->x[6] * trig_temp5 * trig_temp4 + this->x[7] * trig_temp4 * trig_temp2;
      hx[4] = -this->x[6] * (trig_temp3 * trig_temp2 - trig_temp5 * trig_temp0 * trig_temp1) + this->x[7] * (trig_temp3 * trig_temp5 + trig_temp0 * trig_temp2 * trig_temp1) + this->x[8] * trig_temp4 * trig_temp0;
      hx[5] = this->x[6] * (trig_temp0 * trig_temp2 + trig_temp3 * trig_temp5 * trig_temp1) - this->x[7] * (trig_temp5 * trig_temp0 - trig_temp3 * trig_temp2 * trig_temp1) + this->x[8] * trig_temp3 * trig_temp4;
      hx[6] = this->x[9];
      hx[7] = this->x[10];
      hx[8] = this->x[11];
      hx[9] = this->x[12];
      hx[10] = this->x[13];
      hx[11] = this->x[14];

      /* Measurement Jacobian */
      H[0][0] = 1.0;
      H[1][1] = 1.0;
      H[2][2] = 1.0;
      H[3][6] = trig_temp5 * trig_temp4;
      H[3][7] = trig_temp4 * trig_temp2;
      H[3][8] = -trig_temp1;
      H[3][10] = -this->x[8] * trig_temp4 - this->x[6] * trig_temp5 * trig_temp1 - this->x[7] * trig_temp2 * trig_temp1;
      H[3][11] = trig_temp4 * (this->x[7] * trig_temp5 - this->x[6] * trig_temp2);
      H[4][6] = -trig_temp3 * trig_temp2 + trig_temp5 * trig_temp0 * trig_temp1;
      H[4][7] = trig_temp3 * trig_temp5 + trig_temp0 * trig_temp2 * trig_temp1;
      H[4][8] = trig_temp4 * trig_temp0;
      H[4][9] = this->x[8] * trig_temp3 * trig_temp4 - this->x[7] * trig_temp5 * trig_temp0 + this->x[6] * trig_temp0 * trig_temp2 + this->x[6] * trig_temp3 * trig_temp5 * trig_temp1 + this->x[7] * trig_temp3 * trig_temp2 * trig_temp1;
      H[4][10] = trig_temp0 * (-this->x[8] * trig_temp1 + this->x[6] * trig_temp5 * trig_temp4 + this->x[7] * trig_temp4 * trig_temp2);
      H[4][11] = -this->x[6] * trig_temp3 * trig_temp5 - this->x[7] * trig_temp3 * trig_temp2 + this->x[7] * trig_temp5 * trig_temp0 * trig_temp1 - this->x[6] * trig_temp0 * trig_temp2 * trig_temp1;
      H[5][6] = trig_temp0 * trig_temp2 + trig_temp3 * trig_temp5 * trig_temp1;
      H[5][7] = -trig_temp5 * trig_temp0 + trig_temp3 * trig_temp2 * trig_temp1;
      H[5][8] = trig_temp3 * trig_temp4;
      H[5][9] = -this->x[7] * trig_temp3 * trig_temp5 + this->x[6] * trig_temp3 * trig_temp2 - this->x[8] * trig_temp4 * trig_temp0 - this->x[6] * trig_temp5 * trig_temp0 * trig_temp1 - this->x[7] * trig_temp0 * trig_temp2 * trig_temp1;
      H[5][10] = trig_temp3 * (-this->x[8] * trig_temp1 + this->x[6] * trig_temp5 * trig_temp4 + this->x[7] * trig_temp4 * trig_temp2);
      H[5][11] = this->x[6] * trig_temp5 * trig_temp0 + this->x[7] * trig_temp0 * trig_temp2 + this->x[7] * trig_temp3 * trig_temp5 * trig_temp1 - this->x[6] * trig_temp3 * trig_temp2 * trig_temp1;
      H[6][9] = 1.0;
      H[7][10] = 1.0;
      H[8][11] = 1.0;
      H[9][12] = 1.0;
      H[10][13] = 1.0;
      H[11][14] = 1.0;
    };

};

KalmanFilter EKF; // GLOBAL OBJECT 'EKF'

