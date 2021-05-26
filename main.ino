#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <SPI.h>  // SPI library for communicating with the IMU:

// Include custom EKF class
#include "kalmanfilter.hpp"

// Variables for computation time inspection
unsigned elapsed_serial, elapsed_imu, elapsed_ekf, elapsed_control = 0;
elapsedMicros loop_clock; // Teensy specific clock-type


//-------------------------------------
// SETUP
//-------------------------------------
void setup() {

  Serial.begin(2000000); // opens serial port and sets baud rate
  while (!Serial) delay(10); // wait for Serial initialization

  // Sensor code here (examples)
  IMU.init(); // state com with IMU
  RPi.init(); // start com with Raspberry Pi

};




//-------------------------------------
// LOOP
//-------------------------------------
void loop() {
  loop_clock = 0; // Reset timer

  // Sensor Update
  // -----------------------------
  
  elapsed_serial = RPi.update();
  elapsed_imu = IMU.update();


  
  // EKF Update
  // -----------------------------

  // 1) Data acquisition from IMU
  if (IMU.update_received()) {
    // Insert IMU measurements into the EKF
    for (byte i = 0; i < 3; i++) {
      EKF.Z[i + 3] = IMU.acceleration[i];
      EKF.Z[i + 6] = IMU.euler[i];
      EKF.Z[i + 9] = IMU.omega[i];
      EKF.Z_update[i + 3] = true;
      EKF.Z_update[i + 6] = true;
      EKF.Z_update[i + 9] = true;
    }
  }

  // 2) Data acquisition from an Intel T265 Camera
  if (RPi.update_received()) {
    // Insert POSE measurements into the EKF
    for (byte i = 0; i < 3; i++) {
      EKF.Z[i] = RPi.T265_pose[i];
      EKF.Z_update[i] = true;
    }
  }

  // Perform the EKF step
  elapsed_ekf = EKF.update();

  // Fetch estimated states (double) structures
  double* x_outer = EKF.getXref(0); // [x,y,z,dx,dy,dz,ddx,ddy,ddz]
  double* x_inner = EKF.getXref(9); // [phi,theta,psi,p,q,r]
  double* x_dist = EKF.getXref(15); // [Fx,Fy,Fz,Tx,Ty,Tz] - Parameters estimated by the EKF
  
  
  
  // CONTROL
  // -----------------------------

  // Refresh control loops
  elapsed_control = Controller.update(x_outer, x_inner); 
  // Example usage for a function defined within the Controller class:
  // Controller::update(double x_outer[6], double x_inner[6], double x_dist[6]) {
  //
  //    /* control code here */
  //    ...
  //    ...
  //
  // }
  
  
  print_loop_timers();

};

// Prints execution time for tasks in the main loop
void print_loop_timers() {
  if (elapsed_serial + elapsed_imu + elapsed_ekf + elapsed_control > 0) {
    Serial.print("Serial:");
    Serial.print(elapsed_serial);
    Serial.print(",IMU:");
    Serial.print(elapsed_imu);
    Serial.print(",EKF:");
    Serial.print(elapsed_ekf);
    Serial.print(",Control:");
    Serial.print(elapsed_control);
    Serial.print(",LOOP:");
    Serial.println(loop_clock);
  }
};
