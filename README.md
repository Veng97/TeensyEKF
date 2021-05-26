# TeensyEKF
TeensyEKF is a lightweight C/C++ implementation of the Extended Kalman Filter that supports sensors with varying update rates.

The EKF implementation was made with the intention of supporting _fast_ sensor updates on microcontrollers. The core code is essentially build around another implementation "TinyEKF" (https://github.com/simondlevy/TinyEKF). There are only small differences in how the user interacts with the library, and more guides may be found at the previously mentioned repository. Both implementations are efficient for microcontrollers in that they use static (compile-time) memory allocation (no "new" or "malloc").

The difference between the implementations is that 'TeensyEKF' uses an iterative method for including each sensor reading. This means we can perform the update step without need of computing the 'inverse matrix' part of the Kalman gain step i.e. (P * H * P + R)^-1. In order for this method to work, we make the assumption that the measurement noise is uncorrelated, and that each sensor noise component is described solely by its corresponding diagonal term.

![image](https://user-images.githubusercontent.com/40239379/119681705-96ad7080-be42-11eb-9e74-1ce96becbbdc.png)
Image from: http://www.anuncommonlab.com/articles/how-kalman-filters-work/part2.html

The 'Main.ino' file shows an example sketch of how the sensor updates would be included, and the 'kalmanfilter.hpp' file shows how to define a custom EKF for the specific implementation. In the example a 21 state model is used for a combined Inertial Navigation System and Attitude Heading and Reference System, including also the disturbance forces and torques. For the interested reader, and to give an indication of the computational performance the models are shown below:


The model contains the following states and inputs:
![image](https://user-images.githubusercontent.com/40239379/119683535-02dca400-be44-11eb-9d8b-18b9890e375a.png)

The dynamic prediction models for the individual states are derived as follows:
![image](https://user-images.githubusercontent.com/40239379/119683144-b5603700-be43-11eb-9eaa-6db330413de4.png)
![image](https://user-images.githubusercontent.com/40239379/119683418-ec364d00-be43-11eb-9c96-cc9086ae4174.png)
![image](https://user-images.githubusercontent.com/40239379/119683001-9cf01c80-be43-11eb-9ae9-fc97a5dcddbf.png)

Note that covariance terms are defined assuming integrated noise (and the setQ_emperical() function is derived)
![image](https://user-images.githubusercontent.com/40239379/119683362-dfb1f480-be43-11eb-808d-db5a5c875c9c.png)
