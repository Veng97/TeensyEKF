# TeensyEKF
TeensyEKF is a lightweight C/C++ implementation of the Extended Kalman Filter that supports sensors with varying update rates.

The EKF implementation was made with the intention of supporting _fast_ sensor updates on microcontrollers. The core code is essentially build around another implementation "TinyEKF" (https://github.com/simondlevy/TinyEKF). There are only small differences in how the user interacts with the library, and more guides may be found at the previously mentioned repository. Both implementations are efficient for microcontrollers in that they use static (compile-time) memory allocation (no "new" or "malloc").

The difference between the implementations is that 'TeensyEKF' uses an iterative method for including each sensor reading. This means we can perform the update step without need of computing the 'inverse matrix' part of the Kalman gain step i.e. \((P*H*P + R)^-1\)

<img src="https://render.githubusercontent.com/render/math?math=e^{i \pi} = -1">



In order to make it practical for running on Arduino, STM32, and other microcontrollers, it uses static (compile-time) memory allocation (no "new" or "malloc"). The examples folder includes an Arduino example of sensor fusion. The extras/python folder includes an abstract Python class that you can use to prototype your EKF before implementing it in C/C++. The extrasc/c folder contains a "pure C" example from the literature.

Arduino users can simply install or drag the whole TinyEKF folder into their Arduino libraries folder. The examples/SensorFusion folder contains a little sensor fusion example using a BMP180 barometer and LM35 temperature sensor. I have run this example on an Arduino Uno and a Teensy 3.2. The BMP180, being an I^2C sensor, should be connected to pins 4 (SDA) and 5 (SCL) of the Uno, or pins 18 (SDA) and 19 (SCL) of the Teensy. For other Arduino boards, consult the documentation on the Wire library. The analog output from the LM35 should go to the A0 pin of your Arduino or Teensy.
