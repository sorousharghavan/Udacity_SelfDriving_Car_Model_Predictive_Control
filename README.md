# CarND-Controls-MPC
## Self-Driving Car Engineer Nanodegree Program
## _By: Soroush Arghavan_

---

_Build instructions at the bottom_

# Car is now able to reach 101 MPH!

## Lap Stats
_Lap time: 31.4 s_
_Top speed: 101.06 MPH_

![Race Lap](./lap.gif)

## The Model

The model that is used to represent the vehicle takes into account the coordinates _x_ and _y_, the heading psi and velocity _v_. The actuators that are used to control the vehicle are steering angle and throttle. Steering angle is a value between -25 and 25 (scaled to [-1, 1] for the simulator) and the throttle lies within [-1, 1]. For each actuation parameter an error is defined as the deviation from the center (cte) and deviation from the heading (epsi) of the road. To update the values of the prediction for each step, each parameter at step _t+1_ is calculated based on the state of the vehicle at the current step _t_ as follows:

```c++
AD<double> f0 = coeffs[0] + coeffs[1] * x0 + coeffs[2] * CppAD::pow(x0, 2); //+ coeffs[3] * CppAD::pow(x0, 3);
AD<double> psides0 = CppAD::atan(coeffs[1] + 2 * coeffs[2] * x0); //+ 3 * coeffs[3] * CppAD::pow(x0, 2));

// Recall the equations for the model:
// x_[t] = x[t-1] + v[t-1] * cos(psi[t-1]) * dt
// y_[t] = y[t-1] + v[t-1] * sin(psi[t-1]) * dt
// psi_[t] = psi[t-1] + v[t-1] / Lf * delta[t-1] * dt
// v_[t] = v[t-1] + a[t-1] * dt
// cte[t] = f(x[t-1]) - y[t-1] + v[t-1] * sin(epsi[t-1]) * dt
// epsi[t] = psi[t] - psides[t-1] + v[t-1] * delta[t-1] / Lf * dt
fg[1 + x_start + t] = x1 - (x0 + v0 * CppAD::cos(psi0) * dt);
fg[1 + y_start + t] = y1 - (y0 + v0 * CppAD::sin(psi0) * dt);
fg[1 + psi_start + t] = psi1 - (psi0 - v0 * delta / Lf * dt);
fg[1 + v_start + t] = v1 - (v0 + a * dt);
fg[1 + cte_start + t] =
  cte1 - ((f0 - y0) + (v0 * CppAD::sin(epsi0) * dt));
fg[1 + epsi_start + t] =
  epsi1 - ((psi0 - psides0) - v0 * delta0 / Lf * dt);
```

The cost of the predictions are evaluated to allow the vehicle drive as smoothly and quickly as possible reaching a top speed of more than 101 MPH. To achieve this, certain considerations have been taken into account such as penalizing braking and defining an inverse cost relationship between vehicle speed and road curvature which can be found in MPC.cpp from line 119.

## Timestep Length and Elapsed Duration (N & dt)

A variety of timesteps and durations were chosen for testing. It was found that if _N * dt_ is too low (less than 0.5 second), the vehicle would overcompensate cte at low speeds and oscillate about the cte to instability. If _N * dt_ is too high, the curve would deviate from the road as the prediction period is too long at high speeds and would cause the vehicle to "miss" the turn. Furthermore, if N is too high, the computational power needed for each update would introduce a lag that makes the car unstable. Also, if N is too low, the MPC curve becomes too volatile for the car to be able to drive. If dt is too high, the actuation steps would become too slow for the vehicle to use. Also, if dt is too low, the model would try to overcompensate again at low velocities and stop the car from reaching stable speeds.

For the purpose of having a stable race car around the track, a total prediction horizon of 0.84 second with 12 steps of length 0.07 second have been chosen. For a regular cruise speed of 40-60, other parameters such as 10 steps of 0.1 second proved to be most efficient. With said settings at cruise speed, the car is capable of finishing a lap without any braking at set speed.

## Polynomial Fitting and MPC Preprocessing

Second order polynomial fitting is chosen for the purpose of having a race car. Third order polynomials were proven to result in spline and curvy shapes in the fit and therefore adding unnecessary actuation resulting in a harsh ride. Second order polynomials result in smoother steering due to their profile and constant curvature. Furthermore, the curvature of the fit is also used for perfecting the cost function to define areas of the track in which the car should drive more cautiously.

As for preprocessing, the state data from the vehicle is transformed into the vehicle coordinates before curve fitting. This results in _x_, _y_ and psi being zero and making computations much simpler.

## Model Predictive Control with Latency

In order to compensate for the latency of the data coming into the model, in the kinematics models, the actuation data are introduced into the model with a delay of one timestep. This can be found in lines 84-95 of MPC.cpp

## Dependencies

* cmake >= 3.5
 * All OSes: [click here for installation instructions](https://cmake.org/install/)
* make >= 4.1
  * Linux: make is installed by default on most Linux distros
  * Mac: [install Xcode command line tools to get make](https://developer.apple.com/xcode/features/)
  * Windows: [Click here for installation instructions](http://gnuwin32.sourceforge.net/packages/make.htm)
* gcc/g++ >= 5.4
  * Linux: gcc / g++ is installed by default on most Linux distros
  * Mac: same deal as make - [install Xcode command line tools]((https://developer.apple.com/xcode/features/)
  * Windows: recommend using [MinGW](http://www.mingw.org/)
* [uWebSockets](https://github.com/uWebSockets/uWebSockets)
  * Run either `install-mac.sh` or `install-ubuntu.sh`.
  * If you install from source, checkout to commit `e94b6e1`, i.e.
    ```
    git clone https://github.com/uWebSockets/uWebSockets 
    cd uWebSockets
    git checkout e94b6e1
    ```
* Fortran Compiler
  * Mac: `brew install gcc` (might not be required)
  * Linux: `sudo apt-get install gfortran`. Additionall you have also have to install gcc and g++, `sudo apt-get install gcc g++`. Look in [this Dockerfile](https://github.com/udacity/CarND-MPC-Quizzes/blob/master/Dockerfile) for more info.
* [Ipopt](https://projects.coin-or.org/Ipopt)
  * Mac: `brew install ipopt`
       +  Some Mac users have experienced the following error:
       ```
       Listening to port 4567
       Connected!!!
       mpc(4561,0x7ffff1eed3c0) malloc: *** error for object 0x7f911e007600: incorrect checksum for freed object
       - object was probably modified after being freed.
       *** set a breakpoint in malloc_error_break to debug
       ```
       This error has been resolved by updrading ipopt with
       ```brew upgrade ipopt --with-openblas```
       per this [forum post](https://discussions.udacity.com/t/incorrect-checksum-for-freed-object/313433/19).
  * Linux
    * You will need a version of Ipopt 3.12.1 or higher. The version available through `apt-get` is 3.11.x. If you can get that version to work great but if not there's a script `install_ipopt.sh` that will install Ipopt. You just need to download the source from the Ipopt [releases page](https://www.coin-or.org/download/source/Ipopt/) or the [Github releases](https://github.com/coin-or/Ipopt/releases) page.
    * Then call `install_ipopt.sh` with the source directory as the first argument, ex: `sudo bash install_ipopt.sh Ipopt-3.12.1`. 
  * Windows: TODO. If you can use the Linux subsystem and follow the Linux instructions.
* [CppAD](https://www.coin-or.org/CppAD/)
  * Mac: `brew install cppad`
  * Linux `sudo apt-get install cppad` or equivalent.
  * Windows: TODO. If you can use the Linux subsystem and follow the Linux instructions.
* [Eigen](http://eigen.tuxfamily.org/index.php?title=Main_Page). This is already part of the repo so you shouldn't have to worry about it.
* Simulator. You can download these from the [releases tab](https://github.com/udacity/self-driving-car-sim/releases).
* Not a dependency but read the [DATA.md](./DATA.md) for a description of the data sent back from the simulator.

## Basic Build Instructions

1. Clone this repo.
2. Make a build directory: `mkdir build && cd build`
3. Compile: `cmake .. && make`
4. Run it: `./mpc`.
