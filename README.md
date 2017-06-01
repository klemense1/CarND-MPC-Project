# CarND-Controls-MPC
[![Udacity - Self-Driving Car NanoDegree](https://s3.amazonaws.com/udacity-sdc/github/shield-carnd.svg)](http://www.udacity.com/drive)
---

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
* [uWebSockets](https://github.com/uWebSockets/uWebSockets) == 0.14, but the master branch will probably work just fine
  * Follow the instructions in the [uWebSockets README](https://github.com/uWebSockets/uWebSockets/blob/master/README.md) to get setup for your platform. You can download the zip of the appropriate version from the [releases page](https://github.com/uWebSockets/uWebSockets/releases). Here's a link to the [v0.14 zip](https://github.com/uWebSockets/uWebSockets/archive/v0.14.0.zip).
  * If you have MacOS and have [Homebrew](https://brew.sh/) installed you can just run the ./install-mac.sh script to install this.
* [Ipopt](https://projects.coin-or.org/Ipopt)
  * Mac: `brew install ipopt --with-openblas`
  * Linux
    * You will need a version of Ipopt 3.12.1 or higher. The version available through `apt-get` is 3.11.x. If you can get that version to work great but if not there's a script `install_ipopt.sh` that will install Ipopt. You just need to download the source from the Ipopt [releases page](https://www.coin-or.org/download/source/Ipopt/) or the [Github releases](https://github.com/coin-or/Ipopt/releases) page.
    * Then call `install_ipopt.sh` with the source directory as the first argument, ex: `bash install_ipopt.sh Ipopt-3.12.1`.
  * Windows: TODO. If you can use the Linux subsystem and follow the Linux instructions.
* [CppAD](https://www.coin-or.org/CppAD/)
  * Mac: `brew install cppad`
  * Linux `sudo apt-get install cppad` or equivalent.
  * Windows: TODO. If you can use the Linux subsystem and follow the Linux instructions.
* [Eigen](http://eigen.tuxfamily.org/index.php?title=Main_Page). This is already part of the repo so you shouldn't have to worry about it.
* Simulator. You can download these from the [releases tab](https://github.com/udacity/CarND-MPC-Project/releases).



## Basic Build Instructions


1. Clone this repo.
2. Make a build directory: `mkdir build && cd build`
3. Compile: `cmake .. && make`
4. Run it: `./mpc`.

## Implementation

#### State
As a state vector, I used x, y, psi (orientation) and v (velocity). To achieve a control that stays on the street with the right orientation, I have to minimize the following errors
* Cross Track Error cte1 = (f0 - y0) + v0 * sin(epsi0) * dt;
* Orientation Error epsi1 = (psi0 - psides0) + v0 * delta0 * dt / Lf;

Thus, the state vector is extended by the two errors cte and epsi. By predicting how the error is going to evolve to the next time step given a certain actuation, one can find the actuation minimizing the error over time.

#### Actuator
As a actuator, I used theta (steering angle) and a (acceleration combining throttle and breaking)

#### Model
I am using a kinematic bycicle model assuming both axes as one wheel and beeing able to only steer the front wheel.

The discrete update equations for the state vector are the following:
* x1 = x0 + v0 * cos(psi0) * dt
* y1 = y0 + v0 * sin(psi0) * dt
* psi1 = psi0 + v0 * dt delta0 / Lf
* v1 = v0 + a0 * dt

The road is given as xy-coordinates and is fitted as a 3rd order polynom to allow for an evaluation at the points we are predicting in the MPC.

#### Cost function
Costs and constraints (e.g. kinematic constraints for steering or throttle) are stored in the vector fg.

The cost functions contains
* cte
* epsi
* speed difference to reference speed
* use of actuation (steering and throttle)
* change of actuation (change of steering and throttle)

#### Horizon
The time horizon should only be a second or so, as the environment will change too much behind that to allow for a meaningfull prediction. I tried a couple of combinations betweeen number of points (N in the code) and time step (dt in the code) and finally settled with N=20 and dt=0.03. I assume the model to be quite robust to work with slightly different parameters as well.

#### Simulator
The simulator returns the state of the vehicle in global coordinates. I then transform them to a car-reference system, as it makes it easier to calculate the errors epsi and cte.

To allow for an actuation delay (which is modelled as a 100ms sleep of the thread), I am measuring the time the opimization and actuation takes (elapsed_sec in the code) and solve the actuation for the current state vector
* x = v*elapsed_sec
* y = 0
* psi = 0
* v = v
* cte = f  + v * sin(epsi) * mpc.dt_
* epsi = psi_des + v * delta * mpc.dt_ / mpc.Lf_

The mpc problem is then solved.
