# CarND-Controls-MPC
**Liang Zhang**

Self-Driving Car Engineer Nanodegree Program

---
## The Model: state, actuators and update equations

I used a kinematic (bicycle) model for this project, which ignores tire forces, gravity and mass. This model has four states: `x`, `y`, `psi` and `v`, where `x` and `y` denote car position in 2D coordinates, `psi` denotes heading direction and `v` is velocity. The two actuators are steering wheel `delta` and throttle `a` (an unified actuator for throttle pedal and break pedal). The detailed update equation from `t` to `t+1` is the following. Note `Lf` is the distance between the center of mass of the vehicle and the front wheels.

```
      x_[t+1] = x[t] + v[t] * cos(psi[t]) * dt
      y_[t+1] = y[t] + v[t] * sin(psi[t]) * dt
      psi_[t+1] = psi[t] + v[t] / Lf * delta[t] * dt
      v_[t+1] = v[t] + a[t] * dt
```
In addition, I used two errors as states based on the kinematic model. These two errors are cross-track error `cte` and orientation error `epsi`. The two new update equations are in the following.   
```
      cte[t+1] = y[t] - f(x[t]) + v[t] * sin(epsi[t]) * dt
      epsi[t+1] = psi[t] - psides[t] + v[t] * delta[t] / Lf * dt
```

## Polynomial Fitting and MPC Preprocessing
The waypoints is first transformed from map coordinates to vehicle coordinates (line 97-102 in main.cpp) for displaying and calculating CTE and orientation errors. Then, the waypoints is fitted by a third order polynomial using polyfit function in main.cpp.

The steer value is divided by deg2rad(25) to scale down in between [-1, 1]. Note that the simulater uses an opposite sign for tuning left and tuning right compared to the kinematic model. Hence, the steer value will multiply -1 in the update equations.

## Timestep Length and Elapsed Duration (N & dt)
The timestep length and elapsed duration I tried are in the following table.

#| N | dt
---:|:---:|:---
1   | 20 | 0.05
2   | 10 | 0.1
3   | 8  | 0.125
4   | 5  | 0.2

Firstly, I chose 10 seconds as the duration for prediction. I found it is reasonable for speed from 40 mph to 70 mph. For timestep size dt, I used 0.05, 0.1, 0.125 and 0.2. With a small dt, the accuracy of the predicted trajectory could be very high. But I found for a optimization problem, a high resolution could lead to amplification of modelling error and measurement error. Hence, I take N = 5 and dt = 0.2 for reference speed = 70 mph. I found in this case, a large dt leads to a much more stable (robust) solution. 

## Tune the cost function

In MPC.cpp line 45 - 59, I use the following cost function. 
```
    for (size_t t = 0; t < N; t++) {
      fg[0] += 10 * CppAD::pow(vars[cte_start + t], 2);
      fg[0] += 10 * CppAD::pow(vars[epsi_start + t], 2);
      fg[0] += CppAD::pow(vars[v_start + t] - ref_v, 2);
    }
    // Minimize the use of actuators.
    for (size_t t= 0; t < N - 1; t++) {
      fg[0] += 300 * CppAD::pow(vars[delta_start + t], 2);
      fg[0] += CppAD::pow(vars[a_start + t], 2);
    }
    // Minimize the value gap between sequential actuations.
    for (size_t t= 0; t < N - 2; t++) {
      fg[0] += 5000 * CppAD::pow(vars[delta_start + t + 1] - vars[delta_start + t], 2);
      fg[0] += CppAD::pow(vars[a_start + t + 1] - vars[a_start + t], 2);
    }
```
I found that the weight for cross track error and orientation error should not be too large. Because the system is too sensitive in this case. A large weight for the gap between sequential steering is needed to smooth the steering actuator. 

## Model Predictive Control with Latency
I use the vehicle model to compute the 'latency state' starting from the current state for the duration of the latency. Then the 'latency state' is used as the input state for MPC. See the detailed update equations in the following (line 120-135 in main.cpp).
```
        double pred_px = 0.0 + v *  dt; 
        const double pred_py = 0.0; 
        double pred_psi = 0.0 - v *  delta / Lf * dt;
        double pred_v = v + a * dt;
        double pred_cte = cte + v * sin(epsi) * dt;
        double pred_epsi = epsi- v * delta / Lf * dt;
```

## Dependencies

* cmake >= 3.5
 * All OSes: [click here for installation instructions](https://cmake.org/install/)
* make >= 4.1(mac, linux), 3.81(Windows)
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
    Some function signatures have changed in v0.14.x. See [this PR](https://github.com/udacity/CarND-MPC-Project/pull/3) for more details.

* **Ipopt and CppAD:** Please refer to [this document](https://github.com/udacity/CarND-MPC-Project/blob/master/install_Ipopt_CppAD.md) for installation instructions.
* [Eigen](http://eigen.tuxfamily.org/index.php?title=Main_Page). This is already part of the repo so you shouldn't have to worry about it.
* Simulator. You can download these from the [releases tab](https://github.com/udacity/self-driving-car-sim/releases).
* Not a dependency but read the [DATA.md](./DATA.md) for a description of the data sent back from the simulator.


## Basic Build Instructions

1. Clone this repo.
2. Make a build directory: `mkdir build && cd build`
3. Compile: `cmake .. && make`
4. Run it: `./mpc`.

## Tips

1. It's recommended to test the MPC on basic examples to see if your implementation behaves as desired. One possible example
is the vehicle starting offset of a straight line (reference). If the MPC implementation is correct, after some number of timesteps
(not too many) it should find and track the reference line.
2. The `lake_track_waypoints.csv` file has the waypoints of the lake track. You could use this to fit polynomials and points and see of how well your model tracks curve. NOTE: This file might be not completely in sync with the simulator so your solution should NOT depend on it.
3. For visualization this C++ [matplotlib wrapper](https://github.com/lava/matplotlib-cpp) could be helpful.)
4.  Tips for setting up your environment are available [here](https://classroom.udacity.com/nanodegrees/nd013/parts/40f38239-66b6-46ec-ae68-03afd8a601c8/modules/0949fca6-b379-42af-a919-ee50aa304e6a/lessons/f758c44c-5e40-4e01-93b5-1a82aa4e044f/concepts/23d376c7-0195-4276-bdf0-e02f1f3c665d)
5. **VM Latency:** Some students have reported differences in behavior using VM's ostensibly a result of latency.  Please let us know if issues arise as a result of a VM environment.