# psr-HI-kinematics
This repository provides Python programs to calculate HI kinematic distances for pulsars, based on the Reid et al. (2019) Galactic rotation curve and input radial velocities. Please see Romero-Ruiz & Ocker (2026) for details and cite our research note if you use this software.

## Usage
Command-line: `python h1_kinematic_distance.py`
The script asks for users to input the source name, Galactic coordinates, velocity bounds, and velocity error (typically assumed to be about 7 km/s). 

The script accepts input velocities at the tangent point. Velocities are automatically scaled down by a factor of 1.6 for sightlines through Perseus (100 <= l <= 150 deg). 

An optional plot is output showing distance vs. radial velocity for the sightline.

## Dependencies
The command-line program requires `numpy` and the `kinematic_functions` module (which should be installed in the same directory as `h1_kinematic_distance.py`).




