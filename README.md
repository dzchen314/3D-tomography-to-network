# 3D-tomography-to-network
This project imports tomographic slices of hydrogel particles, separates and labels them using watershed segmentation, and tracks contacts and deformation to obtain force networks and fabric information.
This code is written in Matlab. For questions, contact me at ender314@gmail.com

File "raw_input_image.jpg" shows a typical optical slice of our system, with noise and line artifacts from scattering.
The first part of the Matlab code imports these images and filters them separately using FFT filters, with end result being a binarized image "filtered_binary_image.jpg". Optionally, a non-local means filter works well to smooth the image noise. Such filters can be found in other Matlab repositories.

The code outputs particle files for each reconstructed particle, including particle center of mass and surface voxels. The latter portions of the code uses this information to calculate particle contacts (overlaps) and corresponding contact networks (fabric tensor) and stress tensors from contact areas using contact mechanics.

To obtain particle trajectories (time-series data), one can use particle tracking code, such as the code developed by John Crocker, David Grier, and Eric Weeks: http://www.physics.emory.edu/faculty/weeks//idl/
