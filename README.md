# 3D-tomography-to-network
This project imports tomographic slices of hydrogel particles, separates and labels them, and tracks contacts and deformation to obtain force networks and fabric information.
This code is written in Matlab. For questions, contact me at ender314@gmail.com

file "raw_input_image.jpg" shows a typical optical slice of our system, with noise and line artifacts from scattering.
The first part of the Matlab code imports these images and filteres them separately using FFT filters, with end result being a binarized image "filtered_binary_image.jpg". Optionally, a non-local means filter works well to smooth the image noise. Such filters can be found in other Matlab repositories.
