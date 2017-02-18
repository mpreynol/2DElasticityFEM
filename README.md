# 2DElasticityFEM is a 2D elasticity FEM program written in MATLAB. 
Original Written Dec 2016 by Mathew Reynolds.

The supplied nodal (NN) and element arrays (NEL) 
follows a standard format. Included is a Grid Square that creates a 2D square mesh or Gred Rectangle that creates a 2D 
rectangular mesh. 

The boundary conditions are imposed using arrays BE and BN. If an exact integration function handle isn't passed to the 
boundary integration loop that linear intepolation is performed based on nodal values assigned during BN array creation. 

The latest commit includes a script that calls the FE program with succesively smaller mesh sizes to verify convergence. 

Note that the FEM mesh includes a quadlinear 2D element. Theoretically the program should adapt is a higher order shape function
object is used inplace of "quadlinear". This is not tested...

Included in this repo is an analtyical derivation of one element quadlinear elasticity element that was used to verify the program - 
and its gaussian quadraute scheme.
