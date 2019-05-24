We provide a reference implementation of the RedMax algorithm, written in object-oriented MATLAB (2018b). This code is not designed for performance but is rather designed for pedagogical purposes. To run the code, go to the directory containing testRedMax.m and type:

```
>> testRedMax(1,0)
```
This will show a swinging chain with alternating revolute/fixed joints. The first two arguments modify the integrator type and the scene ID:

	(1) `itype`: Integrator type
		- `1` Use the recursive O(n) algorithm [Kim and Pollard 2011; Kim 2012] with ode45.
		- `2` Use RedMax with ode45. This gives numerically the same solution as the recursive O(n) algorithm.
		- `3` Use RedMax with linearly implicit Euler. This is the most feature-rich option.
	(2) `sceneID`: There are many preset scenes. The first few are:
		- `0` Simple serial chain
		- `1` Different revolute axes
		- `2` Branching
		- `3` Spherical joint
		- `4` Loop
		- `5` ···

The full list of scenes is in testRedMax.m. The other arguments are used to control what gets displayed:

	(3) `drawScene`: Whether to draw the scene. Thiscansignificantly speed up the program, since drawing in MATLAB is very slow.
	(4) `plotH`: Whether to compute and plot the energy over time.

To run all the test cases, copy and paste the following into the command window:

```
clear; clc;
for itype = 1 : 3
  for sceneID = 0 : 33
    testRedMax(itype,sceneID,false,false);
  end
end
```
