We provide a reference implementation of the RedMax algorithm, written in object-oriented MATLAB (2018b). This code is not designed for performance but is rather designed for pedagogical purposes. To run the code, go to the directory containing testRedMax.m and type:

```
>> testRedMax(0)
```
This will show a swinging chain with alternating revolute/fixed joints (Fig. 1a). The argument switches the scene to be simulated:

  0. Simple serial chain
  1. Different revolute axes
  2. Branching
