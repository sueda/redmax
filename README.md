# REDMAX: Efficient & Flexible Approach for Articulated Dynamics

**NEW**: We now have an _analytically differentiable_ version! See the the [notes](notes.pdf) and the [reference MATLAB implementation](matlab-diff).


### Contents

- `notes.pdf`: An extensive writeup with details on:
  - Maximal coordinates
  - Reduced coordinates
  - Analytical derivatives
  - Implicit integration
  - Adjoint method
- `matlab-diff`: Object-oriented MATLAB implementation of differentiable redmax
  - Fully implicit time integration: BDF1 and BDF2
  - Parameter optimization with the adjoint method
  - Frictional contact with the ground [[Geilinger et al. 2020]](https://arxiv.org/pdf/2007.00987.pdf)
- `matlab-simple`: Simpler object-oriented MATLAB implementation for getting started
- `matlab`: Object-oriented MATLAB implementation with many features, including:
  - Recursive hybrid dynamics (Featherstone's algorithm) for comparison
  - Time integration using `ode45` or `euler`
  - Frictional dynamics with Bilateral Staggered Projections
  - Spline curve and surface joints [[Lee and Terzopoulos 2008]](http://web.cs.ucla.edu/~dt//papers/siggraph08/siggraph08.pdf)
- `c++`: C++ implementation including Projected Block Jacobi Preconditioner


### Citation

**ACM Transactions on Graphics, 38 (4) 104:1-104:10 (SIGGRAPH), 2019.**

[Ying Wang](http://www.yingwang.io/), [Nicholas J. Weidner](http://weidnern.github.io/), [Margaret A. Baxter](https://www.linkedin.com/in/baxter-margareta/), [Yura Hwang](http://yurahwang.com/), [Danny M. Kaufman](http://dannykaufman.io), [Shinjiro Sueda](http://faculty.cs.tamu.edu/sueda/)

	@article{Wang2019,
	  author = {Wang, Ying and Weidner, Nicholas J. and Baxter, Margaret A. and Hwang, Yura and Kaufman, Danny M. and Sueda, Shinjiro},
	  title = {\textsc{RedMax}: Efficient \& Flexible Approach for Articulated Dynamics},
	  year = {2019},
	  issue_date = {July 2019},
	  publisher = {ACM},
	  address = {New York, NY, USA},
	  volume = {38},
	  number = {4},
	  issn = {0730-0301},
	  url = {https://doi.org/10.1145/3306346.3322952},
	  doi = {10.1145/3306346.3322952},
	  journal = {{ACM} Trans.\ Graph.},
	  month = jul,
	  articleno = {104},
	  numpages = {10},
	  keywords = {friction, rigid body dynamics, physical simulation, constraints, contact}
	}

