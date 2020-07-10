# REDMAX: Efficient & Flexible Approach for Articulated Dynamics

### ACM Transactions on Graphics, 38 (4) 104:1-104:10 (SIGGRAPH), 2019.

[Ying Wang](http://www.yingwang.io/), [Nicholas J. Weidner](http://weidnern.github.io/), [Margaret A. Baxter](https://www.linkedin.com/in/baxter-margareta/), [Yura Hwang](http://yurahwang.com/), [Danny M. Kaufman](http://dannykaufman.io), [Shinjiro Sueda](http://faculty.cs.tamu.edu/sueda/)

**NEW**: We now have a differentiable version! See the `matlab-diff` folder.

- `c++`: C++ implementation including Projected Block Jacobi Preconditioner
- `matlab`: Object-oriented MATLAB implementation with many features, including:
  - Recursive hybrid dynamics (Featherstone's algorithm) for comparison
  - Time integration using `ode45` or `euler`
  - Frictional dynamics with Bilateral Staggered Projections
  - Spline curve and surface joints [Lee and Terzopoulos 2008]
- `matlab-simple`: Simpler object-oriented MATLAB implementation for getting started
- `matlab-diff`: Object-oriented MATLAB implementation of differentiable redmax
  - Fully implicit time integration: BDF1 and BDF2
- `notes.pdf`: An extensive writeup with details on maximal coordinates, reduced coordinates, and their various derivatives

### Citation

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

