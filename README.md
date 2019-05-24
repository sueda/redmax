# REDMAX: Efficient & Flexible Approach for Articulated Dynamics

#### Ying Wang, Nicholas J. Weidner, Margaret A. Baxter, Yura Hwang, Danny M. Kaufman, Shinjiro Sueda

### ACM Transactions on Graphics, 38 (4) 104:1-104:10 (SIGGRAPH), 2019.

Source folders:
- `c++`: C++ implementation including Projected Block Jacobi Preconditioner
- `matlab-simple`: Simpler object-oriented MATLAB implementation for getting started
- `matlab`: Object-oriented MATLAB implementation with many features, including:
  - Recursive hybrid dynamics (Featherstone's algorithm) for comparison
  - Time integration using `ode45` or `euler`
  - Frictional dynamics with Bilateral Staggered Projections
  - Spline curve and surface joints [Lee and Terzopoulos 2008]
