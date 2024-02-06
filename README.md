# Background viscosity dependent Riemann solver for the MHD-Burgers model


## DESCRIPTION:
A family of exact Riemann solvers for the MHD-Burgers model that incorporates the crossing of characteristic speeds and non-classical shock waves.

The MHD-Burgers model
```math
\begin{pmatrix} u\\v\end{pmatrix}_t + \begin{pmatrix} \frac{1}{2}u\\\frac{1}{2}v\end{pmatrix}_x = \begin{pmatrix} 0\\0\end{pmatrix}(\#eq:ideal)
```
is non-strictly hyperbolic. The corresponding system with viscosity reads

```math
\begin{pmatrix} u\\v\end{pmatrix}_t + \begin{pmatrix} \frac{1}{2}u\\\frac{1}{2}v\end{pmatrix}_x = B\begin{pmatrix} u\\v\end{pmatrix}_{xx},\#eq:visc)
```
where $B$ is any symmetric and positive definite $2\times 2$ matrix.

The function `RiemannSolverMHD.m` evaluates the unique(!) solution to the Riemann problem
```math
u(x,0)=\begin{cases}(u_l,v_l),&x<0\\(u_r,v_r),&x>0\end{cases},\quad (u_l,v_l),(u_r,v_r)\in\R^2

```
for \@ref(eq:ideal) such that each shock wave has a visous profile with respect to \@ref(eq:visc).

## Author
+ [Valentin Pellhammer](http://www.math.uni-konstanz.de/~pellhammer/)  
 Department of Mathematics and Statistics,  
 University of Konstanz,  
 78457 Konstanz, Germany

