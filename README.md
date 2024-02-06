# Background viscosity dependent Riemann solver for the MHD-Burgers model


## Description:
A family of exact Riemann solvers for the MHD-Burgers model that incorporates the crossing of characteristic speeds and non-classical shock waves.
It serves as a builing block for numerical schemes like the Godunov method, the random choice method (Glimm's scheme) or numerical front tracking.

The MHD-Burgers model
```math
\begin{pmatrix} u\\v\end{pmatrix}_t + \begin{pmatrix} \frac{1}{2}u^2\\\frac{1}{2}cv^2\end{pmatrix}_x = \begin{pmatrix} 0\\0\end{pmatrix},\quad c>0 \tag{1}
```
is non-strictly hyperbolic. The corresponding system with viscosity reads

```math
\begin{pmatrix} u\\v\end{pmatrix}_t + \begin{pmatrix} \frac{1}{2}u^2\\\frac{1}{2}cv^2\end{pmatrix}_x = B\begin{pmatrix} u\\v\end{pmatrix}_{xx},\tag{2}
```
where $B\in\mathbb{R}^{2\times 2}$ is any symmetric and positive definite matrix.

The function `RiemannSolverMHD.m` evaluates the unique(!) solution to the Riemann problem
```math
u(x,0)=\begin{cases}(u_l,v_l),&x<0\\(u_r,v_r),&x>0\end{cases},\qquad ((u_l,v_l),(u_r,v_r))\in(\mathbb{R}^2)^2

```
for (1) such that each shock wave has a visous profile with respect to (2). The solution is sensitively dependent on the viscosity matrix $B$ and consists of 
Lax-shocks, rarefaction waves, and non-classical (undercompressive) shock waves. 


## Example







## Author
+ [Valentin Pellhammer](http://www.math.uni-konstanz.de/~pellhammer/)  
 Department of Mathematics and Statistics,  
 University of Konstanz,  
 78457 Konstanz, Germany
