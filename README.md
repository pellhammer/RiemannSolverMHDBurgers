# Background viscosity dependent Riemann solver for the MHD-Burgers model


## DESCRIPTION:
An exact Riemann solver for the MHD-Burgers model
```math
\begin{pmatrix} u\\v\end{pmatrix}_t + \begin{pmatrix} u\\v\end{pmatrix}_t
'''
such that each shock wave has a viscous profile with respect to the parabolic regularization
```math
\begin{pmatrix} u\\v\end{pmatrix}_t + \begin{pmatrix} u\\v\end{pmatrix}_t = B\begin{pmatrix} u\\v\end{pmatrix}_{xx}
'''
where $B$ is a symmetric and positive definite $2\times 2$ matrix.


## Author
+ [Valentin Pellhammer](http://www.math.uni-konstanz.de/~pellhammer/)  
 Department of Mathematics and Statistics,  
 University of Konstanz,  
 78457 Konstanz, Germany

