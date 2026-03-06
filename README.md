# Stiff 9x9 System (Implicit Euler + Newton)

Numerical solution of a stiff 9-dimensional initial value problem on $[0,100]$, using implicit Euler and Newton iterations with LU factorization.

## Problem Statement

Consider the IVP:

```math
y_1'(t) = 10^7\,(y_2 - 2y_1) - e^{y_1}
```

```math
y_i'(t) = 10^7\,(y_{i+1} - 2y_i + y_{i-1}) - e^{y_i}, \quad i=2,\dots,8
```

```math
y_9'(t) = 10^7\,(y_8 - 2y_9) - e^{y_9}
```

with initial data:

```math
y_i(0)=\sin(0.1\pi i), \quad i=1,\dots,9.
```

Goal: approximate

```math
y(100)=(y_1(100),\dots,y_9(100))
```

and estimate the approximation error by comparing two time steps.

## Numerical Method

Because the system is stiff, the project uses:

- **Implicit Euler** for time integration
- **Newton method** at each time step to solve the nonlinear implicit update
- **LU factorization with pivoting** for linear solves inside Newton

For the implicit step, solve:

```math
F(Y)=Y-hf(Y)-y^n=0
```

At each Newton iteration:

```math
J(Y^{(k)})\,\delta^{(k)}=-F(Y^{(k)}), \qquad Y^{(k+1)}=Y^{(k)}+\delta^{(k)}.
```

## Error Estimate

Two integrations are performed:

- $h_1=0.1$
- $h_2=0.05$

Then:

```math
\text{diff}=y^{(h_1)}(100)-y^{(h_2)}(100), \qquad E_2=\|\text{diff}\|_2.
```

## Repository Structure

- `main2.cpp`  
  Runs the two integrations, computes error norm, prints points $(i/10, y_i(100))$.
- `funzioni2.h`  
  Constants, declarations, and global data for implicit stepping.
- `funzioni2.cpp`  
  LU decomposition, linear solver, Newton iteration, residual/Jacobian, time integration.

## Build

```bash
g++ -std=c++17 -O2 main2.cpp funzioni2.cpp -o stiff
```

## Run

```bash
./stiff
```

## Numerical Output (current implementation)

With the current code and tolerances:

- Estimated error norm between $h_1=0.1$ and $h_2=0.05$: `0.0000000000` (to printed precision)

Printed polyline points $(i/10, y_i(100))$:

- (0.1, -0.0000004500)
- (0.2, -0.0000008000)
- (0.3, -0.0000010500)
- (0.4, -0.0000012000)
- (0.5, -0.0000012500)
- (0.6, -0.0000012000)
- (0.7, -0.0000010500)
- (0.8, -0.0000008000)
- (0.9, -0.0000004500)

## Notes

- The Jacobian is assembled from the residual of the implicit Euler step.
- This project is focused on numerical stability and solver workflow for stiff systems.
