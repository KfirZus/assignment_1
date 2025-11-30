# Ex1 – Polynomial Functions in Java

This project implements a set of numerical and polynomial operations as part of the Intro to Computer Science course (Ariel University).

The code provides:
- Polynomial evaluation
- Polynomial addition and multiplication
- Polynomial equality
- Numerical derivative
- Root finding (bisection method)
- Area between two polynomials
- Curve length approximation
- Polynomial interpolation from points
- Parsing and printing polynomial string format

---

##  Features

### 1. **Evaluate Polynomial**
`f(double[] p, double x)`  
Computes the value of a polynomial at point `x`.

Polynomial representation:
p[0] = c0 (constant)
p[1] = c1 (coefficient of x)
p[2] = c2 (coefficient of x^2)

---

### 2. **Root Finding**
`root_rec(p, x1, x2, eps)`  
Finds a root of the polynomial using the **bisection method**.  
Requires sign difference between f(x1) and f(x2).

---

### 3. **Find Intersection Between Two Polynomials**
`sameValue(p1, p2, x1, x2, eps)`  
Finds a point where `p1(x) = p2(x)` numerically.

---

### 4. **Area Between Polynomials**
`area(p1, p2, x1, x2, N)`  
Approximates area between curves using the trapezoid method.

---

### 5. **Curve Length**
`length(p, x1, x2, segments)`  
Approximates arc length over a segment.

---

### 6. **Derivative**
`derivative(p)`  
Returns the polynomial derivative.

---

### 7. **Add & Multiply**
`add(p1, p2)`  
`mul(p1, p2)`  
Basic polynomial operations.

---

### 8. **Polynomial From Points**
`PolynomFromPoints(xx, yy)`  
Finds a polynomial (degree ≤ 2) passing exactly through given points.

---

### 9. **String <-> Polynomial Conversion**
- `poly(double[])` → Convert polynomial to a human-readable string
- `getPolynomFromString(String)` → Parse the string format back to coefficients

Example:
"-1.0x^2 + 3.0x + 2.0" <--> {2.0, 3.0, -1.0}

---

## Author

Kfir Zusman (כפיר זוסמן) 
207450503
Ariel University — Intro to CS


