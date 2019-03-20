# Numerical

Numerical computing tools and functions in Swift.

### Special functions

A collection of functions with many applications in stats/ml and the sciences.

* Error Function

  * `invErf()` - Inverse of the Error function (`erf`)

  * `invErfC()` - Inverse of the complement of the Error function (`erfc`)

* Gamma Function

  * `p_gamma(a:x:)` - The regularized lower incomplete gamma function P(a,x)

  * `q_gamma(a:x:)` - The regularized upper incomplete gamma function Q(a,x)

  * `p_gamma_deriv(a:x:)` - The derivative of P(a,x), P'(a,x)

  * `inv_p_gamma(a:p:)` - The inverse of P such that P(a,x) = p

* Beta Function

  * `beta(a:b:)` - The beta function, Beta(a,b) = Gamma(a) * Gamma(b) / Gamma(a + b)

  * `lbeta(a:b:)` - The log of the beta function

  * `beta_reg(x:a:b)` - The regularized incomplete beta function, I_x(a,b)

  * `beta_reg_deriv(x:a:b:)` - The derivative of I, I'_x(a,b)

  * `inv_beta_reg(p:a:b:)` - The inverse of I such that I_x(a,b) = p

### Root finding

Functions to find the root of a function of interest. These come in two main types: methods where you know the derivative
of your function and methods where you do not. This initial release only includes Halley's method, which is good functions where
you can cheaply get both the first and second derivatives.

### Series

Some tools around evaluating convergent series. Included is a general way to test if a sequence is converging
and then functions specifically for sums and product series where the terms are recursively defined.

### Notes

Care is taken to cite sources for algorithms and also to spell out derivations in the comments. While the source algorithms
are generally performant this code has not been optimized for performance.

