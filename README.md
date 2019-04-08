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

Functions to find the root of a function of interest.

* `root(guess:xmin?:xmax?:epsilon:method:f)` - Finds a root of f without any derivatives. It first brackets the root and then finds it within that interval. Multiple methods to find the root are available, including: Brent's (default), Dekker's, Ridders'.

* `root(guess:xmin?:xmax?:xtol:f:f1)` - Finds a root of f using its first derivative. This is the Newton-Raphson method.

* `root(guess:xmin?:xmax?:xtol:f:f1:f2)` - Finds a root of f using its first and second derivatives. This uses Halley's method.

* `root(guess:xmin?:xmax?:xtol:f:f1:f2f1)` - Same as above but takes the ratio of the second derivative to the first. For use in cases where the ratio is cheaper.

### Series

Some tools around evaluating convergent series. Included is a general way to test if a sequence is converging
and then functions specifically for sums and product series where the terms are recursively defined.

### Notes

Care is taken to cite sources for algorithms and also to spell out derivations in the comments. While the source algorithms
are generally performant this code has not been optimized for performance.

