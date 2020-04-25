# Numerical

Numerical computing tools and functions in Swift.

### Basic functions

Fills in a few missing basic functions.

* `expm1mx(_:)` - Calculate e^x - 1 - x with precision even when x is small

* `log1pmx(_:)` - Calculate log(1 + x) - x with precision even when x is small

* `xmsin(_:)` - Calculate x - sin(x) with precision even when x is small

Provides some standard summation methods

  * `sum_naive()` - Simple sequential summation. This method is subject to accumulation of rounding and truncation errors. Implemented as extension on Sequence where Element is FloatingPoint.
  
  * `sum_pairwise()` - Pairwise summation. As fast as naive sum but accumulates error much more slowly. Implemented as extension on Collection where Element is FloatingPoint and Index is Int.
  
  * `sum_kahan()` - Kahan's compensated sum. Estimates the error after adding each term and tries to correct for it. More accurate in many cases than pairwise but slower. Implemented as extension on Sequence where Element is FloatingPoint.
  
  * `sum_kbn()` - Kahan-Babuška-Neumaier sum. Improves on Kahan by compensating for rounding error in either the sum or the addend. The same number of operations as Kahan sum and more accurate so KBN should be preferred. Implemented as extension on Sequence where Element is FloatingPoint.
  
### Special functions

A collection of functions with many applications in stats/ml and the sciences.

* Error Function

  * `invErf(_:)` - Inverse of the Error function (`erf`)

  * `invErfC(_:)` - Inverse of the complement of the Error function (`erfc`)

* Gamma Function

  * `p_gamma(a:x:)` - The regularized lower incomplete gamma function Pₐ(x) = γ(a,x) / Γ(a)

  * `q_gamma(a:x:)` - The regularized upper incomplete gamma function Qₐ(x) = Γ(a,x) / Γ(a)

  * `p_gamma_deriv(a:x:)` - The derivative of P(a,x), P'(a,x)

  * `inv_p_gamma(a:p:)` - The inverse of P such that P(a,x) = p

* Beta Function

  * `beta(a:b:)` - The beta function, Beta(a,b) = Gamma(a) * Gamma(b) / Gamma(a + b)

  * `lbeta(a:b:)` - The log of the beta function

  * `beta_reg(x:a:b)` - The regularized incomplete beta function, I_x(a,b)

  * `beta_reg_deriv(x:a:b:)` - The derivative of I, Iʹ_x(a,b)

  * `inv_beta_reg(p:a:b:)` - The inverse of I such that I_x(a,b) = p

* Marcum Q Function

  * `marcum(µ:x:y:)` - The Marcum Q function, Qᵤ(x,y) = e^(-x) Σᵢ₌₀... xⁱ / i! Qᵤ₊ᵢ(y). Output is a `Probability` value that can store values close to 0 or 1 with precision.

  * `marcum_deriv(µ:x:y:)` - The derivative of the Marcum Q function with respect to y.

  * `inv_marcum(µ:x:p:)` - The inverse of the Marcum Q function with respect to y such that Qᵤ(x,y) = p, where p is a `Probability` value. 

### Root finding

Functions to find the root of a function of interest.

* `root(guess:xmin?:xmax?:tolerance:method:f)` - Finds a root of f without any derivatives. It first brackets the root and then finds it within that interval. Multiple methods to find the root are available, including: Brent's (default), Dekker's, Ridders', and TOMS Algo 748.

* `root(guess:xmin?:xmax?:xtol:f:f1)` - Finds a root of f using its first derivative. This is the Newton-Raphson method.

* `root(guess:xmin?:xmax?:xtol:f:f1:f2)` - Finds a root of f using its first and second derivatives. This uses Halley's method.

* `root(guess:xmin?:xmax?:xtol:f:f1:f2f1)` - Same as above but takes the ratio of the second derivative to the first. For use in cases where the ratio is cheaper to calculate.

### Tools

* Series - Some tools around evaluating convergent series. Included is a general way to test if a sequence is converging and then functions specifically for sums and product series where the terms are recursively defined.

* Polynomials - Evaluation of polynomials using Horner's method. Also includes ratios of polynomials.

* Chebyshev polynomials - Evaluation of Chebyshev polynomials, including rescaling the interval.

* Continued Fractions - Evaluation of continued fractions using the modified Lentz method.

* Quadrature - Numerical integration of functions by quadrature. Currently only implements the Trapezoidal Rule.

### Accuracy

Extensive accuracy measurement and testing to make sure accuracy doesn't change from the expected levels. Accuracy is measured in terms of Log Relative Error.

#### Summation
| Case | Kahan | Kahan-Babuška-Neumaier | Naive | Pairwise |
| --- | ---: | ---: | ---: | ---: |
| Easy | 15.0 | 15.0 | 12.9 | 15.0 |
| Hard | 5.6 | 15.0 | 4.1 | 4.2 |
| Peters | -0.0 | 15.0 | -0.0 | -0.0 |

Easy is 10k random doubles from [0,1] and their negatives shuffled into an array. Hard is 10k random doubles drawn from [0,1] multiplied by random powers of ten between -10 and 10. Peters is the small pathological case [1.0,1e100,1.0,-1e100].

#### Root Finding
| Case | ln 5 | ∛3647963 |
| --- | ---: | ---: |
| Bisection | 15.0 | 14.5 |
| Brent | 15.0 | 15.0 |
| Dekker | 15.0 | 15.0 |
| Halley | 15.0 | 15.0 |
| Newton | 15.0 | 15.0 |
| Ridder | 15.0 | 15.0 |
| Secant | 15.0 | 4.8 |
| TOMS 748 | 15.0 | 15.0 |

#### Error Function
| Case | f | f⁻¹ |
| --- | ---: | ---: |
| 1.2345 | 15.0 | 15.0 |

#### Regularized Incomplete Gamma Function, reference values from Mathematica
| Case | f | f⁻¹ |
| --- | ---: | ---: |
| (a:4,x:0.7) | 15.0 | 14.7 |

#### Regularized Incomplete Gamma Function, 'GammaCHI', Gil, Tegura, and Temme 2015, Table 1
| Case | f | f⁻¹ |
| --- | ---: | ---: |
| (a:1e-13,x:0.01) | 12.8 | 12.2 |
| (a:1e-13,x:6.310e-15) | 4.3 | 2.8 |
| (a:1e-13,x:7.110e-7) | 4.0 | 2.9 |
| (a:1e-249,x:0.01) | 15.0 | 0.0 |
| (a:1e-249,x:6.310e-15) | 4.3 | 0.0 |
| (a:1e-249,x:7.110e-7) | 4.0 | 0.0 |

#### Regularized Incomplete Beta Function
| Case | f | f⁻¹ |
| --- | ---: | ---: |
| (a:3,b:5,x:0.4) | 15.0 | 15.0 |

#### Marcum Q Function, reference values from Mathematica
| Case | f | f⁻¹ |
| --- | ---: | ---: |
| µ:11.5,x:15.3,y:23 | 13.3 | 14.0 |
| µ:11.5,x:15.3,y:29 | 11.0 | 11.4 |
| µ:25,x:35,y:49 | 13.8 | 14.6 |
| µ:25,x:35,y:65 | 14.0 | 14.7 |

#### Marcum Q large µ. µ = 8192, y = 1.05µ, and x is a fraction of µ. 'Computation of the Marcum Q-function', Gil, Segura, Temme 2013, Table 6.1
| Case | f | f⁻¹ |
| --- | ---: | ---: |
| x:0.01µ | 7.0 | 9.5 |
| x:0.02µ | 7.9 | 10.3 |
| x:0.03µ | 9.2 | 11.5 |
| x:0.04µ | 9.6 | 11.8 |
| x:0.05µ | 10.0 | 11.8 |
| x:0.06µ | 9.5 | 11.6 |
| x:0.07µ | 9.3 | 11.5 |
| x:0.08µ | 8.1 | 10.5 |
| x:0.09µ | 7.2 | 9.7 |
| x:0.10µ | 6.5 | 9.1 |

#### Marcum Q far left tail, 'GammaCHI', Gil, Tegura, and Temme 2015, Table 3
| Case | f | f⁻¹ |
| --- | ---: | ---: |
| µ:1,x:75,y:0.5 | 13.0 | 13.0 |
| µ:10,x:100,y:1 | 13.0 | 13.0 |
| µ:2,x:100,y:2 | 12.3 | 13.0 |
| µ:5,x:150,y:30 | 13.0 | 13.0 |

#### Marcum Q inverse, tail values. µ: 15.3, x: 11.5, p varies. Reference values from Mathematica
| Case | f⁻¹ |
| --- | ---: |
| p:1e-10 | 6.0 |
| p:1e-20 | 6.0 |
| p:1e-30 | 5.0 |
| q:1e-10 | 6.0 |
| q:1e-20 | 6.0 |
| q:1e-30 | 5.4 |

#### Continued Fraction
| Case | LRE |
| --- | ---: |
| e | 15.0 |
| π | 13.1 |
| √2 | 15.0 |

#### Chebyshev, e^x test case from Approximation Theory Approximation Practice, Trefethen
| Case | exp(x) |
| --- | ---: |
| -0.5 | 14.8 |
| 0.5 | 15.0 |
