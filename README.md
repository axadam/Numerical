# Numerical

Numerical computing tools and functions in Swift.

### Basic functions

Fills in a few missing basic functions.

* `expm1mx(_:)` - Calculate e^x - 1 - x with precision even when x is small

* `log1pmx(_:)` - Calculate log(1 + x) - x with precision even when x is small

* `xmsin(_:)` - Calculate x - sin(x) with precision even when x is small

* `signum` - property on `FloatingPoint` that gives a value of same sign and magnitude 1 (or zero if the magnitude is 0)

* `absmax(_:_:)` - Returns whichever of the two arguments has a larger magnitude regardless of sign. This is a generic wrapper of FloatingPoint's `maximumMagnitude` static function.

### Special functions

A collection of functions with many applications in stats/ml and the sciences.

* Error Function

  * `invErf(_:)` - Inverse of the Error function (`erf`)

  * `invErfC(_:)` - Inverse of the complement of the Error function (`erfc`)

* Gamma Function

  * `gammaReg(a:x:)` - The regularized incomplete gamma function. This returns a `Probability` value that can store values close to 0 or 1 with precision.
  
  * `pGamma(a:x:)` - The regularized lower incomplete gamma function Pₐ(x) = γ(a,x) / Γ(a). This is the lower tail of `gammaReg(a:x:)`.

  * `qGamma(a:x:)` - The regularized upper incomplete gamma function Qₐ(x) = Γ(a,x) / Γ(a). This is the upper tail of `gammaReg(a:x:)`.

  * `pGammaDeriv(a:x:)` - The derivative of P(a,x), P'(a,x)

  * `invGammaReg(a:pq:)` - The inverse of the regularized incomplete gamma function such that `gamma_reg(a,x)` = p, where p is a `Probability` value.

  * `invPGamma(a:p:)` - The inverse of P such that P(a,x) = p

  * `invQGamma(a:q:)` - The inverse of Q such that Q(a,x) = q

* Beta Function

  * `beta(a:b:)` - The beta function, Beta(a,b) = Gamma(a) * Gamma(b) / Gamma(a + b)

  * `lbeta(a:b:)` - The log of the beta function

  * `betaReg(x:a:b)` - The regularized incomplete beta function, I_x(a,b). This returns a `Probability` value that can store values close to 0 or 1 with precision.

  * `betaRegDeriv(x:a:b:)` - The derivative of I, Iʹ_x(a,b)

  * `invBetaReg(p:a:b:)` - The inverse of I such that I_x(a,b) = p, where p is a `Probability` value.

* Marcum Q Function

  * `marcum(mu:x:y:)` - The Marcum Q function, Qᵤ(x,y) = e^(-x) Σᵢ₌₀... xⁱ / i! Qᵤ₊ᵢ(y). Output is a `Probability` value that can store values close to 0 or 1 with precision.

  * `marcumDeriv(mu:x:y:)` - The derivative of the Marcum Q function with respect to y.

  * `invMarcum(mu:x:p:)` - The inverse of the Marcum Q function with respect to y such that Qᵤ(x,y) = p, where p is a `Probability` value. 

### Root finding

Functions to find the root of a function of interest.

* `root(guess:xmin?:xmax?:tolerance:method:intercept:f)` - Finds a root of f without any derivatives. It first brackets the root and then finds it within that interval. Multiple methods to find the root are available, including: Secant, Brent's (default), Dekker's, Ridders', and TOMS Algo 748. A non-zero intercept may optionally be specified which then informs the tolerance of the exit condition. This is useful when using the root to invert a function.

* `root(guess:xmin?:xmax?:tolerance:f:f1)` - Finds a root of f using its first derivative. This is the Newton-Raphson method.

* `root(guess:xmin?:xmax?:tolerance:f:f1:f2)` - Finds a root of f using its first and second derivatives. This uses Halley's method.

* `root(guess:xmin?:xmax?:tolerance:f:f1:f2f1)` - Same as above but takes the ratio of the second derivative to the first. For use in cases where the ratio is cheaper to calculate.

### Quadrature

Numerically integrate a function on a closed interval.

* `integrate(range:maxIter:method:f:)` - Integrates the function f on the specified range using the specified method.

  * `trapezoidal` - Trapezoidal Rule: a simple method that is is especially efficient when integrating a periodic function over its period or when integrating a peak function. In other cases it will usually not be the most efficient method.

  * `romberg` - (default) Romberg's Method: this method layers Richardson extrapolation on top of the Trapezoidal Rule to achieve more accuracy for the same number of function evaluations in many cases.

### Types

A few types with use cases in numerical computing.

* `Probability` - A value representing a probability between zero and one. The underlying value may be stored as either the probability itself or as its complement. This allows us to store numbers either very close to zero or very close to one without loss of precision.

* `CountedFunction` - An object to wrap a function and count how many times it is called. Can be useful if you are evaluating an expensive function and have a computation budget. In Swift 5.2 or later the object is callable.

* `IterativeValue` - A protocol for types that need to store a value and how many iterations or evaluations went into computing it.

### Tools

* More flexible stopping criteria for sequences - extensions on Sequence that can prove powerful for defining converence and other stopping criteria for numerical methods

  * `first(where:)` - a variant of the Standard Library function where the closure also has access to the preceding element.
  
  * `until(minIter:maxIter:_)` - finds the first element after minIter and before maxIter satisfying the predicate, or the last element it processed otherwise. It returns an `UntilValue` that tells how many elements it saw and how it exited. Available with the closure seeing only the last element or also seeing the preceding element.

* Approximate Equality - methods for determing if two `FloatingPoint` numbers are close to each other in a way useful for numerical computations. This introduces an `EqualityTarget` type which lets you specify arguments related to determining the scale for relative tolerance, and an `EqualityTolerance` type which lets you specify tolerance both relatively and absolutely. Several preset tolerances are provided as static members to represent best practices of different strictness.

  * `isApprox(_:tolerance:)` - extension on `FloatingPoint` that takes an `EqualityTarget` and an `EqualityTolerance`. Default tolerance is matching about half the digits supported by your type.
  
  * `isWithinULP(of:n:)` - extension on `BinaryFloatingPoint` that lets you specify tolerance as the number of discrete floating point values away another number is. This one has no special handling for zero and should not be used there.

* Summation - Provides some standard summation methods

  * `sumNaive()` - Simple sequential summation. This method is subject to accumulation of rounding and truncation errors. Implemented as extension on Sequence where Element is FloatingPoint.
  
  * `sumPairwise()` - Pairwise summation. As fast as naive sum but accumulates error much more slowly. Implemented as extension on Collection where Element is FloatingPoint and Index is Int.
  
  * `sumKahan()` - Kahan's compensated sum. Estimates the error after adding each term and tries to correct for it. More accurate in many cases than pairwise but slower. Implemented as extension on Sequence where Element is FloatingPoint.
  
  * `sumKBN()` - Kahan-Babuška-Neumaier sum. Improves on Kahan by compensating for rounding error in either the sum or the addend. The same number of operations as Kahan sum and more accurate so KBN should be preferred. Implemented as extension on Sequence where Element is FloatingPoint.

* Series - Some tools around evaluating convergent series and products.

  * `series(indices:initialSum:initialState:maxIter:tolerance:update:)` - Evaluates the potentially infinite series truncated either at its `maxIter` term or when it converges within `tolerance`. Terms are defined by the closure `update` which has access to the index and optionally some State as of the previous term.

  * `product(indices:initialProduct:initialState:maxIter:tolerance:update:)` - Evaluates the potentially infinite product truncated either at its `maxIter` term or when it converges within `tolerance`. Terms are defined by the closure `update` which has access to the index and optionally some State as of the previous term.

* Polynomials - Evaluation of polynomials using Horner's method. Also includes ratios of polynomials.

  * `polynomial(coeffs:z:)` - evaluate the polynomial with coefficients `coeffs` at point `z`
  
  * `polynomialRatio(num:denom:z:)` - evaluate the ratio of polynomials defined by `num` and `denom` at point `z`

* Chebyshev polynomials - Evaluation of Chebyshev polynomials, including rescaling the interval.

  * `chebyshev(coeffs:z:interval:maxTerms:)` - evaluates the Chebyshev polynomial defined by `coeffs` at point `z`. Optionally an interval other than [-1,1] can be specified and/or a maximum number of terms.

* Continued Fractions - Evaluation of continued fractions using the modified Lentz method.

  * `continuedFraction(b0:a:b:maxIter:)` - evaluate a continued fraction of the form b₀ + a₁ / (b₁ +) a₂ / (b₂ +) ... The terms aᵢ are defined by the closure `a` and the terms bᵢ by the closure `b` for i>0. The initial term b₀ is specified by `b0`. A maximum number of levels may be set by `maxIter`.

### Accuracy

Extensive accuracy measurement and testing to make sure accuracy doesn't change from the expected levels. Accuracy is measured in terms of Log Relative Error. Bolded results have achieved full accuracy relative to the reference value or machine precision of the data type.

#### Summation
| Case | Easy | Hard | Peters |
| --- | ---: | ---: | ---: |
| Kahan | __15.0__ | 5.6 | -0.0 |
| Kahan-Babuška-Neumaier | __15.0__ | __15.0__ | __15.0__ |
| Naive | 12.9 | 4.1 | -0.0 |
| Pairwise | __15.0__ | 4.2 | -0.0 |

Easy is 10k random doubles from [0,1] and their negatives shuffled into an array. Hard is 10k random doubles drawn from [0,1] multiplied by random powers of ten between -10 and 10. Peters is the small pathological case [1.0,1e100,1.0,-1e100].

#### Basic Functions Near Zero
| Case | e^x - 1 - x | log(1 + x) - x | x - sin(x) |
| --- | ---: | ---: | ---: |
| 1e-1 | __15.0__ | 14.6 | __15.0__ |
| 1e-10 | 10.5 | 10.5 | __15.0__ |
| 3e-1 | __15.0__ | __15.0__ | __15.0__ |
| 3e-10 | 10.0 | 10.0 | __15.0__ |
| 9e-1 | __15.0__ | __15.0__ | __15.0__ |

#### Root Finding
| Case | ln 5; x₀=2 | sin(x) - x/2; x₀=2 | x³ - 2x - 5; x₀=2 | ∛3647963; x₀=364 |
| --- | ---: | ---: | ---: | ---: |
| Bisection | __15.0__  (3/50*) | __15.0__  (3/49) | __15.0__  (2/50) | __15.0__  (3/50*) |
| Brent | __15.0__  (3/8) | __15.0__  (3/9) | __15.0__  (2/6) | __15.0__  (3/11) |
| Dekker | __15.0__  (3/9) | __15.0__  (3/10) | __15.0__  (2/6) | __15.0__  (3/10) |
| Halley | __15.0__  (5) | __15.0__  (5) | __15.0__  (5) | __15.0__  (6) |
| Newton | __15.0__  (7) | __15.0__  (6) | __15.0__  (6) | __15.0__  (9) |
| Ridders | __15.0__  (3/8) | __15.0__  (3/12) | __15.0__  (2/16) | __15.0__  (3/14) |
| Secant | __15.0__  (3/11) | 0.0  (3/5†) | __15.0__  (2/6) | 7.8  (3/30*) |
| TOMS 748 | __15.0__  (3/11) | __15.0__  (3/8) | __15.0__  (2/7) | __15.0__  (3/11) |

(number of function evaluations in parentheses. for bracketing methods first number in parentheses is how many evaluations to bracket. * indicates method didn't converge. † indicates converged to a different root.)

#### Error Function
| Case | f | f⁻¹ |
| --- | ---: | ---: |
| 1.2345 | __15.0__ | __15.0__ |

#### Regularized Incomplete Gamma Function, reference values from Mathematica
| Case | f | f⁻¹ |
| --- | ---: | ---: |
| (a:4,x:0.7) | __15.0__ | 14.7 |

#### Regularized Incomplete Gamma Function, 'GammaCHI', Gil, Tegura, and Temme 2015, Table 1
| Case | f | f⁻¹ |
| --- | ---: | ---: |
| (a:1e-13,x:0.01) | 12.8 | 12.2 |
| (a:1e-13,x:6.310e-15) | 4.3 | 2.8 |
| (a:1e-13,x:7.110e-7) | 4.0 | 2.9 |
| (a:1e-249,x:0.01) | __15.0__ | 0.0 |
| (a:1e-249,x:6.310e-15) | 4.3 | 0.0 |
| (a:1e-249,x:7.110e-7) | 4.0 | 0.0 |

#### Regularized Incomplete Beta Function
| Case | f | f⁻¹ |
| --- | ---: | ---: |
| (a:3,b:5,x:0.4) | __15.0__ | __15.0__ |

#### Marcum Q Function, reference values from Mathematica
| Case | f | f⁻¹ |
| --- | ---: | ---: |
| µ:11.5,x:15.3,y:23 | 13.3 | 14.0 |
| µ:11.5,x:15.3,y:29 | __15.0__ | __15.0__ |
| µ:25,x:35,y:49 | 13.8 | __15.0__ |
| µ:25,x:35,y:65 | 13.9 | __15.0__ |

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
| µ:1,x:75,y:0.5 | __13.0__ | __13.0__ |
| µ:10,x:100,y:1 | __13.0__ | __13.0__ |
| µ:2,x:100,y:2 | 12.3 | __13.0__ |
| µ:5,x:150,y:30 | __13.0__ | __13.0__ |

#### Marcum Q inverse, tail values. µ: 15.3, x: 11.5, p varies. Reference values from Mathematica
| Case | f⁻¹ |
| --- | ---: |
| p:1e-10 | __6.0__ |
| p:1e-20 | __6.0__ |
| p:1e-30 | __5.0__ |
| q:1e-10 | __6.0__ |
| q:1e-20 | __6.0__ |
| q:1e-30 | 5.4 |

#### Continued Fraction
| Case | LRE |
| --- | ---: |
| e | __15.0__  (16) |
| π | 10.1  (1000*) |
| √2 | __15.0__  (20) |

(number of terms in parentheses. * indicates method didn't converge.)

#### Chebyshev, e^x test case from Approximation Theory Approximation Practice, Trefethen
| Case | exp(x) |
| --- | ---: |
| -0.5 | 14.8 |
| 0.5 | 15.0 |

#### Horner's Method of polynomial evaluation
| Case | LRE |
| --- | ---: |
| -19 + 7x - 4x² + 6x³; x=3 | __15.0__ |

#### Infinite Series
| Case | LRE |
| --- | ---: |
| Chudnovsky algorithm for π | __15.0__  (3) |
| Newton's arcsin series for π | __15.0__  (23) |
| e = Σ i=0... 1 / i! | __15.0__  (19) |
| ln 2 = 2 Σ i=0... 3⁻²ⁱ⁻¹ / (2i + 1) | __15.0__  (16) |

(number of terms in parentheses. * indicates method didn't converge.)

#### Infinite Product
| Case | LRE |
| --- | ---: |
| Viète's formula for 2/π | __15.0__  (26) |
| Wallis product for π/2 | 2.6  (100*) |
| Wallis product truncated at n=100 | __15.0__  (100*) |

(number of terms in parentheses. * indicates method didn't converge.)

#### Quadrature
| Case | Romberg | Trapezoidal |
| --- | ---: | ---: |
| 1/x; [1,100] | __15.0__  (16385) | 9.8  (1048577*) |
| e^(-x²) / √π; [-10,10] | 14.7  (1025*) | __15.0__  (129) |
| e^cos(θ); [0,2π] | __15.0__  (1025) | __15.0__  (33) |
| x; [0,5000] | __15.0__  (9) | __15.0__  (9) |
| x³; [0,1] | __15.0__  (9) | 6.0  (1025*) |
| π = 4 ∫0...1 1 / (1 + t²) dt; [0,1] | __15.0__  (257) | 7.3  (1025*) |
| √(1 - 0.36sin²θ) / √(2π); [0,2π] | __15.0__  (1025) | __15.0__  (65) |

(number of function evaluations in parentheses. * indicates method didn't converge.)
