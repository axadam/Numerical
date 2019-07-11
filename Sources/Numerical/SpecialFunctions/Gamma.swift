//
//  Gamma.swift
//  Numerical
//
//  Created by Adam Roberts on 9/25/18.
//

import Foundation

/// Regularized Incomplete Gamma Function (lower), P(a,x)
///
/// P(a,x) = 𝛾(a,x) / 𝛤(a), 𝛾(a,x) = ∫0..x e^(-t) t^(a-1) dt, a > 0
///
/// For small x (less than a + 1) use the series approach, for large
/// use the continued fraction approach. NR also advises using quadrature
/// for large a but that is not yet implemented.
///
/// Numerical Receipes §6.2
public func p_gamma(_ a: Double, _ x: Double) -> Double {
    switch (a,x) {
    case (...0,_): return .nan
    case (_,..<0): return .nan
    case (_,0): return 0
    case (_,..<(a + 1)): return p_gamma_series(a: a, x: x)
    case (_,_): return 1 - q_gamma_frac(a: a, x: x)
    }
}

/// Regularized Incomplete Gamma Function (upper), Q(a,x)
///
/// Q(a,x) = 𝛤(a,x) / 𝛤(a), 𝛤(a,x) = ∫x..∞ e^(-t) t^(a-1) dt, a > 0
///
/// Implemented simply as complement of lower: Q(a,x) = 1 - P(a,x)
public func q_gamma(_ a: Double, _ x: Double) -> Double {
    switch (a,x) {
    case (...0,_): return .nan
    case (_,..<0): return .nan
    case (_,0): return 1
    case (_,..<(a + 1)): return 1 - p_gamma_series(a: a, x: x)
    case (_,_): return q_gamma_frac(a: a, x: x)
    }
}

/// Series approximation of P(a,x)
///
/// 𝛾(a,x) = e^(-x) x^a Σ0..∞ 𝛤(a) / 𝛤(a + 1 + n) x^n
///
/// Compute the denominator recursively: 𝛤(z + 1) = z 𝛤(z)
///
/// For initial term: 𝛤(a) / 𝛤(a + 1) = 1 / a
///
/// Numerical Receipes §6.2
func p_gamma_series(a: Double, x: Double) -> Double {
    let prefix = exp(a * log(x) - x - lgamma(a))
    let first = 1 / a
    let sum = recursiveSum(indices: 1..., sum0: first, state0: first, update: { i, state in
        let ap = a + Double(i)
        let state1 = state * x / ap
        return (state1, state1)
    }, until: { a, b in a.0 == b.0 }, max_iter: 1_000_000)
    return prefix * sum
}

/// Continued fraction approximation of Q(a,x)
///
/// Q(a,x) = e^(-x) x^a 1 / (1 + x - a -) 1 (1 - a) / (3 + x - a -)  2 (2 - a) / (5 + x - a -)
///
/// This is the even part of the following (converges faster):
///
/// Q(a,x) = e^(-x) x^a 1 / (x +) (1 - a) / (1 +) 1 / (x +) (2 - a) / (2 +) 2 / (x +)
///
/// Numerical Receipes §6.2
func q_gamma_frac(a: Double, x: Double) -> Double {
    let prefix = exp(a * log(x) - x - lgamma(a))
    let b = 1 + x - a
    let c = 1 / Double.leastNormalMagnitude
    let d = 1 / b
    let frac = recursiveProduct(indices: 1..., product0: d, state0: (b: b, c: c, d: d), update: { i, state in
        let an = Double(-i) * (Double(i) - a)
        let (b0, c0, d0) = state
        let b1 = b0 + 2
        let c1 = max(Double.leastNormalMagnitude, b1 + an / c0)
        let d1 = 1 / max(Double.leastNormalMagnitude, b1 + an * d0)
        return (c1 * d1, (b: b1, c: c1, d: d1))
    }, until: { a, b in abs(b.1 - 1) < 1e-15 })
    return prefix * frac
}

/// Derivative of regularized lower incomplete gamma function, P
///
/// e^-x * x^(a-1) / Γ(a) = e^(-x + (a-1) * log(x) - logΓ(a))
public func p_gamma_deriv(a: Double, x: Double) -> Double {
    return exp(-x + (a - 1) * log(x) - lgamma(a))
}

/// Inverse of the lower regularized incomplete gamma P(a,x) function.
/// Gives x such that P(a,x) = p.
///
/// Start with approximation and then use Halley's method to find root of P(a,x) - p.
public func inv_p_gamma(_ a: Double, _ p: Double) -> Double {
    switch (a, p) {
    // handle domain edges
    case (...0,_): return .nan
    case (_,..<0): return .nan
    case (_,1.0.nextUp...): return .nan
    case (_,0): return 0
    case (_,1): return .infinity
        
    // closed form solution when a is 1, quantile is -log(1 - p)
    // only valid when 1 - p doesn't lose precision
    case (1,1e-3...): return -log(1 - p)
        
    // normal case
    case (_,_):
        // initial guess
        let guess = invertGuess(a: a, p: p, q: 1 - p)
        
        // Halley method
        // p_gamma'(x) = e^-x * x^(a-1) / Γ(a)
        // p_gamma''(x) = e^-x (a - x - 1) x^(a-2) / Γ(a)
        // p_gamma''(x) / p_gamma'(x) = e^-x (a - x - 1) x^(a-2) / e^-x x^(a-1)
        //                            = (a - x - 1) / x = (a-1)/x - 1
        let a1 = a - 1
        let lna1 = log(a1)
        let gln: Double = lgamma(a)
        let afac = exp(a1 * (lna1 - 1) - gln)
        let x = rootSecondOrder(guess: guess,
                        xmin: 0,
                        maxIter: 11,
                        f: { x in p_gamma(a, x) - p },
                        f1: { x in
                            switch a {
                            case ...1:
                                return exp( -x + a1 * log(x) - gln)
                            case _:
                                return afac * exp( -(x - a1) + a1 * (log(x) - lna1))
                            }
                        },
                        f2f1: { x in a1 / x - 1 })
        return x
    }
}

/// Inverse of the upper regularized incomplete gamma Q(a,x) function.
/// Gives x such that Q(a,x) = q.
///
/// Start with approximation and then use Halley's method to find root of Q(a,x) - q.
public func inv_q_gamma(_ a: Double, _ q: Double) -> Double {
    switch (a, q) {
    // handle domain edges
    case (...0,_): return .nan
    case (_,..<0): return .nan
    case (_,0): return .infinity
    case (_,1): return 0
    case (_,1...): return .nan
        
    // close form solution when a is 1, quantile is -log(q)
    case (1,_): return -log(q)
        
    // normal case
    case (_,_):
        // initial guess
        let guess = invertGuess(a: a, p: 1 - q, q: q)
        
        // Halley method
        // q_gamma'(x) = -e^-x * x^(a-1) / Γ(a)
        // q_gamma''(x) = -e^-x (a - x - 1) x^(a-2) / Γ(a)
        // q_gamma''(x) / q_gamma'(x) = -e^-x (a - x - 1) x^(a-2) / -e^-x x^(a-1)
        //                            = (a - x - 1) / x = (a-1)/x - 1
        let a1 = a - 1
        let lna1 = log(a1)
        let gln: Double = lgamma(a)
        let afac = exp(a1 * (lna1 - 1) - gln)
        let x = rootSecondOrder(guess: guess,
                                xmin: 0,
                                maxIter: 11,
                                f: { x in q_gamma(a, x) - q },
                                f1: { x in
                                    switch a {
                                    case ...1:
                                        return -exp( -x + a1 * log(x) - gln)
                                    case _:
                                        return -afac * exp( -(x - a1) + a1 * (log(x) - lna1))
                                    }
                                },
                                f2f1: { x in a1 / x - 1 })
        return x
    }
}

/// Provide initial guess for inverse P and Q regularized incomplete gamma functions
///
/// Primarily based on the method describe in "EFFICIENT AND ACCURATE ALGORITHMS FOR THE
/// COMPUTATION AND INVERSION OF THE INCOMPLETE GAMMA FUNCTION RATIOS", Gil, Segura,
/// Temme 2013. Also falls back in one case on an approximation from A & S.
public func invertGuess(a: Double, p: Double, q: Double) -> Double {
    let r = exp( (log(p) + lgamma(1 + a)) / a )
    switch (a,r,q) {
        
    // If a is 1 then we have the closed form -log(q). This works everywhere
    // q is well defined (i.e. not when p is very small). Could probably expand to
    // a small region around 1.
    case (1,_,..<0.999):
        return -log(q)
        
    // When p = q = 1/2 we have an expansion from Temme. Could probably expand
    // to a small region around 1/2
    //
    // x₀ = a (1 - 1/3 a⁻¹ + 8 / 405 a⁻² + 184 / 25515 a⁻³ + 2248 / 344425 a⁻⁴ + ...)
    //
    // Asymptotic Inversion of Incomplete Gamma Functions, Temme 1992, Eq. 6.2
    case (_,_,0.5):
        return a - 1/3 + (8/405)/a + (184/25515)/(a^^2) + (2248/3444525)/(a^^3)
        
    // When p is close to zero and a is relatively small we have an asymptotic expansion
    //
    // x = r + i=2... cᵢ rⁱ,
    // r = (pΓ(a + 1))^(1/a)
    //
    // EFFICIENT AND ACCURATE ALGORITHMS FOR THE COMPUTATION AND INVERSION OF THE INCOMPLETE
    // GAMMA FUNCTION RATIOS, Gil, Segura, Temme 2013, Eq. 3.2 and 3.3
    case (_,..<(0.2 * (1 + a)),_):
        let c2 = 1 / (a + 1)
        let c3 = (3 * a + 5) / (2 * (a+1)^^2 * (a + 2))
        let c4 = (8 * a^^2 + 33 * a + 31) / (3 * (a + 1)^^3 * (a + 2) * (a + 3))
        let c5 = (125 * a^^4 + 1179 * a^^3 + 3971 * a^^2 + 5661 * a + 2888) / (24 * (a + 1)^^4 * (a + 2)^^2 * (a + 3) * (a + 4))
        return r + c2 * r^^2 + c3 * r^^3 + c4 * r^^4 + c5 * r^^5
        
    // When q is close to zero and a is relatively small we have an asymptotic expansion
    //
    // x ~ x₀ - L + b Σ i=1... dᵢ / x₀ⁱ,
    //
    // EFFICIENT AND ACCURATE ALGORITHMS FOR THE COMPUTATION AND INVERSION OF THE INCOMPLETE
    // GAMMA FUNCTION RATIOS, Gil, Segura, Temme 2013, Eq. 2.5 and 3.5
    case (..<10,_,..<(exp(-a / 2) / tgamma(a + 1))):
        let η = eta(a, q)
        let λ = lambda(η)
        let x₀ = a * λ
        let b = 1 - a
        let L = log(x₀)

        let d₁ = L - 1
        let d₂ = (1/2) * (2 + 3 * b - 2 * b * L - 2 * L + L^^2)
        let d₃ = (1/6) * (24 * b * L - 11 * b^^2 - 24 * b - 6 * L^^2 + 12 * L - 12 - 9 * b * L^^2 + 6 * b^^2 * L + 2 * L^^3)
        let d₄ = (1/12) * (72 + 36 * L^^2 + 3 * L^^4 - 72 * L + 162 * b - 168 * b * L - 12 * L^^3 + 25 * b^^3 - 22 * b * L^^3 + 36 * b^^2 * L^^2 - 12 * b^^3 * L + 84 * b * L^^2 + 120 * b^^2 - 114 * b^^2 * L)
        return x₀ - L + b * evaluate_polynomial(poly: [0.0,d₁,d₂,d₃,d₄], z: 1 / x₀)
        
    // When a < 1 the following starting point leads to inexpensive iteration
    // by Halley's method
    //
    // x₀ = (pΓ(a + 1))^(1/a)
    //
    // EFFICIENT AND ACCURATE ALGORITHMS FOR THE COMPUTATION AND INVERSION OF THE INCOMPLETE
    // GAMMA FUNCTION RATIOS, Gil, Segura, Temme 2013, Eq. 3.8
    case (..<1,_,_):
        return pow(p * tgamma(a + 1), 1/a)
        
    // When a is large we have an asymptotic expansion. It depends on q so it is not
    // well defined when q is not (i.e. when p is very small)
    //
    // η(a,q) = η₀(a,q) + ε₁(η₀) / a + ε₂(η₀) / a² + ε₃(η₀) / a³
    //
    // EFFICIENT AND ACCURATE ALGORITHMS FOR THE COMPUTATION AND INVERSION OF THE INCOMPLETE
    // GAMMA FUNCTION RATIOS, Gil, Segura, Temme 2013, Eq. 3.11 and 3.12
    case (_,_,..<0.999):
        let η₀ = eta0(a: a, q: q)
        
        // Use temme 1992 method to get epsilons
        let (ε₁,ε₂,ε₃) = epsilon(η₀: η₀)
        
        let η = evaluate_polynomial(poly: [η₀,ε₁,ε₂,ε₃], z: 1 / a)
        let λ = lambda(η)
        
        return λ * a
        
    // When a is large we have an alternative method from A & S. In practice
    // we use this when p is very small and the previous method won't work
    //
    // Q(p) = 2 P⁻¹(ν / 2, p) = ν  ( 1 - 2/(9ν)  + xp √(2/(9ν))  )³
    //      = 2 P⁻¹(a    , p) = 2a ( 1 - 2/(18a) + xp √(2/(18a)) )³
    //          P⁻¹(a    , p) =  a ( 1 - 1/(9a)  + xp √(1/(9a))  )³
    //
    // Handbook of Mathematical Functions, §26.4.17
    case (_,_,_):
        let xp = p < 0.5 ? qapprox(p: p) : -qapprox(p: q)
        return fmax(1e-3, a * (1 - 1/(9 * a) - xp / (3 * sqrt(a)))^^3)
    }
}

/// Calculate first three εᵢ from η₀
///
/// Method based on Temme 1992 section 5
func epsilon(η₀: Double) -> (Double, Double, Double) {
    switch η₀ {
    case -0.3...0.3:
        let coef1: [Double] = [-1/3, 1/36, 1/1620, -7/6480, 5/18144, -11/382725, -101/16329600]
        let ε₁ = evaluate_polynomial(poly: coef1, z: η₀)
        let coef2: [Double] = [-7.0/405, -7/2592, 533/204120, -1579/2099520, 109/1749600, 10217/251942400]
        let ε₂ = evaluate_polynomial(poly: coef2, z: η₀)
        let coef3: [Double] = [449/102060, -63149/20995200, 29233/36741600, 346793/5290790400, -18442139/130947062400]
        let ε₃ = evaluate_polynomial(poly: coef3, z: η₀)
        return (ε₁,ε₂,ε₃)
    case ..<1000:
        let λ₀ = lambda(η₀)
        let µ  = λ₀ - 1
        
        // temme 1992 eq 3.6
        let f = η₀ / µ
        
        // Temme 2013 eq 3.13
        let ε₁ = log(f) / η₀
        
        // Temme 1992 section 5
        let ε₂ = (12 / η₀^^2 - 12 * f^^2 / η₀^^2 - 12 * f / η₀ - 12 * f^^2 * ε₁ / η₀ - 12 * f * ε₁ - 1 - 6 * ε₁^^2) / (12 * η₀)
        let ε₃ = (-30 / η₀^^4 + 12 * f^^2 * ε₁ / η₀^^3 + 12 * f * ε₁ / η₀^^2 + 24 * f^^2 * ε₁ / η₀ + 6 * ε₁^^3 / η₀ - 12 * f^^2 / η₀^^4 + 60 * f^^3 * ε₁ / η₀^^2 + 31 * f^^2 / η₀^^2 + 72 * f^^3 / η₀^^3 + 42 * f^^4 / η₀^^4 + 18 * f^^3 * ε₁^^2 / η₀ + 6 * f^^2 * ε₁^^2 + 36 * f^^4 * ε₁ / η₀^^3 + 12 * f * ε₁^^2 / η₀ + 12 * f^^2 * ε₁^^2 / η₀^^2 - 12 * ε₁ / η₀^^3 + ε₁ / η₀ + f / η₀ - 12 * f / η₀^^3 + 12 * f^^4 * ε₁^^2 / η₀^^2) / (12 * η₀)
        
        return (ε₁,ε₂,ε₃)
    case _:
        let λ₀ = lambda(η₀)
        let µ  = λ₀ - 1
        
        // temme 1992 eq 3.6
        let f = η₀ / µ

        let ε₁ = log(f) / η₀
        let ε₂ = -1 / (12 * η₀)
        let ε₃ = ε₁ / (12 * η₀^^2)

        return (ε₁,ε₂,ε₃)
    }
}

/// Gamma star function from Temme
///
/// Γ∗(a) = Γ(a) / (√(2π/a) a^a e^(-a))
///
/// = Γ(a) / ( √(2π) e^( -x + ( x - 0.5 ) * log(x) ) )
///
/// or if a >> 0 then Stirling series:
///
/// = ∼ 1 + 1/12a−1+ 1/288a−2 +...
///
/// "EFFICIENT AND ACCURATE ALGORITHMS FOR THE COMPUTATION AND INVERSION OF THE INCOMPLETE
/// GAMMA FUNCTION RATIOS", Gil, Segura, Temme 2012, Eq. 2.5, 2.7
public func gammastar(_ a: Double) -> Double {
    return tgamma(a) / ( sqrt(2 * .pi) * exp((a - 0.5) * log(a) - a))
}

/// Find η from a and q. For relatively small a and q
///
/// q = x^a e^-x / Γ(a + 1)
///
/// = e^(-1/2 a η²) / √(2πa) Γ∗(a)
///
/// q √(2πa) Γ∗(a) = e^(-1/2 a η²)
///
/// -1/2 a η² = log(q √(2πa) Γ∗(a))
///
/// η = √(-2 log(q √(2πa) Γ∗(a)) / a)
///
/// "EFFICIENT AND ACCURATE ALGORITHMS FOR THE COMPUTATION AND INVERSION OF THE INCOMPLETE
/// GAMMA FUNCTION RATIOS", Gil, Segura, Temme 2012, Eq. 2.5
public func eta(_ a: Double, _ q: Double) -> Double {
    return sqrt( -2 * log(q * sqrt(2 * .pi) * gammastar(a)) / a )
}

/// Find η₀ from a and q. Works on wide range of values
///
/// 1/2 erfc(η₀ √(a/2)) = q
///           η₀ √(a/2) = erfc⁻¹(2q)
///                  η₀ = erfc⁻¹(2q) / √(a/2)
///
/// temme 1992, Eq 3.2
public func eta0(a: Double, q: Double) -> Double {
    return invErfC(2 * q) / sqrt(a / 2)
}

/// Finds λ for a given η
///
/// Use Lambert W to solve the following for λ
///
/// η² / 2 = λ - 1 - log(λ)
///
/// -η² / 2 - 1 = log(λ) - λ
///
/// e^(-η² / 2 - 1) = λ e^(-λ)
///
/// -e^(-η² / 2 - 1) = -λ e^(-λ)
///
/// -W[-e^(-η² / 2 - 1)] = λ
///
/// temme 2013 Eq. 2.6
public func lambda(_ η: Double) -> Double {
    let s = 0.5 * η^^2
    let λ: Double = {
        switch η {
        case 0:
            return 1.0
        case ..<(-1):
            // Taylor series of the principle branch of the Lambert W function
            // near 0 with argument e^(-1 - η² / 2)
            //
            // W(x) = x - x² + 3/2 x³ - 8/3 x⁴ + 125/24 x⁵ + ...
            let coef: [Double] = [0, 1, -1, 3.0/2, -8.0/3, 125.0/24]
            return evaluate_polynomial(poly: coef, z: exp(-1 - s) )
        case ..<1:
            // Expansion when η is near zero
            //
            // λ = 1 + η + 1/3 η² + 1/36 η³ - 1/270 η⁴ + 1/4320 η⁵
            //
            // temme 1992, below Eq. 6.1
            // This is also the expansion of the Lambert W function's W₋₁ branch
            let coef: [Double] = [1, 1, 1/3, 1/36, -1/270, 1/4320]
            return evaluate_polynomial(poly: coef, z: η)
        case _:
            // Expansion of the principle branch of the Lambert W function for large values
            // with argument e^(η² / 2 + 1)
            let L₁ = 1 + s
            let L₂ = log(L₁)
            let a₁ = 1.0
            let a₂ = (2 - L₂) / 2
            let a₃ = (6 - 9 * L₂ + 2 * L₂^^2) / 6
            let a₄ = -(-12 + 36 * L₂ - 22 * L₂^^2 + 3 * L₂^^3) / 12
            let a₅ = (60 - 300 * L₂ + 350 * L₂^^2 - 125 * L₂^^3 + 12 * L₂^^4) / 60
            let a₆ = -(-120 + 900 * L₂ - 1700 * L₂^^2 + 1125 * L₂^^3 - 274 * L₂^^4 + 20 * L₂^^5) / 120
            return L₁ + L₂ * evaluate_polynomial(poly: [1,a₁,a₂,a₃,a₅,a₄,a₆], z: 1 / L₁)
        }
    }()
    
    // temme suggests iterating from here for -3.5 < η < -0.03 or 0.03 < η < 40
    // η² / 2 = λ - 1 - log(λ)
    // η² / 2 + log(λ) = λ - 1
    // (η² / 2 + log(λ)) / (λ - 1) = 1
    // λ₁ = λ₀ (η² / 2 + log(λ₀)) / (λ₀ - 1)
    switch η {
    case (-3.5...(-0.03)),(0.03...40):
        let λʹ = recursiveSequence(indices: 1..., initialState: λ, maxIter: 100, update: { i, λ₀ in
            let λ₁ = λ₀ * (s + log(λ₀)) / (λ₀ - 1)
            return λ₁
        }, until: { a, b in b / a - 1 < 1e-8 })
        return λʹ ?? λ
    case _:
        return λ
    }
}

/// Finds η for a given λ
///
/// Finds the root of the following with the sign of λ - 1
///
/// η² / 2 = λ - 1 - log(λ)
///
/// η = s √(2 (λ - 1 - log(λ))), s = sign(λ - 1)
public func eta(_ λ: Double) -> Double {
    return (λ - 1).signum * sqrt(2 * (λ - 1 - log(λ)))
}
