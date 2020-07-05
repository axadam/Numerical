//
//  Gamma.swift
//  Numerical
//
//  Created by Adam Roberts on 9/25/18.
//

import Foundation

// MARK: Regularized Gamma

/// Regularized Incomplete Gamma Function
///
/// This function gives both the upper and lower regularized gamma functions as a
/// `Probability` value.
///
/// Lower regularized gamma function:
///
/// P(a,x) = 𝛾(a,x) / 𝛤(a),
///
/// 𝛾(a,x) = ∫0..x e^(-t) t^(a-1) dt, a > 0
///
/// Upper regularized gamma function:
///
/// Q(a,x) = 𝛤(a,x) / 𝛤(a),
///
/// 𝛤(a,x) = ∫x..∞ e^(-t) t^(a-1) dt, a > 0
///
/// We split up the domains of computation into four areas according to Temme.
///
/// EFFICIENT AND ACCURATE ALGORITHMS FOR THE COMPUTATION AND INVERSION OF THE INCOMPLETE
/// GAMMA FUNCTION RATIOS, Gil, Segura, Temme 2013, Section 2
public func gamma_reg(_ a: Double, _ x: Double) -> Probability {
    let α = x >= 0.5 ? x : log(0.5) / log(0.5 * x)
    switch (a,x) {
    case (...0,_): return .nan
    case (_,..<0): return .nan
    case (_,0): return .p(0)
    case (12...,(0.3*a)...(2.35*a)):
        let pq = pq_gamma_uniform_asymptotic(a: a, x: x, isLower: a > α)
        return a > α ? .p(pq) : .q(pq)
    case (α...,     _):
        let p = p_gamma_series(a: a, x: x)
        return .p(p)
    case (   _,..<1.5):
        let q = q_gamma_series(a: a, x: x)
        return .q(q)
    case (   _,     _):
        let q = q_gamma_frac(a: a, x: x)
        return .q(q)
    }
}
/// Regularized Incomplete Gamma Function (lower), P(a,x)
///
/// P(a,x) = 𝛾(a,x) / 𝛤(a),
///
/// 𝛾(a,x) = ∫0..x e^(-t) t^(a-1) dt, a > 0
///
/// Split up the domains of computation into four areas according to Temme.
///
/// EFFICIENT AND ACCURATE ALGORITHMS FOR THE COMPUTATION AND INVERSION OF THE INCOMPLETE
/// GAMMA FUNCTION RATIOS, Gil, Segura, Temme 2013, Section 2
public func p_gamma(_ a: Double, _ x: Double) -> Double {
    return gamma_reg(a,x).p
}

/// Regularized Incomplete Gamma Function (upper), Q(a,x)
///
/// Q(a,x) = 𝛤(a,x) / 𝛤(a),
///
/// 𝛤(a,x) = ∫x..∞ e^(-t) t^(a-1) dt, a > 0
///
/// Split up the domains of computation into four areas according to Temme.
///
/// EFFICIENT AND ACCURATE ALGORITHMS FOR THE COMPUTATION AND INVERSION OF THE INCOMPLETE
/// GAMMA FUNCTION RATIOS, Gil, Segura, Temme 2013, Section 2
public func q_gamma(_ a: Double, _ x: Double) -> Double {
    return gamma_reg(a,x).q
}

// MARK: Derivative

/// Derivative of regularized lower incomplete gamma function, P
///
/// e^-x * x^(a-1) / Γ(a) = e^(-x + (a-1) * log(x) - logΓ(a))
public func p_gamma_deriv(a: Double, x: Double) -> Double {
    switch a {
    case 0.5:
        return exp(-x) / (sqrt(x) * sqrt(.pi))
    case 1:
        return exp(-x)
    case _:
        return exp(-x + (a - 1) * log(x) - lgamma(a))
    }
}

// MARK: Inverse

/// Inverse regularized gamma function
///
/// Calculates x such that P(a,x) = p and Q(a,x) = q
///
/// Takes a `Probability` value as argument allowing either very small p or q.
///
/// Start with an approximation and then use Halley's method to find the exact value.
public func inv_gamma_reg(_ a: Double, _ pq: Probability) -> Double {
    switch (a, pq.p, pq.q) {
    // handle domain edges
    case (...0,_,_): return .nan
    case (_,..<0,_): return .nan
    case (_,1.0.nextUp...,_): return .nan
    case (_,0,_): return 0
    case (_,_,0): return .infinity
        
    // closed form solution when a is 1, quantile is -log(q)
    // only valid when q doesn't lose precision (p isn't too small)
    case (1,1e-3...,_): return -log(pq.q)
        
    // normal case
    case (_,_,_):
        // initial guess
        let guess = invertGuess(a: a, p: pq.p, q: pq.q)
        
        // Halley method
        // Derivatives of the lower regularized gamma. Negate for upper.
        // Pʹ(a,x) = e^-x * x^(a-1) / Γ(a)
        // Pʺ(a,x) = e^-x (a - x - 1) x^(a-2) / Γ(a)
        // Pʺ(a,x) / Pʹ(a,x) = e^-x (a - x - 1) x^(a-2) / e^-x x^(a-1)
        //                            = (a - x - 1) / x = (a-1)/x - 1
        let a1 = a - 1
        let lna1 = log(a1)
        let gln: Double = lgamma(a)
        let afac = exp(a1 * (lna1 - 1) - gln)
        let x = root(guess: guess,
                        xmin: 0,
                        maxIter: 11,
                        f: { x in gamma_reg(a, x).difference(pq) },
                        f1: { x in
                            switch a {
                            case ...1:
                                return exp( -x + a1 * log(x) - gln)
                            case _:
                                return afac * exp( -(x - a1) + a1 * (log(x) - lna1))
                            }
                        },
                        f2f1: { x in a1 / x - 1 }).value
        return x
    }

}
/// Inverse of the lower regularized incomplete gamma P(a,x) function.
/// Gives x such that P(a,x) = p.
///
/// Start with approximation and then use Halley's method to find root of P(a,x) - p.
public func inv_p_gamma(_ a: Double, _ p: Double) -> Double {
    return inv_gamma_reg(a, .p(p))
}

/// Inverse of the upper regularized incomplete gamma Q(a,x) function.
/// Gives x such that Q(a,x) = q.
///
/// Start with approximation and then use Halley's method to find root of Q(a,x) - q.
public func inv_q_gamma(_ a: Double, _ q: Double) -> Double {
    return inv_gamma_reg(a, .q(q))
}

// MARK: Implementation

/// Series approximation of P(a,x)
///
/// 𝛾(a,x) = e^(-x) x^a Σ0..∞ 𝛤(a) / 𝛤(a + 1 + n) x^n
///
/// Compute the denominator recursively: 𝛤(z + 1) = z 𝛤(z)
///
/// For initial term: 𝛤(a) / 𝛤(a + 1) = 1 / a
///
/// Numerical Receipes §6.2
fileprivate func p_gamma_series(a: Double, x: Double) -> Double {
    let prefix = exp(a * log(x) - x - lgamma(a))
    let first = 1 / a
    let sum = series(indices: 1..., initialSum: first, initialState: first) { i, state in
        let ap = a + Double(i)
        let state1 = state * x / ap
        return (state1, state1)
    }
    return prefix * sum.value
}

/// Taylor series approximation of Q(a,x)
///
/// Q(a,x) = u + v,
///
/// u = 1 - x^a / Γ(a + 1)
///
///   = 1 - 1 / Γ(a + 1) + (1 - x^a) / Γ(a + 1),
///
/// v = x^a / Γ(a + 1) [1 - Γ(a + 1) 𝛾*(a,x)],
///
/// 𝛾*(a,x) = x^(-a) / Γ(a) 𝛾(a,x)
///
/// A Computational Procedure for Incomplete Gamma Functions, Gautschi 1979, Section 4.1
///
/// "EFFICIENT AND ACCURATE ALGORITHMS FOR THE COMPUTATION AND INVERSION OF THE INCOMPLETE
/// GAMMA FUNCTION RATIOS", Gil, Segura, Temme 2013, Section 2.3
fileprivate func q_gamma_series(a: Double, x: Double) -> Double {
    // u₁ = 1 - 1 / Γ(a + 1)
    let u₁ = -inverse_gamma_p1m1(a)
    
    /// 1 / Γ(a + 1)
    let Γ⁻¹a1 = 1 - u₁
    
    // u₂ = (1 - x^a) / Γ(a + 1)
    //    = -(e^(a log x) - 1) / Γ(a + 1)
    let lnx = log(x)
    let u₂ = -expm1(a * lnx) * Γ⁻¹a1
    
    /// u = 1 - 1 / Γ(a + 1) + (1 - x^a) / Γ(a + 1)
    let u = u₁ + u₂
    
    // v = -x^a Σ i=0... (-x)ⁱ / ((a + i) i!) / Γ(a)
    //   = x^(a+1) / (a+1) Σ i=0... tᵢ / Γ(a),
    // tᵢ = (a + 1) (-x)ⁱ / ((a + i + 1) (i + 1)!)
    //    = -(a + i) x / ((a + i + 1)(i + 1)) tᵢ₋₁,
    //    = -pᵢ tᵢ₋₁ / qᵢ, t₀ = 1
    // pᵢ = (a + i) x = pᵢ₋₁ + x, p₀ = ax
    // qᵢ = (a + i + 1) (i + 1) = qᵢ₋₁ + rᵢ₋₁, q₀ = a + 1
    // rᵢ = a + 2i + 3 = rᵢ₋₁ + 2, r₀ = a + 3
    //
    // A Computational Procedure for Incomplete Gamma Functions, Gautschi 1979, Eq 4.10
    let Σtᵢ = series(indices: 1..., initialSum: 1.0, initialState: (a * x,a + 1,a + 3,1.0)) { i, prev in
        let (pᵢ₋₁, qᵢ₋₁, rᵢ₋₁, tᵢ₋₁) = prev
        let pᵢ = pᵢ₋₁ + x
        let qᵢ = qᵢ₋₁ + rᵢ₋₁
        let rᵢ = rᵢ₋₁ + 2
        let tᵢ = -pᵢ * tᵢ₋₁ / qᵢ
        return (tᵢ, (pᵢ,qᵢ,rᵢ,tᵢ))
    }
    
    /// v = 1 / Γ(a) x^(a + 1) / (a + 1) Σtᵢ
    let v = a * Γ⁻¹a1 * exp((a + 1) * lnx) * Σtᵢ.value / (a + 1)
    
    return u + v
}

/// Continued fraction approximation of Q(a,x)
///
/// Q(a,x) = e^(-x) x^a / 𝛤(a) * 1 / (1 + x - a -) 1 (1 - a) / (3 + x - a -)  2 (2 - a) / (5 + x - a -)
///
/// This is the even part of the following (converges faster):
///
/// Q(a,x) = e^(-x) x^a / 𝛤(a) * 1 / (x +) (1 - a) / (1 +) 1 / (x +) (2 - a) / (2 +) 2 / (x +)
///
/// Numerical Receipes §6.2
fileprivate func q_gamma_frac(a: Double, x: Double) -> Double {
    let prefix = exp(a * log(x) - x - lgamma(a))
    let frac = continued_fraction(
        b0: 0,
        a: { iInt in let i = Double(iInt); return iInt == 1 ? 1 : (i - 1) * (a - (i - 1)) },
        b: { 1 + x - a + 2 * Double($0 - 1) })
    return prefix * frac.value
}

/// Series repesentation of Q(a,x) or P(a,x) when a and x are large
///
/// Q(a,x) = 1/2 erfc(η √(a/2)) + Rₐ(η),
///
/// P(a,x) = 1/2 erfc(-η √(a/2)) - Rₐ(η),
///
/// Rₐ(η) = e^(-1/2 η²a) / √(2πa) Sₐ(η),
///
/// Sₐ(η) = a / (a + β₁) Σi=0... βᵢηⁱ
///
// EFFICIENT AND ACCURATE ALGORITHMS FOR THE COMPUTATION AND INVERSION OF THE INCOMPLETE
// GAMMA FUNCTION RATIOS, Gil, Segura, Temme 2013, Section 2.5
fileprivate func pq_gamma_uniform_asymptotic(a: Double, x: Double, isLower: Bool = true) -> Double {
    /// Sign depends on which tail we want
    let sgn = isLower ? -1.0 : 1.0
    
    /// µ = λ - 1, λ = x / a
    let µ = (x - a) / a
    
    /// 1/2 η² = µ - log(1 + µ), Temme 1979 Eq. 1.3
    let hη² = -log1pmx(µ)
    
    /// η = s √(2 (µ - log(1 + µ)), s = sign(µ)
    let η = µ.signum * sqrt(2 * hη²)
    
    /// u = 1/2 erfc(√(a/2) η)
    let u: Double = 0.5 * erfc(sgn * η * sqrt(a / 2.0))

    /// prefix = e^(-1/2 η²a) / √2πa
    let Rprefix = exp(-hη² * a) / sqrt(2 * .pi * a)
    
    /// βᵢ = 1/a (i + 2) βᵢ₊₂ + dᵢ₊₁
    let β = C.temme_d.enumerated().reversed().scan((βᵢ₊₁: 0.0, βᵢ₊₂: 0.0)) { prev, term in
        let (βᵢ₊₁,βᵢ₊₂) = prev
        let (n   ,dᵢ₊₁) = term
        let i = Double(n - 1)
        let βᵢ = (i + 2) * βᵢ₊₂ / a + dᵢ₊₁
        return (βᵢ₊₁: βᵢ, βᵢ₊₂: βᵢ₊₁)
        }.dropFirst().dropLast().map { $0.βᵢ₊₁ }.reversed()
    
    /// S = a / (a + β₁) Σi=0... βᵢηⁱ
    let S = a / (a + Array(β)[1]) * evaluate_polynomial(poly: Array(β), z: η)
    
    let v = sgn * Rprefix * S
    return u + v
}

/// Provide initial guess for inverse P and Q regularized incomplete gamma functions
///
/// Primarily based on the method describe in "EFFICIENT AND ACCURATE ALGORITHMS FOR THE
/// COMPUTATION AND INVERSION OF THE INCOMPLETE GAMMA FUNCTION RATIOS", Gil, Segura,
/// Temme 2013. Also falls back in one case on an approximation from A & S.
fileprivate func invertGuess(a: Double, p: Double, q: Double) -> Double {
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
        return a - 1.0/3.0 + (8.0/405.0)/a + (184.0/25515.0)/(a^^2) + (2248.0/3444525.0)/(a^^3)
        
    // When p is close to zero and a is relatively small we have an asymptotic expansion
    //
    // x = r + i=2... cᵢ rⁱ,
    // r = (pΓ(a + 1))^(1/a)
    //
    // EFFICIENT AND ACCURATE ALGORITHMS FOR THE COMPUTATION AND INVERSION OF THE INCOMPLETE
    // GAMMA FUNCTION RATIOS, Gil, Segura, Temme 2013, Eq. 3.2 and 3.3
    case (_,..<(0.2 * (1.0 + a)),_):
        let c2 = 1.0 / (a + 1.0)
        let c3 = (3.0 * a + 5.0) / (2.0 * (a+1.0)^^2 * (a + 2.0))
        let c4 = (8.0 * a^^2 + 33.0 * a + 31.0) / (3.0 * (a + 1.0)^^3 * (a + 2.0) * (a + 3.0))
        let c5 = (125.0 * a^^4 + 1179.0 * a^^3 + 3971.0 * a^^2 + 5661.0 * a + 2888.0) / (24.0 * (a + 1.0)^^4 * (a + 2.0)^^2 * (a + 3.0) * (a + 4.0))
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

        let d₁ = L - 1.0
        let d₂ = (1.0/2.0) * (2.0 + 3.0 * b - 2.0 * b * L - 2.0 * L + L^^2)
        let d₃ = (1.0/6.0) * (24.0 * b * L - 11.0 * b^^2 - 24.0 * b - 6.0 * L^^2 + 12.0 * L - 12.0 - 9.0 * b * L^^2 + 6.0 * b^^2 * L + 2.0 * L^^3)
        let d₄ = (1.0/12.0) * (72.0 + 36.0 * L^^2 + 3.0 * L^^4 - 72.0 * L + 162.0 * b - 168.0 * b * L - 12.0 * L^^3 + 25.0 * b^^3 - 22.0 * b * L^^3 + 36.0 * b^^2 * L^^2 - 12.0 * b^^3 * L + 84.0 * b * L^^2 + 120.0 * b^^2 - 114.0 * b^^2 * L)
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
        return fmax(1e-3, a * (1.0 - 1.0/(9.0 * a) - xp / (3.0 * sqrt(a)))^^3)
    }
}

/// Calculate first three εᵢ from η₀
///
/// Method based on Temme 1992 section 5
fileprivate func epsilon(η₀: Double) -> (Double, Double, Double) {
    switch η₀ {
    case -0.3...0.3:
        let coef1: [Double] = [-1.0/3.0, 1.0/36.0, 1.0/1620.0, -7.0/6480.0, 5.0/18144.0, -11.0/382725.0, -101.0/16329600.0]
        let ε₁ = evaluate_polynomial(poly: coef1, z: η₀)
        let coef2: [Double] = [-7.0/405.0, -7.0/2592.0, 533.0/204120.0, -1579.0/2099520.0, 109.0/1749600.0, 10217.0/251942400.0]
        let ε₂ = evaluate_polynomial(poly: coef2, z: η₀)
        let coef3: [Double] = [449.0/102060.0, -63149.0/20995200.0, 29233.0/36741600.0, 346793.0/5290790400.0, -18442139.0/130947062400.0]
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
        let ε₂ = (12.0 / η₀^^2 - 12.0 * f^^2 / η₀^^2 - 12.0 * f / η₀ - 12.0 * f^^2 * ε₁ / η₀ - 12.0 * f * ε₁ - 1.0 - 6.0 * ε₁^^2) / (12.0 * η₀)
        let ε₃ = (-30.0 / η₀^^4 + 12.0 * f^^2 * ε₁ / η₀^^3 + 12.0 * f * ε₁ / η₀^^2 + 24.0 * f^^2 * ε₁ / η₀ + 6.0 * ε₁^^3 / η₀ - 12.0 * f^^2 / η₀^^4 + 60.0 * f^^3 * ε₁ / η₀^^2 + 31.0 * f^^2 / η₀^^2 + 72.0 * f^^3 / η₀^^3 + 42.0 * f^^4 / η₀^^4 + 18.0 * f^^3 * ε₁^^2 / η₀ + 6.0 * f^^2 * ε₁^^2 + 36.0 * f^^4 * ε₁ / η₀^^3 + 12.0 * f * ε₁^^2 / η₀ + 12.0 * f^^2 * ε₁^^2 / η₀^^2 - 12.0 * ε₁ / η₀^^3 + ε₁ / η₀ + f / η₀ - 12.0 * f / η₀^^3 + 12.0 * f^^4 * ε₁^^2 / η₀^^2) / (12.0 * η₀)
        
        return (ε₁,ε₂,ε₃)
    case _:
        let λ₀ = lambda(η₀)
        let µ  = λ₀ - 1
        
        // temme 1992 eq 3.6
        let f = η₀ / µ

        let ε₁ = log(f) / η₀
        let ε₂ = -1 / (12.0 * η₀)
        let ε₃ = ε₁ / (12.0 * η₀^^2)

        return (ε₁,ε₂,ε₃)
    }
}

/// Gamma star function from Temme
///
/// Γ∗(a) = Γ(a) / (√(2π/a) (a/e)^a), a > 0
///
/// = Γ(a) / ( √(2π) e^( -a + ( a - 0.5 ) * log(a) ) )
///
/// if a >> 0 we use the Stirling series:
///
/// = ∼ 1 + 1/12a−1+ 1/288a−2 +...
///
/// Γ∗ tries to capture just the correction term in the Stirling series for Γ:
///
/// Γ(a) = √(2π/a) (a/e)^a Σi...N-1 (cᵢ / aⁱ)
///
/// "EFFICIENT AND ACCURATE ALGORITHMS FOR THE COMPUTATION AND INVERSION OF THE INCOMPLETE
/// GAMMA FUNCTION RATIOS", Gil, Segura, Temme 2013, Eq. 2.5, 2.7
fileprivate func gammastar(_ a: Double) -> Double {
    switch a {
    case ...3:
        return tgamma(a) / ( sqrt(2 * .pi) * exp((a - 0.5) * log(a) - a))
    case    _:
        return evaluate_polynomial(poly: C.stirling, z: 1 / a)
    }
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
/// GAMMA FUNCTION RATIOS", Gil, Segura, Temme 2013, Eq. 2.4
fileprivate func eta(_ a: Double, _ q: Double) -> Double {
    return sqrt( -2.0 * log(q * sqrt(2.0 * .pi) * gammastar(a)) / a )
}

/// Find η₀ from a and q. Works on wide range of values
///
/// 1/2 erfc(η₀ √(a/2)) = q
///
/// η₀ √(a/2) = erfc⁻¹(2q)
///
/// η₀ = erfc⁻¹(2q) / √(a/2)
///
/// temme 1992, Eq 3.2
fileprivate func eta0(a: Double, q: Double) -> Double {
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
fileprivate func lambda(_ η: Double) -> Double {
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
            let a₃ = (6.0 - 9.0 * L₂ + 2.0 * L₂^^2) / 6.0
            let a₄ = -(-12.0 + 36.0 * L₂ - 22.0 * L₂^^2 + 3.0 * L₂^^3) / 12.0
            let a₅ = (60.0 - 300.0 * L₂ + 350.0 * L₂^^2 - 125.0 * L₂^^3 + 12.0 * L₂^^4) / 60.0
            let a₆ = -(-120.0 + 900.0 * L₂ - 1700.0 * L₂^^2 + 1125.0 * L₂^^3 - 274.0 * L₂^^4 + 20.0 * L₂^^5) / 120.0
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
        let λʹ = sequence(first: λ) { λ₀ in
            let λ₁ = λ₀ * (s + log(λ₀)) / (λ₀ - 1)
            return λ₁
        }.until(maxIter: 100) { a, b in b.isApprox(.maybeZero(a), tolerance: .strict) }
        return λʹ?.result ?? λ
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
fileprivate func eta(λ: Double) -> Double {
    return (λ - 1).signum * sqrt(2 * (λ - 1 - log(λ)))
}

/// Finds η for a given µ
///
/// Finds the root of the following with sign of µ:
///
/// η = s √(2 (µ - log(1 + µ)), s = sign(µ)
///
/// THE ASYMPTOTIC EXPANSION OF THE INCOMPLETE GAMMA FUNCTIONS, Temme 1979, Eq. 1.3
fileprivate func eta(µ: Double) -> Double {
    return µ.signum * sqrt(2 * (-log1pmx(µ)))
}

/// Calculates 1 / Γ(x + 1) - 1, uses an expansion when x is small
///
/// 1 / Γ(x) = Σi=1... aᵢxⁱ
///
/// 1 / Γ(x + 1) - 1 = -1 + Σi=1... aᵢ₊₁xⁱ
///
/// Concerning two series for the gamma function, JW Wrench 1967, Eq. 22
fileprivate func inverse_gamma_p1m1(_ x: Double) -> Double {
    switch x {
    case ..<1.5: return evaluate_polynomial(poly: C.wrench, z: x)
    case      _: return 1 / tgamma(x + 1) - 1
    }
}

// MARK: Coefficients

/// Coefficient vectors
///
/// This is fine while we are using Double but needs more thought if we
/// want to go generic. In particular, note that literals don't currently
/// work as expected for types other than Float or Double.
fileprivate struct C {
    /// Stirling series for Γ(a)
    ///
    /// Provides the cᵢ in
    ///
    /// Γ(a) = √(2π/a) (a/e)^a Σi...N-1 (cᵢ / aⁱ)
    ///
    /// Note that this series is not convergent so more terms start to hurt at
    /// some point (where depends on a).
    ///
    /// Concerning two series for the gamma function, JW Wrench 1967, Table 2
    static let stirling: [Double] = [
         1,
         0.08333_33333_33333_33333_33333_33333_33333_33333_33333_33333,
         0.00347_22222_22222_22222_22222_22222_22222_22222_22222_22222,
        -0.00268_13271_60493_88888_88888_88888_88888_88888_88888_88888,
        -0.00022_94720_93621_39917_69547_32510_28806_58444_44444_44444,
         0.00078_40392_21720_06662_74740_34881_44228_88496_96257_10366,
         0.00006_97281_37583_65857_77429_39882_85757_83308_29359_63594,
        -0.00059_21664_37353_69388_28648_36225_60440_11873_91585_19680,
        -0.00005_17179_09082_60592_19337_05784_30020_58822_81785_34534,
         0.00083_94987_20672_08727_99933_57516_76498_34451_98182_11159,
         0.00007_20489_54160_20010_55908_57193_02250_15052_06345_17380,
        -0.00191_44384_98565_47752_65008_98858_32852_25448_76893_57895,
        -0.00016_25162_62783_91581_68986_35123_98027_09981_05872_59193,
         0.00640_33628_33808_06979_48236_38090_26579_58304_01893_93280,
         0.00054_01647_67892_60451_51804_67508_57024_17355_47254_41598,
        -0.02952_78809_45699_12050_54406_51054_69382_44465_65482_82544,
        -0.00248_17436_00264_99773_09156_58368_74346_43239_75168_04723,
         0.17954_01170_61234_85610_76994_07722_22633_05309_12823_38692,
         0.01505_61130_40026_42441_23842_21877_13112_72602_59815_45541,
        -1.39180_10932_65337_48139_91477_63542_27314_93580_45617_72646,
        -0.11654_62765_99463_20085_07340_36907_14796_96789_37334_38371,
    ]

    /// Taylor expansion of 1 / Γ(1 + x) - 1, x < 1.5
    ///
    /// This is modified from the original series for 1 / Γ(x) in two ways: (1)
    /// we remove the first coefficient, thererby dividing the whole series by x
    /// and making it a series for 1 / Γ(1 + x), and (2) we subtract 1 from the
    /// first (constant) term to make it 1 / Γ(1 + x) - 1.
    ///
    /// Concerning two series for the gamma function, JW Wrench 1967, Table 5
    static let wrench: [Double] = [
         0,
         0.57721_56649_01532_86060_65120_90082_4,
        -0.65587_80715_20253_88107_70195_15145_4,
        -0.04200_26350_34095_23552_90039_34875_4,
         0.16653_86113_82291_48950_17007_95102_1,
        -0.04219_77345_55544_33674_82083_01289_2,
        -0.00962_19715_27876_97356_21149_21672_3,
         0.00721_89432_46663_09954_23950_10340_5,
        -0.00116_51675_91859_06511_21139_71084_0,
        -0.00021_52416_74114_95097_28157_29963_1,
         0.00012_80502_82388_11618_61531_98626_3,
        -0.00002_01348_54780_78823_86556_89391_4,
        -0.00000_12504_93482_14267_06573_45359_5,
         0.00000_11330_27231_98169_58823_74128_9,
        -0.00000_02056_33841_69776_07103_45015_9,
         0.00000_00061_16095_10448_14158_17863_4,
         0.00000_00050_02007_64446_92229_30056_2,
        -0.00000_00011_81274_57048_70201_44588_3,
         0.00000_00001_04342_67116_91100_51048_8,
         0.00000_00000_07782_26343_99050_71253_7,
        -0.00000_00000_03696_80561_86422_05708_2,
         0.00000_00000_00510_03702_87454_47597_9,
        -0.00000_00000_00020_58326_05356_65067_9,
        -0.00000_00000_00005_34812_25394_23018_0,
         0.00000_00000_00001_22677_86282_38260_9,
        -0.00000_00000_00000_11812_59301_69745_6,
         0.00000_00000_00000_00118_66922_54751_7,
         0.00000_00000_00000_00141_23806_55318_0,
        -0.00000_00000_00000_00022_98745_68443_6,
         0.00000_00000_00000_00001_71440_63219_3,
         0.00000_00000_00000_00000_01337_35173_1,
        -0.00000_00000_00000_00000_02054_23355_1,
         0.00000_00000_00000_00000_00273_60300_6,
        -0.00000_00000_00000_00000_00017_32356_4,
        -0.00000_00000_00000_00000_00000_23606_0,
         0.00000_00000_00000_00000_00000_18650_0,
        -0.00000_00000_00000_00000_00000_02218_0,
         0.00000_00000_00000_00000_00000_00129_9,
         0.00000_00000_00000_00000_00000_00001_2,
        -0.00000_00000_00000_00000_00000_00001_1,
         0.00000_00000_00000_00000_00000_00000_1
    ]
    
    /// Temme's d coefficients used in the uniform asymptotic expansion
    /// of the incomplete gamma function with large a and x near a.
    ///
    /// They are defined as the coefficients in the following expansion:
    ///
    /// η / (λ - 1) = Σi=0... dᵢ ηⁱ,
    ///
    /// d₀ = -1 / 3, dᵢ = (i + 2) αᵢ₊₂,
    ///
    /// Where the αᵢ are from the expansion for the Lambert W's -1 branch. Temme
    /// says we only need 25 terms when a > 12
    ///
    /// THE ASYMPTOTIC EXPANSION OF THE INCOMPLETE GAMMA FUNCTIONS, Temme 1979,
    /// Eq. 3.8 and following
    ///
    /// Listings of the Lambert W coefficients available as OEIS A005447/A005446
    static let temme_d: [Double] = [
         1,
        -3.33333_33333_33333_33333_33333_33333e-1,
         8.33333_33333_33333_33333_33333_33333e-2,
        -1.48148_14814_81481_48148_14814_81481e-2,
         1.15740_74074_07407_40740_74074_07407e-3,
         3.52733_68606_70194_00352_73368_60670e-4,
        -1.78755_14403_29218_10699_58847_73663e-4,
         3.91926_31785_22437_78169_70409_56300e-5,
        -2.18544_85106_79992_16147_36429_55124e-6,
        -1.85406_22107_15159_96070_17988_36230e-6,
         8.29671_13409_53086_00501_62421_31664e-7,
        -1.76659_52736_82607_93043_60054_24574e-7,
         6.70785_35434_01498_58036_93971_00296e-9,
         1.02618_09784_24030_80425_73957_32273e-8,
        -4.38203_60184_53353_18655_29746_22447e-9,
         9.14769_95822_36790_23418_24881_76331e-10,
        -2.55141_93994_94624_97668_77953_79939e-11,
        -5.83077_21325_50425_06746_40894_50400e-11,
         2.43619_48020_66741_62436_94069_67078e-11,
        -5.02766_92801_14175_58909_05498_59257e-12,
         1.10043_92031_95613_47708_37417_44972e-13,
         3.37176_32624_00985_37882_76988_41692e-13,
        -1.39238_87224_18162_06591_93661_84895e-13,
         2.85348_93807_04744_32039_66909_90528e-14,
        -5.13911_18342_42572_61899_06458_03004e-16,
        -1.97522_88294_34944_28353_96240_15807e-15,
         8.09952_11567_04561_33407_11566_87025e-16,
    ]
}
