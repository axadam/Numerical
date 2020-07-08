//
//  Sums.swift
//  Numerical
//
//  Created by Adam Roberts on 7/8/20.
//

public extension Sequence where Element: FloatingPoint {
    /// Sum a sequence of numbers
    ///
    /// Straightforward calculation of:
    ///
    /// Σ1...n aᵢ
    ///
    /// This method is subject to accumulation of rounding and truncation errors.
    func sum_naive() -> Element { return reduce(Element.zero, +) }
    
    /// Kahan Summation
    ///
    /// A compensated sum. After adding each term we check for error in low order digits
    /// and compensate. S2 is an estimate of the error of the last time S was rounded and
    /// truncated. The Neumaier variant, `sum_kbn`, is more accurate at the same cost
    /// and should be preferred.
    ///
    /// Sᵢ = Sᵢ₋₁ + (S2ᵢ₋₁ + Yᵢ),
    ///
    /// S2ᵢ = (Sᵢ₋₁ - Sᵢ) + (S2ᵢ₋₁ + Yᵢ)
    ///
    /// "Further remarks on reducing truncation errors", Kahan, 1965
    func sum_kahan() -> Element {
        return reduce((s: Element.zero, s2: Element.zero)) { accum, yᵢ in
            let (sᵢ₋₁,s2ᵢ₋₁) = accum
            let s2py = s2ᵢ₋₁ + yᵢ
            let sᵢ = sᵢ₋₁ + s2py
            let s2ᵢ = (sᵢ₋₁ - sᵢ) + s2py
            return (sᵢ,s2ᵢ)
        }.s
    }

    /// Kahan-Babuška-Neumaier Sum
    ///
    /// A variant of compensated sum that attends to the fact that sometimes the loss of
    /// precision is on the side of the sum and sometimes on the side of the term. Compensation
    /// is applied at the end.
    ///
    /// s₀ = 0, w₀ = 0
    ///
    /// sᵢ = aᵢ + sᵢ₋₁,
    ///
    /// wᵢ = wᵢ₋₁ + (aᵢ + (sᵢ₋₁ - sᵢ)), |aᵢ| ≤ |sᵢ₋₁|
    ///
    /// wᵢ = wᵢ₋₁ + (sᵢ₋₁ + (aᵢ - sᵢ)),  |aᵢ| > |sᵢ₋₁|
    ///
    /// Σ1...n aᵢ = s_n + w_n
    ///
    /// "Rundungsfehleranalyse einiger Verfahren zur Summation endlicher Summen", Neumaier, 1974, Eq. 2.IV
    func sum_kbn() -> Element {
        let (s,w) = reduce((s: Element.zero,w: Element.zero)) { accum, aᵢ in
            let (sᵢ₋₁,wᵢ₋₁) = accum
            let sᵢ = sᵢ₋₁ + aᵢ
            let Δ = abs(aᵢ) <= abs(sᵢ₋₁) ? (sᵢ₋₁ - sᵢ) + aᵢ : (aᵢ - sᵢ) + sᵢ₋₁
            let wᵢ = wᵢ₋₁ + Δ
            return (s: sᵢ, w: wᵢ)
        }
        return s + w
    }
}

public extension Collection where Index == Int, Element: FloatingPoint {
    /// Pairwise sum
    ///
    /// A divide and conquer algorithm where we divide the collection in two and then recursively
    /// apply the algorithm to each half. For collections below a certain size we revert to naive
    /// summation to avoid recurisve overhead.
    ///
    /// Pairwise summation is a special case of superblock accumulation where the number of
    /// blocks at each step, b, is equal to 2.
    func sum_pairwise(N: Int = 100) -> Element {
        switch count {
        case ..<N:
            return sum_naive()
        case    _:
            // Integer division to find index of mid-point
            let m = count / 2 + startIndex
            return self[startIndex...m].sum_pairwise(N: N) + self[(m+1)...(endIndex - 1)].sum_pairwise(N: N)
        }
    }
}
