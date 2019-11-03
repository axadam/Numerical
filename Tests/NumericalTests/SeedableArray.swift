//
//  SeedableArray.swift
//  NumericalTests
//
//  Created by Adam Roberts on 11/3/19.
//

import Foundation

/// A single draw of a random number
///
/// Uses `drand48()` which can be seeded. Optionally accepts a range, otherwise by default
/// will draw from [0,1].
func draw(range: ClosedRange<Double> = 0.0...1.0) -> Double {
    return range.lowerBound + drand48() * (range.upperBound - range.lowerBound)
}

/// A reproducible array of random numbers
///
/// Requires a seed to give reproducibiity. The intended use is for producing test cases. Should not
/// be used for cryptography. Can choose whether to include negative numbers and the range of
/// exponents. By default it gives an array drawn from [-1,1].
func randomArray(seed: Int, n: UInt, exponentRange: ClosedRange<Double> = 0.0...0.0, signed: Bool = true) -> [Double] {
    srand48(seed)
    
    let c: [Double] = (0..<n).map { _ in
        let sign = signed ? (drand48() > 0.5 ? 1.0 : -1.0) : 1.0
        let significand = drand48()
        let exponent = pow(10.0, draw(range: exponentRange).rounded())
        return sign * significand * exponent
    }
    
    return c
}

extension Array {
    /// Reproducible shuffle
    ///
    /// Requires a seed to give reproducibility. The intended use is for producing test cases.
    /// This simple implementation has some non-zero risk of collisions causing ambiguity in
    /// the sort. Risk should be low for small (less than 1M element) arrays.
    func shuffled(seed: Int) -> Array<Element> {
        let rand = randomArray(seed: seed, n: UInt(count))
        return zip(rand,self).sorted { $0.0 < $1.0 }.map { $0.1 }
    }
}
