//
//  ApproximateEquality.swift
//  Numerical
//
//  Created by Adam Roberts on 6/20/20.
//

import Foundation

/// A target value for approximate equality checking.
///
/// This value lets you specify if you know that the target value is zero. The  preferred method of
/// checking equality of floating point values is relative to the magnitude of the values, but that approach
/// does not work when comparing to zero.
///
/// Additional control is given by two optional parameters. One allows a scale to be specified when the
/// target is known to be zero. The other allows the target to be declared a trusted value so that allowable
/// error will be calculated relative to it.
public enum EqualityTarget<T: FloatingPoint> {
    /// The comparison target is zero
    ///
    /// If the scale of the quantities that led to zero is known that informs the internal tolerance calculation.
    ///
    /// This most often applies when you get zero by subtracting a calculation and a target value, such as
    /// `f(x) - c = 0`. Specify `c` for scale to define the relative tolerance. When possible it is better
    /// to rewrite this as `f(x) = c` and use the `maybeZero` target option.
    ///
    /// If `scaleRelativeTo` is zero or omitted the `absoluteForZero` tolerance is used.
    case zero(scaleRelativeTo: T = 0)
    
    /// The comparison target is a number not known to be zero.
    ///
    /// This will calculate error relative to the target value if `trusted` is set `true`, or to the larger of
    /// self and the target value by default.
    case maybeZero(_ x: T, trusted: Bool = false)
}

/// Tolerance for checking equality of floating point values
///
/// The tolerance value consists of a relative and an absolute tolerance for use when scale is known
/// and an additional aboslute tolerance that is used when scale is unknown.
public struct EqualityTolerance<T: FloatingPoint> {
    /// The allowable relative tolerance. In the hybrid allowable error equation `ε ﹤ rx + a` this is `r`
    public let relative: T
    
    /// The allowable absolute tolerance. In the hybrid allowable error equation `ε ﹤ rx + a` this is `a`
    public let absolute: T
    
    /// The allowable absolute tolerance in the case where we have no scale to define relative tolerance. 
    public let absoluteForZero: T
    
    public init(relative: T = T.ulpOfOne.squareRoot(), absolute: T = 0, absoluteForZero: T = T.ulpOfOne.squareRoot() ) {
        self.relative = relative
        self.absolute = absolute
        self.absoluteForZero = absoluteForZero
    }
}

public extension EqualityTolerance {
    static var standard: EqualityTolerance { EqualityTolerance() }
    static var strict: EqualityTolerance { EqualityTolerance(relative: 2 * T.ulpOfOne, absolute: 0, absoluteForZero: 8 * T.ulpOfOne) }
}

public extension FloatingPoint {
    func isApprox(_ other: EqualityTarget<Self>, tolerance: EqualityTolerance<Self> = .standard) -> Bool {
        switch other {
        case .zero(scaleRelativeTo: let s):
            let tol = s == 0 ? tolerance.absoluteForZero : tolerance.relative * s + tolerance.absolute
            return abs(self) < tol
        case .maybeZero(let x, trusted: let trust):
            let scale = trust ? abs(x) : max(abs(self),abs(x))
            return abs(self - x) < tolerance.relative * scale + tolerance.absolute
        }
    }
}

public extension BinaryFloatingPoint {
    func isWithinULP(of reference: Self, n: Int = 1024) -> Bool {
        switch (self, reference) {
        case (reference, _): return true
        case (..<Self(0),Self(0).nextUp), (Self(0).nextUp,..<Self(0)): return false
        case (_,_): return ulpDist(a: self, b: reference) < Self(n)
        }
    }
}

/// Distance between binary floating point values a and b expressed as the number of discrete floating
/// point values separating them
///
/// Floating point numbers have finite precision and as such represent discrete points. One way to express
/// the distance between two floating point numbers is as the number of these discrete points there are between
/// them.
///
/// The space between two adjacent floating point values is termed the Unit in the Last Place, or "ULP".
/// Binary floating point values within an interval bounded by powers of two, [2ⁿ,2ⁿ⁺¹), will all be separated
/// by the same size ULP. Such an interval is called a "binade" . Every binade has the same number of
/// discrete values, so there as the same number of binary floating point values between 1 and 2 as between
/// 1024 and 2048.
///
/// This distance method works by simply finding the distance between the two points and dividing by the size
/// of the ULP. If the interval [a,b] spans more than one binade then we are counting different sized ULPs and the
/// distance is pieced together. It is the distance from a to the top of it's binade, the distance of any complete
/// binades in between the top of a's binade and bottom of b's, and then the distance from the bottom of b's
/// binade to b.
public func ulpDist<T: BinaryFloatingPoint>(a: T, b: T) -> T {
    // if they are equal we are done
    if a == b { return 0 }

    // if they are opposite signs say they are far... we could return distance
    // from a to -0 and then 0 to b and sum them but this is rarely useful
    if a.sign != b.sign { return T.greatestFiniteMagnitude }
    
    // now we can assume same sign and if they are negative we can flip positive
    if a.sign == .minus { return ulpDist(a: -a, b: -b) }
    
    // make sure b is bigger than a
    if a > b { return ulpDist(a: b, b: a) }
    
    // the next exponent bigger than the smaller number
    let nextBinade: T = scalbn(T(1.0), Int(a.exponent + 1))
    
    // if b is not more than the next binade we can directly find the distance
    if b <= nextBinade { return (b - a) / a.ulp }
  
    // otherwise we are finding a distance that spans multiple binades and so
    // consists of multiple sizes of ulp.
    
    // a to next binade
    let aToNextBinade = (nextBinade - a) / a.ulp
    
    // distance between next binade and b's binade is the difference in
    // exponents times the number of values per exponent
    let upperToBinade = T(b.exponent - (a.exponent + 1)) * scalbn(T(1),T.significandBitCount)

    // b's binade to b
    let bsBinadeToB = (b - b.binade) / b.ulp
    
    // sum up the segments to get the total distance
    let total = aToNextBinade + upperToBinade + bsBinadeToB
    
    return floor(total)
}
