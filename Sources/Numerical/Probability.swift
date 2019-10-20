//
//  Probability.swift
//  Numerical
//
//  Created by Adam Roberts on 10/6/19.
//

/// A value representing a probability between zero and one.
///
/// The underlying value may be stored as either the probability itself or
/// as its complement. This allows us to store numbers either very close
/// to zero or very close to one without loss of precision.
///
/// This value can be initialized based on either the probability or the
/// complement, and both are available as properties.
public struct Probability {
    
    /// The underlying value used when initializing.
    /// May be either the probability or its complement.
    private let value: Double
    
    /// Whether the underlying value is the probability of the complement
    private let isComplement: Bool
    
    public init(value: Double, isComplement: Bool) {
        self.value = value
        self.isComplement = isComplement
    }
}

public extension Probability {
    /// The probability, p
    var p: Double { return isComplement ? 1 - value : value }
    
    /// The complementary probability, q
    var q: Double { return isComplement ? value : 1 - value }
    
    init(p: Double) {
        self.init(value: p, isComplement: false)
    }
    
    init(q: Double) {
        self.init(value: q, isComplement: true)
    }
    
    static func p(_ v: Double) -> Probability { return Probability(p: v) }
    static func q(_ v: Double) -> Probability { return Probability(q: v) }

    static let nan = Probability(value: .nan, isComplement: false)
    
    /// Difference between this and another probability.
    ///
    /// Result not defined to be a probability because it may not be between zero and one.
    func difference(_ rhs: Probability) -> Double {
        switch (p < 0.5, rhs.p < 0.5) {
        case (false,false): return rhs.q - q
        case (    _,    _): return p - rhs.p
        }
    }
}

extension Probability: Codable, Hashable, Equatable {}

extension Probability: Comparable {
    public static func < (lhs: Probability, rhs: Probability) -> Bool {
        switch (lhs.p < 0.5, rhs.p < 0.5) {
        case (false,false): return rhs.q < lhs.q
        case (    _,    _): return lhs.p < rhs.p
        }
    }
}

extension Probability: CustomStringConvertible {
    public var description: String { return "p: \(p), q: \(q)" }
}

