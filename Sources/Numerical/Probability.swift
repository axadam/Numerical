//
//  Probability.swift
//  Numerical
//
//  Created by Adam Roberts on 10/6/19.
//

public struct Probability {
    private let value: Double
    private let isComplement: Bool
    
    public init(value: Double, isComplement: Bool) {
        self.value = value
        self.isComplement = isComplement
    }
}

public extension Probability {
    var p: Double { return isComplement ? 1 - value : value }
    var q: Double { return isComplement ? value : 1 - value }
    
    init(p: Double) {
        value = p
        isComplement = false
    }
    
    init(q: Double) {
        value = q
        isComplement = true
    }
    
    static func p(_ v: Double) -> Probability { return Probability(p: v) }
    static func q(_ v: Double) -> Probability { return Probability(q: v) }

    static let nan = Probability(value: .nan, isComplement: false)
    
    var pair: (Double, Double) {
        switch isComplement {
        case false: return (p, 1 - p)
        case true:  return (1 - q, q)
        }
    }
    
    static func -(lhs: Probability, rhs: Probability) -> Double {
        switch (lhs.p < 0.5, rhs.p < 0.5) {
        case (false,false): return rhs.q - lhs.q
        case (    _,    _): return lhs.p - rhs.p
        }
    }
}

extension Probability: CustomStringConvertible {
    public var description: String { return "p: \(p), q: \(q)" }
}

