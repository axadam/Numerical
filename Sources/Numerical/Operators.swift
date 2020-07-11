//
//  Operators.swift
//  Numerical
//
//  Created by Adam Roberts on 12/29/18.
//

precedencegroup ExponentiationPrecedence {
    associativity: right
    higherThan: MultiplicationPrecedence
}

infix operator ^^: ExponentiationPrecedence

public func ^^<T: FloatingPoint>(x: T, y: Int) -> T {
    if y == 0 { return 1 }
    let temp = x^^(y/2)
    switch (y,y % 2) {
    case (_,0):
        return temp * temp
    case (..<0,_):
        return temp * temp / x
    case (_,_):
        return temp * temp * x
    }
}
