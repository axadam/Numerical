//
//  RootFindingTests.swift
//  NumericalTests
//
//  Created by Adam Roberts on 4/2/19.
//

import XCTest
@testable import Numerical

final class RootFindingTests: XCTestCase {
    let f = { (x: Double) in exp(x) - 5 }
    let f1 = { (x: Double) in exp(x) }
    let f2f1 = { (x: Double) in 1.0 }
    
    func testNewton() {
        let a = root(guess: 2, f: f, f1: f1)
        let b = f(a)
        XCTAssertEqual(b, 0, accuracy: 1e-10)
    }
    
    func testHalley() {
        let a = rootSecondOrder(guess: 2, f: f, f1: f1, f2f1: f2f1)
        let b = f(a)
        XCTAssertEqual(b, 0, accuracy: 1e-10)
    }
    
    func testBracket() {
        let b = bracket(f: f, guess: 2)
        let fa = f(b!.a)
        let fb = f(b!.b)
        XCTAssert(fa * fb < 0)
    }
    
    func testBisection() {
        let a = root(guess: 2, method: bisectionRoot, f: f)
        let b = f(a)
        XCTAssertEqual(b, 0, accuracy: 1e-10)
    }
    
    func testSecant() {
        let a = root(guess: 2, method: secantRoot, f: f)
        let b = f(a)
        XCTAssertEqual(b, 0, accuracy: 1e-10)
    }
    
    func testDekker() {
        let a = root(guess: 2, method: dekkerRoot, f: f)
        let b = f(a)
        XCTAssertEqual(b, 0, accuracy: 1e-10)
    }
    
    func testRidder() {
        let a = root(guess: 2, method: ridderRoot, f: f)
        let b = f(a)
        XCTAssertEqual(b, 0, accuracy: 1e-10)
    }
    
    func testBrent() {
        let a = root(guess: 2, method: brentRoot, f: f)
        let b = f(a)
        XCTAssertEqual(b, 0, accuracy: 1e-10)
    }
    
}
