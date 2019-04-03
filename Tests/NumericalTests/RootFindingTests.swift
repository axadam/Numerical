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
    
    func testHalley() {
        let a = halley(guess: 2, f: f, f1: f1, f2f1: f2f1)
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
        let a = root(f: f, guess: 2, method: bisectionRoot)
        let b = f(a)
        XCTAssertEqual(b, 0, accuracy: 1e-10)
    }
    
    func testSecant() {
        let a = root(f: f, guess: 2, method: secantRoot)
        let b = f(a)
        XCTAssertEqual(b, 0, accuracy: 1e-10)
    }
    
    func testDekker() {
        let a = root(f: f, guess: 2, method: dekkerRoot)
        let b = f(a)
        XCTAssertEqual(b, 0, accuracy: 1e-10)
    }
    
    func testRidder() {
        let a = root(f: f, guess: 2, method: ridderRoot)
        let b = f(a)
        XCTAssertEqual(b, 0, accuracy: 1e-10)
    }
    
    func testBrent() {
        let a = root(f: f, guess: 2, method: brentRoot)
        let b = f(a)
        XCTAssertEqual(b, 0, accuracy: 1e-10)
    }
    
    //    brent 0.0
    //    brent 0.13567309812608463
    //    brent 0.13567309812608463
    //    brent 0.9240200047821446
    //    brent 2.0115564010052522
    //    brent 1.4684092311870962
    //    brent 1.5824402670778779
    //    brent 1.6099169644814735
    //    brent 1.60943141722178
    //    brent 1.6094379108784505
    //    brent 1.6094379124341005
    //    brent 1.6094379124341005
    
    func testBrentCPP() {
        let a = brentCPP(f: f, x1: 0, x2: 5)
        let b = f(a)
        XCTAssertEqual(b, 0, accuracy: 1e-10)
    }

}
