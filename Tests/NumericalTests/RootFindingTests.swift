//
//  RootFindingTests.swift
//  NumericalTests
//
//  Created by Adam Roberts on 4/2/19.
//

import XCTest
import LogRelativeError
@testable import Numerical

final class RootFindingTests: XCTestCase {
    override class func tearDown() {
        print(resultStore.md())
        super.tearDown()
    }

    static let resultStore: ResultStore = ResultStore()
    var rs: ResultStore { return type(of: self).resultStore }

    let t = "Root Finding"
    let fl = "ln 5"
    let fc = "∛3647963"
    
    // Case 1: ln 5
    let f = { (x: Double) in exp(x) - 5 }
    let f1 = { (x: Double) in exp(x) }
    let f2f1 = { (x: Double) in 1.0 }
    let guessA = 2.0
    // ln(5)
    let c = "1.609437912434100374600759333226187639525601354268517721912"
    
    // Case 2: ∛3647963
    let g = { (x: Double) in x * x * x - 3647963.0 }
    let gʹ = { (x: Double) in 3 * x * x }
    let gʺ = { (x: Double) in 6 * x }
    let guessB = 364.0
    // cbrt(3647963)
    let d = "153.9395248014257325366856084389202348736844263890449292507"
    
    func testNewton() {
        let a = root(guess: guessA, f: f, f1: f1)
        AssertLRE(a,c,resultStore: rs, table: t, testCase: "Newton", field: fl)
        
        let b = root(guess: guessB, f: g, f1: gʹ)
        AssertLRE(b,d,resultStore: rs, table: t, testCase: "Newton", field: fc)
    }
    
    func testHalley() {
        let a = rootSecondOrder(guess: guessA, f: f, f1: f1, f2f1: f2f1)
        AssertLRE(a,c,resultStore: rs, table: t, testCase: "Halley", field: fl)
        
        let b = rootSecondOrder(guess: guessB, f: g, f1: gʹ, f2: gʺ)
        AssertLRE(b,d,resultStore: rs, table: t, testCase: "Halley", field: fc)
    }
    
    func testBracket() {
        let b = bracket(f: f, guess: 2)
        let fa = f(b!.a)
        let fb = f(b!.b)
        XCTAssert(fa * fb < 0)
    }
    
    func testBisection() {
        let a = root(guess: guessA, method: bisectionRoot, f: f)
        AssertLRE(a,c,digits: 14.9,resultStore: rs, table: t, testCase: "Bisection", field: fl)

        let b = root(guess: guessB, method: bisectionRoot, f: g)
        AssertLRE(b,d,digits: 14.5,resultStore: rs, table: t, testCase: "Bisection", field: fc)
    }
    
    func testSecant() {
        let a = root(guess: guessA, method: secantRoot, f: f)
        AssertLRE(a,c,resultStore: rs, table: t, testCase: "Secant", field: fl)
        
        let b = root(guess: guessB, method: secantRoot, f: g)
        AssertLRE(b,d,digits: 4.8,resultStore: rs, table: t, testCase: "Secant", field: fc)
    }
    
    func testDekker() {
        let a = root(guess: guessA, method: dekkerRoot, f: f)
        AssertLRE(a,c,resultStore: rs, table: t, testCase: "Dekker", field: fl)
        
        let b = root(guess: guessB, method: dekkerRoot, f: g)
        AssertLRE(b,d,resultStore: rs, table: t, testCase: "Dekker", field: fc)
    }
    
    func testRidder() {
        let a = root(guess: guessA, method: riddersRoot, f: f)
        AssertLRE(a,c,resultStore: rs, table: t, testCase: "Ridder", field: fl)
        
        let b = root(guess: guessB, method: riddersRoot, f: g)
        AssertLRE(b,d,resultStore: rs, table: t, testCase: "Ridder", field: fc)
    }
    
    func testBrent() {
        let a = root(guess: guessA, method: brentRoot, f: f)
        AssertLRE(a,c,resultStore: rs, table: t, testCase: "Brent", field: fl)
        
        let b = root(guess: guessB, method: brentRoot, f: g)
        AssertLRE(b,d,resultStore: rs, table: t, testCase: "Brent", field: fc)
    }
    
    func test748() {
        let a = root(guess: guessA, method: toms748Root, f: f)
        AssertLRE(a,c,resultStore: rs, table: t, testCase: "TOMS 748", field: fl)
        
        let b = root(guess: guessB, method: toms748Root, f: g)
        AssertLRE(b,d,resultStore: rs, table: t, testCase: "TOMS 748", field: fc)
    }
}
