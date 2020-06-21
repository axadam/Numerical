//
//  RootFindingTests.swift
//  NumericalTests
//
//  Created by Adam Roberts on 4/2/19.
//

import XCTest
import LogRelativeError
@testable import Numerical

struct RootTC {
    let name: String
    let f: (Double) -> Double
    let fʹ: (Double) -> Double
    let fʺ: (Double) -> Double
    let fʺfʹ: (Double) -> Double
    let intercept: Double
    let guess: Double
    let answer: String
}

final class RootFindingTests: XCTestCase {
    override class func tearDown() {
        print(resultStore.md())
        super.tearDown()
    }

    static let resultStore: ResultStore = ResultStore()
    var rs: ResultStore { return type(of: self).resultStore }

    let t = "Root Finding"
    
    /// Calculate the natural logarithm as the inverse of the exponential function
    let ln5 = RootTC(
        name: "ln 5; x₀=2"
        , f: { (x: Double) in exp(x) }
        , fʹ: { (x: Double) in exp(x) }
        , fʺ: { (x: Double) in exp(x) }
        , fʺfʹ: { (x: Double) in 1.0 }
        , intercept: 5
        , guess: 2.0
        , answer: "1.609437912434100374600759333226187639525601354268517721912")

    /// Calculate the cube root as the inverse of x³
    let bigCube = RootTC(
        name: "∛3647963; x₀=364"
        , f: { (x: Double) in x * x * x }
        , fʹ: { (x: Double) in 3 * x * x }
        , fʺ: { (x: Double) in 6 * x }
        , fʺfʹ: { (x: Double) in 2 / x }
        , intercept: 3647963.0
        , guess: 364.0
        , answer: "153.9395248014257325366856084389202348736844263890449292507")

    /// Newton's first test case for his method, a cubic equation
    let newt = RootTC(
        name: "x³ - 2x - 5; x₀=2"
        , f: { (x: Double) in x * x * x - 2 * x - 5 }
        , fʹ: { (x: Double) in 3 * x * x - 2 }
        , fʺ: { (x: Double) in 6 * x }
        , fʺfʹ: { (x: Double) in 6 * x / (3 * x * x - 2) }
        , intercept: 0.0
        , guess: 2.0
        , answer: "2.094551481542326591482386540579303")

    /// TOMS Algorithm 748 (1995), Table I, Case 1
    let toms1 = RootTC(
        name: "sin(x) - x/2; x₀=2"
        , f: { (x: Double) in sin(x) - x / 2 }
        , fʹ: { (x: Double) in cos(x) - 1 / 2 }
        , fʺ: { (x: Double) in -sin(x) }
        , fʺfʹ: { (x: Double) in -sin(x) / (cos(x) - 1 / 2) }
        , intercept: 0.0
        , guess: 2.0
        , answer: "1.895494267033980947144035738093601691751")

    func testNewton() {
        let tc = "Newton"
        
        let a = root(guess: ln5.guess, f: { self.ln5.f($0) - self.ln5.intercept }, f1: ln5.fʹ)
        AssertLRE(a.value,ln5.answer,resultStore: rs, table: t, testCase: tc, field: ln5.name, annotation: a.note)
        
        let b = root(guess: bigCube.guess, f: { self.bigCube.f($0) - self.bigCube.intercept }, f1: bigCube.fʹ)
        AssertLRE(b.value,bigCube.answer,resultStore: rs, table: t, testCase: tc, field: bigCube.name, annotation: b.note)
        
        let c = root(guess: newt.guess, f: newt.f, f1: newt.fʹ)
        AssertLRE(c.value,newt.answer,resultStore: rs, table: t, testCase: tc, field: newt.name, annotation: c.note)
        
        let d = root(guess: toms1.guess, f: toms1.f, f1: toms1.fʹ)
        AssertLRE(d.value,toms1.answer,resultStore: rs, table: t, testCase: tc, field: toms1.name, annotation: d.note)
    }
    
    func testHalley() {
        let tc = "Halley"
        
        let a = root(guess: ln5.guess, f: { self.ln5.f($0) - self.ln5.intercept }, f1: ln5.fʹ, f2f1: ln5.fʺfʹ)
        AssertLRE(a.value,ln5.answer,resultStore: rs, table: t, testCase: tc, field: ln5.name, annotation: a.note)

        let b = root(guess: bigCube.guess, f: { self.bigCube.f($0) - self.bigCube.intercept }, f1: bigCube.fʹ, f2: bigCube.fʺ)
        AssertLRE(b.value,bigCube.answer,resultStore: rs, table: t, testCase: tc, field: bigCube.name, annotation: b.note)

        let c = root(guess: newt.guess, f: newt.f, f1: newt.fʹ, f2: newt.fʺ)
        AssertLRE(c.value,newt.answer,resultStore: rs, table: t, testCase: tc, field: newt.name, annotation: c.note)

        let d = root(guess: toms1.guess, f: toms1.f, f1: toms1.fʹ, f2: toms1.fʺ)
        AssertLRE(d.value,toms1.answer,resultStore: rs, table: t, testCase: tc, field: toms1.name, annotation: d.note)
    }
    
    func testBracket() {
        let g = { self.ln5.f($0) - self.ln5.intercept }
        let b = bracket(f: g , guess: 2)
        guard case let .bracket(_,e) = b else { XCTFail(); return }
        let fa = g(e.a)
        let fb = g(e.b)
        XCTAssert(fa * fb < 0)
    }
    
    func testBisection() {
        let tc = "Bisection"
        let method = bisectionRoot
        
        let a = root(guess: ln5.guess, method: method, intercept: ln5.intercept, f: ln5.f)
        AssertLRE(a.value,ln5.answer,digits: 14.9,resultStore: rs, table: t, testCase: tc, field: ln5.name, annotation: a.note)

        let b = root(guess: bigCube.guess, method: method, intercept: bigCube.intercept, f: bigCube.f)
        AssertLRE(b.value,bigCube.answer,resultStore: rs, table: t, testCase: tc, field: bigCube.name, annotation: b.note)

        let c = root(guess: newt.guess, method: method, f: newt.f)
        AssertLRE(c.value,newt.answer,resultStore: rs, table: t, testCase: tc, field: newt.name, annotation: c.note)

        let d = root(guess: toms1.guess, method: method, f: toms1.f)
        AssertLRE(d.value,toms1.answer,resultStore: rs, table: t, testCase: tc, field: toms1.name, annotation: d.note)
    }
    
    func testSecant() {
        let tc = "Secant"
        let method = secantRoot
        
        let a = root(guess: ln5.guess, method: method, intercept: ln5.intercept, f: ln5.f)
        AssertLRE(a.value,ln5.answer,resultStore: rs, table: t, testCase: tc, field: ln5.name, annotation: a.note)

        let b = root(guess: bigCube.guess, method: method, intercept: bigCube.intercept, f: bigCube.f)
        AssertLRE(b.value,bigCube.answer,digits: 7.8,resultStore: rs, table: t, testCase: tc, field: bigCube.name, annotation: b.note)

        let c = root(guess: newt.guess, method: method, f: newt.f)
        AssertLRE(c.value,newt.answer,resultStore: rs, table: t, testCase: tc, field: newt.name, annotation: c.note)

        let d = root(guess: toms1.guess, method: method, f: toms1.f)
        AssertLRE(d.value,toms1.answer,digits: 0.0, resultStore: rs, table: t, testCase: tc, field: toms1.name, annotation: d.note + "†")
    }
    
    func testDekker() {
        let tc = "Dekker"
        let method = dekkerRoot
        
        let a = root(guess: ln5.guess, method: method, intercept: ln5.intercept, f: ln5.f)
        AssertLRE(a.value,ln5.answer,resultStore: rs, table: t, testCase: tc, field: ln5.name, annotation: a.note)

        let b = root(guess: bigCube.guess, method: method, intercept: bigCube.intercept, f: bigCube.f)
        AssertLRE(b.value,bigCube.answer,resultStore: rs, table: t, testCase: tc, field: bigCube.name, annotation: b.note)

        let c = root(guess: newt.guess, method: method, f: newt.f)
        AssertLRE(c.value,newt.answer,resultStore: rs, table: t, testCase: tc, field: newt.name, annotation: c.note)

        let d = root(guess: toms1.guess, method: method, f: toms1.f)
        AssertLRE(d.value,toms1.answer,resultStore: rs, table: t, testCase: tc, field: toms1.name, annotation: d.note)
    }
    
    func testRidders() {
        let tc = "Ridders"
        let method = riddersRoot
        
        let a = root(guess: ln5.guess, method: method, intercept: ln5.intercept, f: ln5.f)
        AssertLRE(a.value,ln5.answer,resultStore: rs, table: t, testCase: tc, field: ln5.name, annotation: a.note)

        let b = root(guess: bigCube.guess, method: method, intercept: bigCube.intercept, f: bigCube.f)
        AssertLRE(b.value,bigCube.answer,resultStore: rs, table: t, testCase: tc, field: bigCube.name, annotation: b.note)

        let c = root(guess: newt.guess, method: method, f: newt.f)
        AssertLRE(c.value,newt.answer,resultStore: rs, table: t, testCase: tc, field: newt.name, annotation: c.note)

        let d = root(guess: toms1.guess, method: method, f: toms1.f)
        AssertLRE(d.value,toms1.answer,resultStore: rs, table: t, testCase: tc, field: toms1.name, annotation: d.note)
    }
    
    func testBrent() {
        let tc = "Brent"
        let method = brentRoot
        
        let a = root(guess: ln5.guess, method: method, intercept: ln5.intercept, f: ln5.f)
        AssertLRE(a.value,ln5.answer,resultStore: rs, table: t, testCase: tc, field: ln5.name, annotation: a.note)

        let b = root(guess: bigCube.guess, method: method, intercept: bigCube.intercept, f: bigCube.f)
        AssertLRE(b.value,bigCube.answer,resultStore: rs, table: t, testCase: tc, field: bigCube.name, annotation: b.note)

        let c = root(guess: newt.guess, method: method, f: newt.f)
        AssertLRE(c.value,newt.answer,resultStore: rs, table: t, testCase: tc, field: newt.name, annotation: c.note)

        let d = root(guess: toms1.guess, method: method, f: toms1.f)
        AssertLRE(d.value,toms1.answer,resultStore: rs, table: t, testCase: tc, field: toms1.name, annotation: d.note)
    }
    
    func test748() {
        let tc = "TOMS 748"
        let method = toms748Root
        
        let a = root(guess: ln5.guess, method: method, intercept: ln5.intercept, f: ln5.f)
        AssertLRE(a.value,ln5.answer,resultStore: rs, table: t, testCase: tc, field: ln5.name, annotation: a.note)

        let b = root(guess: bigCube.guess, method: method, intercept: bigCube.intercept, f: bigCube.f)
        AssertLRE(b.value,bigCube.answer,resultStore: rs, table: t, testCase: tc, field: bigCube.name, annotation: b.note)

        let c = root(guess: newt.guess, method: method, f: newt.f)
        AssertLRE(c.value,newt.answer,resultStore: rs, table: t, testCase: tc, field: newt.name, annotation: c.note)

        let d = root(guess: toms1.guess, method: method, f: toms1.f)
        AssertLRE(d.value,toms1.answer,resultStore: rs, table: t, testCase: tc, field: toms1.name, annotation: d.note)
    }
}

fileprivate extension BracketAndRootResult {
    var note: String {
        return "\(bracketEvals)/\(rootEvals)\(converged ? "" : "*")"
    }
}

fileprivate extension RootResult {
    var note: String {
        return "\(evals)\(converged ? "" : "*")"
    }
}
