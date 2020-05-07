//
//  QuadratureTests.swift
//  NumericalTests
//
//  Created by Adam Roberts on 4/24/20.
//

import Foundation


import XCTest
import LogRelativeError
@testable import Numerical

final class Quadrature: XCTestCase {
    override class func tearDown() {
        print(resultStore.md())
        super.tearDown()
    }

    static let resultStore: ResultStore = ResultStore()
    var rs: ResultStore { return type(of: self).resultStore }

    let cube = { (x: Double) in x * x * x }
    let recip = { (x: Double) in 1 / x }
    let iden = { (x: Double) in x }
    // Poisson's example: perimeter of ellipse of 1/π by 0.6/π
    let pois = { (x: Double) in sqrt(1 - 0.36 * sin(x) * sin(x)) / (2 * Double.pi) }
    let ecos = { (x: Double) in exp(cos(x)) }
    let ex2  = { (x: Double) in exp(-x * x) / Double.pi.squareRoot() }
    let pi   = { (t: Double) in 1 / (1 + t * t) }
    
    let t = "Quadrature"
    
    let tcube = "x³; [0,1]"
    let trecip = "1/x; [1,100]"
    let tiden = "x; [0,5000]"
    let tpois = "√(1 - 0.36sin²θ) / √(2π); [0,2π]"
    let tecos = "e^cos(θ); [0,2π]"
    let tex2  = "e^(-x²) / √π; [-10,10]" // truncating [-∞,∞] which would be 1
    let tpi   = "π = 4 ∫0...1 1 / (1 + t²) dt; [0,1]"
    
    let vcube = "0.25" // exact
    let vrecip = "4.605170185988091368035982909368728415202202977257545952066"
    let viden = "12500000" // exact
    let vpois = "0.902779927772193884716139509461678950653101811833151562088"
    let vecos = "7.95492652101284"
    let vex2  = "0.999999999999999999999999999999999999999999997911512416237"
    let vpi   = "3.141592653589793238462643383279502884197169399375105820974"
    
    func testTrapezoidal() {
        let f = "Trapezoidal"
        
        let c = integrate(range: 0...1, method: trapezoidal, f: cube)
        let r = integrate(range: 1...100, maxIter: 20, method: trapezoidal, f: recip)
        let i = integrate(range: 0...5000, method: trapezoidal, f: iden)
        let p = integrate(range: 0...(2 * Double.pi), method: trapezoidal, f: pois)
        let e = integrate(range: 0...(2 * Double.pi), method: trapezoidal, f: ecos)
        let x = integrate(range: -10...10, method: trapezoidal, f: ex2)
        let π4 = integrate(range: 0...1, method: trapezoidal, f: pi)
        
        AssertLRE(c.value, vcube, exact: true, digits: 6.0, resultStore: rs, table: t, testCase: tcube, field: f, annotation: "\(c.evals)\(c.converged ? "" : "*")")
        AssertLRE(r.value, vrecip, digits: 9.7, resultStore: rs, table: t, testCase: trecip, field: f, annotation: "\(r.evals)\(r.converged ? "" : "*")")
        AssertLRE(i.value, viden, exact: true, resultStore: rs, table: t, testCase: tiden, field: f, annotation: "\(i.evals)\(i.converged ? "" : "*")")
        AssertLRE(p.value, vpois, resultStore: rs, table: t, testCase: tpois, field: f, annotation: "\(p.evals)\(p.converged ? "" : "*")")
        AssertLRE(e.value, vecos, resultStore: rs, table: t, testCase: tecos, field: f, annotation: "\(e.evals)\(e.converged ? "" : "*")")
        AssertLRE(x.value, vex2, resultStore: rs, table: t, testCase: tex2, field: f, annotation: "\(x.evals)\(x.converged ? "" : "*")")
        AssertLRE(4 * π4.value, vpi, digits: 7.2, resultStore: rs, table: t, testCase: tpi, field: f, annotation: "\(π4.evals)\(π4.converged ? "" : "*")")
    }
    
    func testRomberg() {
        let f = "Romberg"
        
        let c = integrate(range: 0...1, method: romberg, f: cube)
        let r = integrate(range: 1...100, maxIter: 20, method: romberg, f: recip)
        let i = integrate(range: 0...5000, method: romberg, f: iden)
        let p = integrate(range: 0...(2 * Double.pi), method: romberg, f: pois)
        let e = integrate(range: 0...(2 * Double.pi), method: romberg, f: ecos)
        let x = integrate(range: -10...10, method: romberg, f: ex2)
        let π4 = integrate(range: 0...1, method: romberg, f: pi)

        AssertLRE(c.value, vcube, exact: true, resultStore: rs, table: t, testCase: tcube, field: f, annotation: "\(c.evals)\(c.converged ? "" : "*")")
        AssertLRE(r.value, vrecip, digits: 14.9, resultStore: rs, table: t, testCase: trecip, field: f, annotation: "\(r.evals)\(r.converged ? "" : "*")")
        AssertLRE(i.value, viden, exact: true, resultStore: rs, table: t, testCase: tiden, field: f, annotation: "\(i.evals)\(i.converged ? "" : "*")")
        AssertLRE(p.value, vpois, resultStore: rs, table: t, testCase: tpois, field: f, annotation: "\(p.evals)\(p.converged ? "" : "*")")
        AssertLRE(e.value, vecos, resultStore: rs, table: t, testCase: tecos, field: f, annotation: "\(e.evals)\(e.converged ? "" : "*")")
        AssertLRE(x.value, vex2, digits: 14.6, resultStore: rs, table: t, testCase: tex2, field: f, annotation: "\(x.evals)\(x.converged ? "" : "*")")
        AssertLRE(4 * π4.value, vpi, resultStore: rs, table: t, testCase: tpi, field: f, annotation: "\(π4.evals)\(π4.converged ? "" : "*")")
    }
}
