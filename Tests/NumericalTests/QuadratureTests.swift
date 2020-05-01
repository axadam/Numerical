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
    
    let t = "Quadrature"
    
    let fcube = "x³, [0,1]"
    let frecip = "1/x, [1,100]"
    let fiden = "x, [0,5000]"
    let fpois = "√(1 - 0.36sin²θ) / √(2π), [0,2π]"
    let fecos = "e^cos(θ), [0,2π]"
    let fex2  = "e^(-x²) / √π, [-10,10]" // truncating [-∞,∞] which would be 1
    
    let vcube = "0.25" // exact
    let vrecip = "4.605170185988091368035982909368728415202202977257545952066"
    let viden = "12500000" // exact
    let vpois = "0.902779927772193884716139509461678950653101811833151562088"
    let vecos = "7.95492652101284"
    let vex2  = "0.999999999999999999999999999999999999999999997911512416237"
    
    func testTrapezoidal() {
        let tc = "Trapezoidal"
        
        let c = integrate(range: 0...1, maxIter: 100, method: trapezoidal, f: cube)
        let r = integrate(range: 1...100, maxIter: 1000, method: trapezoidal, f: recip)
        let i = integrate(range: 0...5000, method: trapezoidal, f: iden)
        let p = integrate(range: 0...(2 * Double.pi), method: trapezoidal, f: pois)
        let e = integrate(range: 0...(2 * Double.pi), method: trapezoidal, f: ecos)
        let x = integrate(range: -10...10, method: trapezoidal, f: ex2)
        
        AssertLRE(c, vcube, exact: true, digits: 10.8, resultStore: rs, table: t, testCase: tc, field: fcube)
        AssertLRE(r, vrecip, digits: 10.9, resultStore: rs, table: t, testCase: tc, field: frecip)
        AssertLRE(i, viden, exact: true, resultStore: rs, table: t, testCase: tc, field: fiden)
        AssertLRE(p, vpois, resultStore: rs, table: t, testCase: tc, field: fpois)
        AssertLRE(e, vecos, resultStore: rs, table: t, testCase: tc, field: fecos)
        AssertLRE(x, vex2, resultStore: rs, table: t, testCase: tc, field: fex2)
    }
    
    func testRomberg() {
        let tc = "Romberg"
        
        let c = integrate(range: 0...1, maxIter: 100, method: romberg, f: cube)
        let r = integrate(range: 1...100, maxIter: 1000, method: romberg, f: recip)
        let i = integrate(range: 0...5000, method: romberg, f: iden)
        let p = integrate(range: 0...(2 * Double.pi), method: romberg, f: pois)
        let e = integrate(range: 0...(2 * Double.pi), method: romberg, f: ecos)
        let x = integrate(range: -10...10, method: romberg, f: ex2)
        
        AssertLRE(c, vcube, exact: true, resultStore: rs, table: t, testCase: tc, field: fcube)
        AssertLRE(r, vrecip, digits: 14.9, resultStore: rs, table: t, testCase: tc, field: frecip)
        AssertLRE(i, viden, exact: true, resultStore: rs, table: t, testCase: tc, field: fiden)
        AssertLRE(p, vpois, resultStore: rs, table: t, testCase: tc, field: fpois)
        AssertLRE(e, vecos, resultStore: rs, table: t, testCase: tc, field: fecos)
        AssertLRE(x, vex2, digits: 14.6, resultStore: rs, table: t, testCase: tc, field: fex2)
    }
}
