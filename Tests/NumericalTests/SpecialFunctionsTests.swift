//
//  SpecialFunctionsTests.swift
//  NumericalTests
//
//  Created by Adam Roberts on 3/20/19.
//

import XCTest
import LogRelativeError
@testable import Numerical

final class SpecialFunctionsTests: XCTestCase {
    override class func tearDown() {
        print(resultStore.md())
        super.tearDown()
    }

    static let resultStore: ResultStore = ResultStore()
    var rs: ResultStore { return type(of: self).resultStore }

    let f = "f"
    let f1 = "f⁻¹"
    let te = "Error Function"
    let tg = "Regularized Incomplete Gamma Function, 'GammaCHI', Gil, Tegura, and Temme 2015, Table 1"
    let tg1 = "Regularized Incomplete Gamma Function, reference values from Mathematica"
    let tb = "Regularized Incomplete Beta Function"
    let tm = "Marcum Q Function, reference values from Mathematica"
    
    func testInvErf() {
        let a = erf(1.2345)
        AssertLRE(a, "0.9191623964135658", resultStore: rs, table: te, testCase: "1.2345", field: f)
        let b = invErf(a)
        AssertLRE(b, "1.2345", exact: true, resultStore: rs, table: te, testCase: "1.2345", field: f1)
    }
    
    func testIncompleteGamma() {
        let a = q_gamma(4, 0.7)
        // Q(4,7) = 0.9942465424077004166914...
        AssertLRE(a, "0.9942465424077004166914", resultStore: rs, table: tg1, testCase: "(a:4,x:0.7)", field: f)

        // "GammaCHI: a package for the inversion and computation of the gamma
        // and chi-square cumulative distribution functions (central and noncentral)",
        // Gil, Tegura, and Temme 2015, Table 1
        AssertLRE(q_gamma(1e-249, 6.310e-15), "3.212101109661167e-248", digits: 4.3, resultStore: rs, table: tg, testCase: "(a:1e-249,x:6.310e-15)", field: f)
        AssertLRE(q_gamma(1e-249, 7.110e-7 ), "1.3580785912009393e-248", digits: 3.9, resultStore: rs, table: tg, testCase: "(a:1e-249,x:7.110e-7)", field: f)
        AssertLRE(q_gamma(1e-249, 0.01     ), "4.0379295765381135e-249", resultStore: rs, table: tg, testCase: "(a:1e-249,x:0.01)", field: f)
        AssertLRE(q_gamma(1e-13, 6.310e-15 ), "3.212101109660651e-12", digits: 4.3, resultStore: rs, table: tg, testCase: "(a:1e-13,x:6.310e-15)", field: f)
        AssertLRE(q_gamma(1e-13, 7.110e-7  ), "1.358078591200848e-12", digits: 3.9, resultStore: rs, table: tg, testCase: "(a:1e-13,x:7.110e-7)", field: f)
        AssertLRE(q_gamma(1e-13, 0.01      ), "4.0379295765380405e-13", digits: 12.7, resultStore: rs, table: tg, testCase: "(a:1e-13,x:0.01)", field: f)
    }
    
    func testInverseIncompleteGamma() {
        let a = inv_p_gamma(4, 1 - 0.9942465424077004166914)
        AssertLRE(a, "0.7", exact: true, digits: 14.7, resultStore: rs, table: tg1, testCase: "(a:4,x:0.7)", field: f1)
        
        AssertLRE(inv_q_gamma(1e-249, 3.212101109661167e-248), "6.310e-15", exact: true, digits: 0.0, resultStore: rs, table: tg, testCase: "(a:1e-249,x:6.310e-15)", field: f1)
        AssertLRE(inv_q_gamma(1e-249, 1.3580785912009393e-248), "7.110e-7", exact: true, digits: 0.0, resultStore: rs, table: tg, testCase: "(a:1e-249,x:7.110e-7)", field: f1)
        AssertLRE(inv_q_gamma(1e-249, 4.0379295765381135e-249), "0.01", exact: true, digits: 0.0, resultStore: rs, table: tg, testCase: "(a:1e-249,x:0.01)", field: f1)
        AssertLRE(inv_q_gamma(1e-13, 3.212101109660651e-12), "6.310e-15", exact: true, digits: 2.8, resultStore: rs, table: tg, testCase: "(a:1e-13,x:6.310e-15)", field: f1)
        AssertLRE(inv_q_gamma(1e-13, 1.358078591200848e-12), "7.110e-7", exact: true, digits: 2.8, resultStore: rs, table: tg, testCase: "(a:1e-13,x:7.110e-7)", field: f1)
        AssertLRE(inv_q_gamma(1e-13, 4.0379295765380405e-13), "0.01", exact: true, digits: 12.1, resultStore: rs, table: tg, testCase: "(a:1e-13,x:0.01)", field: f1)
    }
    
    func testIncompleteBeta() {
        let a = beta_reg(x: 0.4, a: 3, b: 5)
        // I₀ٖ₄(3,5) = 0.580096, exact
        AssertLRE(a, "0.580096", exact: true, digits: 14.8, resultStore: rs, table: tb, testCase: "(a:3,b:5,x:0.4)", field: f)
    }
    
    func testInverseIncompleteBeta() {
        let a = inv_beta_reg(p: 0.580096, a: 3, b: 5)
        AssertLRE(a, "0.4", exact: true, resultStore: rs, table: tb, testCase: "(a:3,b:5,x:0.4)", field: f1)
    }

    /// Marcum Q asymptotic large µ test cases
    ///
    /// Ten tests of the large µ case for the Marcum Q function. The varied dimension
    /// is the non-centrality parameter.
    ///
    /// "Computation of the Marcum Q-function", Gil, Segura, Temme 2013, Table 6.1
    func testMarcumBigMu() {
        let t = "Marcum Q large µ. µ = 8192, y = 1.05µ, and x is a fraction of µ. 'Computation of the Marcum Q-function', Gil, Segura, Temme 2013, Table 6.1"
        let µ: Double = 8192
        let y = 1.05 * µ
        func x(_ x: Double) -> Double { return x * 8192 }
        AssertLRE(marcum(µ: µ, x: x(0.01), y: y).q, "1.9845278031193e-4", digits: 7.0, resultStore: rs, table: t, testCase: "x:0.01µ", field: f)
        AssertLRE(marcum(µ: µ, x: x(0.02), y: y).q, "4.138241872117e-3", digits: 7.9, resultStore: rs, table: t, testCase: "x:0.02µ", field: f)
        AssertLRE(marcum(µ: µ, x: x(0.03), y: y).q, "0.04000364971081", digits: 9.1, resultStore: rs, table: t, testCase: "x:0.03µ", field: f)
        AssertLRE(marcum(µ: µ, x: x(0.04), y: y).q, "0.191650654805848", digits: 9.6, resultStore: rs, table: t, testCase: "x:0.04µ", field: f)
        AssertLRE(marcum(µ: µ, x: x(0.05), y: y).q, "0.498535453743169", digits: 9.9, resultStore: rs, table: t, testCase: "x:0.05µ", field: f)
        AssertLRE(marcum(µ: µ, x: x(0.06), y: y).p, "0.1964796269915073", digits: 9.5, resultStore: rs, table: t, testCase: "x:0.06µ", field: f)
        AssertLRE(marcum(µ: µ, x: x(0.07), y: y).p, "0.04434265824612003", digits: 9.2, resultStore: rs, table: t, testCase: "x:0.07µ", field: f)
        AssertLRE(marcum(µ: µ, x: x(0.08), y: y).p, "5.526239087335513e-3", digits: 8.0, resultStore: rs, table: t, testCase: "x:0.08µ", field: f)
        AssertLRE(marcum(µ: µ, x: x(0.09), y: y).p, "3.750276163593746e-4", digits: 7.1, resultStore: rs, table: t, testCase: "x:0.09µ", field: f)
        AssertLRE(marcum(µ: µ, x: x(0.10), y: y).p, "1.386276448162126e-5", digits: 6.5, resultStore: rs, table: t, testCase: "x:0.10µ", field: f)

        let ystr = "\(y)"
        AssertLRE(inv_marcum(µ: µ, x: x(0.01), p: .q(1.9845278031193e-4)), ystr, exact: true, digits: 9.5, resultStore: rs, table: t, testCase: "x:0.01µ", field: f1)
        AssertLRE(inv_marcum(µ: µ, x: x(0.02), p: .q(4.138241872117e-3)), ystr, exact: true, digits: 10.3, resultStore: rs, table: t, testCase: "x:0.02µ", field: f1)
        AssertLRE(inv_marcum(µ: µ, x: x(0.03), p: .q(0.04000364971081)), ystr, exact: true, digits: 11.4, resultStore: rs, table: t, testCase: "x:0.03µ", field: f1)
        AssertLRE(inv_marcum(µ: µ, x: x(0.04), p: .q(0.191650654805848)), ystr, exact: true, digits: 11.7, resultStore: rs, table: t, testCase: "x:0.04µ", field: f1)
        AssertLRE(inv_marcum(µ: µ, x: x(0.05), p: .q(0.498535453743169)), ystr, exact: true, digits: 11.8, resultStore: rs, table: t, testCase: "x:0.05µ", field: f1)
        AssertLRE(inv_marcum(µ: µ, x: x(0.06), p: .p(0.1964796269915073)), ystr, exact: true, digits: 11.6, resultStore: rs, table: t, testCase: "x:0.06µ", field: f1)
        AssertLRE(inv_marcum(µ: µ, x: x(0.07), p: .p(0.04434265824612003)), ystr, exact: true, digits: 11.5, resultStore: rs, table: t, testCase: "x:0.07µ", field: f1)
        AssertLRE(inv_marcum(µ: µ, x: x(0.08), p: .p(5.526239087335513e-3)), ystr, exact: true, digits: 10.4, resultStore: rs, table: t, testCase: "x:0.08µ", field: f1)
        AssertLRE(inv_marcum(µ: µ, x: x(0.09), p: .p(3.750276163593746e-4)), ystr, exact: true, digits: 9.7, resultStore: rs, table: t, testCase: "x:0.09µ", field: f1)
        AssertLRE(inv_marcum(µ: µ, x: x(0.10), p: .p(1.386276448162126e-5)), ystr, exact: true, digits: 9.1, resultStore: rs, table: t, testCase: "x:0.10µ", field: f1)
    }
    
    /// Marcum Q far tail test cases
    ///
    /// Four tests on the far left tail. The first three use the quadrature method and
    /// the last uses the asymptotic expansion for large ξ.
    ///
    /// "GammaCHI: a package for the inversion and computation of the gamma
    /// and chi-square cumulative distribution functions (central and noncentral)",
    /// Gil, Tegura, and Temme 2015, Table 3
    func testMarcumExtremeTail() {
        let t = "Marcum Q far left tail, 'GammaCHI', Gil, Tegura, and Temme 2015, Table 3"
        AssertLRE(marcum(µ: 5 , x: 150, y: 30 ).p, "1.215915354045e-23", resultStore: rs, table: t, testCase: "µ:5,x:150,y:30", field: f)
        AssertLRE(marcum(µ: 1 , x: 75,  y: 0.5).p, "3.287840255874e-30", resultStore: rs, table: t, testCase: "µ:1,x:75,y:0.5", field: f)
        AssertLRE(marcum(µ: 2 , x: 100, y: 2  ).p, "1.557081489535e-35", digits: 12.3, resultStore: rs, table: t, testCase: "µ:2,x:100,y:2", field: f)
        AssertLRE(marcum(µ: 10, x: 100, y: 1  ).p, "5.152185145235e-48", resultStore: rs, table: t, testCase: "µ:10,x:100,y:1", field: f)

        AssertLRE(inv_marcum(µ: 5, x: 150, p: .p(1.215915354045e-23)), "30.00000000000", resultStore: rs, table: t, testCase: "µ:5,x:150,y:30", field: f1)
        AssertLRE(inv_marcum(µ: 1, x: 75, p: .p(3.287840255874e-30)), "0.5000000000000", resultStore: rs, table: t, testCase: "µ:1,x:75,y:0.5", field: f1)
        AssertLRE(inv_marcum(µ: 2, x: 100, p: .p(1.557081489535e-35)), "2.000000000000", resultStore: rs, table: t, testCase: "µ:2,x:100,y:2", field: f1)
        AssertLRE(inv_marcum(µ: 10, x: 100, p: .p(5.152185145235e-48)), "1.000000000000", resultStore: rs, table: t, testCase: "µ:10,x:100,y:1", field: f1)

    }
    
    /// Marcum Q test cases for various methods
    ///
    /// Method used by each test case is noted in comments. Reference values from
    /// Mathematica using the syntax: N[MarcumQ[µ,sqrt(2 * x),0,sqrt(2 * y)]]
    func testMarcum() {
        // p series
        AssertLRE(marcum(µ: 11.5, x: 15.3, y: 23).p, "0.2948691834572695", digits: 13.3, resultStore: rs, table: tm, testCase: "µ:11.5,x:15.3,y:23", field: f)
        AssertLRE(inv_marcum(µ: 11.5, x: 15.3, p: .p(0.2948691834572695)), "23", exact: true, digits: 13.9, resultStore: rs, table: tm, testCase: "µ:11.5,x:15.3,y:23", field: f1)

        // q series
        AssertLRE(marcum(µ: 11.5, x: 15.3, y: 29).p, "0.6555891257392535", digits: 11.0, resultStore: rs, table: tm, testCase: "µ:11.5,x:15.3,y:29", field: f)
        AssertLRE(inv_marcum(µ: 11.5, x: 15.3, p: .p(0.6555891257392535)), "29", exact: true, digits: 11.4, resultStore: rs, table: tm, testCase: "µ:11.5,x:15.3,y:29", field: f1)

        // p recursion
        AssertLRE(marcum(µ: 25, x: 35, y: 49).p, "0.1258610027087132", digits: 13.7, resultStore: rs, table: tm, testCase: "µ:25,x:35,y:49", field: f)
        AssertLRE(inv_marcum(µ: 25, x: 35, p: .p(0.1258610027087132)), "49", exact: true, resultStore: rs, table: tm, testCase: "µ:25,x:35,y:49", field: f1)

        // q recursion
        AssertLRE(marcum(µ: 25, x: 35, y: 65).p, "0.7079055833201373", digits: 14.0, resultStore: rs, table: tm, testCase: "µ:25,x:35,y:65", field: f)
        AssertLRE(inv_marcum(µ: 25, x: 35, p: .p(0.7079055833201373)), "65", exact: true, digits: 14.5, resultStore: rs, table: tm, testCase: "µ:25,x:35,y:65", field: f1)
    }
    
    /// More inverse Marcum Q cases on the tails
    ///
    /// Reference from Mathematica using the following formulation:
    /// Solve[-1e-30 + MarcumQ[15.3, Sqrt[2 * 11.5], 0, Sqrt[2 * x]] == 0, {x}]
    func testInverseMarcum() {
        let t = "Marcum Q inverse, tail values. µ: 15.3, x: 11.5, p varies. Reference values from Mathematica"
        AssertLRE(inv_marcum(µ: 15.3, x: 11.5, p: .p(1e-10)), "3.26354", resultStore: rs, table: t, testCase: "p:1e-10", field: f1)
        AssertLRE(inv_marcum(µ: 15.3, x: 11.5, p: .p(1e-20)), "0.690816", resultStore: rs, table: t, testCase: "p:1e-20", field: f1)
        AssertLRE(inv_marcum(µ: 15.3, x: 11.5, p: .p(1e-30)), "0.15206", resultStore: rs, table: t, testCase: "p:1e-30", field: f1)

        AssertLRE(inv_marcum(µ: 15.3, x: 11.5, p: .q(1e-10)), "83.7808", resultStore: rs, table: t, testCase: "q:1e-10", field: f1)
        AssertLRE(inv_marcum(µ: 15.3, x: 11.5, p: .q(1e-20)), "122.571", resultStore: rs, table: t, testCase: "q:1e-20", field: f1)
        AssertLRE(inv_marcum(µ: 15.3, x: 11.5, p: .q(1e-30)), "157.415", digits: 5.3, resultStore: rs, table: t, testCase: "q:1e-30", field: f1)
    }

    static var allTests = [
        ("testInvErf", testInvErf),
    ]
}
