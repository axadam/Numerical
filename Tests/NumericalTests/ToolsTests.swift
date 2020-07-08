//
//  ToolsTests.swift
//  NumericalTests
//
//  Created by Adam Roberts on 4/19/20.
//

import Foundation

import XCTest
import LogRelativeError
@testable import Numerical

final class ToolsTests: XCTestCase {
    override class func tearDown() {
        print(resultStore.md())
        super.tearDown()
    }

    static let resultStore: ResultStore = ResultStore()
    var rs: ResultStore { return type(of: self).resultStore }

    func testContinuedFraction() {
        let t = "Continued Fraction"
        let f = "LRE"
        
        let r2 = continuedFraction(b0: 1, a: { _ in 1 }, b: { _ in 2 })
        AssertLRE(r2.value, "1.414213562373095048801688724209698078569671875376948073176", resultStore: rs, table: t, testCase: "√2", field: f, annotation: "\(r2.iterations)\(r2.converged ? "" : "*")")
        
        let e  = continuedFraction(b0: 2, a: { $0 == 1 ? 1.0 : Double($0 - 1) }, b: { Double($0) })
        AssertLRE(e.value,  "2.718281828459045235360287471352662497757247093699959574966", resultStore: rs, table: t, testCase: "e", field: f, annotation: "\(e.iterations)\(e.converged ? "" : "*")")

        let pi = continuedFraction(b0: 3, a: { pow(Double(2 * $0 - 1),2) }, b: { _ in 6 }, maxIter: 1000)
        AssertLRE(pi.value, "3.141592653589793238462643383279502884197169399375105820974", digits: 10.1, resultStore: rs, table: t, testCase: "π", field: f, annotation: "\(pi.iterations)\(pi.converged ? "" : "*")")
    }
    
    func testChebyshev() {
        let t = "Chebyshev, e^x test case from Approximation Theory Approximation Practice, Trefethen"
        let f = "exp(x)"
        let coeffsExp = [
        1.266065877752008,
        1.130318207984970,
        0.271495339534077,
        0.044336849848664,
        0.005474240442094,
        0.000542926311914,
        0.000044977322954,
        0.000003198436462,
        0.000000199212481,
        0.000000011036772,
        0.000000000550590,
        0.000000000024980,
        0.000000000001039,
        0.000000000000040,
        0.000000000000001
        ]
        let e = { (x: Double) -> Double in
            chebyshev(coeffs: coeffsExp, z: x)
        }
        AssertLRE(e(0.5), "\(exp(0.5))", digits: 14.9, resultStore: rs, table: t, testCase: "0.5", field: f)
        AssertLRE(e(-0.5), "\(exp(-0.5))", digits: 14.7, resultStore: rs, table: t, testCase: "-0.5", field: f)
    }
    
    func testHorners() {
        let t = "Horner's Method of polynomial evaluation"
        let coefs: [Double] = [-19,7,-4,6]
        let r = polynomial(coeffs: coefs, z: 3)
        AssertLRE(r, "128", exact: true, resultStore: rs, table: t, testCase: "-19 + 7x - 4x² + 6x³; x=3", field: "LRE")
    }
    
    func testProduct() {
        let t = "Infinite Product"
        let f = "LRE"
        
        // Viète's formula: 2 / π = √2/2 √(2 + √2)/2 √(2 + √(2 + √2)) / 2 ...
        let v = product(indices: 1..., initialState: 2.0.squareRoot()) { i, prev in
            return (prev / 2,(2 + prev).squareRoot())
        }
        AssertLRE(v.value, "0.636619772367581343075535053490057448137838582961825794990", resultStore: rs, table: t, testCase: "Viète's formula for 2/π", field: f, annotation: "\(v.iterations)\(v.converged ? "" : "*")")
        
        // Wallis product: π / 2 = ∏ n=1... 4n² / (4n² - 1)
        let w = product(indices: 1..., initialState: ()) { i, _ in
            let t = Double(4 * i * i) / Double(4 * i * i - 1)
            return (t,())
        }
        AssertLRE(w.value, "1.570796326794896619231321691639751442098584699687552910487", digits: 2.6, resultStore: rs, table: t, testCase: "Wallis product for π/2", field: f, annotation: "\(w.iterations)\(w.converged ? "" : "*")")

        // Verification of Wallis product truncated at n=100 since it is so slow
        // "An investigation of the comparative efficiency of the different methods in which π is calculated", Nouri Al-Othman, 2013, §8.1, Table 1, n=100
        AssertLRE(w.value, "1.566893745314081", resultStore: rs, table: t, testCase: "Wallis product truncated at n=100", field: f, annotation: "\(w.iterations)\(w.converged ? "" : "*")")
    }
    
    func testSeries() {
        let t = "Infinite Series"
        let f = "LRE"
        
        // ln 2 = 2 Σ i=0... 3⁻²ⁱ⁻¹ / (2i + 1)
        let l = series(indices: 0..., initialState: ()) { i, _ in
            let j = 2 * i + 1
            let t = 1 / (3.0^^j * Double(j))
            return (t, ())
        }
        AssertLRE(2 * l.value, "0.693147180559945309417232121458176568075500134360255254120", resultStore: rs, table: t, testCase: "ln 2 = 2 Σ i=0... 3⁻²ⁱ⁻¹ / (2i + 1)", field: f, annotation: "\(l.iterations)\(l.converged ? "" : "*")")
        
        // Newton's arcsin series to find pi
        // π = 6 sin⁻¹ 1/2 = Σ i=0... 3 C(2i, i) / (16ⁱ (2i + 1))
        let p = series(indices: 0..., initialState: ()) { i, _ in
            let c = i == 0 ? 1 : zip((1...i),((i+1)...(2*i))).reduce(1) { $0 * Double($1.1) / Double($1.0) }
            let t = 3 * c / 16.0^^i / Double(2 * i + 1)
            return (t,())
        }
        AssertLRE(p.value, "3.141592653589793238462643383279502884197169399375105820974", resultStore: rs, table: t, testCase: "Newton's arcsin series for π", field: f, annotation: "\(p.iterations)\(p.converged ? "" : "*")")
        
        // Chudnovsky algorithm for pi
        // 426880 √10005 / π = Σ i=0... Mᵢ Lᵢ / Xᵢ,
        // Mᵢ₊₁ = Mᵢ (Kᵢ^3 - 16Kᵢ) / (i + 1)^3, M₀ = 1
        // Lᵢ₊₁ = Lᵢ + 545140134, L₀ = 13591409
        // Kᵢ₊₁ = Kᵢ + 12, K₀ = 6
        // Xᵢ₊₁ = Xᵢ (-262537412640768000), X₀ = 1
        let c = series(indices: 0..., initialState: (Kᵢ: 6.0, Lᵢ: 13591409.0, Mᵢ: 1.0, Xᵢ: 1.0)) { i, prev in
            let (Kᵢ,Lᵢ,Mᵢ,Xᵢ) = prev
            let t = Mᵢ * Lᵢ / Xᵢ
            let Kᵢ₊₁ = Kᵢ + 12
            let Lᵢ₊₁ = Lᵢ + 545140134
            let Mᵢ₊₁ = Mᵢ * (Kᵢ^^3 - 16.0 * Kᵢ) / Double(i + 1)^^3
            let Xᵢ₊₁ = Xᵢ * (-262537412640768000)
            return (t, (Kᵢ₊₁,Lᵢ₊₁,Mᵢ₊₁,Xᵢ₊₁))
        }
        AssertLRE(426880 * 10005.0.squareRoot() / c.value, "3.141592653589793238462643383279502884197169399375105820974", resultStore: rs, table: t, testCase: "Chudnovsky algorithm for π", field: f, annotation: "\(c.iterations)\(c.converged ? "" : "*")")
        
        // e = Σ i=0... 1 / i!
        let e = series(indices: 0..., initialState: 1.0) { i, prev in
            let fac = i == 0 ? 1 : prev * Double(i)
            return (1 / fac, fac)
        }
        AssertLRE(e.value,"2.718281828459045235360287471352662497757247093699959574966", resultStore: rs, table: t, testCase: "e = Σ i=0... 1 / i!", field: f, annotation: "\(e.iterations)\(e.converged ? "" : "*")")
    }

    let tSum = "Summation"
    
    let hardSum: [Double] = { () -> [Double] in
        let rand = randomArray(seed: 123, n: 10000, exponentRange: -10...10)
        return (rand.map { [$0,-$0] }.flatMap() { $0 } + [1.0]).shuffled(seed: 456)
    }()

    let easySum: [Double] = { () -> [Double] in
        let rand = randomArray(seed: 123, n: 10000, exponentRange: -0...0)
        return (rand.map { [$0,-$0] }.flatMap() { $0 } + [1.0]).shuffled(seed: 456)
    }()
    
    let petersSum: [Double] = [1.0, 1e100, 1.0, -1e100]

    func testNaive() {
        let f = "Naive"
        AssertLRE(easySum.sum_naive(), "1.0", exact: true, digits: 12.9, resultStore: rs, table: tSum, testCase: f, field: "Easy")
        AssertLRE(hardSum.sum_naive(), "1.0", exact: true, digits: 4.1, resultStore: rs, table: tSum, testCase: f, field: "Hard")
        AssertLRE(petersSum.sum_naive(), "2.0", exact: true, digits: 0.0, resultStore: rs, table: tSum, testCase: f, field: "Peters")
    }
    
    func testPairwise() {
        let f = "Pairwise"
        AssertLRE(easySum.sum_pairwise(), "1.0", exact: true, resultStore: rs, table: tSum, testCase: f, field: "Easy")
        AssertLRE(hardSum.sum_pairwise(), "1.0", exact: true, digits: 4.2, resultStore: rs, table: tSum, testCase: f, field: "Hard")
        AssertLRE(petersSum.sum_pairwise(), "2.0", exact: true, digits: 0.0, resultStore: rs, table: tSum, testCase: f, field: "Peters")
    }
    
    func testKahan() {
        let f = "Kahan"
        AssertLRE(easySum.sum_kahan(), "1.0", exact: true, resultStore: rs, table: tSum, testCase: f, field: "Easy")
        AssertLRE(hardSum.sum_kahan(), "1.0", exact: true, digits: 5.5, resultStore: rs, table: tSum, testCase: f, field: "Hard")
        AssertLRE(petersSum.sum_kahan(), "2.0", exact: true, digits: 0.0, resultStore: rs, table: tSum, testCase: f, field: "Peters")
    }

    func testKBN() {
        let f = "Kahan-Babuška-Neumaier"
        AssertLRE(easySum.sum_kbn(), "1.0", exact: true, resultStore: rs, table: tSum, testCase: f, field: "Easy")
        AssertLRE(hardSum.sum_kbn(), "1.0", exact: true, resultStore: rs, table: tSum, testCase: f, field: "Hard")
        AssertLRE(petersSum.sum_kbn(), "2.0", exact: true, resultStore: rs, table: tSum, testCase: f, field: "Peters")
    }
}
