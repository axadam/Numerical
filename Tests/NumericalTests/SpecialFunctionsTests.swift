import XCTest
@testable import Numerical

final class SpecialFunctionsTests: XCTestCase {
    func testInvErf() {
        let a = erf(1.2345)
        let b = invErf(a)
        XCTAssertEqual(b, 1.2345, accuracy: 1e-10)
    }
    
    func testIncompleteGamma() {
        let a = q_gamma(4, 0.7)
        // Q(4,7) = 0.9942465424077004166914...
        XCTAssertEqual(a, 0.9942465424077004166914, accuracy: 1e-10)
    }
    
    func testInverseIncompleteGamma() {
        let a = inv_p_gamma(4, 1 - 0.9942465424077004166914)
        XCTAssertEqual(a, 0.7, accuracy: 1e-10)
        
        // "GammaCHI: a package for the inversion and computation of the gamma
        // and chi-square cumulative distribution functions (central and noncentral)",
        // Gil, Tegura, and Temme 2015, Table 1
        XCTAssertEqual(q_gamma(1e-249, 6.310e-15), 3.212101109661167e-248 , accuracy: 1e-251)
        XCTAssertEqual(q_gamma(1e-249, 7.110e-7 ), 1.3580785912009393e-248, accuracy: 1e-251)
        XCTAssertEqual(q_gamma(1e-249, 0.01     ), 4.0379295765381135e-249, accuracy: 1e-251)
        XCTAssertEqual(q_gamma(1e-13, 6.310e-15 ), 3.212101109660651e-12  , accuracy: 1e-14)
        XCTAssertEqual(q_gamma(1e-13, 7.110e-7  ), 1.358078591200848e-12  , accuracy: 1e-14)
        XCTAssertEqual(q_gamma(1e-13, 0.01      ), 4.0379295765380405e-13 , accuracy: 1e-25)
    }
    
    func testIncompleteBeta() {
        let a = beta_reg(x: 0.4, a: 3, b: 5)
        // I₀ٖ₄(3,5) = 0.580096
        XCTAssertEqual(a, 0.580096, accuracy: 1e-10)
    }
    
    func testInverseIncompleteBeta() {
        let a = inv_beta_reg(p: 0.580096, a: 3, b: 5)
        XCTAssertEqual(a, 0.4, accuracy: 1e-10)
    }

    /// Marcum Q asymptotic large µ test cases
    ///
    /// Ten tests of the large µ case for the Marcum Q function. The varied dimension
    /// is the non-centrality parameter.
    ///
    /// "Computation of the Marcum Q-function", Gil, Segura, Temme 2013, Table 6.1
    func testMarcumBigMu() {
        let µ: Double = 8192
        let y = 1.05 * µ
        func x(_ x: Double) -> Double { return x * 8192 }
        XCTAssertEqual(marcum(µ: µ, x: x(0.01), y: y).q, 1.9845278031193e-4  , accuracy: 1e-10)
        XCTAssertEqual(marcum(µ: µ, x: x(0.02), y: y).q, 4.138241872117e-3   , accuracy: 1e-10)
        XCTAssertEqual(marcum(µ: µ, x: x(0.03), y: y).q, 0.04000364971081    , accuracy: 1e-10)
        XCTAssertEqual(marcum(µ: µ, x: x(0.04), y: y).q, 0.191650654805848   , accuracy: 1e-10)
        XCTAssertEqual(marcum(µ: µ, x: x(0.05), y: y).q, 0.498535453743169   , accuracy: 1e-10)
        XCTAssertEqual(marcum(µ: µ, x: x(0.06), y: y).p, 0.1964796269915073  , accuracy: 1e-10)
        XCTAssertEqual(marcum(µ: µ, x: x(0.07), y: y).p, 0.04434265824612003 , accuracy: 1e-10)
        XCTAssertEqual(marcum(µ: µ, x: x(0.08), y: y).p, 5.526239087335513e-3, accuracy: 1e-10)
        XCTAssertEqual(marcum(µ: µ, x: x(0.09), y: y).p, 3.750276163593746e-4, accuracy: 1e-10)
        XCTAssertEqual(marcum(µ: µ, x: x(0.10), y: y).p, 1.386276448162126e-5, accuracy: 1e-10)
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
        XCTAssertEqual(marcum(µ: 5 , x: 150, y: 30 ).p, 1.215915354045e-23, accuracy: 1e-35)
        XCTAssertEqual(marcum(µ: 1 , x: 75,  y: 0.5).p, 3.287840255874e-30, accuracy: 1e-42)
        XCTAssertEqual(marcum(µ: 2 , x: 100, y: 2  ).p, 1.557081489535e-35, accuracy: 1e-45)
        XCTAssertEqual(marcum(µ: 10, x: 100, y: 1  ).p, 5.152185145235e-48, accuracy: 1e-60)
    }
    
    /// Marcum Q test cases for various methods
    ///
    /// Method used by each test case is noted in comments. Reference values from
    /// Mathematica using the syntax: N[MarcumQ[µ,sqrt(2 * x),0,sqrt(2 * y)]]
    func testMarcum() {
        // p series
        XCTAssertEqual(marcum(µ: 11.5, x: 15.3, y: 23).p, 0.2948691834572695, accuracy: 1e-12)
        
        // q series
        XCTAssertEqual(marcum(µ: 11.5, x: 15.3, y: 29).p, 0.6555891257392535, accuracy: 1e-11)
        
        // p recursion
        XCTAssertEqual(marcum(µ: 25, x: 35, y: 49).p, 0.1258610027087132, accuracy: 1e-11)

        // q recursion
        XCTAssertEqual(marcum(µ: 25, x: 35, y: 65).p, 0.7079055833201373, accuracy: 1e-11)
    }

    
    static var allTests = [
        ("testInvErf", testInvErf),
    ]
}
