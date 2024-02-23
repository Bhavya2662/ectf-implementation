use curv::elliptic::curves::{secp256_k1::Secp256k1, Point, Scalar};
use curv::arithmetic::traits::Group;
use curv::cryptographic_primitives::hashing::hash_sha256;

// Input: P1 = (x1,y1) E EC(Fp) from P, P2 = (x2,y2) Є EC(Fp) from V.
// Output: P and V output s1 and s2 such that s1+s2 =x where (x,y) = P1+P2 in EC.
pub fn ectf_protocol(p1: &Point<Secp256k1>, p2: &Point<Secp256k1>) -> (Scalar<Secp256k1>, Scalar<Secp256k1>) {
    // P (and V) sample p; <- Zp for i ∈ {1,2} respectively.
    let (p1, p2) = (Scalar::random(), Scalar::random());

    // P and V run a1, α2 := MtA( (−x1, P1),(p2,x2)).
    let (a1, a2) = mta(&(-p1.x, p1), &(p2, p2.x));

    // P computes δ1 =-x1* p1+α1 and V computes δ2 =x2 * p2+ α2.
    let (d1, d2) = (-p1.x * p1 + a1, p2.x * p2 + a2);

    // P (and V) reveal δ1 (and δ2) to each other and compute δ = δ1 + δ2 .
    let d = d1 + d2;

    // P (and V) compute ηi = pi* δ-1 for i = {1,2} respectively.
    let (n1, n2) = (p1 * d.invert(), p2 * d.invert());

    //P and V run B1, B2 :=MtA( (-y1, N1), (N2,y2 )).
    let (b1, b2) = mta(&(-p1.y, n1), &(n2, p2.y));

    // P computes λ1 =-y1 * n1+ B1 and V computes 1 =y2*n2+ B2. They run Y1, Y2 =MtA(A1, A2).
    let (l1, l2) = (-p1.y * n1 + b1, p2.y * n2 + b2);
    let (y1, y2) = mta(&l1, &l2);

    // P (and V) computes si =2; + λ;2 - xi for i = {1,2} respectively.
    let (s1, s2) = (n1 + l1.square() - p1.x, n2 + l2.square() - p2.x);

    // P outputs s1 and V outputs s2.
    (s1, s2)
}

fn main() {
    // Generate two points on the secp256k1 curve.
    let p1 = Point::generator();
    let p2 = Point::random();

    // Run the ECtF protocol.
    let (s1, s2) = ectf_protocol(&p1, &p2);

    // Check if s1 + s2 = x where (x, y) = p1 + p2.
    let p3 = p1 + p2;
    assert_eq!(s1 + s2, p3.x);

    println!("ECtF protocol completed successfully.");
}
