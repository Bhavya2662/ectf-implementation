// use curv::elliptic::curves::{secp256_k1::Secp256k1, Point, Scalar};

// use curv::cryptographic_primitives::hashing::hash_sha256;
#![allow(unused_imports)]
#![allow(non_snake_case)]
#![allow(unused_variables)]
#![allow(unused_parens)]
#![allow(dead_code)]
use curv::arithmetic::Integer;
use curv::arithmetic::NumberTests;
use paillier::*;
use curv::BigInt;
use curv::arithmetic::BasicOps;
use curv::arithmetic::Modulo;
use curv::arithmetic::traits::Converter;
use curv::arithmetic::Zero;
use curv::arithmetic::Samplable;
use curv::elliptic::curves::{secp256_k1::Secp256k1, Point, Scalar};
use rand::Rng;


fn mta(a: &Scalar<Secp256k1>, b: &Scalar<Secp256k1>) -> (Scalar<Secp256k1>, Scalar<Secp256k1>){
    let zp = BigInt::from_str_radix("115792089237316195423570985008687907853269984665640564039457584007908834671663", 10).unwrap();
    let order = Scalar::<Secp256k1>::group_order();

    let (ek, dk) = Paillier::keypair().keys();

    // Alice's input
    let alice_input= a;
    // Bob's input
    let bob_input= b;

    // Alice computes cA = EncryptA(a)
    let c_alice = Paillier::encrypt(&ek, RawPlaintext::from(&alice_input.to_bigint()));

    // Bob selects Beta Tag <– Z(N) -> where n is ek.n 
    // let res = &ek.n - ((order-BigInt::from(1))*(order-BigInt::from(1)));
    // dbg!(&res);
    let beta_tag: BigInt = BigInt::sample_below(&ek.n);
    dbg!(&beta_tag);

    // Compute Encrypt(BetaTag) using key of A
    let enc_betatag = Paillier::encrypt(&ek, RawPlaintext::from(&beta_tag));

    // Compute cB = b * cA + EncryptA(BetaTag) = EncryptA(ab+Tag)
    // Compute b * CA
    let b_mul_c_alice = Paillier::mul(&ek, RawPlaintext::from(&bob_input.to_bigint()), c_alice.clone());

    // Compute cB
    let c_bob = Paillier::add(&ek, b_mul_c_alice.clone(), enc_betatag.clone());

    // Bob sets additive share Beta = -BetaTag mod n
    let beta = Scalar::<Secp256k1>::from(&(&BigInt::zero() - &beta_tag));

    // Alice decrypts = dec(cB)
    let dec_alice = Paillier::decrypt(&dk, &c_bob);

    // Alice sets alpha = dec_alice mod n
    let alpha = Scalar::<Secp256k1>::from(BigInt::from(dec_alice.clone()));
    // let left = ((&alpha + &beta).to_bigint());
    // dbg!(&left);
    // let right =( (alice_input * bob_input).to_bigint());
    // dbg!(&right);
    // assert_eq!(left, right, "Verification failed: Left side ({}) is not equal to right side ({})", left, right);
    (alpha, beta)
}

fn mta_2(
    alice_secret: &Scalar<Secp256k1>,
    r1: &BigInt,
    r2: &BigInt,
    bob_secret: &Scalar<Secp256k1>, 
) -> (Scalar<Secp256k1>, Scalar<Secp256k1>) {
    let zp = BigInt::from_str_radix("115792089237316195423570985008687907853269984665640564039457584007908834671663", 10).unwrap();

    let order = Scalar::<Secp256k1>::group_order();

    // Generate Paillier key pair
    let (ek, dk) = Paillier::keypair().keys();

    // Alice's secret shares
    let a1 = &alice_secret.to_bigint();
    // Bob's secret shares
    let a2 = &bob_secret.to_bigint();
    
    // 1. Alice computes CA = EncryptA(a1)
    let c_a = Paillier::encrypt(&ek, RawPlaintext::from(&a1.clone()));
    // 2. Alice computes CR1 = EncryptA(r1)
    let c_r1 = Paillier::encrypt(&ek, RawPlaintext::from(r1.clone()));

    // 4. Bob selects β' <- ZN
    // let res = &ek.n - ((order-BigInt::from(1))*(order-BigInt::from(1)));
    let beta_prime = BigInt::sample_below(&ek.n);

    // 5. Bob computes CB = (r2 * CA) + (a2 * CR1) + EncryptA(β')
    let r2_mul_c_alice =  Paillier::mul(&ek, RawPlaintext::from(r2.clone()), c_a.clone());
    let a2_mul_c_r1 = Paillier::mul(&ek, RawPlaintext::from(a2.clone()), c_r1.clone());
    let enc_beta_prime = Paillier::encrypt(&ek, RawPlaintext::from(&beta_prime));


    let c_bob = Paillier::add(&ek, r2_mul_c_alice.clone(), a2_mul_c_r1.clone());
    let c_bob = Paillier::add(&ek, c_bob.clone(), enc_beta_prime.clone());

    // 6. Bob sets additive share δ2 = -β′ mod q
    let delta_2 = (BigInt::zero() - &beta_prime) ;

    // 8. Alice decrypts α' = dec(CB)
    let dec_alice = Paillier::decrypt(&dk, &c_bob);

    // 9. Alice sets δ 1 = α' mod q
    let delta_1 = BigInt::from(dec_alice);


    // let left = &delta_1 + &delta_2;
    // dbg!(&left);
    // let right = (a1*r2) + (a2*r1);
    // dbg!(&right);
    // assert_eq!(left, right, "Verification failed: Left side ({}) is not equal to right side ({})", left, right);

    // println!("MTA Test Satisfied");
    (delta_1.into(), delta_2.into())
}
fn ectf_protocol(
    p1: &Point<Secp256k1>,
    p2: &Point<Secp256k1>,
) -> (Scalar<Secp256k1>, Scalar<Secp256k1>) {
    let x1 = p1.x_coord().unwrap();
    let y1 = p1.y_coord().unwrap();
    let x2 = p2.x_coord().unwrap();
    let y2 = p2.y_coord().unwrap();

    // Alice and Bob sample random values
    let zp = BigInt::from_str_radix("115792089237316195423570985008687907853269984665640564039457584007908834671663", 10).unwrap();
    dbg!(&zp); 
    let order = Scalar::<Secp256k1>::group_order();

    let rho_1 = BigInt::sample_below(order);
    dbg!(&rho_1); 
    let rho_2 = BigInt::sample_below(order);
    dbg!(&rho_2);

    // Alice and Bob run MTA
    let (alpha_1, alpha_2) = mta_2(&Scalar::<Secp256k1>::from(-x1.clone()), &rho_1, &rho_2, &Scalar::<Secp256k1>::from(x2.clone()));
    dbg!(&alpha_1);
    dbg!(&alpha_2);
    // let aa = (-x1.clone() * rho_2.clone()) + (x2.clone() * rho_1.clone());
    // let aa = BigInt::from(2)*
    // Check 1 ----------------------------------------------------------------------------
    // assert_eq!((alpha_1.clone()+alpha_2.clone()).to_bigint(), ((-x1.clone() * rho_2.clone()) + (x2.clone() * rho_1.clone()))%order);
    // let alpha_sum = (alpha_1.clone().to_bigint() + alpha_2.clone().to_bigint()) % order.clone();
    // let mut sum_x1r2_x2r1 = ((-x1.clone() * rho_2.clone()) % order.clone()) + ((x2.clone()*rho_1.clone()% order.clone()));
    
    // // sum_x1r2_x2r1 = -sum_x1r2_x2r1.clone();

    // dbg!(&alpha_sum);
    // dbg!(&sum_x1r2_x2r1);
    // assert_eq!(alpha_sum.clone(), sum_x1r2_x2r1.clone(),"Check 1 Failed");
    dbg!((alpha_1.clone()+alpha_2.clone()).to_bigint());
    let right_ch1 = (alpha_1.clone()+alpha_2.clone()).to_bigint();
    let mut left_ch1 = ((-x1.clone() * rho_2.clone()) + (x2.clone() * rho_1.clone()))%order;

    if BigInt::is_negative(&left_ch1) {
        left_ch1 = order - (left_ch1.abs() % order);
    }

    dbg!(&left_ch1);
    dbg!(order - ( left_ch1.abs() % order));
    dbg!(((-x1.clone() * rho_2.clone()) + (x2.clone() * rho_1.clone()))%order);

    // Check 1 ----------------------------------------------------------------------------
    assert_eq!(right_ch1, left_ch1);
    println!("check 1 passed");

    dbg!((alpha_1.clone()+alpha_2.clone()).to_bigint());
    dbg!(((-x1.clone() * rho_2.clone()) + (x2.clone() * rho_1.clone()))%order);
    // Compute delta values
    let delta_1 = (-x1.clone() * rho_1.clone()) + alpha_1.clone().to_bigint();
    let delta_2 = (x2.clone() * rho_2.clone() )+ alpha_2.clone().to_bigint();
    dbg!(&delta_1);
    dbg!(-(5%4));
    dbg!(&delta_2);

    let delta = &delta_1 + &delta_2;
    dbg!(&delta);
    // let gcd = delta.gcd(&zp); 
    // if gcd != BigInt::from(1) {
    //     panic!("GCD of delta and zp is not 1, modular inverse does not exist.");
    // }
    let delta_inv = BigInt::mod_inv(&delta, order).unwrap();
    dbg!(&delta_inv);
    

    // Compute eta values
    let eta_1 = rho_1*delta_inv.clone();
    let eta_2 = rho_2*delta_inv.clone();
    dbg!(&eta_1);
    dbg!(&eta_2);
    // assert_eq!((eta_1.clone() + eta_2.clone())%order, (y2.clone() - y1.clone()) / (x2.clone()-x1.clone()));
    // Check 2------------------------------------------------------------------------------
    assert_eq!((eta_1.clone() + eta_2.clone())%order, (BigInt::mod_inv(&(x2.clone() - x1.clone()), order).unwrap())%order );
    // Run MTA for beta values
    let (beta_1, beta_2) = mta_2(&Scalar::<Secp256k1>::from(-y1.clone()), &eta_1, &eta_2, &Scalar::<Secp256k1>::from(y2.clone()));
    dbg!(&beta_1);
    dbg!(&beta_2);
    // Compute lambda values
    let mut lambda_1 = ((-y1.clone() * eta_1) + beta_1.clone().to_bigint());
    // let mut lambda_1_pos = (y2.clone()-y1.clone())*BigInt::mod_inv(&(x2.clone() - x1.clone()), order).unwrap()%order;

    if BigInt::is_negative(&lambda_1) {
        lambda_1 = order - (lambda_1.abs() % order);
    }
    let mut lambda_2 =( (y2.clone() * eta_2 )+ beta_2.clone().to_bigint())%order;
    if BigInt::is_negative(&lambda_2) {
        lambda_2 = order - (lambda_2.abs() % order);
    }
    dbg!(&lambda_1);
    dbg!(&lambda_2);
    let lambda = lambda_1.clone() + lambda_2.clone();
    let slope = (y2.clone()-y1.clone())*BigInt::mod_inv(&(x2.clone() - x1.clone()), order).unwrap()%order;
    // Check 3----------------------------------------------------------------------------
    // assert_eq!((lambda_1.clone()+ lambda_2.clone())%order, ((y2.clone()-y1.clone())*BigInt::mod_inv(&(x2.clone() - x1.clone()), order).unwrap())%order );
    let right_ch3 = (lambda_1.clone()+ lambda_2.clone())%order;
    let mut left_ch3 = (y2.clone()-y1.clone())*BigInt::mod_inv(&(x2.clone() - x1.clone()), order).unwrap()%order;

    if BigInt::is_negative(&left_ch3) {
        left_ch3 = order - (left_ch3.abs() % order);
    }

    dbg!(&left_ch3);
    dbg!(order - ( left_ch3.abs() % order));
    dbg!((y2.clone()-y1.clone())*BigInt::mod_inv(&(x2.clone() - x1.clone()), order).unwrap()%order);
    assert_eq!(right_ch1, left_ch1);

    // Run MTA for gamma values
    let (gamma_1, gamma_2) = mta(&Scalar::<Secp256k1>::from(lambda_1.clone()), &Scalar::<Secp256k1>::from(lambda_2.clone()));
    dbg!(&gamma_1);
    dbg!(&gamma_2);
    let left = ((&gamma_1 + &gamma_2).to_bigint());
    dbg!(&left);
    let right =( (&Scalar::<Secp256k1>::from(lambda_1.clone()) * &Scalar::<Secp256k1>::from(lambda_2.clone())).to_bigint());
    dbg!(&right);
    assert_eq!(left, right, "Verification failed: Left side ({}) is not equal to right side ({})", left, right);

    let mut s = (slope.clone() * slope.clone() - x1.clone() - x2.clone())%order;
    if BigInt::is_negative(&s) {
        s = order - (s.abs() % order);
    }
    dbg!(&s);

    // Compute final output s values
    let mut s1 = ((BigInt::from(2) * gamma_1.to_bigint()) +( (lambda_1.clone()*lambda_1.clone()) - x1.clone()))%order;
    if BigInt::is_negative(&s1) {
        s1 = order - (s1.abs() % order);
    }
    let mut s2 = ((BigInt::from(2) * gamma_2.to_bigint()) + ((lambda_2.clone()*lambda_2.clone()) - x2.clone()))%order;
    if BigInt::is_negative(&s2) {
        s2 = order - (s2.abs() % order);
    }
    dbg!(&s1);
    dbg!(&s2);
    assert_eq!((s1.clone() + s2.clone())%order, s);

    // Check 4----------------------------------------------------------------------------
    assert_eq!((s1.clone() + s2.clone())%order, ((lambda.clone()*lambda.clone() )- x1.clone() - x2.clone())%order);
    (s1.into(), s2.into())
}
fn main() {
    let scalar_1 = Scalar::<Secp256k1>::random();
    let scalar_2 = Scalar::<Secp256k1>::random();
    let g_1 = Point::<Secp256k1>::generator();
    let g_2 = Point::<Secp256k1>::generator();
// mod p
    let zp = BigInt::from_str_radix("115792089237316195423570985008687907853269984665640564039457584007908834671663", 10).unwrap();
    let order = Scalar::<Secp256k1>::group_order();

    let p1: Point<Secp256k1> = scalar_1 * g_1;
    let p2: Point<Secp256k1> = scalar_2 * g_2;
    // let p1_x = p1.x_coord().unwrap();
    // let p1_y = p1.y_coord().unwrap();
    // dbg!(&p1_x%order);
    // dbg!(&p1_x);
    // dbg!(&p1_y%order);
    // dbg!(&p1_y);

    // Generate two points on the secp256k1 curve.
    // let p1 = Point::<Secp256k1>::base_point2(); // user
    // let p2 = Point::<Secp256k1>::generator().to_point(); // oracle
    // dbg!(&p1);
    // dbg!(&p2);
    // let x1 = p1.x_coord().unwrap();
    // dbg!(&x1);
    // let zp = BigInt::from_str_radix("115792089237316195423570985008687907853269984665640564039457584007908834671663", 10).unwrap();
    // dbg!(&x1%&zp);
    // Run the ECtF protocol.
    let (s1, s2) = ectf_protocol(&p1, &p2);

    // Check if s1 + s2 = x where (x, y) = p1 + p2.
    let p3 = p1 + p2;

    assert_eq!(((s1 + s2).to_bigint())%order, p3.x_coord().unwrap());

    println!("ECtF protocol completed successfully.");
}
