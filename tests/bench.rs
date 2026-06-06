#![feature(test)]
#![allow(non_snake_case)]

use modular_inversion_gcd_study::{
    mod_inverse_A, mod_inverse_A2, mod_inverse_C, mod_inverse_C2, mod_inverse_D, mod_inverse_D2,
};

extern crate test;
use test::Bencher;

#[bench]
fn bench_mod_inverse_A(b: &mut Bencher) {
    let modulus = 2_559_031_471;
    let unit = 1_956_912_061;

    b.iter(|| _ = mod_inverse_A(modulus, unit));
}

#[bench]
fn bench_mod_inverse_A2(b: &mut Bencher) {
    let modulus = 2_559_031_471;
    let unit = 1_956_912_061;

    b.iter(|| _ = mod_inverse_A2(modulus, unit));
}

#[bench]
fn bench_mod_inverse_C(b: &mut Bencher) {
    let modulus = 2_559_031_471;
    let unit = 1_956_912_061;

    b.iter(|| _ = mod_inverse_C(modulus, unit));
}

#[bench]
fn bench_mod_inverse_C2(b: &mut Bencher) {
    let modulus = 2_559_031_471;
    let unit = 1_956_912_061;

    b.iter(|| _ = mod_inverse_C2(modulus, unit));
}

#[bench]
fn bench_mod_inverse_D(b: &mut Bencher) {
    let modulus = 2_559_031_471;
    let unit = 1_956_912_061;

    b.iter(|| _ = mod_inverse_D(modulus, unit));
}

#[bench]
fn bench_mod_inverse_D2(b: &mut Bencher) {
    let modulus = 2_559_031_471;
    let unit = 1_956_912_061;

    b.iter(|| _ = mod_inverse_D2(modulus, unit));
}

// rustup default nightly-x86_64-unknown-linux-gnu
// rustup default stable-x86_64-unknown-linux-gnu
//
// cargo fmt && cargo bench --test bench
// cargo fmt && cargo test
