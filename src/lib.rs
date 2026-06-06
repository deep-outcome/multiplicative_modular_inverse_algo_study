#![allow(dead_code)]
//! If `gcd(left,right) =1`, then linear combination, linear Diophantine equation,
//! `l ⋅left +r ⋅right =1` expresses respective modular multiplicative inverses:
//!
//! ⏵ `l ⋅left  ≡1 (mod right)`, `l ≡left⁻¹ (mod right)`<br/>
//! ⏵ `r ⋅right ≡1 (mod left)`,  `r ≡right⁻¹ (mod left)`
//!
//! Inverses are called Bézout's coefficients and equation Bézout's identity.
//!
//! Note: Bézout's coefficients and Bézout's identity are not limited to `gcd(left,right) =1`,
//! yet both apply to all gcds.
use core::convert::TryInto;
type InverseRes = Result<Option<usize>, &'static str>;

const fn val_rels(modulus: usize, unit: usize) -> Result<(), &'static str> {
    return if modulus & 1 == 0 {
        Err("Odd modulus required.")
    } else if modulus == 1 {
        Err("Modulus is unity.")
    } else if modulus <= unit {
        Err("Unit must be less than modulus.")
    } else {
        Ok(())
    };
}

const fn val_invs(
    modulus: usize,
    unit: usize,
    g: usize,
    h: usize,
    a: usize,
    b: usize,
) -> Result<(), &'static str> {
    let modulus = modulus as u128;
    let unit = unit as u128;
    let g = g as u128;
    let h = h as u128;
    let a = a as u128;
    let b = b as u128;

    let a_aux = (g * unit) % modulus;
    let b_aux = (h * unit) % modulus;

    return if a != a_aux {
        Err("`a` invariant failed.")
    } else if b != b_aux {
        Err("`b` invariant failed.")
    } else {
        Ok(())
    };
}

// inverse is know in iteration when either of gcd operand reaches 1
const fn val_coeffs(g: usize, h: usize) -> Result<(), &'static str> {
    if g == h {
        Ok(())
    } else {
        Err("Coefficients must equal.")
    }
}

#[allow(non_snake_case)]
pub fn mod_inverse_A(modulus: usize, unit: usize) -> InverseRes {
    #[cfg(test)]
    #[cfg(feature = "extended-tests")]
    if let Err(e) = val_rels(modulus, unit) {
        return Err(e);
    }

    let mut a = unit as isize;
    let mut b = modulus as isize;

    let mut g = 1;
    let mut h = 0;

    let mut swap;

    while a != 0 {
        if a & 1 == 0 {
            a >>= 1;

            if g & 1 == 1 {
                g += modulus as isize;
            }

            g >>= 1;
        } else {
            if a < b {
                swap = a;
                a = b;
                b = swap;

                swap = g;
                g = h;
                h = swap;
            }

            a -= b;
            a >>= 1;

            g = g - h;
            if g & 1 == 1 {
                g += modulus as isize;
            }

            g >>= 1
        }
    }

    assert_eq!(true, TryInto::<usize>::try_into(h).is_ok());
    let res = if b == 1 { Some(h as usize) } else { None };
    Ok(res)
}

#[allow(non_snake_case)]
pub const fn mod_inverse_B(modulus: usize, unit: usize) -> InverseRes {
    #[cfg(test)]
    #[cfg(feature = "extended-tests")]
    if let Err(e) = val_rels(modulus, unit) {
        return Err(e);
    }

    let mut a = unit;
    let mut b = modulus;

    let mut g = 1;
    let mut h = 0;

    let mut swap;

    // e.g. 51 ≡2⁻¹ (mod 101)
    let unit_2_inverse = (modulus + 1) / 2;

    while a != 0 {
        if a & 1 == 0 {
            let ctz = a.trailing_zeros();
            a >>= ctz;
            g = g * (unit_2_inverse.pow(ctz)) % modulus;
        }
        if a < b {
            swap = a;
            a = b;
            b = swap;

            swap = g;
            g = h;
            h = swap;
        }

        a -= b;
        a >>= 1;

        #[cfg(test)]
        #[cfg(feature = "extended-tests")]
        if a == 0 && b == 1 {
            if let Err(e) = val_coeffs(g, h) {
                return Err(e);
            }
        }

        g = sub_halve_mod(g, h, modulus);

        #[cfg(test)]
        #[cfg(feature = "extended-tests")]
        if let Err(e) = val_invs(modulus, unit, g, h, a, b) {
            return Err(e);
        }
    }

    let res = if b == 1 { Some(h) } else { None };
    Ok(res)
}

#[allow(non_snake_case)]
pub const fn mod_inverse_C(modulus: usize, unit: usize) -> InverseRes {
    #[cfg(test)]
    #[cfg(feature = "extended-tests")]
    if let Err(e) = val_rels(modulus, unit) {
        return Err(e);
    }

    let mut a = unit;
    let mut b = modulus;

    let mut g = 1;
    let mut h = 0;

    let mut swap;

    while a != 0 {
        if a & 1 == 0 {
            let mut ctz = a.trailing_zeros();
            a >>= ctz;

            while ctz > 0 {
                ctz -= 1;

                if g & 1 == 1 {
                    g += modulus;
                }

                g >>= 1;
            }
        }
        if a < b {
            swap = a;
            a = b;
            b = swap;

            swap = g;
            g = h;
            h = swap;
        }

        a -= b;
        a >>= 1;

        #[cfg(test)]
        #[cfg(feature = "extended-tests")]
        if a == 0 && b == 1 {
            if let Err(e) = val_coeffs(g, h) {
                return Err(e);
            }
        }

        g = sub_halve_mod(g, h, modulus);

        #[cfg(test)]
        #[cfg(feature = "extended-tests")]
        if let Err(e) = val_invs(modulus, unit, g, h, a, b) {
            return Err(e);
        }
    }

    let res = if b == 1 { Some(h) } else { None };
    Ok(res)
}

#[allow(non_snake_case)]
pub const fn mod_inverse_C2(modulus: usize, unit: usize) -> InverseRes {
    #[cfg(test)]
    #[cfg(feature = "extended-tests")]
    if let Err(e) = val_rels(modulus, unit) {
        return Err(e);
    }

    let mut a = unit;
    let mut b = modulus;

    let mut g = 1;
    let mut h = 0;

    let mut swap;

    let ok = loop {
        if a & 1 == 0 {
            let mut ctz = a.trailing_zeros();
            a >>= ctz;

            while ctz > 0 {
                ctz -= 1;

                if g & 1 == 1 {
                    g += modulus;
                }

                g >>= 1;
            }
        }
        if a < b {
            swap = a;
            a = b;
            b = swap;

            swap = g;
            g = h;
            h = swap;
        }

        a -= b;
        a >>= 1;

        if a == 0 {
            break None;
        }

        g = sub_halve_mod(g, h, modulus);

        #[cfg(test)]
        #[cfg(feature = "extended-tests")]
        if let Err(e) = val_invs(modulus, unit, g, h, a, b) {
            return Err(e);
        }

        if a == 1 {
            break Some(g);
        }
    };

    Ok(ok)
}

#[allow(non_snake_case)]
pub const fn mod_inverse_D(modulus: usize, unit: usize) -> InverseRes {
    #[cfg(test)]
    #[cfg(feature = "extended-tests")]
    if let Err(e) = val_rels(modulus, unit) {
        return Err(e);
    }

    let mut a = unit;
    let mut b = modulus;

    let mut g = 1;
    let mut h = 0;

    loop {
        if a & 1 == 0 {
            let mut ctz = a.trailing_zeros();
            a >>= ctz;

            while ctz > 0 {
                ctz -= 1;

                if g & 1 == 1 {
                    g += modulus;
                }

                g = g >> 1;
            }
        } else if b & 1 == 0 {
            let mut ctz = b.trailing_zeros();
            b >>= ctz;

            while ctz > 0 {
                ctz -= 1;
                if h & 1 == 1 {
                    h += modulus;
                }

                h = h >> 1;
            }
        }

        if a < b {
            b -= a;
            b >>= 1;
            h = sub_halve_mod(h, g, modulus);
        } else {
            a -= b;

            if a == 0 {
                break;
            }

            a >>= 1;
            g = sub_halve_mod(g, h, modulus);
        }

        #[cfg(test)]
        #[cfg(feature = "extended-tests")]
        if let Err(e) = val_invs(modulus, unit, g, h, a, b) {
            return Err(e);
        }
    }

    let res = if b == 1 {
        #[cfg(test)]
        #[cfg(feature = "extended-tests")]
        if let Err(e) = val_coeffs(g, h) {
            return Err(e);
        }
        Some(g)
    } else {
        None
    };

    return Ok(res);
}

#[allow(non_snake_case)]
pub fn mod_inverse_D2(modulus: usize, unit: usize) -> InverseRes {
    #[cfg(test)]
    #[cfg(feature = "extended-tests")]
    if let Err(e) = val_rels(modulus, unit) {
        return Err(e);
    }

    let mut a = unit;
    let mut b = modulus;

    let mut g = 1;
    let mut h = 0;

    loop {
        if a & 1 == 0 {
            let mut ctz = a.trailing_zeros();
            a >>= ctz;

            while ctz > 0 {
                ctz -= 1;

                if g & 1 == 1 {
                    g += modulus;
                }

                g = g >> 1;
            }
        } else if b & 1 == 0 {
            let mut ctz = b.trailing_zeros();
            b >>= ctz;

            while ctz > 0 {
                ctz -= 1;
                if h & 1 == 1 {
                    h += modulus;
                }

                h = h >> 1;
            }
        }

        if a < b {
            b -= a;
            b >>= 1;
            h = sub_halve_mod(h, g, modulus);

            #[cfg(test)]
            #[cfg(feature = "extended-tests")]
            if let Err(e) = val_invs(modulus, unit, g, h, a, b) {
                return Err(e);
            }

            if b == 1 {
                return Ok(Some(h));
            }
        } else {
            a -= b;

            if a == 0 {
                break;
            }

            a >>= 1;
            g = sub_halve_mod(g, h, modulus);

            #[cfg(test)]
            #[cfg(feature = "extended-tests")]
            if let Err(e) = val_invs(modulus, unit, g, h, a, b) {
                return Err(e);
            }

            if a == 1 {
                return Ok(Some(g));
            }
        }
    }

    return Ok(None);
}

const fn sub_halve_mod(mut l: usize, r: usize, m: usize) -> usize {
    if l < r {
        l += m;
    }

    l -= r;

    if l & 1 == 0 { l >> 1 } else { (l + m) >> 1 }
}

#[cfg(test)]
mod tests_of_units {

    mod val_rels {
        use crate::val_rels;

        #[test]
        fn non_odd_modulus() {
            for v in [0, 2] {
                let err = val_rels(v, 0);
                assert_eq!(Err("Odd modulus required."), err);
            }
        }

        #[test]
        fn unity_modulus() {
            let err = val_rels(1, 0);
            assert_eq!(Err("Modulus is unity."), err);
        }

        #[test]
        fn too_large_unit() {
            for u in [11, 12] {
                let err = val_rels(11, u);
                assert_eq!(Err("Unit must be less than modulus."), err);
            }
        }
    }

    mod val_invs {
        use crate::val_invs;

        #[test]
        fn ok() {
            let m = 35;
            let u = 24;
            let g = 2;
            let h = 4;
            let a = 13;
            let b = 26;

            let ok = val_invs(m, u, g, h, a, b);
            assert_eq!(Ok(()), ok);
        }

        #[test]
        fn invariant_failed() {
            let m = 35;
            let u = 24;
            let g = 2;
            let h = 4;

            for (a, b, c) in [(12, 26, 'a'), (13, 25, 'b')] {
                let msg = format!("`{c}` invariant failed.");
                let err = val_invs(m, u, g, h, a, b);
                assert_eq!(Err(msg.as_str()), err);
            }
        }
    }

    mod val_coeffs {
        use crate::val_coeffs;

        #[test]
        fn ok() {
            assert_eq!(Ok(()), val_coeffs(3, 3));
        }

        #[test]
        fn nok() {
            let err = Err("Coefficients must equal.");
            assert_eq!(err, val_coeffs(0, 1));
        }
    }

    mod mod_inverse {

        use crate::{
            InverseRes, mod_inverse_A, mod_inverse_B, mod_inverse_C, mod_inverse_C2, mod_inverse_D,
            mod_inverse_D2,
        };

        type InvCom = fn(usize, usize) -> InverseRes;
        type TestInf = (&'static str, InvCom);

        const A: TestInf = ("A", mod_inverse_A);
        const B: TestInf = ("B", mod_inverse_B);
        const C: TestInf = ("C", mod_inverse_C);
        const C2: TestInf = ("C2", mod_inverse_C2);
        const D: TestInf = ("D", mod_inverse_D);
        const D2: TestInf = ("D2", mod_inverse_D2);

        fn val(unit: u128, modulus: u128, inverse: InverseRes) {
            assert_eq!(
                false,
                inverse.is_err(),
                "Precondition failure, {:?}",
                inverse
            );
            assert_ne!(Ok(None), inverse, "Non-invertable value");

            let inverse = inverse.unwrap().unwrap() as u128;
            println!("inverse: {inverse}");
            let u_summand = unit * inverse;

            let mut res = eq_m_summand(u_summand + 1, modulus);
            assert_eq!(false, res);
            res = eq_m_summand(u_summand - 1, modulus);
            assert_eq!(true, res);

            fn eq_m_summand(sum: u128, m: u128) -> bool {
                let ratio = sum / m;
                let m_summand = ratio * m;
                sum == m_summand
            }
        }

        fn test(unit: usize, modulus: usize, f: &[TestInf]) {
            for f in f {
                println!("function|{}", f.0);
                let inverse = f.1(modulus, unit);
                val(unit as u128, modulus as u128, inverse);
            }
        }

        #[test]
        fn basic_test() {
            let modulus = 17;
            let unit = 11;
            let f = [A, B, C, C2, D, D2];

            test(unit, modulus, &f);
        }

        #[test]
        fn basic_test2() {
            let modulus = 59;
            let unit = 58;
            let f = [B, C, C2, D, D2];

            test(unit, modulus, &f);
        }

        #[test]
        fn non_invertible() {
            let modulus = 33;
            let unit = 11;
            let f = [A, B, C, C2, D, D2];

            for f in f {
                println!("function|{}", f.0);
                let res = f.1(modulus, unit);
                assert_eq!(Ok(None), res);
            }
        }

        #[test]
        #[cfg(feature = "extended-tests")]
        fn val_rels() {
            let modulus = 59;
            let unit = 60;
            let f = [A, B, C, C2, D, D2];
            for f in f {
                println!("function|{}", f.0);
                let res = f.1(modulus, unit);
                assert_eq!(true, res.is_err());
            }
        }

        #[test]
        fn coprime_primes() {
            let modulus = 1_299_709;
            let unit = 56_999;
            let f = [A, B, C, C2, D, D2];

            test(unit, modulus, &f);
        }

        #[test]
        fn coprime_odd() {
            let modulus = 2_559_031_471; // 150531263ᵖ ⋅17ᵖ
            let unit = 1_956_912_061; // 150531697ᵖ ⋅13ᵖ
            let f = [A, C, C2, D, D2];

            test(unit, modulus, &f);
        }

        #[test]
        fn coprime_mixed_a() {
            let modulus = 52_535_230_703; // 150530747ᵖ ⋅349ᵖ
            let unit = 46_664_522_890; // 150530719ᵖ ⋅310ᶜ
            let f = [A, C, C2, D, D2];

            test(unit, modulus, &f);
        }

        #[test]
        fn coprime_mixed_b() {
            let modulus = 19_209_934_347; // 56666473ᵖ ⋅113ᵖ ⋅3ᵖ
            let unit = 10_993_312_058; // 56666557ᵖ ⋅2ᵖ ⋅97ᵖ
            let f = [A, C, C2, D, D2];

            test(unit, modulus, &f);
        }
    }

    mod sub_halve_mod {
        use crate::sub_halve_mod;

        #[test]
        fn correction_negative_number() {
            let res = sub_halve_mod(3, 7, 12);
            assert_eq!(4, res);
        }

        #[test]
        fn correction_odd_halving() {
            let res = sub_halve_mod(8, 1, 3);
            assert_eq!(5, res);
        }

        #[test]
        fn halving() {
            let res = sub_halve_mod(13, 3, usize::MAX);
            assert_eq!(5, res);
        }
    }
}

// cargo fmt & cargo test --features extended-tests
