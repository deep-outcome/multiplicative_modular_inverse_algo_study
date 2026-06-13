# Multiplicative Modular Inverse Algorithm Study

## Multiplicative Modular Inverse

If `gcd(left,right) =1`, then linear combination, linear Diophantine equation, 
`l ⋅left +r ⋅right =1` expresses respective modular multiplicative inverses:

⏵ `l ⋅left  ≡1 (mod right)`, `l ≡left⁻¹ (mod right)`<br/>
⏵ `r ⋅right ≡1 (mod left)`,  `r ≡right⁻¹ (mod left)`

Inverses are called Bézout's coefficients and equation Bézout's identity.

Note: Bézout’s coefficients and Bézout’s identity are not limited to `gcd(left,right) =1`, yet both apply to all gcds.

## Algorithm in Pseudocode 

**Requires:** Modulus `m ≥3, m mod 2 =1` and unit `0 ≤u <m`. </br>
**Ensures:** `¹⁄ᵤ mod m` if `gcd(u,m) =1`, or `0`

```
a 🠈u, g 🠈1, b 🠈m, h 🠈0
while a ≠0 do
  if a =0 mod 2 then
    a 🠈a /2
    g 🠈g /2 mod m
 else
   if a <b then
     (a,g,b,h) 🠈(b,h,a,g)
   a 🠈(a −b) /2
   g 🠈(g −h) /2 mod m
if b =1 then
  return h
return 0
```

Notes
- Loop invariants: `a =g ⋅u mod m` and `b =h ⋅u mod m`.

## Simple Performance Comparison

Statistically inaccurate, bound to one computation scenario simple performance overview, based on one run of performance comparison of various algorithm implementations outlined in table below. However, 
results are stable accross runs.

```rust
#[bench]
fn method(b: &mut Bencher) {
    let modulus = 2_559_031_471;
    let unit = 1_956_912_061;

    b.iter(|| _ = method(modulus, unit));
}
```


|     Method     | Description                                           |   Mean   | Deviation |
|----------------|-------------------------------------------------------|----------|-----------|
| mod_inverse_A  | Basic, avoids modular corrections                     | 61.25 ns | ± 0.45    |
| mod_inverse_A2 | CTZ, avoids modular corrections                       | 30.84 ns | ± 0.25    |
| mod_inverse_C  | CTZ, modular corrections                              | 35.78 ns | ± 0.45    |
| mod_inverse_C2 | CTZ, modular corrections, early escape                | 39.52 ns | ± 1.10    |
| mod_inverse_D  | CTZ, modular corrections, no swap                     | 60.48 ns | ± 2.16    |
| mod_inverse_D2 | CTZ, modular corrections, no swap, early escape       | 52.91 ns | ± 0.66    |


<small>
Configuration:
<ul>
<li>Operating System: openSUSE Leap 16.0</li>
<li>Kernel Version: 6.12.0-160000.33-default (64-bit)</li>
<li>Processor: 16 × AMD Ryzen 7 3800X 8-Core Processor</li>
<li>RAM: DDR4 2133 MT/s</li>
</ul>
</small>

See [Stein's Greatest Common Divisor Algorithm Study](https://github.com/deep-outcome/stein_gcd_algo_study#simple-performance-comparison) to compare with other Stein's Extended Greatest Common Divisor algorithm performance.

More or less interesting times are gotten when `fat` [Link Time Optimizations](https://doc.rust-lang.org/cargo/reference/profiles.html#lto) is used.


|     Method     |     Mean      | Deviation |
|----------------|---------------|-----------|
| mod_inverse_A  | 62.33 ns/iter | ± 1.38    |
| mod_inverse_A2 | 30.83 ns/iter | ± 1.01    |
| mod_inverse_C  | 35.75 ns/iter | ± 0.41    |
| mod_inverse_C2 | 39.54 ns/iter | ± 0.29    |
| mod_inverse_D  | 29.96 ns/iter | ± 5.06    |
| mod_inverse_D2 | 24.85 ns/iter | ± 7.32    |


## References

For further scholastic reference on topic of Greates Common Divisor and Modular Inverse computations, especially with focus on optimization, see:
- [Optimized Binary GCD for Modular Inversion](./optimized-binary-gcd-for-modular-inversion.pdf) by Thomas Pornin
- [Fast constant-time
gcd computation and modular inversion](./fast-constant-time-gcd-computation-and-modular-inversion.pdf) by Daniel J. Bernstein and Bo-Yin Yang
