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
