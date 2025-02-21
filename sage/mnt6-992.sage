# MNT6_992
k = 6
u = -0x642f8c925a0278b003be024b9cf438c83f203c80536acd5ec712fb1f57601dcfcba5fd496eae61c6584ffee6e9e8cd591dd6c71140985f6e87042ce193ba
D = 95718723
c = 1
a = -3
b = 0xa0c214d66abeed117834b3812f966d30b7ce0cef1dc8aec978c8da94cbcf67a2d99c2c1428f9ffb86b280c13144154fb1be8a9e1f4fd886271c3816e25f623c01e344d24440dd2873f6b207862dab186fde6c075b5a0dccdcd86c7656865d75ca94f63a091cf6fc5d538342b81b871bda4b63fac1a9ed8b36ccd906
pnbits = 992
rnbits = 992
p = 4*u**2 + 1
r = 4*u**2 - 2*u + 1
t = 2*u + 1
y = 0x914c26752f9ff10b6a9060342744c1971bd11ee0a5a9f3dbff85fd462a6807cdbe69fef2c1fc2731306a180a5490999b3e810db101c100f78eee893a3
assert t**2-4*p == -D*y**2
assert p+1-t == r*c
(p-1) % (2^4 * 3^2 * 5^10) == 0
(r-1) % (2^2 * 3 * 5^5 * 13) == 0
Fp = GF(p)
E0 = EllipticCurve([Fp(a),Fp(b)])

# there is only one isogeny here
phi0 =  E0.isogenies_prime_degree(3)[0]
E1 = phi0.codomain()