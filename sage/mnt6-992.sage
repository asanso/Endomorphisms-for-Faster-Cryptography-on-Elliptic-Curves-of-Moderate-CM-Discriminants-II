def multi_scalar_mul(P, k1, endo, k2):
    return k1*P + k2*endo(P)

def fast_scalar_mul(n,P):
    beta = vector([n,0])*N_inv
    b = vector([int(beta[0]), int(beta[1])]) * N
    k1 = n-b[0]
    k2 = -b[1]
    return  multi_scalar_mul(P,k1, full_end, k2)

def projective_maps_optimized_simple(phi,Fp, neg,k):
    rX,sXY = phi
    Fpx = Fp['x']
    x = Fpx.gen()
    FpX = Fpx.fraction_field()
    X = FpX.gen()
    Fpxz = Fp['x', 'z']
    FpXZ = FractionField(Fpxz)
    X, Z = FpXZ.gens()
    
    psi1 = rX.numerator()
    if neg:
        psi3 = -rX.denominator().sqrt()
    else:
        psi3 = rX.denominator().sqrt()

    sX = sXY(y=1)
    assert psi3^3 == sX.denominator()
    psi2 = sX.numerator()
    psi1XZ = psi1(x=X/Z)
    psi2XZ = psi2(x=X/Z)
    psi3XZ = psi3(x=X/Z)

    a = psi1XZ*psi3XZ *Z^(k+1)
    b = psi2XZ *Z^k
    c = psi3XZ^3 *Z^k
    return a,b,c

def projective_maps_optimized(phi,Fp, neg):
    rX,sXY = phi
    Fpx = Fp['x']
    x = Fpx.gen()
    FpX = Fpx.fraction_field()
    X = FpX.gen()
    Fpxz = Fp['x', 'z']
    FpXZ = FractionField(Fpxz)
    X, Z = FpXZ.gens()
    
    psi1 = rX.numerator()
    if neg:
        psi3 = -rX.denominator().sqrt()
    else:
        psi3 = rX.denominator().sqrt()

    sX = sXY(y=1)

    assert psi3^3 == (sX.denominator()/10314424798490535546171949056)
    psi2 = sX.numerator()
    psi1XZ = psi1(x=X/Z)
    psi2XZ = psi2(x=X/Z)
    psi3XZ = psi3(x=X/Z)

    a = psi1XZ*psi3XZ *Z^19
    b = psi2XZ *Z^18
    c = psi3XZ^3 *Z^18
    return a,b,c

def end_composition_optimized(P):
    x0 = P[0]
    y0 = P[1]
    z0 = 1 

    #1st isogeny
    x1 = a0(x0,z0)  
    y1 = y0 *b0(x0,z0)
    z1 = z0*c0(x0,z0)
    #2nd isogeny
    x2 = a1(x1,z1)  
    y2 = y1 *b1(x1,z1)
    z2 = z1*c1(x1,z1)
    #3rd isogeny
    x3 = a2(x2,z2)  
    y3 = y2 *b2(x2,z2)
    z3 = z2*c2(x2,z2)
    #4th isogeny
    x4 = a3(x3,z3)  
    y4 = y3 *b3(x3,z3)
    z4 = z3*c3(x3,z3)
    #5th isogeny
    x5 = a4(x4,z4)  
    y5 = y4 *b4(x4,z4)
    z5 = z4*c4(x4,z4)
    #6th isogeny
    x6 = a5(x5,z5)  
    y6 = y5 *b5(x5,z5)
    z6 = z5*c5(x5,z5)
    #7th isogeny
    x7 = a6(x6,z6)  
    y7 = y6 *b6(x6,z6)
    z7 = z6*c6(x6,z6)
    return  isoX(x7, y7), isoY(x7,y7), z7

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

phi0 =  E0.isogenies_prime_degree(3)[0]
E1 = phi0.codomain()
phi1 = E1.isogenies_prime_degree(13)[0]
E2 = phi1.codomain()
phi2 = E2.isogenies_prime_degree(13)[0]
E3 = phi2.codomain()
phi3 = E3.isogenies_prime_degree(13)[0]
E4 = phi3.codomain()
phi4 = E4.isogenies_prime_degree(17)[10]
E5 = phi4.codomain()
phi5 = E5.isogenies_prime_degree(23)[1]
E6 = phi5.codomain()
phi6 = E6.isogenies_prime_degree(41)[0]
assert phi6.codomain().j_invariant() == E0.j_invariant()

# Computing eigenvalue 
end =(phi6*phi5*phi4*phi3*phi2*phi1*phi0)
iso =end.codomain().isomorphism_to(E0)
full_end = (iso*end)
trace = full_end.trace()

norm = 3*13^3*17*23*41
Fr = GF(r)
R.<x> = PolynomialRing(Fr)
poly = x^2 - trace*x + norm
roots = poly.roots()
P = E0.random_point()
Q = full_end(P)
eigen = roots[1][0]

assert Q == eigen*P

# GLV

M = Matrix([[int(-eigen),1], [int(r),0]])
#print(M)
N = M.LLL()
N_inv = N**-1

n = ZZ.random_element(r)
S1 = n*P
S2 = fast_scalar_mul(n,P)
assert S1 == S2

a0,b0,c0 = projective_maps_optimized_simple(phi0,Fp, True, 3)
a1,b1,c1 = projective_maps_optimized(phi1,Fp, False)
a2,b2,c2 = projective_maps_optimized(phi2,Fp, True)
a3,b3,c3 = projective_maps_optimized(phi3,Fp, True)
a4,b4,c4 = projective_maps_optimized_simple(phi4,Fp, True, 24)
a5,b5,c5 = projective_maps_optimized_simple(phi5,Fp, True, 33)
a6,b6,c6 = projective_maps_optimized_simple(phi6,Fp, False, 60)

#isoX = iso.rational_maps()[0]
#isoY = iso.rational_maps()[1]


