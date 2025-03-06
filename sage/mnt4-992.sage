# MNT4_992
k = 4
u = 0xc85f1924b404f160077c049739e871907e407900a6d59abd8e25f63eaec03b9f974bfa92dd5cc38cb09ffdcdd3d19ab23bad8e228130bedd0e0859c32774
D = 95718723
c = 1
a = -3
b = 0x1c517e4f1632c3879c949ad49791cd8d9a6fc4bba403a3e69053e909ecbd42dbafb80100e0e46c53e85a07d869637e75ca6c3371de1bb090517613b30782f01489562fe913cd858d6671ab2a9d7ceb5c3ce8ea426c577ebb16b3c87a501e4bef0df989d7ec1037f32e852df87dce54696805764e1a072268d650a1ea
pnbits = 992
rnbits = 992
p = u**2 + u + 1
r = u**2 + 1
t = u + 1
y = 0x914c26752f9ff10b6a9060342744c1971bd11ee0a5a9f3dbff85fd462a6807cdbe69fef2c1fc2731306a180a5490999b3e810db101c100f78eee893a3
assert t**2-4*p == -D*y**2
assert p+1-t == r*c
(p-1) % (2**2 * 3 * 5**5 * 13) == 0
(r-1) % (2**4 * 3**2 * 5**10) == 0
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
