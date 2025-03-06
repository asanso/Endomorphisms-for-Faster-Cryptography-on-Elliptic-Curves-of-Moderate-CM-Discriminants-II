def multi_scalar_mul(P, k1, endo, k2):
    return k1*P + k2*endo(P)

def fast_scalar_mul(n,P):
    beta = vector([n,0])*N_inv
    b = vector([int(beta[0]), int(beta[1])]) * N
    k1 = n-b[0]
    k2 = -b[1]
    return  multi_scalar_mul(P,k1, full_end, k2)


def projective_maps_optimized(phi,Fp, neg,n,k):
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
    assert psi3^3 == (sX.denominator()/n)
    psi2 = sX.numerator()/n
    psi1XZ = psi1(x=X/Z)
    psi2XZ = psi2(x=X/Z)
    psi3XZ = psi3(x=X/Z)

    a = psi1XZ*psi3XZ *Z^k
    b = psi2XZ *Z^(k-1)
    c = psi3XZ^3 *Z^(k-1)
    return a,b,c

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


x=589042076226215548287689476330582191609793749744861871604900539741454676076257591997245318667180479210398140701841229910422936321486475345296314; 
p=x^2-x+1
Fp = GF(p)
a=327109485262376726136172643805932252878964212240064991046803990546570692113955204504295853421320565971168014448160518214189382593729022922934619191479697131508401422054977766988346990626462981083698419529403033388974770759484439350971607246476622466769281101195803851961334206387433757192;
b=279050255428736347835787201657572365478299307911044416076073175912778177975625581983202449094627882457913470267629624709093466687587530683922179493820089554667794629254883792790600467080376496365458891552310274342079659991659178694655626276630944082371547597372703391340583858884507859049;

h=76182347360104295691389198837779736649246768246453788622413701157110735036917351906437908803769129014598039561668077363397378012546375286033150471993610
r=4554474620279221025376588013428792331658085875800439965576470444142136125074436360422816311318800204172284882642774195680976231192452677

E0 = EllipticCurve(GF(p), [a,b])
phi0 =  E0.isogenies_prime_degree(5)[0]    
E1 = phi0.codomain()
phi1 = E1.isogenies_prime_degree(5)[0]
E2 = phi1.codomain()
phi2 = E2.isogenies_prime_degree(5)[1]
E3 = phi2.codomain()
phi3 = E3.isogenies_prime_degree(5)[0]
E4 = phi3.codomain()
phi4 = E4.isogenies_prime_degree(5)[1]
E5 = phi4.codomain()
phi5 = E5.isogenies_prime_degree(5)[0]
E6 = phi5.codomain()        
phi6 = E6.isogenies_prime_degree(5)[0]
E7 = phi6.codomain()
phi7 = E7.isogenies_prime_degree(7)[0]
E8 = phi7.codomain()
phi8 = E8.isogenies_prime_degree(7)[0]
E9 = phi8.codomain()
phi9 = E9.isogenies_prime_degree(23)[1]
assert phi9.codomain().j_invariant() == E0.j_invariant()

# Computing eigenvalue 
end =(phi9*phi8*phi7*phi6*phi5*phi4*phi3*phi2*phi1*phi0)
iso =end.codomain().isomorphism_to(E0)
full_end = (iso*end)
trace = full_end.trace()

norm = 5^7*7^2*23
Fr = GF(r)
R.<x> = PolynomialRing(Fr)
poly = x^2 - trace*x + norm
roots = poly.roots()
P = h*E0.random_point()
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

a0,b0,c0 = projective_maps_optimized(phi0,Fp, True, 2176782336, 7)
a1,b1,c1 = projective_maps_optimized(phi1,Fp, False, 2176782336, 7)
a2,b2,c2 = projective_maps_optimized(phi2,Fp, False, 2176782336, 7)
a3,b3,c3 = projective_maps_optimized(phi3,Fp, False, 2176782336, 7)
a4,b4,c4 = projective_maps_optimized(phi4,Fp, False, 2176782336, 7)
a5,b5,c5 = projective_maps_optimized(phi5,Fp, False, 2176782336, 7)
a6,b6,c6 = projective_maps_optimized(phi6,Fp, False, 2176782336, 7)
a7,b7,c7 = projective_maps_optimized(phi7,Fp, True, 101559956668416, 10)
a8,b8,c8 = projective_maps_optimized(phi8,Fp, True, 101559956668416, 10)
q9,b9,c9 = projective_maps_optimized_simple(phi9,Fp, False, 33)
