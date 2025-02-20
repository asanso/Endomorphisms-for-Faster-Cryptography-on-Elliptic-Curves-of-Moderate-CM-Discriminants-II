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

#t =E0.trace_of_frobenius()
#disc = t^2 - 4*p
#print(factor(disc))
#print()
#list of curves on the crater
#FpX = Fp['X']
#X = FpX.gen()
#D = fundamental_discriminant(-D)
#H = FpX(hilbert_class_polynomial(D))
#print(H.roots())
#print()