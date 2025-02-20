from sage.all import ZZ, FiniteField

class WeierstrassCurve():
    def __init__(self, p, a, b, r, cofactor):
        self.p = p
        self.Fp = FiniteField(p)
        self.a = self.Fp(a)
        self.b = self.Fp(b)
        self.r = r
        self.cofactor = cofactor

    def __str__(self):
        a = ZZ(self.a)
        b = ZZ(self.b)
        p = ZZ(self.p)
        return f"Elliptic curve in Weierstrass form y² = x³ + {a}x + {b} over GF({p})"

    def random_point(self):
        x = self.Fp.random_element()
        rhs = x**3 + self.a * x + self.b
        while not rhs.is_square():
            x = self.Fp.random_element()
            rhs = x**3 + self.a * x + self.b
        y = rhs.sqrt()
        z = self.Fp.random_element()  # Random projective coordinate
        x *= z
        y *= z
        return PointWeierstrass(x, y, z, self)

    def point_of_order_r(self):
        P = self.random_point().clear_cofactor()
        while P.is_zero():
            P = self.random_point().clear_cofactor()
        assert P.scalar_mul(self.r).is_zero()
        return P


class PointWeierstrass():
    def __init__(self, X, Y, Z, curve):
        self.X = X  # Projective coordinate X
        self.Y = Y  # Projective coordinate Y
        self.Z = Z  # Projective coordinate Z
        self.curve = curve

    def __str__(self):
        if self.is_zero():
            return "Point at infinity"
        return f"x:{self.X / self.Z}\ny:{self.Y / self.Z}\nz:{self.Z}"

    def __eq__(self, other):
        if self.is_zero() and other.is_zero():
            return True
        if self.is_zero() or other.is_zero():
            return False
        return self.X * other.Z == other.X * self.Z and self.Y * other.Z == other.Y * self.Z

    def is_zero(self):
        return self.Z.is_zero()

    def neg(self):
        return PointWeierstrass(self.X, -self.Y, self.Z, self.curve)

    def on_curve(self):
        X, Y, Z = self.X, self.Y, self.Z
        a, b = self.curve.a, self.curve.b
        return Y**2 * Z == X**3 + a * X * Z**2 + b * Z**3

    def add(self, other):
        if self.is_zero():
            return other
        if other.is_zero():
            return self
        if self.X * other.Z == other.X * self.Z and self.Y * other.Z == -other.Y * self.Z:
            return PointWeierstrass(self.curve.Fp(0), self.curve.Fp(1), self.curve.Fp(0), self.curve)  # Point at infinity

        X1, Y1, Z1 = self.X, self.Y, self.Z
        X2, Y2, Z2 = other.X, other.Y, other.Z

        # Projective point addition formulas
        U1 = Y2 * Z1
        U2 = Y1 * Z2
        V1 = X2 * Z1
        V2 = X1 * Z2

        if V1 == V2:
            if U1 == U2:
                return self.double()  # Same point
            else:
                return PointWeierstrass(self.curve.Fp(0), self.curve.Fp(1), self.curve.Fp(0), self.curve)  # Point at infinity

        U = U1 - U2
        V = V1 - V2
        W = Z1 * Z2
        A = U**2 * W - V**3 - 2 * V**2 * V2
        X3 = V * A
        Y3 = U * (V**2 * V2 - A) - V**3 * U2
        Z3 = V**3 * W

        return PointWeierstrass(X3, Y3, Z3, self.curve)

    def double(self):
        if self.is_zero():
            return self

        X, Y, Z = self.X, self.Y, self.Z
        a = self.curve.a

        # Projective point doubling formulas
        W = 3 * X**2 + a * Z**2
        S = Y * Z
        B = X * Y * S
        H = W**2 - 8 * B
        X3 = 2 * H * S
        Y3 = W * (4 * B - H) - 8 * Y**2 * S**2
        Z3 = 8 * S**3

        return PointWeierstrass(X3, Y3, Z3, self.curve)

    def scalar_mul(self, n):
        """
        Perform scalar multiplication using the double-and-add method.
        Handles negative scalars by negating the point.
        """
        if n == 0 or self.is_zero():
            # Multiplying by 0 or starting with the point at infinity
            return PointWeierstrass(self.curve.Fp(0), self.curve.Fp(1), self.curve.Fp(0), self.curve)  # Point at infinity
        
        if n < 0:
            # Handle negative scalars by negating the point
            n = -n
            P = self.neg()
        else:
            P = self

        # Initialize result as the point at infinity
        R = PointWeierstrass(self.curve.Fp(0), self.curve.Fp(1), self.curve.Fp(0), self.curve)

        # Perform the double-and-add algorithm
        while n > 0:
            if n & 1:  # If the current bit is 1
                R = R.add(P)
            P = P.double()  # Always double the point
            n >>= 1  # Shift to the next bit

        return R


    def multi_scalar_mul(self, k1, other, k2):
        P = self
        if k1 < 0:
            k1 = -k1
            P = P.neg()
        if k2 < 0:
            k2 = -k2
            other = other.neg()
        PplusOther = P.add(other)
        bits_k1 = ZZ(k1).bits()
        bits_k2 = ZZ(k2).bits()
        while len(bits_k1) < len(bits_k2):
            bits_k1.append(0)
        while len(bits_k2) < len(bits_k1):
            bits_k2.append(0)
        R = PointWeierstrass(self.curve.Fp(0), self.curve.Fp(1), self.curve.Fp(0), self.curve)  # Point at infinity
        for i in range(len(bits_k1) - 1, -1, -1):
            R = R.double()
            if bits_k1[i] == 1 and bits_k2[i] == 0:
                R = R.add(self)
            if bits_k1[i] == 0 and bits_k2[i] == 1:
                R = R.add(other)
            if bits_k1[i] == 1 and bits_k2[i] == 1:
                R = R.add(PplusOther)
        return R

    def clear_cofactor(self):
        return self.scalar_mul(self.curve.cofactor)
