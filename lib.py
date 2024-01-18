import time
import timeit

def _ror(val: int, n: int, m: int):
    save = val & ((1 << n) - 1)
    val >>= n
    val |= save << (m - n)
    return val


def _rol(val: int, n: int, m: int):
    save = (val & (((1 << n) - 1) << (m - n))) >> (m - n)
    val <<= n
    val &= (1 << m) - 1
    val |= save
    return val


class GF2NBElement:
    m = 419
    mask = (1 << m) - 1
    mul_matrix = None

    def __init__(self, val, radix: int = 16):
        if isinstance(val, int):
            self.val = val
        if isinstance(val, str):
            self.val = int(val, radix)

        # cut
        self.val = self.val & (1 << self.m) - 1

    @classmethod
    def from_other(cls, gfe) -> "GF2NBElement":
        new = cls(0)
        new.val = gfe.val
        new.m = gfe.m
        return new

    @classmethod
    def ONE(cls) -> "GF2NBElement":
        return cls((1 << cls.m) - 1)

    @classmethod
    def ZERO(cls) -> "GF2NBElement":
        return cls(0)

    def assert_same_field(self, other: "GF2NBElement"):
        if self.m != other.m:
            raise ValueError("elements over different fields")

    def __str__(self) -> str:
        return f'{self.val:X}'

    def __repr__(self):
        return f'{self.val:X}'

    def __format__(self, format_spec):
        fs = f'{{:{format_spec}}}'
        return fs.format(self.val)

    def __eq__(self, other):
        return self.val == other.val

    def __xor__(self, other: "GF2NBElement"):
        self.assert_same_field(other)
        return self.val ^ other.val

    def __add__(self, other):
        return self.add(other)

    def __mul__(self, other):
        return self.mul(other)

    def sqr(self) -> "GF2NBElement":
        new = self.from_other(self)
        new.val = _ror(new.val, 1, new.m)
        return new

    def add(self, other: "GF2NBElement") -> "GF2NBElement":
        self.assert_same_field(other)
        new = self.from_other(self)
        new.val ^= other.val
        return new

    def trace(self):
        t = self.val.bit_count() % 2
        if t == 1:
            return self.ONE()
        elif t == 0:
            return self.ZERO()
        else:
            raise RuntimeError("implementation is incorrect")

    def mul(self, other) -> "GF2NBElement":
        self.assert_same_field(other)
        if self.mul_matrix is None:
            self.mul_matrix = self._gen_mm()

        m = self.mul_v_m_u(self.val, other.val)

        return self.__class__(m)

    def pow(self, e: int):
        if e > 2 ** self.m - 1:
            raise ValueError("exponent > 2**m -1")

        if self == GF2NBElement.ZERO():
            return GF2NBElement.ZERO()

        if e == 0:
            return GF2NBElement.ONE()

        c = GF2NBElement.ONE()
        for i in range(e.bit_length(), -1, -1):
            c = c.sqr()
            if e & (1 << i) != 0:
                c = c * self

        return c

    def inverse(self):
        return self.ito()

    def ito(self):
        b = self
        k = 1
        f = self.m - 1
        # for i: = t-1 -> 0
        for i in range(f.bit_length() - 2, -1, -1):
            t = b
            for j in range(k):
                t = t.sqr()
            b = t * b
            k *= 2
            if f & (1 << i) != 0:
                b = b.sqr() * self
                k += 1

        b = b.sqr()
        return b

    def _gen_mm(self):
        p = 2 * self.m + 1
        rows = list()
        for i in range(self.m):
            rows.append(list())
            for j in range(self.m):
                if (pow(2, i, p) + pow(2, j, p)) % p == 1 \
                        or (pow(2, i, p) - pow(2, j, p)) % p == 1 \
                        or (-pow(2, i, p) + pow(2, j, p)) % p == 1 \
                        or (-pow(2, i, p) - pow(2, j, p)) % p == 1:
                    r = 1
                else:
                    r = 0

                rows[i].append(r)

        return rows

    def generate_mm(self):
        if self.mul_matrix is None:
            self.mul_matrix = self._gen_mm()

    def mul_v_m_u(self, u: int, v: int):
        r = 0

        #0
        for _ in range(self.m):
            r = (r << 1) & self.mask
            d = 0

            #1
            for j in range(self.m):
                d <<= 1
                #2
                for i in range(self.m):
                    mat = self.mul_matrix[self.m - 1 - i][j]
                    # opt1 (works)
                    if mat == 0:
                        continue
                    d ^= (((u & (1 << i)) >> i) & mat)

            c = 0
            #3
            for i in range(self.m):
                # opt 2
                # c ^= ((d & (1 << i)) >> i) & ((v & (1 << i)) >> i)
                c ^= (d & v & (1 << i)) >> i

            r |= c
            u = _rol(u, 1, self.m)
            v = _rol(v, 1, self.m)

        return r
