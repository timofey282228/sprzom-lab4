from lib import *

a = GF2NBElement(input("a: "))
b = GF2NBElement(input("b: "))
n = int(input("n: "), 16)

a_add_b = a + b
a_mul_b = a * b
a_sqr = a.sqr()
a_inv = a.inverse()
a_trace = a.trace()

print(
f"""
a = {a:X} ({a:b})
a = {b:X} ({b:b})
a + b = {a_add_b:X} ({a_add_b:b})
a * b = {a_mul_b:X} ({a_mul_b:b})
a^2 = {a_sqr:X} ({a_sqr:b})
a_inv = {a_inv:X} ({a_inv:b})
a_trace = {a_trace:X} ({a_trace:b})
"""
)

a_pow_n = a.pow(n)
print(f"a^n = {a_pow_n:x} ({a_pow_n:b})")