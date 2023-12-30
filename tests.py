import unittest
from lib import *


class AllTests(unittest.TestCase):
    a = GF2NBElement(
        "02416715E64D7DF197CD3699B8FC37CD1813E5FEB8C05641FCFD20B9DFC9CF479111BFD12B47FFFC5ADD656F1C997A1119EE6523A6")
    b = GF2NBElement(
        "060BA7A69C89F783A2907C4641D05F92D777D085946D5961707C682A1E0D4D91484FDDF58F14663CF69183492A63C28074CBBC22B0")

    c = GF2NBElement(
        "060FD6E96E43D932117D566E6E369187993D64D674E6F242788E8BA326F2743130266678F4419987DC77F2B53DA63DB40581978B05"
    )
    n = 0x123

    def test_sqr(self):
        r = GF2NBElement(
            "0120B38AF326BEF8CBE69B4CDC7E1BE68C09F2FF5C602B20FE7E905CEFE4E7A3C888DFE895A3FFFE2D6EB2B78E4CBD088CF73291D3")

        self.assertEqual(self.a.sqr(), r)

    def test_add(self):
        r = GF2NBElement(
            "044AC0B37AC48A72355D4ADFF92C685FCF64357B2CAD0F208C814893C1C482D6D95E6224A45399C0AC4CE62636FAB8916D25D90116")

        self.assertEqual(self.a + self.b, r)

    def test_mul(self):
        r = GF2NBElement(
            "07791F2710C184BF5A9FCD4C8FEA582E390E95447D943C06488B716AD82965BD0C75C50EC0D747BB30F173DF4C0F7AA116FCB898C0")

        f = self.a * self.b
        print(f"{f:x}")

        self.assertEqual(r, f)

    def test_trace(self):
        r = GF2NBElement.ZERO()

        self.assertEqual(self.a.trace(), r)

    def test_inverse(self):
        r = GF2NBElement(
            "01F1EBC0F68A3874AF25977E3E1B9524D31EBCCAB45D9FEA4194B3120E760B3FEB3064334321D605463F956A053BBC97D44D9A3BC6")

        self.assertEqual(r, self.a.inverse())

    def test_pow(self):
        r = GF2NBElement(
            "0676107F7754906E3B58C11FA81ECA479AC84C8FD860EFD299FFBD8ECE8328BD0F0AF14DF4DCA2648991682603B5449C2A60F810E4")

        self.assertEqual(r, self.a.pow(self.n))


    def test_eq1(self):
        self.assertEqual(
            (self.a + self.b) * self.c, (self.b * self.c + self.c * self.a)
        )


if __name__ == '__main__':
    unittest.main()