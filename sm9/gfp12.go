package sm9

// For details of the algorithms used, see "Multiplication and Squaring on
// Pairing-Friendly Fields, Devegili et al.
// http://eprint.iacr.org/2006/471.pdf.

import (
	"math/big"
)

// gfP12 implements the field of size p¹² as a quadratic extension of gfP6
// where ω²=τ.
type gfP12 struct {
	x, y gfP6 // value is xω + y
}

var gfP12Gen *gfP12 = &gfP12{
	x: gfP6{
		x: gfP2{
			x: gfP{0xeb2aeaa2823d010c, 0xe192c39d7c3e6440, 0x68411e843fea2a9b, 0x5f23b1ce3ac438e7},
			y: gfP{0x065c1ad6d376db4f, 0xe2447d6d5edfdda6, 0x0d4eba5c8c017781, 0x61ebca2110d736bf},
		},
		y: gfP2{
			x: gfP{0xc219536a54552cae, 0xc4e4ad66027f8f55, 0xff31b23d5bc78184, 0x3b0fc03d5711c93d},
			y: gfP{0x290e1c8bdb9441aa, 0x074e1694c800c130, 0xfa196a2583564700, 0x254eb32dea84e64d},
		},
		z: gfP2{
			x: gfP{0x24fb5abe38626c9c, 0xd32d71f71d7bd3de, 0x671d686fd9c9271d, 0xa3eec3cd6a795be8},
			y: gfP{0x7b9c733c1f964b52, 0x9b988c0c238fb05e, 0xe546ccb8d6e1f9b8, 0xb101d668bfbf8ac8},
		},
	},
	y: gfP6{
		x: gfP2{
			x: gfP{0x487ab1a6229d91f3, 0x7e2a3e36c6c822c7, 0x282c24f00c10930f, 0x2efe33f18332bb77},
			y: gfP{0x346965f4dc5b5813, 0xed43ed38c0ce33e6, 0x9ba7630e295a5ce7, 0xa6db7142e0ca24ae},
		},
		y: gfP2{
			x: gfP{0xfea0bce10965b32b, 0x441e074b4573390c, 0xe9d6067a4cf3c571, 0x9ee43c7e3740bcd8},
			y: gfP{0x0e06727b47ee6118, 0xb01ab631f2f10a18, 0xb0ebd9852fc780ef, 0xaa07010f9d42787c},
		},
		z: gfP2{
			x: gfP{0xbe7381e2bce90a00, 0x2a72158dbf514e31, 0x44e199bee3498d4d, 0x6a5fed210720de58},
			y: gfP{0xb55d63ee8d7a8468, 0x9ef5d413e3176666, 0x796c802ec3f1370b, 0xa0f422c35d7b6262},
		},
	},
}

func (e *gfP12) String() string {
	return "(" + e.x.String() + "," + e.y.String() + ")"
}

func (e *gfP12) Set(a *gfP12) *gfP12 {
	e.x.Set(&a.x)
	e.y.Set(&a.y)
	return e
}

func (e *gfP12) SetZero() *gfP12 {
	e.x.SetZero()
	e.y.SetZero()
	return e
}

func (e *gfP12) SetOne() *gfP12 {
	e.x.SetZero()
	e.y.SetOne()
	return e
}

func (e *gfP12) IsZero() bool {
	return e.x.IsZero() && e.y.IsZero()
}

func (e *gfP12) IsOne() bool {
	return e.x.IsZero() && e.y.IsOne()
}

func (e *gfP12) Conjugate(a *gfP12) *gfP12 {
	e.x.Neg(&a.x)
	e.y.Set(&a.y)
	return e
}

func (e *gfP12) Neg(a *gfP12) *gfP12 {
	e.x.Neg(&a.x)
	e.y.Neg(&a.y)
	return e
}

// Frobenius computes (xω+y)^p = x^p ω·ξ^((p-1)/6) + y^p
func (e *gfP12) Frobenius(a *gfP12) *gfP12 {
	e.x.Frobenius(&a.x)
	e.y.Frobenius(&a.y)
	e.x.MulScalar(&e.x, xiToPMinus1Over6)
	return e
}

// FrobeniusP2 computes (xω+y)^p² = x^p² ω·ξ^((p²-1)/6) + y^p²
func (e *gfP12) FrobeniusP2(a *gfP12) *gfP12 {
	e.x.FrobeniusP2(&a.x)
	e.x.MulGFP(&e.x, xiToPSquaredMinus1Over6)
	e.y.FrobeniusP2(&a.y)
	return e
}

func (e *gfP12) FrobeniusP4(a *gfP12) *gfP12 {
	e.x.FrobeniusP4(&a.x)
	e.x.MulGFP(&e.x, xiToPSquaredMinus1Over3)
	e.y.FrobeniusP4(&a.y)
	return e
}

func (e *gfP12) Add(a, b *gfP12) *gfP12 {
	e.x.Add(&a.x, &b.x)
	e.y.Add(&a.y, &b.y)
	return e
}

func (e *gfP12) Sub(a, b *gfP12) *gfP12 {
	e.x.Sub(&a.x, &b.x)
	e.y.Sub(&a.y, &b.y)
	return e
}

func (e *gfP12) Mul1(a, b *gfP12) *gfP12 {
	tx := (&gfP6{}).Mul(&a.x, &b.y)
	t := (&gfP6{}).Mul(&b.x, &a.y)
	tx.Add(tx, t)

	ty := (&gfP6{}).Mul(&a.y, &b.y)
	t.Mul(&a.x, &b.x).MulTau(t)

	e.x.Set(tx)
	e.y.Add(ty, t)
	return e
}
func (e *gfP12) Mul(a, b *gfP12) *gfP12 { //lazy reduction;1113version;

	t0, t1, t2, t3 := &gfP6{}, &gfP6{}, &gfP6{}, &gfP6{}
	t0.Add(&a.x, &a.y)
	t1.Add(&b.x, &b.y)
	t2.Mul(&a.x, &b.x)
	t3.Mul(&a.y, &b.y)

	t0.Mul(t0, t1)
	t0.Sub(t0, t2)
	t0.Sub(t0, t3)

	t2.MulTau(t2)
	t3.Add(t3, t2)

	e.x.Set(t0)
	e.y.Set(t3)
	return e
}

func (e *gfP12) MulScalar(a *gfP12, b *gfP6) *gfP12 {
	e.x.Mul(&e.x, b)
	e.y.Mul(&e.y, b)
	return e
}

func (c *gfP12) Exp(a *gfP12, power *big.Int) *gfP12 {
	sum := (&gfP12{}).SetOne()
	t := &gfP12{}

	for i := power.BitLen() - 1; i >= 0; i-- {
		t.Square(sum)
		if power.Bit(i) != 0 {
			sum.Mul(t, a)
		} else {
			sum.Set(t)
		}
	}

	c.Set(sum)
	return c
}

func (e *gfP12) Square(a *gfP12) *gfP12 {
	// Complex squaring algorithm
	v0 := (&gfP6{}).Mul(&a.x, &a.y)

	t := (&gfP6{}).MulTau(&a.x)
	t.Add(&a.y, t)
	ty := (&gfP6{}).Add(&a.x, &a.y)
	ty.Mul(ty, t).Sub(ty, v0)
	t.MulTau(v0)
	ty.Sub(ty, t)

	e.x.Add(v0, v0)
	e.y.Set(ty)
	return e
}

func (e *gfP12) Square1(a *gfP12) *gfP12 { //not faster;1115
	v0 := (&gfP6{}).Square(&a.x)
	v1 := (&gfP6{}).Square(&a.y)
	v2 := (&gfP6{}).Add(&a.x, &a.y)
	v2.Square(v2)

	e.x.Sub(v2, v0)
	e.x.Sub(v2, v1)

	v0.MulTau(v0)
	v0.Add(v0, v1)
	e.y.Set(v0)

	return e
}

func (e *gfP12) Invert(a *gfP12) *gfP12 {
	// See "Implementing cryptographic pairings", M. Scott, section 3.2.
	// ftp://136.206.11.249/pub/crypto/pairings.pdf
	t1, t2 := &gfP6{}, &gfP6{}

	t1.Square(&a.x)
	t2.Square(&a.y)
	t1.MulTau(t1).Sub(t2, t1)
	t2.Invert(t1)

	e.x.Neg(&a.x)
	e.y.Set(&a.y)
	e.MulScalar(e, t2)
	return e
}
