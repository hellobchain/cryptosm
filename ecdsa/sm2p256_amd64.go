//go:build amd64
// +build amd64

package ecdsa

import (
	"bytes"
	"crypto/elliptic"
	"encoding/binary"
	"errors"
	"github.com/wsw365904/cryptosm/sm3"
	"math/big"
	"sync"
)

type (
	sm2Curve struct {
		*elliptic.CurveParams
	}

	sm2Point struct {
		xyz [12]uint64
	}
)

var (
	sm2P256             sm2Curve
	sm2Precomputed      *[43][32 * 8]uint64
	sm2precomputeOnce   sync.Once
	defaultZaBeforeByte []byte
	zaBeforeByte        []byte
)

const (
	SM2CurveName = "SM2-P-256"
)

func initSM2() {
	sm2P256.CurveParams = &elliptic.CurveParams{Name: SM2CurveName}
	sm2P256.P, _ = new(big.Int).SetString("FFFFFFFEFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF00000000FFFFFFFFFFFFFFFF", 16)
	sm2P256.N, _ = new(big.Int).SetString("FFFFFFFEFFFFFFFFFFFFFFFFFFFFFFFF7203DF6B21C6052B53BBF40939D54123", 16)
	sm2P256.B, _ = new(big.Int).SetString("28E9FA9E9D9F5E344D5A9E4BCF6509A7F39789F515AB8F92DDBCBD414D940E93", 16)
	sm2P256.Gx, _ = new(big.Int).SetString("32C4AE2C1F1981195F9904466A39C9948FE30BBFF2660BE1715A4589334C74C7", 16)
	sm2P256.Gy, _ = new(big.Int).SetString("BC3736A2F4F6779C59BDCEE36B692153D0A9877CC62A474002DF32E52139F0A0", 16)
	sm2P256.BitSize = 256
	zaBeforeByte = getZBefore(nil)
	defaultZaBeforeByte = getZBefore([]byte(defaultUid))
	return
}

func big2Bytes(big *big.Int) []byte {
	r := make([]byte, 32)
	bigBytes := big.Bytes()
	copy(r[32-len(bigBytes):], bigBytes)
	return r
}

const defaultUid = "1234567812345678"

func msgHash(za, msg []byte) *big.Int {
	eHashed := make([][]byte, 2)
	eHashed[0] = za
	eHashed[1] = msg
	return new(big.Int).SetBytes(sm3.Sm3Sum(bytes.Join(eHashed, nil))[:32])
}

func getE(pub *PublicKey, uid []byte, msg []byte) (*big.Int, error) {
	za, err := getZ(pub, uid)
	if err != nil {
		return nil, err
	}
	return msgHash(za, msg), nil
}

func getZBefore(uidValue []byte) []byte {
	uidValueLen := len(uidValue)
	var entl []byte
	var zHashedLen int
	if uidValueLen != 0 {
		zHashedLen = 6
		entl = make([]byte, 2)
		binary.BigEndian.PutUint16(entl, uint16(uidValueLen*8))
	} else {
		zHashedLen = 4
		entl = nil
	}
	a := big2Bytes(new(big.Int).Sub(sm2P256.Params().P, new(big.Int).SetInt64(3))) // a = p-3
	b := big2Bytes(sm2P256.Params().B)
	xG := big2Bytes(sm2P256.Params().Gx)
	yG := big2Bytes(sm2P256.Params().Gy)
	zHashed := make([][]byte, zHashedLen)
	if zHashedLen == 4 {
		zHashed[0] = a
		zHashed[1] = b
		zHashed[2] = xG
		zHashed[3] = yG
	} else {
		zHashed[0] = entl
		zHashed[1] = uidValue
		zHashed[2] = a
		zHashed[3] = b
		zHashed[4] = xG
		zHashed[5] = yG
	}
	return bytes.Join(zHashed, nil)
}

// Z = H256(ENTLA || IDA || a || b || xG || yG || xA || yA)
func getZ(pub *PublicKey, uid []byte) ([]byte, error) {
	uidLen := len(uid)
	var zBeforeByte []byte
	if uidLen == 0 {
		zBeforeByte = defaultZaBeforeByte
	} else if uidLen >= 8192 {
		return []byte{}, errors.New("SM2: uid too large")
	} else {
		entl := make([]byte, 2)
		binary.BigEndian.PutUint16(entl, uint16(uidLen*8))
		zBeforeByteTmp := make([][]byte, 3)
		zBeforeByteTmp[0] = entl
		zBeforeByteTmp[1] = uid
		zBeforeByteTmp[2] = zaBeforeByte
		zBeforeByte = bytes.Join(zBeforeByteTmp, nil)
	}
	x := big2Bytes(pub.X)
	y := big2Bytes(pub.Y)
	zHashed := make([][]byte, 3)
	zHashed[0] = zBeforeByte
	zHashed[1] = x
	zHashed[2] = y
	return sm3.Sm3Sum(bytes.Join(zHashed, nil)), nil
}

func (curve sm2Curve) Params() *elliptic.CurveParams {
	return curve.CurveParams
}

// Functions implemented in sm2_asm_amd64.s
// Montgomery multiplication modulo sm2
//go:noescape
func sm2Mul(res, in1, in2 []uint64)

// Montgomery square modulo sm2, repeated n times (n >= 1)
//go:noescape
func sm2Sqr(res, in []uint64, n int)

// Montgomery multiplication by 1
//go:noescape
func sm2FromMont(res, in []uint64)

// iff cond == 1  val <- -val
//go:noescape
func sm2NegCond(val []uint64, cond int)

// if cond == 0 res <- b; else res <- a
//go:noescape
func sm2MovCond(res, a, b []uint64, cond int)

// Endianness swap
//go:noescape
func sm2BigToLittle(res []uint64, in []byte)

//go:noescape
func sm2LittleToBig(res []byte, in []uint64)

// Constant time table access
//go:noescape
func sm2Select(point, table []uint64, idx int)

//go:noescape
func sm2SelectBase(point, table []uint64, idx int)

// Montgomery multiplication modulo Ord(G)
//go:noescape
func sm2OrdMul(res, in1, in2 []uint64)

// Montgomery square modulo Ord(G), repeated n times
//go:noescape
func sm2OrdSqr(res, in []uint64, n int)

// Point add with in2 being affine point
// If sign == 1 -> in2 = -in2
// If sel == 0 -> res = in1
// if zero == 0 -> res = in2
//go:noescape
func sm2PointAddAffineAsm(res, in1, in2 []uint64, sign, sel, zero int)

// Point add. Returns one if the two input points were equal and zero
// otherwise. (Note that, due to the way that the equations work out, some
// representations of ∞ are considered equal to everything by this function.)
//go:noescape
func sm2PointAddAsm(res, in1, in2 []uint64) int

// Point double
//go:noescape
func sm2PointDoubleAsm(res, in []uint64)

func (curve sm2Curve) Inverse(k *big.Int) *big.Int {
	if k.Sign() < 0 {
		// This should never happen.
		k = new(big.Int).Neg(k)
	}

	if k.Cmp(sm2P256.N) >= 0 {
		// This should never happen.
		k = new(big.Int).Mod(k, sm2P256.N)
	}

	// table will store precomputed powers of x. The four words at index
	// 4×i store x^(i+1).
	var table [4 * 15]uint64

	x := make([]uint64, 4)
	fromBig(x[:], k)
	// This code operates in the Montgomery domain where R = 2^256 mod n
	// and n is the order of the scalar field. (See initSM2 for the
	// value.) Elements in the Montgomery domain take the form a×R and
	// multiplication of x and y in the calculates (x × y × R^-1) mod n. RR
	// is R×R mod n thus the Montgomery multiplication x and RR gives x×R,
	// i.e. converts x into the Montgomery domain.
	RR := []uint64{0x901192af7c114f20, 0x3464504ade6fa2fa, 0x620fc84c3affe0d4, 0x1eb5e412a22b3d3b}
	sm2OrdMul(table[:4], x, RR)

	// Prepare the table, no need in constant time access, because the
	// power is not a secret. (Entry 0 is never used.)
	for i := 2; i < 16; i += 2 {
		sm2OrdSqr(table[4*(i-1):], table[4*((i/2)-1):], 1)
		sm2OrdMul(table[4*i:], table[4*(i-1):], table[:4])
	}

	x[0] = table[4*14+0] // f
	x[1] = table[4*14+1]
	x[2] = table[4*14+2]
	x[3] = table[4*14+3]

	sm2OrdSqr(x, x, 4)
	sm2OrdMul(x, x, table[4*14:4*14+4]) // ff
	t := make([]uint64, 4, 4)
	t[0] = x[0]
	t[1] = x[1]
	t[2] = x[2]
	t[3] = x[3]

	sm2OrdSqr(x, x, 8)
	sm2OrdMul(x, x, t) // ffff

	sm2OrdSqr(x, x, 8)
	sm2OrdMul(x, x, t) // ffffff
	t[0] = x[0]
	t[1] = x[1]
	t[2] = x[2]
	t[3] = x[3]

	sm2OrdSqr(x, x, 4)
	sm2OrdMul(x, x, table[4*14:4*14+4]) // fffffff

	sm2OrdSqr(x, x, 4)
	sm2OrdMul(x, x, table[4*13:4*13+4]) // fffffffe

	for i := 0; i < 4; i++ {
		sm2OrdSqr(x, x, 24)
		sm2OrdMul(x, x, t) // fffffffe_ffffff fffffffe_ffffffff_ffff fffffffe_ffffffff_ffffffff_ff fffffffe_ffffffff_ffffffff_ffffffff
	}

	// Remaining 32 windows
	expLo := [32]byte{0x7, 0x2, 0x0, 0x3, 0xd, 0xf, 0x6, 0xb, 0x2, 0x1, 0xc, 0x6, 0x0, 0x5, 0x2, 0xb, 0x5, 0x3, 0xb, 0xb, 0xf, 0x4, 0x0, 0x9, 0x3, 0x9, 0xd, 0x5, 0x4, 0x1, 0x2, 0x1}
	for i := 0; i < 32; i++ {
		sm2OrdSqr(x, x, 4)
		if expLo[i] != 0 {
			sm2OrdMul(x, x, table[4*(expLo[i]-1):])
		}
	}

	// Multiplying by one in the Montgomery domain converts a Montgomery
	// value out of the domain.
	one := []uint64{1, 0, 0, 0}
	sm2OrdMul(x, x, one)

	xOut := make([]byte, 32)
	sm2LittleToBig(xOut, x)
	return new(big.Int).SetBytes(xOut)
}

// fromBig converts a *big.Int into a format used by this code.
func fromBig(out []uint64, big *big.Int) {
	for i := range out {
		out[i] = 0
	}

	for i, v := range big.Bits() {
		out[i] = uint64(v)
	}
}

// sm2GetScalar endian-swaps the big-endian scalar value from in and writes it
// to out. If the scalar is equal or greater than the order of the group, it's
// reduced modulo that order.
func sm2GetScalar(out []uint64, in []byte) {
	n := new(big.Int).SetBytes(in)

	if n.Cmp(sm2P256.N) >= 0 {
		n.Mod(n, sm2P256.N)
	}
	fromBig(out, n)
}

// sm2Mul operates in a Montgomery domain with R = 2^256 mod p, where p is the
// underlying field of the curve. (See initSM2 for the value.) Thus rr here is
// R×R mod p. See comment in Inverse about how this is used.
var sm2rr = []uint64{0x0000000200000003, 0x00000002ffffffff, 0x0000000100000001, 0x0000000400000002}

func maybeReduceModP(in *big.Int) *big.Int {
	if in.Cmp(sm2P256.P) < 0 {
		return in
	}
	return new(big.Int).Mod(in, sm2P256.P)
}

func (curve sm2Curve) CombinedMult(bigX, bigY *big.Int, baseScalar, scalar []byte) (x, y *big.Int) {
	scalarReversed := make([]uint64, 4)
	var r1, r2 sm2Point
	sm2GetScalar(scalarReversed, baseScalar)
	r1IsInfinity := scalarIsZero(scalarReversed)
	r1.sm2BaseMult(scalarReversed)

	sm2GetScalar(scalarReversed, scalar)
	r2IsInfinity := scalarIsZero(scalarReversed)
	fromBig(r2.xyz[0:4], maybeReduceModP(bigX))
	fromBig(r2.xyz[4:8], maybeReduceModP(bigY))
	sm2Mul(r2.xyz[0:4], r2.xyz[0:4], sm2rr[:])
	sm2Mul(r2.xyz[4:8], r2.xyz[4:8], sm2rr[:])

	// This sets r2's Z value to 1, in the Montgomery domain.
	r2.xyz[8] = 0x0000000000000001
	r2.xyz[9] = 0x00000000ffffffff
	r2.xyz[10] = 0x0000000000000000
	r2.xyz[11] = 0x0000000100000000

	r2.sm2ScalarMult(scalarReversed)

	var sum, double sm2Point
	pointsEqual := sm2PointAddAsm(sum.xyz[:], r1.xyz[:], r2.xyz[:])
	sm2PointDoubleAsm(double.xyz[:], r1.xyz[:])
	sum.CopyConditional(&double, pointsEqual)
	sum.CopyConditional(&r1, r2IsInfinity)
	sum.CopyConditional(&r2, r1IsInfinity)

	return sum.sm2PointToAffine()
}

func (curve sm2Curve) ScalarBaseMult(scalar []byte) (x, y *big.Int) {
	scalarReversed := make([]uint64, 4)
	sm2GetScalar(scalarReversed, scalar)

	var r sm2Point
	r.sm2BaseMult(scalarReversed)
	return r.sm2PointToAffine()
}

func (curve sm2Curve) ScalarMult(bigX, bigY *big.Int, scalar []byte) (x, y *big.Int) {
	scalarReversed := make([]uint64, 4)
	sm2GetScalar(scalarReversed, scalar)

	var r sm2Point
	fromBig(r.xyz[0:4], maybeReduceModP(bigX))
	fromBig(r.xyz[4:8], maybeReduceModP(bigY))
	sm2Mul(r.xyz[0:4], r.xyz[0:4], sm2rr[:])
	sm2Mul(r.xyz[4:8], r.xyz[4:8], sm2rr[:])
	// This sets r2's Z value to 1, in the Montgomery domain.
	r.xyz[8] = 0x0000000000000001
	r.xyz[9] = 0x00000000ffffffff
	r.xyz[10] = 0x0000000000000000
	r.xyz[11] = 0x0000000100000000

	r.sm2ScalarMult(scalarReversed)
	return r.sm2PointToAffine()
}

// uint64IsZero returns 1 if x is zero and zero otherwise.
func uint64IsZero(x uint64) int {
	x = ^x
	x &= x >> 32
	x &= x >> 16
	x &= x >> 8
	x &= x >> 4
	x &= x >> 2
	x &= x >> 1
	return int(x & 1)
}

// scalarIsZero returns 1 if scalar represents the zero value, and zero
// otherwise.
func scalarIsZero(scalar []uint64) int {
	return uint64IsZero(scalar[0] | scalar[1] | scalar[2] | scalar[3])
}

func (p *sm2Point) sm2PointToAffine() (x, y *big.Int) {
	zInv := make([]uint64, 4)
	zInvSq := make([]uint64, 4)
	sm2Inverse(zInv, p.xyz[8:12])
	sm2Sqr(zInvSq, zInv, 1)
	sm2Mul(zInv, zInv, zInvSq)

	sm2Mul(zInvSq, p.xyz[0:4], zInvSq)
	sm2Mul(zInv, p.xyz[4:8], zInv)

	sm2FromMont(zInvSq, zInvSq)
	sm2FromMont(zInv, zInv)

	xOut := make([]byte, 32)
	yOut := make([]byte, 32)
	sm2LittleToBig(xOut, zInvSq)
	sm2LittleToBig(yOut, zInv)

	return new(big.Int).SetBytes(xOut), new(big.Int).SetBytes(yOut)
}

// CopyConditional copies overwrites p with src if v == 1, and leaves p
// unchanged if v == 0.
func (p *sm2Point) CopyConditional(src *sm2Point, v int) {
	pMask := uint64(v) - 1
	srcMask := ^pMask

	for i, n := range p.xyz {
		p.xyz[i] = (n & pMask) | (src.xyz[i] & srcMask)
	}
}

// sm2Inverse sets out to in^-1 mod p.
func sm2Inverse(out, in []uint64) {
	var stack [6 * 4]uint64
	p2 := stack[4*0 : 4*0+4]
	p4 := stack[4*1 : 4*1+4]
	p8 := stack[4*2 : 4*2+4]
	p16 := stack[4*3 : 4*3+4]
	p32 := stack[4*4 : 4*4+4]

	sm2Sqr(out, in, 1)
	sm2Mul(p2, out, in) // 3*p

	sm2Sqr(out, p2, 2)
	sm2Mul(p4, out, p2) // f*p

	sm2Sqr(out, p4, 4)
	sm2Mul(p8, out, p4) // ff*p

	sm2Sqr(out, p8, 8)
	sm2Mul(p16, out, p8) // ffff*p

	sm2Sqr(out, p16, 8)
	sm2Mul(out, out, p8) //ffffff*p

	sm2Sqr(out, out, 4)
	sm2Mul(out, out, p4) // fffffff*p

	sm2Sqr(out, out, 2)
	sm2Mul(out, out, p2) // fffffff*p

	sm2Sqr(out, out, 1)
	sm2Mul(out, out, in)

	sm2Sqr(out, out, 1) //fffffffe*p

	sm2Mul(p32, out, in) // ffffffff*p

	for j := 0; j < 4; j++ {
		sm2Sqr(out, out, 32)
		sm2Mul(out, out, p32)
	}

	sm2Sqr(out, out, 64)
	sm2Mul(out, out, p32)

	sm2Sqr(out, out, 16)
	sm2Mul(out, out, p16)

	sm2Sqr(out, out, 8)
	sm2Mul(out, out, p8)

	sm2Sqr(out, out, 4)
	sm2Mul(out, out, p4)

	sm2Sqr(out, out, 2)
	sm2Mul(out, out, p2)

	sm2Sqr(out, out, 2)
	sm2Mul(out, out, in)
}

func (p *sm2Point) sm2StorePoint(r *[16 * 4 * 3]uint64, index int) {
	copy(r[index*12:], p.xyz[:])
}

func boothW5(in uint) (int, int) {
	var s uint = ^((in >> 5) - 1)
	var d uint = (1 << 6) - in - 1
	d = (d & s) | (in & (^s))
	d = (d >> 1) + (d & 1)
	return int(d), int(s & 1)
}
func boothW6(in uint) (int, int) {
	var s uint = ^((in >> 6) - 1)
	var d uint = (1 << 7) - in - 1
	d = (d & s) | (in & (^s))
	d = (d >> 1) + (d & 1)
	return int(d), int(s & 1)
}

func boothW7(in uint) (int, int) {
	var s uint = ^((in >> 7) - 1)
	var d uint = (1 << 8) - in - 1
	d = (d & s) | (in & (^s))
	d = (d >> 1) + (d & 1)
	return int(d), int(s & 1)
}

func sm2InitTable() {
	sm2Precomputed = new([43][32 * 8]uint64)

	basePoint := []uint64{
		0x61328990f418029e, 0x3e7981eddca6c050, 0xd6a1ed99ac24c3c3, 0x91167a5ee1c13b05,
		0xc1354e593c2d0ddd, 0xc1f5e5788d3295fa, 0x8d4cfb066e2a48f8, 0x63cd65d481d735bd,
		0x0000000000000001, 0x00000000ffffffff, 0x0000000000000000, 0x0000000100000000,
	}
	t1 := make([]uint64, 12)
	t2 := make([]uint64, 12)
	copy(t2, basePoint)

	zInv := make([]uint64, 4)
	zInvSq := make([]uint64, 4)
	for j := 0; j < 32; j++ {
		copy(t1, t2)
		for i := 0; i < 43; i++ {
			// The window size is 6 so we need to double 6 times.
			if i != 0 {
				for k := 0; k < 6; k++ {
					sm2PointDoubleAsm(t1, t1)
				}
			}
			// Convert the point to affine form. (Its values are
			// still in Montgomery form however.)
			sm2Inverse(zInv, t1[8:12])
			sm2Sqr(zInvSq, zInv, 1)
			sm2Mul(zInv, zInv, zInvSq)

			sm2Mul(t1[:4], t1[:4], zInvSq)
			sm2Mul(t1[4:8], t1[4:8], zInv)

			copy(t1[8:12], basePoint[8:12])
			// Update the table entry
			copy(sm2Precomputed[i][j*8:], t1[:8])
		}
		if j == 0 {
			sm2PointDoubleAsm(t2, basePoint)
		} else {
			sm2PointAddAsm(t2, t2, basePoint)
		}
	}
}

func (p *sm2Point) sm2BaseMult(scalar []uint64) {
	sm2precomputeOnce.Do(sm2InitTable)

	wvalue := (scalar[0] << 1) & 0x7f
	sel, sign := boothW6(uint(wvalue))
	sm2SelectBase(p.xyz[0:8], sm2Precomputed[0][0:], sel)
	sm2NegCond(p.xyz[4:8], sign)

	// (This is one, in the Montgomery domain.)
	p.xyz[8] = 0x0000000000000001
	p.xyz[9] = 0x00000000ffffffff
	p.xyz[10] = 0x0000000000000000
	p.xyz[11] = 0x0000000100000000

	var t0 sm2Point
	// (This is one, in the Montgomery domain.)
	t0.xyz[8] = 0x0000000000000001
	t0.xyz[9] = 0x00000000ffffffff
	t0.xyz[10] = 0x0000000000000000
	t0.xyz[11] = 0x0000000100000000

	index := uint(5)
	zero := sel

	for i := 1; i < 43; i++ {
		if index < 192 {
			wvalue = ((scalar[index/64] >> (index % 64)) + (scalar[index/64+1] << (64 - (index % 64)))) & 0x7f
		} else {
			wvalue = (scalar[index/64] >> (index % 64)) & 0x7f
		}
		index += 6
		sel, sign = boothW6(uint(wvalue))
		sm2SelectBase(t0.xyz[0:8], sm2Precomputed[i][0:], sel)
		sm2PointAddAffineAsm(p.xyz[0:12], p.xyz[0:12], t0.xyz[0:8], sign, sel, zero)
		zero |= sel
	}
}

func (p *sm2Point) sm2ScalarMult(scalar []uint64) {
	// precomp is a table of precomputed points that stores powers of p
	// from p^1 to p^16.
	var precomp [16 * 4 * 3]uint64
	var t0, t1, t2, t3 sm2Point

	// Prepare the table
	p.sm2StorePoint(&precomp, 0) // 1

	sm2PointDoubleAsm(t0.xyz[:], p.xyz[:])
	sm2PointDoubleAsm(t1.xyz[:], t0.xyz[:])
	sm2PointDoubleAsm(t2.xyz[:], t1.xyz[:])
	sm2PointDoubleAsm(t3.xyz[:], t2.xyz[:])
	t0.sm2StorePoint(&precomp, 1)  // 2
	t1.sm2StorePoint(&precomp, 3)  // 4
	t2.sm2StorePoint(&precomp, 7)  // 8
	t3.sm2StorePoint(&precomp, 15) // 16

	sm2PointAddAsm(t0.xyz[:], t0.xyz[:], p.xyz[:])
	sm2PointAddAsm(t1.xyz[:], t1.xyz[:], p.xyz[:])
	sm2PointAddAsm(t2.xyz[:], t2.xyz[:], p.xyz[:])
	t0.sm2StorePoint(&precomp, 2) // 3
	t1.sm2StorePoint(&precomp, 4) // 5
	t2.sm2StorePoint(&precomp, 8) // 9

	sm2PointDoubleAsm(t0.xyz[:], t0.xyz[:])
	sm2PointDoubleAsm(t1.xyz[:], t1.xyz[:])
	t0.sm2StorePoint(&precomp, 5) // 6
	t1.sm2StorePoint(&precomp, 9) // 10

	sm2PointAddAsm(t2.xyz[:], t0.xyz[:], p.xyz[:])
	sm2PointAddAsm(t1.xyz[:], t1.xyz[:], p.xyz[:])
	t2.sm2StorePoint(&precomp, 6)  // 7
	t1.sm2StorePoint(&precomp, 10) // 11

	sm2PointDoubleAsm(t0.xyz[:], t0.xyz[:])
	sm2PointDoubleAsm(t2.xyz[:], t2.xyz[:])
	t0.sm2StorePoint(&precomp, 11) // 12
	t2.sm2StorePoint(&precomp, 13) // 14

	sm2PointAddAsm(t0.xyz[:], t0.xyz[:], p.xyz[:])
	sm2PointAddAsm(t2.xyz[:], t2.xyz[:], p.xyz[:])
	t0.sm2StorePoint(&precomp, 12) // 13
	t2.sm2StorePoint(&precomp, 14) // 15

	// Start scanning the window from top bit
	index := uint(254)
	var sel, sign int

	wvalue := (scalar[index/64] >> (index % 64)) & 0x3f
	sel, _ = boothW5(uint(wvalue))

	sm2Select(p.xyz[0:12], precomp[0:], sel)
	zero := sel

	for index > 4 {
		index -= 5
		sm2PointDoubleAsm(p.xyz[:], p.xyz[:])
		sm2PointDoubleAsm(p.xyz[:], p.xyz[:])
		sm2PointDoubleAsm(p.xyz[:], p.xyz[:])
		sm2PointDoubleAsm(p.xyz[:], p.xyz[:])
		sm2PointDoubleAsm(p.xyz[:], p.xyz[:])

		if index < 192 {
			wvalue = ((scalar[index/64] >> (index % 64)) + (scalar[index/64+1] << (64 - (index % 64)))) & 0x3f
		} else {
			wvalue = (scalar[index/64] >> (index % 64)) & 0x3f
		}

		sel, sign = boothW5(uint(wvalue))

		sm2Select(t0.xyz[0:], precomp[0:], sel)
		sm2NegCond(t0.xyz[4:8], sign)
		sm2PointAddAsm(t1.xyz[:], p.xyz[:], t0.xyz[:])
		sm2MovCond(t1.xyz[0:12], t1.xyz[0:12], p.xyz[0:12], sel)
		sm2MovCond(p.xyz[0:12], t1.xyz[0:12], t0.xyz[0:12], zero)
		zero |= sel
	}

	sm2PointDoubleAsm(p.xyz[:], p.xyz[:])
	sm2PointDoubleAsm(p.xyz[:], p.xyz[:])
	sm2PointDoubleAsm(p.xyz[:], p.xyz[:])
	sm2PointDoubleAsm(p.xyz[:], p.xyz[:])
	sm2PointDoubleAsm(p.xyz[:], p.xyz[:])

	wvalue = (scalar[0] << 1) & 0x3f
	sel, sign = boothW5(uint(wvalue))

	sm2Select(t0.xyz[0:], precomp[0:], sel)
	sm2NegCond(t0.xyz[4:8], sign)
	sm2PointAddAsm(t1.xyz[:], p.xyz[:], t0.xyz[:])
	sm2MovCond(t1.xyz[0:12], t1.xyz[0:12], p.xyz[0:12], sel)
	sm2MovCond(p.xyz[0:12], t1.xyz[0:12], t0.xyz[0:12], zero)
}

var initOnce sync.Once

// SM2 returns a Curve which implements SM2
// The cryptographic operations are implemented using constant-time algorithms.
func SM2() elliptic.Curve {
	initOnce.Do(initSM2)
	return sm2P256
}
