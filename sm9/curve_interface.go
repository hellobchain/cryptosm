package sm9

import (
	"math/big"
)

// type impl struct{}

/*
type ECCInternal interface {
	AffineToPoint(x, y *big.Int) (point *ECCInternalPoint)
	PointToAffine(point *ECCInternalPoint) (x, y *big.Int)
	//JacobianAdd(p1, p2 *ECCInternalPoint) (r *ECCInternalPoint, equal int)
	JacobianAdd(p1, p2 *ECCInternalPoint) (r *ECCInternalPoint)
	JacobianDouble(p1 *ECCInternalPoint) (r *ECCInternalPoint)
	PointNegCondition(p1 *ECCInternalPoint, c int)
	FieldMul(res, in1, in2 []uint64)
	MontgomaryR() []uint64
	MontgomaryR2() []uint64
	ModSqrtP(a *big.Int) *big.Int
	ModInverseP(a []uint64) []uint64
	ModInverseOrder(a *big.Int) *big.Int
}

type ECCInternalPoint struct {
	XYZ [12]uint64
}
*/
// func New() supercurve.JacobianIF { //name is not good
// 	return &impl{}
// }

func InternalP2CurveP(xyz *[12]uint64) (out *curvePoint) {
	out = &curvePoint{}
	copy(out.x[:], xyz[0:4])
	copy(out.y[:], xyz[4:8])
	copy(out.z[:], xyz[8:12])
	return
}

func InternalP2CurveAffineP(xyz *[12]uint64) (out *curvePointAffine) {
	out = &curvePointAffine{}
	copy(out.x[:], xyz[0:4])
	copy(out.y[:], xyz[4:8])
	return
}

func CurveP2InternalP(point *curvePoint) (out *[12]uint64) {
	out = &[12]uint64{}
	copy(out[0:4], point.x[:])
	copy(out[4:8], point.y[:])
	copy(out[8:12], point.z[:])
	return
}

func AffineToPoint(x, y *big.Int) (xyz *[12]uint64) {

	xyz = &[12]uint64{}

	t := new(big.Int).Mul(x, mogo_bigInt)
	t.Mod(t, p)
	fromBig(xyz[0:4], t) // change x to montgomery format

	t.Mul(y, mogo_bigInt)
	t.Mod(t, p)
	fromBig(xyz[4:8], t) // change y to montgomery format

	xyz[8] = 0x1a9064d81caeba83
	xyz[9] = 0xde0d6cb4e5851124
	xyz[10] = 0x29fc54b00a7138ba
	xyz[11] = 0x49bffffffd5c590e

	return
}

func PointToAffine(xyz *[12]uint64) (x, y *big.Int) {

	ret1 := make([]byte, 32)
	ret2 := make([]byte, 32)

	temp := InternalP2CurveP(xyz)
	temp.MakeAffine()

	montDecode(&temp.x, &temp.x)
	montDecode(&temp.y, &temp.y)

	temp.x.Marshal(ret1)
	temp.y.Marshal(ret2)

	x1 := new(big.Int).SetBytes(ret1)
	y1 := new(big.Int).SetBytes(ret2)

	return x1, y1
}

func JacobianAdd(p1, p2 *[12]uint64) (r1 *[12]uint64) {

	t1 := InternalP2CurveP(p1)
	t2 := InternalP2CurveP(p2)
	t1.Add(t1, t2)

	r1 = CurveP2InternalP(t1)
	return
}

func JacobianAddAffine(p1, p2 *[12]uint64) (r1 *[12]uint64) {

	t1 := InternalP2CurveP(p1)
	t2 := InternalP2CurveAffineP(p2)
	t1.Add_JA(t1, t2)

	r1 = CurveP2InternalP(t1)
	return
}

func JacobianDouble(p1 *[12]uint64) (r1 *[12]uint64) {
	t1 := InternalP2CurveP(p1)
	t1.Double(t1)
	r1 = CurveP2InternalP(t1)
	return
}

func PointNegCondition(p1 *[12]uint64, condition int) { //if c!=0, p1 <--(-p1)
	if condition == 0 {
		return
	}
	t1 := InternalP2CurveP(p1)
	t1.Neg(t1)
	copy(p1[0:4], t1.x[:])
	copy(p1[4:8], t1.y[:])
	copy(p1[8:12], t1.z[:])
	return
}

func isOdd(a *big.Int) bool {
	return a.Bit(0) == 1
}

// decompressPoint decompresses a point on the given curve given the X point and the solution to use.
func DecompressPointX4Y(x *big.Int, ybit bool) *big.Int {
	// x³ - 3x + b
	x3 := new(big.Int).Mul(x, x)
	x3.Mul(x3, x)

	//	threeX := new(big.Int).Lsh(x, 1)
	//	threeX.Add(threeX, x)

	//	x3.Sub(x3, threeX)
	x3.Add(x3, big.NewInt(5))
	x3.Mod(x3, Sm9_p)

	//y := new(big.Int).ModSqrt(x3, curve.Params().P)
	//y := new(big.Int)
	//y:=new(big.Int).ModSqrt(x3, SM9().Params().P)
	y := ModSqrtP(x3)

	if y == nil {
		return nil
	}

	if ybit != isOdd(y) {
		y.Sub(Sm9_p, y)
	}

	return y
}

// Compress compress a point to a byte string
func CompressP(x, y *big.Int) []byte {
	if x.Sign() == 0 && y.Sign() == 0 {
		return []byte{0x00}
	}
	//byteLen := (sm9inter.SM9().Params().BitSize+7)/8
	byteLen := 32

	ret := make([]byte, 1+byteLen)
	if y.Bit(0) == 0 {
		ret[0] = 0x02
	} else {
		ret[0] = 0x03
	}

	xBytes := x.Bytes()
	copy(ret[byteLen-len(xBytes)+1:], xBytes)
	return ret
}

// Decompress decompress a byte string to a point
func DecompressP(data []byte) (x, y *big.Int) {
	if data == nil {
		return nil, nil
	}
	if len(data) == 0 {
		return nil, nil
	}
	//byteLen := (sm9inter.SM9().Params().BitSize+7)/8
	byteLen := 32
	switch data[0] {
	case 0x00:
		if len(data) == 1 {
			return new(big.Int), new(big.Int)
		}
		return nil, nil
	case 0x02, 0x03:
		{
			if len(data) != 1+byteLen {
				return nil, nil
			}
			x = new(big.Int).SetBytes(data[1:])

			// x³ + b
			x3 := new(big.Int).Mul(x, x)
			x3.Mul(x3, x)

			x3.Add(x3, big.NewInt(5))
			x3.Mod(x3, Sm9_p)

			y := ModSqrtP(x3)

			if y == nil {
				return nil, nil
			}
			if y.Bit(0) != uint(data[0]&0x01) {
				y.Sub(Sm9_p, y)
			}
			return x, y
		}
	case 0x04:
		{
			x = new(big.Int).SetBytes(data[1 : 1+byteLen])
			y = new(big.Int).SetBytes(data[1+byteLen:])
			return x, y
		}
	default:
		return nil, nil
	}
}

func FieldMul(res, in1, in2 []uint64) { //res,in1,in2:mogo format
	var temp1 = &gfP{in1[0], in1[1], in1[2], in1[3]}
	var temp2 = &gfP{in2[0], in2[1], in2[2], in2[3]}
	gfpMul(temp2, temp2, temp1)
	for i := 0; i < 4; i++ {
		res[i] = temp2[i]
	}
}

func MontgomaryR() []uint64 {
	return []uint64{0x1a9064d81caeba83, 0xde0d6cb4e5851124, 0x29fc54b00a7138ba, 0x49bffffffd5c590e}
}

func MontgomaryR2() []uint64 {
	return []uint64{0x27dea312b417e2d2, 0x88f8105fae1a5d3f, 0xe479b522d6706e7b, 0x2ea795a656f62fbd}
	//2ea795a656f62fbd e479b522d6706e7b 88f8105fae1a5d3f 27dea312b417e2d2
}

/*
func (c *impl) ModSqrtP(a *big.Int) *big.Int {
	var temp = big.NewInt(1)
	temp.ModSqrt(a, p)
	return temp
}
*/

func ModSqrtP(a *big.Int) *big.Int { //a>=0,a!=2;//ALG see:SM_1_P29

	var temp = big.NewInt(1)
	var one = big.NewInt(1)
	var fone = big.NewInt(-1)
	fone.Mod(fone, p)
	temp.Exp(a, const2uplus1, p) //g^(2u+1)

	if temp.Cmp(one) == 0 {
		temp.Exp(a, constuplus1, p)
	} else if temp.Cmp(fone) == 0 {
		one.Add(a, a) //2g
		one.Mod(one, p)
		fone.Add(one, one) //4g
		fone.Mod(fone, p)
		fone.Exp(fone, constu, p)
		fone.Mul(fone, one)
		fone.Mod(fone, p)
		temp.Set(fone)
	} else {
		return a
	}
	return temp
}

//func montEncode(c, a *gfP) { gfpMul(c, a, r2) }
func ModInverseP(a []uint64) []uint64 { //a: mogo format;return mogo format too;
	temp := &gfP{a[0], a[1], a[2], a[3]}
	temp.Invert(temp)

	X := make([]uint64, 4)
	for i := 0; i < 4; i++ {
		X[i] = temp[i]
	}
	return X
}

func ModInverseOrder(a *big.Int) *big.Int {
	var b = big.NewInt(1)
	b.ModInverse(a, Order)
	return b
}

func InitBigTable(Bx, By *big.Int) *[43][32 * 8]uint64 {
	var precomputed [43][32 * 8]uint64

	btpoint := AffineToPoint(Bx, By)

	basepoint := &curvePoint{}
	for i := 0; i < 4; i++ {
		basepoint.x[i] = btpoint[i]
		basepoint.y[i] = btpoint[i+4]
		basepoint.z[i] = btpoint[i+8]
	}

	t1 := new(G1).Set(&G1{basepoint})
	t2 := new(G1).Set(&G1{basepoint})
	t3 := new(G1).Set(&G1{basepoint})
	var count int
	for j := 0; j < 32; j++ {
		t1.Set(t2)

		for i := 0; i < 43; i++ {
			// The window size is 6 so we need to double 6 times.
			if i != 0 {
				for k := 0; k < 6; k++ {
					//p256PointDoubleAsm(t1, t1)
					t3.p.Double(t1.p)
					t1.Set(t3)
				}
			}

			t1.p.MakeAffine()
			t1.p.z[0] = 0x1a9064d81caeba83
			t1.p.z[1] = 0xde0d6cb4e5851124
			t1.p.z[2] = 0x29fc54b00a7138ba
			t1.p.z[3] = 0x49bffffffd5c590e

			if t1.p.IsInfinity() {
				for count = 0; count < 4; count++ {
					precomputed[i][j*8+count] = 0
					precomputed[i][j*8+count+4] = 0
				}
			}
			for count = 0; count < 4; count++ {
				precomputed[i][j*8+count] = t1.p.x[count]
				precomputed[i][j*8+count+4] = t1.p.y[count]
			}

		}
		if j == 0 {
			//p256PointDoubleAsm(t2, basePoint)
			t2.p.Double(basepoint)
		} else { //t2:3-->32G
			t2.p.Add(t2.p, basepoint)
		}
	}
	return &precomputed
}
