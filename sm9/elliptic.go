package sm9

import (
	"crypto/elliptic"
	"errors"
	"fmt"
	"math/big"
	"sync"
)

type Curve struct {
	*elliptic.CurveParams
}

var curve Curve
var initOnce sync.Once

func initSM9() {
	curve.CurveParams = &elliptic.CurveParams{}
	curve.CurveParams.P, _ = new(big.Int).SetString("B640000002A3A6F1D603AB4FF58EC74521F2934B1A7AEEDBE56F9B27E351457D", 16)
	curve.CurveParams.N, _ = new(big.Int).SetString("B640000002A3A6F1D603AB4FF58EC74449F2934B18EA8BEEE56EE19CD69ECF25", 16)
	curve.CurveParams.B = new(big.Int).SetInt64(5)
	curve.CurveParams.Gx, _ = new(big.Int).SetString("93DE051D62BF718FF5ED0704487D01D6E1E4086909DC3280E8C4E4817C66DDDD", 16)
	curve.CurveParams.Gy, _ = new(big.Int).SetString("21FE8DDA4F21E607631065125C395BBC1C1C00CBFA6024350C464CD70A3EA616", 16)
	curve.CurveParams.BitSize = 256
	curve.CurveParams.Name = "SM9"
}

// SM9 return the elliptic.Curve interface of SM9 curve
func SM9() *Curve {
	initOnce.Do(initSM9)
	return &curve
}

// Params returns the parameters for the curve.
func (curve *Curve) Params() *elliptic.CurveParams {
	return curve.CurveParams
}

// IsOnCurve reports whether the given (x,y) lies on the curve.
func (curve *Curve) IsOnCurve(x, y *big.Int) bool {

	// y² = x³ + b
	y2 := new(big.Int).Mul(y, y)
	y2.Mod(y2, curve.P)

	x3 := new(big.Int).Mul(x, x)
	x3.Mul(x3, x)

	x3.Add(x3, curve.B)
	x3.Mod(x3, curve.P)

	return x3.Cmp(y2) == 0
}

// Add returns the sum of (x1,y1) and (x2,y2)
func (curve *Curve) Add(x1, y1, x2, y2 *big.Int) (x, y *big.Int) {
	p1 := BigToG1(x1, y1)
	p2 := BigToG1(x2, y2)
	rp := new(G1).Add(p1, p2)
	return G1ToBig(rp)
}

// Double returns 2*(x,y)
func (curve *Curve) Double(x1, y1 *big.Int) (x, y *big.Int) {
	return curve.Add(x1, y1, x1, y1)
}

// ScalarMult returns k*(Bx,By) where k is a number in big-endian form.
func (curve *Curve) ScalarMult(x1, y1 *big.Int, k []byte) (x, y *big.Int) {
	p1 := BigToG1(x1, y1)
	bigK := new(big.Int).SetBytes(k)
	rp := new(G1).ScalarMult(p1, bigK)
	return G1ToBig(rp)
}

// ScalarBaseMult returns k*G, where G is the base point of the group
// and k is an integer in big-endian form.
func (curve *Curve) ScalarBaseMult(k []byte) (x, y *big.Int) {
	bigK := new(big.Int).SetBytes(k)
	rp := new(G1).ScalarBaseMult(bigK)
	return G1ToBig(rp)
}

//Neg is (x, -y)
func (curve *Curve) Neg(x1, y1 *big.Int) (x, y *big.Int) {
	return new(big.Int).Set(x1), new(big.Int).Sub(curve.Params().P, y1)
}

//CombinedMult do baseScalar*G + scalar*(X,Y)
func (curve *Curve) CombinedMult(bigX, bigY *big.Int, baseScalar, scalar []byte) (x, y *big.Int) {
	x1, y1 := curve.ScalarBaseMult(baseScalar)
	x2, y2 := curve.ScalarMult(bigX, bigY, scalar)
	return curve.Add(x1, y1, x2, y2)
}

func BigToG1(x, y *big.Int) *G1 {
	m := make([]byte, 64)
	xBytes := x.Bytes()
	yBytes := y.Bytes()
	copy(m[32-len(xBytes):32], xBytes)
	copy(m[64-len(yBytes):64], yBytes)
	r := new(G1)
	_, success := r.Unmarshal(m)
	if success != nil {
		fmt.Printf("error in big int to G1")
		return nil
	}
	return r
}

func G1ToBig(g1 *G1) (x, y *big.Int) {
	m := g1.Marshal()
	return new(big.Int).SetBytes(m[0:32]), new(big.Int).SetBytes(m[32:64])
}

func (curve *Curve) Compress(x, y *big.Int) []byte {
	return CompressP(x, y)
}

func (curve *Curve) Decompress(in []byte) (x, y *big.Int, err error) {
	x, y = DecompressP(in)
	if x == nil || y == nil {
		return nil, nil, errors.New("decompress fail")
	}
	return x, y, nil
}

func BytesToG2(in []byte) *G2 {
	point := new(G2)
	if _, success := point.Unmarshal(in); success != nil {
		fmt.Printf("point is not on curve G2")
		return nil
	}
	return point
}

func G2ToBytes(point *G2) []byte {
	return point.Marshal()
}

func BytesToGt(in []byte) *GT {
	point := new(GT)
	if _, success := point.Unmarshal(in); success != nil {
		fmt.Printf("point is not on curve GT")
		return nil
	}
	return point
}

func GtToBytes(point *GT) []byte {
	return point.Marshal()
}
