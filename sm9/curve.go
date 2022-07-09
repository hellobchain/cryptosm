package sm9

import (
	"math/big"
	"sync"
)

/*****************************************************************************/
// curvePoint implements the elliptic curve y²=x³+3. Points are kept in Jacobian
// form and t=z² when valid. G₁ is the set of points of this curve on GF(p).
// sm9:y²=x³+5
type curvePoint struct {
	x, y, z, t gfP
}
type curvePointAffine struct {
	x, y gfP
}

var curveB = newGFp(5) //sm9:b

// curveGen is the generator of G₁.
var curveGen = &curvePoint{
	x: *G1x,
	y: *G1y,
	z: *newGFp(1),
	t: *newGFp(1),
}

func (c *curvePoint) String() string {
	c.MakeAffine()
	x, y := &gfP{}, &gfP{}
	montDecode(x, &c.x)
	montDecode(y, &c.y)
	return "(" + x.String() + ", " + y.String() + ")"
}

func (c *curvePoint) Set(a *curvePoint) {
	c.x.Set(&a.x)
	c.y.Set(&a.y)
	c.z.Set(&a.z)
	c.t.Set(&a.t)
}

// IsOnCurve returns true iff c is on the curve.
func (c *curvePoint) IsOnCurve() bool {
	c.MakeAffine()
	if c.IsInfinity() {
		return true
	}

	y2, x3 := &gfP{}, &gfP{}
	gfpMul(y2, &c.y, &c.y)
	gfpMul(x3, &c.x, &c.x)
	gfpMul(x3, x3, &c.x)
	gfpAdd(x3, x3, curveB)

	return *y2 == *x3
}

func (c *curvePoint) SetInfinity() {
	c.x = gfP{0}
	c.y = *newGFp(1)
	c.z = gfP{0}
	c.t = gfP{0}
}

func (c *curvePoint) IsInfinity() bool {
	return c.z == gfP{0}
}

//simon:add Add_JA for jacobi+affine 190103
func (c *curvePoint) Add_JA(a *curvePoint, b *curvePointAffine) {
	//hyperelliptic.org/EFD/g1p/auto-shortw-jacobian-0.html#addition-madd-2007-bl
	if a.IsInfinity() {
		//c.Set(b)
		c.x.Set(&b.x)
		c.y.Set(&b.y)
		c.z = *newGFp(1)
		c.t = *newGFp(1)
		return
	}

	z12, u2, s2, temp, H, HH, I, J, r, V := &gfP{}, &gfP{}, &gfP{}, &gfP{}, &gfP{}, &gfP{}, &gfP{}, &gfP{}, &gfP{}, &gfP{}
	gfpMul(z12, &a.z, &a.z) //z1z1 = z1^2
	gfpMul(temp, &a.z, z12) //z1*z1z1
	gfpMul(u2, &b.x, z12)
	gfpMul(s2, &b.y, temp) //s2 = y2*z1*z1z1
	gfpSub(H, u2, &a.x)
	gfpMul(HH, H, H)
	gfpAdd(I, HH, HH)
	gfpAdd(I, I, I) //I = 4*H
	gfpMul(J, H, I)
	gfpSub(r, s2, &a.y)
	gfpAdd(r, r, r)    //r = 2*(S2-Y1)
	gfpMul(V, &a.x, I) //V = X1*I
	gfpMul(temp, r, r)
	gfpSub(temp, temp, J)
	gfpSub(temp, temp, V)
	gfpSub(&c.x, temp, V)

	gfpSub(temp, V, &c.x)
	gfpMul(temp, r, temp) //temp = r*(V-X3)
	gfpMul(u2, &a.y, J)
	gfpSub(temp, temp, u2)
	gfpSub(&c.y, temp, u2) //Y3
	gfpAdd(temp, &a.z, H)
	gfpMul(temp, temp, temp)
	gfpSub(&c.z, temp, z12)
	gfpSub(&c.z, &c.z, HH)

}

func (c *curvePoint) Add(a, b *curvePoint) {
	if a.IsInfinity() {
		c.Set(b)
		return
	}
	if b.IsInfinity() {
		c.Set(a)
		return
	}

	// See http://hyperelliptic.org/EFD/g1p/auto-code/shortw/jacobian-0/addition/add-2007-bl.op3

	// Normalize the points by replacing a = [x1:y1:z1] and b = [x2:y2:z2]
	// by [u1:s1:z1·z2] and [u2:s2:z1·z2]
	// where u1 = x1·z2², s1 = y1·z2³ and u1 = x2·z1², s2 = y2·z1³
	z12, z22 := &gfP{}, &gfP{}
	gfpMul(z12, &a.z, &a.z)
	gfpMul(z22, &b.z, &b.z)

	u1, u2 := &gfP{}, &gfP{}
	gfpMul(u1, &a.x, z22)
	gfpMul(u2, &b.x, z12)

	t, s1 := &gfP{}, &gfP{}
	gfpMul(t, &b.z, z22)
	gfpMul(s1, &a.y, t)

	s2 := &gfP{}
	gfpMul(t, &a.z, z12)
	gfpMul(s2, &b.y, t)

	// Compute x = (2h)²(s²-u1-u2)
	// where s = (s2-s1)/(u2-u1) is the slope of the line through
	// (u1,s1) and (u2,s2). The extra factor 2h = 2(u2-u1) comes from the value of z below.
	// This is also:
	// 4(s2-s1)² - 4h²(u1+u2) = 4(s2-s1)² - 4h³ - 4h²(2u1)
	//                        = r² - j - 2v
	// with the notations below.
	h := &gfP{}
	gfpSub(h, u2, u1)
	xEqual := *h == gfP{0}

	gfpAdd(t, h, h)
	// i = 4h²
	i := &gfP{}
	gfpMul(i, t, t)
	// j = 4h³
	j := &gfP{}
	gfpMul(j, h, i)

	gfpSub(t, s2, s1)
	yEqual := *t == gfP{0}
	if xEqual && yEqual {
		c.Double(a)
		return
	}
	r := &gfP{}
	gfpAdd(r, t, t)

	v := &gfP{}
	gfpMul(v, u1, i)

	// t4 = 4(s2-s1)²
	t4, t6 := &gfP{}, &gfP{}
	gfpMul(t4, r, r)
	gfpAdd(t, v, v)
	gfpSub(t6, t4, j)

	gfpSub(&c.x, t6, t)

	// Set y = -(2h)³(s1 + s*(x/4h²-u1))
	// This is also
	// y = - 2·s1·j - (s2-s1)(2x - 2i·u1) = r(v-x) - 2·s1·j
	gfpSub(t, v, &c.x) // t7
	gfpMul(t4, s1, j)  // t8
	gfpAdd(t6, t4, t4) // t9
	gfpMul(t4, r, t)   // t10
	gfpSub(&c.y, t4, t6)

	// Set z = 2(u2-u1)·z1·z2 = 2h·z1·z2
	gfpAdd(t, &a.z, &b.z) // t11
	gfpMul(t4, t, t)      // t12
	gfpSub(t, t4, z12)    // t13
	gfpSub(t4, t, z22)    // t14
	gfpMul(&c.z, t4, h)
}

func (c *curvePoint) Double(a *curvePoint) { //simon: modify the case:a!=c in 20190107
	// See http://hyperelliptic.org/EFD/g1p/auto-code/shortw/jacobian-0/doubling/dbl-2009-l.op3
	if a.IsInfinity() {
		c.Set(a)
		return
	}
	temp := &curvePoint{}
	temp.Set(c)

	A, B, C := &gfP{}, &gfP{}, &gfP{}
	gfpMul(A, &a.x, &a.x)
	gfpMul(B, &a.y, &a.y)
	gfpMul(C, B, B)

	t, t2 := &gfP{}, &gfP{}
	gfpAdd(t, &a.x, B)
	gfpMul(t2, t, t)
	gfpSub(t, t2, A)
	gfpSub(t2, t, C)

	d, e, f := &gfP{}, &gfP{}, &gfP{}
	gfpAdd(d, t2, t2)
	gfpAdd(t, A, A)
	gfpAdd(e, t, A)
	gfpMul(f, e, e)

	gfpAdd(t, d, d)
	gfpSub(&temp.x, f, t)

	gfpAdd(t, C, C)
	gfpAdd(t2, t, t)
	gfpAdd(t, t2, t2)
	gfpSub(&temp.y, d, &temp.x)
	gfpMul(t2, e, &temp.y)
	gfpSub(&temp.y, t2, t)

	gfpMul(t, &a.y, &a.z)
	gfpAdd(&temp.z, t, t)
	c.Set(temp)
}

/******************************************/

func (c *curvePoint) MakeAffine() {
	if c.z == *newGFp(1) { //mogo format 1
		return
	} else if c.z == *newGFp(0) {
		c.x = gfP{0}
		c.y = *newGFp(1)
		c.t = gfP{0}
		return
	}

	zInv := &gfP{}
	zInv.Invert(&c.z)

	t, zInv2 := &gfP{}, &gfP{}
	gfpMul(t, &c.y, zInv)
	gfpMul(zInv2, zInv, zInv)

	gfpMul(&c.x, &c.x, zInv2)
	gfpMul(&c.y, t, zInv2)

	c.z = *newGFp(1)
	c.t = *newGFp(1)
}

func (c *curvePoint) Neg(a *curvePoint) {
	c.x.Set(&a.x)
	gfpNeg(&c.y, &a.y)
	c.z.Set(&a.z)
	c.t = gfP{0}
}

// fromBig converts a *big.Int into a format used by this code.
func fromBig(out []uint64, big *big.Int) { //little endian
	for i := range out {
		out[i] = 0
	}

	for i, v := range big.Bits() {
		out[i] = uint64(v)
	}
}

/*********************************************/
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

func (c *curvePoint) sm9CombinedMult(H *curvePoint, baseScalar, scalar *big.Int) {
	scalarReversed := make([]uint64, 4)
	var r1, r2 curvePoint
	//p256GetScalar(scalarReversed, baseScalar)
	fromBig(scalarReversed, baseScalar)
	r1IsInfinity := scalarIsZero(scalarReversed)
	if r1IsInfinity == 1 {
		r1.SetInfinity()
	} else {
		r1.sm9BaseMult(scalarReversed)
	}

	//p256GetScalar(scalarReversed, scalar)
	fromBig(scalarReversed, scalar)
	r2IsInfinity := scalarIsZero(scalarReversed)
	if r2IsInfinity == 1 {
		c.Set(&r1)
		return
	}

	r2.Mul(H, scalar)

	c.Add(&r2, &r1)
}

//simon:rewrite k*baseG on 1214
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

var p256Precomputed1 *[43][32 * 8]uint64
var precomputeOnce sync.Once
var ginfinite = &curvePoint{
	x: *newGFp(0),
	y: *newGFp(1),
	z: *newGFp(0),
	t: *newGFp(0),
}

/*****************************************************/
func initTable() {
	//fmt.Println("enter JAmode initable test:")
	p256Precomputed1 = new([43][32 * 8]uint64)

	t1 := new(G1).Set(&G1{curveGen})
	t2 := new(G1).Set(&G1{curveGen})

	var count int
	for j := 0; j < 32; j++ {
		t1.Set(t2)
		for i := 0; i < 43; i++ {
			// The window size is 6 so we need to double 6 times.
			if i != 0 {
				for k := 0; k < 6; k++ {
					t1.p.Double(t1.p)
				}
			}
			// Convert the point to affine form. (Its values are
			// still in Montgomery form however.)
			zInv, zInvSq := &gfP{}, &gfP{}
			zInv.Invert(&t1.p.z)
			gfpMul(zInvSq, zInv, zInv) //z^2
			gfpMul(zInv, zInv, zInvSq) //z^3
			gfpMul(&t1.p.x, &t1.p.x, zInvSq)
			gfpMul(&t1.p.y, &t1.p.y, zInv)
			t1.p.z = *newGFp(1)
			for count = 0; count < 4; count++ {
				p256Precomputed1[i][j*8+count] = t1.p.x[count]
				p256Precomputed1[i][j*8+count+4] = t1.p.y[count]
			}
		}
		if j == 0 {
			t2.p.Double(curveGen) //2G
		} else { //t2:3G-->32G
			for count = 0; count < 4; count++ {
				t2.p.x[count] = g32[j+1][count]
				t2.p.y[count] = g32[j+1][count+4]
				t2.p.z[count] = g32[j+1][count+8]
			}
		}
	}
}

func (c *curvePoint) sm9BaseMult(scalar []uint64) {
	precomputeOnce.Do(initTable)
	var count int

	t1 := new(G1).Set(&G1{curveGen})
	t2 := new(G1).Set(&G1{curveGen})
	//t3 := new(G1).Set(&G1{curveGen})
	t4 := &curvePointAffine{}

	p := make([]uint64, 12)
	wvalue := (scalar[0] << 1) & 0x7f
	sel, sign := boothW6(uint(wvalue))

	if sel > 0 {
		copy(p, p256Precomputed1[0][(sel-1)*8:(sel-1)*8+8])
		p[8] = 0x1a9064d81caeba83
		p[9] = 0xde0d6cb4e5851124
		p[10] = 0x29fc54b00a7138ba
		p[11] = 0x49bffffffd5c590e
	}

	//p256NegCond(p[4:8], sign)

	for count = 0; count < 4; count++ {
		t1.p.x[count] = p[count]
		t1.p.y[count] = p[count+4]
		t1.p.z[count] = p[count+8]
	}
	if sign == 1 {
		t1.p.Neg(t1.p)
	}

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

		if sel > 0 {
			copy(p, p256Precomputed1[i][(sel-1)*8:(sel-1)*8+8])
		}

		for count = 0; count < 4; count++ {
			t2.p.x[count] = p[count]
			t2.p.y[count] = p[count+4]
			//t2.p.z[count] = p[count+8]
		}
		if sign == 1 {
			gfpNeg(&t2.p.y, &t2.p.y)
		}
		//p256PointAddAffineAsm(p.xyz[0:12], p.xyz[0:12], t0.xyz[0:8], sign, sel, zero)
		if sel != 0 {
			if zero == 0 {
				t1.p.x.Set(&t2.p.x)
				t1.p.y.Set(&t2.p.y)
				t1.p.z = *newGFp(1)

			}
			if zero != 0 {
				t4.x = t2.p.x
				t4.y = t2.p.y
				t1.p.Add_JA(t1.p, t4)
			}
		}
		zero |= sel
	}
	c.Set(t1.p)
}

func p256pre_store(precomp *[16][4 * 3]uint64, index int, in *curvePoint) {
	for count := 0; count < 4; count++ {
		precomp[index][count] = in.x[count]
		precomp[index][4+count] = in.y[count]
		precomp[index][8+count] = in.z[count]
	}
}

func p256Select1(output *curvePoint, precomp *[16][4 * 3]uint64, sel int) {
	var count int
	temp := new(G1).Set(&G1{curveGen})

	if sel > 0 {
		//copy(p, precomp[0][(sel-1)*12:(sel-1)*12+12])
		for count = 0; count < 4; count++ {
			temp.p.x[count] = precomp[sel-1][count]
			temp.p.y[count] = precomp[sel-1][count+4]
			temp.p.z[count] = precomp[sel-1][count+8]
		}
	}
	if sel == 0 {
		temp.p.Set(ginfinite)
	}
	output.Set(temp.p)
}

func splitK_bigint(k *big.Int) (*big.Int, *big.Int, int, int) {
	c1, c2 := new(big.Int), new(big.Int)
	tmp1, tmp2 := new(big.Int), new(big.Int)
	k1, k2 := new(big.Int), new(big.Int)

	// c1 = round(b2 * k / n) from step 4.
	// Rounding isn't really necessary and costs too much, hence skipped
	c1.Mul(b2, k)
	c1.Div(c1, Order)
	// c2 = round(b1 * k / n) from step 4 (sign reversed to optimize one step)
	// Rounding isn't really necessary and costs too much, hence skipped
	c2.Mul(b1, k)
	c2.Div(c2, Order)
	// k1 = k - c1 * a1 - c2 * a2 from step 5 (note c2's sign is reversed)
	tmp1.Mul(c1, a1)
	tmp2.Mul(c2, a2)
	k1.Sub(k, tmp1)
	k1.Add(k1, tmp2)
	// k2 = - c1 * b1 - c2 * b2 from step 5 (note c2's sign is reversed)
	tmp1.Mul(c1, b1)
	tmp2.Mul(c2, b2)
	k2.Sub(tmp2, tmp1)

	return k1, k2, k1.Sign(), k2.Sign()
}

////////

//apply k = k1+k2*beta
func (c *curvePoint) Mul(a *curvePoint, k *big.Int) {
	// precomp is a table of precomputed points that stores powers of p
	// from p^1 to p^16.
	//var precomp [16 * 4 * 3]uint64
	//k = k1+k2*beta
	var precomp1, precomp2 [16][4 * 3]uint64

	ta, enda := &curvePoint{}, &curvePoint{}

	ta.Set(a)

	gfpMul(&enda.x, beta, &a.x)
	enda.y.Set(&a.y)
	enda.z.Set(&a.z)
	enda.t.Set(&a.t)

	k1, k2, signk1, signk2 := splitK_bigint(k)
	//fmt.Println("k1,k2:",k1,"--",k2)

	scalar1 := make([]uint64, 2)
	fromBig(scalar1, k1)

	scalar2 := make([]uint64, 2)
	fromBig(scalar2, k2)

	if signk1 == -1 {
		ta.Neg(ta)
	}
	if signk2 == -1 {
		enda.Neg(enda)
	}

	t0, t1, t2, t3, temp, temp2 := &curvePoint{}, &curvePoint{}, &curvePoint{}, &curvePoint{}, &curvePoint{}, &curvePoint{}

	//----------------- Prepare the table of point a---------------------------
	in := &curvePoint{}
	in.Set(ta)
	p256pre_store(&precomp1, 0, in) //store1

	t0.Double(in) //2G
	t1.Double(t0) //4G
	t2.Double(t1) //8G
	t3.Double(t2) //16G
	p256pre_store(&precomp1, 1, t0)
	p256pre_store(&precomp1, 3, t1)
	p256pre_store(&precomp1, 7, t2)
	p256pre_store(&precomp1, 15, t3)

	temp.Add(in, t0)  //3G
	temp2.Add(in, t1) //5G
	t3.Add(in, t2)    //9G
	p256pre_store(&precomp1, 2, temp)
	p256pre_store(&precomp1, 4, temp2)
	p256pre_store(&precomp1, 8, t3)

	//temp.Add(t0, t1) //6G
	//t3.Add(t0, t2)   //10G
	temp.Double(temp)   //6G
	temp2.Double(temp2) //10G
	p256pre_store(&precomp1, 5, temp)
	p256pre_store(&precomp1, 9, temp2)
	t1.Add(in, temp)  //7G
	t2.Add(in, temp2) //11G
	p256pre_store(&precomp1, 6, t1)
	p256pre_store(&precomp1, 10, t2)
	temp.Double(temp) //12G
	t3.Double(t1)     //14G

	p256pre_store(&precomp1, 11, temp)
	p256pre_store(&precomp1, 13, t3)
	t2.Add(temp, in) //13G
	t1.Add(t3, in)   //15G
	p256pre_store(&precomp1, 12, t2)
	p256pre_store(&precomp1, 14, t1)

	//----------------- Prepare the table of point enda---------------------------
	in.Set(enda)
	p256pre_store(&precomp2, 0, in) //store1

	t0.Double(in) //2G
	t1.Double(t0) //4G
	t2.Double(t1) //8G
	t3.Double(t2) //16G
	p256pre_store(&precomp2, 1, t0)
	p256pre_store(&precomp2, 3, t1)
	p256pre_store(&precomp2, 7, t2)
	p256pre_store(&precomp2, 15, t3)

	temp.Add(in, t0)  //3G
	temp2.Add(in, t1) //5G
	t3.Add(in, t2)    //9G
	p256pre_store(&precomp2, 2, temp)
	p256pre_store(&precomp2, 4, temp2)
	p256pre_store(&precomp2, 8, t3)

	//temp.Add(t0, t1) //6G
	//t3.Add(t0, t2)   //10G
	temp.Double(temp)   //6G
	temp2.Double(temp2) //10G
	p256pre_store(&precomp2, 5, temp)
	p256pre_store(&precomp2, 9, temp2)
	t1.Add(in, temp)  //7G
	t2.Add(in, temp2) //11G
	p256pre_store(&precomp2, 6, t1)
	p256pre_store(&precomp2, 10, t2)
	temp.Double(temp) //12G
	t3.Double(t1)     //14G

	p256pre_store(&precomp2, 11, temp)
	p256pre_store(&precomp2, 13, t3)
	t2.Add(temp, in) //13G
	t1.Add(t3, in)   //15G
	p256pre_store(&precomp2, 12, t2)
	p256pre_store(&precomp2, 14, t1)

	// Start scanning the window from top bit
	index := uint(124)
	var sel, sign int

	wvalue := (scalar1[index/64] >> (index % 64)) & 0x3f
	sel, _ = boothW5(uint(wvalue))

	p256Select1(temp, &precomp1, sel)

	wvalue = (scalar2[index/64] >> (index % 64)) & 0x3f
	sel, _ = boothW5(uint(wvalue))

	p256Select1(temp2, &precomp2, sel)

	t0.Add(temp, temp2)

	for index > 4 {
		index -= 5

		t0.Double(t0)
		t0.Double(t0)
		t0.Double(t0)
		t0.Double(t0)
		t0.Double(t0) //t0 <--2^5t0

		if index < 64 {
			wvalue = ((scalar1[index/64] >> (index % 64)) + (scalar1[index/64+1] << (64 - (index % 64)))) & 0x3f
		} else {
			wvalue = (scalar1[index/64] >> (index % 64)) & 0x3f
		}

		sel, sign = boothW5(uint(wvalue))

		p256Select1(t1, &precomp1, sel)

		if sign == 1 {
			t1.Neg(t1)
		}

		t0.Add(t1, t0) //

		if index < 64 {
			wvalue = ((scalar2[index/64] >> (index % 64)) + (scalar2[index/64+1] << (64 - (index % 64)))) & 0x3f
		} else {
			wvalue = (scalar2[index/64] >> (index % 64)) & 0x3f
		}

		sel, sign = boothW5(uint(wvalue))

		p256Select1(t1, &precomp2, sel)

		if sign == 1 {
			t1.Neg(t1)
		}

		t0.Add(t1, t0) //
	}

	t0.Double(t0)
	t0.Double(t0)
	t0.Double(t0)
	t0.Double(t0)
	t0.Double(t0) //t0 <--2^5t0

	wvalue = (scalar1[0] << 1) & 0x3f
	sel, sign = boothW5(uint(wvalue))

	p256Select1(t1, &precomp1, sel)

	if sign == 1 {
		t1.Neg(t1)
	}

	t0.Add(t1, t0) //

	wvalue = (scalar2[0] << 1) & 0x3f
	sel, sign = boothW5(uint(wvalue))

	p256Select1(t1, &precomp2, sel)
	if sign == 1 {
		t1.Neg(t1)
	}

	t0.Add(t1, t0) //

	c.Set(t0)
}

var g32 = [33][12]uint64{ //G-->32G jacobi format; Note 33 is needed.
	{0x22e935e29860501b, 0xa946fd5e0073282c, 0xefd0cec817a649be, 0x5129787c869140b5, 0xee779649eb87f7c7, 0x15563cbdec30a576, 0x326353912824efbf, 0x7215717763c39828, 0x1a9064d81caeba83, 0xde0d6cb4e5851124, 0x29fc54b00a7138ba, 0x49bffffffd5c590e},
	{0x4b1315758013ba8a, 0xc1766d14e5c40522, 0xc2934a22bcdb3562, 0x2b25f38ad2488ddf, 0xd1a4e6f3617676af, 0x35e083458214cbb, 0x5b50be741c809b69, 0x2adc9853c8e95baf, 0xf77f916bf3beaa11, 0x8b9e630bde65c11, 0x8ec2fbd25abb1839, 0x2deae2eec4e3895e},
	{0x55ca89163d966593, 0xeb61bdf266fad7fe, 0x930cc3c2953ed25b, 0x3b3d7c9b40e5f7cb, 0xa4885d1647a56f19, 0x2f01a68b2c78f57, 0xa16c9eaed4f2fe97, 0x6f296770cf1d9955, 0x72c7b851f4713332, 0xf0397231f1dde474, 0x7cde4cec21a47ff0, 0x9500c2c507d155},
	{0xc7c5b606a9ad4cb7, 0x5752ce48f2fad90a, 0xfc61137b83e48ff2, 0x3dd2479279089fdb, 0x3bf2b296c4960b41, 0x1d54b3b919fb33cb, 0xd621f25e6658a1e0, 0x3dd928a89728075a, 0x1133c5a4eb41e1cc, 0x16b04348e11b283b, 0xbc7b33abaeb1d7f, 0xa518ec079964e1c3},
	{0xacd063dc8b82dc22, 0x3d405779552b1de8, 0x6f341ce0ddc31bb2, 0x21312b65e9e3a053, 0x39627269dd774051, 0x6f94cff3d4ea855c, 0xe546254135ea4710, 0x8a0337f3a56691c5, 0xe1073def6c46f62, 0x7b12c98e22de92cf, 0x15246128f053640f, 0x1c2956aebb8a4ab7},
	{0x80cb385aabd18359, 0xff47a036aea57692, 0xa963c091620c3d7c, 0x2d77216e1a972c0e, 0xfbcf3f7dc3750bf9, 0xb1605dfac3410f6d, 0xb8833403c9f5ed1f, 0x4fa582d15715862f, 0xc4c4be420820573b, 0x954e4b238620e81a, 0xc3b5906a27ee1e08, 0x2d0fc12d78ab5513},
	{0x36a711d0bb47e05a, 0x361c09f5f50d6836, 0x4784f00d78c872a6, 0x3d23c50c347ca9bb, 0x8762c6b78f7f74f0, 0x39fe396d30f24a78, 0x47aba03fceea989c, 0x49b1356bdf203d86, 0x6e40555271b00723, 0x38ae38dad39edee9, 0xda56bf13f955d372, 0x645cce4f4a9bf1d1},
	{0x44949441c9015e1a, 0xdfd9ba223220330b, 0xde7d2586b318fa09, 0x89d49975bf5c116d, 0xbcca9dfb0d844882, 0xb3a4a2054f5afaa, 0x1edf67cd2a2ce709, 0xa041d227906ab04, 0xb53545ee77804ad7, 0x9664cbea81f4b6d0, 0xcf7eadbf7e93cc93, 0x5fda931291b8a1e6},
	{0xe17d57a9f8fa63df, 0x618dca119f98e568, 0x4f03420981d1c310, 0x173574c676ad3235, 0xddd501c345ff8d28, 0x2ea3f364852f7062, 0x47b699ef07ce21cd, 0xade77ee93cd104de, 0x9acfb3264074ab43, 0xcfe98048753944a, 0x771498a60cafc103, 0x4cdf2529c48c2b85},
	{0xb8e5135d7fc89827, 0x4706d5241a131d11, 0xff19d594bbeeca63, 0x7fadc28a3b848268, 0x3b5504b424eb48c5, 0xe70b2949b0eab29e, 0xadf9f583bb53f58, 0x9268cac3c058cdb1, 0x968146aec2c9d241, 0x2f667c336af0a384, 0x9945fe61b1b11739, 0x5142f477baf28187},
	{0xf79b4e207359c3e8, 0x859d8e82097d1fd3, 0xfba85b4a4459a850, 0x390390dfaff45f6d, 0x80e9b1ba48f2bb58, 0x5ece495ca45e7e7a, 0x786586a390b1d640, 0x8795d5c9b9feafec, 0x5e7e2c0cc578326e, 0xa9d08d3933ebe06b, 0x6c75774f9526257e, 0x93afe2137c381d30},
	{0xda5052bbe4803d9, 0x84eff24e72f42559, 0x5019b9e6bcdf64dd, 0x7b472f076a967384, 0x2180741bd3933ef3, 0xdcfe902a9279892e, 0x156d34b6c8ae16a5, 0x2e7b08518d99c479, 0x99030644615bf8be, 0xcc262957d7035f92, 0xb354a75e2e1e39d4, 0xb8dbc7bc3ff1a28},
	{0x3a89d14e5474b63a, 0xd474f4a02c381de2, 0x2ecec43d29079790, 0x1ac52eb0e852f59d, 0x222fcdc5f99c013a, 0xec4808ad0211a424, 0x9ec192af1957a349, 0x14ade8b69c8fa5e3, 0x6193a5e40ca25a9a, 0x442fe8e0083fd4b1, 0xb65da4567e8f60fb, 0xf231d00a2f0cf9d},
	{0xebbb083102317e6, 0xc3a365142d701e78, 0xe8acfec088e7c9bc, 0x190ffb5b0c84ba0c, 0xd1d39181e913efc5, 0x3010b6a3ee6977d3, 0x128f3b80c3cd7846, 0x832b85ae9d053b75, 0x2064057515c012f2, 0x6d165fc40a0227d7, 0x7677e5f116b3a3e5, 0x9cc276699c36b48c},
	{0x67b219218d5966a, 0x290d3bcdb14df588, 0x6278338020e507cb, 0x5abbb37037e95c22, 0x67a63e58bb700364, 0x4d61f457c17b6586, 0xac0feae3076f8c90, 0x5d42f13850e9d4f0, 0xf87dc37f29d39cfb, 0x11727eda3ba868e6, 0xd8ac59434c393565, 0x8f8b85755bac5659},
	{0xfbd81f4322ebff5d, 0x799911b6e64d777a, 0x2f7d30872bed9d14, 0x132d426d0b88f62e, 0x9b17a424bda60838, 0x280e5b4595ba5fd8, 0xf95d2889227d62f4, 0x667c25ef6d65ec2d, 0x95682ebdf3f17c81, 0x952da4be3ddcf492, 0xdbc9d30fe785ba6a, 0xa5c7207807bfa985},
	{0x5ef1fd0d7b3f69f4, 0xe98ecc5bc1566b52, 0x187b3662171f4c10, 0x9155fe311d25379f, 0x9851a7885be9207, 0xd6eb2f9120fb527, 0x367be4efa57b68c3, 0x86866e3c1133897c, 0xa8006df7cd1b6956, 0xcf5be8db2a3de350, 0xd7e2abd10214889b, 0xaff383a649bd1daa},
	{0xca6600fadce0a128, 0x3d5108c2ec21c71, 0xc63b4530b873c71, 0x60b742c35ad735e7, 0x1f5221c3e4d21f9d, 0xa6647b6879b87efc, 0xd4a776b20d79837e, 0xc1d507ca6278782, 0x22d59a03943c3f3d, 0x4e079856bb6e6fc3, 0xd022bf114b40bb5a, 0x916cb328382084fd},
	{0x6e05798d846bb754, 0xefa37d0d5e222e91, 0x5ea28eda3588d2e1, 0xc65563e49331a47, 0xd2463bfcdf8b7701, 0x53b2d737ff08132b, 0x2c10bca4e525dd92, 0x9bed1a5c60f946cd, 0x5e31196168a522ea, 0xc5291ad8c678cc20, 0x11717ee8455808e, 0x5e05d7230ca51bb8},
	{0x27003dd5e61db81e, 0x7a099325722a64f2, 0x7c428d7ec89d2133, 0x71cb37af34a38363, 0x8212188c3a6f2b12, 0xfdfa1fb914925958, 0x277ea3132e3aac1f, 0xa4c33887779aa3f9, 0xe7a60cfbb9d63c6f, 0x177cc2f8f9e281d7, 0x3eb743a5b0d377b2, 0x7ab0feac2f30cf5e},
	{0xe694524bd8c75a36, 0x4a4f4900a90f5faa, 0xeb6ebc036c8ab765, 0x22476b4e25509186, 0xd039f0c073e7b2d6, 0x801b571f6d2aa8b2, 0x8ba7d77f1a4927b7, 0xafb3bc9f5e7ba8a7, 0x8b84e7ef2eb45490, 0xa5cdc51dea6b2a28, 0x233ea9a9ca4f110a, 0x6f559bfddcb5596f},
	{0x799d8b2f3ea582a5, 0x13276fc39e5099cd, 0xe2c17dc08283530a, 0x36fb60b15f86d241, 0x33b3fa8b87ff9e3a, 0x9d41e09c7eb8dcb0, 0x95fcd658df9eb864, 0x64dfa75f9e52e5c4, 0x4e9ac90fc078f04d, 0xe8de2306ddd45921, 0x8c24a059aa2ffe33, 0x8c995ac97a0fe063},
	{0x38eca76390bb113a, 0x2f0b4e54ee6356f5, 0xe3aaca2ffa8d6996, 0x627f1afbc2ee2214, 0x9061b4d4ea3363bb, 0x94d47bdd148e9be4, 0x5994c59683bfe915, 0x567396336d3a2689, 0xe85a21665b54a106, 0x831529fbb411936b, 0x5c04e2e23ebdc807, 0xd7a32537cce3d22},
	{0x408dc8873f0aa53c, 0xd31293fef7b352b4, 0x232d9654ef09e4cb, 0x3ff6ef4312f65ba3, 0xa7afc6f996cd707e, 0xb8291cd8d54d0a12, 0x3719569dc08d80c6, 0x1285435f338cd01c, 0xec437cebc2740e1d, 0x227542827e3c9493, 0xe7a9ef1b5335e092, 0xde97620de838602},
	{0xb113f6e7813cefbf, 0x68d5f1014691f30e, 0x59e60ba7d4d6031b, 0x8ffcfd3896976532, 0xe4ac614a60d2de5f, 0x1f63288f500690f6, 0x2bd59886d99ca3a1, 0xafc1a645ed2aa9f0, 0x54f6a9c0fccb7fcc, 0x884f244f369ead5b, 0x97ec1d2a44d15f98, 0x1c59e630b0b95602},
	{0x9a832df67c1fb660, 0x39b474439d65e5f5, 0x4278a075e13d6660, 0x229d0cb37a2a6db9, 0x634c7db232f8c5ca, 0x5de621787aefde5c, 0xd371d496dad5a127, 0x36e690172d6fdbbf, 0x20e243fcca373cbc, 0x9efd3b7bde0d3c80, 0x8e2cb568d414cfc6, 0x1eaebb6ab29a5410},
	{0x36c373d9497dde0b, 0x1725701832f5cb7d, 0xb820e4f58f90ee81, 0x7ec3f29a480739e2, 0x4a70394b1d5ef22a, 0xd3aec7524385a703, 0x9043abe590c46836, 0xa9b1ef78feb1864, 0xc44e6fd54b2fb679, 0x58311e62cddac409, 0x4e3a0f2387b1feae, 0xcebc79ac74372fc},
	{0x487699a0cb9d15fc, 0x44394aa6554cd8a2, 0x9c20a7e124f3c8e8, 0x850f94a4da1d28c0, 0xb0c04241c55a1878, 0x685315db6f20dcc8, 0x539c0076e5e1b7ae, 0x8e2ae9770752bcf9, 0x3222716fad7be4de, 0xbee702f16e4ab488, 0xe71cdda5d919015f, 0x3055bfaca8e6487d},
	{0x19dbf18c98f6266f, 0x15d4b9069e1afc4d, 0x5c6c7bb2e0508df6, 0xa2390c76b6ef863f, 0xf157b84106a1d0c2, 0x4e3799701a85d4bc, 0x439a35f9ff8cc096, 0x2dffd18a55c386c5, 0x4bcb8806fe9681f3, 0xd9570902282f5700, 0x7f2aefe1f42fe0fb, 0x77262a4ad9b62855},
	{0x18e33b64dd0917c6, 0xc52702e49e8a0df9, 0x2ce22ca6a5f9b55c, 0xa6cd01f5361b165c, 0x9460528369141c4f, 0x6241262634305079, 0xf9a56860d9746a4, 0x5e7f283f3ac834e1, 0xcb209ef3f4caa077, 0x47aa21c6b089d4ef, 0x85358abdb53abc7d, 0xc94d13e3916f2ac},
	{0xfc3a1c6ba3c00435, 0x90f330ae7f04a58c, 0x7a4a21aa35cf3083, 0x3eb91e6fa015b8bd, 0xf6d3efa9727e5a24, 0xdb006b99294e8a42, 0x5856b38bb14e0b62, 0x5b44c5d4c85ce28, 0x57b144b4133647c5, 0xbce93e0f3833f943, 0x736cc1f130361d9d, 0x8b694b38d769a82b},
	{0x87fbeec59a535668, 0xddedaa36aa5710b4, 0xabba13eeb88c0569, 0x874b6507e1cb3f32, 0x2eb38c7694e903d3, 0xf218742f104a2503, 0xf1975832862bfb34, 0x461d9fec046a16de, 0x4e57f11be9880caa, 0x79c4d66549fdc53, 0x86010bd6bebeaf74, 0x7ddd7c66ca954114}}
