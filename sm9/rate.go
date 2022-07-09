package sm9

func lineFunctionAdd(r, p *twistPoint, q *curvePoint, r2 *gfP2) (a, b, c *gfP2, rOut *twistPoint) {
	// See the mixed addition algorithm from "Faster Computation of the
	// Tate Pairing", http://arxiv.org/pdf/0904.0854v3.pdf
	B := (&gfP2{}).Mul(&p.x, &r.t)
	D := (&gfP2{}).Add(&p.y, &r.z)
	D.Square(D)
	D.Sub(D, r2)
	D.Sub(D, &r.t)
	D.Mul(D, &r.t)

	H := (&gfP2{}).Sub(B, &r.x)
	I := (&gfP2{}).Square(H)

	E := (&gfP2{}).Add(I, I)
	E.Add(E, E)

	J := (&gfP2{}).Mul(H, E)

	L1 := (&gfP2{}).Sub(D, &r.y)
	L1.Sub(L1, &r.y)

	V := (&gfP2{}).Mul(&r.x, E)

	rOut = &twistPoint{} //siomn:2R part
	rOut.x.Square(L1)
	rOut.x.Sub(&rOut.x, J)
	rOut.x.Sub(&rOut.x, V)
	rOut.x.Sub(&rOut.x, V)

	rOut.z.Add(&r.z, H)
	rOut.z.Square(&rOut.z)
	rOut.z.Sub(&rOut.z, &r.t)
	rOut.z.Sub(&rOut.z, I)

	t := (&gfP2{}).Sub(V, &rOut.x)
	t.Mul(t, L1)
	t2 := (&gfP2{}).Mul(&r.y, J)
	t2.Add(t2, t2)
	rOut.y.Sub(t, t2)
	rOut.t.Square(&rOut.z) //R=P+Q end;

	//cf version
	t.Add(&p.y, &rOut.z).Square(t).Sub(t, r2).Sub(t, &rOut.t)

	t2.Mul(L1, &p.x)
	t2.Add(t2, t2)
	a = (&gfP2{}).Sub(t2, t)

	c = (&gfP2{}).MulScalar(&rOut.z, &q.y)
	c.Add(c, c)
	c.MulXi(c) //simon add1030;first test success;

	b = (&gfP2{}).Neg(L1)
	b.MulScalar(b, &q.x).Add(b, b)

	return
}

func lineFunctionDouble(r *twistPoint, q *curvePoint) (a, b, c *gfP2, rOut *twistPoint) {
	// See the doubling algorithm for a=0 from "Faster Computation of the
	// Tate Pairing", http://arxiv.org/pdf/0904.0854v3.pdf

	A := (&gfP2{}).Square(&r.x)
	B := (&gfP2{}).Square(&r.y)
	C := (&gfP2{}).Square(B)

	D := (&gfP2{}).Add(&r.x, B)
	D.Square(D)
	D.Sub(D, A)
	D.Sub(D, C)
	D.Add(D, D)

	E := (&gfP2{}).Add(A, A)
	E.Add(E, A) //E = 3A

	G := (&gfP2{}).Square(E) //G = E^2

	rOut = &twistPoint{} //X3 = G-2D
	rOut.x.Sub(G, D)
	rOut.x.Sub(&rOut.x, D)

	rOut.z.Add(&r.y, &r.z) //Z3 = (Y1+Z1)^2-B-T1;
	rOut.z.Square(&rOut.z)
	rOut.z.Sub(&rOut.z, B)
	rOut.z.Sub(&rOut.z, &r.t)

	rOut.y.Sub(D, &rOut.x) //Y3 = E(D-X3)-8C
	rOut.y.Mul(&rOut.y, E)
	t := (&gfP2{}).Add(C, C)
	t.Add(t, t)
	t.Add(t, t)
	rOut.y.Sub(&rOut.y, t) //t =8C

	rOut.t.Square(&rOut.z) // T3 = Z3^2; R = 2R end

	//bn256 version//1030first test success;speed up 25%
	t.Mul(E, &r.t) //
	t.Add(t, t)
	b = &gfP2{}
	b.SetZero()
	b.Sub(b, t)
	b.MulScalar(b, &q.x) //b = -2E*T1*q.x;

	a = &gfP2{}
	a.Add(&r.x, E) //a = (X1 +E)^2 -A-G-4B; // a(W^3)
	a.Square(a)
	a.Sub(a, A)
	a.Sub(a, G)
	t.Add(B, B)
	t.Add(t, t)
	a.Sub(a, t)

	c = &gfP2{}
	c.Mul(&rOut.z, &r.t)
	c.Add(c, c)
	c.MulScalar(c, &q.y) //c = 2*Z3*T1*q.y
	c.MulXi(c)

	return
}

//simon: 0928change this mulLine;
//0929version
func mulLine(ret *gfP12, a, b, c *gfP2) {
	a2 := &gfP6{}
	a3 := &gfP6{}
	a12 := &gfP12{}

	a2.z.SetZero()
	a2.y.Set(a)
	a2.x.Set(b)

	//c.MulXi(c, pool) //(0,0,c)
	a3.z.Set(c)
	a3.y.SetZero()
	a3.x.SetZero()

	a12.x.Set(a2) //test
	a12.y.Set(a3)
	ret.Mul(ret, a12)

}

// sixuPlus2NAF is 6u+2 in non-adjacent form.//NAF is checked again in 0918;but the pair result is wrong!
//var sixuPlus2NAF1 = []int8{0, -1, 0, 0, 0, 0, 1, 0, 1, 0, 0, -1, 0, -1, 0,
//	0, 0, -1, 0, -1, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
//	0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1}

var sixuPlus2NAF = [66]int8{0, 1, 1, 1, 1, 1, 0, 0, 1, 0, 0, 1, 1, 0, 1, 1,
	1, 0, 1, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0,
	0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
	0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1}

// miller implements the Miller loop for calculating the Optimal Ate pairing.
// See algorithm 1 from http://cryptojedi.org/papers/dclxvi-20100714.pdf
func miller(q *twistPoint, p *curvePoint) *gfP12 {
	ret := &gfP12{} //ret = f
	ret.SetOne()
	//fmt.Println("In miller:")

	aAffine := &twistPoint{}
	aAffine.Set(q)
	aAffine.MakeAffine()

	bAffine := &curvePoint{}
	bAffine.Set(p)
	bAffine.MakeAffine()

	minusA := &twistPoint{}
	minusA.Neg(aAffine)

	r := &twistPoint{} // R<--Q
	r.Set(aAffine)

	r2 := (&gfP2{})
	r2.Square(&aAffine.y) //r2 <-- Q.y^2
	//fmt.Println("loop:")
	for i := len(sixuPlus2NAF) - 2; i >= 0; i-- {
		a, b, c, newR := lineFunctionDouble(r, bAffine) // lR,R(P); R <--2R

		if i != len(sixuPlus2NAF)-1 {
			ret.Square(ret)
		}

		mulLine(ret, a, b, c)

		r = newR

		switch sixuPlus2NAF[i] {
		case 1:
			a, b, c, newR = lineFunctionAdd(r, aAffine, bAffine, r2)
		case -1:
			a, b, c, newR = lineFunctionAdd(r, minusA, bAffine, r2)
		default:
			continue
		}

		mulLine(ret, a, b, c)
		r = newR
	}

	// In order to calculate Q1 we have to convert q from the sextic twist
	// to the full GF(p^12) group, apply the Frobenius there, and convert
	// back.
	//
	// The twist isomorphism is (x', y') -> (xω², yω³). If we consider just
	// x for a moment, then after applying the Frobenius, we have x̄ω^(2p)
	// where x̄ is the conjugate of x. If we are going to apply the inverse
	// isomorphism we need a value with a single coefficient of ω² so we
	// rewrite this as x̄ω^(2p-2)ω². ξ⁶ = ω and, due to the construction of
	// p, 2p-2 is a multiple of six. Therefore we can rewrite as
	// x̄ξ^((p-1)/3)ω² and applying the inverse isomorphism eliminates the
	// ω².
	//
	// A similar argument can be made for the y value.

	q1 := &twistPoint{}
	q1.x.Conjugate(&aAffine.x) //q1.x = 共轭(q.x)
	q1.x.Mul(&q1.x, xiToPMinus1Over3i)
	q1.y.Conjugate(&aAffine.y)
	q1.y.Mul(&q1.y, xiToPMinus1Over2i)
	q1.z.SetOne()
	q1.t.SetOne() //simon:0926 test pass;

	// For Q2 we are applying the p² Frobenius. The two conjugations cancel
	// out and we are left only with the factors from the isomorphism. In
	// the case of x, we end up with a pure number which is why
	// xiToPSquaredMinus1Over3 is ∈ GF(p). With y we get a factor of -1. We
	// ignore this to end up with -Q2.

	minusQ2 := &twistPoint{}
	minusQ2.x.MulScalar(&aAffine.x, xiToPSquaredMinus1Over3i)
	minusQ2.y.Set(&aAffine.y)
	minusQ2.z.SetOne()
	minusQ2.t.SetOne()

	r2.Square(&q1.y)
	a, b, c, newR := lineFunctionAdd(r, q1, bAffine, r2) //step11
	mulLine(ret, a, b, c)
	r = newR

	r2.Square(&minusQ2.y)
	a, b, c, newR = lineFunctionAdd(r, minusQ2, bAffine, r2) //step12
	mulLine(ret, a, b, c)

	r = newR

	return ret
}

// finalExponentiation computes the (p¹²-1)/Order-th power of an element of
// GF(p¹²) to obtain an element of GT (steps 13-15 of algorithm 1 from
// http://cryptojedi.org/papers/dclxvi-20100714.pdf)
func finalExponentiation(in *gfP12) *gfP12 {
	t1 := &gfP12{}

	// This is the p^6-Frobenius
	t1.x.Neg(&in.x)
	t1.y.Set(&in.y)

	inv := &gfP12{}
	inv.Invert(in)
	t1.Mul(t1, inv) //f^(p^6-1);//simon check ok;0930

	t2 := (&gfP12{}).FrobeniusP2(t1)
	t1.Mul(t1, t2) //simon:(f^(p^6-1))^(p^2+1) 0920

	//simon1031version:
	st0 := (&gfP12{}).Frobenius(t1)
	sx0 := (&gfP12{}).Frobenius(st0)
	sx1 := (&gfP12{}).Mul(t1, st0)

	sx0 = sx0.Mul(sx0, sx1).Frobenius(sx0)
	sx1 = sx1.Invert(t1) //simon check ok 1031;

	sx4 := (&gfP12{})
	//sx4.Exp(t1, tneg)
	sx4 = sx4.Exp(t1, t)
	sx4.Invert(sx4)

	sx3 := (&gfP12{}).Frobenius(sx4)
	//sx2 := (&gfP12{}).Exp(sx4, tneg)
	sx2 := (&gfP12{}).Exp(sx4, t)
	sx2.Invert(sx2)
	sx5 := (&gfP12{}).Invert(sx2)

	//st0.Exp(sx2, tneg)
	st0.Exp(sx2, t)
	st0.Invert(st0)

	sx2.Frobenius(sx2)
	sx2inv := (&gfP12{}).Invert(sx2)
	sx4.Mul(sx4, sx2inv)
	sx2.Frobenius(sx2)

	t1.Frobenius(st0)
	st0.Mul(st0, t1)

	st0.Mul(st0, st0)
	st0.Mul(st0, sx4)
	st0.Mul(st0, sx5)
	t1.Mul(sx3, sx5)
	t1.Mul(t1, st0)

	st0.Mul(st0, sx2)
	t1.Mul(t1, t1)
	t1.Mul(t1, st0)

	t1.Mul(t1, t1)
	st0.Mul(t1, sx1)
	t1.Mul(t1, sx0)

	st0.Mul(st0, st0)
	st0.Mul(st0, t1)
	//1031version end;

	return st0
}

func optimalAte(a *twistPoint, b *curvePoint) *gfP12 {
	e := miller(a, b)
	//t1 := time.Now()
	ret := finalExponentiation(e)
	//elapsed := time.Since(t1)
	//fmt.Println("final spend time: ", elapsed)

	if a.IsInfinity() || b.IsInfinity() {
		ret.SetOne()
	}

	return ret
}
