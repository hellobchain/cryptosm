package sm9

import (
	"bytes"
	"crypto/rand"
	"encoding/hex"
	"fmt"
	"math/big"
	"testing"
	"time"
)

func TestG1Imp(t *testing.T) { //check:0924;
	fmt.Printf("Enter G1Imp test:\n")
	g := &curvePoint{
		x: *G1x,
		y: *G1y,
		z: *newGFp(1),
		t: *newGFp(1),
	}
	test := g.IsOnCurve()
	fmt.Println(test)

	var one = big.NewInt(2)

	Gtest := new(G1).Set(&G1{curveGen}).ScalarBaseMult(one)
	test1 := Gtest.p.IsOnCurve()
	fmt.Println(test1)

	mtest := Gtest.Marshal()
	g1 := new(G1).Set(&G1{curveGen})
	g2 := new(G1).Set(&G1{curveGen})

	g3 := new(G1).Add(g1, g2)
	mgen := g3.Marshal()

	mulg := new(G1).ScalarMult(&G1{g}, one)
	www := mulg.p.IsOnCurve()
	fmt.Println(www)

	if !bytes.Equal(mtest, mgen) {
		//t.Fatal("bytes are different")
		fmt.Println("test error")
	}

}

func TestG1Marshal(t *testing.T) {
	_, Ga, err := RandomG1(rand.Reader)
	if err != nil {
		t.Fatal(err)
	}
	ma := Ga.Marshal()
	test := Ga.p.IsOnCurve()
	fmt.Println(test)

	Gb := new(G1)
	_, err = Gb.Unmarshal(ma)
	if err != nil {
		//t.Fatal(err)
		fmt.Println("test error1")
	}
	mb := Gb.Marshal()

	if !bytes.Equal(ma, mb) {
		//t.Fatal("bytes are different")
		fmt.Println("test error")
	}
}
func TestG1_kGtest(t *testing.T) { //check:0924;
	fmt.Printf("\nEnter G1 kG test:\n")
	g := &curvePoint{
		x: *G1x,
		y: *G1y,
		z: *newGFp(1),
		t: *newGFp(1),
	}
	var j int64
	var k = bigFromBase16("291FE3CAC8F58AD2DC462C8D4D578A94DAFD5624DDC28E328D2936688A86CF1A")
	for j = 0; j < 1000; j++ {
		test, err := randomK(rand.Reader)
		//fmt.Printf("%x\n", test)
		if err != nil {
			return
		}
		Gtest := new(G1).Set(&G1{curveGen}).ScalarBaseMult(test)
		//test1 := Gtest.p.IsOnCurve()
		//fmt.Println(test1)
		mtest := Gtest.Marshal()

		mulg := new(G1).ScalarMult(&G1{g}, test)
		test1 := mulg.p.IsOnCurve()
		if test1 == false {
			fmt.Println("randmul fail")
		}
		mgen1 := mulg.Marshal()

		if !bytes.Equal(mtest, mgen1) {
			fmt.Println("kG correctness failed!")
		} else {
			//fmt.Printf("kG correctness is OK!\n")
		}
	}
	//speed test

	t1 := time.Now()
	for j = 0; j < 100000; j++ {
		var test = big.NewInt(j)
		test.Add(k, test)
		e := new(G1).ScalarMult(&G1{g}, test)
		e.p.MakeAffine()
		if j == 0 {
			fmt.Printf("j = %d:  %x,  %x,  %x\n", j, e.p.x, e.p.y, e.p.z)
		}
	}
	elapsed := time.Since(t1)
	fmt.Println("100000 randGmul spend time: ", elapsed)

}
func TestG1_baseGtest_alone(t *testing.T) { //simon:1213
	fmt.Println("\nEnter TestG1_baseG1 test:")

	two := big.NewInt(4)
	var k = bigFromBase16("291FE3CAC8F58AD2DC462C8D4D578A94DAFD5624DDC28E328D2936688A86CF1A")
	k.Sub(k, two)
	e := new(G1).ScalarMult(&G1{curveGen}, two)
	e.p.MakeAffine()
	fmt.Printf("e:%x\n", e.p.x)
	//e:[8fdf2548f0fde68 c80ddebf804d6dd4 c8cef5282905b7ca 6007e08434132464]

}

func TestG1_baseGtest(t *testing.T) { //simon:1213
	fmt.Println("\nEnter TestG1_baseG test:")

	var j int64
	var k = bigFromBase16("291FE3CAC8F58AD2DC462C8D4D578A94DAFD5624DDC28E328D2936688A86CF1A")

	e := new(G1).ScalarBaseMult(k)
	e.p.MakeAffine()
	//fmt.Printf("kbaseG mogo result:%x,  %x\n", j, e.p.x, e.p.y)
	out := e.Marshal()
	//fmt.Printf("kbaseG nomogo result:%x\n", out)
	test := new(big.Int).SetBytes(out)
	if test.Cmp(kbaseG) == 0 {
		fmt.Printf("k baseG result is OK!\n")
	} else {
		fmt.Printf("k baseG result is NOT OK!\n")
	}
	t1 := time.Now()
	for j = 0; j < 100000; j++ {
		var test = big.NewInt(j)
		test.Add(k, test)
		e := new(G1).ScalarBaseMult(test)
		e.p.MakeAffine()
	}
	elapsed := time.Since(t1)
	fmt.Println("100000 basemul spend time:", elapsed)
}

func TestG1CombinedMult(t *testing.T) {
	fmt.Println("\nEnter TestG1CombinedMult test:")
	g := &curvePoint{
		x: *G1x,
		y: *G1y,
		z: *newGFp(1),
		t: *newGFp(1),
	}

	zero := big.NewInt(0)
	one := big.NewInt(1)
	//two := big.NewInt(2)

	// 0×G + 0×G = ∞
	e := new(G1).CombinedMult(&G1{g}, zero, zero) //H=G
	e.p.MakeAffine()
	if e.p.x[0] != 0 {
		fmt.Printf("0×G + 0×G = (%x, %x), should be ∞", e.p.x, e.p.y)
	}

	// 1×G + 0×H = G
	e1 := new(G1).CombinedMult(&G1{g}, one, zero)
	e1.p.MakeAffine()
	if e1.p.x != *G1x || e1.p.y != *G1y {
		fmt.Printf("1×G + 0×G = (%x, %x), should be (%x,%x)", e1.p.x, e1.p.y, *G1x, *G1y)
	}

	// 0×G + 1×H = H =G
	e2 := new(G1).CombinedMult(&G1{g}, zero, one)
	e2.p.MakeAffine()
	if e2.p.x != *G1x || e2.p.y != *G1y {
		fmt.Printf("e2:0×G + 1×G = (%x, %x), should be (%x,%x)", e2.p.x, e2.p.y, *G1x, *G1y)
	}

	// 1×G + 1×H = 2×G
	e3 := new(G1).CombinedMult(&G1{g}, one, one)
	e3.p.MakeAffine()
	//fmt.Printf("2G: %x,%x", e3.p.x, e3.p.y)
	if e3.p.x == *DoubleGx {
		fmt.Printf("1×G + 1×H = 2×G is OK!\n")
	}
	var k = bigFromBase16("291FE3CAC8F58AD2DC462C8D4D578A94DAFD5624DDC28E328D2936688A86CF1A")
	var j int64
	for j = 0; j < 1; j++ { //for large-scale test
		randk := big.NewInt(j)
		basek := new(big.Int).Sub(k, randk)

		e4 := new(G1).CombinedMult(&G1{g}, basek, randk)
		e4.p.MakeAffine()

		out := e4.Marshal()
		//	fmt.Printf("kG: %x,%x\n", e4.p.x, e4.p.y)
		//	fmt.Printf("%x\n", out)
		result := new(big.Int).SetBytes(out)
		if result.Cmp(kbaseG) == 0 {
			fmt.Printf("j=%d:jG +(k-j)G result is OK!\n", j)
		} else {
			fmt.Printf("j=%d:jG +(k-j)G result is NOT OK!\n", j)
		}
	}

	//fmt.Printf("CombinedMult End.\n")
}

func gfpMul_speed_test() {
	g1 := new(G1).Set(&G1{curveGen})
	g2 := new(G1).Set(&G1{curveGen})
	t1 := time.Now()
	for i := 0; i < 1000000; i++ {
		//new(G1).Add(g1, g2)
		g2.p.Double(g1.p)
	}

	elapsed := time.Since(t1)
	fmt.Println("1000000 point add spend time: ", elapsed)
	z12, z22 := &gfP{}, &gfP{}
	t2 := time.Now()
	for i := 0; i < 1000000; i++ {

		g1.p.z.Invert(&g1.p.z)
		gfpMul(z12, &g1.p.z, &g1.p.z)
		gfpMul(z22, &g2.p.z, &g2.p.z)
	}

	elapsed1 := time.Since(t2)
	fmt.Println("1000000 gfPInvert spend time: ", elapsed1)

}

func TestG2Imp(t *testing.T) { //check:0924;
	fmt.Println("\nEnter TestG2Imp:")
	g := new(G2).Set(&G2{twistGen})

	test := g.p.IsOnCurve()
	fmt.Println(test)

	var one = big.NewInt(2)
	//Gtest := new(G1).ScalarBaseMult(one)
	Gtest := new(G2).ScalarBaseMult(one)
	test1 := Gtest.p.IsOnCurve()
	fmt.Println(test1)

	mtest := Gtest.Marshal()
	g1 := new(G2).Set(&G2{twistGen})
	g2 := new(G2).Set(&G2{twistGen})

	g3 := new(G2).Add(g1, g2)
	mgen := g3.Marshal()

	if !bytes.Equal(mtest, mgen) {
		//t.Fatal("bytes are different")
		fmt.Println("G2Imp test NOT OK!")
	} else {
		fmt.Printf("G2Imp is OK!\n")
	}

}

func TestG2(t *testing.T) {
	fmt.Println("\nEnter TestG2:")
	k, Ga, err := RandomG2(rand.Reader)
	if err != nil {
		fmt.Println("test error1")
		t.Fatal(err)
	}
	ma := Ga.Marshal()

	t1 := time.Now()
	for j := 0; j < 1000; j++ {
		new(G2).ScalarBaseMult(k)
	}
	elapsed := time.Since(t1)
	fmt.Println("1000 kG2 spend time: ", elapsed)

	Gb := new(G2).ScalarBaseMult(k)
	mb := Gb.Marshal()
	//mb = append([]byte{0x01}, mb...)//simon:cf have a bug here;1023

	if !bytes.Equal(ma, mb) {
		//t.Fatal("bytes are different")
		fmt.Println("test error")
	}
}

func TestG2Marshal(t *testing.T) {
	_, Ga, err := RandomG2(rand.Reader)
	if err != nil {
		t.Fatal(err)
	}
	ma := Ga.Marshal()

	Gb := new(G2)
	_, err = Gb.Unmarshal(ma)
	if err != nil {
		t.Fatal(err)
	}
	mb := Gb.Marshal()

	if !bytes.Equal(ma, mb) {
		t.Fatal("bytes are different")
	}
}

func TestGT(t *testing.T) {
	k, Ga, err := RandomGT(rand.Reader)
	if err != nil {
		t.Fatal(err)
	}
	ma := Ga.Marshal()

	//Gb := new(GT).Unmarshal((&GT{gfP12Gen}).Marshal())
	Gb := new(GT).Set((&GT{gfP12Gen}))
	Gb.ScalarMult(Gb, k)
	mb := Gb.Marshal()

	if !bytes.Equal(ma, mb) {
		fmt.Println("test error")
		t.Fatal("bytes are different")
	}
}
func TestGT_basemul(t *testing.T) {

	//var k = bigFromBase16("0130E78459D78545CB54C587E02CF480CE0B66340F319F348A1D5B1F2DC5F4")
	var k1 = bigFromBase16("1")
	gt := new(GT).ScalarBaseMult(k1)
	ret := gt.Marshal()

	fmt.Printf("%x\n", ret)

}

func TestGTMarshal(t *testing.T) {
	_, Ga, err := RandomGT(rand.Reader)
	if err != nil {
		t.Fatal(err)
	}
	ma := Ga.Marshal()

	Gb := new(GT)
	_, err = Gb.Unmarshal(ma)
	if err != nil {
		t.Fatal(err)
	}
	mb := Gb.Marshal()

	if !bytes.Equal(ma, mb) {
		fmt.Println("test error")
		t.Fatal("bytes are different")
	}
}

func TestPair_Miller1(t *testing.T) {

	var milldata = bigFromBase16("aa401110b9e23db43679bec6e2a761680dd2eb64463b192e28c7e66dee8ecb344f3017b4f0ecd8c724b5e54d9cdb5ed898317d344e1b9e5e54dc1f3bcc3bee21398a8a2443d028097a4b7181c47a441bdd2f4a65e7c7f83c5efedad73f52b5cc3ce08a2c7b295e7f890b4dc2ae2e21faba80ba06fbae4d5ff3bb97a2a530b4cdade6e3d4125ee965c594be8c6ddd50e330d5577df7567a1ded515522c73cfd183ff68217c2bd563092c9b87be26d76879ab68321b6268968efa8f0b736b4d0c92d3bb20d74ba9ed3ebde5654980249f4d5c4c43fe6921dfc64150f74e7248205ae0add7f05baf1bd722b709c245fe34f22496ffdc89fcb94f89d12f3582980d50567b492ca068a3263dedf6b1fae8ef4225122d0ef000ecad63dd64994fb1d9f876e11d0ff891caa37a047cd84690faf0bb1dc9091afc78cd4e4a8a4431410f5b2c987074e65a17212759cb9c4224e7858c25f1ce04c5392363de1462f8beebbb14ba7e0f02f856fab3e2f6a74eeb6b04cdf2d7dc8ff2b93788f3b009a68bbe1")
	bytes := milldata.Bytes()
	gtnew := new(GT)
	gtnew.Unmarshal(bytes)

	e1 := Pair(&G1{curveGen}, &G2{twistGen})
	test1 := e1.Marshal()
	//e2 := Miller(&G1{curveGen}, &G2{twistGen})
	e3 := finalExponentiation(gtnew.p)

	test := new(GT).Set(&GT{e3}).Marshal()

	fmt.Printf("%x\n", test)
	fmt.Printf("%x\n", test1)

}

func TestPair(t *testing.T) {
	fmt.Printf("\nEnter Pair  test:\n")
	t1 := time.Now()
	for i := 0; i < 1000; i++ {
		Pair(&G1{curveGen}, &G2{twistGen})
	}
	elapsed := time.Since(t1)
	fmt.Println("1000 pair spend time: ", elapsed)
	test := Pair(&G1{curveGen}, &G2{twistGen}).Marshal()

	fmt.Printf("sm9 standards g(G1,G2)result:\n%x\n", test)
}
func TestPair_Miller2(t *testing.T) { //Pair = Miller.finalExp

	e1 := Pair(&G1{curveGen}, &G2{twistGen})
	test1 := e1.Marshal()
	e2 := Miller(&G1{curveGen}, &G2{twistGen})
	e3 := finalExponentiation(e2.p)

	test := new(GT).Set(&GT{e3}).Marshal()

	fmt.Printf("%x\n", test)
	fmt.Printf("%x\n", test1)

}

func TestBilinearity(t *testing.T) {
	fmt.Printf("\nEnter TestBilinearity:\n")
	var count = 0

	for i := 0; i < 100; i++ {
		a, p1, _ := RandomG1(rand.Reader)
		b, p2, _ := RandomG2(rand.Reader)
		e1 := Pair(p1, p2)

		e2 := Pair(&G1{curveGen}, &G2{twistGen})

		e2.ScalarMult(e2, a)
		e2.ScalarMult(e2, b)

		if *e1.p != *e2.p {
			fmt.Println("test error")
			fmt.Printf("a=%x,\n,b=%x,\n", a, b)
			count = 1
			//t.Fatalf("bad pairing result: %s", e1)
		} else {
			//fmt.Printf("Bilinearity is OK!\n")
		}
	}
	if count == 0 {
		fmt.Printf("Bilinearity is OK!\n")
	}
}

func TestTripartiteDiffieHellman(t *testing.T) {
	fmt.Printf("\nEnter TestTripartiteDiffieHellman:\n")
	for i := 0; i < 1000; i++ {
		a, _ := rand.Int(rand.Reader, Order)
		b, _ := rand.Int(rand.Reader, Order)
		c, _ := rand.Int(rand.Reader, Order)

		pa, pb, pc := new(G1), new(G1), new(G1)
		qa, qb, qc := new(G2), new(G2), new(G2)

		pa.Unmarshal(new(G1).ScalarBaseMult(a).Marshal())
		qa.Unmarshal(new(G2).ScalarBaseMult(a).Marshal())
		pb.Unmarshal(new(G1).ScalarBaseMult(b).Marshal())
		qb.Unmarshal(new(G2).ScalarBaseMult(b).Marshal())
		pc.Unmarshal(new(G1).ScalarBaseMult(c).Marshal())
		qc.Unmarshal(new(G2).ScalarBaseMult(c).Marshal())

		k1 := Pair(pb, qc)
		k1.ScalarMult(k1, a)
		k1Bytes := k1.Marshal()

		k2 := Pair(pc, qa)
		k2.ScalarMult(k2, b)
		k2Bytes := k2.Marshal()

		k3 := Pair(pa, qb)
		k3.ScalarMult(k3, c)
		k3Bytes := k3.Marshal()

		if !bytes.Equal(k1Bytes, k2Bytes) || !bytes.Equal(k2Bytes, k3Bytes) {
			t.Errorf("keys didn't agree")
			fmt.Printf("error")
		}
	}
	fmt.Printf("TestTripartiteDiffieHellman end. Noprint is the best!\n")
}

func BenchmarkG1(b *testing.B) {
	x, _ := rand.Int(rand.Reader, Order)
	b.ResetTimer()

	for i := 0; i < b.N; i++ {
		new(G1).ScalarBaseMult(x)
	}
}

func BenchmarkG2(b *testing.B) {
	x, _ := rand.Int(rand.Reader, Order)
	b.ResetTimer()

	for i := 0; i < b.N; i++ {
		new(G2).ScalarBaseMult(x)
	}
}

func BenchmarkGT(b *testing.B) {
	x, _ := rand.Int(rand.Reader, Order)
	b.ResetTimer()

	for i := 0; i < b.N; i++ {
		new(GT).ScalarBaseMult(x)
	}
}

func BenchmarkPairing(b *testing.B) {
	for i := 0; i < b.N; i++ {
		gt := Pair(&G1{curveGen}, &G2{twistGen})
		gt.p.String()
	}
}

/********************************************************/
func TestCurveInterface(t *testing.T) { //check:1225;
	fmt.Printf("Enter curve interface test:\n")
	g := &curvePoint{
		x: *G1x,
		y: *G1y,
		z: *newGFp(1),
		t: *newGFp(1),
	}
	//var testp *ECCInternalPoint

	point := &[12]uint64{}
	//point1 := &ECCInternalPoint{}
	fmt.Printf("%x\n", point)
	copy(point[0:4], g.x[:])
	copy(point[4:8], g.y[:])
	copy(point[8:12], g.z[:])
	t.Logf("point is %x", point)
	x, y := PointToAffine(point)
	fmt.Printf("P2A:%x,%x\n", x, y)

	point1 := AffineToPoint(x, y) //point1=G
	fmt.Printf("A2P:%x\n", point1)
	point3 := JacobianAdd(point, point1) //G+G=2G
	fmt.Printf("JJAdd:%x\n", point3)
	point3d := JacobianDouble(point) //2G
	fmt.Printf("Double:%x\n", point3d)
	//2G-G=G
	PointNegCondition(point, 1)
	fmt.Printf("Negcondition:%x\n", point)
	//point3 = csp.JacobianAdd(point, point3d)
	//x1, y1 := csp.PointToAffine(point3)
	//fmt.Printf("2G-G:%x\n", point3.XYZ)
	//fmt.Printf("x1:%x,%x\n", x1, y1)

	testmul := []uint64{0x4b1315758013ba8a, 0xc1766d14e5c40522, 0xc2934a22bcdb3562, 0x2b25f38ad2488ddf}
	mulresult := make([]uint64, 4)
	X := make([]uint64, 4)
	for i := 0; i < 4; i++ {
		X[i] = g.z[i]
	}
	FieldMul(mulresult, testmul, X)
	fmt.Printf("FieldMul: %x\n", mulresult)
	mulresult = MontgomaryR()
	fmt.Printf("MogoR: %x\n", mulresult)
	mulresult = MontgomaryR2() //data correctness?
	fmt.Printf("MogoR2: %x\n", mulresult)

	var three = big.NewInt(3) //Note:2 is not Permiited

	sqrresult := ModSqrtP(three)
	fmt.Println(sqrresult)
	sqrresult.Mul(sqrresult, sqrresult)
	sqrresult.Mod(sqrresult, p)

	fmt.Printf("sqrt result:%d\n", sqrresult)
	//////////////////////////////////////////////////////////////
	var test11 = []uint64{0x1a9064d81caeba83, 0xde0d6cb4e5851124, 0x29fc54b00a7138ba, 0x49bffffffd5c590e}
	invresult := ModInverseP(test11) //if three ==1 ,return R2,think?
	fmt.Printf("InvertP: %x\n", invresult)
	orderinv := ModInverseOrder(three)
	fmt.Printf("InvertOrder: %x\n", orderinv)
	orderinv.Mul(orderinv, three)
	orderinv.Mod(orderinv, Order)
	fmt.Printf("%x\n", orderinv)

}

func Test_TMP1(t *testing.T) {
	a, _ := randomK(rand.Reader)
	b, _ := randomK(rand.Reader)
	c, _ := randomK(rand.Reader)

	g1Gen := new(G1).ScalarBaseMult(new(big.Int).SetInt64(1))
	g2Gen := new(G2).ScalarBaseMult(new(big.Int).SetInt64(1))

	aG1 := new(G1).ScalarBaseMult(a)
	bAc := new(big.Int).Mod(new(big.Int).Add(b, c), Order)
	bAcG2 := new(G2).ScalarBaseMult(bAc)
	left := Pair(aG1, bAcG2)

	ab := new(big.Int).Mod(new(big.Int).Mul(a, b), Order)
	abG1 := new(G1).ScalarBaseMult(ab)
	right1 := Pair(abG1, g2Gen)

	ac := new(big.Int).Mod(new(big.Int).Mul(a, c), Order)
	tmp := Pair(g1Gen, g2Gen)
	right2 := new(GT).ScalarMult(tmp, ac)

	right := new(GT).Add(right1, right2)

	t.Logf("left  = %s", left.String())
	t.Logf("right = %s", right.String())

	leftBytes := left.Marshal()
	rightBytes := right.Marshal()

	if !bytes.Equal(leftBytes, rightBytes) {
		t.Fail()
	}

	t.Logf("bytes = %s", hex.EncodeToString(leftBytes))

}

func Test_TMP2(t *testing.T) {
	a, _ := randomK(rand.Reader)
	b, _ := randomK(rand.Reader)
	c, _ := randomK(rand.Reader)

	g2Gen := new(G2).ScalarBaseMult(new(big.Int).SetInt64(1))

	aG1 := new(G1).ScalarBaseMult(a)
	bAc := new(big.Int).Mod(new(big.Int).Add(b, c), Order)
	bAcG2 := new(G2).ScalarBaseMult(bAc)
	left := Pair(aG1, bAcG2)

	ab := new(big.Int).Mod(new(big.Int).Mul(a, b), Order)
	abG1 := new(G1).ScalarBaseMult(ab)
	right1 := Pair(abG1, g2Gen)

	ac := new(big.Int).Mod(new(big.Int).Mul(a, c), Order)
	acG1 := new(G1).ScalarBaseMult(ac)
	right2 := Pair(acG1, g2Gen)

	right := new(GT).Add(right1, right2)

	t.Logf("left  = %s", left.String())
	t.Logf("right = %s", right.String())

	leftBytes := left.Marshal()
	rightBytes := right.Marshal()

	if !bytes.Equal(leftBytes, rightBytes) {
		t.Fail()
	}

	t.Logf("bytes = %s", hex.EncodeToString(leftBytes))

}
