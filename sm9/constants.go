package sm9

import (
	"math/big"
)

func bigFromBase10(s string) *big.Int {
	n, _ := new(big.Int).SetString(s, 10)
	return n
}

//simon add
func bigFromBase16(s string) *big.Int {
	n, _ := new(big.Int).SetString(s, 16)
	return n
}

//sm9 para t:
var t = bigFromBase16("600000000058F98A")

//6*t+2=2400000000215d93c

// p is a prime over which we form a basic field: 36u⁴+36u³+24u²+6u+1.//u = t
var Sm9_p = bigFromBase16("B640000002A3A6F1D603AB4FF58EC74521F2934B1A7AEEDBE56F9B27E351457D") //checked
var p = bigFromBase16("B640000002A3A6F1D603AB4FF58EC74521F2934B1A7AEEDBE56F9B27E351457D")     //checked
var constu = bigFromBase16("16c80000005474de3ac07569feb1d8e8a43e5269634f5ddb7cadf364fc6a28af")
var constuplus1 = bigFromBase16("16c80000005474de3ac07569feb1d8e8a43e5269634f5ddb7cadf364fc6a28b0")
var const2uplus1 = bigFromBase16("2d90000000a8e9bc7580ead3fd63b1d1487ca4d2c69ebbb6f95be6c9f8d4515f")
var mogo_bigInt = bigFromBase16("49bffffffd5c590e29fc54b00a7138bade0d6cb4e58511241a9064d81caeba83")     //R mod p
var mogo_sqr_bigInt = bigFromBase16("2ea795a656f62fbde479b522d6706e7b88f8105fae1a5d3f27dea312b417e2d2") //R*R mod p
// Order is the number of elements in both G₁ and G₂: 36u⁴+36u³+18u²+6u+1.
var Order = bigFromBase16("B640000002A3A6F1D603AB4FF58EC74449F2934B18EA8BEEE56EE19CD69ECF25")

// xiToPMinus1Over6 is ξ^((p-1)/6) where ξ = i; //SM9 //checked;
var xiToPMinus1Over6 = &gfP2{gfP{0x0, 0x0, 0x0, 0x0}, gfP{0x1a98dfbd4575299f, 0x9ec8547b245c54fd, 0xf51f5eac13df846c, 0x9ef74015d5a16393}}

// xiToPMinus1Over3 is ξ^((p-1)/3) where ξ = i.
var xiToPMinus1Over3 = &gfP2{gfP{0x0, 0x0, 0x0, 0x0}, gfP{0xb626197dce4736ca, 0x08296b3557ed0186, 0x9c705db2fd91512a, 0x1c753e748601c992}}

// xiToPMinus1Over2 is ξ^((p-1)/2) where ξ = i.
var xiToPMinus1Over2 = &gfP2{gfP{0x0, 0x0, 0x0, 0x0}, gfP{0x39b4ef0f3ee72529, 0xdb043bf508582782, 0xb8554ab054ac91e3, 0x9848eec25498cab5}}

/**********************************************************/ //0919:add xi^-1
// xiToPMinus1Over3i is (ξ^-1)^((p-1)/3) where ξ = i.(1,0) //checked.
var xiToPMinus1Over3i = &gfP2{gfP{0x0, 0x0, 0x0, 0x0}, gfP{0x646a4b5a4e6783b9, 0xd5e4017f8d980f9d, 0x8d8bf6fd0cdfe790, 0x2d4ac18b775a8f7b}}

// xiToPMinus1Over2i is (ξ^-1)^((p-1)/2) where ξ = i. //checked.
var xiToPMinus1Over2i = &gfP2{gfP{0x0, 0x0, 0x0, 0x0}, gfP{0xabbaac18a46a2054, 0x46ee57561222c759, 0x1dae609fa0e23561, 0x1df7113dae0adc3c}}

// xiToPSquaredMinus1Over3i is (ξ^-1)^((p²-1)/3) where ξ = i. //checked.
var xiToPSquaredMinus1Over3i = &gfP{0x2f4981aa150a0eb3, 0x19c92815c28ded55, 0x39934d9cf7fd761b, 0x99cac18b7ca1dd5f}

// xiToPSquaredMinus1Over3 is ξ^((p²-1)/3) where ξ = i.
var xiToPSquaredMinus1Over3 = &gfP{0x81054fcd94e9c1c4, 0x4c0e91cb8ce2df3e, 0x4877b452e8aedfb4, 0x88f53e748b491776}

// xiTo2PSquaredMinus2Over3 is ξ^((2p²-2)/3) where ξ = i (a cubic root of unity, mod p).
var xiTo2PSquaredMinus2Over3 = &gfP{0x2f4981aa150a0eb3, 0x19c92815c28ded55, 0x39934d9cf7fd761b, 0x99cac18b7ca1dd5f}

// xiToPSquaredMinus1Over6 is ξ^((1p²-1)/6) where ξ = i+3 (a cubic root of -1, mod p).
var xiToPSquaredMinus1Over6 = &gfP{0xb626197dce4736ca, 0x08296b3557ed0186, 0x9c705db2fd91512a, 0x1c753e748601c992}

// xiTo2PMinus2Over3 is ξ^((2p-2)/3) where ξ = i+3.
var xiTo2PMinus2Over3 = &gfP2{gfP{0x0, 0x0, 0x0, 0x0}, gfP{0x81054fcd94e9c1c4, 0x4c0e91cb8ce2df3e, 0x4877b452e8aedfb4, 0x88f53e748b491776}}

// p2 is p, represented as little-endian 64-bit words.
var p2 = [4]uint64{0xe56f9b27e351457d, 0x21f2934b1a7aeedb, 0xd603ab4ff58ec745, 0xb640000002a3a6f1}

// np is the negative inverse of p, mod 2^256.
var np = [4]uint64{0x892bc42c2f2ee42b, 0x181ae39613c8dbaf, 0x966a4b291522b137, 0xafd2bac5558a13b3}

// rN1 is R^-1 where R = 2^256 mod p.
var rN1 = &gfP{0x0a1c7970e5df544d, 0xe74504e9a96b56cc, 0xcda02d92d4d62924, 0x7d2bc576fdf597d1}

// r2 is R^2 where R = 2^256 mod p.
var r2 = &gfP{0x27dea312b417e2d2, 0x88f8105fae1a5d3f, 0xe479b522d6706e7b, 0x2ea795a656f62fbd}

// r3 is R^3 where R = 2^256 mod p.
var r3 = &gfP{0x130257769df5827e, 0x36920fc0837ec76e, 0xcbec24519c22a142, 0x219be84a7c687090}

//need mogo format //somin 1023change and check;
//var G1x = &gfP{0xe8c4e4817c66dddd, 0xe1e4086909dc3280, 0xf5ed0704487d01d6, 0x93de051d62bf718f}
var G1x = &gfP{0x22e935e29860501b, 0xa946fd5e0073282c, 0xefd0cec817a649be, 0x5129787c869140b5}
var G1y = &gfP{0xee779649eb87f7c7, 0x15563cbdec30a576, 0x326353912824efbf, 0x7215717763c39828}

var DoubleGx = &gfP{0x8fdf2548f0fde68, 0xc80ddebf804d6dd4, 0xc8cef5282905b7ca, 0x6007e08434132464}
var kbaseG = bigFromBase16("a5702f05cf1315305e2d6eb64b0deb923db1a0bcf0caff90523ac8754aa6982078559a844411f9825c109f5ee3f52d720dd01785392a727bb1556952b2b013d3")
var mogo = &gfP{0x1a9064d81caeba83, 0xde0d6cb4e5851124, 0x29fc54b00a7138ba, 0x49bffffffd5c590e}
var mogo_sqr = &gfP{0x27dea312b417e2d2, 0x88f8105fae1a5d3f, 0xe479b522d6706e7b, 0x2ea795a656f62fbd}

/*****************************homocrypt para*************************/
//beta:beta^3 =1 mod p
//(1,beta,beta^2)
var beta = &gfP{0x81054fcd94e9c1c4, 0x4c0e91cb8ce2df3e, 0x4877b452e8aedfb4, 0x88f53e748b491776}
var a1 = bigFromBase16("c000000000b1f315")
var b1 = bigFromBase16("-d8000000019062edc000b98b0d64696c")
var a2 = bigFromBase16("d8000000019062ee8000b98b0e165c81")
var b2 = bigFromBase16("c000000000b1f315")

//beta^2
var beta2 = &gfP{0x2f4981aa150a0eb3, 0x19c92815c28ded55, 0x39934d9cf7fd761b, 0x99cac18b7ca1dd5f}
var a1plus = bigFromBase16("d8000000019062edc000b98b0d64696c")
var b1plus = bigFromBase16("-c000000000b1f315")
var a2plus = bigFromBase16("c000000000b1f315")
var b2plus = bigFromBase16("d8000000019062ee8000b98b0e165c81")

/*****************************homo para*************************/
