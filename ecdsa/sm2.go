/*
Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

	http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.
*/

package ecdsa

// reference to ecdsa
import (
	"crypto/cipher"
	"crypto/elliptic"
	"math/big"
)

func sm2SignGeneric(priv *PrivateKey, csprng *cipher.StreamReader, c elliptic.Curve, hash []byte) (r, s *big.Int, err error) {
	N := c.Params().N
	if N.Sign() == 0 {
		return nil, nil, errZeroParam
	}
	var k *big.Int
	e, err := getE(&priv.PublicKey, nil, hash)
	if err != nil {
		return nil, nil, err
	}
	for { // 调整算法细节以实现SM2
		for {
			k, err = randFieldElement(c, csprng)
			if err != nil {
				r = nil
				return
			}
			r, _ = priv.Curve.ScalarBaseMult(k.Bytes())
			r.Add(r, e)
			r.Mod(r, N)
			if r.Sign() != 0 {
				break
			}
			if t := new(big.Int).Add(r, k); t.Cmp(N) == 0 {
				break
			}
		}
		rD := new(big.Int).Mul(priv.D, r)   //r*D
		s = new(big.Int).Sub(k, rD)         //k-r*D
		d1 := new(big.Int).Add(priv.D, one) //D+1
		//(D+1)^-1
		var d1Inv *big.Int
		if in, ok := priv.Curve.(invertible); ok {
			d1Inv = in.Inverse(d1)
		} else {
			d1Inv = fermatInverse(d1, N) // N != 0
		}
		s.Mul(s, d1Inv) //s*d1Inv
		s.Mod(s, N)     //(s*d1Inv)mod N
		if s.Sign() != 0 {
			break
		}
	}
	return
}

func sm2VerifyGeneric(pub *PublicKey, c elliptic.Curve, hash []byte, r, s *big.Int) bool {
	N := c.Params().N

	if r.Sign() <= 0 || s.Sign() <= 0 {
		return false
	}
	if r.Cmp(N) >= 0 || s.Cmp(N) >= 0 {
		return false
	}

	// 调整算法细节以实现SM2
	t := new(big.Int).Add(r, s)
	t.Mod(t, N)
	if t.Sign() == 0 {
		return false
	}
	var x *big.Int
	if opt, ok := c.(combinedMult); ok {
		x, _ = opt.CombinedMult(pub.X, pub.Y, s.Bytes(), t.Bytes())
	} else {
		x1, y1 := c.ScalarBaseMult(s.Bytes())
		x2, y2 := c.ScalarMult(pub.X, pub.Y, t.Bytes())
		x, _ = c.Add(x1, y1, x2, y2)
	}

	e, err := getE(pub, nil, hash)
	if err != nil {
		return false
	}
	x.Add(x, e)
	x.Mod(x, N)
	return x.Cmp(r) == 0
}
