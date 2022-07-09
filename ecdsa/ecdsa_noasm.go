// Copyright 2020 The Go Authors. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

//go:build !s390x
// +build !s390x

package ecdsa

import (
	"crypto/cipher"
	"crypto/elliptic"
	"math/big"
)

func sign(priv *PrivateKey, csprng *cipher.StreamReader, c elliptic.Curve, hash []byte) (r, s *big.Int, err error) {
	if c.Params().Name == SM2CurveName {
		return sm2SignGeneric(priv, csprng, c, hash)
	}
	return signGeneric(priv, csprng, c, hash)
}

func verify(pub *PublicKey, c elliptic.Curve, hash []byte, r, s *big.Int) bool {
	if c.Params().Name == SM2CurveName {
		return sm2VerifyGeneric(pub, c, hash, r, s)
	}
	return verifyGeneric(pub, c, hash, r, s)
}
