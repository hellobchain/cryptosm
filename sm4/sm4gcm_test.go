package sm4

import (
	"bytes"
	"testing"
)

func TestSM4GCM(t *testing.T) {
	key := []byte("1234567890abcdef")
	data := []byte{0x01, 0x23, 0x45, 0x67, 0x89, 0xab, 0xcd, 0xef, 0xfe, 0xdc, 0xba, 0x98, 0x76, 0x54, 0x32, 0x10}
	IV := make([]byte, BlockSize)
	testA := [][]byte{ // the length of the A can be random
		{},
		{0x01, 0x23, 0x45, 0x67, 0x89},
		{0x01, 0x23, 0x45, 0x67, 0x89, 0xab, 0xcd, 0xef, 0xfe, 0xdc, 0xba, 0x98, 0x76, 0x54, 0x32, 0x10},
	}
	for _, A := range testA {
		gcmMsg, T, err := SM4GCM(key, IV, data, A, true)
		if err != nil {
			t.Errorf("sm4 enc error:%s", err)
		}
		t.Logf("gcmMsg = %x T=%x\n", gcmMsg, T)
		gcmDec, T_, err := SM4GCM(key, IV, gcmMsg, A, false)
		if err != nil {
			t.Errorf("sm4 dec error:%s", err)
		}
		t.Logf("gcmDec = %x T_=%x\n", gcmDec, T_)
		if bytes.Compare(T, T_) == 0 {
			t.Log("authentication succeed")
		} else {
			t.Error("authentication failed")
		}
		//Failed Test : if we input the different A , that will be a falied result.
		A = []byte{0x01, 0x23, 0x45, 0x67, 0x89, 0xab, 0xcd}
		gcmDec, T_, err = SM4GCM(key, IV, gcmMsg, A, false)
		if err != nil {
			t.Errorf("sm4 dec error:%s", err)
		}
		t.Logf("gcmDec = %x T=%x T_=%x\n", gcmDec, T, T_)
		if bytes.Compare(T, T_) != 0 {
			t.Log("authentication should fail")
		} else {
			t.Error("authentication should not succeed")
		}
	}

}
