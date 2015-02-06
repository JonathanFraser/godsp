package godsp

import "testing"
import "math"

type FFTData struct {
	input     []complex128
	output    []complex128
	tolerance float64
}

var FFTtest = []FFTData{
	{[]complex128{}, []complex128{}, 1e-12},
	{[]complex128{1}, []complex128{1}, 1e-12},
	{[]complex128{1, 0}, []complex128{1, 1}, 1e-12},
	{[]complex128{1, 0, 0}, []complex128{1, 1, 1}, 1e-12},
	{[]complex128{1, 0, 0, 0}, []complex128{1, 1, 1, 1}, 1e-12},
}

var IFFTtest = []FFTData{
	{[]complex128{}, []complex128{}, 1e-12},
	{[]complex128{1}, []complex128{1}, 1e-12},
	{[]complex128{1, 1}, []complex128{1, 0}, 1e-12},
	{[]complex128{1, 1, 1}, []complex128{1, 0, 0}, 1e-12},
	{[]complex128{1, 1, 1, 1}, []complex128{1, 0, 0, 0}, 1e-12},
}

func verify(expected, received []complex128, tol float64) bool {
	if len(expected) != len(received) {
		return false
	}
	good := true
	for i, v := range expected {
		diff := v - received[i]
		good = good || math.Abs(real(diff)) < tol
		good = good || math.Abs(imag(diff)) < tol
	}
	return good
}

func TestFFT(t *testing.T) {
	for i, v := range FFTtest {
		out, err := FFT(v.input)
		if err != nil {
			t.Error("FFT erroneously resulted in: " + err.Error())
			return
		}
		if !verify(v.output, out, v.tolerance) {
			t.Error("FFT failed tolerance check item:", i, "\n")
			return
		}
	}
	t.Log("FFT tests Passed\n")
}

func TestIFFT(t *testing.T) {
	for i, v := range IFFTtest {
		out, err := IFFT(v.input)
		if err != nil {
			t.Error("IFFT erroneously resulted in: " + err.Error())
			return
		}
		if !verify(v.output, out, v.tolerance) {
			t.Error("IFFT failed tolerance check item:", i, "\n")
			return
		}
	}
	t.Log("FFT tests Passed\n")
}

func BenchmarkFFT(t *testing.B) {
	t.StopTimer()
	data := make([]complex128, 997) //this is a prime
	data[0] = complex(1, 0)
	t.StartTimer()
	for i := 0; i < t.N; i++ {
		FFT(data)
	}
}
