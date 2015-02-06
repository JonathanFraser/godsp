package godsp

import "math"
import "math/cmplx"

func FFT(data []complex128) ([]complex128, error) {

	//could special case a few more multiples here like 4 which is particularily efficient
	//Also could exploit the concurrent nature of go a little more to get implicit threading
	//could also check for perfect squares because bluestein can be done very efficient there
	if len(data) == 0 {
		return make([]complex128, 0), nil
	}

	if len(data) == 1 {
		return []complex128{data[0]}, nil
	}

	if len(data)%2 == 0 {
		return FFTmod2(data)
	}

	return FFTBluestein(data)
}

func IFFT(data []complex128) ([]complex128, error) {
	swapped := make([]complex128, len(data))
	for i, v := range data {
		swapped[i] = complex(imag(v), real(v))
	}
	ret, err := FFT(swapped)
	if err != nil {
		return nil, err
	}
	for i, v := range ret {
		swapped[i] = complex(imag(v)/float64(len(data)), real(v)/float64(len(data)))
	}
	return swapped, nil
}

//Not very efficient because it does a bunch of copies
//but it does the job
func FFTmod2(data []complex128) ([]complex128, error) {
	ret := make([]complex128, len(data))

	evenPart := make([]complex128, len(data)/2)
	oddPart := make([]complex128, len(data)/2)
	twFactors := make([]complex128, len(data)/2)

	for i := range evenPart {
		evenPart[i] = data[2*i]
		oddPart[i] = data[2*i+1]
		twFactors[i] = cmplx.Exp(complex(0, 2*math.Pi*float64(i)/float64(len(data))))
	}

	evenFFT, err := FFT(evenPart)
	if err != nil {
		return nil, err
	}

	oddFFT, err := FFT(oddPart)
	if err != nil {
		return nil, err
	}

	for i := range evenFFT {
		ret[i] = evenFFT[i] + twFactors[i]*oddFFT[i]
		ret[i+len(data)/2] = evenFFT[i] - twFactors[i]*oddFFT[i]
	}
	return ret, nil
}
func nextPower(value int) int {
	ret := 1
	for ret < value {
		ret = ret * 2
	}
	return ret
}

func FFTBluestein(data []complex128) ([]complex128, error) {
	PadLen := nextPower(2*len(data) - 1)
	chirpedData := make([]complex128, PadLen)
	chirp := make([]complex128, len(data))

	B := make([]complex128, PadLen)

	for i, v := range data {
		chirp[i] = cmplx.Exp(complex(0, -1*math.Pi*float64(i*i)/float64(len(data))))
		B[i] = cmplx.Conj(chirp[i])
		if i != 0 {
			B[PadLen-i] = -1 * cmplx.Conj(chirp[i])
		}
		chirpedData[i] = chirp[i] * v
	}

	retFFT, err := FFT(chirpedData)
	if err != nil {
		return nil, err
	}

	//This can be precomputed and should be... but I'm lazy
	retChirp, err := FFT(B)
	if err != nil {
		return nil, err
	}

	for i, v := range retFFT {
		retFFT[i] = v * retChirp[i]
	}

	retConv, err := IFFT(retFFT)
	if err != nil {
		return nil, err
	}

	for i, v := range retConv[0:len(data)] {
		retConv[i] = v * chirp[i]
	}

	return retConv[0:len(data)], nil
}
