package godsp

import "math"
import "math/cmplx"

func FFT(data []complex128) error {

	//could special case a few more multiples here like 4 which is particularily efficient
	//Also could exploit the concurrent nature of go a little more to get implicit threading
	//could also check for perfect squares because bluestein can be done very efficient there
	if len(data) == 0 {
		return nil
	}

	if len(data) == 1 {
		return nil
	}

	if len(data)%2 == 0 {
		return FFTmod2(data)
	}

	return FFTBluestein(data)
}

func IFFT(data []complex128) error {
	swapped := make([]complex128, len(data))
	for i, v := range data {
		swapped[i] = complex(imag(v), real(v))
	}
	err := FFT(swapped)
	if err != nil {
		return err
	}
	for i, v := range swapped {
		data[i] = complex(imag(v)/float64(len(data)), real(v)/float64(len(data)))
	}
	return nil
}

//Not very efficient because it does a bunch of copies
//but it does the job
func FFTmod2(data []complex128) error {
	evenPart := make([]complex128, len(data)/2)
	oddPart := make([]complex128, len(data)/2)
	twFactors := make([]complex128, len(data)/2)

	for i := range evenPart {
		evenPart[i] = data[2*i]
		oddPart[i] = data[2*i+1]
		twFactors[i] = cmplx.Exp(complex(0, 2*math.Pi*float64(i)/float64(len(data))))
	}

	err := FFT(evenPart)
	if err != nil {
		return err
	}

	err = FFT(oddPart)
	if err != nil {
		return err
	}

	for i := range evenPart {
		data[i] = evenPart[i] + twFactors[i]*oddPart[i]
		data[i+len(data)/2] = evenPart[i] - twFactors[i]*oddPart[i]
	}
	return nil
}
func nextPower(value int) int {
	ret := 1
	for ret < value {
		ret = ret * 2
	}
	return ret
}

func FFTBluestein(data []complex128) error {
	PadLen := nextPower(2*len(data) - 1)
	chirpedData := make([]complex128, PadLen)

	chirp := make([]complex128, len(data))

	wrappedChirp := make([]complex128, PadLen)

	//Prepare chirps
	for i, v := range data {
		chirp[i] = cmplx.Exp(complex(0, -1*math.Pi*float64(i*i)/float64(len(data))))
		wrappedChirp[i] = cmplx.Conj(chirp[i])
		if i != 0 {
			wrappedChirp[PadLen-i] = -1 * cmplx.Conj(chirp[i])
		}
		chirpedData[i] = chirp[i] * v
	}

	err := FFT(chirpedData)
	if err != nil {
		return err
	}

	//This can be precomputed and should be... but I'm lazy
	err = FFT(wrappedChirp)
	if err != nil {
		return err
	}

	for i, v := range wrappedChirp {
		chirpedData[i] = v * wrappedChirp[i]
	}

	err = IFFT(chirpedData)
	if err != nil {
		return err
	}

	for i, v := range chirpedData[0:len(data)] {
		data[i] = v * chirp[i]
	}

	return nil
}
