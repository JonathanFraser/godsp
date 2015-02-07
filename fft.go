package godsp

import "math"
import "math/cmplx"
import "errors"

func isPowerOf2(n int) bool {
	return ((n != 0) && !((n & (n - 1)) != 0))
}

func FFT(data []complex128) error {
	//could special case a few more multiples here like 4 which is particularily efficient
	//Also could exploit the concurrent nature of go a little more to get implicit threading
	//could also check for perfect squares because bluestein can be done very efficient there
	n := len(data)
	if n == 0 {
		return nil
	}

	if n == 1 {
		return nil
	}

	if isPowerOf2(n) {
		return FFTradix2(data)
	}

	if n%4 == 0 {
		return FFTmod4(data)
	}

	if n%2 == 0 {
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

func bitReverse(x, levels int) int {
	result := 0
	for i := 0; i < levels; i++ {
		result = (result << 1) | (x & 1)
		x = x >> 1
	}
	return result
}

func FFTradix2(data []complex128) error {
	// Variables
	n := len(data)

	// Compute levels = floor(log2(n))
	temp := n
	levels := 0
	compare := 1
	for temp > 1 {
		levels++
		temp = temp >> 1
		compare = compare << 1
	}

	if compare != n {
		return errors.New("Not a Power of 2") // n is not a power of 2
	}

	// Trignometric tables
	cmplx_table := make([]complex128, n/2)

	for i := range cmplx_table {
		cmplx_table[i] = cmplx.Conj(cmplx.Exp(complex(0, 2*math.Pi*float64(i)/float64(n))))
	}

	// Bit-reversed addressing permutation
	for i := range data {
		j := bitReverse(i, levels)
		if j > i {
			data[i], data[j] = data[j], data[i]
		}
	}

	// Cooley-Tukey decimation-in-time radix-2 FFT
	for size := 2; size <= n; size *= 2 {
		halfsize := size / 2
		tablestep := n / size

		for i := 0; i < n; i += size {
			k := 0
			for j := i; j < i+halfsize; j++ {
				t := data[j+halfsize] * cmplx_table[k]
				data[j+halfsize] = data[j] + t
				data[j] += t
				k += tablestep
			}
		}
	}

	return nil
}

func FFTmod4(data []complex128) error {
	//TODO: implement a recursing multiple of 4 implentation
	//useful because twiddle factors are all 1's and j's
	return FFTmod2(data)
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
