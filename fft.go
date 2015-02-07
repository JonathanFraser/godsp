package godsp

import "math"
import "math/cmplx"

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

	if n > 1000 { //tune this number against the goroutine overhead
		if n%4 == 0 {
			return FFTmod4(data)
		}

		if n%2 == 0 {
			return FFTmod2(data)
		}
	}

	if isPowerOf2(n) {
		return FFTradix2(data)
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

func getLevels(size int) int {
	// Compute levels = floor(log2(n))
	levels := 0
	for size > 1 {
		levels++
		size = size >> 1
	}

	return levels
}

func FFTradix2(data []complex128) error {
	// Variables
	n := len(data)

	levels := getLevels(n)

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

func FFTmod2(data []complex128) error {
	FirstPart := data[:len(data)/2]
	SecondPart := data[len(data)/2:]

	for i := range FirstPart {
		tw := cmplx.Exp(complex(0, 2*math.Pi*float64(i)/float64(len(data))))
		t := tw * (FirstPart[i] - SecondPart[i])
		FirstPart[i] = FirstPart[i] + SecondPart[i]
		SecondPart[i] = t
	}

	//Parallelize the decomp here, mod2 should only be called for very large data
	evenResult := make(chan error)
	oddResult := make(chan error)
	go func() { evenResult <- FFT(FirstPart) }()
	go func() { oddResult <- FFT(SecondPart) }()

	if err := <-evenResult; err != nil {
		return err
	}

	if err := <-oddResult; err != nil {
		return err
	}
	levels := getLevels(len(data))

	// Bit-reversed addressing permutation
	for i := range data {
		j := bitReverse(i, levels)
		if j > i {
			data[i], data[j] = data[j], data[i]
		}
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
	//This would save alot fo memory and increase speed by a factor
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
