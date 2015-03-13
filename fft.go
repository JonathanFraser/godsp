package godsp

import "errors"
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

	if n > 1e6 { //tune this number against the goroutine overhead
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

	// Bit-reversed addressing permutation
	for i := range data {
		j := bitReverse(i, levels)
		if j > i {
			data[i], data[j] = data[j], data[i]
		}
	}

	factors := make([]complex128, n)
	for i := range factors {
		arg := -2.0 * math.Pi * float64(i) / float64(n)
		val := complex(math.Cos(arg), math.Sin(arg))
		factors[i] = val
	}

	// Cooley-Tukey decimation-in-time radix-2 FFT
	for size := 2; size <= n; size *= 2 {
		halfsize := size / 2
		tablestep := n / size
		for i := 0; i < n; i += size {
			k := 0
			for j := i; j < i+halfsize; j++ {
				t := data[j+halfsize] * factors[k]
				data[j+halfsize] = data[j] + t
				data[j] += t
				k += tablestep
			}
		}
	}

	return nil
}

func getMLevel(input int) (int, error) {
	out := 0
	for input != 1 {
		if input%4 != 0 {
			return 0, errors.New("input not a power of 4")
		}
		input = input >> 2
		out++
	}
	return out, nil
}

func FFTRadix4(data []complex128) error {
	//public static void FFTR4(double[] X, double[] Y, int N, int M) {
	M, err := getMLevel(len(data))
	if err != nil {
		return err
	}

	// N = 4 ^ M
	// N = 1 << (M+M);
	N := len(data)
	N1 := 0
	N2 := len(data)
	I1 := 0
	I2 := 0
	I3 := 0

	twiddles := make([]complex128, N)
	for i := range twiddles {
		arg := 2 * math.Pi * float64(i) / float64(N)
		twiddles[i] = complex(math.Cos(arg), math.Sin(arg))
	}

	for K := 0; K < M; K++ {
		N1 = N2
		N2 = N2 / 4

		for J := 0; J < N2; J++ {
			//Should be pre-calculated for optimization
			ind := J << uint(K)
			Tw1 := twiddles[ind%N]
			Tw2 := twiddles[(2*ind)%N]
			Tw3 := twiddles[(3*ind)%N]

			CO1 := real(Tw1)
			CO2 := real(Tw2)
			CO3 := real(Tw3)
			SI1 := imag(Tw1)
			SI2 := imag(Tw2)
			SI3 := imag(Tw3)

			for I := J; I < N; I += N1 {
				I1 = I + N2
				I2 = I1 + N2
				I3 = I2 + N2
				C1 := data[I] + data[I2]
				C3 := data[I] - data[I2]
				C2 := data[I1] + data[I3]
				C4 := data[I1] - data[I3]
				data[I] = C1 + C2

				//These adds and subtracts are probably efficient enough
				C2 = C1 - C2
				C1 = complex(real(C3)-imag(C4), imag(C3)+real(C4))
				C3 = complex(real(C3)+imag(C4), imag(C3)-real(C4))

				//These mults should be able to be converted into the complex domain
				data[I1] = complex(CO1*real(C3)+SI1*imag(C3), CO1*imag(C3)-SI1*real(C3))
				data[I2] = complex(CO2*real(C2)+SI2*imag(C2), CO2*imag(C2)-SI2*real(C2))
				data[I3] = complex(CO3*real(C1)+SI3*imag(C1), CO3*imag(C1)-SI3*real(C1))
			}
		}
	}

	// Radix-4 bit-reverse
	J := 0
	N2 = N >> 2
	for I := 0; I < N-1; I++ {
		if I < J {
			data[I], data[J] = data[J], data[I]
		}
		N1 = N2
		for J >= 3*N1 {
			J -= 3 * N1
			N1 >>= 2
		}
		J += N1
	}
	return nil
}

func FFTmod4(data []complex128) error {
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
