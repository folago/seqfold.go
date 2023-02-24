package seqfold

import (
	"math"

	"golang.org/x/exp/constraints"
)

func max[T constraints.Ordered](a, b T) T {
	if a > b {
		return a
	}
	return b
}

func min[T constraints.Ordered](a, b T) T {
	if a < b {
		return a
	}
	return b
}

func abs(x int) int {
	if x < 0 {
		return -x
	}
	return x
}

func Reverse(s string) string {
	runes := []rune(s)
	for i, j := 0, len(runes)-1; i < j; i, j = i+1, j-1 {
		runes[i], runes[j] = runes[j], runes[i]
	}
	return string(runes)
}

// RoundFloat runds at the given decimal place
func RoundFloat(number float64, decimalPlace int) float64 {
	// Calculate the 10 to the power of decimal place
	temp := math.Pow(10, float64(decimalPlace))
	// Multiply floating-point number with 10**decimalPlace and round it
	// Divide the rounded number with 10**decimalPlace to get decimal place rounding
	return math.Round(number*temp) / temp
}

func isDNA(seq string) bool {
	for _, base := range seq {
		switch base {
		case 'A', 'C', 'T', 'G':
			continue
		default:
			return false
		}
	}
	return true
}

func isRNA(seq string) bool {
	for _, base := range seq {
		switch base {
		case 'A', 'C', 'U', 'G':
			continue
		default:
			return false
		}
	}
	return true
}

func buildDefaultCaches(seq string) ([][]Struct, [][]Struct) {
	var (
		n       = len(seq)
		v_cache = make([][]Struct, n)
		w_cache = make([][]Struct, n)
		row     = make([]Struct, n)
	)
	for i := 0; i < n; i++ {
		// IJs are nil in the STRUCT_DEFAULT, so copying this is fine
		row[i] = STRUCT_DEFAULT
	}
	for j := 0; j < n; j++ {
		v_cache[j] = make([]Struct, n)
		// IJs are nil in the STRUCT_DEFAULT, so copying this is fine
		copy(v_cache[j], row)

		w_cache[j] = make([]Struct, n)
		copy(w_cache[j], row)
	}
	return v_cache, w_cache
}
