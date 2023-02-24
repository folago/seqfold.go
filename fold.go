package seqfold

import (
	"fmt"
	"math"
	"strings"
)

// Fold the DNA sequence and return the lowest free energy score.
//
// Based on the approach described in:
// Zuker and Stiegler, 1981
// https://www.ncbi.nlm.nih.gov/pmc/articles/PMC326673/pdf/nar00394-0137.pdf
//
// If the sequence is 50 or more bp long, "isolated" matching bp
// are ignored in V(i,j). This is based on an approach described in:
// Mathews, Sabina, Zuker and Turner, 1999
// https://www.ncbi.nlm.nih.gov/pubmed/10329189
//
// Args:
//     seq: The sequence to Fold
// Keyword Args:
//     temp: The temperature the Fold takes place in, in Celcius
// Returns:
//     List[Struct]: A list of structures. Stacks, bulges, hairpins, etc.
func Fold(seq string, temp float64) []Struct {
	v_cache, w_cache := cache(seq, temp)
	n := len(seq)

	// get the minimum free energy structure out of the cache
	return traceback(0, n-1, v_cache, w_cache)
}

// Fold the sequence and return just the delta G of the structure
//
// Args:
//     seq: The sequence to fold
// Keyword Args:
//     temp: The temperature to fold at
// Returns:
//     float: The minimum free energy of the folded sequence
func dg(seq string, temp float64) float64 {
	structs := Fold(seq, temp)
	sumE := 0.0
	for _, s := range structs {
		sumE += s.E
	}
	return RoundFloat(sumE, 2)
}

// Fold a nucleic acid sequence and return the estimated dg of each (i,j) pairing.
//
// Args:
//     seq: The nucleic acid sequence to fold
// Keyword Args:
//     temp: The temperature to fold at
// Returns:
//     Cache: A 2D matrix where each (i, j) pairing corresponds to the
//         minimum free energy between i and j
func dg_cache(seq string, temp float64) [][]float64 {
	_, w_cache := cache(seq, temp)

	cache := make([][]float64, len(seq))
	for i := range cache {
		cache[i] = make([]float64, len(seq))
	}

	for i, w_row := range w_cache {
		for j, w_struct := range w_row {
			cache[i][j] = w_struct.E
		}
	}

	return cache
}

// Get the dot bracket notation for a secondary structure.
//
// Args:
//     structs: A list of structs, usually from the fold function
// Returns:
//     str: The dot bracket notation of the secondary structure
func DotBracket(structs []Struct) string {
	// log structure with dot-bracket notation
	// result = ["."] * (max(j for s in structs for _, j in s.ij) + 1)
	maxj := 0
	for _, s := range structs {
		for _, ij := range s.IJ {
			if ij.J > maxj {
				maxj = ij.J
			}
		}
	}
	maxj += 1
	result := make([]byte, maxj)
	for i := range result {
		result[i] = '.'
	}
	for _, s := range structs {
		if len(s.IJ) == 1 {
			ij := s.IJ[0]
			result[ij.I] = '('
			result[ij.J] = ')'
		}
	}
	return string(result)
}

// Create caches for the w_cache and v_cache
//
// The Structs is useful for gathering many possible energies
// between a series of (i,j) combinations.
// Args:
//     seq: The sequence to fold
// Keyword Args:
//     temp: The temperature to fold at
// Returns:
//     (Structs, Structs): The w_cache and the v_cache for traversal later
func cache(seq string, temp float64) ([][]Struct, [][]Struct) {
	// if it's a SeqRecord, gather just the seq
	//    if "SeqRecord" in str(type(seq)):
	//       seq = str(seq.seq)  # type: ignore

	seq = strings.ToUpper(seq)
	temp = temp + 273.15 // kelvin

	// figure out whether it's DNA or RNA, choose energy map
	var emap Energies
	switch {
	case isDNA(seq):
		emap = DNA_ENERGIES
	case isRNA(seq):
		emap = RNA_ENERGIES
	default:
		// TODO handle errors
		panic("not DNA nor RNA")
	}

	var (
		n       = len(seq)
		v_cache = make([][]Struct, n)
		w_cache = make([][]Struct, n)
		row     = make([]Struct, n)
	)
	for i := 0; i < n; i++ {
		row[i] = STRUCT_DEFAULT
	}
	for j := 0; j < n; j++ {
		v_cache[j] = make([]Struct, n)
		copy(v_cache[j], row)

		w_cache[j] = make([]Struct, n)
		copy(w_cache[j], row)
	}

	// fill the cache
	w(seq, 0, n-1, temp, v_cache, w_cache, emap)

	return v_cache, w_cache
}

// Calculate the free energy associated with a bulge.
//    seq: The full folding DNA sequence
//    i: The start index of the bulge
//    i1: The index to the right of i
//    j: The end index of the bulge
//    j1: The index to the left of j
//    loop: The sequence of the bulge
//    temp: Temperature in Kelvin
//    emap: Map to DNA/RNA energies
// Returns:
//    float: The increment in free energy from the bulge
func bulge(seq string, i, i1, j, j1 int, temp float64, emap Energies) float64 {
	loop_len := max(i1-i-1, j-j1-1)
	if loop_len <= 0 {
		// raise RuntimeError
		// TODO handle errors
		panic("loop_len <= 0")
	}

	var dG float64

	// add penalty based on size
	if en, ok := emap.BULGE_LOOPS[loop_len]; ok {
		d_h, d_s := en.H, en.S
		dG = d_g(d_h, d_s, temp)
	} else {
		// it's too large for pre-calculated list, extrapolate
		en := emap.BULGE_LOOPS[30]
		d_h, d_s := en.H, en.S
		dG = d_g(d_h, d_s, temp)
		dG = j_s(loop_len, 30, dG, temp)
	}

	if loop_len == 1 {
		// if len 1, include the delta G of intervening NN (SantaLucia 2004)
		pair := pair(seq, i, i1, j, j1)
		// TODO handle errors
		// assert pair in emap.NN
		if _, ok := emap.NN[pair]; !ok {
			panic("assert pair not in emap.NN")
		}
		dG += stack(seq, i, i1, j, j1, temp, emap)
	}

	// penalize AT terminal bonds
	for _, k := range []int{i, i1, j, j1} {
		if seq[k] == 'A' {
			dG += 0.5
		}
	}

	return dG
}

// Find and return the lowest free energy structure in Sij subsequence
//
// Figure 2B in Zuker and Stiegler, 1981
// Args:
//     seq: The sequence being folded
//     i: The start index
//     j: The end index (inclusive)
//     temp: The temperature in Kelvin
//     v_cache: Free energy cache for if i and j bp
//     w_cache: Free energy cache for lowest energy structure from i to j. 0 otherwise
// Returns:
//     float: The free energy for the subsequence from i to j
func w(seq string, i, j int, temp float64, v_cache, w_cache [][]Struct, emap Energies) Struct {
	if !w_cache[i][j].Equal(STRUCT_DEFAULT) {
		return w_cache[i][j]
	}

	if j-i < 4 {
		w_cache[i][j] = STRUCT_NULL
		return w_cache[i][j]
	}

	w1 := w(seq, i+1, j, temp, v_cache, w_cache, emap)
	w2 := w(seq, i, j-1, temp, v_cache, w_cache, emap)
	w3 := v(seq, i, j, temp, v_cache, w_cache, emap)

	w4 := STRUCT_NULL
	for k := i + 1; k < j-1; k++ {
		w4_test := multi_branch(seq, i, k, j, temp, v_cache, w_cache, emap, false)

		if w4_test.Valid() && w4_test.E < w4.E {
			w4 = w4_test
		}
	}

	wret := min_struct([]Struct{w1, w2, w3, w4})
	w_cache[i][j] = wret
	return wret
}

// Find, store and return the minimum free energy of the structure between i and j
//
// If i and j don't bp, store and return INF.
// See: Figure 2B of Zuker, 1981
// Args:
//     seq: The sequence being folded
//     i: The start index
//     j: The end index (inclusive)
//     temp: The temperature in Kelvin
//     v_cache: Free energy cache for if i and j bp. INF otherwise
//     w_cache: Free energy cache for lowest energy structure from i to j. 0 otherwise
//     emap: Energy map for DNA/RNA
// Returns:
//    float: The minimum energy folding structure possible between i and j on seq
func v(seq string, i, j int, temp float64, v_cache, w_cache [][]Struct, emap Energies) Struct {
	if !v_cache[i][j].Equal(STRUCT_DEFAULT) {
		return v_cache[i][j]
	}

	// the ends must basepair for V(i,j)
	if emap.COMPLEMENT[seq[i]] != seq[j] {
		v_cache[i][j] = STRUCT_NULL
		return v_cache[i][j]
	}
	// if the basepair is isolated, and the seq large, penalize at 1,600 kcal/mol
	// heuristic for speeding this up
	// from https://www.ncbi.nlm.nih.gov/pubmed/10329189
	isolated_outer := true
	if i > 0 && j < len(seq)-1 {
		isolated_outer = emap.COMPLEMENT[seq[i-1]] != seq[j+1]
	}
	isolated_inner := emap.COMPLEMENT[seq[i+1]] != seq[j-1]

	if isolated_outer && isolated_inner {
		v_cache[i][j] = Struct{E: 1600}
		return v_cache[i][j]
	}

	// E1 = FH(i, j); hairpin
	p := pair(seq, i, i+1, j, j-1)
	e1 := Struct{E: hairpin(seq, i, j, temp, emap), Desc: "HAIRPIN:" + p}
	if j-i == 4 { // small hairpin; 4bp
		v_cache[i][j] = e1
		w_cache[i][j] = e1
		return v_cache[i][j]
	}

	// E2 = min{FL(i, j, i', j') + V(i', j')}, i<i'<j'<j
	// stacking region or bulge or interior loop; Figure 2A(2)
	// j-i=d>4; various pairs i',j' for j'-i'<d
	n := len(seq)
	e2 := Struct{E: math.Inf(1)}
	for i1 := i + 1; i1 < j-4; i1++ {
		for j1 := i1 + 4; j1 < j; j1++ {
			// i1 and j1 must match
			if emap.COMPLEMENT[seq[i1]] != seq[j1] {
				continue
			}

			p := pair(seq, i, i1, j, j1)
			pair_left := pair(seq, i, i+1, j, j-1)
			pair_right := pair(seq, i1-1, i1, j1+1, j1)
			_, plin := emap.NN[pair_left]
			_, prin := emap.NN[pair_right]
			pair_inner := plin || prin

			stck := i1 == i+1 && j1 == j-1
			bulge_left := i1 > i+1
			bulge_right := j1 < j-1

			var (
				e2_test      = math.Inf(1)
				e2_test_type = ""
			)
			switch {
			case stck:
				// it's a neighboring/stacking pair in a helix
				e2_test = stack(seq, i, i1, j, j1, temp, emap)
				e2_test_type = fmt.Sprintf("STACK:%s", p)

				if i > 0 && j == n-1 || i == 0 && j < n-1 {
					// there's a dangling end
					e2_test_type = fmt.Sprintf("STACK_DE:%s", p)
				}
			case bulge_left && bulge_right && !pair_inner:
				// it's an interior loop
				e2_test = internal_loop(seq, i, i1, j, j1, temp, emap)
				e2_test_type = fmt.Sprintf("INTERIOR_LOOP:%d/%d", i1-i, j-j1)

				if i1-i == 2 && j-j1 == 2 {
					loop_left := seq[i : i1+1]
					loop_right := seq[j1 : j+1]
					// technically an interior loop of 1. really 1bp mismatch
					e2_test_type = fmt.Sprintf("STACK:%s/%s", loop_left, Reverse(loop_right))
				}
			case bulge_left && !bulge_right:
				// it's a bulge on the left side
				e2_test = bulge(seq, i, i1, j, j1, temp, emap)
				e2_test_type = fmt.Sprintf("BULGE:%d", i1-i)
			case !bulge_left && bulge_right:
				// it's a bulge on the right side
				e2_test = bulge(seq, i, i1, j, j1, temp, emap)
				e2_test_type = fmt.Sprintf("BULGE:%d", j-j1)
			default:
				// it's basically a hairpin, only outside bp match
				continue
			}

			// add V(i', j')
			e2_test += v(seq, i1, j1, temp, v_cache, w_cache, emap).E
			if e2_test != math.Inf(-1) && e2_test < e2.E {
				e2 = Struct{E: e2_test, Desc: e2_test_type, IJ: []IJ{{i1, j1}}}
			}
		}
	}

	// E3 = min{W(i+1,i') + W(i'+1,j-1)}, i+1<i'<j-2
	e3 := STRUCT_NULL
	if !isolated_outer || i == 0 || j == len(seq)-1 {
		for k := i + 1; k < j-1; k++ {
			e3_test := multi_branch(seq, i, k, j, temp, v_cache, w_cache, emap, true)

			if e3_test.Valid() && e3_test.E < e3.E {
				e3 = e3_test
			}
		}
	}

	e := min_struct([]Struct{e1, e2, e3})
	v_cache[i][j] = e
	return e
}

func add_branch(s Struct, branches *[]IJ, seq string, i, j int, temp float64, v_cache, w_cache [][]Struct, emap Energies) {
	if !s.Valid() || len(s.IJ) == 0 {
		return
	}
	if len(s.IJ) == 1 {
		*branches = append(*branches, s.IJ[0])
		return
	}
	for _, ij := range s.IJ {
		str := w(seq, ij.I, ij.J, temp, v_cache, w_cache, emap)
		add_branch(str, branches, seq, i, j, temp, v_cache, w_cache, emap)
	}
}

// Calculate a multi-branch energy penalty using a linear formula.
//
// From Jaeger, Turner, and Zuker, 1989.
// Found to be better than logarithmic in Ward, et al. 2017
// Args:
//     seq: The sequence being folded
//     i: The left starting index
//     k: The mid-point in the search
//     j: The right ending index
//     temp: Folding temp
//     v_cache: Structs of energies where V(i,j) bond
//     w_cache: Structs of min energy of substructures between W(i,j)
//     helix: Whether this multibranch is enclosed by a helix
//     emap: Map to DNA/RNA energies
// Keyword Args:
//     helix: Whether V(i, j) bond with one another in a helix
// Returns:
//     Struct: A multi-branch structure
func multi_branch(seq string, i, k, j int, temp float64, v_cache, w_cache [][]Struct, emap Energies, helix bool) Struct {
	var left, right Struct
	if helix {
		left = w(seq, i+1, k, temp, v_cache, w_cache, emap)
		right = w(seq, k+1, j-1, temp, v_cache, w_cache, emap)
	} else {
		left = w(seq, i, k, temp, v_cache, w_cache, emap)
		right = w(seq, k+1, j, temp, v_cache, w_cache, emap)
	}

	if !left.Valid() || !right.Valid() {
		return STRUCT_NULL
	}

	// gather all branches of this multi-branch structure
	var branches []IJ

	// in python this was a recursive closure, in Go this is not possible so
	// we pull it out and pass all the parameters
	add_branch(left, &branches, seq, i, j, temp, v_cache, w_cache, emap)
	add_branch(right, &branches, seq, i, j, temp, v_cache, w_cache, emap)

	// this isn't multi-branched
	if len(branches) < 2 {
		return STRUCT_NULL
	}

	// if there's a helix, i,j counts as well
	if helix {
		branches = append(branches, IJ{i, j})
	}

	// count up unpaired bp and asymmetry
	branches_count := len(branches)
	unpaired := 0
	e_sum := 0.0
	ij := IJ{i, j}
	for index, ij2 := range branches {
		i2, j2 := ij2.I, ij2.J
		ij1 := branches[abs((index-1)%len(branches))]
		j1 := ij1.J
		ij3 := branches[abs((index+1)%len(branches))]
		i3, j3 := ij3.I, ij3.J

		// add energy from unpaired bp to the right
		// of the helix as though it was a dangling end
		// if there's only one bp, it goes to whichever
		// helix (this or the next) has the more favorable energy
		unpaired_left := 0
		unpaired_right := 0
		de := 0.0
		if index == len(branches)-1 && !helix {
			// TODO what here?
			// pass
		} else if ij3 == ij {
			unpaired_left = i2 - j1 - 1
			unpaired_right = j3 - j2 - 1

			if unpaired_left != 0 && unpaired_right != 0 {
				de = stack(seq, i2-1, i2, j2+1, j2, temp, emap)
			} else if unpaired_right != 0 {
				de = stack(seq, -1, i2, j2+1, j2, temp, emap)
				if unpaired_right == 1 {
					de = min(stack(seq, i3, -1, j3, j3-1, temp, emap), de)
				}
			}
		} else if ij2 == ij {
			unpaired_left = j2 - j1 - 1
			unpaired_right = i3 - i2 - 1

			if unpaired_left != 0 && unpaired_right != 0 {
				de = stack(seq, i2-1, i2, j2+1, j2, temp, emap)
			} else if unpaired_right != 0 {
				de = stack(seq, i2, i2+1, j2, -1, temp, emap)
				if unpaired_right == 1 {
					de = min(stack(seq, i3-1, i3, -1, j3, temp, emap), de)
				}
			}
		} else {
			unpaired_left = i2 - j1 - 1
			unpaired_right = i3 - j2 - 1

			if unpaired_left != 0 && unpaired_right != 0 {
				de = stack(seq, i2-1, i2, j2+1, j2, temp, emap)
			} else if unpaired_right != 0 {
				de = stack(seq, -1, i2, j2+1, j2, temp, emap)
				if unpaired_right == 1 {
					de = min(stack(seq, i2-1, i2, j2+1, j2, temp, emap), de)
				}
			}
		}
		e_sum += de
		unpaired += unpaired_right
		// assert unpaired_right >= 0
		// TODO handle assertion
		if unpaired_right < 0 {
			panic("unpaired_right >= 0")
		}

		if ij2 != ij { // add energy
			e_sum += w(seq, i2, j2, temp, v_cache, w_cache, emap).E
		}

	}

	// assert unpaired >= 0
	// TODO handle assertion
	if unpaired < 0 {
		panic("assert unpaired >= 0")
	}

	// penalty for unmatched bp and multi-branch
	a, b, c, d := emap.MULTIBRANCH.A, emap.MULTIBRANCH.B, emap.MULTIBRANCH.C, emap.MULTIBRANCH.D
	e_multibranch := a + b*float64(len(branches)) + c*float64(unpaired)

	if unpaired == 0 {
		e_multibranch = a + d
	}

	// energy of min-energy neighbors
	e := e_multibranch + e_sum

	// pointer to next structures
	if helix {
		// branches.pop()
		branches = branches[:len(branches)-1]
	}

	return Struct{E: e, Desc: fmt.Sprintf("BIFURCATION:%dn/%dh", unpaired, branches_count), IJ: branches}
}

// Calculate the free energy of an internal loop.
//
// The first and last bp of both left and right sequences
// are not themselves parts of the loop, but are the terminal
// bp on either side of it. They are needed for when there's
// a single internal looping bp (where just the mismatching
// free energies are used)
// Note that both left and right sequences are in 5' to 3' direction
// This is adapted from the "Internal Loops" section of SantaLucia/Hicks, 2004
// Args:
//     seq: The sequence we're folding
//     i: The index of the start of structure on left side
//     i1: The index to the right of i
//     j: The index of the end of structure on right side
//     j1: The index to the left of j
//     temp: Temperature in Kelvin
//     emap: Dictionary mapping to energies for DNA/RNA
// Returns:
//    float: The free energy associated with the internal loop
func internal_loop(seq string, i, i1, j, j1 int, temp float64, emap Energies) float64 {
	loop_left := i1 - i - 1
	loop_right := j - j1 - 1
	loop_len := loop_left + loop_right

	if loop_left < 1 || loop_right < 1 {
		// raise RuntimeError
		// TODO error handling
		panic("innternal_oop")
	}

	// single bp mismatch, sum up the two single mismatch pairs
	if loop_left == 1 && loop_right == 1 {
		mm_left := stack(seq, i, i1, j, j1, temp, emap)
		mm_right := stack(seq, i1-1, i1, j1+1, j1, temp, emap)
		return mm_left + mm_right
	}
	var d_h, d_s, dG float64
	// apply a penalty based on loop size
	if en, ok := emap.INTERNAL_LOOPS[loop_len]; ok {
		d_h, d_s = en.H, en.S
		dG = d_g(d_h, d_s, temp)
	} else {
		// it's too large an internal loop, extrapolate
		en := emap.INTERNAL_LOOPS[30]
		d_h, d_s = en.H, en.S
		dG = d_g(d_h, d_s, temp)
		dG = j_s(loop_len, 30, dG, temp)
	}

	// apply an asymmetry penalty
	loop_asymmetry := math.Abs(float64(loop_left - loop_right))
	dG += 0.3 * loop_asymmetry

	// apply penalty based on the mismatching pairs on either side of the loop
	pair_left_mm := pair(seq, i, i+1, j, j-1)
	en := emap.TERMINAL_MM[pair_left_mm]
	d_h, d_s = en.H, en.S
	dG += d_g(d_h, d_s, temp)

	pair_right_mm := pair(seq, i1-1, i1, j1+1, j1)
	en = emap.TERMINAL_MM[pair_right_mm]
	d_h, d_s = en.H, en.S
	dG += d_g(d_h, d_s, temp)

	return dG
}

// Return the struct with the lowest free energy that isn't -inf (undef)
// Args:
//    structs: Structures being compared
// Returns:
//    struct: The min free energy structure
func min_struct(structs []Struct) Struct {
	s := STRUCT_NULL
	for _, str := range structs {
		if str.E != math.Inf(-1) && str.E < s.E {
			s = str
		}
	}
	return s
}

// Get the free energy for a stack.
//
// Using the indexes i and j, check whether it's at the end of
// the sequence or internal. Then check whether it's a match
// or mismatch, and return.
// Two edge-cases are terminal mismatches and dangling ends.
// The energy of a dangling end is added to the energy of a pair
// where i XOR j is at the sequence's end.
// Args:
//     seq: The full folding sequence
//     i: The start index on left side of the pair/stack
//     i1: The index to the right of i
//     j: The end index on right side of the pair/stack
//     j1: The index to the left of j
//     temp: Temperature in Kelvin
// Returns:
//     float: The free energy of the NN pairing
func stack(seq string, i, i1, j, j1 int, temp float64, emap Energies) float64 {
	// if any(x >= len(seq) for x in [i, i1, j, j1]):
	//    return 0.0
	for _, x := range []int{i, i1, j, j1} {
		if x >= len(seq) {
			return 0
		}
	}

	pair := pair(seq, i, i1, j, j1)
	// if any(x == -1 for x in [i, i1, j, j1]):
	for _, x := range []int{i, i1, j, j1} {
		if x == -1 {
			// it's a dangling end
			en := emap.DE[pair]
			d_h, d_s := en.H, en.S
			return d_g(d_h, d_s, temp)
		}
	}

	if i > 0 && j < len(seq)-1 {
		// it's internal
		en, ok := emap.NN[pair]
		if !ok {
			en = emap.INTERNAL_MM[pair]
		}
		d_h, d_s := en.H, en.S
		return d_g(d_h, d_s, temp)
	}
	if i == 0 && j == len(seq)-1 {
		// it's terminal
		en, ok := emap.NN[pair]
		if !ok {
			en = emap.INTERNAL_MM[pair]
		}
		d_h, d_s := en.H, en.S
		return d_g(d_h, d_s, temp)
	}

	if i > 0 && j == len(seq)-1 {
		// it's dangling on left
		en, ok := emap.NN[pair]
		if !ok {
			en = emap.INTERNAL_MM[pair]
		}
		d_h, d_s := en.H, en.S
		dG := d_g(d_h, d_s, temp)

		pair_de := fmt.Sprintf("%c%c/.%c", seq[i-1], seq[i], seq[j])
		if en, ok := emap.DE[pair_de]; ok {
			d_h, d_s := en.H, en.S
			dG += d_g(d_h, d_s, temp)
		}
		return dG
	}

	if i == 0 && j < len(seq)-1 {
		// it's dangling on right
		en, ok := emap.NN[pair]
		if !ok {
			en = emap.INTERNAL_MM[pair]
		}
		d_h, d_s := en.H, en.S
		dG := d_g(d_h, d_s, temp)

		pair_de := fmt.Sprintf(".%c/%c%c", +seq[i], seq[j+1], seq[j])
		if en, ok := emap.DE[pair_de]; ok {
			d_h, d_s := en.H, en.S
			dG += d_g(d_h, d_s, temp)
		}
	}
	return 0
}

// Find the free energy given delta h, s and temp
// Args:
//    d_h: The enthalpy increment in kcal / mol
//    d_s: The entropy increment in cal / mol
//    temp: The temperature in Kelvin
// Returns:
//    The free energy increment in kcal / (mol x K)
func d_g(d_h, d_s, temp float64) float64 {
	return d_h - temp*(d_s/1000.0)
}

// Estimate the free energy of length query_len based on one of length known_len.
//
// The Jacobson-Stockmayer entry extrapolation formula is used
// for bulges, hairpins, etc that fall outside the 30nt upper limit
// for pre-calculated free-energies. See SantaLucia and Hicks (2004).
// Args:
//     query_len: Length of element without known free energy value
//     known_len: Length of element with known free energy value (d_g_x)
//     d_g_x: The free energy of the element known_len
//     temp: Temperature in Kelvin
// Returns:
//     float: The free energy for a structure of length query_len
func j_s(query_len, known_len int, d_g_x, temp float64) float64 {
	gas_constant := 1.9872e-3
	return d_g_x + 2.44*gas_constant*temp*math.Log(float64(query_len)/float64(known_len))
}

// pair Returns a stack representation, a key for the NN maps
// Args:
//    s: Sequence being folded
//    i: leftmost index
//    i1: index to right of i
//    j: rightmost index
//    j1: index to left of j
// Returns:
//    str: string representation of the pair
func pair(s string, i, i1, j, j1 int) string {
	ss := []rune(s)

	//return (
	//    (s[i] if i >= 0 else ".")
	//    + (s[i1] if i1 >= 0 else ".")
	//    + "/"
	//    + (s[j] if j >= 0 else ".")
	//    + (s[j1] if j1 >= 0 else ".")
	//)
	ret := []rune{'.', '.', '/', '.', '.'}
	if i >= 0 {
		ret[0] = ss[i]
	}
	if i1 >= 0 {
		ret[1] = ss[i1]
	}
	if j >= 0 {
		ret[3] = ss[j]
	}
	if j1 >= 0 {
		ret[4] = ss[j1]
	}
	return string(ret)
}

// Calculate the free energy of a hairpin.
// Args:
//    seq: The sequence we're folding
//    i: The index of start of hairpin
//    j: The index of end of hairpin
//    temp: Temperature in Kelvin
//    emap: Map of energies
// Returns:
//    float: The free energy increment from the hairpin structure
func hairpin(seq string, i, j int, temp float64, emap Energies) float64 {
	if j-i < 4 {
		return math.Inf(1)
	}

	hairpinSeq := seq[i : j+1]
	hairpin_len := len(hairpinSeq) - 2
	pair := pair(seq, i, i+1, j, j-1)

	if emap.COMPLEMENT[hairpinSeq[0]] != hairpinSeq[len(hairpinSeq)-1] {
		// not known terminal pair, nothing to close "hairpin"

		// TODO errors
		// raise RuntimeError
		panic("emap.COMPLEMENT[hairpin[0]] != hairpin[-1]")
	}

	dG := 0.0
	if emap.TRI_TETRA_LOOPS != nil {
		if en, ok := emap.TRI_TETRA_LOOPS[hairpinSeq]; ok {
			// it's a pre-known hairpin with known value
			d_h, d_s := en.H, en.S
			dG = d_g(d_h, d_s, temp)
		}
	}

	// add penalty based on size
	if en, ok := emap.HAIRPIN_LOOPS[hairpin_len]; ok {
		d_h, d_s := en.H, en.S
		dG += d_g(d_h, d_s, temp)
	} else {
		// it's too large, extrapolate
		en := emap.HAIRPIN_LOOPS[30]
		d_h, d_s := en.H, en.S
		d_g_inc := d_g(d_h, d_s, temp)
		dG += j_s(hairpin_len, 30, d_g_inc, temp)
	}

	// add penalty for a terminal mismatch
	en, ok := emap.TERMINAL_MM[pair]
	if hairpin_len > 3 && ok {
		d_h, d_s := en.H, en.S
		dG += d_g(d_h, d_s, temp)
	}

	// add penalty if length 3 and AT closing, formula 8 from SantaLucia, 2004
	if hairpin_len == 3 && (hairpinSeq[0] == 'A' || hairpinSeq[len(hairpinSeq)-1] == 'A') {
		dG += 0.5 // convert to entropy
	}

	return dG
}

// Traceback thru the V(i,j) and W(i,j) caches to find the structure
// For each step, get to the lowest energy W(i,j) within that block
// Store the structure in W(i,j)
// Inc i and j
// If the next structure is viable according to V(i,j), store as well
// Repeat
// Args:
//    i: The leftmost index to start searching in
//    j: The rightmost index to start searching in
//    v_cache: Energies where i and j bond
//    w_cache: Energies/sub-structures between or with i and j
// Returns:
//    A list of Structs in the final secondary structure
func traceback(i, j int, v_cache, w_cache [][]Struct) []Struct {
	// move i,j down-left to start coordinates
	s := w_cache[i][j]
	if !strings.Contains(s.Desc, "HAIRPIN") {
		for w_cache[i+1][j].Equal(s) {
			i += 1
		}
		for w_cache[i][j-1].Equal(s) {
			j -= 1
		}
	}

	structs := []Struct{}
	// ij := IJ{I: i, J: j}
	for {
		s = v_cache[i][j]

		// structs = append(structs, Struct{E: s.E, Desc: s.Desc, IJs: []IJ{ij}})
		structs = append(structs, Struct{E: s.E, Desc: s.Desc, IJ: []IJ{{I: i, J: j}}})

		// it's a hairpin, end of structure
		if len(s.IJ) == 0 {
			// set the energy of everything relative to the hairpin
			return trackback_energy(structs)
		}

		// it's a stack, bulge, etc
		// there's another single structure beyond this
		if len(s.IJ) == 1 {
			i, j = s.IJ[0].I, s.IJ[0].J
			// ij = s.IJs[0]
			continue
		}

		// it's a multibranch
		e_sum := 0.0
		structs = trackback_energy(structs)
		branches := []Struct{}
		for _, ij1 := range s.IJ {
			i1, j1 := ij1.I, ij1.J
			tb := traceback(i1, j1, v_cache, w_cache)
			if len(tb) > 0 && len(tb[0].IJ) > 0 {
				ij2 := tb[0].IJ[0]
				i2, j2 := ij2.I, ij2.J
				e_sum += w_cache[i2][j2].E
				branches = append(branches, tb...)
			}
		}

		last := structs[len(structs)-1]
		newlast := Struct{
			E:    RoundFloat(last.E-e_sum, 1),
			Desc: last.Desc,
			IJ:   make([]IJ, len(last.IJ)),
		}
		copy(newlast.IJ, last.IJ)

		structs[len(structs)-1] = newlast
		return append(structs, branches...)
	}

	// this is unreachable also in the python code
	// return trackback_energy(structs)
}

// Add energy to each structure, based on how it's W(i,j) differs from the one after
// Args:
//    structs: The structures for whom energy is being calculated
// Returns:
//    List[Struct]: Structures in the folded DNA with energy
func trackback_energy(structs []Struct) []Struct {
	structs_e := []Struct{}
	for index, str := range structs {
		e_next := 0.0
		if index < len(structs)-1 {
			e_next = structs[index+1].E
		}
		structs_e = append(structs_e, Struct{RoundFloat(str.E-e_next, 1), str.Desc, str.IJ})
	}
	return structs_e
}
