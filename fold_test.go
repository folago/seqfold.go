package seqfold

import (
	"math"
	"strings"
	"testing"

	"github.com/stretchr/testify/assert"
)

func TestFold(t *testing.T) {
	t.Run("FoldCache", func(t *testing.T) {
		seq := "ATGGATTTAGATAGAT"
		cache := dg_cache(seq, 37.0)
		seq_dg := dg(seq, 37.0)

		assert.InDelta(t, seq_dg, cache[0][len(seq)-1], 1)
	})
	t.Run("FoldDNA", func(t *testing.T) {
		// unafold's estimates for free energy estimates of DNA oligos
		unafold_dgs := map[string]float64{
			"GGGAGGTCGTTACATCTGGGTAACACCGGTACTGATCCGGTGACCTCCC":                         -10.94, // three branched structure
			"GGGAGGTCGCTCCAGCTGGGAGGAGCGTTGGGGGTATATACCCCCAACACCGGTACTGATCCGGTGACCTCCC": -23.4,  // four branched structure
			"CGCAGGGAUACCCGCG":                         -3.8,
			"TAGCTCAGCTGGGAGAGCGCCTGCTTTGCACGCAGGAGGT": -6.85,
			"GGGGGCATAGCTCAGCTGGGAGAGCGCCTGCTTTGCACGCAGGAGGTCTGCGGTTCGATCCCGCGCGCTCCCACCA": -15.50,
			"TGAGACGGAAGGGGATGATTGTCCCCTTCCGTCTCA":                                         -18.10,
			"ACCCCCTCCTTCCTTGGATCAAGGGGCTCAA":                                              -3.65,
		}

		for seq, ufold := range unafold_dgs {
			d := dg(seq, 37.0)
			// accepting a 60% difference
			delta := math.Abs(0.6 * min(d, ufold))
			assert.InDeltaf(t, d, ufold, delta, seq)
		}
	})
	t.Run("FoldRNA", func(t *testing.T) {
		// unafold's estimates for free energy estimates of RNA oligos
		// most tests available at https://github.com/jaswindersingh2/SPOT-RNA/blob/master/sample_inputs/batch_seq.fasta
		unafold_dgs := map[string]float64{
			"ACCCCCUCCUUCCUUGGAUCAAGGGGCUCAA":        -9.5,
			"AAGGGGUUGGUCGCCUCGACUAAGCGGCUUGGAAUUCC": -10.1,
			"UUGGAGUACACAACCUGUACACUCUUUC":           -4.3,
			"AGGGAAAAUCCC":                           -3.3,
			"GCUUACGAGCAAGUUAAGCAAC":                 -4.6,
			"UGGGAGGUCGUCUAACGGUAGGACGGCGGACUCUGGAUCCGCUGGUGGAGGUUCGAGUCCUCCCCUCCCAGCCA":   -32.8,
			"GGGCGAUGAGGCCCGCCCAAACUGCCCUGAAAAGGGCUGAUGGCCUCUACUG":                         -20.7,
			"GGGGGCAUAGCUCAGCUGGGAGAGCGCCUGCUUUGCACGCAGGAGGUCUGCGGUUCGAUCCCGCGCGCUCCCACCA": -31.4,
		}

		for seq, ufold := range unafold_dgs {
			d := dg(seq, 37.0)

			// accepting a 5% difference
			delta := math.Abs(0.5 * min(d, ufold))
			assert.InDeltaf(t, d, ufold, delta, seq)
		}
	})
	t.Run("DotBracket", func(t *testing.T) {
		// t.Skip()
		seq := "GGGAGGTCGTTACATCTGGGTAACACCGGTACTGATCCGGTGACCTCCC"
		structs := Fold(seq, 37.0)

		assert.Equal(t, "((((((((.((((......))))..((((.......)))).))))))))", DotBracket(structs))
	})
	t.Run("Multibranch", func(t *testing.T) {
		seq := "GGGAGGTCGTTACATCTGGGTAACACCGGTACTGATCCGGTGACCTCCC" // three branch

		structs := Fold(seq, 37.0)

		found := false
		foundIJ := IJ{7, 41}
		for _, s := range structs {
			if strings.Contains(s.Desc, "BIFURCATION") {
				for _, ij := range s.IJ {
					if ij == foundIJ {
						found = true
					}
				}
			}
		}
		assert.True(t, found, "not found a  BIFURCATION with (7, 41) in ij")
	})
	t.Run("Pair", func(t *testing.T) {
		seq := "ATGGAATAGTG"
		assert.Equal(t, pair(seq, 0, 1, 9, 10), "AT/TG")
	})
	t.Run("Stack", func(t *testing.T) {
		seq := "GCUCAGCUGGGAGAGC"
		temp := 310.15

		e := stack(seq, 1, 2, 14, 13, temp, RNA_ENERGIES)
		assert.InDelta(t, e, -2.1, 0.1)
	})
	t.Run("Bulge", func(t *testing.T) {
		// mock bulge of CAT on one side and AG on other
		// from pg 429 of SantaLucia, 2004
		seq := "ACCCCCATCCTTCCTTGAGTCAAGGGGCTCAA"

		pair_dg := bulge(seq, 5, 7, 18, 17, 310.15, DNA_ENERGIES)
		assert.InDelta(t, pair_dg, 3.22, 0.4)
	})
	t.Run("Hairpin", func(t *testing.T) {
		// hairpin = "CCTTGG"
		seq := "ACCCCCTCCTTCCTTGGATCAAGGGGCTCAA"
		i := 11
		j := 16
		temp := 310.15
		hairpin_dg := hairpin(seq, i, j, temp, DNA_ENERGIES)
		// this differs from Unafold
		assert.InDelta(t, hairpin_dg, 4.3, 1.0)

		// from page 428 of SantaLucia, 2004
		// hairpin = "CGCAAG"
		seq = "ACCCGCAAGCCCTCCTTCCTTGGATCAAGGGGCTCAA"
		i = 3
		j = 8
		hairpin_dg = hairpin(seq, i, j, temp, DNA_ENERGIES)
		assert.InDelta(t, hairpin_dg, 0.67, 0.1)

		seq = "CUUUGCACG"
		i = 0
		j = 8
		hairpin_dg = hairpin(seq, i, j, temp, RNA_ENERGIES)
		assert.InDelta(t, hairpin_dg, 4.5, 0.2)
	})
	t.Run("InternalLoop", func(t *testing.T) {
		seq := "ACCCCCTCCTTCCTTGGATCAAGGGGCTCAA"
		i := 6
		j := 21
		temp := 310.15
		dg := internal_loop(seq, i, i+4, j, j-4, temp, DNA_ENERGIES)
		assert.InDelta(t, dg, 3.5, 0.1)
	})
	t.Run("W", func(t *testing.T) {
		seq := "GCUCAGCUGGGAGAGC"
		i := 0
		j := 15
		temp := 310.15
		v_cache, w_cache := buildDefaultCaches(seq)
		struc := w(seq, i, j, temp, v_cache, w_cache, RNA_ENERGIES)
		assert.InDelta(t, struc.E, -3.8, 0.2)

		seq = "CCUGCUUUGCACGCAGG"
		i = 0
		j = 16
		temp = 310.15
		v_cache, w_cache = buildDefaultCaches(seq)
		struc = w(seq, i, j, temp, v_cache, w_cache, RNA_ENERGIES)
		assert.InDelta(t, struc.E, -6.4, 0.2)

		seq = "GCGGUUCGAUCCCGC"
		i = 0
		j = 14
		v_cache, w_cache = buildDefaultCaches(seq)
		struc = w(seq, i, j, temp, v_cache, w_cache, RNA_ENERGIES)
		assert.InDelta(t, struc.E, -4.2, 0.2)
	})
}
