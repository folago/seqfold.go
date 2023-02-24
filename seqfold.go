package seqfold

import (
	"fmt"
	"math"
)

type Comp map[byte]byte

type MultiBranch struct {
	A, B, C, D float64
}

// Energy holds two energies, enthaply and entropy
// SantaLucia & Hicks (2004), Annu. Rev. Biophys. Biomol. Struct 33: 415-440
type Energy struct {
	// enthapl
	H float64
	// entropy
	S float64
}

type BpEnergy map[string]Energy

type LoopEnergy map[int]Energy

type Energies struct {
	BULGE_LOOPS     LoopEnergy
	COMPLEMENT      Comp
	DE              BpEnergy
	HAIRPIN_LOOPS   LoopEnergy
	MULTIBRANCH     MultiBranch
	INTERNAL_LOOPS  LoopEnergy
	INTERNAL_MM     BpEnergy
	NN              BpEnergy
	TERMINAL_MM     BpEnergy
	TRI_TETRA_LOOPS BpEnergy
}

type IJ struct {
	I, J int
}

// A single structure with a free energy, description, and inward children
type Struct struct {
	E    float64
	Desc string
	IJ   []IJ
}

func (s Struct) Equal(other Struct) bool {
	if len(s.IJ) != len(other.IJ) {
		return false
	}
	for i, val := range s.IJ {
		if val != other.IJ[i] {
			return false
		}
	}
	return s.E == other.E
}

func (s Struct) Valid() bool {
	return s.E != math.Inf(1) && s.E != math.Inf(-1)
}

func (s Struct) String() string {
	i, j := " ", " "
	if len(s.IJ) > 0 {
		i, j = fmt.Sprint(s.IJ[0].I), fmt.Sprint(s.IJ[0].J)
	}
	return fmt.Sprintf("(%s, %s) % 6.2f %s", i, j, s.E, s.Desc)
}

// MultiString returns all the fields as strings, this is useful to output the
// result in a tabular fashion using the same format string.
func (s Struct) MultiString() (string, string, string, string) {
	i, j := "", ""
	if len(s.IJ) > 0 {
		i, j = fmt.Sprint(s.IJ[0].I), fmt.Sprint(s.IJ[0].J)
	}
	return i, j, fmt.Sprintf("%6.2f", s.E), s.Desc
}

// STRUCT_DEFAULT = Struct(-math.inf)
var STRUCT_DEFAULT = Struct{
	Desc: "",
	E:    math.Inf(-1),
}

// STRUCT_NULL = Struct(math.inf)
var STRUCT_NULL = Struct{
	Desc: "",
	E:    math.Inf(1),
}
