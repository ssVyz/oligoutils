
//// Oligoutils
//// version 0.2 - 260903
////

package oligoutils

import "bufio"
import "os"
import "io"
import "fmt"
import "strings"
import "slices"


type Seqr struct {
	Header string
	Seq string
}

type posBase struct {
	bases []byte
}


const Iupac = "ATGCURYWVBMHDSKNI"


var complements = map[string]string{
	"A": "T",
	"T": "A",
	"G": "C",
	"C": "G",
	"U": "A",

	"R": "Y",
	"Y": "R",
	"S": "S",
	"W": "W",
	"K": "M",
	"M": "K",

	"B": "V",
	"V": "B",
	"D": "H",
	"H": "D",

	"N": "N",
	"I": "N",
}


var matches = map[string][]string{
	"A": {"T"},
	"T": {"A"},
	"G": {"C"},
	"C": {"G"},
	"U": {"A"},

	"R": {"T", "C"},
	"Y": {"G", "A"},
	"S": {"C", "G"},
	"W": {"T", "A"},
	"K": {"C", "A"},
	"M": {"T", "G"},

	"B": {"G", "C", "A"},
	"D": {"T", "C", "A"},
	"H": {"T", "G", "A"},
	"V": {"T", "G", "C"},

	"N": {"A", "T", "G", "C"},
	"I": {"A", "T", "G", "C"},
}

var identical = map[string][]string{
	"A": {"A"},
	"T": {"T"},
	"G": {"G"},
	"C": {"C"},
	"U": {"T"},

	"R": {"A", "G"},
	"Y": {"C", "T"},
	"S": {"C", "G"},
	"W": {"A", "T"},
	"K": {"G", "T"},
	"M": {"A", "C"},

	"B": {"C", "G", "T"},
	"D": {"A", "G", "T"},
	"H": {"A", "C", "T"},
	"V": {"A", "C", "G"},

	"N": {"A", "C", "G", "T"},
	"I": {"A", "C", "G", "T"},
}


/// Oligo ops area. Functions that act on entire oligos
// reverse complement an oligo
func MakeReverseComplement(oligo string) (string, error) {
	oligoLength := len(oligo)
	idx := oligoLength - 1
	var result string

	for i := 0; i < oligoLength; i++ {
		currentPos := oligo[idx]
		if !isValidBase(currentPos) {return "", fmt.Errorf("Invalid base: %v", currentPos)}
		result = result + complements[string(currentPos)]
		idx--
	}
	return result, nil
} 

// Oligo comparer. Checks if two oligos are identical, IUPAC aware for query
func OligoMatch(query string, template string) bool {
	if len(query) != len(template) {
		return false
	}
	for i := 0; i < len(query); i++ {
		if !isIdentical(query[i], template[i]) {return false;}
	}
	return true
}

// Returns false if this Seqr contains a non-cannonical base
func IsCanonOligo(oli Seqr) bool {
	if len(oli.Seq) == 0 {return false;}
	result := true
	for i := 0; i < len(oli.Seq); i++ {
		if !isCanonBase(oli.Seq[i]) {result = false;}
	}
	return result
}



/// Per position area. Functions that compare or act on individual bases.
// is this position identical? IUPAC aware for query, i.e. bas
func isIdentical(bas byte, template byte) bool {
	if isValidBase(bas) != true || isValidBase(template) != true {
		return false
	}
	basConv := strings.ToUpper(string(bas))
	templConv := strings.ToUpper(string(template))
	matching := false

	lst, ok := identical[basConv]
	if !ok {return false;}

	for _, item := range lst {
		if item == templConv {
			matching = true
		}
	}
	return matching
}

// Is this a valid base pairing? both as bytes
func isComplementMatch(bas byte, template byte) bool {
	if isValidBase(bas) != true || isValidBase(template) != true {
		return false
	}
	basConv := strings.ToUpper(string(bas))
	templConv := strings.ToUpper(string(template))
	matching := false

	lst, ok := matches[basConv]
	if !ok {return false;}

	for _, item := range lst {
		if item == templConv {
			matching = true
		}
	}
	return matching
}

// checks if a byte character represents a valid IUPAC base
func isValidBase(b byte) bool {
	isValid := false
	bString := strings.ToUpper(string(b))
	for i := 0; i < len(Iupac); i++ {
		if bString == string(Iupac[i]) {
			isValid = true
		}
	}
	return isValid
}

// Check if this byte is a cannonical base
func isCanonBase(bas byte) bool {
	switch bas {
	case 'A', 'T', 'G', 'C':
		return true
	}
	return false
}

func collapseIupac(bases []byte) (byte, error) {
	if len(bases) == 0 || len(bases) > 4 {return 0, fmt.Errorf("invalid base list");}

	// make string representation of the slice in alphabetical order
	var stringSlice = []string{}
	for _, base := range bases {
		if base == byte('A') {stringSlice = append(stringSlice, "A");}
	}
	for _, base := range bases {
		if base == byte('C') {stringSlice = append(stringSlice, "C");}
	}
	for _, base := range bases {
		if base == byte('G') {stringSlice = append(stringSlice, "G");}
	}
	for _, base := range bases {
		if base == byte('T') {stringSlice = append(stringSlice, "T");}
	}

	// Find the entry that matches the stringSlice. Uses slices.Equal to compare two slices
	var result byte = 0
	for key, slce := range identical {
		if slices.Equal(slce, stringSlice) && key != "U" && key != "I" {
			result = key[0]
		}
	}

	if result == 0 {
		return 0, fmt.Errorf("IUPAC character not found");
	} else {
		return result, nil
	}
}



/// File ops. Reading and writing of fasta files
// Load a fasta file and turn it as a slice of Seqr
func ParseFasta(path string) ([]Seqr, error) {
	// initialize the result slice
	var res []Seqr

	//file opening and assigning reader
	f, err := os.Open(path)
	if err != nil {
		return res, fmt.Errorf("Error while trying to open a file")
	}
	defer f.Close()

	r := bufio.NewReader(f)

	// Go into the Seqr building loop

	firstPass := true
	buildHeader := true
	var currentRec = Seqr{}

	for true {
		// read a byte and handle EOF or errors
		c, err := r.ReadByte()
		if err == io.EOF {
			if firstPass == false && buildHeader == false {
				res = append(res, currentRec)
			}
			break
		} else if err != nil {
			return []Seqr{}, fmt.Errorf("Error while reading input: %w", err)
		}
		
		// Handle first pass > symbol OR a new sequence with reset
		if c == '>' && firstPass == true {
			firstPass = false
			continue
		} else if c != '>' && firstPass == true {
			return []Seqr{}, fmt.Errorf("fasta formatting, initial position is not >")
		} else if c == '>' && firstPass == false {
			res = append(res, currentRec)
			currentRec = Seqr{}
			buildHeader = true
			continue
		}

		// Handle backslash n and or backslash r delimiting the header
		if (c == '\n' || c == '\r') && buildHeader == true {
			buildHeader = false
			continue
		}

		// handle any letter or character
		if buildHeader == true {
			currentRec.Header = currentRec.Header + string(c)
			continue
		} else if buildHeader == false {
			if isValidBase(c) == true {
				currentRec.Seq = currentRec.Seq + strings.ToUpper(string(c))
			}
			continue
		}
	}
	return res, nil
}


func WriteFasta(path string, sl []Seqr) (error) {
	if len(sl) == 0 {
		return fmt.Errorf("The provided sequence list is empty\n")
	}

	f, err := os.Create(path)
	if err != nil {
		return fmt.Errorf("Error while trying to create file: %v\n", path)
	}
	defer f.Close()

	w := bufio.NewWriter(f)
	defer w.Flush()

	// write Seqr entries one by one
	for i := 0; i < len(sl); i++ {
		_, err = fmt.Fprintf(w, ">%v\n%v\n\n", sl[i].Header, sl[i].Seq)
		if err != nil {
			return fmt.Errorf("Error while writing Sequence entry %v.\n", i)
		}
	}
	return nil
}


/// List area. Functions that act on entire lists
// takes a Seqr list and returns a clean version without ambiguity letters
func CleanSeqList(sl []Seqr) ([]Seqr, int) {
	if len(sl) == 0 {return []Seqr{}, 0;}
	
	result := []Seqr{}
	eliminated := 0

	for i := 0; i < len(sl); i++ {
		if IsCanonOligo(sl[i]) {
			result = append(result, sl[i])
		} else {
			eliminated++
		}
	}
	return result, eliminated
}

// Checks if all Seqr have the same length. Returns the length if yes. Otherwise returns an error
func SeqListLength(sl []Seqr) (int, error) {
	if len(sl) == 0 {return 0, fmt.Errorf("Can not find consensus length, Seq list is empty");}
	first := len(sl[0].Seq)

	for i := 0; i < len(sl); i++ {
		if len(sl[i].Seq) != first {
			return 0, fmt.Errorf("Not all sequences are identical in length")
		}
	}
	return first, nil
}

// Collapse a seqr list into a consensus sequence
func MakeConsensus(sl []Seqr) (string, error) {
	leng, err := SeqListLength(sl)
	if err != nil {
		return "", err
	}

	var resultString string = ""
	var resultSlce = []byte{}
	var positions = []posBase{}

	// populate the positions slice with the variants for each position
	for i := 0; i < leng; i++ {
		posBuild := posBase{}
		for j := 0; j < len(sl); j++ {
			if sl[j].Seq[i] == 'A' {
				posBuild.bases = append(posBuild.bases, 'A')
				break
			}
		}
		for j := 0; j < len(sl); j++ {
			if sl[j].Seq[i] == 'C' {
				posBuild.bases = append(posBuild.bases, 'C')
				break
			}
		}
		for j := 0; j < len(sl); j++ {
			if sl[j].Seq[i] == 'G' {
				posBuild.bases = append(posBuild.bases, 'G')
				break
			}
		}
		for j := 0; j < len(sl); j++ {
			if sl[j].Seq[i] == 'T' {
				posBuild.bases = append(posBuild.bases, 'T')
				break
			}
		}

		positions = append(positions, posBuild)
	}

	//collapse the positions into iupacs
	for i := 0; i < len(positions); i++ {
		var newBase byte
		newBase, err := collapseIupac(positions[i].bases)
		if err != nil {
			return "", err
		}

		if newBase != 0 {
			resultSlce = append(resultSlce, newBase)
		}
	}

	resultString = string(resultSlce)
	return resultString, nil

}

