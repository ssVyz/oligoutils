package oligoutils

import "testing"
import "fmt"

func TestParseError(t *testing.T) {

	_, err := ParseFasta("not_exist.txt")
	if err != nil {
		fmt.Println(err)
	}
}

func TestIsValidBase(t *testing.T) {

	if isValidBase('A') != true {
		t.Fatalf("Valid for A is incorrect")
	}
	if isValidBase('Z') != false {
		t.Fatalf("Valid for A is incorrect")
	}
	if isValidBase('_') != false {
		t.Fatalf("Valid for A is incorrect")
	}

}

func TestParseFasta(t *testing.T) {

	content, _ := ParseFasta("example/exa1.fa")

	fmt.Printf("Reading first entry name: %v, and first entry sequence: %v\nLength of the slice: %v\n", content[0].Header, content[0].Seq, len(content))
}

func TestParseFastaMalformed(t *testing.T) {
	_, err := ParseFasta("example/malformed.txt")
	if err != nil {
		fmt.Println(err)
	} else {
		t.Fatalf("malformed did not result in an error")
	}
}

func TestBaseMatching(t *testing.T) {

	bases := "ART"
	template := "TCT"

	fmt.Println("Matching three bases:")
	fmt.Printf("input base %v, template base %v, matching: %v\n", string(bases[0]), string(template[0]), isComplementMatch(bases[0], template[0]))
	fmt.Printf("input base %v, template base %v, matching: %v\n", string(bases[1]), string(template[1]), isComplementMatch(bases[1], template[1]))
	fmt.Printf("input base %v, template base %v, matching: %v\n", string(bases[2]), string(template[2]), isComplementMatch(bases[2], template[2]))

}

func TestBaseIdentical(t *testing.T) {

	bases := "WRT"
	template := "AGA"

	fmt.Println("Matching three bases:")
	fmt.Printf("input base %v, template base %v, matching: %v\n", string(bases[0]), string(template[0]), isIdentical(bases[0], template[0]))
	fmt.Printf("input base %v, template base %v, matching: %v\n", string(bases[1]), string(template[1]), isIdentical(bases[1], template[1]))
	fmt.Printf("input base %v, template base %v, matching: %v\n", string(bases[2]), string(template[2]), isIdentical(bases[2], template[2]))

}

func TestOligoMatch(t *testing.T) {

	oligo1 := "WGTTRCCTGA"
	oligo2 := "AGTTGCCTGA"
	oligo3 := "TGTTCCCTGA"
	oligo4 := "AATGTCA"
	oligo5 := "AATGTCA"

	fmt.Println("Matching 3 oligo combinations:")
	fmt.Printf("Oligo1 %v, Oligo2 %v, matching: %v\n", oligo1, oligo2, OligoMatch(oligo1, oligo2))
	fmt.Printf("Oligo1 %v, Oligo2 %v, matching: %v\n", oligo1, oligo3, OligoMatch(oligo1, oligo3))
	fmt.Printf("Oligo1 %v, Oligo2 %v, matching: %v\n", oligo4, oligo5, OligoMatch(oligo4, oligo5))

}

func TestMakeReverseOligo(t *testing.T) {

	oligo1 := "ATAGTGCCCRGVIGTG"
	oligo2 := "ATAGTTC-CTA"

	res1, err1 := MakeReverseComplement(oligo1)
	res2, err2 := MakeReverseComplement(oligo2)

	fmt.Println("Testing reverse complement:")
	fmt.Printf("Oligo %v, reverse complement %v, error %v\n", oligo1, res1, err1)
	fmt.Printf("Oligo %v, reverse complement %v, error %v\n", oligo2, res2, err2)
}

func TestIsCanonBase(t *testing.T) {

	oligo1 := "CURRYWURST"

	for i := 0; i < len(oligo1); i++ {
		fmt.Printf("Result for %v is %v\n", oligo1[i], isCanonBase(oligo1[i]))
	}
}

func TestIsCanonOligo(t *testing.T) {
	seq1 := Seqr{Header: "Testseq1", Seq: "AAAGGGTACCCATTGTC",}
	seq2 := Seqr{Header: "Testseq2", Seq: "AAAGRGTACCCATTGTC",}

	fmt.Printf("Testing canon oligo with %v, result: %v\n", seq1.Header, IsCanonOligo(seq1))
	fmt.Printf("Testing canon oligo with %v, result: %v\n", seq2.Header, IsCanonOligo(seq2))
}

func TestCleanList(t *testing.T) {
	seq3 := Seqr{Header: "Testseq1", Seq: "AAAGGGTACCCATTGTC",}
	seq4 := Seqr{Header: "Testseq2", Seq: "AAAGRGTACCCATTGTC",}
	seq5 := Seqr{Header: "Testseq2", Seq: "AAAGAGTACCCATTGTC",}

	lst := []Seqr{}
	lst = append(lst, seq3)
	lst = append(lst, seq4)
	lst = append(lst, seq5)

	cleanList, elim := CleanSeqList(lst)

	fmt.Printf("Length of input list %v, length of output list %v, eliminated: %v\n", len(lst), len(cleanList), elim) 
}
