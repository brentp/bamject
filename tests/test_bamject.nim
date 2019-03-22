# vim: ts=2 sw=2 sts=2 sr et ai si
import bamject
import hts/bam
import strutils
import sequtils
import unittest

# met refgetter for testing.
type RG = object
  sequence: string
  offset: int

proc get(r:RG, region:string): string =
  var pos = region.split(':')
  var se = pos[1].split('-')
  return r.sequence[parseInt(se[0])-1-r.offset..<parseInt(se[1])-r.offset]

suite "bamjection":
  test "test simple snp injection":
    var r = Read(chrom:"1", start: 1, sequence: "AAAAA", cigar:"5M".toCigar)
    var m = Mutation(chrom:"1", start: 3, alt: "C", reference:"A")
    var reference = RG(sequence:"AAAAAA")
    discard m.mutate(r, reference)
    check r.sequence == "AACAA"

  test "test simple insertion injection":
    var r = Read(chrom:"1", start: 0, sequence: "AAAAA", cigar:"5M".toCigar)
    var m = Mutation(chrom:"1", start: 2, alt: "ATC", reference:"A")
    var reference = RG(sequence:"AAAAAA")
    discard m.mutate(r, reference)
    doAssert r.sequence == "AAATC"
    check $r.cigar == "3M2I"

    r = Read(chrom:"1", start: 0, sequence: "AAAAA", cigar:"5M".toCigar)
    m = Mutation(chrom:"1", start: 3, alt: "ATC", reference:"A")
    discard m.mutate(r, reference)
    doAssert r.sequence == "AAAAT"
    check $r.cigar == "4M1I"

    r = Read(chrom:"1", start: 4, sequence: "AAAAAC", cigar: "5M1I".toCigar)
    m = Mutation(chrom:"1", start: 3, alt: "ATC", reference:"A")
    discard m.mutate(r, reference)
    check r.sequence == "TCAAAA"
    check $r.cigar == "2I4M"

    r = Read(chrom:"1", start: 4, sequence: "AAAAAC", cigar: "4M1I1M".toCigar)
    m = Mutation(chrom:"1", start: 6, alt:    "ATC", reference:"A")
                                            #3M2I1M
    discard m.mutate(r, reference)
    check r.sequence == "AAATCA"
    check $r.cigar == "3M2I1M"

  test "test simple deletion injection":
    var r = Read(chrom:"1", start: 0, sequence: "AAAAACT", cigar:"7M".toCigar)
    var m = Mutation(chrom:"1", start: 2, alt:     "A", reference:"AAA")
    var reference = RG(sequence:"AAAAACTGG")
    discard m.mutate(r, reference)
    check r.sequence == "AAACTGG"
    check $r.cigar == "3M2D4M"

  test "long insertion":
    var m = Mutation(chrom: "10", start: 43607600, alt: "TGGAGTGTGA", reference: "T")
    var r = Read(chrom: "10", start: 43607570-1, cigar: "98M39S".toCigar, sequence: "CCCCTGTCCTGTGCAGTCAGCAAGAGACGGCTGGAGTGTGAGGAGTGTGGCGGCCTGGGCTCCCCAACAGGCAGGTGTGAGTGGAGGCAAGGAGATGGAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTTTAGCG")
    var reference = RG(sequence:"TGGAGTGTGAGGAGTGTGGCGGCCTGGGCTCCCCAACAGGCAGGTGTGAGTGGAGGCAAGGAGATGGCAAAGGTAAGCCCTGGAAACGCCCAAGGGAGGC", offset: 43607600-1)
    #echo $r.cigar
    #echo r.sequence
    discard m.mutate(r, reference)
    #echo r.sequence
    #echo $r.cigar
    check $r.cigar == "32M9I66M30S"

  test "another long insertion":
    var m = Mutation(chrom: "10", start: 43607600, alt: "TGGAGTGTGA", reference: "T")
    var r = Read(chrom: "10", start: 43607578-1, cigar: "41S96M".toCigar, sequence: "CTTGGCCCATACACTCTTTCCCTACACGACGCTCTTCCGATCTGTGCAGTCAGCAAGAGACGGCTGGAGTGTGAGGAGTGTGGCGGCCTGGGCTCCCCAACAGGCAGGTGTGAGTGGAGGCAAGGAGATGGCAAAGG")
    var reference = RG(sequence:"TGGAGTGTGAGGAGTGTGGCGGCCTGGGCTCCCCAACAGGCAGGTGTGAGTGGAGGCAAGGAGATGGCAAAGGTAAGCCCTGGAAACGCCCAAGGGAGGC", offset: 43607600-1)
    #echo $r.cigar
    #echo r.sequence
    discard m.mutate(r, reference)
    #echo r.sequence
    check $r.cigar == "41S24M9I63M"

  test "at end":
    var m = Mutation(chrom: "10", start: 89624240, alt: "CA", reference: "C")
    var r = Read(chrom: "10", start: 89624148-1, cigar: "2S135M".toCigar, sequence: "CTCATCCTGCAGAAGAAGCCCCGCCACCAGCAGCTTCTGCCATCTCTCTCCTCCTTTTTCTTCAGCCACAGGCTCCCAGACATGACAGCCATCATCAAAGAGATCGTTAGCAGAAACAAAAGGAGATATCAAGAGGA")
    var reference = RG(sequence:"")
    discard m.mutate(r, reference)
    check $r.cigar == "2S94M1I40M"
  test "missed soft":
    var m = Mutation(chrom: "10", start: 89624240, alt: "C", reference: "CAA")
    var r = Read(chrom: "10", start: 89624242-1, cigar: "51S86M".toCigar, sequence: "NNNNNACACGNGCAGTTGACACTCTTTCCCTACACGACGCTCTTCCGATCTAAAGAGATCGTTAGCAGAAACAAAAGGAGATATCAAGAGGATGGATTCGACTTAGACTTGACCTGTATCCATTTCTGCGGCTGCTC")
    var reference = RG(sequence:"TCAAAGAGATCGTTAGCAGAAACAAAAGGAGATATCAAGAGGATGGATTCGACTTAGACTTGACCTGTATCCATTTCTGCGGCTGCTCCTCTTTACCTTTCTGTCACTCTCTTAGAACGTGGGAGTAGACGGATGCGAAAA", offset:89624240-1)
    echo $r.cigar
    discard m.mutate(r, reference)
    # TODO: make sure this is right.
    check $r.cigar == "51S2D86M"


  test "another case":
    var m = Mutation(chrom: "10", start: 89685265, alt: "T", reference: "TAAGG")
    var reference = RG(sequence:"TCAAAGAGATCGTTAGCAGAAACAAAAGGAGATATCAAGAGGATGGATTCGACTTAGACTTGACCTGTATCCATTTCTGCGGCTGCTCCTCTTTACCTTTCTGTCACTCTCTTAGAACGTGGGAGTAGACGGATGCGAAAA", offset:89685265-1)
    var r = Read(chrom:"10", start:89685158-1, cigar:"96M1I40M".toCigar, sequence:"ATTAGGAAAAAGAAAATCTGTCTTTTGGTTTTTCTTGATAGTATTAATGTAATTTCAAATGTTAGCTCATTTTTGTTAATGGTGGCTTTTTGTTTGTTTTGTTTTGTTTTAAGGTTTTTGGATTCAAAGCATAAAAA")
    echo m
    echo $r.cigar
    discard m.mutate(r, reference)
    echo $r.cigar

  test "casexx":
    var m = Mutation(chrom: "11", start: 533847, alt: "GGTCCCGCATGGCGCTGTACTC", reference: "G")
    var r = Read(chrom: "11", start: 533787-1, cigar: "43M1D94M".toCigar, sequence:"AAAGACTTGGTGTTGTTGATGGCAAACACACACAGGAAGCCCTCCCGGTGCGCATGTACTGGTCCCGCATGGCGCTGTACTCCTCCTGGCCGGCGGTATCCAGGATGTCCAACAGGCACGTCTCCCCATCAATGACC")
    var reference = RG(sequence:"")
    discard m.mutate(r, reference)
    echo $r.cigar

  test "snp in complex cigar":
    var m = Mutation(chrom: "13", start: 28602339, alt: "C", reference: "G")
    var r = Read(chrom: "13", start: 28602200-1, cigar: "27M4D110M".toCigar, sequence: "TGACTGAGAAAAGACAAAGAATTAAAAAGAGAGAGAGAGAGAGAGCAAACATTCTCTTTGTCATCAAGCTACAGTCTTTTTGATGAGGTGATTTTCGTGGAAGTGGGTTACCTGACAGTGTGCACGCCCCCAGCAGG")

    discard m.mutate_snp(r)
    echo $r.cigar

  test "end insertion":
    var m = Mutation(chrom: "10", start: 89624240, alt: "CA", reference: "C")
    var r = Read(chrom: "10", start: 89624104-1, cigar: "96M1D41M".toCigar, sequence:"GCAGAGCGAGGGGCATCAGCTACCGCCAAGTCCAGAGCCATTTCCATCCTGCAGAAGAAGCCCCGCCACCAGCAGCTTCTGCCATCTCTCTCCTCCTTTTCTTCAGCCACAGGCTCCCAGACATGACAGCCATCATC")
    var reference = RG(sequence:"")
    discard m.mutate(r, reference)
    echo $r.cigar

  test "another end insertion":
    var m = Mutation(chrom: "10", start: 89624240, alt: "CA", reference: "C")
    var r = Read(chrom: "10", start: 89624113-1, cigar:"36M1D96M5S".toCigar, sequence:"GGGGCATCAGCTACCGCCAAGTCCAGAGCCATTTCCTCCTGCAGAAGAAGCCCCGCCACCAGCAGCTTCTGCCATCTCTCTCCTCCTTTTTCTTCAGCCACAGGCTCCCAGACATGACAGCCATCATCAAAGNNNNN")
    var reference = RG(sequence:"")
    discard m.mutate(r, reference)
    echo $r.cigar

  test "anoter deletion":
    var r = Read(chrom: "10", start: 89685129-1, cigar:"117M4D20M".toCigar, sequence:"TATTTATTTTCACTTAAATGGTATTTGAGATTAGGAAAAAGAAAATCTGTCTTTTGGTTTTTCTTGATAGTATTAATGTAATTTCAAATGTTAGCTCATTTTTGTTAATGGTGGCTTTTTGTTTGTTTTGTTTTAAG")
    var m = Mutation(chrom: "10", start: 89685265, alt: "T", reference: "TAAGG")
    var reference = RG(sequence:"TCAAAGAGATCGTTAGCAGAAACAAAAGGAGATATCAAGAGGATGGATTCGACTTAGACTTGACCTGTATCCATTTCTGCGGCTGCTCCTCTTTACCTTTCTGTCACTCTCTTAGAACGTGGGAGTAGACGGATGCGAAAA", offset:89685265-1)
    discard m.mutate(r, reference)
    echo $r.cigar

  test "looong insertion":
    var m = Mutation(chrom: "11", start: 533847, alt: "GGTCCCGCATGGCGCTGTACTC", reference: "G")
    var r = Read(chrom:"11", start: 533734-1, cigar: "121M16S".toCigar, sequence:"GTGGGCTCCCGGGCCAGCCTCACGGGGTTCACCTGTACTGGTGGATGTCCTCAAAAGACTTGGTGTTGTTGATGGCAAACACACACAGGAAGCCCTCCCCGGTGCGCATGTACTGGTCCCGNNNNNNNNNNNNNNNN")
    var reference = RG(sequence:"")
    discard m.mutate(r, reference)
    echo $r.cigar
