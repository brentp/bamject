# vim: ts=2 sw=2 sts=2 sr et ai si
import bamject
import hts/bam
import strutils
import unittest

# met refgetter for testing.
type RG = object
  sequence: string

proc get(r:RG, region:string): string =
  var pos = region.split(':')
  var se = pos[1].split('-')
  return r.sequence[parseInt(se[0])-1..<parseInt(se[1])]

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
