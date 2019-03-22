import hts
import lapper
import random
import strformat
import algorithm
import sequtils
import argparse
import ./cigar
import strutils
import tables

type RefGetter* = concept i
  ## An object/tuple must implement these 2 methods to use this module
  get(i, string) is string

type Mutation = object
  chrom: string
  start: int
  alt: string
  reference: string

# satisfy nim-lapper interface.
proc start(m:Mutation): int {.inline.} =
  return m.start
proc stop(m:Mutation): int {.inline.} =
  return m.start + m.reference.len

type Read = object
  chrom: string
  start: int
  sequence: string
  cigar: seq[CigarElement]


proc is_snp(m:Mutation): bool {.inline.} =
  return m.alt.len == 1 and m.reference.len == 1

proc is_insertion(m:Mutation): bool {.inline.} =
  return m.alt.len > m.reference.len

proc is_deletion(m:Mutation): bool {.inline.} =
  return m.reference.len > m.alt.len

proc mutate_snp(m:Mutation, r:var Read): bool =
  doAssert m.is_snp
  var offset = m.start - r.start
  r.sequence[offset] = m.alt[0]
  return true

proc mutate_insertion(m:Mutation, r:var Read, reference:RefGetter): bool =
  doAssert m.is_insertion
  var offset = m.start - r.start + 1
  var alt = m.alt[1..m.alt.high]
  var L = r.sequence.len
  if offset < 0:
    alt = alt[abs(offset)..alt.high]
    offset = 0

  var sav = r.sequence
  r.sequence = r.sequence[0..<offset] & alt
  if r.sequence.len < L:
    r.sequence &= sav[offset..sav.high]
  r.sequence = r.sequence[0..<L]
  doAssert r.sequence.len == L
  var coff = 0

  var cigarCopy = newSeqOfCap[CigarElement](r.cigar.len)
  var added = false

  for i, op in r.cigar:
    if op.consumes.query:
      coff += op.len
    if coff < offset or added:
      if i == r.cigar.high:
        # truncate last cigar
        var n = L - (coff - op.len)
        cigarCopy.add(($n & $op.op).toCigar)
      else:
        cigarCopy.add(op)
      if coff >= L: break
      continue
    if not op.consumes.reference:
      continue
    added = true
    #echo $op.op
    var before = coff - op.len + offset
    if before > 0:
      cigarCopy.add(($before & $op.op).toCigar)
    var left = op.len - (L - op.len + before)

    if left + before < op.len:
      cigarCopy.add(($(alt.len) & 'I').toCigar)
      var n = if i == r.cigar.len - 1:
        op.len - alt.len
      else:
        op.len - before
      if coff - op.len + before + n + alt.len > L:
        n = L - (coff - op.len + before + alt.len)
      cigarCopy.add(($n & $op.op).toCigar)
      break
    else:
      cigarCopy.add(($left & 'I').toCigar)

  r.cigar = cigarCopy

proc mutate_deletion(m:Mutation, r:var Read, reference:RefGetter): bool =
  doAssert m.is_deletion
  var offset = m.start - r.start + 1
  var ref_allele = m.reference[1..m.reference.high]
  var L = r.sequence.len
  if offset < 0:
    return
  r.sequence = r.sequence[0..<offset] & r.sequence[offset + ref_allele.len  .. r.sequence.high]
  var dropped = ref_allele.len
  var seqadd = reference.get(&"{r.chrom}:{m.start+L-dropped+1}-{m.start+L}")
  r.sequence &= seqadd
  doAssert r.sequence.len == L

  var coff = 0

  var cigarCopy = newSeqOfCap[CigarElement](r.cigar.len)
  var added = false
  for i, op in r.cigar:
    if op.consumes.query:
      coff += op.len
    if coff < offset or added:
      cigarCopy.add(op)
      continue

    if not op.consumes.reference: continue
    added = true

    var before = coff - op.len + offset
    if before > 0:
      cigarCopy.add(($before & $op.op).toCigar)
    cigarCopy.add(($ref_allele.len & 'D').toCigar)
    var left = op.len - before
    cigarCopy.add(($left & $op.op).toCigar)

  r.cigar = cigarCopy

proc mutate(m:Mutation, r:var Read, reference:RefGetter): bool =
  ## mutate the read, adjusting the sequence and the cigar and the start if needed.
  doAssert m.chrom == r.chrom
  if m.is_snp:
    return m.mutate_snp(r)
  if m.is_insertion:
    return m.mutate_insertion(r, reference)
  if m.is_deletion:
    return m.mutate_deletion(r, reference)
  quit "wtf:" & $r

proc mutate_alignment(m:Mutation, aln:Record, reference:RefGetter): string =
  var iread = Read(chrom: $aln.chrom, start: aln.start, cigar: toSeq(aln.cigar))
  discard aln.sequence(iread.sequence)
  discard m.mutate(iread, reference)

  var t = aln.tostring().split('\t')
  t[5] = $iread.cigar
  t[9] = iread.sequence
  return t.join("\t")

proc readMutations(ivcf:VCF): TableRef[string, Lapper[Mutation]] =
  var ivs = newTable[string, seq[Mutation]]()
  for v in ivcf:
    if $v.CHROM notin ivs:
      ivs[$v.CHROM] = newSeqOfCap[Mutation](1000)
    ivs[$v.CHROM].add(Mutation(chrom: $v.CHROM, start: v.start, reference: $v.REF, alt: v.ALT[0]))
  result = newTable[string, Lapper[Mutation]]()
  for k, v in ivs.mpairs:
    result[k] = lapify(v)

proc main*() =
  var p = newParser("bamject"):
    help("bamject: inject variants into a bam file and write sam to stdout")
    option("-v", "--vcf", help="required VCF containing sites to mutate")
    option("-b", "--bam", help="required bam file to mutate")
    option("-f", "--fasta", help="required reference fasta to which the bam is aligned")
    option("-a", "--af", help="allele frequency of each injected variant", default="0.05")

  let opts = p.parse()
  if opts.fasta == "":
    echo p.help
    quit "fasta is required"
  if opts.vcf == "":
    echo p.help
    quit "vcf is required"
  if opts.bam == "":
    echo p.help
    quit "bam is required"

  var
    ibam:Bam
    fai:Fai
    ivcf:VCF
  var af = parseFloat(opts.af)
  if not open(ibam, opts.bam, threads=3, index=true, fai=opts.fasta):
    quit "couldn't open bam file"
  if not open(ivcf, opts.vcf, threads=2):
    quit "couldn't open vcf file"
  if not open(fai, opts.fasta):
    quit "couldn't open fasta file"

  var mutationsByChrom = readMutations(ivcf)
  ivcf.close()
  randomize(40)

  echo $ibam.hdr

  var res = newSeqOfCap[Mutation](2)
  for aln in ibam:
    if ($aln.chrom in mutationsByChrom) and mutationsByChrom[$aln.chrom].find(aln.start, aln.stop, res) and random(1'f) < af:
      echo res[0].mutate_alignment(aln, fai)
      continue
    echo aln.tostring()


when isMainModule:
  import unittest

  main()

#[
proc test() =

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
]#
