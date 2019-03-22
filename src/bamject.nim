# vim: ts=2 sw=2 sts=2 sr et ai si
import hts
import lapper
import random
import strformat
import algorithm
import sequtils
import argparse
import ./cigar
export cigar
import strutils
import tables

type RefGetter* = concept i
  ## An object/tuple must implement these 2 methods to use this module
  get(i, string) is string

type Mutation* = object
  chrom*: string
  start*: int
  alt*: string
  reference*: string

# satisfy nim-lapper interface.
proc start*(m:Mutation): int {.inline.} =
  return m.start
proc stop*(m:Mutation): int {.inline.} =
  return m.start + m.reference.len

type Read* = object
  chrom*: string
  start*: int
  sequence*: string
  cigar*: seq[CigarElement]


proc is_snp(m:Mutation): bool {.inline.} =
  return m.alt.len == 1 and m.reference.len == 1

proc is_insertion(m:Mutation): bool {.inline.} =
  return m.alt.len > m.reference.len

proc is_deletion(m:Mutation): bool {.inline.} =
  return m.reference.len > m.alt.len

proc mutate_snp*(m:Mutation, r:var Read): bool =
  doAssert m.is_snp
  var offset = m.start - r.start
  var coff = 0
  for op in r.cigar:
    #echo "op:", $op, " consumes reference:", op.consumes.reference, " consumes query:", op.consumes.query
    if op.consumes.query:
      coff += op.len.int
    else:
      offset -= op.len.int
    if coff >= offset: break

  #echo offset

  #if offset == r.sequence.len: return false
  r.sequence[offset] = m.alt[0]
  return true

proc mutate_insertion*(m:Mutation, r:var Read, reference:RefGetter): bool =
  doAssert m.is_insertion
  var offset = m.start - r.start + 1
  var alt = m.alt[1..m.alt.high]
  var L = r.sequence.len
  if offset < 0:
    alt = alt[abs(offset)..alt.high]
    offset = 0
  var coff = 0
  for op in r.cigar:
    if not op.consumes.query:
      offset -= op.len
    coff += op.len
    if coff >= offset: break
  if offset == r.sequence.len: return false

  var sav = r.sequence
  r.sequence = r.sequence[0..<offset] & alt
  if r.sequence.len < L:
    r.sequence &= sav[offset..sav.high]
  r.sequence = r.sequence[0..<L]
  doAssert r.sequence.len == L
  coff = 0

  var cigarCopy = newSeqOfCap[CigarElement](r.cigar.len)
  var added = false

  for i, op in r.cigar:
    when defined(testing):
      echo "op:", op, " cigar:", $cigarCopy
    if op.consumes.query:
      coff += op.len
    if (op.consumes.reference and coff < offset) or added:
      if i == r.cigar.high:
        # truncate last cigar
        var n = L - (coff - op.len)
        #echo "n:", n, "L:", L, "op:", op, " coff:", coff, " offset:", coff
        #echo "cigar:", $cigarCopy
        #echo "nup:", n
        cigarCopy.add(($n & $op.op).toCigar)
      else:
        cigarCopy.add(op)
      if coff >= L:
        break
      continue
    if not (op.consumes.query and op.consumes.reference):
      ## hack to adjust for softclips before the event.
      if not added:
        coff -= op.len
        L -= op.len
      cigarCopy.add(op)
      continue
    added = true
    # added above to allow continue. subtract op.len here for
    # more intuitive maths.
    coff -= op.len
    var before = offset - coff
    #echo "coff:", coff, " before:", before
    if before > 0:
      cigarCopy.add(($before & $op.op).toCigar)
      if op.consumes.query:
        coff += before
    # now add the event.
    # TODO: handle event at end which must be truncated.
    if coff + alt.len <= L:
      cigarCopy.add(($alt.len & 'I').toCigar)
    else:
      cigarCopy.add(($(L - coff) & 'I').toCigar)
    coff += alt.len
    if coff >= L:
      break

    if before + alt.len < op.len:
      var rem = op.len - before #- alt.len
      when defined(testing):
        echo "coff:", coff, " rem:", rem, " alt.len:", " before:", before, " current:", $cigarCopy
      if rem + coff > L:
        rem = L - coff
      if rem + coff - alt.len > op.len:
        #echo "adj"
        rem = op.len - coff
      if i == r.cigar.high:
        rem = L - coff
      if rem <= 0:
        continue
      cigarCopy.add(($rem & $op.op).toCigar)
      coff += rem
    else:
      var rem = op.len - before
      if rem + coff > L:
        rem = L - coff

      cigarCopy.add(($rem & $op.op).toCigar)
      coff += rem
    added = true
    #echo "aa", $r.cigar

    if coff >= L:
      break
    #echo $cigarCopy

  r.cigar = cigarCopy

proc mutate_deletion*(m:Mutation, r:var Read, reference:RefGetter): bool =
  doAssert m.is_deletion
  var offset = m.start - r.start + 1
  var ref_allele = m.reference[1..m.reference.high]
  var L = r.sequence.len
  if offset < 0:
    return
  #echo "offset:", offset
  var coff = 0
  for op in r.cigar:
    if op.consumes.reference and not op.consumes.query:
      #echo "adj:", op
      offset -= op.len
    coff += op.len
    if offset > coff: break
  if offset < 0 or offset > r.sequence.len: return
    #elif op.consumes.query and not op.consumes.reference:
    #  offset += op.len

  when defined(testing):
    echo "offset:", offset
  var dropped = ref_allele.len
  if offset + ref_allele.len < r.sequence.len:
    #echo r.sequence[0..<offset] #$& r.sequence[offset + ref_allele.len  .. r.sequence.high]
    r.sequence = r.sequence[0..<offset] & r.sequence[offset + ref_allele.len  .. r.sequence.high]
  else:
    dropped = L - offset
    r.sequence = r.sequence[0..<offset] #& r.sequence[offset + ref_allele.len  .. r.sequence.high]

  var seqadd = reference.get(&"{r.chrom}:{m.start+L-dropped+1}-{m.start+L}")
  r.sequence &= seqadd
  doAssert r.sequence.len == L

  coff = 0

  var cigarCopy = newSeqOfCap[CigarElement](r.cigar.len)
  var added = false
  for i, op in r.cigar:
    when defined(testing):
      echo "op:", op, " coff:", coff
    if op.consumes.query:
      coff += op.len
    if (op.consumes.reference and coff < offset) or added:
      cigarCopy.add(op)
      continue

    if not op.consumes.reference:
      when defined(testing):
        echo "adjusting by:", op.len
      coff -= op.len
      L -= op.len
      cigarCopy.add(op)
      continue
    added = true

    coff -= op.len
    var before = offset - coff

    if before > 0:
      cigarCopy.add(($before & $op.op).toCigar)
      coff += before
    cigarCopy.add(($ref_allele.len & 'D').toCigar)
    coff += ref_allele.len
    var left = op.len - before
    when defined(testing):
      echo "before:", before, " coff:", coff, " offset:", offset, " op.len:", op.len, " left:", left
      echo $cigarCopy
    cigarCopy.add(($left & $op.op).toCigar)
    coff += left

  r.cigar = cigarCopy

proc mutate_general*(m:Mutation, r:var Read, reference:RefGetter): bool =
  ## we view a mutation as a loss of the reference and a gain of the alternate
  # so we first remove the reference and then we add the alternate.

  var
    query_offset = 0
    ref_offset = r.start
  var newCigar = newSeqOfCap[CigarElement](r.cigar.len+2)
  var L = r.sequence.len
  var added = false

  for op in r.cigar:
    if op.consumes.query:
      query_offset += op.len
    if op.consumes.reference:
      ref_offset += op.len

    if ref_offset < m.start or added:
      newCigar.add(op)
      continue

    var overshot = ref_offset - m.start #- delta

    # remove the allele from the reference
    r.sequence = r.sequence[0..<(query_offset-overshot)] & m.alt & r.sequence[(query_offset-overshot+m.reference.len)..r.sequence.high]
    added = true

    # don't need to update cigar string for snp
    if m.is_snp:
      newCigar.add(op)
      continue

    var back = op.len - overshot
    # add 1 because the alt is e.g. 3 bases to indicate a 2-base deletion
    var before = op.len - overshot + 1

    newCigar.add(($before & $op.op).toCigar)
    if r.sequence.len > L:
      r.sequence = r.sequence[0..<L]
      newCigar.add(($(L - query_offset + op.len - before) & 'I').toCigar)

  r.cigar = newCigar


proc mutate*(m:Mutation, r:var Read, reference:RefGetter): bool =
  ## mutate the read, adjusting the sequence and the cigar and the start if needed.
  doAssert m.chrom == r.chrom

  discard m.mutate_general(r, reference)
  return

  #[
  if m.is_snp:
    return m.mutate_snp(r)
  if m.is_insertion:
    return m.mutate_insertion(r, reference)
  if m.is_deletion:
    return m.mutate_deletion(r, reference)
  quit "wtf:" & $r
  ]#

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
    if len(v.ALT) == 0:
      stderr.write_line "[bamject] skipping variant with empty alternate allele:" & v.tostring()[0..<60]
      continue
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
    obam:Bam
  var af = parseFloat(opts.af)
  if not open(ibam, opts.bam, threads=3, index=true, fai=opts.fasta):
    quit "couldn't open bam file"
  if not open(ivcf, opts.vcf, threads=2):
    quit "couldn't open vcf file"
  if not open(fai, opts.fasta):
    quit "couldn't open fasta file"
  if not open(obam, "bamject.bam", mode="wb"):
    quit "couldn't open output bam file"

  var mutationsByChrom = readMutations(ivcf)
  ivcf.close()

  #randomize(40)
  obam.write_header(ibam.hdr)

  var res = newSeqOfCap[Mutation](2)
  for aln in ibam:
    if ($aln.chrom in mutationsByChrom) and mutationsByChrom[$aln.chrom].find(aln.start, aln.stop, res) and random(1'f) < af:
      var before = aln.tostring
      var after: string
      try:
        for m in res:
          if not m.is_snp: continue
          after = m.mutate_alignment(aln, fai)
          aln.from_string(after)
      except:
        stderr.write_line "mutation:" & $res
        stderr.write_line "before:" & before
        stderr.write_line "after:" & after
        raise

    obam.write(aln)
  obam.close()
  ibam.close()

when isMainModule:
  main()
