import hts/bam
import strutils

const BAM_CIGAR_SHIFT = 4'u32
const BAM_CIGAR_STR = "MIDNSHP=XB"

var cigar_tab : array[128, int]
for i in 0..<128:
    cigar_tab[i] = -1
for i, c in BAM_CIGAR_STR:
    cigar_tab[c.int] = i

type Cigar* = seq[CigarElement]

proc `$`*(c:Cigar): string =
  result = ""
  for e in c:
    result &= $e

proc tocigar*(cs:bam.Cigar): Cigar {.inline.} =
  result = newSeqOfCap[CigarElement](cs.len)
  for e in cs:
    result.add(e)

proc tocigar*(cs:string): Cigar {.inline.} =
  if cs.len > 4:
    result = newSeqOfCap[CigarElement](2)
  var off = 0
  while off < cs.len:
    var i = 0
    while cs[off + i].isdigit:
      i += 1
    var num = parseInt(cs[off..<off+i])
    var ops = cs[off+i]
    off += i + 1
    var el:uint32 = num.uint32 shl BAM_CIGAR_SHIFT
    if cigar_tab[ops.int] == -1:
      quit "unknown cigar op from " & cs & ": " & $ops
    el = el or cigar_tab[ops.int].uint32
    result.add(cast[CigarElement](el))
