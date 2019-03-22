###  bamject: inject variants into a bam file

```
Usage:
  bamject [options]

Options:
  -v, --vcf=VCF              required VCF containing sites to mutate
  -b, --bam=BAM              required bam file to mutate
  -f, --fasta=FASTA          required reference fasta to which the bam is aligned
  -a, --af=AF                allele frequency of each injected variant (default: 0.05)
  -h, --help                 Show this help

E.g.:
```

```
bamject -b $bam -f $fasta -v mut.vcf --af 0.05 | samtools view -o mutated.bam
```
