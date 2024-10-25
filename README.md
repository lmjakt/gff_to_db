# gff files to gene databases

This repository contains a small C++ program that extracts information about
genes, transcripts, exons and coding sequences. It's designed to:

1. Be correct; i.e. to sanity check parent-child relationships.
2. To make as few assumptions as possible about the encoding of the gene
   information.
2. To convert exon to gene mapping to exon to transcript to gene maps.
   This is done in order to create a uniform structure that can more
   easily be queried.
3. To remove annotation redundancy information. Currently annotation for protein
   products in gff files are often bound to all exons and all CDS regions
   associated with a transcript or a gene. The program tries to collapse this
   to individual transcripts.
4. To remove repeated exon and CDS information. That is, an exon should only be
   defined once, regardless of the number of transcripts it belongs to (and the
   same is true for CDS regions).
5. To output a set of files for each gff that use numeric identifiers to map the
   relationships between genes, transcripts, exons and CDS regions, in a similar
   way to how Ensembl represents gene information.
6. To hopefully be reasonably fast. Right now it seems to take about 20 seconds to
   process a file with about 1 million features. So processing a thousand such files
   should take around 6 hours of CPU time. I don't know if this is fast or slow,
   but I should compare it to other tools.
   
The current version (2024-10-25) is experimental. That means that it outputs redundant
information in the output files in order to more easily check that the output
is correct. This behaviour should be optional.

## Output

The program outputs eight tab separated files:

1. `region` Names and lengths (if known) of sequences mentioned in the gff file.
2. `gene` Names and locations of genes.
3. `rna` Names and locations of transcripts. Includes a gene parent identifier.
4. `tr_annotation` Annotation of transcripts; collapsed from the attributes fields
   of exons, CDS and transcript features. One transcript can contain several entries
   of several fields. At this time the `name`, `product` and `description` fields are
   extracted. This table should have a companion one for different types of attributes
   and potentially a bitwise flag defining the origin of the attribute (i.e., CDS,
   exon or transcript).
5. `exon` Identifiers (numeric) and locations of exons features.
6. `cds` Identifiers (numeric) and locations of CDS features.
7. `tr_exon` Map of transcript identifiers to exon identifiers (both numeric).
8. `tr_cds` Map of transcript identifiers to CDS identifiers as for exons.

There should also be a file giving attribute type identifiers in order to avoid repeating
the names of types in the tr_annotation file. One can also consider to collapse the
annotation fields as these are likely to have a fair amount of redundancy.

## Limitations

### Incomplete data extraction

Only features identifed (after conversion to lower-case) as: 

1. `cds`
2. `exon`
3. `mrna`
4. `lnc_rna`
5. `transcript`
6. `rna`
7. `gene`
8. `region`

are extracted. This is reasonble for _my_ current requirements, but it means that
other types of information is not obtained.

### Correctness

I have not yet tested if the program behaves correctly. This means it currently outputs
more data fields than necessary.

### Exons lacking for CDS fields (orphan CDS regions)

Some annotation pipelines only output CDS regions without mapping these to exons. In a sense
this is probably technically correct if the annotation is not supported by RNA sequencing;
however, it makes automating the analyses of large numbers of genomes inconvenient. It is
also inconsistent, since in many cases it seems clear that exons are simply defined from
CDS regions. Unfortunately, I have not come up with a method for defining exons for orphan
CDS regions; defining a rule that works for almost all CDS-regions is fairly trivial but
there are some gff files that I have come across that have very strange definitions of exons,
and I have not yet decided how to deal with these. Ideally, they should just be flagged
as `bad`, and ignored.

