(C) 2008 to 2014 NovoCraft Technologies Sdn Bhd

Freely downloaded versions of Novocraft programs are licensed for educational use
and for use in not-for-profit organisation for internal use only. Some features 
are disabled in the free versions.

For-profit organisations without paid up licenses can only use this software for
the purpose of evaluating it.

For commercial licenses and support contracts contact sales@novocraft.com

Release Novoalign V3.02.07 & NovoalignCS V1.05.00, Novosort (V1.03.01), Novomethyl (V1.01)
------------------------------------------------------------------------------------------
Novoalign
      1. Earlier versions of Novoalign had a 15bp indel limitation on the read that anchors a pair 
         alignment. This limitation was causing problems with long reads that have large overlaps 
         such that both reads of the pair have indels longer than 15bp. This change removes the 
         15bp indel limit.

Release  Novosort (V1.03.01)
----------------------------
      1. Fix: Novosort --markduplicates -i creates an invalid BAM index.
      2. Fix: Seg Fault at termination if marking duplicates on presorted BAM files.

Release Novoalign V3.02.06 & NovoalignCS V1.05.00, Novosort (V1.03.00), Novomethyl (V1.01)
---------------------------------------------------------------------
Novoalign
      1. Added support for single line tabbed separated read files. (-F TSV)
      2. Added new option for processing reads from a BAM file. -F BAM will 
         align both single end and paired end reads in the BAM input.

Novoalign & CS
      1. The log file has been changed to include % of reads and to add new 
         classifications. It should be self explanatory but may break any scripts 
         that process the log file.
      2. Added a --tags option to control the SAM tag attributes. 
      3. If single end alignment quality for a read is greater than paired end alignment 
         quality the SAM QUAL attribute is set to the single end quality. However, in 
         some cases where a read is very short (e.g. due to low quality bases) or aligns 
         to highly repetitive sequence the single end quality may not be accurate. This 
         change stops the QUAL attribute being set to the single end quality in these cases.
      4. Soft Clipping was adjusted so that bases with quality of 2 or bases aligned against 
         a reference N get a score of zero as per Novoalign V2. This allows a read alignment
         to extend into a region of N's if the --softclip option is set to a non-zero value.
      5. Fix: novoalign: src/findamate.cc:18: Assertion `tgtlen <= ncols - 32' failed. This
         was most likely to have happen with circular genomes.
NovoalignCS
      1. Added support for XSQ format read files in NovoalignCS and NovoalignCSMPI. 
         This will require HDF5 to be installed if using MAC OSX or if building
         your own NovoalignCSMPI.
           1. Download the source code of HDF5 1.8.8 from HDF5 web site
           2. Open terminal and extract the tar 
           3. ./configure –prefix=/usr/local –enable-cxx 
           4. make
           5. make check
           6. sudo make install
Novosort
      1. Added an option to remove or mark duplicate reads. More details are in the 
         Novosort manual and white paper.

Novomethyl
      1. Fix: In diploid consensus report the counts of converted & total cytosines 
         are wrong for heterozygous bases.

Release Novoalign V3.02.05 & NovoalignCS V1.04.06, Novosort (V1.02.02)
---------------------------------------------------------------------
Novoalign
      1. Fix: Assertion `tgtlen <= ncols - 32' failed.

Release Novoalign V3.02.04 & NovoalignCS V1.04.05, Novosort (V1.02.02)
---------------------------------------------------------------------
Novoalign
      1. Fix: A Signal 8 (Floating Point Divide by Zero) could occur when a fixed 
         threshold was used with hard clipping option on paired end reads. (e.g. -H -t 320)
      2. Fix: In V3.02.03 paired end mode the adapter trimming option -a was 
         mistakenly changed so that if only one adapter was specified the second adapter
         would default. In all previous versions the second adapter would be set equal 
         to the first adapter, this has been restored.

Release Novoalign V3.02.03 & NovoalignCS V1.04.05, Novosort (V1.02.02)
---------------------------------------------------------------------
Novoalign & CS
      1. Fix: When reporting multiple alignments for paired end reads it was possible
         that the same alignment reported as secondary with multiple alternate mate alignments
         had different values for the SM tag (Single End alignment quality)
      2. Runtime has been reduce slightly as a result of performance tuning.
NovoalignMPI & CSMPI
      1. This release allows you to compile and link your own MPI programs. Necessary 
         files are in folder buildMPI together with instructions. This has been tested 
         on several Linux releases with MPICH2, MPICH and Cray MPICH. It has not been 
         tested with other MPIs such as OpenMPI or IntelMPI.

Release Novoalign V3.02.02 & NovoalignCS V1.04.04, Novosort (V1.02.02)
---------------------------------------------------------------------
Novoalign & CS
      1. Fix: CPU time of long running jobs may get reported as a negative number.
      2. Fix: Single end reads that are satellite repeats may report duplicate alignment 
         locations when using -r All.
      3. The -R option increases the score range for printing multiple alignments to 
         a read. This change increases sensitivity for high scoring alignments when 
         -R setting is >60.
      4. Allow alignment score threshold to be set with -t A,B format option when 
         using -r Exhaustive, thus allowing threshold to be a function of read length.
         Previously -r Exhaustive required a specific threshold for all reads.
      5. Allow use of -R option with -r Exhaustive. The -t setting is used to determine 
         the set of seeds for the alignment process and then all alignments with alignment 
         score range [0,'-t' + '-R'] setting are reported.
         Sensitivity for alignments with score in the range ('-t', '-t' + '-R'] is
         reduced compared to alignments with score [0,'-t'].
         In earlier releases the -R option had no effect on -r Exhaustive alignments.
         e.g. -r Exhaustive 1000 -t 150 -R 60 will report alignments with a score less than 210 
         but sensitivity is reduced for alignments with a score >150.
         This change is designed to work with LSC (http://www.stanford.edu/~kinfai/LSC) 
         when using Illumina reads to correct Pacbio reads. 
      6. The actual alignment score threshold used by Novoalign may be less than the 
         threshold set with the -t option as Novoalign calculates an upper limit for 
         threshold based on read length and base qualities. This upper limit is applied 
         to stop nonsense alignments like deleting the whole read or aligning the whole 
         read to a stretch of N's. This change enforces the threshold upper limit for 
         -r Exhaustive alignments.
      7. Fix: Amplicon primer trimming for amplicon alignments from -ve strand were being counted 
         as individual primer matches rather than as full primer pair match.
      8. Upgraded memory allocator to tcmalloc 2.1
      9. Reduced memory usuage during long runs by about 0.5Gbytes.
     10. Fix: A rare buffer underrun could cause various errors such as Seg Fault 
         or Checksum errors in malloc.
     11. In MD tag of SAM alignments a multibase insert was shown like 20^A0^T0^C10. While 
         this is valid, the 0^ between inserted bases is not required and has been removed 
         so the example insert will now be shown as 20^ATC10
     12. Fix: Novoalign in mate pair mode with long & short insert sizes may fail an assert
         if maximum fragment length allowed by the insert size settings exceeds 32767bp. The fix 
         increases the limit to 65000bp and will report an error if the insert size is set higher 
         than this.
         
NovoalignCS
      1. Increased alignment score range in the first alignment stage. This will identify 
         more alignments for a read and increase the accuracy of the alignment quality 
         score.

Novoalign[CS]MPI
      1. Fix: Some run stats in the log file were only showing values for one slave process. 
         Most noticeably the CPU time used.
      2. Fix: NovoalignMPI crashes with an assert failure in quality calibration if the 
         first 2000 reads all have alignment quality < 70
         
Novosort
      1. Fix: If the first alignment in a sort buffer was out of order and remaining 
         records were all in order then the buffer would not be sorted. This problem was 
         observered when sorting BAM files with two alignments.
      2. Fix: Novosort may seg fault when adding RG tags or merging multiple bam files 
         due to a buffer overrun.
      
Release Novoalign V3.02.00 & NovoalignCS V1.04.02, Novosort (V1.02.01)
---------------------------------------------------------------------
Novoalign/MPI
      1. Fix: In V3.01.02, Paired End adapter trimming (-a option) involves alignment
         of the two reads of a pair against each other and the adapter sequences. In 
         this alignment process, the gap penalty was excessively high relative to 
         the mismatch penalties effectively disallowing gapped alignments.
      2. Fix: In Paired End adapter trimming (V3 Only), if the adapter trimming 
         alignment was gapped and -r None was used (default) then it was possible 
         the read was reported as a multi mapper when in fact it had a unique 
         alignment location. This problem was evident in V3 releases prior to V3.01.02 
         but masked in V3.01.02 due to the high gap penalty used in adapter trimming. 
         This only affects reads where DNA template molecule was shorter than the read 
         length.
      3. The adapter trimming option -a now allows degenerate bases (N) in the adapter
         sequence.
      4. Add option for soft clipping amplicon primers from alignments,  
          --amplicons amplicons.BED
      
      Bed File Format
          chrom        Name of the chromosome 
          chromStart   Start position of the amplicon (includes primer bases)
          chromEnd     End position of the amplicon
          name         Amplicon name if any
          score        ignored
          strand       + or -, ignored for now.
          thickStart   Start of amplicon excluding primer
          thickEnd     end of amplicon excluding primer
          itemRgb      ignored
          
      Example
          chr2	29083861	29084059	AMP.1	100	-	29083881	29084039
          chr2	29085075	29085273	AMP.2	100	-	29085095	29085254
          chr2	29089969	29090233	AMP.3	100	-	29089989	29090214
          chr2	29091056	29091241	AMP.4	100	-	29091076	29091220

      At end of run the counts of amplicon clipping events is printed to stderr 
      e.g.
        #	Amplicon	Count	SE5	SE3
        #	AMP.1        371	  0	  0
        #	AMP.2        190      0   0
 
      There are three counters, first is the number of hits where both reads of 
      pair aligned to primers of same amplicons. Next two counts are where read1 
      & read2 of pair aligned to different amplicons or perhaps one read of the 
      pair failed to align.
      
Novoalign/MPI/CS/CSMPI
    1. Fix: In SAM format proper pair alignments, the single end quality tags SM: and
       AM: could report incorrect quality scores. This fix improves the accuracy 
       of the single end quality scores (SM: tag). The change also affects some 
       proper pair alignment qualities by a small amount (typically +-1).
       
NovoalignCS/CSMPI
    1. Allow Colour Space Reads without the Primer nucleotide. In colour space
       reads we normally see the last base of the primer in the read sequence 
       followed by colours. This gives a clue to the translation of the first colour. 
       In some XSQ to CSFASTQ translations the primer base is missing meaning that 
       the read could start with any base. If NovoalignCS finds the primer base is 
       missing it assumes a first base of 'G' and reduces the quality of the first 
       colour to 2 hence removing any penalty for an initial colour error.
       Note. It is better to use a conversion tool that includes the first base and 
       colour in the csfastq file as it maximises specificity of alignments.
       
Novosort
    1. Fix: V1.02.00 of Novosort would seg fault if an input BAM file had a missing 
       EOF block or other error.
    2. Fix: If novosort program was renamed then it would not find the novosort.lic file.


Release Novoalign V3.01.02 & NovoalignCS V1.04.02, Novosort (V1.02.00)
---------------------------------------------------------------------
Novoalign*
      1. Honour --nonc option in NovoalignMPI & NovoalignCSMPI. This results in 
         slightly faster alignments as there is no need to synchronise fragment 
         length and quality calibration tables across slaves.
      2. Fix: If novoalign V3 was run with a low threshold (-t option) it was possible for 
         it to get in an infinite loop. 
      3. Add an option to only process a fixed % of reads. 
             Format -# 9.9%
         Example with -# 1% Novoalign will align every 100th read, skipping intervening
         reads.
      4. Fix: When aligning with a low -t setting some reads with a lot of low quality 
         bases were being aligned even though their best possible alignment score was 
         greater than the -t setting. 
      5. Fix: When aligning RNA reads with a reference that contains exon junction 
         sequences and with three penalty format of -v option then memory utilsation 
         may grow excessively 
      6. Fix: Paired end adapter trimming required higher identity for long reads than 
         for shorter reads and may not have been sensitive enough for 250bp MiSeq reads. 
         The fix requires approx >90% identity for all read lengths. The change reduced 
         mismatch penalties and match rewards in the adapter trimming function..
         
Release Novoalign V3.01.01 & NovoalignCS V1.04.01, Novosort (V1.02.00)
---------------------------------------------------------------------
Novoalign*
      1. Fix. Novoalign reports errors and stops for fastq files with zero length reads & 
         qualities.
         These can occur as result of external adapter trimming programs and were acceptable 
         in earlier version of Novoalign.       
      2. Fix. Novoalign was nonconcordant (ie. rerun produces different results) if 
         quality calibration option (-k) was used.
      3. Fix. If reference ambiguous nucleotide of BDH or V is aligned to a base with quality of < 9
         then an incorrect and excessively high mismatch penalty is applied. This would usually 
         cause the alignment to exceed the score threshold. Ambiguous codes for 2 nucleotides have no issues.
      
Novosort (V1.2.00)
---------------------------------------------------------------------
         Memory allocator was changed to Google's tcmalloc. This reduces memory overhead for 
         threads and allows Novosort to be run with less memory. Minimum memory (-m) with 16 threads
         is 0.16 * SQRT(B) where B is the size of input BAM file in Gbytes. 
         e.g For 100Gbyte unsorted compressed BAM the -m should be >= 0.16 * 100^0.5 = 1.6G
         Using less memory than this creates a large number of temporary sorted segments 
         and then a large number of threads to merge the temporary segments. Additional 
         per thread memory can cause actual memory used to exceed the -m value.
         Using more memory than above is perfectly okay. Best performance is best when the 
         uncompressed BAM file fits in -m memory.
         Actual memory used will be about 500M higher than the -m setting.

Release Novoalign V3.01.00 & NovoalignCS V1.04.00, Novosort (V1.1.00)
---------------------------------------------------------------------
Novoalign
      1. Fix: Version 2 single end adapter trimming could trim as little as one base 
         from a read. Version 3 was trimming a minimum of 2bp matching adapter sequence. 
         This change restores V2 behaviour.
      2. Softclipping algorithm has been changed to allow for a score bonus for 
         alignments extending to the end of a read. This gives some control over soft 
         clipping versus earlier versions. The option to enable softclipping with a 
         end of read bonus is --softclip <99> where <99> is the bonus, typically setting 
         is 40. The current -o Softclip option is equivalent to --softclip 0 
         and the -o FullNW option to --softclip 9999. Old format options are still available.
         Softclipping typically improves specificity of downstream SNP calling and 
         large indel detection with programs like PINDEL but it can reduce SNP sensitivity 
         in regions with low read coverage and is most noticeable with protocols that 
         produce variable read depth such as sequence capture projects. Use of a value 
         greater than mismatch penalty (30) and less than a single base indel (46) penalty 
         is suggested.
      3. Fix: When aligning mate pairs with both long & short fragment lengths Novoalign
         could crash with an assert failure or a seg fault.
Novoalign & CS
      1. Fix: When setting alignment threshold to a fixed value (e.g. -t 200) then 
         reads with very low quality bases and alignment scores greater than the -t 
         setting could get reported.
         
Release Novoalign V3.00.05 & NovoalignCS V1.03.05, Novosort (V1.0.05)
---------------------------------------------------------------------
Novoalign
      1. Usually paired end alignment quality is higher than the single end alignment
         quality of each read, and in paired end mode earlier versions of Novoalign reported  
         the alignment quality of the pair for each read. However, there are cases where  
         the quality of an alignment is higher for an individual read than for the pair. 
         eg. One read of pair aligns uniquely and the other has multiple nearby mappings 
         due to a tandem repeat.
         We now report the maximum of paired or single end quality for each read.
         This change has also improved the accuracy of single end quality score (SM & 
         AM tag) for paired end reads.
NovoalignCS
      1. Fix(V3): If a read alignment overlapped the ends of a sequence it could be reported 
         as not aligning or as an insert. This change will report more alignments with 
         inserts at the end of the reference seqeunces. This mainly affects circular genomes.
      2. When generating base qualities from a colour space alignment cap the base quality at 60.
NovoalignMPI
      1. Fix: NovoalignMPI V3 would go into an infinite loop at EOF when processing reads 
         from a BAM file.
      2. Fix: NovoalignMPI V3 incorrectly calculated fragment length penalties when aligning 
         mate pairs reads with two fragment length ranges specified.  
      3. Fix: NovoalignMPI V3 fragment length counts were inconsistent with Novoalign. 
      4. Fix(V3): NovoalignMPI results were not identical to Novoalign results.     
         
Release Novoalign V3.00.04 & NovoalignCS V1.03.04, Novosort (V1.0.04)
---------------------------------------------------------------------
Novoalign
      1. Fix: Version 3 would not align reads with large numbers of low quality (#) 
         bases due to calculation of alignment threshold. The % of low quality bases 
         depended on the gap extend penalty.
      2. Fix: Version 3 could give low alignment quality to reads that had many low 
         quality bases even if read was uniquely aligned.
      3. Fix: If a read aligned at the end of sequence and include extra bases beyond
         the end of the sequence it was possible that it would be reported as an insert
         at the beginning of the next reference sequence.
NovoalignCS
      1. Change to how _F3 _F5 _R3 _R5 are trimmed from read names. In previous 
         versions read name was first trimmed to make the name unique and then any
         trailing underscore was removed. This could leave an _F on the end of the 
         read names. We now trim until unique and then remove trailing F, R and/or 
         underscore.
      2. Fix: If a read aligned at the end of sequence and include extra bases beyond
         the end of the sequence it was possible that it would be reported as an insert
         at the beginning of the next reference sequence.
Novoutil
      1. Novoutil IUPAC normally takes a vcf SNP file and a fasta sequence
         file and adds the SNPs to the fasta file as IUPAC ambiguous codes. 
         It can optionally add homozygous indels to the fasta file which is 
         useful for correcting assemblies. This change adds the ability to 
         remap bed file coordinates based on the coordinate changes from the 
         homozygous indels so they match the new fasta file.
          
Release Novoalign V3.00.03 & NovoalignCS V1.03.03, Novosort (V1.0.03)
---------------------------------------------------------------------
Novoalign & CS
      1. Fix: If -i option set a mean fragment length of zero (e.g. -i 0,100) or 
      minimum length of 0 (e.g. -i 0-500) then the alignment process fails.
NovoalignCS
      1. Fix (V3.00 Bug): If a read aligned close to the start of a chromosome 
         it could be aligned as though the full read was inserted at the end of
         the previous chromosome. 
         
Release Novoalign V3.00.02 & NovoalignCS V1.03.02, Novosort (V1.0.02)
---------------------------------------------------------------------
Novoalign & CS
      1. More informative error messages if there are problems processing the 
         read sequence files. Issues such as different length sequence and quality 
         strings will now report the affected read header.
      2. After processing a fastq read if next line does not start with '@' then 
         we report the problem and proceed to process the next read line. Novoalign 
         will only stop after 10 such errors.
      3. Spaces are now accepted as valid qualities in Colour Space fastq files.
         
Release Novoalign V3.00.00 & NovoalignCS V1.03.00, Novosort (V1.0.02)
---------------------------------------------------------------------
Novoalign
      1. Increase maximum read length by implementing score scaling during SIMD 
         alignment phase. In earlier versions a mismatch at a good quality base
         scored 30 points limiting alignments to 8 mismatches due to score limit
         of 255 in the byte arithmetic of the SIMD routines, and hence limiting 
         ability to align reads with more than 8 mismatches (low quality bases 
         score less than 30 for a mismatch so it was possible to have more than 
         8 mismatches if some were at low quality bases).
         This release allows scoring to be scaled by a factor of 1/2, 1/3,
         to 1/6 in SIMD routines allowing more mismatches and longer indels per read.
         Scaling factor is determined automatically based on the alignment score 
         threshold for the read.
         
         Maximum read length has now been increased to 950bp with a score threshold 
         of 1500, allowing 50 mismatches at high quality bases or indel of 250bp.
         
      2. Remove -Q option
      3. Remove -r [0.99] option
      4. In SAM report format remove custom tag ZN as it had the same value as NH tag
      5. Increased maximum length gap that can be aligned in single end mode from 
         15bp to ~50% of read length with upper limit of 250bp.
      6. In SAM report format remove AS tag as it has same value as UQ tag.
      7. In paired end mode the reported alignment scores (UQ tag) no longer includes
         the fragment length penalty.
      8. Quality calibration is now done on full Needleman-Wunsch alignment rather 
         than using the softclipped Smith-Waterman alignment. This means that mismatches
         in first & last few bases of read are more likely to be included in quality 
         calibration and hence lead to slightly lower quality for these bases. The previous
         version artifically inflated quality of bases near the ends of the reads as 
         mismatches were clipped off the alignment.
      9. Adapter trimming for single end reads now allows indels in the adapter alignment. 
         This is particularly helpful for 454 reads.
      10. Increase in match reward from 6 to 8 when soft clipping alignments.
      11. Default alignment theshold is now -t A,4.5 and results in a threshold of (N-A)*4.5 
          where A =log4(Reference Genome Length) and N is read length.
      12. Corrections for TLEN attribute in SAM format files. Was out by 1 in version 2. Also
          we now use mapping location after soft clipping. 
      13. Reruns of Novoalign with identical data and parameters now produce identical 
          results. Previously multi-threading could lead to slight differences in results 
          between reruns of the same data. NovoalignMPI results should also be identical
          to Novoalign results.
          
NovoalignCS
       1. Code restructuring has reduced run time.
       2. Default alignment theshold is now -t A,4.5 where A =log4(Reference Genome Length). 
          This slightly lower than the previous default and slightly reduces the number of 
          reads aligned. You can get back the extra alignments by setting -t 20,5.
          
Novosort
      1. Default memory reduced to 50% of RAM to allow for more disk cache.
      2. Fix: When sorting alignments made to a reference with a very large number of
         sequences (eg. RNA Exon junction sequences) novosort could be extremely slow.
      3. Name sort order was changed to keep secondary alignment pairs together and order the
         reads closer to the original Novoalign report order. 
         Sort order is now..
         
         For Proper pairs:
         Name, Primary/secondary flag, Hit Index (or Alignment Location if no HI tag), Read1/Read2
         
         For non-Proper pairs:
         Name, Read1/Read2, Primary/secondary flag, Hit Index, Alignment location

Novoutil
      1. Fix: Function extractsv was not working.
         
Novobarcode
      1. Support for Casava V1.8 file format with index read in the header.
         EG.    @EAS139:136:FC706VJ:2:5:1000:12850 1:Y:18:ATCACG

Release Novoalign V2.08.05 & NovoalignCS V1.02.05, Novosort (V1.0.01)
--------------------------------------------------------------------
NovoalignCS
      1. Fix: Some csfastq files have quality codes of space (quality -1), if this
         occurs as the last quality novoalignCS stops and reports quality string is 
         shorter than the colour sequence.
      2. Fix: Some base qualities genertaed for aligned read are being rejected by 
         some downstream tools. Cap generated base qualities at 60.
      
NovoalignMPI
      1. NovoalignMPI would skip ~5% of reads when input was a BAM file

Release Novoalign V2.08.04 & NovoalignCS V1.02.04, Novosort (V1.0.01)
--------------------------------------------------------------------
Novoalign
      1. Fix: When adapter trimming in paired end mode and one read of pair is 
         very low quality causing search to drop back to single end mode then 
         the hard clipping of the adapter may not get reported in the CIGAR field.
      2. Remove ZN: tag from SAM format reports.

Novobarcode
      1. Changed algorithm for paired end classification when tag read is on 5' of both
         reads. The algorithm now aligns both tag reads and combines scores to find 
         tag with lowest total alignment score. Previously each tag read was treated 
         separately and both tag reads had to align to the same tag.

Release Novoalign V2.08.03 & NovoalignCS V1.02.03, Novosort (V1.0.01)
--------------------------------------------------------------------
Novoalign
      1. Fix: Some alignment jobs were running very slowly. In this case it was caused by 
         use of quality calibration with reads contained many bases with a quality of 0. 
      2. The hard clipping option -H now clips N's even if there quality is > 2.
      3. Added a new option for limiting time spent trying to align reads which are 
         primarily homopolymeric. In many cases these reads are artifacts from the 
         sequencing process and they can take a long time to align as they seed to 
         many locations in the genome and when they do align they usually have very 
         low algnment quality.
         Option format --HLimit <F>
         Where F is typically in range 5-15 and limits threshold to F*nH where nH 
         is number of mismatches required to align to a perfect homopolymer. As 
         mismatches typically score 30 a value of 10 would allow 1/3 the number 
         of mismatches as there are bases differing from a homopolymer.
         We suggest using --HLimit 8 fro Bi-sulphite alignments.

Novoalign & CS 
      1. Allow the alignment threshold to be specified as a function of read length in form 
         threshold = (L - A) * B where L is read length (sum of pairs) and A&B are specified
         on the -t option in format  -t A,B
         B can be fractional and should always be <= the gap extend penalty.
         This particularly useful for ION Torrent reads if you want to use a threshold less
         than the default.
      2. Existing option -q <9> allows alignment quality to be reported as a float 
         with specified decimal places. In previous releases this only worked 
         with Native and Pairwise report formats. We now extend this to SAM format
         by adding the tag ZQ:f:9.99 with alignment quality to the specified number
         of decimal places.
      3. When checking that two headers of a pair match (hamming distance) the check is now
         stopped at the first space or tab character, ignoring any characters after this.
      
Novoutil iupac
      1. When updating a reference from a vcf file and the option to apply indels is selected,
         we now apply heterozygous inserts as well as all homozygous indels. This makes 
         the reference heterozygous to a deletion.
      2. Added support for multi-sample vcf files.
         a) Added an option to choose sample -s <samplenumber> with the alt alleles for forming 
            IUPAC code selected based on the GT tag of sample. If GT tag is ./. then SNP is not
            applied (not applicable to this sample, effectively GT of 0/0). Example If GT tag is
            2/2 then only the second alt allele is used to form the iupac code (genotype).
         b) If sample number is not specified then zygosity is determined 
            from the AF attribute rather than any of the GT attributes. AF=1.0 means genotype
            is formed from alt alleles only.
         c) If file is single sample and there is no AF tag, then GT tag of first and only 
            sample determines the genotype.
      3. Fix: In previous versions any lower case bases codes were converted to upper case 
         and counted as modified bases (SNPs). Case is now preserved.

Release Novosort (V1.0.01)
--------------------------
Novosort
      1. Fix: Name order sort failed if sort required merge of temporary work files.

Release Novoalign V2.08.02 & NovoalignCS V1.02.02, Novosort (V1.0)
-------------------------------------------------------------------
Novoalign
      1. Fix: GATK complains if the CIGAR field contains a multibase substitution
         encoded as adjacent insert & delete operations. Multi base substitutions are 
         now encoded in the CIGAR as M's. This does not change the dynamic programming 
         algorithm and scoring. Internally a long substitution may be scored as 
         adjacent inserts and deletes but then shown in the CIGAR as mismatches plus 
         an indel to make up any length difference.
      2. Removed input file buffering when reading from pipes. This allows Novoalign 
         to be used as a service.
      3. Fix: If a paired read sequence contains many invalid IUB NA codes then Novolaign
         performance is reduced.
         Rather than rejecting read sequences with invalid base codes Novoalign treats them 
         as mismatches to all bases. This fix will QC any read with more than 8 invalid 
         base codes.
NovolaignMPI
      1. Fix: if option --hdrhd off was used on NovoalignMPI the reads were processed 
         as single end rather than paired end.
NovoalignCS
      1. Fix: When using option --rOQ to report pre-calibrated colour qualities the 
         qualities on -ve strand alignments were reversed relative to the CS & CQ tags. 
         The fix reports CS, CQ & OQ in the same direction as the original read.
Novoutil bgzf
      1. Fix: Novoutil is completing with odd exit status values. Exit status 
         is now zero unless there was an error.
Novoutil iupac
      1. Fix: If quality attribute was '.' then SNPs were rejected as low quality. We 
         now accept these SNPs
      2. Fix: If genotype couldn't be determined from GT, AC or AC1 attribute SNPs 
         were assumed to be homozygous for the alt allele. We now assume heterozygous 
         for reference & alt allele. Note. This only has an affect when -g option is used.
      3. Fix: SNPs would not match chromosome names if the case was different. We now 
         do a case insensitive compare of names.
Novosort
      1. Fix: In some situations the sort was not "stable" and could change the order 
         of two alignnments with the same alignment location. A stable sort is one
         where if two records have the same sort key their relative positions are 
         not changed. This was likely to happen if the input BAM file size was 
         between 1 & 2 times the amount of RAM allocated for the sort.
      2. Fix: The -a option was not being recognised.
      3. Fix: If an @RG record was given as an option and the input BAM file already 
         had an @RG record then the replacement RG tag was only substituted on alignments 
         that had an RG:Z: tag. Alignments without RG tags stayed that way.
      4. Added option for name sort
      5. Added option to index the sorted bam file

Release Novoalign V2.08.01 & NovoalignCS V1.02.01, Novosort (beta3)
-------------------------------------------------------------------
Novoalign
      1. *** Novoalign now attempts to place indels in the most 5' position of reference
         given alternative alignments with the same score. Previously indels were 
         placed in the most 3' position.
      2. Fix: When aligning reads from a BAM file, Novoalign would stop with an
         error message if read2 of a pair was before read1 of the pair. The fix 
         accepts either ordering.
      3. Fix: Novoalign stops with a Seg Fault when -m option is used.
      4. Softclipping was changed so that if bases at the end of a read aligned
         to N's in the reference they will not be soft clipped. This facilitates 
         calling the consensus sequence at these locations.
      5. Fix: If a folder name was given instead of a read file then Novoalign assumed prb
         file format and entered an infinite loop processing zero length records
Novosort
      1. Report an error if creation of a temporary file fails. This may happen if the 
         temporary folder does not exist.
      2. Allow addition or replacement of @RG record in all files with option
             -r "@RG\tID:..."  Defines an @RG record to add or replace existing @RGs
      3. Added an "--assumedsorted" option that assumes all input files are already
         sorted even if there is no @HD record showing coordinate sorted. Files are
         passed directly to the merge phase. No check is made on the order of the files.
      4. When merging bam files it is not necessary that all the bam files have the
         same @SQ entries or the same ordering of the @SQ records. The output order
         is based on the order in which the @SQ records are first seen.
      5. Added an option to set the compression level for temporary files. This allows trade off 
         between CPU & IO time for temporary files. Default is level 3 (previously 1)
         [--tmpcompression |-x][0-9]    Set compression level for temporary BAM in range 0-9.
                                        Defaults to 3. Higher compression levels reduce temporary
                                        file size and hence IO time at the expense of CPU time.
                                        Suggested range is 0 to 3.

Novoutil iupac
      1. Add option ( -g ) to substitute genotype for the reference base. Without 
         the -g option ambiguous bases are formed by combining the reference base 
         and the genotype.
      2. Add option -i to substitute homozygous indel calls into the new reference. 
         This is useful for reference guided assembly.
      3. Add option -q 99 that sets a lower quality limit (default 30) on usable VCF 
         data lines.
      4. If the REF column of VCF file shows an N then IUPAC code is formed from the 
         ALT alleles.
 

Release Novosort (beta2)
------------------------
      1. Fix problem with reading from stdin
      2. Fix problem where we could exceed the Linux limit on open files per process
         when sorting very large files or when buffer memory was limited. The maximum 
         number of files opened is now less than two times the number of CPU cores.

Release Novoalign V2.07.18 & NovoalignCS V1.01.18
-------------------------------------------------
Novoalign*
      1. Add option --rOQ
         For SAM report with quality calibration, write original base 
         qualities as OQ:Z: tag.
      2. Add option --rNMOri  
         For SAM format report, if a read is unmapped then report the original 
         read and qualities before any hard clipping or quality calibration.
         This facilitates the extraction of unmapped reads from the SAM report
         in their original form.
      3. Fix: When a fasta file was used as input and the file name did not have a 
         '.fa' suffix, then Novoalign fails if there is no '.qual' file.
      4. Fix: Alignment quality could be underestimated in the case where there 
         were multiple possible alignments (e.g. Align with a deletion or align 
         with multiple mismatches) at the same locus with similar alignment 
         scores and most base qualities were <= 10.
      5. When using a BAM file for input and in paired mode any single end reads will 
         be skipped. Previously it caused Novoalign to stop with an appropriate 
         error message.
Novoindex
      1. Fix: Checksum errors were occuring on the novoindex files due to use of 
         unitialised data in unused fields of the index.
         The checksum function on Novoindex files was added in V2.07.15 to detect
         corruption of index files due to file system errors. 
Novoutil
      1. Fix: Novoutil bgzf would not find the novoalign.lic file if it was in the 
         same folder as novoutil and folder path was specified when starting novoutil. 
         e.g. ~/bin/novoutil bgzf ...
      2. Fix: Novoutil bgzf was not writing the bgzf EOF marker which meant some BAM
         readers gave a warning about truncated files.
Novosort
      1. Beta release of a new multithreaded BAM coordinate sort/merge function
         Benefits:
         * Reduced wall time for sort/merge of BAM files
         * Can include strand in sort key
         * Picard like handling of @RG & @PG records

         novosort [-c threads] [-s] [-t tmpfolder] [-m memory[G|M|K]] [-compressionlevel] input bam files...  >sorted.bam

         Examples:
             novosort in.bam >sorted.bam
             novosort -c 4 -s -t ./ -m 4G -1 in1.bam in2.bam >sorted.bam
         This function requires a novoalign or novosort license to operate in multithreaded mode.


Release Novoalign V2.07.17 & NovoalignCS V1.01.17
-------------------------------------------------
novoutil
      1. Fix: iupac function was not writing the modified reference sequence
      2. Fix: iupac function, if multiple vcf SNP records are present for the same locus then
         IUPAC code includes all possible SNPs. Previously on firts SNP for a locus was used.
novoalignCS
      1. Will now operate in single thread mode without a license file.

Release Novoalign V2.07.16 & NovoalignCS V1.01.16
-------------------------------------------------
novoutil
      1. Add utility function iupac to merge VCF file of SNPs into a fasta reference sequence as IUPAC ambiguous codes.
      2. Add multi-threaded utility function bgzf to compress a file into BAM, tabix, or gzip format.

Change History
Release Novoalign V2.07.15 & NovoalignCS V1.01.15
-------------------------------------------------
Novoalign*
      1. FIX: When using unaligned BAM input and SAM output the first character was being 
         erased from read headers.
      2. Limit alignment quality to a maximum value of 70 except for -r Exhaustive.
      3. Added an option to lock the reference genome index in RAM. Use option --LockIdx. 
         This only applies when using a memory mapped index.
      4. For Paired end reads change the default insert size option from -i 250,30 to -i 250,50
      5. Fix Picard validation error where mate alignment locations didn't match the mate. 

NovoalignMPI
      1. FIX: Changes in V2.07.14 for unaligned BAM support caused  errors when passing 
         @RG record to slave processes.

NovoalignCS
      1. Fix: Seg Fault could occur if using quality calibration (-k or -K option). The 
         error was using a wrong index when counting counting mismatches to colours coded 
         as periods '.' and was also getting the position of all colour errors off by one 
         base. The correction improves the accuracy of mismatch counting and the quality 
         calibration function and also changes the quality of base calls near colour errors
         and hence may improve SNP calls especially in low coverage projects. Alignment 
         location and base calls are not affected.
Novoindex
      1. Added a checksum attribute to the index file. This is used to validate correct 
         save/load of the index.
novo2sam.pl
      1. Added conversion for alignment score tag AS:i:99

Change History
Release Novoalign V2.07.14 & NovoalignCS V1.01.14
-------------------------------------------------
Novoalign*
      1. To avoid Picard validation errors adapter trimming will always leave at
         least 1bp in a read.
      2. FIX: In paired end mode with all fragments exactly the same length (usually 
         simulated data) it was possible that floating point errors cause the 
         square root of a -ve number in the calculation of the standard deviation of 
         the fragment lengths.
      3. Fixed issue with iterative alignment that was doing unnecessary iterations 
         for one read when the other read of the pair failed to align with maximum 
         possible alignment score or had been flagged as a low quality read. 
         The problem was evident when aligning pairs against incomplete genomes with 
         many contigs. This change has also improved the accuracy of alignment quality.
      4. Change the default gap extend penalty to 6. Set -x15 to get previous behaviour.
         Using the new default of -x6 may be slower than using -x 15.
      5. Added support for unaligned BAM input. Use option -F BAMSE or -F BAMPE.
         For paired end, mates are expected to be adjacent. If there is an @RG record on
         the BAM file it will be used for the report. Any @RG on -o SAM option overrides 
         the @RG in the BAM file. Not supported for Colour Space reads.
      6. Fix: In SAM report format the NM tag was counting a mismatch between a '.' in
         the read and a 'N' in the reference. This could result in picard ValidateSamFile 
         errors.
      7. Fix: Novoalign uses memory mapping to load the index file and used MAP_POPULATE 
         option to force loading of the index into RAM. Some older Linux Kernals do not  
         support MAP_POPULATE with result that the index pages were not loaded at startup. 
         This could cause slow operation while index pages gradually faulted into memory.  
         We now touch each page at startup to force loading of the index.
      8. Improved run time performance of Bi-Seq strand specific alignments when run with 
         option -b2.

Novobarcode
      1. Add option (--GZIP) to write demuxed reads in gzip format.

Release Novoalign V2.07.13 & NovoalignCS V1.01.13
-------------------------------------------------
Novoalign
      1. Correct documentation for --hdrhd option.
      2. Fix. Reporting stops at first read with Illumina Low Quality Flag set 
         when -o Sync option is used. Applies to CASAVA 1.8 files only.

Release Novoalign V2.07.12 & NovoalignCS V1.01.12
-------------------------------------------------
Novoalign
      1. Fix Seg Fault that occurred if -Q option was used with single end reads.
      2. Add option to collect statistics on homopolymer run length errors.
            --hpstats <filename>
         These can be charted using included R script IONTorrent.R
            IONTorrent.R -f <filename> -r indelcharts.pdf
         The data collection and script will work for 454 Paired & Illumina reads.
      3. Add option to check that paired end read files are in the same order.
         This involves checking that the headers of the two reads in a pair match within
         a specified Hamming Distance.

         --hdrhd [9|off] Controls checking of identity between headers in paired end reads. 
                        Sets the Hamming Distance or disables the check. Default is a 
                        Hamming Distance of not more than 1. Processing will stop with 
                        appropriate error messages if Hamming Distance exceeds the limit.

         The default of 1 should be valid for standard Illumina headers. 
      4. In Paired end mode, previously reads shorter than log4(reference length) would 
         not be aligned, we now honour the -l setting, allowing reads as short as the 
         index k+s-4.

NovoalignMPI
      1. Fix: If the number of reads input to an MPI run was insufficient to transfer 
         at least one buffer of reads to each slave then the MPI process would hang 
         after processing all reads. Buffers are 16Kbyte.
      2. Compile and link using MPICH2 V1.3.2p1
Novoindex
      1. Added an option to control the number of threads used for indexing. Option is -t 9
Novobarcode
      1. Add option (--QSEQ_OUT) to force qseq format output if input is qseq and
         index tags are embedded in reads. If not used then reads are converted to fastq format.
      2. Add option -d folder that sets base folder name for demux'd read files. The folder
         should exist.

Release Novoalign V2.07.11 & NovoalignCS V1.01.11
-------------------------------------------------
Novoalign*
      1. Add an option for Native report format to report bases that match 
         IUB ambiguous codes in the reference sequence. Option is -o IUBMatch
         >Read1	L	CTGTAGTAAAATTAAATTAATTATAAAAT	.	U	32	150	>NC_003663.2	1	F	.	137	R	8R>A 15V>A 23N>A
         This will simplify task of calling SNPs at these locations.
      2. Fix. The manual stated that -h -1 -1 turned off the homopolymer filter 
         however from V2.07.09 it actually set the same as -h 1 1. Code has 
         been corrected so that -h -1 -1 does disable the filter.
      3. In SAM report format add tags NH,IH & HI for reads with multiple alignments.
      4. Port to Solaris 10 is now available on request.
      5. Added support for Illumina Casava V1.8 format reads (-F ILM1.8). These have same
         quality encoding as Sanger Fastq (-F STDFQ) format reads. If V1.8 reads 
         are detected Novoalign will check the header and will (by default) skip reads 
         where is_filtered = 'Y'.
      6. When processing QSEQ format files if Illumina Low Quality indicator is 
         set then by default we will now skip the reads. A count of skipped reads will appear 
         in the log but the reads will not be aligned and will not be written to the SAM report. 
      7. For QSEQ and ILM1.8 format read files you can now specify how reads flagged as 
         Low Quality will be treated. Options are:
             --ILQ_SKIP       Read is not aligned and not written to output report. (Default)
             --ILQ_QC         Read is flagged as QC using Novoalign status and SAM flag bits.
             --ILQ_USE        Quality flag is ignored and read is aligned as per any other read.
         Example:
                  novoalign -d ....  -f ..._sequence.txt -F ILM1.8 --ILQ_USE ....
         If -F option is not specified then Casava 1.8 files are recognised by the header 
         matching regular expression '@*:*:*:*:*:*:* *:[YN]:*:*'.
      8. Fix "Bus Error" that could occur if the gap extend was set to extremely low 
         values of 0 or 1.
      9. In small RNA mode (-m option) and Native report format, changes implemented in 
         V2.07.03 were causing the report to have additional "\t.\t.\t." after location of 
         reverse complement alignment. This has now been removed.
     10. In  small RNA mode (-m option) we now allow a G/U match in the alignment for 
         the hairpin structure. This should improve identification of potential precursor 
         sites.
     11. Fix a malloc exception that could occur when using -v option with multiple penalties.
     12. Fix to paired end adapter trimming in presence of sequencing errors. In some 
         cases adapter was not being recognised and not trimmed from the read, usually 
         resulting in failure to find an alignment.
novoalignMPI
      1. Fix: "src/tcmalloc.cc:387, Attempt to free invalid pointer:" that 
         occured when aligning Illumina mate pair libraries.
Novobarcode
      1. Fix: Demux of qseq.txt files with 3' tags creates incorrectly formatted fastq files.
      2. Allow demultiplex of paired end QSEQ files where the tag read is in a third file
      3. For QSEQ files, reads with Illumina low quality flag set are not classified.
      4. Casava V1.8 format files are automatically recognised and is_filtered == 'Y' reads 
         are not classified.
      5. Add option --NC_OFF which inhibits writing of unclassified reads to the NC folder.
      
Change History
Release Novoalign V2.07.10 & NovoalignCS V1.01.10
-------------------------------------------------
Novoalign*
      1. When reporting mapping locations, we now truncate sequence headers at first 
         whitespace. Prior versions truncated at the first space.
      2. When reporting multiple alignments per read in SAM format, one character was 
         being dropped from next header in the CC tag. This caused CC tag to be null
         if the reference sequence header was only 1 character long.

Change History
Release Novoalign V2.07.09 & NovoalignCS V1.01.09
-------------------------------------------------
novoindex
      1.  Fix problem with building indexes with more than 8Gbp of reference sequence.
novoalign*
      1.  Add option to append text to every read header in SAM and Native report formats.
          '-o Header text'   will append text to every read header. Some users requested 
          this to add slide/lane data to read headers in colour space.
      2.  For Bisulphite alignments change default homopolymer filter setting to 120 
          (i.e. default is -h 120)
novoalignMPI
      1.  Fix: When '-K Filename' was being used to write quality calibration file 
          the calibration data for read 2 of pairs was incorrect.
novoalignCS
      1.  Fix a Seg Fault that could occur if read lengths were less than the 
          index k-mer size. Short reads resulted from trimming adapter sequences.
novomethyl
      1.  Fix interpretation of mpileup base qualities when base string contains 
          a $. This had caused bases and qualities to become misaligned. Novomethyl
          should probably be rerun to correct results.
novo2sam.pl
      1.  Add capability to convert Bi-Seq alignments.

Change History
Release Novoalign V2.07.08 & NovoalignCS V1.01.08
-------------------------------------------------
novo2sam.pl
      1. Fix for single end conversion, previously was setting paired end flags.
novoalignCS
      1. Fix Assertion failed: (tgt[tposn] != '\0') which could occur if read 
         aligned at the end of a reference sequence.

Change History
Release Novoalign V2.07.07 & NovoalignCS V1.01.07
-------------------------------------------------
Novoalign*
      1. If a read is both soft clipped and hard clipped on the 5' end of alignment
         earlier versions would put the soft clipping in the cigar before the hard 
         clipping. The hard clipping should be first. This is fixed in this release.
      2. Performance improvements for paired end read by optimisation of search process.
         a) During iterative search Novoalign gradually increases error tolerance for 
            each read of a pair until a mapping is found. The choice of which read
            to increase next has been optimised to reduce the cost of each iteration.
         b) In some cases Novoalign aligns a pair but because the mate alignment has a
            very high alignment score it reports a mapping for one read and No Match
            to the mate. Occasionally this test failed to report pair alignments
            that were within the search space. This has been corrected and you may 
            see more pair alignments (0.05% in tests on 48bp reads) where one read of
            the pair has a high (+200) alignment score.
Novoalign*MPI
      1. For MPI versions, when running multi-threaded with exactly one thread per
         CPU core (default) we now set processor affinity for each thread to force
         specific CPU per thread. This overcomes problems with Linux CFS scheduler
         where several cores may be idle while running a single multi-threaded job.
      2. Increased the message buffer size from 2K to 16K bytes.
      3. Fix problem with license checking for Bisulphite mode.
      4. Add option, --mmapoff, to disable memory mapping of the index file. This 
         forces each instance of NovoalignMPI to loads its own copy of the index 
         and can help performance on servers with NUMA memory subsystems.
NovoalignCS
      1. When calling colour space alignments, if the quality of a base call is <3
         we now put an 'N' in the nucleotide read sequence rather than the reference base.
Novobarcode
      1. Added support for files in qseq format. Qseq files are supported for input 
         however demuxed files are written in fastq format.
      2. Treat carriage return characters in the tag file as end of line. This allows
         tag files edited on MS Windows systems with CR/LF line delimiters to be used.
      3. Changed default for -t option from 30 to 30*Distance/2. Where Distance is 
         read from the tag file.
      4. Strip ".gz" or ".bz2" off filenames when creating output files. Compressed 
         files are supported for input (with license) but output is always uncompressed.
      5. Add option for adapter trimming to be used with 3' index tags. Option -a [adapter sequence]
         Adapter sequence is appended to tags before checking reads against tags.
         Tag plus adapter are then trimmed form the read. Requires 3' tags and single end reads.


Change History
Release Novoalign V2.07.06 & NovoalignCS V1.01.06
-------------------------------------------------
Novoalign
     	1. For Bi-Seq alignments in SAM report format we add tag ZB:Z:mode where 
         mode is either CT or GA and indicates which mode/strand the read was aligned 
         in. Reads aligned in CT mode are usually from the 5'-3' strand 
         of the chromosome. A simple methylation analysis pipeline could be 
         constructed by splitting alignments into two SAM files using the ZB 
         tag and then running samtools pileup on each file. Non-methylated 
         cytosines should show up as 'T' on the CT alignments and as A on
         the GA alignments. More details on our wiki at http://tinyurl.com/Novocraft-BiSeq
      2. Change to automatic fastq format detection to allow -F STDFQ even if 
         quality range looks like Illumina fastq format. We recommend always using 
         the -F option to ensure correct interpretation of quality values.
      3. For Bi-Seq alignments if a read aligns in the same location and direction and 
         with same score in both GA & CT modes (typically there are no unmethylated 
         cytosines) then we choose randomly whether to report as CT or GA aligned. 
         Earlier versions would have biased reporting to the CT mode. Note. This has 
         no effect if your protocol has preserved strand and you are using the -b2 
         option which should be the case if you are using the Illumina protocol.
      4. In SAM report format when reporting multiple alignments per read and 
         one read of a pair is unaligned then the mate location is now shown as the 
         primary alignment location.
      5. Support for read files compressed with bzip2. i.e. *.bz2 files.
      6. Change to alignment seeding to allow drop back to a single seed location if
         other seed locations have too many low quality bases. More mismatches are 
         allowed in the single seed so that sensitivity is not affected.
      7. Add option --3Prime that enables reporting of 3' locus of alignments. This 
         will appear in SAM files as tag Z3:i:9999
Novoalign*
      1. When running multi-threaded with exactly one thread per CPU core (default)
         we now set processor affinity for each thread to force specific CPU per 
         thread. This overcomes problems with Linux CFS scheduler where one or two 
         cores may be idle while running a single multi-threaded job.
Novomethyl
      Beta release of methylation status caller. Please refer to our wiki at 
      http://tinyurl.com/Novocraft-BiSeq


Release Novoalign V2.07.05 & NovoalignCS V1.01.05
-------------------------------------------------

In this release we have addressed some Picard validation problems. One issue is that
alignments can be hard clipped by -H or -a option and also soft clipped as a result
of Smith-Waterman alignment. In Picard V1.35 and earlier this was flagged as an error.
Picard V1.36 should fix this.

novoalign*
      1. In command line processing we now validate that readgroup record 
         includes  ID: tag. Applies to -o SAM ["@RG.."] option.
	  2. Fix for Picard validation failure when -Q option was used. Alignments 
         were being reported as mate was mapped but then it's mate was reported as 
         unmapped if its alignment quality was below the -Q reporting limit. The 
         -Q option is really redundant for use with SAM format and we may remove 
         it in a future release.
      3. Correct problem where insert size was +1bp for non proper pairs that 
         aligned on same chromosome.
      4. Fix Picard validation problem where mate alignment location (MRNM) may
         differ by 1bp from actual mate alignment.
      5. In miRNA mode and SAM report format (option -m -oSAM) add custom tags ZH:i:
         & ZL:i: for hairpin score and alignment location.

novo2sam.pl
      1. Corrections to CIGAR & MD fields for case where alignments are softclipped.

novobarcode
      1. If read folder path included ..'s then the process of creating folders and files
         for demuxed reads could fail. This change creates tag folders in current
         folder and places readfiles directly into these, ignoring the folder of the
         original read files.

Release Novoalign V2.07.04 & NovoalignCS V1.01.04
-----------------------------------------------
novoalign*
      1. When using -H option to trim low quality bases from reads if a read would be clipped to length <= 1
         then we don't clip, leaving the read intact. It will still get a QC alignment status.
      2. In SAM report format, if due to adapter trimming a read is clipped to length 1bp and the base
         has a quality code of "*" then we convert the quality code to ")". This avoids an ambiguity 
         in the SAM specifications.
      3. Fastq file format detection now defaults to Illumina coding ("A" + q) unless lower quality 
         values indicate Solexa or Sanger formats. It is still advisable to use -F option to ensure 
         correct decoding of quality values.
      4. In command line processing we now validate that readgroup record starts with @RG. 
         Applies to -o SAM [readgroup] option
      5. Fixed a problem in paired end mode SAM report format where occasionally two reads of a pair 
         would have different fragment lengths reported.
	  6. Fix a seg fault in bisulphite mode that occurred if index k+s was >= 21
      7. Fix problem introduced in V2.06.10 where iterative alignment search was extended too far 
         for paired reads that did not align as a proper pair. This will restore runtime performance 
         to V2.06.09 level or better.
      8. Added command line option --Q2Off to disable treatment of Q=2 bases as Illumina 
         "Read Segment Quality Control Indicator". Setting Q2 off will treat Q=2 bases as normal 
         bases with a quality of 2.  When off Q=2 bases are included in quality calibration. By 
         default it is off in NovoalignCS.
      9. When Paired end adapter trimming was used  with -H option and with mixed mate pair/paired end reads,
         earlier versions attempted to trim both reads to half the fragment length. This could remove 
         too many bases if one read had been hard clipped to less than half the fragment length. 
         Trimming now removes just enough bases to remove any overlap between the two reads.

novobarcode
      1. Fix to correctly strip index tag from  Illumina fastq format files. This applies when 
         index tag is part of read and is has no affect when index tag is in the header.
      2. Change so that tag folders and files are only created for tags that have reads.
      3. Fix to use base qualities when tag is in reads.

Release Novoalign V2.07.03 & NovoalignCS V1.01.03
-----------------------------------------------
novoalign*
	1. Fix *** glibc detected *** novoalign: munmap_chunk(): invalid pointer: 0x0000000000d7f620 ***
         which happens at end of some mate pair runs and was caused by a one byte buffer overrun. 
         Results were not affected by the fault.
      2. Added basic validity checking on novoindex file before attempting any alignments. If supplied 
         index file is not recognised as a novoindex file the program will now stop with an appropriate 
         message. In earlier versions it usually crashed with a seg fault.
      3. Reduced memory utilisation by changing to a memory allocator (tcalloc) with a thread cache.
         This greatly reduces the memory growth seen in NovoalignCS runs.
      4. Fix for abort:
            terminate called after throwing an instance of 'std::length_error'
 		what():  vector::_M_fill_insert
novoalign
      1. miRNA mode and Hair Pin changes.
          a) Hair Pin score is now part of sort key when reporting alignments. This will affect alignment
             reporting order for -r Random, All & Exhaustive
          b) Mate location fields in report are now used to report alignment of reverse complement of the
             read giving an indication of orientation of precursor molecule.

novoalignMPI & novoalignCSMPI
      1. Fixed a seg fault that could occur on startup in single end mode

isnovoindex
      1. Change return status
            -1 - not a novoindex
             0 - Normal nucleotide space index
             1 - Bisulphite nucleotide index
             2 - Colour space index
             3 - (reserved for Colour space bisulphite index)
         previously status 0 was for any valid novoindex.

Release Novoalign V2.07.02 & NovoalignCS V1.01.02
-----------------------------------------------
 novoalign
      1. Novoalign does not accept sequence files in SCARF format. In previous versions these files 
         could be misinterpreted as PRB format files. This change should report a file format error 
         for SCARF files.
      2. Runtime performance improvements, with neutral affect on sensitivity and specificity, of around
         40% for Bi-seq alignments when not using -b2 option. 
    
novoalign*
	1. Fix Picard validation errors
	   a) When one read of pair is unmapped set MRNM & MPOS in mapped read to location of mapped read.

novo2sam.pl
      1. Convert softclipped alignments from Native to SAM format


Release Novoalign V2.07.01 & NovoalignCS V1.01.01
-----------------------------------------------
novoalign*
	1. Allow named pipes to be used as read sequence files. This disables file 
         format checking and requires that the -F option is used to specify the 
         sequence file format.

Release Novoalign V2.07 & NovoalignCS V1.01
-----------------------------------------------
novoalign 
	1. Extended support for Illumina Mate pair libraries
         a) Can have mixed mate pair & paired end alignments by specifying long
            and short fragment lengths on the -i option.
            novoalign -d ....    -i MP 2500,500 250,60 
         b) When both short & long fragment lengths have been specified Novoalign
            can detect read pairs where circularisation junction is inside one of 
            the reads and split the read at the junction. The alignment to longer 
            portion of the split read will be reported.
         c) Paired end, short fragment, adapter trimming will trim adapter and 
            shorten reads (hard clip) to remove overlap. This allows (b) to work 
            on short fragments.

Release Novoalign V2.06.11 & NovoalignCS V1.0.11
-----------------------------------------------
novoalign *
	1. Increased alignment match reward from 6 to 12 for improved snp and indel calling 
         with softclipped alignments.
	2. Polyclonal filter is applied after trimming low quality bases. This helps stop QC
         filtering of reads with trailing qualities of #

Change History
Release Novoalign V2.06.10 & NovoalignCS V1.0.10
-----------------------------------------------
novoalign*
	1. Changes to -i option to allow mate pairs and paired end reads. Mate pairs and
         paired end reads have different orientations when aligned on the genome. Paired 
         end reads are on opposite strands in +- orientation. Illumina mate pairs are 
         reverse of paired end with 5' read on -ve strand and 3' read on +ve strand. 
         ABI Solid mate pairs are in same orientation, either ++ or --. Changes to -i 
         option allow expected orientation to be set and ensure proper operation of paired 
         end alignment.
         Insert lengths can be specifed as a range rather than as a mean and standard deviation.
         Option format -i 999-9999 specifies a range of fragment lengths
         Option formats -i  99 99 or -i 99,999 specify a mean and standard deviation.
         When specified as a range, fragment length penalties are not applied, which may reduce
         the ability to resolve alignments near tandem and satellite repeats.
         The option is now defined as -i [PE|MP|++|+-|-+] 99[-|,| ]99
         Examples:
         -i 250 50     Default to paired end Illumina (+-) or Mate Pair ABI (++) with 250bp insert
         -i PE 250,50  Uses paired end orientation (+-) with 250bp insert
         -i MP 2000,200  Uses mate pair orientation (Illumina -+, ABI ++)with 2000bp insert
         -i ++ 200-1000  Uses ++/-- orientation with fragment lengths from 200bp to 1000bp
         Proper setting of orientation is important. If in doubt about mean fragment length 
         and standard deviation err on the high side.
      2. Changes to -p option to use the fraction of bases rather than number of bases below
         quality limit for full read filter option.
         This makes it easier to set the when the two reads of a pair have different lengths.
         Example  -p 7,10 0.2,10  would flag as low quality any read with 7 or bases 
         in first 25 below Q10 or more than 20% of bases in whole read below Q12
	3. Changes to MD tag in SAM reports so that a multibase deletion is reported with a
         count of matching bases and a caret between each deleted base. Example a 3 base delete
         of ACT was reported as ^ACT and is now reported as ^A0^C0^T
      4. Change to ISIZE field in SAM report so it is based on unclipped positions of the
         read alignments. This pertains to softclipped bases, not hard clipped bases. Bases hard
         clipped by functions such as paired end adapter trimming are not included in ISIZE.
      5. Reporting option (-o Softclip) to softclip alignments back to best local alignment
         is now on by default for SAM report format. This is to avoid problems with GATK toolkit
	   and also improves SNP calling with sampttols.pl varfilter.
         You can turn off soft clipping with the option -o FullNW.
      6. Illumina use a base quality of 2 to indicate a problem with calling of a read.
         "The Read Segment Quality Control Indicator: At the ends of some reads, quality scores
          are unreliable. Illumina has an algorithm for identifying these unreliable runs of 
          quality scores, and we use a special indicator to flag these portions of reads with a 
          quality score of 2, encoded as a "B", is used as a special indicator. A quality score
          of 2 does not imply a specific error rate, but rather implies that the marked region 
          of the read should not be used for downstream analysis. Some reads will end with a 
          run of B (or Q2) basecalls, but there will never be an isolated Q2 basecall."
         Q=2 is Perr = 0.63 and in earlier releases were scored as 4 for a match and 7 for 
         a mismatch so they contributed slightly to alignment scores.
         The net effect was a slight increase in SNP noise. In this release Q=2 bases are treated
         as N's.
         We also add an option -H that will hard clip trailing bases with Q<=2 from reads.
      6. Implemented Hard clipping of trailing bases (option -H to enable) with Q <=2 in response
         to Illumina use of this as a special indicator.
      7. Added option -# 999[K|M]  that sets a limit on the number of reads to be processed.
         e.g. -# 1k will stop after 1000 reads or pairs have been processed. This option is 
         useful for doing a test run of a few alignments.

NovoalignMPI
      1. License checking is now only done on the master process. It is not necessary to 
         have the license file available on the slave servers.

novoalignCS
	1. Fix - If first read in a csfastq format file had a '.' in colour string then
         file type may be set incorrectly causing either the program to stop or the
         colour qualities to be offset from the base positions.
	2. In Native report format the strand of the mate alignment was inverted.
      3. Correct a problem in alignment of long (> 7bp) deletions in paired end reads.
      4. When one read of pair fails quality control checks we may still align it by 
         anchoring it to the alignment of the other read of the pair. We now only report
         the alignment to the low quality read of the pair if it's alignment score 
         satisfies quality limits (expected Q>30)
      5. Changed default settings of polyclonal filter to -p 7,10 0.3,10. This should improve 
         performance with minimal impact on the number of alignments.
      6. Changed called base qualities to be based on colour qualities rather than alignment
         penalties. This will change quality of last base by 4 or 5 points.
      7. Implemented soft clipping of alignments. This trims alignments back to best 
         local alignment.

novoalignCSMPI
      1. Initial Release of MPI version of NovoalignCS

Release Novoalign V2.06.09 & NovoalignCS V1.0.9
-----------------------------------------------
novoalignCS
      1. Performance improvements for mate-pairs (order 20-30%)
novoalign
	1. Fix for difference between code generated by GCC 4.4.1 vs earlier versions of 
         GCC and the Sun compilers. GCC compiled versions of Novoalign since Nov 2009 
         may (P~0.1) have incorrect alignment of inserts over 7bp in length. 
         The effect was most noticeable on alignments that extend passed the end of sequence. 
         In Sun compiled version, the extra bases are shown as inserts after the last base 
         in the sequence. In GCC compiled versions the inserts are position 1 bp before 
         the end of the sequence. It may have also caused some alignments to long 
         inserts to be missed. It did not show up as a problem in simulated reads 
         designed to test indel detection.

Change History
Release Novoalign V2.06.08 & NovoalignCS V1.0.8
-----------------------------------------------
novoindex
	1. Novoindex may cause 1 or 2bp to be added to the length of reference sequence
         in SAM format @SQ records. This could happen to the last fasta sequence in a 
         file if multiple fasta sequences files were used with novoindex. The 
         problem would not happen if all sequences were in a single file. This 
         is now fixed and affected indexes should be rebuilt.

novoalign*
      1. It was possible that for paired end reads if only one read of a pair aligned 
         that the read was not reported. This affected about 1/1000000 reads. 

novoalignCS
	1. Fixed a problem in SIMD alignment routine that was affecting alignment of reads
         with deletions relative to the reference.
      2. Colour Mismatch counts (SAM tag CM:i) were +1 on some alignments.

Release Novoalign V2.06.07 & NovoalignCS V1.0.7
-----------------------------------------------
novoalign
Change History
novoindex
	1. Skip any space characters in the reference sequence. Previously they
         were included as not ACG or T. All other characters in the reference 
         sequence that are not valid IUB NA codes are included in the index as 
         not ACG or T
novoalign
	1. Adapter trimming can trim read sequences to zero length if the read 
         starts with an adapter or if all the bases in the read are phred 
         quality 1 or 2 (equivalent to N). When reporting in SAM format earlier 
         versions would ouput an empty sequence which could cause Picard validation 
         to fail. In this version if a read is trimmed to zero length we 
         output a '*' in the sequence field.
	2. Fix to a Seg Fault error introduced in V2.05.33
novoalignCS
	1. Enable multithreading, quality calibration and Gzipped files whenever 
         colour space is enabled.
      2. Apply alignment penalties for match to IUB Ambiguous codes in the 
         reference sequence: N=6, V,H,D,B = 5, remainder 3.
	3. CSFASTQ files with primer base qualitywere being processed incorrectly 
         with result that alignment was off by 1bp. This is now corrected.
      4. Fix to a Seg Fault error.
      5. In SAM report format the CQ tag should always have a quality for 
         primer base even if input read file didn't. In NovoalignCS the primer 
         base quality is always '!'.

Release Novoalign V2.06.06 & NovoalignCS V1.0.6
-----------------------------------------------
novoalign
	1. Paired end adapter trimming could result in alignments with
         invalid insert at the end opposite the trimming. Problem was introduced in
         V2.05.33 and partially fixed in V2.00.05.
         It is still possible to get inserts called at the end of read: 1) if the
         read overlaps the end of a reference sequence, 2) where calling a longer
         indel has a lower penalty than calling mismatches. 
         Some programs in samTools/Picard may not handle this. Consider using
         option -oSoftclip to trim alignments back to best local alignment.

Release Novoalign V2.06.05 & NovoalignCS V1.0.5
-----------------------------------------------
novoalignCS
	1. Allow colour space fastq with/without primer base quality value. The first read 
         is checked to see if quality string is the same length as the read or 1bp shorter.
         If the same length then first quality is treated as primer quality.
	2. Fix: In V1.00.04 a csfastq read were being aligned with an extra base.

Change History
Release Novoalign V2.06.04 & NovoalignCS V1.0.4
-----------------------------------------------
novoalignCS
      1. Allow reads up to 100bp
      2. SAM format RNAME truncation now at first white space rather than first space.
	3. Don't use first colour in a seed.
	4. Fix seg fault in seed selection routine.
novoalign
	1. Paired end adapter trimming could result in alignments with
         invalid insert at end opposite the trimming. Problem was introduced in V2.05.33


Release Novoalign V2.06.03 & NovoalignCS V1.0.3
-----------------------------------------------
novoalignCS
	1.  NovoalignCS passed Valgrind
	2.  Corrections to csfastq format to allow for the primer base quality. 
      3.  New SIMD alignment routine for mate pairs may improve performance.
      4.  When parsing read headers to match pairs the third numeric field need not 
          be terminated by an underscore.
novoalign*
	1.  When using -s option to progressively trim reads and realign, the final 
          result is now either an aligned read or a read with 'NM' status. Previously 
          it was 'QC' status when a read was trimmed below the quality control limits.
      2.  Changes to RNAME processing in SAM format reports
          Single End - from Read header truncated at first space.
          Paired End - from Read header truncated at first space or first character where 
          headers of the pair differ and then a single trailing underscore '_' or 
          slash '/' is removed if it exists.

Release Novoalign V2.06.02 & NovoalignCS V1.0.2
-----------------------------------------------
novoalign*
	1.  All aligners, fixed an intermittent seg fault.
	2.  Change default for Polyclonal filter to 8bp below Q10 in first 20bp.

Release Novoalign V2.06.00 & NovoalignCS V1.0.0
-----------------------------------------------
novoalignCS
	1. Initial release of aligner for ABI Solid Colour Space reads.
	2. Change default for Polyclonal filter to 10bp below Q10 in first 20bp.

novoalign
	1. Fragment length penalties could be set incorrectly if standard deviation was set to less than ~8% of the mean.
         This resulted from an underflow in some calculations. The problem is fixed in this release.
         The correction may increase the alignment score of some pairs especially for fragments shorter than the mean 
         by 10 or more standard deviations.
	2. Increased the resolution of fragment length histogram to support long insert libraries. There may be a small 
	   change (1-2) in alignment scores of paired end reads. 
	3. Performance improvement of 30-50% for normal DNA/RNA reads (Unfortunately it doesn't help Bi-Seq Alignments)
	4. Addition of a filter for Polyclonal reads. This filter flags as QC any reads that have more than 7bp of the first 
         20bp with a phred quality below 10. This filter is based on the paper "Filtering error from SOLiD Output" by 
         Ariella Sasson and Todd P. Michael. Command line option -p can be used to change settings or to turn off the filter.
         This filter is more applicable to ABI Solid reads than Illumina reads and defaults will usually only filter a 
         handful of reads from an Illumina lane.
         However if you are using Novoalign to align unaligned reads from your normal aligner, then its 
         worth considering using this filter to remove low quality unaligned reads as they will just slow down alignment
         process.
novoindex
	1. Fix a problem with code that determines memory on MAC OSX computers. Memory estimate is used to when defaulting
	   k&s to determine the maximum possible settings that will fit in memory.

novo2sam.pl
      1. Correction to convert alignments with status of 'R' if these had been reported in the novoalign run.

Release 2.05.33 - 21 Apr 2010
----------------------------------
novoalign
	1. Performance improvement. Paired end processing can find duplicate alignments, this
         change improves performance of the code to filter duplicates.
      2. With paired end short fragment detection, adapter is trimmed from reads and 
         then the reads aligned as a normal pair, in presence of tandem repeats it was possible
         that the two sides could align as a longer proper fragment. This restricts the two
         reads of a short fragment to align to the same location.
		Reference ---------repeat copy1----------------------repeatcopy2-----

		Valid Alignment    ---Read1--->
                               <--Read2----
            And also                                              ---Read1--->
                                                                  <--Read2----

            Invalid            ---Read1--->.......................<--Read2----  (not possible any more for short fragments)

Release 2.05.32 - 30 March 2010 (Updated 21st Apr 2010)
----------------------------------
novoalign
	1. Change to default adapter sequence for paired end adapter trimming. The read1 
         default adapter sequence is now "AGATCGGAAGAGCGGTTCAGCAGGAATGCCGAG". The change
         will have no effect as Novoalign only uses the first 12 bp and these haven't changed.
         Please refer to our Wiki for a description of paired end adapter trimming.
      2. In 2.05.31, the correction for 'Seg Fault' error was to adjust alignment threshold 
         when reads had adapter trimmed from them. This change also had unintended effect of 
         reducing threshold even when the read wasn't trimmed. This has been corrected.
	3. Fixed a problem in SAM report format where the MATE_STRAND bit (0x0020) could be set
         incorrectly for read pairs that failed to align as a proper pair.

Release 2.05.31 - 19 March 2010
----------------------------------
novoalign
	1. SAM format, fix problem with @RG tags where the colon in the RG:Z: tag was duplicated.
	2. Fix problem in dinucleotide filter when used with filter scores that allowed gaps longer than
         7bp in the alignment. Previous versions may have filtered out some dinucleotide reads
         whose alignmnent score was higher than the specified threshold for the filter.
      3. Fix problem with use of -rAll and -rExhaustive to report multiple alignments per read
         where occasionally multiple alignments were reported for the same location
         but with differing edits. One alignment might have an indel and another mismatches.
         This fix may also improve quality scores for some alignments.
	4. Native report format - Added reporting of the number of bases soft clipped from 3' end
         of alignment for soft trimmed alignments.
      5. Fix assert failure when using soft clipping of alignments.
      6. Fix seg fault when using paired end adapter triming and read trimmed to less than index k-mer length.
      7. In SAM format correct the MD field so that it conforms to the specified regular expression  "[0-9]+(([ACGTN]|\^[ACGTN]+)[0-9]+)*"
         This required addition of 0 counts of matched bases between mismatches and at the end of the MD tag.
	   Earlier MD:Z:31C30T3^CA3TC6TCT new MD:Z:31C30T3^CA3T0C6T0C0T0 
	8. Allow more threads to be created using -c option than there are CPUs on the system. Previous versions
         limited threads to sysconf(_SC_NPROCESSORS_ONLN) which may be incorrect on some systems running hyper
         threading.
      9. Allow MS Windows format text files as input (i.e. With CR/LF line separators)

novoindex
	1. Force s=2 for reference genomes >4Gbp. Minimum s = INT((ReferenceGenome Size)/4^16);

Change History
Release 2.05.25 - 22 Feb 2010
----------------------------------
novoalign
	1. In SAM report fix problem where unaligned reads were being flagged with QC status when
         they should have been given NM status (Novoalign custom ZS:Z:... tag)
	2. Soft Trimming (-o Soft) of alignments was added. Alignment locations are still determined
         using a global Needleman-Wunsch alignment of the read, this option trims alignments back 
         to the best local alignment. This helps remove spurious SNP & Indel calls near the ends 
         of the reads and can improve specificity of SNP and micro Indels.
novo2maq
	1. Fix a problem where the list of sequence headers printed using -r option contained spurious entries.
      2. Fixed a problem where novo2maq wouldn't write a map file with uncompressed file size >4GB (Note. This problem 
         also exists in MAQ 0.7.1  novo2maq and eland2maq)

Change History
Release 2.05.24 - 5th Feb 2010
----------------------------------
novoalign
	1. When printing multiple alignments for a read with -r All or Exhaustive options and 
         limiting the number of alignments reported (e.g. -r Exhaustive 100) then the alignments 
         reported are now selected randomly.
      2. The MAC OS X version is now compiled and linked for MAC OSX 10.5 rather than 10.6. It 
         should run on 10.5 & 10.6 systems.

Release 2.05.23 - Not Released
----------------------------------
novoalign
	1. Performance Tuning

Release 2.05.22 - 22nd Jan 2010
----------------------------------
novoalign
	1. In SAM report format, paired end reads, -rNone, in the case where alignments did
         not form a proper pair, one read aligned uniquely, and the other read of pair aligned 
         to more than one location we now set flag 0x8, "Mate is not mapped" on the read that 
         aligned uniquely as we will not report an alignmnet location for the mate.

Release 2.05.21 - 22nd Jan 2010
----------------------------------
novoalign
	1. Add extra tags, ZS & ZN, to SAM format report line for Alignment Status
         and count of alignments found.
      2. Fixed Seg Fault that occurred occasionally for reads that failed to align 
         as a proper pair, one read of the pair aligned to multiple locations, and 
         random reporting of multiple hits was selected.
      3. Fixed a problem with -rNone and paired end in the case where: alignments did
         not form a proper pair; one read of pair aligned uniquely; the other read 
         aligned to more than one location, and
         -RNone was selected. In this situation Novoalign was reporting a mate location 
         for the read that was uniquely aligned.

Release 2.05.20 - 10th Jan 2010
----------------------------------
novoalignMPI
	1. Official release of Novoalign MPI version (available on request only)

 
Release 2.05.19 - (not released)
----------------------------------
novoalign
	1. Fixed seg fault in quality calibration when as called base quality was > 40.
	2. Improvements to performance when using quality calibration.
	3. Allow Windows format read files (CR/LF line delimiters)
	4. Refactoring to support upcoming MPI version.
        5. Change quality calibration to factor in reads that failed to align. These are 
           taken as examples of correct quality calls. Quality calibration files are not
           compatible with earlier versions.
novoalignMPI
	1. Beta release of NovoalignMPI

Release 2.05.18 - (not released)
----------------------------------
novoalign
	1. Change to SAM report format to conform with "... when one read of a pair is mapped and the other unmapped
	   the common practice has become to use the RNAME/POS of the mapped mate." in the unmapped SAM record.
	2. Improved multi-threaded performance by moving report formatting from control thread to alignment threads.
	3. Corrected problem in SAM report format where /1 & /2 could be left on read headers if one read of pair had
         less than 16 bits of information. (i.e. approx 8 good quality bp)
      4. Added a variation on Native report format that adds 3 extra columns: Count of bp trimmed from 5'end of read;
         count of bp trimmed from 3' end of read and; a proper pair flag.
         The trimmed bp counts appear before alignment status field and the proper pair flag is after the alignment quality.
         The new format is selected using option -o Extended
      5. Modified SAM format CIGAR to included any bp hard trimmed from the read.


Release 2.05.17 - 5th November 2009
----------------------------------
novoalign
	1. When running multithreaded default output is now asynchronous with read files. If you need 
	   output in sync with reads add option -oSync to the command line.
	2. The homopolymer filter now also filters reads that are dinucleotide repeats, this helps
         reduce run time when aligning against genomes with high dinucleotide repeat content. Default option
         will report perfect dinuc repeats with status of QC. A second threshold specified on the -h option
         applies to dinucleotide repeats.
	3. In SAM report format, tag fields for PG: & RG: are now added to records for reads that didn't align.
         In previous versions they were only added when the read aligned. 
      4. When processing @RG option on for SAM format reports, we now interpret '\t' sequences as tabs
         making it easy to correctly format the @RG record.

Release 2.05.16 - 5th October 2009
----------------------------------
novo2sam.pl
	1. Fix for case where read sequence contains lower case NA codes.
novoalign
	1. Added option in SAM format to define read group, platform and library using an @RG record.
	2. Fixed a problem in reporting of structural variations (not proper pairs) where two 
         identical alignments to the same read could be reported.

Release 2.05.15 - 22nd September 2009
----------------------------------
Novoalign
	1. Fixed problem where 150bp reads were causing a buffer overrun and failing to align.
         150bp is longest read currently handled by Novoalign

Release 2.05.14 - 13th September 2009
----------------------------------
Novoalign
	1. Added support for _qseq.txt format read files.


Release 2.05.13 - 6th September 2009
----------------------------------
Novoalign
	1. Removed sanity tests on values of gap open and extend penalties set using -g & -x options.
         These can now be set as low as you would like but do be careful as it can cause strange results 
         with dynamic programming if gap penalties are less than the mismatch penalties.
         
Release 2.05.12 - 26th August 2009
----------------------------------
Novoalign
	1. Fixed problem in reporting of structural variations. Using paired end reads if two reads aligned
         on different chromsomes or too far apart (when using single SV penalty option) to be a proper fragment 
         then the reads may be reported as unaligned "NM" status.
	2. Fix a problem where a very high default threshold could be calculated if the reference sequence contained
	   IUPAC ambiguous NA codes and no large blocks of N's. This could lead to long alignment times for some reads.

Release 2.05.11 - 11th August 2009
----------------------------------
Novoalign
	1. Fixed problem in paired end adapter trimming where occassionally one base of adapter was left on the reads.

Change History
Release 2.05.10 - 11th August 2009
-------------------------------
Novoalign
	1. Now allows indels of up to 15bp in single end reads. In practice this can only be achieved
	   on long reads and by lowering the gap extension penalty to 10 by adding option -x10
	   Paired end reads still allow indels of 15bp or greater in one read of the pair while the other
         read of the pair is limited to indels <= 7bp.
	2. File format detection now reads up to 20 fastq reads to identify the range of quality values
         and assign a format. However, it is still best to set the format using the -F option rather than
         rely on automatic detection.
	3. Correction to SAM format to reverse complement read and quality strings when alignment is on -ve strand
Novobarcode
	1. Allow length of barcode in Illumina FASTA header to be longer than the barcodes in the tag file.
	2. In previous versions an 'N' in an Illumina format barcode was treated as a mismatch. In this version it is treated
	   as probability P=0.75 of a mismatch and scores 6 in alignment rather than 30. This change
	   improves the ability to classify reads.

Release 2.05.09 - 4th August 2009
-------------------------------
Novoalign
	1. Changes to SAM report format
         a) Set MRNM & MPOS fields in unaligned reads when the mate was aligned.
         b) Truncate @SQ sequence headers at first space to be consistent with alignment records.
	2. Change to  -r option with posterior probability. If option -r 0.01 was used to report alignments
         with posterior probability > 0.01 the test was always done as though reads were single end. 
         This change now uses posterior probability of read pair with paired end reads.

Release 2.05.07 - 4th August 2009
-------------------------------
Novoalign
	1. Performance improvements for large genomes with high repetitive content (e.g. Human and similar) performance
         improvements of 2 to 3 times are seen with some read files.

Release 2.05.04 - 9th July 2009
-------------------------------
Novoalign
	1. Fix - Paired end short fragment adapter stripping was trimming one to many bases from each read.
	2. Correction to read file format identification. In this version any file recognised as a FASTQ 
         format can be changed to another FASTQ format using the the -F option.

Release 2.05.02 - 16th June 2009
-------------------------------
Novoalign
      1. Add a function for stripping adapter off paired end reads where fragment length is less than the read length.
      2. Fix for a problem where occasionally two identical alignments were reported as a repeat alignment. This could 
         happen when there were two alternative alignments within a few bases of each other and rather than report the 
         two different alignments they were both reported as the same.
Novoindex
      1. Fix a problem where occasionally the first line of a reference sequence was included as part of the header.
Novobarcode
      1. Add support for Illumina format index read files with index tag in the read header.

Release 2.04.03 - 9th June 2009
-------------------------------
novoalign
	1. Fix a problem with calibration where Read 1 calibration data was being used
         to calibrate read 2 in paired end mode.
      2. Fix a problem with 3' adapter stripping in paired end mode where the scoring 
         matrices for adapter bases were not always cleared with possible result of 
         false negatives.

Release 2.04.02 - 21st May 2009
-------------------------------
novoalign
	1. Fix Novoalign report where calibrated base qualities were truncated at first N.

Release 2.04.01 - 21st May 2009
-------------------------------
novoalign
	1. Change Novoalign report to show calibrated base qualities when run with -k option.

Release 2.04.00 - 15th May 2009
-------------------------------
novoalign
	1. Addition of a base quality calibration function
	2. A new command line option '-F ILMFQ' for processing read files in Illumina Casava Pipeline 1.3 format.

Release 2.03.12 - 16th Apr 2009
-------------------------------
novoalign
	1. Additional performance improvements for Bi-seq alignments
	2. Allow Paired end alignements to have score as high as 510 (previously limited to 254)

Release 2.03.12 - 16th Apr 2009
-------------------------------
novoalign
	1. Additional performance improvements for long reads and bi-seq alignments
	2. Fix so that structural variations with combined score greater than the threshold are now reported 
         as no alignment found. In previous versons alignments to the two ends would be reported even if 
         alignment scores plus SV penalty was greater than the threshold.
      3. Addition of an optional bisulphite alignment mode where only two alignment searches rather than 
         four are performed per read. In this mode alignments are Forward against CT index and Reverse 
         Complement against the GA index. In normal mode, alignment is Forward and Reverse Complement against 
         both CT & GA indexes.
	4. Novoalign usage now displays the version number.

novoutil
	1. Added an option "n2mhdrs" that will generate a sequence header file for use in novo2maq

novo2maq
	1. Corrections for report format changes in V2.03.01 of novoalign.
      2. Add an option that allows operation without the need to list all reference sequences.
      3. Add an option to report the number of alignments per reference sequence.

novobarcodes
	1. Fix crash when an odd number of index tags was used.


Release 2.03.01  - 17th Mar 2009
-------------------------------
	1. Performance improvements for longer reads.
	2. Change in the way inserts are reported. Old 30+C 30+T 30+G is now 30+CTG

Release 2.02  - 1st March 2009
-------------------------------
	1. Additional options for structural variation penalties. Refer to the manual.

Release 2.01.06  - 17th Feb 2009
-------------------------------
	1. Add an optional penalty for unconverted cytosines at CHH and CHG in bisulfite alignment mode.

Release 2.01.05  - 15th Feb 2009
-------------------------------
	1. Add option to print floating point quality values. This allows more precision in quality values 
         especially with lower values.
	2. Fix for novo2maq when novoalign has been run with the -m option.

Release 2.01.04  - 4th Feb 2009
-------------------------------
	1. Addition of alignment for bisulfite treated reads.

Release 2.0.14  - 21st Jan 2009
-------------------------------
	1. Fix problem in Novoindex when total length of reference sequences are over 4Gbp
	2. Fix a problem in Novoalign when adapter stripping and read trimming are used together.

Release 2.0.13  - 7th Jan 2009
-----------------------------
	1. Improved accuracy of quality scores for paired end reads.

Release 2.0.12  - 3rd Jan 2009
-----------------------------
	1. Fix problem with identification of file formats when used without a license file.

Release 2.0.11  - 20th Dec 2008
-----------------------------
	First release of V2.0, refer to Release notes for changes from V1.0

Release 1.05.03  - 8th Dec 2008
-------------------------------
novopaired
	1. Fix assert() failure that occurred if paired end fragment was shorter than the read lengths and the alignment
	   to -ve strand overlapped the beginning of a reference sequence.
	2. Fix problem with '-r Random' reporting of repeat alignments where occasionally no alignments were reported for a read.

Release 1.05.02  - 13th Nov 2008
-----------------------------
novopaired & novoalign
	1. Fix "Interrupted..4" problem that was occurring on some Intel Xeon CPUs.
	2. Fix novoindex "Interrupted..11" problem that could occur if a reference 
        sequence was shorter than the index word size.

Release 1.05.01  - 3rd Sep 2008
-----------------------------
novopaired & novoalign
	1. Fix segment fault in miRNA mode.

Release 1.05  - 2nd Sep 2008
-----------------------------
novopaired & novoalign
	1. Added multi-thread option in the commercial version.

Change History
Release 1.04  - 14th August 2008
------------------------------
novopaired & novoalign
	1. Added a command line option -e 999 that limits the number of alignments per read. The limit applies to
                the number of alignments with score equal to the best alignment. Once limit is reached no new alignments
                are recorded. Defaults to 100 alignments.
novopaired
	1. Reduces memory usage when a read aligns to a very large number of locations. In previous versions
                excess memory would be used possibly causing out-of-memory fault.

Release 1.03  - 30th July 2008
------------------------------
novopaired
	1. Correction to alignment quality calculation for paired end. In V1.01 some high quality alignments were
	   were getting assigned low qualities. This was evident by high number of true positives with Q0 using
	   simulated reads.

Release 1.02  - 30th July 2008
------------------------------
novopaired & novoalign
	1. Fix an out-of-bounds access. Problem could occurred if a read had no alignments with score better than threshold
                and one or more with score no more than 5 points worse than the threshold.

Release 1.01  - 22nd July 2008
------------------------------
novopaired
	1. Improvement in alignment quality scores and in detection of reads with multiple alignment locations.
	2. Performance improvements for large genome swith high repeat content.
	3. *.fq files are treated as Sanger FASTQ format
	4. Report order is now always the alignment for the read from first file followed by the alignment 
                for the read from the second file of a pair.

Release 1.0  - 8th July 2008
------------------------------
Novoalign & novopaired
	1. Native report format changed. Column order, space delimiter for mismatches &
           addition of pair alignment location. Please refer to the manual.
	2. Added an option for 5' trimming of unaligned reads. Reads are trimmed 
           until they align or until they are too short to produce a quality
           alignment.
	3. Allow prb qualities greater than 30 when converting to Sanger fastq format.
	4. Fix a rare divide by zero in quality calculation routines.
	5. novopaired formatting of prb&seq read identifier is now consistent. See 0.21 changes.


Release 0.22  - 30th June 2008
------------------------------
novoalign
	1. RNA  mode defaults to report All alignments to a read but no longer
           forces this mode.
	2. Calculation of default alignment threshold has been slightly relaxed.

Release 0.21  - 26th June 2008
------------------------------
novoalign & novopaired
	1. Added Sanger FASTQ format quality values to the report line.
	2. Changes to read header formatting when using prb files as input. If 
           prb & seq are used then tabs are changed to _ and the read sequence from the 
           seq record is not included in header.

Release 0.20
-----------------
novopaired
	1. Performance improvements. Fixes a problem where reads with only 5 or 6
           good quality bases were being used for alignment. These reads would never 
           align to a specific location.

Release 0.19
-----------------
Novoalign
	1. Added option to strip adapter sequences from the read
	2. Added an miRNA alignment mode.

Change History
Release 0.18
-----------------
Novopaired & Novoalign
	1. Correction to alignment quality when printing all alignments to a read.
	2. Minor tweak of alignment quality calculations.
	3. Added automatic threshold mode. If no threshold is specified a threshold will be
	    automatically calculated for each read.
	4. For novopaired, allow two read files as input as an alternative to two folders.

Release 0.17
-----------------
Novoalign & paired
	1. Fix problem where reads with no alignment(NM) were reported as low quality (QC)
	2. Changed Native report format to include print of read sequence.
	3. Add option to limit number of repeats reported using posterior probability
	4. Checked and validated input file formats fasta, qual, fastq, solexa fastq and prb.


Release Beta 0.16 1st June 2008
-------------------------------------------
novoindex
	1. Fix problem with auto k&s selection on genome >2Gbp.
	2. Remove memory map to avoid paging problem on SUSE 10.3
novopaired
	1. Fixed problem on alignment for genome >2Gbp.
	2. Fix problem with fasta format file recognition.

Release Beta 0.15
-----------------
novoindex
	No change
novopaired
	1.	Add quality score and related options for reporting.
	2.	Changed Blast format to Pairwise
	3.	Minor changes to Native report format. Now shows which the side each alignment
                originated from.
	4. 	Report alignment score as a positive rather than negative number so consistent
                with Phred.
	5.	Fix for alignment across boundary of two reference sequences.

novoalign
	1.	Minor changes to Native report format
	2.	Changed Blast format to Pairwise
	3. 	Report alignment score as a positive rather than negative number so consistent 
                with Phred.
	4.	Fix for alignment across boundary of two reference sequences.

Beta 0.14
---------
novoindex
	1. Added attributes to index for software version used to construct the index
	2. Added controls when opening index to verify this version of software will work with 
           the index.
novopaired
	1. No change.
novoalign
	1. Added alignment qualities based on posterior alignment probabilities. Field is added 
           to Native and Blast format output.
	2. Added additional program options for quality limit on reporting
	3. Modified repeat detection to be based on quality of first alignment plus number of
           alignments.
	4. Added options for reporting of repeats (refer manual)

Beta 0.13
Novoindex
	1. Change sequence header processing so that headers are truncated at first space
           rather than at 23 letters.
Novopaired
	1. Fix alignment location in Eland & Native format.
Novoalign
	1. Fix problem where a read exactly at the end of a sequence could not be aligned with gaps.

Beta 0.11
Novopaired
	1. First release

Beta 0.8
Novoalign
	1. Added Sanger FASTQ support
	2. Fix for case of DOS format query and quality files.

Beta 0.7
Novoindex
	1.	Added -m option. If set then lowercase sequence is not indexed.

Beta 0.6
Novoindex
	1.	Fix to Blast like output format for gapped alignments
	2.	Limits and controls on gap penalties so that maximum gap length never exceeds 7pb
	3.	Report usage if called with insufficient arguments.

Beta 0.5
novoindex
	1.	Changed command options so k & s are introduced using -k and -s
	2.	k & s are now optional and system will choose reasonable values based
		on the available RAM and the length of the reference sequences.
novoalign - No change

Beta 0.4
novoindex
	1.	Constrain k & s to reasonable values.
novoalign
	1.	Fix problem related to read quality check.

Beta 3.0

novoindex – No change
novoalign
	1.	Added report line to show interpretation of input file formats.
	2.	Added an option -n 99 that can limit length of the reads used for alignments. 
                This enables easier comparison of results with Eland

Beta 0.2
novoindex
	1.	Reduced memory requirements during index construction
	2.	Fixed an indexing issue for step size 3 and above

novoalign
	1.	Fixes to file format detection for prb with no seq file.
	2.	Added some missing newlines in reporting.

Beta 0.0
	Initial Release for evaluation purposes.

This software should be compatible with most modern 64-bit Linux systems running on X64 CPUs.
It was built and tested on 64 bit Open Suse 10.3 on an AMD Athlon X64
