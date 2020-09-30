(C) 2008, 2009, 2010, 2011 NovoCraft Technologies Sdn Bhd

Freely downloaded versions of Novocraft programs are licensed for educational use
and for use in not-for-profit organisation for internal use only.

For-profit organisations without paid up licenses can only use this software for
the purpose of evaluating it.

For commercial licenses and support contracts contact sales@novocraft.com

Documentation for Novomethyl is online at  http://tinyurl.com/Novocraft-BiSeq


NOVOMETHYL  Change History
--------------------------

Release Novomethyl Beta.8.0 9th August 2011
---------------------------------------------
1.    Change calculation of conversion efficiency. It is now the % of converted Cytosines at non CpG sites.

Release Novomethyl Beta.7.0 28th March
---------------------------------------------
1.    Version 6.0 was giving excessive quality to base & methylation state calls. 
      Almost everything had quality of 150 due to max/min test being inverted. 
2.    In consensus format count of methylated and unmethylated cytosines was 
      always off the CT strand, it should be off the GA strand for base calls of G.

Release Novomethyl Beta.6.0 - 24th March 2011
---------------------------------------------
1.    Add a new report format ( -o CONSENSUS )that is similar to samtools 
      pileup format with extra fields for methylation context, quality and 
      percentage.
2.    Reviewed calculation of posterior probabilities with some corrections 
      that to improve accuracy of quality call. Significant difference is
      only for bases with very low depth of reads.
3.    In BED format make offsets 0 based.


Release Novomethyl Beta.5.0
---------------------------
1.    In previous release the option -m (probability of a Cytosine being 
      methylated) was interpreted as P not methylated. This is fixed.
2.    In the run summary the count of Cytosine's not reported because of 
      a SNP was inaccurate. This has now been fixed.
3.    Exit status is now 0 (zero) for normal completion.

Release Novomethyl Beta.4.0
---------------------------
1.    Add mode for reporting percentage methylation. This should be used for 
      heterogeneous samples.
2.    Corrections to calling of Heterozygous cytosines including adding 
      tag (Het) to context in the bed file. In earlier version heterozygous 
      cytosines were never reported.
3.    Reviewed setting of prior probabilities making some small adjustments 
      in particular a double heterozygous SNP (say C to A/T) gets prior Psnp^2.
      Other changes improve accuracy of calls when only 1 or 2 reads cover 
      a locus.
4.    Add a run report to stderr

Release Novomethyl Beta_3.0
---------------------------
      1.  Fix interpretation of mpileup base qualities when base string contains 
          a $. This had caused bases and qualities to become misaligned. Novomethyl
          should probably be rerun to correct results.

Release Novomethyl Beta_2.0
---------------------------
      Beta release of methylation status caller. Please refer to our wiki at 
      http://tinyurl.com/Novocraft-BiSeq

