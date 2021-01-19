def extract_HLA_reads(bam,name,outdir,runshFH,chrNaming):
    # get bam header
    cmd = "samtools view -H %s >%s/%s.hla.sam" % (bam,outdir,name)
    runshFH.write(cmd+'\n')

    if chrNaming == 'no_prefix':
        '''
        Note:
            1) this region is only for hg19
            2) please do not contain chr6_*_hap* in your hg19 fasta file

            #1.chr6_apd_hap1
            #2.chr6_cox_hap2
            #3.chr6_dbb_hap3
            #4.chr6_mann_hap4
            #5.chr6_mcf_hap5
            #6.chr6_qbl_hap6
            #7.chr6_ssto_hap7
        '''
        #region1 = '6:29909037-29913661'
        #region2 = '6:31321649-31324964'
        #region3 = '6:31236526-31239869'
        region1 = '6:29909000-29914000'
        region2 = '6:31321000-31326000'
        region3 = '6:31236000-31241000'
    else:
        #region1 = 'chr6:29909037-29913661'
        #region2 = 'chr6:31321649-31324964'
        #region3 = 'chr6:31236526-31239869'
        region1 = 'chr6:29909000-29914000'
        region2 = 'chr6:31321000-31326000'
        region3 = 'chr6:31236000-31241000'

    cmd = "samtools view %s %s >>%s/%s.hla.sam" % (bam,region1,outdir,name)
    runshFH.write(cmd+'\n')

    cmd = "samtools view %s %s >>%s/%s.hla.sam" % (bam,region2,outdir,name)
    runshFH.write(cmd+'\n')

    cmd = "samtools view %s %s >>%s/%s.hla.sam" % (bam,region3,outdir,name)
    runshFH.write(cmd+'\n')

    # sam to bam
    hla_sam = "%s/%s.hla.sam" % (outdir,name)
    hla_bam = "%s/%s.hla.bam" % (outdir,name)
    cmd = "samtools view -b -o %s %s" % (hla_bam,hla_sam)
    runshFH.write(cmd+'\n')

    # sort by read name
    bam_sort_by_name = "%s/%s.sort_by_name.bam" % (outdir,name)
    cmd = "samtools sort -n %s -o %s" % (hla_bam,bam_sort_by_name)
    runshFH.write(cmd+'\n')

    # bedtools bamtofastq
    fq1 = "%s/%s.chr6region.1.fastq" % (outdir,name)
    fq2 = "%s/%s.chr6region.2.fastq" % (outdir,name)
    cmd = "bedtools bamtofastq -i %s -fq %s -fq2 %s" % (bam_sort_by_name,fq1,fq2)
    runshFH.write(cmd+'\n')

    # Illegal Mate State
    # https://www.biostars.org/p/59521/

    #cmd = "%s -jar %s/SamToFastq.jar I=%s F=%s F2=%s VALIDATION_STRINGENCY=SILENT" % (jre,gatk_dir,hla_sam,fq1,fq2)










    

