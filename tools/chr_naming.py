def chrNaming(fasta):
    infile = open(fasta,'r')
    firstLine = infile.readline()
    infile.close()
    
    if '>chr' in firstLine:
        chrNaming = 'with_prefix'
    else:
        chrNaming = 'no_prefix'

    return(chrNaming)
