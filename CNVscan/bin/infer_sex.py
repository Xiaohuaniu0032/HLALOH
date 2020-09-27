import sys
import argparse

class Options:
    def __init__(self):
        self.parser = argparse.ArgumentParser("Infer Sample Sex")
        self.parser.add_argument('-cov',help=".targetcoverage.cnn file",dest="cov")
        self.parser.add_argument('-o',help='outfile',dest='of')

        args = self.parser.parse_args()

        if args.cov:
            self.cov = args.cov
        else:
            self.parser.error(".targetcoverage.cnn not supplied")

        if args.of:
            self.of = args.of
        else:
            self.parser.error("please specify outfile name")



def chr_naming(cnn_file):
    cnn = open(cnn_file,'r')
    cnn.readline()
    first_line = cnn.readline()
    chrom = first_line.split('\t')[0]
    # if chr1 or 1
    if 'chr' in chrom:
        naming = "with_prefix"
    else:
        naming = "no_prefirx"

    return(naming)
    
        
def infer_sex(infile,outfile):
    auto_cov = []
    chrX_cov = []
    chrY_cov = []

    with open(infile,'r') as cov_H:
        next(cov_H)
        for line in cov_H:
            arr = line.rstrip().split('\t')
            cov = int(float(arr[-2]))
            
            chr_name = chr_naming(infile)
            if chr_name == "with_prefix":
                # chr1
                if arr[0] != "chrX" and arr[0] != "chrY":
                    auto_cov.append(cov)
                if arr[0] == "chrX":
                    chrX_cov.append(cov)
                if arr[0] == "chrY":
                    chrY_cov.append(cov)
            else:
                # 1
                if arr[0] != 'X' and arr[0] != 'Y':
                    auto_cov.append(cov)
                if arr[0] == 'X':
                    chrX_cov.append(cov)
                if arr[0] == 'Y':
                    chrY_cov.append(cov)


    '''
    1.if chrY exists, then infer sex by chrY, and do not consider chrX
    2.if chrY not exists, and chrX exists, then infer sex only by chrX
    3.if chrX and chrY all not exists, then set sex by NA, which means do not correct sex's chr log2ratio

    if chrY exists
        mean_cov_chrY / mean_cov_chrX >= 0.1 ,and mean_cov_chrY >= 100, then sex => Male; else sex => Female
    if chrY not exists and chrX exists:
        mean_cov_chrX / mean_cov_autochr >= 0.67, then sex => Female; else sex => Male
    if chrX and chrY all not exists:
        sex => NA

    '''
    sex = ''

    if len(chrX_cov) == 0:
        sex = 'NA'
        print("can not find chrX, will not infer sex")
        return(sex)

    if len(chrX_cov) == 0 and len(chrY_cov) == 0:
        sex = 'NA'
        print("can not find chrX and chrY, will not infer sex")
        return(sex)

    mean_cov_autochr = int(sum(auto_cov)/len(auto_cov))
    mean_cov_chrX = int(sum(chrX_cov)/len(chrX_cov))
    

    # infer sex
    if len(chrY_cov):
        print("chrY exists!")
        mean_cov_chrY = int(sum(chrY_cov)/len(chrY_cov))
        y_x_ratio = mean_cov_chrY / mean_cov_chrX
        if y_x_ratio >= 0.1 and mean_cov_chrY >= 100:
            sex = 'male'
            print("Y/X ratio is %s (>=0.1), sex is male" % y_x_ratio)
        else:
            sex = 'female'
            print("Y/X ratio is %s (< 0.1), sex is female" % y_x_ratio)
    else:
        # no chrY
        print("chrY not exists!")
        x_autochr_ratio = mean_cov_chrX / mean_cov_autochr
        if x_autochr_ratio >= 0.67:
            sex = 'female'
            print("X/autochr ratio is %s (>=0.67), sex is female" % x_autochr_ratio)
        else:
            sex = 'male'
            print("X/autochr ratio is %s (<0.67), sex is male" % x_autochr_ratio)


    # out cov info
    print("autochr mean cov is: %s" % mean_cov_autochr)
    print("chrX mean cov is: %s" % mean_cov_chrX)

    if len(chrY_cov):
        mean_cov_chrY = int(sum(chrY_cov)/len(chrY_cov))
        print("chrY mean cov is: %s" % mean_cov_chrY)

    print("infered sex is %s" % sex)

    # write to outfile
    of = open(outfile,'w')
    of.write("sex"+'\t'+sex+'\n')
    of.close()

    return(sex)

if __name__ == "__main__":
    options = Options()
    infer_sex(options.cov,options.of)





