/*
 *  19th Aug 2008, Modified by Colin Hercus to include conversion of novoalign and novopaired reports to maq format.
 *  Released under GPL
 */
#include <string.h>
#include <ctype.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <unistd.h>
#include "maqmap.h"
#include "algo.hh"
#include "stdhash.hh"
#include "main.h"
#include "limits.h"

static int DEFAULT_QUAL = 25;
static int log_n[256];
static bool chrlist = true;
static char **name_conv;
static int *cnt_chr;
static int n_hdr = 0;
bool sw_on = false;
bool ref_report = false;
#define QLIM 60

static bool find_n_add(hash_map_char<bit32_t> *hash, char *hdr, bit32_t *i) {
    if(hash->find(hdr, i))
        return true;
    if(chrlist)
        return false;
    hash->insert(hdr, n_hdr);
    if ((n_hdr & 0xff) == 0) {
        name_conv = (char**)realloc(name_conv, sizeof(char*) * (n_hdr + 0x100));
        cnt_chr = (int *)realloc(cnt_chr, sizeof(int) * (n_hdr + 0x100));
        memset( (void *) (cnt_chr + n_hdr) , 0, 0x100 * sizeof(int));
    }
    *i = n_hdr;
    name_conv[n_hdr++] = strdup(hdr);
    return true;
}
static hash_map_char<bit32_t> *read_list(FILE *fp)
{
	hash_map_char<bit32_t> *hash = new hash_map_char<bit32_t>;
    if(fp == NULL)
        return hash;
	char str[1024];
	while (fscanf(fp, "%s", str) != 0) {
		if (feof(fp)) break;
		char *p = (str[0] == '>')? str+1 : str;
		int c;
		bit32_t x;
		assert(!hash->find(p, &x));
		hash->insert(p, n_hdr);
		fscanf(fp, "%s", str);
        p = (str[0] == '>')? str+1 : str;
        if ((n_hdr & 0xff) == 0) {
			name_conv = (char**)realloc(name_conv, sizeof(char*) * (n_hdr + 0x100));
			cnt_chr = (int *)realloc(cnt_chr, sizeof(int) * (n_hdr + 0x100));
            memset( (void *) (cnt_chr + n_hdr) , 0, 0x100 * sizeof(int));
        }
		name_conv[n_hdr++] = strdup(p);
		while ((c = fgetc(fp)) != EOF && c != '\n');
	}
	return hash;
}
static inline int cal_map_qual(int default_qual, bit32_t *count)
{
	if (count[0] == 1) {
		if (count[1] == 0 && count[2] == 0) return 3 * default_qual;
		if (count[1] == 0) return 2 * default_qual - log_n[count[2]];
		return default_qual - log_n[count[1]];
	}
	if (count[1] == 1) {
		if (count[2] == 0) return 2 * default_qual;
		return default_qual - log_n[count[2]];
	}
	if (count[2] == 1) return default_qual - 3;
	return default_qual;
}
static inline int operator < (const maqmap1_t &a, const maqmap1_t &b)
{
	return (a.seqid < b.seqid) || (a.seqid == b.seqid && a.pos < b.pos);
}

/* novo2maq 
  name: read name
  size: the length of the read
  seq: read sequence (see also below)
  seq[MAX_READLEN-1]: single end mapping quality (equals to map_qual if not paired)
  map_qual: the final mapping quality
  alt_qual: the lower quality of the two ends (equals to map_qual if not paired)
  flag: status of the pair
  dist: offset of the mate (zero if not paired)
  info1: mismatches in the 24bp (higher 4 bits) and mismatches (lower 4 bits)
  info2: sum of errors of the best hit
  c[2]: count of all 0- and 1-mismatch hits on the reference
 */
//Count mismatches
static int novomismatches(char *mm) {
    int i = 0;
    while(*mm != '\0') {
        if(*mm == '>')
            i++;
        mm++;
    }
    return i;    
}
// Recover novo indels while making sure to only pick up indels smith-waterman would find.
// Needleman-Wunsch used in novo may report indels too close to the ends of the reads.
static bool gapped(char * mm, char strand, int rdlen, int &readposn, int &indellen) {
    readposn = 0;
    indellen = 0;
    // Bisulphite - Skip mode flag
    if(strncmp(mm,"CT",2) == 0 || strncmp(mm, "GA",2) == 0) {
        mm += 2;
        if(*mm == ' ')
            mm++;
    }
    if(*mm == '\0')
        return false;

    int tposn = 1, qposn = 1, posn = 1;
    int match = 15, mismatch = 30, gapopen = 40, gapextend = 10;
    int maxsc = 0, maxposn = 0;
    bool gapped = false;
    int sc[256], mv[256], qp[256], tp[256];
    sc[0] = qp[0] = tp[0] = '\0';
    mv[0] = '|';
    while(*mm != '\0') {
        int i = 0;
        while(isdigit(*mm))
            i = i* 10 + (*mm++ - '0');
        while(tposn < i) {
            sc[posn] = sc[posn-1] + match;
            mv[posn] = '|';
            qp[posn] = qposn;
            tp[posn] = tposn;
            tposn++;
            qposn++;
            posn++;
        }
        if(sc[tposn-1] > maxsc) {
            maxsc = sc[tposn - 1];
            maxposn = tposn - 1;
        }
        switch(*mm) {
            case '+':
                gapped = true;
                if(mv[posn-1] == '+')
                    sc[posn] = sc[posn-1] - gapextend;
                else
                    sc[posn] = sc[posn-1] - gapopen;
                mm += 2;
                do {
                    if(sc[posn] < 0)
                        sc[posn] = 0;
                    mv[posn] = '+';
                    qp[posn] = qposn;
                    tp[posn] = tposn;
                    qposn++;
                    posn++;
                    char c = *mm++;
                    if(c == ' ' || c == '\0') break;
                    sc[posn] = sc[posn-1] - gapextend;
                } while(true);
                break;
            case '-':
                gapped = true;
                if(mv[posn-1] == '-')
                    sc[posn] = sc[posn-1] - gapextend;
                else
                    sc[posn] = sc[posn-1] - gapopen;
                if(sc[posn] < 0)
                    sc[posn] = 0;
                mv[posn] = '-';
                qp[posn] = qposn;
                tp[posn] = tposn;
                tposn++;
                posn++;
                mm += 3;
                break;
            default:
                sc[posn] = sc[posn-1] - mismatch;
                if(sc[posn] < 0)
                    sc[posn] = 0;
                mv[posn] = '|';
                qp[posn] = qposn;
                tp[posn] = tposn;
                tposn++;
                qposn++;
                posn++;
                mm+= 4;
        }
    }
    if(!gapped)
        return false;
    while(qposn <= rdlen) {
        sc[posn] = sc[posn-1] + match;
        mv[posn] = '|';
        qp[posn] = qposn;
        tp[posn] = tposn;
        tposn++;
        qposn++;
        posn++;
    }
    if(!sw_on) {
        for(readposn = 0; readposn < posn && mv[readposn] == '|'; readposn++);
        for(indellen = 0; readposn + indellen < posn && mv[readposn + indellen] != '|'; indellen ++);
        if(mv[readposn] == '-')
            indellen = -indellen;
        readposn -= 1;
        return true;
    }
    if(sc[posn-1] > maxsc) {
        maxsc = sc[posn - 1];
        maxposn = posn - 1;
    }
    if(maxsc <= 0)
        return false;
    int j = maxposn;
    while(sc[j] > 0 && mv[j] == '|')
        j--;
    if(mv[j] == '|') {
        readposn = tp[j] - qp[j];
        return false;
    }
    indellen = 1;
    while(mv[j-1] == mv[j]) {
        indellen++;
        j--;
    }
    readposn = qp[j] - 1;
    if(mv[j] == '-')
        indellen = -indellen;
    return true;
}

#define min(x,y) ((x)<(y)? (x):(y))

static void novo2maq_core(FILE *fp_list, FILE *fp_novo, gzFile fp)
{
    int miRNAmode = 0;           // In miRNA mode there is an extra score value before alignment location fields.
	hash_map_char<bit32_t> *hash;
	// initialize maqmap_t
	maqmap_t *mm = maq_new_maqmap();
	int max = 0, i, l;
	hash = read_list(fp_list);
	// initialize log_n
	log_n[0] = -1;
	for (i = 1; i != 256; ++i)
		log_n[i] = (int)(3.434 * log((float)i) + 0.5);
	// read the file
    enum {n_readid, n_side, n_seq, n_qual, n_status, n_score, n_quality, N_chr, N_offset, N_strand, N_chr2, N_offset2, N_strand2, N_poly, NLAST};
    #define n_chr (N_chr + miRNAmode)
    #define n_offset (N_offset + miRNAmode)
    #define n_strand (N_strand + miRNAmode)
    #define n_chr2 (N_chr2 + miRNAmode)
    #define n_offset2 (N_offset2 + miRNAmode)
    #define n_strand2 (N_strand2 + miRNAmode)
    #define n_poly (N_poly + miRNAmode)
    char * novo[NLAST+1];
    char * novo2[NLAST+1];
	char str[4096], str2[4096];
	bit8_t tmp_seq[MAX_READLEN];
	bit32_t tmp[4];
	while (fgets( str, 4096,fp_novo ) != NULL) {
		if (feof(fp_novo)) break;
        if(str[0] == '#' || str[0] == '\0') {
            char *mopt;
            // Check for -m option
            if(strncmp(str, "#  novoalign", 11) == 0 && (mopt = strstr(str, " -m")) != NULL)
                miRNAmode = 1;
            continue;
        }
		if (mm->n_mapped_reads + 1 >= (bit64_t)max)      // Make sure room for a pair.
			mm->mapped_reads = (maqmap1_t*)realloc(mm->mapped_reads, sizeof(maqmap1_t) * (max = max == 0 ? 100000 : max * 1.5));
		maqmap1_t *m1 = mm->mapped_reads + mm->n_mapped_reads;
        char *np = str;
        char mtstr[1] = "";
        for(i = 0; i <= NLAST; i++ ) {
            novo[i] = np;
            while(*np != '\t' && *np != '\n' && *np != '\0')
                np++;
            if(*np == '\n' || *np == '\0') break;
            *np++ = '\0';
        }
        *np = '\0';
        for(i++;i <= NLAST; i++ )
            novo[i] = mtstr;
        
        // Read pair
        bool paired =false;
        if(novo[n_side][0] != 'S') {
            paired = true;
            fgets( str2, 4096, fp_novo );
            np = str2;
            for(i = 0; i <= NLAST; i++ ) {
                novo2[i] = np;
                while(*np != '\t' && *np != '\n')
                    np++;
                if(*np == '\n') break;
                *np++ = '\0';
            }
            *np = '\0';
            for(i++;i <= NLAST; i++ ) 
                novo2[i] = mtstr;
        }
        else
            miRNAmode = novo[N_offset][0] == '>' ? 1 : 0;
        // read flag
		if (novo[n_offset][0] != '\0' || (paired && novo2[n_offset][0] != '\0')) { // there is an alignment
            // set name
            strncpy(m1->name, novo[n_readid]+1, MAX_NAMELEN-1);
            m1->name[MAX_NAMELEN - 1] = 0;
            // set seq
            m1->size = l = strlen(novo[n_seq]);
            if(novo[n_qual][0] == '\0')
                for (i = 0; i != l; ++i) {
                    tmp[0] = nst_nt4_table[(int)novo[n_seq][i]];
                    m1->seq[i] = (tmp[0] > 3)? 0 : (tmp[0]<<6 | DEFAULT_QUAL);
                }
            else
                for (i = 0; i != l; ++i) {
                    tmp[0] = nst_nt4_table[(int)novo[n_seq][i]];
                    m1->seq[i] = (tmp[0] > 3)? 0 : (tmp[0]<<6 | ((int)novo[n_qual][i] - 33));
                }
            int n_mis = novomismatches(novo[n_poly]);
            bool aligned = novo[n_chr][0] != '\0';
            if(aligned) {
                aligned = find_n_add(hash, novo[n_chr]+1, &m1->seqid);
                if(aligned) cnt_chr[m1->seqid]++;
            }
            if (novo[n_strand][0] == 'R') {
                for (i = m1->size - 1; i >= 0; --i)
                    tmp_seq[m1->size-i-1] = (m1->seq[i] == 0)? 0 : (0xc0 - (m1->seq[i]&0xc0)) | (m1->seq[i]&0x3f);
                memcpy(m1->seq, tmp_seq, m1->size);
            }

            m1->c[0] = 0;   // 0 mismatches alignments
            m1->c[1] = 0;   // 1 mismatch alignments
            if (n_mis < 2) m1->c[n_mis] = 1;
            m1->flag = 0;   //paired flag  FR, correct pair, 16 good pair, 32 diff chr, 64 one read of pair not mapped, 128 this read not mapped and other read is.
            m1->dist = 0;   // Outer coordinates
            m1->pos = (atoi(novo[n_offset])-1)<<1 | (novo[n_strand][0] == 'F'? 0 : 1);
            m1->info1 = n_mis<<4 | n_mis;
            m1->info2 = atoi(novo[n_score])/3;
            m1->map_qual = m1->seq[MAX_READLEN-1] = min(QLIM, atoi(novo[n_quality]));
            m1->alt_qual = m1->map_qual;  // Lesser of two reads
            if(paired) {
                int dist = 0;
                int flag = 64;
                if (novo[n_offset][0] == '\0') {
                    m1->flag = 192;
                    m1->info2 = 0;
                    m1->map_qual = m1->seq[MAX_READLEN-1] = 0;
                    m1->alt_qual = min(QLIM, atoi(novo2[n_quality]));
                } else if (novo2[n_offset][0] == '\0') {
                    m1->flag = 64;
                    flag = 192;
                    m1->alt_qual = min((int)m1->map_qual, atoi(novo2[n_quality]));
                } else {
                    if((novo[n_chr2][0] == '\0' || strcmp(novo[n_chr2],".") == 0)
                            && strcmp(novo[n_quality], novo2[n_quality]) == 0
                            && (dist = abs( atoi(novo[n_offset]) - atoi(novo[n_offset2]))) < 2000) {
                        if(novo[n_strand][0] == 'R')
                            dist = -(dist + strlen(novo[n_seq]));
                        else
                            dist += strlen(novo2[n_seq]);
                        m1->flag = flag = 2 + 16;
                        m1->dist = dist;
                    } else if(novo[n_chr2][0] != '\0' && strcmp(novo[n_chr2],".") != 0)
                        m1->flag = flag = 32;
                    else if(novo[n_strand][0] == 'R' && novo2[n_strand][0] == 'R')
                        m1->flag = flag = 8;
                    else if(novo[n_strand][0] == 'F' && novo2[n_strand][0] == 'F')
                        m1->flag = flag = 1;
                    else
                        m1->flag = flag = 8; 
                    int gapposn, gaplen;
                    if(gapped(novo[n_poly], novo[n_strand][0], m1->size, gapposn, gaplen)) {
                        //gapped(novo[n_poly], novo[n_strand][0], m1->size, gapposn, gaplen);
                        m1->flag = 130;
                        m1->map_qual = gapposn;
                        m1->seq[MAX_READLEN-1] = gaplen;
                    } else
                        m1->pos = (atoi(novo[n_offset]) + gapposn -1)<<1 | (novo[n_strand][0] == 'F'? 0 : 1);

                }
                if(aligned && m1->flag != 192  && (novo[n_status][0] != 'R' || m1->flag != 130))
                    ++mm->n_mapped_reads;
                if(flag == 192)
                    continue;
                m1 = mm->mapped_reads + mm->n_mapped_reads;
                            // set name
                strncpy(m1->name, novo2[n_readid]+1, MAX_NAMELEN-1);
                m1->name[MAX_NAMELEN - 1] = 0;
                // set seq
                m1->size = l = strlen(novo2[n_seq]);
                if(novo2[n_qual][0] == '\0')
                    for (i = 0; i != l; ++i) {
                        tmp[0] = nst_nt4_table[(int)novo2[n_seq][i]];
                        m1->seq[i] = (tmp[0] > 3)? 0 : (tmp[0]<<6 | DEFAULT_QUAL);
                    }
                else
                    for (i = 0; i != l; ++i) {
                        tmp[0] = nst_nt4_table[(int)novo2[n_seq][i]];
                        m1->seq[i] = (tmp[0] > 3)? 0 : (tmp[0]<<6 | ((int)novo2[n_qual][i] - 33));
                    }
                int n_mis = novomismatches(novo2[n_poly]);
                
                if (novo2[n_chr] == '\0' || !find_n_add(hash, novo2[n_chr]+1, &m1->seqid))
                    continue;                                        

            
                cnt_chr[m1->seqid]++;
                aligned = true; 
                if (novo2[n_strand][0] == 'R') {
                    for (i = m1->size - 1; i >= 0; --i)
                        tmp_seq[m1->size-i-1] = (m1->seq[i] == 0)? 0 : (0xc0 - (m1->seq[i]&0xc0)) | (m1->seq[i]&0x3f);
                    memcpy(m1->seq, tmp_seq, m1->size);
                }

                m1->c[0] = 0;   // 0 mismatches alignments
                m1->c[1] = 0;   // 1 mismatch alignments
                if (n_mis < 2) m1->c[n_mis] = 1;
                m1->flag = flag;   // paired flag  FR, correct pair, 16 good pair, 32 diff chr, 64 one read of pair not mapped, 128 this read not mapped and other read is.
                m1->dist = -dist;   // Outer coordinates
                m1->pos = (atoi(novo2[n_offset])-1)<<1 | (novo2[n_strand][0] == 'F'? 0 : 1);
                m1->info1 = n_mis<<4 | n_mis;
                m1->info2 = atoi(novo2[n_score])/3;
                m1->map_qual = m1->seq[MAX_READLEN-1] = min(QLIM, atoi(novo2[n_quality]));
                if (novo[n_offset][0] != '\0')
                    m1->alt_qual = min(QLIM,min((int)m1->map_qual, atoi(novo[n_quality]))); 
                else
                    m1->alt_qual = m1->map_qual;
                int gapposn, gaplen;
                if(gapped(novo2[n_poly], novo2[n_strand][0], m1->size, gapposn, gaplen)) {
                    //gapped(novo2[n_poly], novo2[n_strand][0], m1->size, gapposn, gaplen);
                    m1->flag = 130;
                    m1->map_qual = gapposn;
                    m1->seq[MAX_READLEN-1] = gaplen;
                } else
                    m1->pos = (atoi(novo2[n_offset]) + gapposn -1)<<1 | (novo2[n_strand][0] == 'F'? 0 : 1);
                if(aligned && (novo[n_status][0] != 'R' || m1->flag != 130))
                    ++mm->n_mapped_reads;
            } else {
                int gapposn, gaplen;
                if(gapped(novo[n_poly], novo[n_strand][0], m1->size, gapposn, gaplen)) {
                    m1->flag = 130;
                    m1->alt_qual = m1->map_qual;
                    m1->map_qual = gapposn;
                    m1->seq[MAX_READLEN-1] = gaplen;
                } else
                    m1->pos = (atoi(novo[n_offset]) + gapposn -1)<<1 | (novo[n_strand][0] == 'F'? 0 : 1);
                if(aligned && (novo[n_status][0] != 'R' || m1->flag != 130))
                    ++mm->n_mapped_reads;
            }
		}
	}
	mm->n_ref = hash->size();
	mm->ref_name = (char**)malloc(sizeof(char*) * mm->n_ref);
	for (i = 0; i != int(hash->size()); ++i)
		mm->ref_name[i] = strdup(name_conv[i]);
	algo_sort(mm->n_mapped_reads, mm->mapped_reads);
	maqmap_write_header(fp, mm);
    unsigned long n_write = UINT_MAX /  sizeof(maqmap1_t);
    unsigned long n = 0;
    while((mm->n_mapped_reads - n) > n_write) {
        gzwrite(fp, mm->mapped_reads + n, sizeof(maqmap1_t) * n_write);
        n += n_write;
    }
	gzwrite(fp, mm->mapped_reads + n, sizeof(maqmap1_t) * (mm->n_mapped_reads - n));
	// free
	maq_delete_maqmap(mm);
    if(ref_report) {
        printf("Count\tHeader\n");
        for(int i = 0; i < n_hdr; i++)
            printf("%d\t%s\n", cnt_chr[i], name_conv[i]);
    }
	for (i = 0; i != int(hash->size()); ++i) free(name_conv[i]);
	free(name_conv);
    free(cnt_chr);
	delete hash;
}

static int novo2maqusage()
{
	fprintf(stderr, "Usage: novo2maq [-s off] [-r] <out.map> <in.list> <in.novo>\n");
	fprintf(stderr, "Options:\n");
	fprintf(stderr, "                -s on|off Turns on or off the Smith-Waterman check of indels. When only indels that"
                    "/n                          would be part of optimum local alignment are converted as MAQ indels. Default is off "
                    "/n                          which will report all indels found by Novoalign.\n");
	fprintf(stderr, "                -r        Produces short report on number of alignments per reference sequence.\n"
                    " Where:\n"
                    "                out.map   is file name for the output MAQ map file.\n"
                    "                in.list   is a list of reference sequence headers to be selected. This file servers to specify\n"
                    "                          reference sequences to be selected for conversion to the MAQ map format and also allows\n"
                    "                          translation of the header sequence. Each line in the file has format:\n"
                    "                              <refheader> <tab> <replacement header>\n"
                    "                          Any reference sequences not listed in this file will not be converted to MAQ map file.\n"
                    "                          Use '-' rather than a filename to specify that all sequences are to be selected and converted.\n"
                    "                in.novo   is the Novoalign report file to be converted. Use '-' to specify that the report is to be read from stdin\n");
    
	return 1;
}

int maq_novo2maq(int argc, char *argv[])
{
	FILE *fp_novo, *fp_list;
	gzFile fp_map;
	int c;
	while ((c = getopt(argc, argv, "rs:")) >= 0) {
		switch (c) {
            case 's': 
                if(strcmp(optarg,"off") == 0)
                    sw_on = false; 
                if(strcmp(optarg,"on") == 0)
                    sw_on = true; 
                break;
            case 'r':
                ref_report = true;
                break;
		}
	}
	if (optind + 2 >= argc) return novo2maqusage();
	fp_novo = (strcmp(argv[optind+2], "-") == 0)? stdin : fopen(argv[optind+2], "r");
    if(strcmp(argv[optind+1], "-") == 0) {
        fp_list = NULL;
        chrlist = false;
    } else
        fp_list = fopen(argv[optind+1], "r");
	fp_map = gzopen(argv[optind], "w");
	assert(fp_novo && (!chrlist || fp_list) && fp_map);
	novo2maq_core(fp_list, fp_novo, fp_map);
	fclose(fp_novo); 
    if(fp_list) fclose(fp_list); 
    gzclose(fp_map);
	return 0;
}

