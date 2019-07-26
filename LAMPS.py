#!/usr/bin/env python
#
#   LAMPS - 2C-ChIP and 5C library processing (Python v2 or v3)
#   Steps:
#       1) maps paired-end sequences to custom BLAST database of defined regions
#       2) provides barcode/primer Quality Checks (QC)
#       3) outputs processed data in standardized format: bedGraph (2C-ChIP) and my5C matrix (5C)
#   See README for more information: https://github.com/BlanchetteLab/LAMPS
#   Author: Christopher JF Cameron
#

from __future__ import print_function

import argparse,glob,itertools,math,matplotlib,multiprocessing,os,re,shutil,subprocess,sys
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

from os import environ
from scipy.stats import spearmanr


#   check Python and sed versions being used
python_2 = True
if sys.version_info > (3, 0):
    python_2 = False
sed_cmd = "sed" if sys.platform in ["linux","linux2"] else "gsed"

#   prevent MatPlotlib view window from opening
plt.ioff()  

#   check if BLAST can be found in the 'PATH' variable
if not re.search(r'.+blast.+',os.environ["PATH"].lower()):
    print("Error - BLAST cannot be found in the 'PATH' environmental variable. Please ensure BLAST is installed")
    sys.exit(-2)

#   check if samtools can be found in the 'PATH' variable for BAM->FASTQ conversion
if not re.search(r'.+samtools.+',os.environ["PATH"].lower()):
    print("Warning - samtools cannot be found in the 'PATH' environmental variable. BAM input will lead to errors")

def get_range(start,stop,step):
    """returns range useing Python version specific function"""
    return xrange(start,stop,step) if python_2 else range(start,stop,step)
    
def build_FASTAs(filepath,file_dict,barcodes=None):
    """parses provided primer file and creates FASTA files for BLASTdb generation"""
    print("\tParsing primer file ... ",end='',file=sys.stderr)
    #   create look-up table of omitted primers
    omitted_primers = set([])
    for key in file_dict.keys():
        for barcode in file_dict[key].keys():
            omitted_primers |= file_dict[key][barcode]["omitted_primers"]
    lookup_dict = {key:None for key in omitted_primers}
    
    #   parse primer file
    primer_dict = {}
    fwd_primers,rvs_primers = [],[]
    min_fwd_length,min_rvs_length = None,None
    with open(filepath,'rU' if python_2 else 'rt') as f:
        for line in f:
            name,num,strand,seq,assembly,chrom,start,end = line.rstrip().split('\t')
            chrom = chrom.lower().replace("chr",'')
            label = ''.join([strand,num,'@',name])
            seq_length = len(seq)
            try:
                primer_dict[label]
                print(''.join(["Error - multiple primers with the same label '",label,"' found in primer file"]))
                sys.exit(-2)
            except KeyError: primer_dict[label] = {"assembly":assembly,"chrom":chrom,"start":int(start),"seq":seq,"length":seq_length}
            try: 
                lookup_dict[num]
                lookup_dict[num] = label
            except KeyError: pass   #   primer isn't omitted
            if strand == 'F':
                fwd_primers.append(label)
                min_fwd_length = seq_length if min_fwd_length == None or seq_length < min_fwd_length else min_fwd_length
            else:
                rvs_primers.append(label)
                min_rvs_length = seq_length if min_rvs_length == None or seq_length < min_rvs_length else min_rvs_length     
    print("done",file=sys.stderr)
    
    print("\tProcessing primers ... ",end='',file=sys.stderr)
    #   sort primers by genomic locus
    fwd_primers = [val[0] for val in sorted([[primer,primer_dict[primer]["chrom"],primer_dict[primer]["start"]] for primer in fwd_primers],key=lambda x:(x[1],x[2]))]
    rvs_primers = [val[0] for val in sorted([[primer,primer_dict[primer]["chrom"],primer_dict[primer]["start"]] for primer in rvs_primers],key=lambda x:(x[1],x[2]))]
    
    #   switch omitted primer names to indices of sorted primers
    primers = fwd_primers+rvs_primers
    for key in file_dict.keys():
        for barcode in file_dict[key].keys():
            file_dict[key][barcode]["omitted_primers"] = [primers.index(lookup_dict[omitted_primer]) for omitted_primer in file_dict[key][barcode]["omitted_primers"]]
    
    #   create similarity matrix for bit scores of primer sequences
    n = len(primers)
    matrix = np.zeros((n,n),dtype=float)
    word_size = min(min_fwd_length,min_rvs_length)
    #   iterate over possible primer pairs: F-F,F-R,R-F,R-R
    tmp_dir = create_directory(os.path.join(primer_dir,"tmp"))
    filepath = os.path.join(os.path.join(tmp_dir,"tmp.fasta"))
    with open(filepath,'wt') as f:
        for primer in primers:
            f.write(''.join(['>',primer,'\n',primer_dict[primer]["seq"],'\n']))
    output = subprocess.check_output(''.join(["blastn -word_size ",str(word_size)," -max_target_seqs ",str(n)," -outfmt \"6 qseqid sseqid bitscore\" -query ",filepath," -subject ",filepath," -dust no"]),shell=True).decode("utf-8").rstrip().split('\n')
    warned = False
    #   parse blastn output
    for i,line in enumerate(output):
        primer_A,primer_B,score = line.split()
        score = float(score)
        if not primer_A == primer_B and score > 0.0 and not warned:
            print("Warning - two or more primers were found to be similar. Please review QC reports. ",end='',file=sys.stderr)
            warned = True
        #insert score into matrix
        index_A,index_B = primers.index(primer_A),primers.index(primer_B)
        matrix[index_A][index_B] = score
        #   R-F and F-R primer pairs are symmetric within similarity matrix
        if (primer_A.startswith('F') and primer_B.startswith('R')) or (primer_A.startswith('R') and primer_B.startswith('F')):
            matrix[index_B][index_A] = score
    #   remove temporary directory
    shutil.rmtree(tmp_dir)
    #   write matrix to storage
    np.savetxt(os.path.join(primer_dir,"primer_similarities.matrix.tsv"),np.hstack((np.reshape(primers,(n,1)),matrix.astype(str))),fmt="%s",delimiter="\t",header="\t".join([""]+primers),comments='')
    
    #   create heatmap of bitwise scores
    plot_heatmap(os.path.join(primer_dir,"primer_similarities.heatmap.png"),np.log(matrix+1.0),primers)
    print("done",file=sys.stderr)

    print("\tWriting paired-primer FASTA to storage ... ",end='',file=sys.stderr)
    primers = fwd_primers+rvs_primers if len(barcodes) == 0 else [''.join([barcode,'-',primer]) for barcode in barcodes for primer in fwd_primers]+rvs_primers
    #   create FASTA of paired-primers including potential T3 and barcode sequences
    with open(os.path.join(primer_dir,"primer_pairs.custom_db.fasta"),'wt') as o:
        for primer_A,primer_B in itertools.product(primers,repeat=2):
            #   split barcode from primer ID if present
            try: barcode_A,primer_A = primer_A.split('-')
            except ValueError: barcode_A,primer_A = '',primer_A
            try: barcode_B,primer_B = primer_B.split('-')
            except ValueError: barcode_B,primer_B = '',primer_B
            barcode_B = ''  #   barcode can never ligate to a primer
                
            seq_A_id = '-'.join([barcode_A,primer_A]) if not barcode_A == '' else primer_A
            seq_B_id = '-'.join([barcode_B,primer_B]) if not barcode_B == '' else primer_B
            seq_A = ''.join([barcode_A,primer_dict[primer_A]["seq"]])
            seq_B = ''.join([barcode_B,primer_dict[primer_B]["seq"]])

            o.write(''.join(['>',seq_A_id,'|',seq_B_id,'\n',seq_A,seq_B,'\n']))
    
    #   create FASTA of individual primers for 'unmappables'
    with open(os.path.join(primer_dir,"short_read.custom_db.fasta"),'wt') as o:
        for primer in primers:
            #   split barcode sequence from primer ID if present
            try: barcode_seq,primer = primer.split('-')
            except ValueError: barcode_seq,primer = '',primer
            if not barcode_seq == '':
                #   fwd primer with barcode
                seq_id = '-'.join([barcode_seq,primer])
                seq = ''.join([barcode_seq,primer_dict[primer]["seq"]])
                o.write(''.join(['>',seq_id,'\n',seq,'\n']))
            #   fwd or rvs primer without barcode
            seq_id,seq = primer,primer_dict[primer]["seq"]
            o.write(''.join(['>',seq_id,'\n',seq,'\n']))
    print("done",file=sys.stderr)
    
    return min_fwd_length,min_rvs_length,primer_dict,fwd_primers+rvs_primers
    
def parse_config(filepath):
    """parses LAMP config file and returns user-defined parameters"""
    file_dict = {}
    barcode_seqs = set([])
    min_barcode_length,norm_factor = None,None
    with open(filepath,'rU' if python_2 else 'rt') as f:
        #   format: dict[filepath][barcode] = {"label":,"omitted":,"norm_factor":}
        #   one line per barcode, multiple lines per file if many barcodes sequenced together
        for line in f:
            #   potential columns - library,barcode_seq,filepath,omitted_primers,batch_num,TAQMAN_factor,dilution_factor,lysate_factor
            line = line.rstrip().split('\t')
            
            library = line[0].replace(' ','_')
            barcode_seq = '' if line[1] in [' ',''] else line[1]
            if os.path.isfile(line[2]):
                filepath = line[2]
            else:
                print(''.join(["Error - sequencing file '",line[2],"' does not exist"]))
                sys.exit(-2)
            try: omitted_primers = set([]) if line[3] in [' ',''] else set(line[3].replace('"','').split(','))
            except IndexError: omitted_primers = set([])   #   no omitted primer column provided
            try: batch_num = int(line[4])
            except IndexError or ValueError: batch_num = None
            norm_factor = np.prod([float(val) for val in line[5:]])
            
            #   handle barcode sequence if present
            if not barcode_seq is '':
                barcode_seqs.add(barcode_seq)
                barcode_length = len(barcode_seq)
                min_barcode_length = barcode_length if min_barcode_length == None or barcode_length < min_barcode_length else min_barcode_length 
            
            try: 
                file_dict[library][barcode_seq]
                print("Warning - multiple declarations of a barcode for the same sequencing run found in config file",end='',file=sys.stderr)
            except KeyError: 
                try: file_dict[filepath][barcode_seq] = {"label":library,"norm_factor":norm_factor,"omitted_primers":omitted_primers,"batch_num":batch_num}
                except KeyError: file_dict[filepath] = {barcode_seq:{"label":library,"norm_factor":norm_factor,"omitted_primers":omitted_primers,"batch_num":batch_num}}
                    
    return file_dict,barcode_seqs,(min_barcode_length if not min_barcode_length is None else 0)

def create_directory(directory):
    """creates directory if it does not exist"""
    if not os.path.exists(directory):
        os.makedirs(directory)
    
    return directory

def adjust_plot_parameters(ax,y_label=None,x_label=None,x_top=False):
    """adjust common Matplotlib attributes"""
    #   adjust axes ticks and labels
    ax.tick_params(axis='x',which="both",bottom=False if x_top else True,top=True if x_top else False,direction="out",length=4,width=0.5,color='k')
    if not x_label == None:
        ax.set_xlabel(x_label)
    ax.tick_params(axis='y',which="both",left=True,right=False,direction="out",length=4,width=0.5,color='k')
    if not y_label == None:
        ax.set_ylabel(y_label)
    
    #   adjust splines
    for axis in ["bottom","left","right","top"] if x_top else ["bottom","left"]:   
        ax.spines[axis].set_visible(True)
        ax.spines[axis].set_color('k')
        ax.spines[axis].set_linewidth(1.)
    if x_top:
        ax.xaxis.set_label_position("top") 
        ax.xaxis.tick_top()
        ax.set_aspect("equal")
    else:
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
    
def BLAST_ligated_products(file_dict,word_size):
    """maps ligated paired-end reads to BLASTdb"""
    #   create custom BLAST database for ligated products
    print("\tBuilding BLAST database ... ",end='',file=sys.stderr)
    with open(''.join([blastdb_dir,"primer_pairs_makeblastdb.log"]),'wt') as o:
        subprocess.check_call(["makeblastdb","-in",os.path.join(primer_dir,"primer_pairs.custom_db.fasta"),"-dbtype","nucl","-out",os.path.join(blastdb_dir,"primer_pairs.custom_db")],stdout=o,stderr=subprocess.STDOUT)
    print("done",file=sys.stderr)
     
    filepaths = []
    read_counts = [[],[]]  #[[total_reads1,total_reads2,...],[mapped_reads1,mapped_reads2,...]]
    keys = file_dict.keys()
    for key in keys:
        basename = '.'.join(os.path.basename(key).split('.')[:-1])
        
        #   check if BAM format provided as input and convert to FASTQ
        if key.lower().endswith(".bam"):
            print("\tConverting '"+basename+"' from BAM to FASTQ ... ",end='',file=sys.stderr)
            output_filepath = os.path.join(mapping_dir,'.'.join([basename,"fastq"]))
            with open(os.path.join(mapping_dir,'.'.join([basename,"BAMtoFASTQ.log"])),'wt') as o:
                subprocess.check_call(''.join(["samtools bam2fq ",key," > ",output_filepath]),shell=True,stdout=o,stderr=subprocess.STDOUT)
            #   update dictionary to point towards FASTQ
            file_dict[output_filepath] = file_dict.pop(key)
            key = output_filepath
            print("done",file=sys.stderr)
            
        print("\tConverting '"+basename+"' from FASTQ to FASTA ... ",end='',file=sys.stderr)
        #   convert fastq to fasta
        output_filepath = os.path.join(mapping_dir,'.'.join([basename,"fasta"]))
        #   commands:
        #       1) extract sequence string from fastq file
        #       2) exclude sequences <min(primer_junction) in length
        #       3) sort and count unique sequences
        #       4) output in fasta format
        with open(os.devnull,'wb') as devnull:
            subprocess.check_call(''.join([sed_cmd," -n '2~4p' ",key," | awk 'length($0)>=",word_size,"' | sort | uniq -c | awk 'BEGIN{SUM=0}{SUM+=1; print \">ID_\"SUM\"_FREQ_\"$1 \"\\n\" $2}' > ",output_filepath]),shell=True,stdout=devnull)

        #   calculate total number of reads to be mapped
        filepaths.append(key)
        read_counts[0].append(int(subprocess.check_output(''.join([sed_cmd," -n '1~2p' ",output_filepath," | awk -F '_' 'BEGIN{sum=0}{sum+=$NF}END{print sum}'"]),shell=True).rstrip()))
        print("done",file=sys.stderr)
        
        print(''.join(["\tBLAST-ing '",basename,"' ... "]),end='',file=sys.stderr)
        #   blast sample queries against database
        query_filepath = output_filepath
        output_filepath = output_filepath.replace(".fasta",".BLAST.tsv")
        BLAST_process = subprocess.Popen(["blastn","-db",os.path.join(blastdb_dir,"primer_pairs.custom_db"),"-word_size",word_size,"-dust","no","-max_target_seqs","1","-num_threads",str(num_threads),"-outfmt",r"6 qseqid qstart qseq sseqid sstart sseq bitscore","-query",query_filepath,"-out",output_filepath],stdout=subprocess.PIPE,stderr=subprocess.STDOUT)
        with open(output_filepath.replace("tsv","log"),'wt') as o:
            o.write(BLAST_process.communicate()[0].decode("utf-8"))
        read_counts[1].append(int(subprocess.check_output(''.join(["awk -F \"\\t\" '{print $1}' ",output_filepath," | sort | uniq | awk -F \"_\" 'BEGIN{sum=0}{sum+=$4}END{print sum}' "]),shell=True).rstrip()))
        print("done",file=sys.stderr)

    #   write mapped and unmapped read counts to file
    with open(os.path.join(mapping_dir,''.join(["BLAST_summary.word_size_",word_size,".tsv"])),'wt') as o:
        o.write('\t'.join(["filepath","total_reads","mapped_reads"])+'\n')
        for i,(mapped,unmapped) in enumerate(zip(*read_counts)):
            assert(unmapped < mapped),''.join(["Error - more reads were mapped than sequenced for library '",filepaths[i],"'"])
            o.write('\t'.join([filepaths[i],str(mapped),str(unmapped)])+'\n')

    print("\tPlotting summary ... ",end='',file=sys.stderr)
    for i,key in enumerate(filepaths):
        output_filepath = os.path.join(mapping_dir,os.path.basename(key).replace(".fastq",".BLAST_summary.bar_plot.png"))
        plot_bar(output_filepath,[read_counts[0][i],read_counts[1][i]],["Total","Mapped"],['b','r'])
    print("done",file=sys.stderr)

def BLAST_short_reads(dictionary,word_size):
    """maps 'unmappables' to short-read BLASTdb"""
    print("\tBuilding short-read BLASTdb ... ",end='',file=sys.stderr)
    with open(os.path.join(blastdb_dir,"short_read.makeblastdb.log"),'wt') as o:
        subprocess.check_call(["makeblastdb","-in",os.path.join(primer_dir,"short_read.custom_db.fasta"),"-dbtype","nucl","-out",os.path.join(blastdb_dir,"short_read.custom_db")],stdout=o,stderr=subprocess.STDOUT)
    print("done",file=sys.stderr)
    
    total_dict = {}
    label_dict = {'u':"Unknown",'F':"FWD",'R':"RVS",'B':"BC+FWD"}
    color_dict={"Total":'k',"Unknown":"grey","FWD":'b',"RVS":'r',"BC+FWD":'w'}
    for key in dictionary.keys():
        #if not "15Min_Set1" in key:
        #    continue
        count_dict = {}
        basename = '.'.join(os.path.basename(key).split('.')[:-1])
        print(''.join(["\tProcessing '",basename,"' reads ... "]),end='',file=sys.stderr)
        filepath = os.path.join(mapping_dir,'.'.join([basename,"BLAST.tsv"]))
        
        #   get mapped IDs
        mapped_ids = set(subprocess.check_output(''.join(["awk -F '\\t' '{print $1}' ",filepath]),shell=True).decode("utf-8").rstrip().split())
        #   get total number of mapped reads
        total_mapped = sum([int(seq_id.split('_')[-1]) for seq_id in mapped_ids])
        
        #   parse out 'unmappables'
        unmapped_seq_dict = {}
        total_unmapped,total = 0,0
        input_filepath = os.path.join(mapping_dir,'.'.join([basename,"fasta"]))
        output_filepath = os.path.join(short_read_dir,'.'.join([basename,"unmapped_reads.fasta"]))
        with open(output_filepath,'wt') as o:
            with open(input_filepath,'rt') as f:
                seq_id = f.readline().rstrip()
                seq = f.readline().rstrip()
                while seq_id and seq:
                    if not seq_id[1:] in mapped_ids:
                        o.write('\n'.join([seq_id,seq,'']))
                        seq_id,freq = [val for val in seq_id.split('_') if val.isdigit()]
                        total_unmapped += int(freq)
                        total += int(freq)
                        unmapped_seq_dict["ID_"+seq_id] = [freq,seq]
                    else:   #   track total reads
                        seq_id,freq = [val for val in seq_id.split('_') if val.isdigit()]
                        total += int(freq)
                    #   read in next pair of lines
                    seq_id = f.readline().rstrip()
                    seq = f.readline().rstrip()
        assert(total_mapped+total_unmapped == total),"Error - mapped and umapped read counts do not sum to the total number of reads"
        print("done",file=sys.stderr)
        
        #   BLAST 'unmappables' against  short-read BLASTdb
        print(''.join(["\tBLAST-ing '",basename,"' short reads ... "]),end='',file=sys.stderr)
        query_filepath = output_filepath
        output_filepath = output_filepath.replace(".unmapped_reads.fasta",".short_reads.BLAST.tsv")
        BLAST_process = subprocess.Popen(["blastn","-db",os.path.join(blastdb_dir,"short_read.custom_db"),"-word_size",word_size,"-dust","no","-max_target_seqs",'1',"-max_hsps",'1',"-num_threads",num_threads,"-outfmt",r"6 qseqid qstart qseq sseqid sstart sseq bitscore","-query",query_filepath,"-out",output_filepath],stdout=subprocess.PIPE,stderr=subprocess.STDOUT)
        with open(output_filepath.replace("tsv","log"),'wt') as o:
            o.write(BLAST_process.communicate()[0].decode("utf-8"))
        
        #   reformat BLAST output
        input_filepath = output_filepath
        output_filepath = input_filepath.replace(".short_reads.BLAST.tsv",".short_reads_summary.tsv")
        with open(output_filepath,'wt') as o:
            #   write header
            o.write('\t'.join(["seq_id","freq","seq","possible_primer"])+'\n')
        with open(output_filepath,'at') as o:
            subprocess.Popen(["awk","-F","\\t|_","{OFS=\"\\t\"; print $1\"_\"$2,$4,$6,$7}",input_filepath],stdout=o,stderr=subprocess.STDOUT).communicate()
        #   add unmapped sequences to output
        mapped_short_reads = subprocess.check_output(''.join(["awk -F '\\t' 'NR>1{print $1}' ",output_filepath]),shell=True).rstrip().split()
        assert(len(mapped_short_reads) == len(set(mapped_short_reads))),"Error - multimapping encountered in short reads library"
        with open(output_filepath,'at') as o:
            for seq_id in set(unmapped_seq_dict.keys())-set(mapped_short_reads):
                o.write('\t'.join([seq_id]+unmapped_seq_dict[seq_id]+["unknown"])+'\n')

        #   store summary of sample in memory
        process_1 = subprocess.Popen(["awk","-F","\\t","NR>1 {print $2,$4}",output_filepath],stdout=subprocess.PIPE,stderr=subprocess.PIPE)
        process_2 = subprocess.Popen(["awk","{a[substr($2,1,1)]+=$1}END{for(i in a) print i,a[i]}"],stdin=process_1.stdout,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
        output = process_2.communicate()[0].decode("utf-8").rstrip().split('\n')
        for i,val in enumerate(output):
            label,count = val.split()
            label = 'B' if label.lower() in ['a','c','g','t'] else label
            try: count_dict[label_dict[label]] += int(count)
            except KeyError: count_dict[label_dict[label]] = int(count)
        
        heights,colors = [],[]
        labels = ["Total"]+sorted(count_dict.keys())
        for label in labels:
            heights.append(sum(count_dict.values()) if label == "Total" else count_dict[label])
            colors.append(color_dict[label])
        plot_bar(output_filepath.replace(".tsv",".bar_plot.png"),heights,labels,colors,title="Short-reads summary")
        
        total_dict[key] = sum(count_dict.values())
        assert(total_dict[key] == total_unmapped),"Error - short reads do not sum to expected count"
        print("done",file=sys.stderr)
        
    return total_dict

def plot_heatmap(output_filepath,matrix,labels):
    """creates heatmap of matrix via MatPlotlib"""
    n = len(labels)
    labels = [val.split('@')[1] for i,val in enumerate(labels)]
    #   create heatmap of bitwise scores
    fig,ax = plt.subplots(figsize=(8,8),facecolor='w')
    im = ax.imshow(matrix,cmap="gist_heat_r",interpolation="none")
    ax.grid(False)
    cb = plt.colorbar(im,fraction=0.046,pad=0.04)
    cb.outline.set_visible(False)
    cb.ax.tick_params(axis='y', direction='out')
    #create ticks
    tick_labels,tick_positions = [],[]
    for i in get_range(0,n,n//10):
        tick_positions.append(i)
        tick_labels.append(labels[i])
    plt.xticks(tick_positions,tick_labels,rotation=-45,ha="right")
    plt.yticks([val for val in tick_positions],tick_labels)
    adjust_plot_parameters(ax,x_top=True)
    plt.tight_layout()
    plt.savefig(output_filepath,format="png",transparent=False,dpi=300.0)
    plt.close()

def plot_bar(output_filepath,y_vals,labels,colors,title="BLAST results of libraries"):
    """creates bar plot of provided heights via MatPlotlib"""
    if len(y_vals) == 2:
        title = title.replace("libraries","library")
    fig,ax = plt.subplots(figsize=(8,8))
    ax.set_title(title,fontweight="bold")
    x_ticks = []
    for i,(height,color) in enumerate(zip(y_vals,colors)):
        ax.bar(i,height,1.,color=color,edgecolor='k')
        x_ticks.append(i+0.5)
    ax.grid(False)
    ax.set_xlim([-1,6])
    ax.set_xticks(x_ticks)
    ax.set_xticklabels(labels,rotation=40,ha="right")
    adjust_plot_parameters(ax,y_label="Read count frequency")
    plt.savefig(output_filepath,format="png",transparent=True,dpi=300.0, bbox_inches="tight")
    plt.close()

def process_mapped_reads(file_dict,primer_dict,sorted_primers,short_read_dict,data_type="2C-ChIP"):
    """process mapped reads and outputs them in user-friendly format"""
    n = len(sorted_primers)
    #   build matrix loci row/column labels
    loci = []
    chrom = set([])
    min_start,max_end = None,None
    for primer in sorted_primers:
        end = primer_dict[primer]["start"]+primer_dict[primer]["length"]
        loci.append('|'.join([primer,primer_dict[primer]["assembly"],''.join(["chr",primer_dict[primer]["chrom"],':',str(primer_dict[primer]["start"]),'-',str(end)])]))
        if data_type == "2C-ChIP":
            chrom.add(primer_dict[primer]["chrom"])
            min_start = primer_dict[primer]["start"] if min_start == None or primer_dict[primer]["start"] < min_start else min_start
            max_end = end if max_end == None or end > max_end else max_end
    if data_type == "2C-ChIP" and len(chrom) > 1:
        print("Warning - multiple chromosomes detected in 2C-ChIP sample. bedGraphs will not load to a specific region",file=sys.stderr)
        region = ''
    elif data_type == "2C-ChIP":
        region = ''.join(["chr",chrom.pop(),':',str(min_start),'-',str(max_end)])
    
    expected_barcodes = [''.join([key,barcode]) for key in file_dict.keys() for barcode in file_dict[key].keys()]
    tensor = np.zeros((len(expected_barcodes),n,n),dtype=np.int32)
    with open(os.path.join(mapping_dir,"read_count_frequency_table.tsv"),'wt') as o:
        o.write('\t'.join(["filepath","total_reads"]+(["on_diagonal","off_diagonal"] if data_type == "2C-ChIP" else ["expected","unexpected"])+["short_reads"])+'\n')
        for key in file_dict.keys():
            total = 0
            cur_barcodes = set([])
            basename = '.'.join(os.path.basename(key).split('.')[:-1])
            print(''.join(["\tParsing '",basename,"' mapped reads ... "]),end='',file=sys.stderr)
            filepath = os.path.join(mapping_dir,'.'.join([basename,"BLAST.tsv"]))
            
            #   check for multi-mappings
            multimap_dict = {}
            output = subprocess.check_output(''.join(["awk -F '\t' '{a[$1]++}END{for(i in a){if(a[i] > 1){print i}}}' ",filepath,"  | xargs -I {} grep {} ",filepath]),shell=True).decode("utf-8").rstrip()
            for line in output.split('\n'):
                try: query_id,query_start,query_seq,subject_id,subject_start,subject_seq,bitscore = line.split()
                except ValueError: continue #   no multi-mappings in library
                try: multimap_dict[query_id].add(subject_id)
                except KeyError: multimap_dict[query_id] = set([subject_id])

            single_primer_count=[0,0]
            invalid_barcode,omitted,multimaps,total,off_diagonal,on_diagonal,expected,unexpected = 0,0,0,0,0,0,0,0
            with open(filepath,'rt') as f:
                for line in f:
                    query_id,query_start,query_seq,subject_id,subject_start,subject_seq,bitscore = line.rstrip().split()
                    freq = int(query_id.split('_')[-1])
                    
                    #   check for multi-mappings
                    try:
                        num_hits = len(multimap_dict[query_id])
                        if num_hits > 1:    #   canonical multi-mapping
                            multimaps += freq
                        elif num_hits == 1:
                            #   multi-mapping, but only to one fragment pair - retain mapping
                            multimap_dict[query_id] = []
                        else: continue  #   multi-mapping, but only to one fragment pair - already retained previous hit, skip
                    except KeyError: pass   #   no multi-mapping
                    
                    total += freq
                    primer_A,primer_B = subject_id.split('|')

                    #   split barcode and primer
                    try: barcode_A,primer_A = primer_A.split('-')
                    except ValueError: barcode_A = ''
                    full_barcode = ''.join([key,barcode_A])

                    #   validate barcodes
                    if (barcode_A == '' and len(expected_barcodes) == 0) or (full_barcode in expected_barcodes):
                        #   skip omitted libraries where all primers are ommitted (still expect mappings for 2C-ChIP products)
                        subject_start = int(subject_start)
                        #ensure mapped sequences map to both primers with at least one nucleotide coverage
                        #   1)  mapped query doesn't start on downstream primer
                        #   2)  query that maps to upstream primer covers downstream primer
                        if ((len(barcode_A)+len(primer_dict[primer_A]["seq"])) > subject_start) and ((subject_start+len(query_seq)) > len(primer_dict[primer_A]["seq"])):   
                            #insert valid mapping into matrix
                            index_A = sorted_primers.index(primer_A)
                            index_B = sorted_primers.index(primer_B)
                            tensor[expected_barcodes.index(full_barcode)][index_A][index_B] += freq
                            cur_barcodes.add(full_barcode)
                            if any(val in file_dict[key][barcode_A]["omitted_primers"] for val in [index_A,index_B]):
                                omitted += freq
                            elif primer_A.split('@')[0][1:] == primer_B.split('@')[0][1:] and primer_A[0] == 'F' and primer_B[0] == 'R' and data_type == "2C-ChIP":
                                on_diagonal += freq
                            elif data_type == "2C-ChIP":
                                off_diagonal += freq
                            elif primer_A[0] == 'F' and primer_B[0] == 'R' and data_type == "5C":
                                expected += freq
                            else:
                                unexpected +=freq
                        else:
                            single_primer_count[0 if primer_A[0] == 'F' else 1] += freq
                    else:
                        try:
                            if file_dict[key][barcode_A]["omitted_primers"] == "all":
                                omitted += freq
                        except KeyError: invalid_barcode += freq
        
            for barcode in cur_barcodes:
                #   write raw matrix of frequency counts to storage
                i = expected_barcodes.index(barcode)
                output_filepath = os.path.join(results_dir,'.'.join([file_dict[key][barcode.split(".fastq")[-1]]["label"],"raw.matrix"]))
                np.savetxt(output_filepath,np.hstack((np.reshape(loci,(n,1)),tensor[i].astype(str))),
                           fmt='%s',delimiter='\t',header='\t'.join(['']+loci),comments='')
                #   plot heatmap of log(values) in matrix
                plot_heatmap(output_filepath.replace('raw.matrix','log_raw.heatmap.png'),np.log(tensor[i]+1.0),sorted_primers)
            vals = [total+short_read_dict[key],expected,sum([unexpected,omitted,multimaps,invalid_barcode,sum(single_primer_count)]),short_read_dict[key]] if data_type == "5C" else [total+short_read_dict[key],on_diagonal,sum([off_diagonal,multimaps,omitted,invalid_barcode,sum(single_primer_count)]),short_read_dict[key]]
            o.write('\t'.join([key]+[str(val) for val in vals])+'\n')
            assert(vals[0] == sum(vals[1:])),''.join(["Error - '",key,"' LP-mapped reads do not sum to total"])
            #   plot bar of total counts
            output_filepath = os.path.join(mapping_dir,'.'.join([basename,"read_count.bar_plot.png"]))
            plot_bar(output_filepath,vals,["Total","Expected","Unexpected","Short reads"] if data_type == "5C" else ["Total","On-diagonal","Off-diagonal","Short reads"],['k','b','r','grey'])
            print("done",file=sys.stderr)

    if data_type == "2C-ChIP":
        print("\tNormalizing data ... ",end='',file=sys.stderr)
        #   determine track batches to be normalized be respective input
        batch_dict = {} #   {number: {"labels":[],"indices":[],"omitted_primers":[],"input":True/False}}
        for key in file_dict.keys():
            for barcode in file_dict[key].keys():
                try: batch_dict[file_dict[key][barcode]["batch_num"]]
                except KeyError: batch_dict[file_dict[key][barcode]["batch_num"]] = {"tracks":[],"input":False}
                #   check for input track
                try:
                    re.search(r'.*input.*',file_dict[key][barcode]["label"].lower()).group(0)
                    batch_dict[file_dict[key][barcode]["batch_num"]]["input"] = True
                    batch_dict[file_dict[key][barcode]["batch_num"]]["tracks"] = [tuple([key,barcode])]+batch_dict[file_dict[key][barcode]["batch_num"]]["tracks"]
                except AttributeError:
                    batch_dict[file_dict[key][barcode]["batch_num"]]["tracks"].append(tuple([key,barcode]))
 
        #   iterate over batches
        offset = None
        for i,val in enumerate(sorted_primers):
            if val.startswith('R'):
                offset = i
                break
        bedGraph_loci,primers = [],[]
        for j in get_range(0,offset,1):
            primer_A,primer_B = sorted_primers[j].split('@')[0],sorted_primers[j+offset].split('@')[0]
            strand_A,primer_A = primer_A[0],primer_A[1:]
            strand_B,primer_B = primer_B[0],primer_B[1:]
            assert(strand_A == 'F' and strand_B == 'R' and primer_A == primer_B),''.join(["Error - FWD and RVS primer labels do not match for row/col (",str(j),',',str(j+offset),')'])
            bedGraph_loci.append('\t'.join(re.split(r':|-',loci[j].replace("omosome_",''))[:2]+[loci[j+offset].split('-')[-1]]))
            if j%10 == 0.0:
                primers.append(sorted_primers[j].split('@')[0][1:])
                
        raw_totals = []
        for key in batch_dict.keys():
            #   warn user if no input track present for batch
            if not batch_dict[key]["input"]:
                print(''.join(["Warning - There was no 'input' track found for batch #",str(key)]),file=sys.stderr)
            
            #   iterate over batch tracks and normalize on-diagonals targets
            data_dict = {}
            input_norm = []
            norm_max_val,raw_max_val = 0.0,0
            for k,(filepath,barcode) in enumerate(batch_dict[key]["tracks"]):
                label = file_dict[filepath][barcode]["label"]
                data_dict[label] = {}
                i = expected_barcodes.index(''.join([filepath,barcode]))
                
                #   extract raw read counts and check primer pairing
                data_dict[label]["raw"] = [0.0 if j in file_dict[filepath][barcode]["omitted_primers"] else float(tensor[i][j][j+offset]) for j in get_range(0,offset,1)]
                raw_totals.append(tuple([label,sum(data_dict[label]["raw"]),file_dict[filepath][barcode]["norm_factor"]]))
                
                #   read count normalize diagonal (reads per million [RPM] for total observed reads associated with barcode for given sequencing run)
                total = np.sum(tensor[i])
                data_dict[label]["RPM"] = [(val/total)*1000000 for val in data_dict[label]["raw"]]
                
                #   normalize for DNA density in sequencing sample
                data_dict[label]["norm"] =[val*file_dict[filepath][barcode]["norm_factor"] for val in data_dict[label]["RPM"]]
                
                #   retain input normalization values if present
                if k == 0 and batch_dict[key]["input"]:
                    input_norm = data_dict[label]["norm"]
                    for i,val in enumerate(input_norm):
                        if np.isclose(val,0.0):
                            chrom,start,end = bedGraph_loci[i].split('\t')
                            print(''.join(["Warning - input value of zero for primer in ",label,". Primer-pair locus: ",chrom,':',start,"-",end]),file=sys.stderr)
                
                #   normalize by input if present
                if batch_dict[key]["input"]:
                    vals = []
                    for i,(val,factor) in enumerate(zip(data_dict[label]["norm"],input_norm)):
                        vals.append(val/factor if factor > 0.0 else 0.0)    #   else case is a catch for underperforming primers
                    data_dict[label]["norm"] = vals

                #   calculate max read counts for any library in batch to set bedgraph y-axis height
                norm_max_val = max(data_dict[label]["norm"]) if max(data_dict[label]["norm"]) > norm_max_val else norm_max_val
                raw_max_val = max(data_dict[label]["raw"]) if max(data_dict[label]["raw"]) > raw_max_val else raw_max_val
            
            n = len(data_dict[label]["norm"])
            x_vals = get_range(0,n,1)
            for filepath,barcode in batch_dict[key]["tracks"]:
                label = file_dict[filepath][barcode]["label"]
                #   write bedGraph files - need to iterate again to get proper 'max_val' in previous loop
                for max_val,data_type in zip([raw_max_val,norm_max_val],["raw","norm"]):
                    output_filepath = os.path.join(results_dir,'.'.join([label,data_type,"bedGraph"]))
                    with open(output_filepath,'wt') as o:
                        #   write header
                        o.write('\n'.join([''.join(["browser position ",region]),"browser hide all","browser pack refGene encodeRegions",''.join(["track type=bedGraph name=\"",label,"\" color=0,0,0 visibility=full priority=20 maxHeightPixels=60 autoScale=off viewLimits=0:",str(max_val)]),'']))
                        for i,(loci,freq) in enumerate(zip(bedGraph_loci,data_dict[label][data_type])):
                            if not i+offset in file_dict[filepath][barcode]["omitted_primers"]:
                                o.write('\t'.join([loci,str(freq)])+'\n')
                
                #   plot line graphs of primer performance
                for data_type in ["RPM","norm"]:
                    fig,ax = plt.subplots(figsize=(8,8),facecolor='w')
                    output_filepath = os.path.join(results_dir,'.'.join([label,data_type,"line_plot.png"]))
                    ax.plot(x_vals,data_dict[label][data_type],label=data_type,color='k' if data_type == "RPM" else 'r')
                    #   set x-ticks
                    xticks = get_range(0,n,10)
                    ax.set_xticks(xticks)
                    ax.set_xticklabels(primers[:len(list(xticks))])
                    ax.grid(color="gray",ls=':',lw=0.5)
                    adjust_plot_parameters(ax,y_label="Raw frequency (in RPM)" if data_type == "RPM" else "Normalized frequency",x_label="Primer pair")
                    ax.margins(0.025)
                    plt.tight_layout()
                    plt.savefig(output_filepath.replace(".tsv",".scatter.png"),format="png",transparent=False,dpi=300.0)
                    plt.close()
        
        #   create summary of F-R raw frequencies vs. norm-factor
        output_filepath = os.path.join(results_dir,"raw_totals_vs_norm_factor.tsv")
        with open(output_filepath,'wt') as f:
            f.write('\t'.join(["library","total_raw_reads","norm_factor"])+'\n')
            for (label,total,norm_val) in raw_totals:
                f.write('\t'.join([label,str(total),str(norm_val)])+'\n')
        #   create scatter of F-R raw frequencies vs. norm-factor
        fig,ax = plt.subplots(figsize=(8,8),facecolor='w')
        x_vals,y_vals = zip(*raw_totals)[1:] if python_2 else list(zip(*raw_totals))[1:]
        ax.plot(x_vals,y_vals,'ok',alpha=0.75,clip_on=False,markersize=5.0)
        #   correlate values
        rho,p_value = spearmanr(x_vals,y_vals)
        ax.set_title(''.join([r'$\rho_{S}$ = ',"{0:.2f} ({1:.2E})".format(rho,p_value)]),fontsize=12)
        ax.grid(None)
        ax.set_ylim([1.0,np.power(10,math.ceil(np.log10(val)))])
        ax.set_yscale("log",nonposy="clip")
        ax.grid(color="gray",ls=':',lw=0.5)
        adjust_plot_parameters(ax,y_label="Normalization factor",x_label="Total reads")
        plt.tight_layout()
        plt.savefig(output_filepath.replace(".tsv",".scatter.png"),format="png",transparent=False,dpi=300.0)
        plt.close()
        print("done",file=sys.stderr)
        
parser = argparse.ArgumentParser()
parser.add_argument("config", help = "path to LAMPS config file", type = str)
parser.add_argument("primers", help = "path to TSV file containing primer information", type = str)
parser.add_argument("type", help = "source of paired-end reads:[2C-ChIP,5C]")
parser.add_argument("output", help = "path to output folder", type = str)
parser.add_argument("--num_cpus", help = "set the number of cpus - default = num_cpus-2", type = str)
parser.add_argument("--word_size", help = "set the minimum required sequence length for processing", type = int)
args = parser.parse_args()

assert(args.type in ["2C-ChIP","5C"]),"Error - invalid type provided. Please specify either '2C-ChIP' or '5C'"

num_threads = str(multiprocessing.cpu_count()-2) if args.num_cpus == None else str(args.num_cpus)

args.output = create_directory(args.output)
primer_dir = create_directory(os.path.join(args.output,"primer_files",''))
mapping_dir = create_directory(os.path.join(args.output,"mapping_files",''))
blastdb_dir = create_directory(os.path.join(mapping_dir,"BLAST",''))
short_read_dir = create_directory(os.path.join(mapping_dir,"short_read_analysis",''))
results_dir = create_directory(os.path.join(args.output,"results",''))

short_read_dict = None

print("Parsing LAMPS config file ... ",end='',file=sys.stderr)
#   parse config file
file_dict,barcode_seqs,min_barcode_length = parse_config(args.config)
print("done",file=sys.stderr)

#   build FASTA file of primer pairs
print("Parsing primer file")
min_fwd_length,min_rvs_length,primer_dict,sorted_primers = build_FASTAs(args.primers,file_dict,barcode_seqs)

#   map paired-end reads to ligated product BLASTdb
print("Mapping to ligated-product BLASTdb")
word_size = args.word_size if not args.word_size == None else str(min((min_barcode_length+min_fwd_length)*2,
                    min_barcode_length+min_fwd_length+min_rvs_length,
                    min_rvs_length+min_rvs_length))
BLAST_ligated_products(file_dict,word_size)

#   map too short or unmappable pair-end reads to short-read BLASTdb
print("Mapping to short-read BLASTdb")
#   first case is unecessary, only present for completeness
word_size = str(min(min_barcode_length+min_fwd_length,min_fwd_length,min_rvs_length))
short_read_dict = BLAST_short_reads(file_dict,word_size)

#   process mapped reads
print("Processing mapped reads")
process_mapped_reads(file_dict,primer_dict,sorted_primers,short_read_dict,args.type)
