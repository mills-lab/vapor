import os
global invert_base
import itertools
import matplotlib.pyplot as plt
import numpy as np
import random
import rpy2
from rpy2.robjects.packages import importr
from rpy2.robjects.vectors import IntVector, FloatVector, StrVector
import scipy
from scipy.cluster.vq import vq, kmeans, whiten
from scipy import stats
from scipy.stats import linregress
from scipy.spatial import distance
from sklearn import cluster
from sklearn import cluster
import sys
import valor_vali.plotting as plotting
global invert_base
invert_base = { 'A' : 'T', 'T' : 'A', 'C' : 'G', 'G' : 'C','N' : 'N','a' : 't', 't' : 'a', 'c' : 'g', 'g' : 'c','n' : 'n'}
global default_flank_length
default_flank_length=500
global default_read_length
default_read_length=4000    #average length of pacbio read
global default_max_sv_test
default_max_sv_test=10000 #when size of a sv block excessed default_max_sv_test, try junctions instead of event
def alt_seq_readin(ref,info,flank_length):
    #eg of info=('abcde/abcde', 'db/db', '1', '14436355', '14436655', '14437075', '14438241', '14438566', '14438940']
    block_hash={}
    for x in range(len(info)-4):
        block_hash[chr(97+x)]=[info[2],int(info[x+3]),int(info[x+4])]
    block_hash['-']=[block_hash[info[0][0]][0],block_hash[info[0][0]][1]-flank_length,block_hash[info[0][0]][1]]
    block_hash['+']=[block_hash[info[0][-1]][0],block_hash[info[0][-1]][2],block_hash[info[0][-1]][2]+flank_length]
    if not info[1].split('/')[0]==info[1].split('/')[1]:
        out=[ref_seq_readin(ref,block_hash['-'][0],block_hash['-'][1],block_hash['-'][2]),ref_seq_readin(ref,block_hash['-'][0],block_hash['-'][1],block_hash['-'][2])]
        rec=-1
        for x in info[1].split('/'):
            rec+=1
            temp_x=[]
            for y in x:
                if not y=='^':temp_x.append(y)
                else:    temp_x[-1]+=y
            for y in temp_x:
                if not '^' in y:
                    out[rec]+=ref_seq_readin(ref,block_hash[y][0],block_hash[y][1],block_hash[y][2])
                else:
                    out[rec]+=ref_seq_readin(ref,block_hash[y[0]][0],block_hash[y[0]][1],block_hash[y[0]][2],'TRUE')
        out[0]+=ref_seq_readin(ref,block_hash['+'][0],block_hash['+'][1],block_hash['+'][2])
        out[1]+=ref_seq_readin(ref,block_hash['+'][0],block_hash['+'][1],block_hash['+'][2])
        return out
    else:
        out=[ref_seq_readin(ref,block_hash['-'][0],block_hash['-'][1],block_hash['-'][2])]
        rec=-1
        for x in [info[1].split('/')[0]]:
            rec+=1
            temp_x=[]
            for y in x:
                if not y=='^':temp_x.append(y)
                else:    temp_x[-1]+=y
            for y in temp_x:
                if not '^' in y:
                    out[rec]+=ref_seq_readin(ref,block_hash[y][0],block_hash[y][1],block_hash[y][2])
                else:
                    out[rec]+=ref_seq_readin(ref,block_hash[y[0]][0],block_hash[y[0]][1],block_hash[y[0]][2],'TRUE')
        out[0]+=ref_seq_readin(ref,block_hash['+'][0],block_hash['+'][1],block_hash['+'][2])
        return [out[0],out[0]]
def bam_in_decide(bam_in,bps):
    if os.path.isfile(bam_in):
        return [bam_in]
    else:
        temp_bam_in=[]
        bam_in_path='/'.join(bam_in.split('/')[:-1])+'/'
        if 'XXX' in bam_in.split('/')[-1]:
            bam_in_keys=bam_in.split('/')[-1].split('XXX')
        elif '*' in bam_in.split('/')[-1]:
            bam_in_keys=bam_in.split('/')[-1].split('*')
        else:
            print 'Error: invalid name for pacbio files !'
        for k1 in os.listdir(bam_in_path):
            if k1.split('.')[-1]==bam_in.split('.')[-1]:
                flag=0
                for y in bam_in_keys:
                    if not y in k1:
                        flag+=1
                if flag==0:
                    temp_bam_in.append(bam_in_path+k1)
        return temp_bam_in
def block_modify(block,chromos):
    #eg of block=['chr16', '34911339', '34913149', 'chr16', '34913149', '34913438']
    out=[]
    for x in block:
        if x in chromos:
            if out==[]:    out.append([x])
            else:
                if not x in out[-1]:    out.append([x])
        else:   out[-1].append(x)
    out_new=[]
    for x in out:
        out_new.append([])
        for y in x:
            if x.count(y)==1:
                out_new[-1].append(y)
    out_new_2=[]
    for x in out_new:
        if len(x)==3:
            out_new_2.append(x)
        else:
            for y in range((len(x)-1)/2):
                out_new_2.append([x[0],x[2*y+1],x[2*y+2]])
    return out_new_2
def bp_to_chr_hash(bps,chromos,flank_length=500):
    #eg of bps=['chr16', '34910548', '34911339', '34913149', '34913438', '36181068', '36181482']
    temp1=[]
    for i in bps:
        if i in chromos:
            temp1.append([i])
        else:
            temp1[-1].append(i)
    out={}
    rec=-1
    for k1 in temp1:
        for k2 in range(len(k1[2:])):
            rec+=1
            out[chr(97+rec)]=[k1[0],k1[k2+1],k1[k2+2]]
    out['+']=[out[sorted(out.keys())[-1]][0],out[sorted(out.keys())[-1]][2],str(int(out[sorted(out.keys())[-1]][2])+flank_length)]
    out['-']=[out['a'][0],str(int(out['a'][1])-flank_length),int(out['a'][1])]
    return out
def block_around_check(alt_allele,ref_allele):
    #eg of alt_allele='abcab'  eg of ref_allele='abcd'
    alt_juncs=[[['-']+letter_split(alt_allele)+['+']][0][j:j+2] for j in range(len(letter_split(alt_allele))+1)]
    ref_juncs=[[['-']+letter_split(ref_allele)+['+']][0][j:j+2] for j in range(len(letter_split(alt_allele))+1)]
    new_juncs=[i for i in alt_juncs if not i in ref_juncs]
    return new_juncs
def block_subsplot(bp_list,chromos):
    #eg of bp_list=['chr1', '9061390', '9061470', 'chr14', '93246136', '93248034']
    out=[]
    for x in bp_list:
        if not x in chromos:    out[-1].append(int(x))
        else:   out.append([x])
    return out
def calcu_log10(x):
    if x==0:
        return 0
    else:
        return np.log10(x)
def calcu_valor_single_read_score_abs_dis_m1(ref_seq,alt_seq,pacbio_read_info,window_size):
    ref_dotdata=dotdata(window_size,pacbio_read_info[0],ref_seq[pacbio_read_info[1]:])
    alt_dotdata=dotdata(window_size,pacbio_read_info[0],alt_seq[pacbio_read_info[1]:])
    if len(ref_dotdata)>0 and len(alt_dotdata)>0:
        [ref_clean_dotdata,ref_kept_segs]=clean_dotdata_m1(ref_dotdata)
        [alt_clean_dotdata,alt_kept_segs]=clean_dotdata_m1(alt_dotdata)
        ref_left=[i for i in ref_dotdata if not list(i) in ref_clean_dotdata]
        alt_left=[i for i in alt_dotdata if not list(i) in alt_clean_dotdata]
        [ref_anti_diag_clean_dotdata,ref_anti_diag_kept_segs]=clean_dotdata_anti_diagnal_m1b(ref_left)
        [alt_anti_diag_clean_dotdata,alt_anti_diag_kept_segs]=clean_dotdata_anti_diagnal_m1b(alt_left)
        ref_clean_dotdata+=ref_anti_diag_clean_dotdata
        alt_clean_dotdata+=alt_anti_diag_clean_dotdata
        if len(ref_clean_dotdata)>0 and len(alt_clean_dotdata)>0:
            return [eu_dis_abs_calcu(ref_clean_dotdata),eu_dis_abs_calcu(alt_clean_dotdata)]
        else:    
            return [1,1]
    else:
        return [1,1]
def calcu_valor_single_read_score_abs_dis_m1b(ref_seq,alt_seq,pacbio_read_info,window_size):
    ref_dotdata=dotdata(window_size,pacbio_read_info[0],ref_seq[pacbio_read_info[1]:])
    alt_dotdata=dotdata(window_size,pacbio_read_info[0],alt_seq[pacbio_read_info[1]:])
    if len(ref_dotdata)>0 and len(alt_dotdata)>0:
        [ref_clean_dotdata,ref_kept_segs]=clean_dotdata_diagnal_m1b(ref_dotdata)
        [alt_clean_dotdata,alt_kept_segs]=clean_dotdata_diagnal_m1b(alt_dotdata)
        ref_left=[i for i in ref_dotdata if not list(i) in ref_clean_dotdata]
        alt_left=[i for i in alt_dotdata if not list(i) in alt_clean_dotdata]
        [ref_anti_diag_clean_dotdata,ref_anti_diag_kept_segs]=clean_dotdata_anti_diagnal_m1b(ref_left)
        [alt_anti_diag_clean_dotdata,alt_anti_diag_kept_segs]=clean_dotdata_anti_diagnal_m1b(alt_left)
        ref_clean_dotdata+=ref_anti_diag_clean_dotdata
        alt_clean_dotdata+=alt_anti_diag_clean_dotdata
        if len(ref_clean_dotdata)>0 and len(alt_clean_dotdata)>0:
            return [eu_dis_abs_calcu(ref_clean_dotdata),eu_dis_abs_calcu(alt_clean_dotdata)]
        else:    
            return [1,1]
    else:
        return [1,1]
def calcu_valor_single_read_score_directed_dis_m1b(ref_seq,alt_seq,pacbio_read_info,window_size):
    ref_dotdata=dotdata(window_size,pacbio_read_info[0],ref_seq[pacbio_read_info[1]:])
    alt_dotdata=dotdata(window_size,pacbio_read_info[0],alt_seq[pacbio_read_info[1]:])
    if len(ref_dotdata)>0 and len(alt_dotdata)>0:
        [ref_clean_dotdata,ref_kept_segs]=clean_dotdata_diagnal_m1b(ref_dotdata)
        [alt_clean_dotdata,alt_kept_segs]=clean_dotdata_diagnal_m1b(alt_dotdata)
        ref_left=[i for i in ref_dotdata if not list(i) in ref_clean_dotdata]
        alt_left=[i for i in alt_dotdata if not list(i) in alt_clean_dotdata]
        [ref_anti_diag_clean_dotdata,ref_anti_diag_kept_segs]=clean_dotdata_anti_diagnal_m1b(ref_left)
        [alt_anti_diag_clean_dotdata,alt_anti_diag_kept_segs]=clean_dotdata_anti_diagnal_m1b(alt_left)
        ref_clean_dotdata+=ref_anti_diag_clean_dotdata
        alt_clean_dotdata+=alt_anti_diag_clean_dotdata
        if len(ref_clean_dotdata)>0 and len(alt_clean_dotdata)>0:
            return [eu_dis_reg_calcu(ref_clean_dotdata),eu_dis_reg_calcu(alt_clean_dotdata)]
        else:    
            return [1,1]
    else:
        return [1,1]
def calcu_valor_single_read_score_within_10Perc_m1b(ref_seq,alt_seq,pacbio_read_info,window_size):
    ref_dotdata=dotdata(window_size,pacbio_read_info[0],ref_seq[pacbio_read_info[1]:])
    alt_dotdata=dotdata(window_size,pacbio_read_info[0],alt_seq[pacbio_read_info[1]:])
    if len(ref_dotdata)>0 and len(alt_dotdata)>0:
        [ref_clean_dotdata,ref_kept_segs]=clean_dotdata_diagnal_m1b(ref_dotdata)
        [alt_clean_dotdata,alt_kept_segs]=clean_dotdata_diagnal_m1b(alt_dotdata)
        ref_left=[i for i in ref_dotdata if not list(i) in ref_clean_dotdata]
        alt_left=[i for i in alt_dotdata if not list(i) in alt_clean_dotdata]
        [ref_anti_diag_clean_dotdata,ref_anti_diag_kept_segs]=clean_dotdata_anti_diagnal_m1b(ref_left)
        [alt_anti_diag_clean_dotdata,alt_anti_diag_kept_segs]=clean_dotdata_anti_diagnal_m1b(alt_left)
        ref_clean_dotdata+=ref_anti_diag_clean_dotdata
        alt_clean_dotdata+=alt_anti_diag_clean_dotdata
        if len(ref_clean_dotdata)>0 and len(alt_clean_dotdata)>0:
            return [eu_dis_dots_within_10perc(alt_clean_dotdata),eu_dis_dots_within_10perc(ref_clean_dotdata)]    #if prediction correct [large, small]
        else:    
            return [1,1]
    else:
        return [1,1]
def calcu_valor_single_read_score_abs_dis_m2(ref_seq,alt_seq,pacbio_read_info,window_size):
    ref_dotdata=dotdata(window_size,pacbio_read_info[0],ref_seq[pacbio_read_info[1]:])
    alt_dotdata=dotdata(window_size,pacbio_read_info[0],alt_seq[pacbio_read_info[1]:])
    if len(ref_dotdata)>0 and len(alt_dotdata)>0:
        [ref_clean_dotdata,ref_kept_segs]=clean_dotdata_m2(ref_dotdata)
        [alt_clean_dotdata,alt_kept_segs]=clean_dotdata_m2(alt_dotdata)
        if len(ref_clean_dotdata)>0 and len(alt_clean_dotdata)>0:
            return [eu_dis_abs_calcu(ref_clean_dotdata),eu_dis_abs_calcu(alt_clean_dotdata)]
        else:    
            return [1,1]
    else:
        return [1,1]
def cigar2alignstart_by_pos(cigar,align_start,start,end):
    #eg cigar2alignstart(pbam[5],int(pbam[3]),bps)
    import re
    pcigar=re.compile(r'''(\d+)([MIDNSHP=X])''')
    cigars=[]
    for m in pcigar.finditer(cigar):
        cigars.append((m.groups()[0],m.groups()[1]))
    read_rec=0
    align_rec=align_start
    cigar_record=''
    for x in cigars:
        if x[1]=='S':
            read_rec+=int(x[0])
        if x[1]=='M':
            read_rec+=int(x[0])
            align_rec+=int(x[0])
        if x[1]=='D':
            align_rec+=int(x[0])
        if x[1]=='I':
            read_rec+=int(x[0])
        cigar_record=x
        if align_rec>start-1: break
    start_dis=int(align_rec)-start
    if cigar_record[1]=='M':
        new_read_rec=read_rec-start_dis
        new_start_dis=0
        return [new_read_rec,new_start_dis]
    else:
        return [read_rec,start_dis]
def chop_pacbio_read_by_pos(bam_in_new,chrom,start,end,flank_length):
    fbam=os.popen(r'''samtools view %s %s:%d-%d'''%(bam_in_new,chrom,start,end))
    out=[]
    for line in fbam:
        pbam=line.strip().split()
        if not pbam[0]=='@': 
            if int(pbam[3])<start+1:
                align_info=cigar2alignstart_by_pos(pbam[5],int(pbam[3]),start,end)
                align_start=align_info[0]
                miss_bp=align_info[1]
                if not miss_bp>flank_length/2:
                    target_read=pbam[9][align_start:]
                    if len(target_read)>end-start-miss_bp:                    
                        out.append([target_read[:end-start-miss_bp],miss_bp,pbam[0]])
    fbam.close()
    return out
def chromos_readin(ref):
    fin=open(ref+'.fai')
    chromos=[]
    for line in fin:
            pin=line.strip().split()
            chromos.append(pin[0])
    fin.close()
    return chromos
def cluster_range_decide(other_cluster):
    out=[]
    for x in other_cluster:
        out.append([])
        out[-1].append([min(x[0]),max(x[0])])
        out[-1].append([min(x[1]),max(x[1])])
    return out
def cluster_size_decide(range_cluster):
    out=[]
    for x in range_cluster:
        size=(x[0][1]-x[0][0])*(x[1][1]-x[1][0])
        out.append(np.sqrt(size))
    return out
def clean_dotdata_m1(ref_dotdata):
    #cluster dots based on their distance to diagnal
    x_axis=[i[0] for i in ref_dotdata]
    y_axis=[i[1] for i in ref_dotdata]
    dis_to_diagnal=[y_axis[i]-x_axis[i] for i in range(len(x_axis))]
    kept_dot_pos=dis_cluster(dis_to_diagnal,dis_cff=10)
    kept_dot=[]
    for x in kept_dot_pos:
        temp=[[x_axis[i],y_axis[i]] for i in x]
        kept_dot_cluster=dis_cluster([i[0] for i in temp],40)
        for y in kept_dot_cluster:
            kept_dot.append([temp[i] for i in y])
    kept_region=[[x[0],x[1]] for x in kept_dot]
    out=[]
    for x in kept_dot:    out+=x
    return [out,kept_region]
def clean_dotdata_diagnal_m1b(ref_dotdata):
    #cluster dots based on their distance to diagnal
    x_axis=[i[0] for i in ref_dotdata]
    y_axis=[i[1] for i in ref_dotdata]
    dis_to_diagnal=[y_axis[i]-x_axis[i] for i in range(len(x_axis))]
    kept_dot_pos=dis_cluster(dis_to_diagnal,dis_cff=10)
    kept_dot=[]
    for x in kept_dot_pos:
        temp=[[x_axis[i],y_axis[i]] for i in x]
        kept_dot+=temp
    kept_region=[[x[0],x[1]] for x in kept_dot]
    return [kept_dot,kept_region]
def clean_dotdata_anti_diagnal_m1b(ref_dotdata):
    #cluster dots based on their distance to diagnal
    if ref_dotdata==[]:    return [[],[]]
    x_axis=[i[0] for i in ref_dotdata]
    y_axis=[i[1] for i in ref_dotdata]
    dis_to_diagnal=[y_axis[i]+x_axis[i] for i in range(len(x_axis))]
    kept_dot_pos=dis_cluster(dis_to_diagnal,dis_cff=10)
    kept_dot=[]
    for x in kept_dot_pos:
        temp=[[x_axis[i],y_axis[i]] for i in x]
        kept_dot+=temp
    kept_region=[[x[0],x[1]] for x in kept_dot]
    return [kept_dot,kept_region]
def clean_dotdata_m2(ref_dotdata):
    #keep dots closest to diagnal
    data_hash={}
    for x in ref_dotdata:
        if not x[0] in data_hash:    data_hash[x[0]]=x[1]
        else:
            if abs(x[1]-x[0])<abs(data_hash[x[0]]-x[0]):    data_hash[x[0]]=x[1]
    out=[[[i,data_hash[i]] for i in sorted(data_hash.keys())],[]]
    return out
def complementary(seq):
    seq2=[]
    for i in seq:
            if i in 'ATGCN':
                    seq2.append('ATGCN'['TACGN'.index(i)])
            elif i in 'atgcn':
                    seq2.append('atgcn'['tacgn'.index(i)])
    return ''.join(seq2)
def compute_bic(kmeans,X):
    """
    Computes the BIC metric for a given clusters
    Parameters:
    -----------------------------------------
    kmeans:  List of clustering object from scikit learn
    X     :  multidimension np array of data points
    Returns:
    -----------------------------------------
    BIC value
    """
    # assign centers and labels
    centers = [kmeans.cluster_centers_]
    labels  = kmeans.labels_
    #number of clusters
    m = kmeans.n_clusters
    # size of the clusters
    n = np.bincount(labels)
    #size of data set
    N, d = X.shape
    #compute variance for all clusters beforehand
    cl_var=[]
    for i in xrange(m):
        if not n[i] - m==0:
            cl_var.append((1.0 / (n[i] - m)) * sum(distance.cdist(X[np.where(labels == i)], [centers[0][i]], 'euclidean')**2))
        else:
            cl_var.append(float(10**20) * sum(distance.cdist(X[np.where(labels == i)], [centers[0][i]], 'euclidean')**2))
    const_term = 0.5 * m * calcu_log10(N)
    BIC = np.sum([n[i] * calcu_log10(n[i]) -
           n[i] * calcu_log10(N) -
         ((n[i] * d) / 2) * calcu_log10(2*np.pi) -
          (n[i] / 2) * calcu_log10(cl_var[i]) -
         ((n[i] - m) / 2) for i in xrange(m)]) - const_term
    return(BIC)
def dotdata(kmerlen,seq1, seq2):
    nth_base = 1
    inversions = True
    hits = kmerhits(seq1, seq2, kmerlen, nth_base, inversions)
    return hits
def dis_cluster(dis_to_diagnal,dis_cff=10):
    #eg of dis_list=[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 683, 854, 0, 683, 854, 0, 683, 854, 0, 683, 854, 0, 683, 0, 683, 0, 683, 0, 683, 0, 341, 512, 683, 0, 341, 512, 683, 0, 341, 512, 683, 0, 341, 512, 683, 1025, 0, 341, 512, 683, 1025, 0, 170, 341, 512, 683, 1025, 0, 170, 341, 683, 854, 1025, 0, 170, 683, 854, 1025, 0, 683, 1025, 0, 0, 0, 0, 0, 0, 0, 0, 512, 0, 512, 0, 512, 1196, 0, 170, 341, 512, 854, 1196, 1367, 0, 170, 341, 512, 854, 1196, 1367, 0]
    sorted_list=sorted(dis_to_diagnal)
    sub_group1=[[sorted_list[0]]]
    for i in sorted_list[1:]:
        if i-sub_group1[-1][-1]<dis_cff:
            sub_group1[-1].append(i)
        else:
            sub_group1.append([i])
    remove_noise_1=[i for i in sub_group1 if len(i)>50]
    return [[i for i in range(len(dis_to_diagnal)) if dis_to_diagnal[i] in j] for j in remove_noise_1]
def dup_block_combine(dup_block,k1_hap,k2_hap):
    #eg of dup_block=['a', 'b'];     k1_hap='abcd'   ; k2_hap='abab'
    all_combines=[]
    for x in range(len(dup_block)):
        all_combines+=[''.join(list(i)) for i in list(itertools.combinations(dup_block,x+1))]
    all_combines=dup_block_combined_qc(all_combines)
    kept_dup=[]
    for x in all_combines[::-1]:
        if k2_hap.count(x)>1:
            kept_dup.append(x)
    return dup_block_kept_qc(kept_dup)[::-1]
def dup_block_combined_qc(all_combines):
    #eg of all_combines=['a', 'b', 'c', 'd', 'ab', 'ac', 'ad', 'bc', 'bd', 'cd', 'abc', 'abd', 'acd', 'bcd', 'abcd']
    out=[]
    for x in all_combines:
        if len(x)==1:   out.append(x)
        else:
            temp=[ord(i) for i in x]
            if interval_dis_calcu_max(temp)>1: continue
            else:   out.append(x)
    return out
def dup_block_kept_qc(kept_dup):
    #eg of kept_dup:
    out=[]
    if len(kept_dup)>0:
        out.append(kept_dup[0])
        for y in kept_dup[1:]:
            flag_y=0 
            for z in out:
                if y in z:  flag_y+=1
            if flag_y==0:   out.append(y)
    return out
def editDistance(seq1,seq2,cost_matrix,move_matrix,r,c):
    [iCost, dCost , mCost]=[1,1,1]
    if cost_matrix[r][c]==float('inf'):
        if r==0 and c==0:
            cost_matrix[r][c]=0
        elif r==0:
            move_matrix[r][c]=0
            cost_matrix[r][c]=editDistance(seq1,seq2,cost_matrix,move_matrix,r,c-1)+iCost
        elif c==0:
            move_matrix[r][c]=1
            cost_matrix[r][c] = editDistance(seq1,seq2,cost_matrix,move_matrix,r-1,c)+dCost
        else:
            iDist=editDistance(seq1,seq2,cost_matrix,move_matrix,r,c-1) + iCost
            dDist=editDistance(seq1,seq2,cost_matrix,move_matrix,r-1,c) + dCost
            if seq1[r-1] == seq2[c-1]:
                mDist=editDistance(seq1,seq2,cost_matrix,move_matrix,r-1,c-1)
            else:
                mDist=editDistance(seq1,seq2,cost_matrix,move_matrix,r-1,c-1)+mCost
            if iDist<dDist:
                if iDist < mDist:  # insertion is optimal
                    move_matrix[r][c] = 0
                    cost_matrix[r][c] = iDist
                else:    #match is optimal
                    move_matrix[r][c] = 2
                    cost_matrix[r][c] = mDist
            else:
                if dDist < mDist:    #deletion is optimal
                    move_matrix[r][c]=1
                    cost_matrix[r][c]=dDist
                else:    #match is optimal
                    move_matrix[r][c]=2
                    cost_matrix[r][c]=mDist
    return cost_matrix[r][c]
def edit_dis_setup(seq1,seq2):
    cost_matrix=np.full((len(seq1)+1,len(seq2)+1),np.inf)
    move_matrix=np.full((len(seq1)+1,len(seq2)+1),-1.0)
    opt_dist=editDistance(seq1, seq2, cost_matrix,move_matrix, len(seq1), len(seq2))
    return opt_dist
def eu_dis_abs_calcu(ref_clean_dotdata):
    #eg of ref_clean_dotdata=[[2380, 471], [2381, 472], [2382, 472], [2383, 472], [2384, 472], [2385, 472], [2386, 472], [2387, 472], [2388, 472], [2389, 473], [2390, 474], [2391, 475], [2402, 486], [2403, 487], [2404, 488], [2405, 489], [2406, 490], [2407, 491], [2408, 492], [2409, 493], [2410, 494], [2411, 495], [2418, 504], [2440, 528], [2441, 529], [2442, 530], [2443, 531], [2454, 542], [2455, 543], [2456, 544], [2473, 563], [2474, 564], [2475, 565], [2476, 566], [2477, 567], [2478, 568], [2479, 569], [2480, 570], [2481, 571], [2482, 572], [2483, 573], [2484, 574], [2485, 575], [2486, 576], [2487, 577], [2488, 578], [2489, 579], [2490, 580], [2500, 589], [2501, 590], [2502, 591], [2503, 592], [2504, 593], [2523, 613], [2524, 614], [2525, 615], [2526, 616], [2527, 617], [2528, 618], [2529, 619], [2530, 620], [2531, 621], [2532, 622], [2533, 623], [2557, 642], [2558, 643], [2559, 644], [2560, 645], [2561, 646], [2562, 647]]
    eu_dis_abs=[abs(i[0]-i[1]) for i in ref_clean_dotdata]
    return np.mean(eu_dis_abs)
def eu_dis_reg_calcu(ref_clean_dotdata):
    #eg of ref_clean_dotdata=[[2380, 471], [2381, 472], [2382, 472], [2383, 472], [2384, 472], [2385, 472], [2386, 472], [2387, 472], [2388, 472], [2389, 473], [2390, 474], [2391, 475], [2402, 486], [2403, 487], [2404, 488], [2405, 489], [2406, 490], [2407, 491], [2408, 492], [2409, 493], [2410, 494], [2411, 495], [2418, 504], [2440, 528], [2441, 529], [2442, 530], [2443, 531], [2454, 542], [2455, 543], [2456, 544], [2473, 563], [2474, 564], [2475, 565], [2476, 566], [2477, 567], [2478, 568], [2479, 569], [2480, 570], [2481, 571], [2482, 572], [2483, 573], [2484, 574], [2485, 575], [2486, 576], [2487, 577], [2488, 578], [2489, 579], [2490, 580], [2500, 589], [2501, 590], [2502, 591], [2503, 592], [2504, 593], [2523, 613], [2524, 614], [2525, 615], [2526, 616], [2527, 617], [2528, 618], [2529, 619], [2530, 620], [2531, 621], [2532, 622], [2533, 623], [2557, 642], [2558, 643], [2559, 644], [2560, 645], [2561, 646], [2562, 647]]
    eu_dis_abs=[i[0]-i[1] for i in ref_clean_dotdata]
    return abs(np.mean(eu_dis_abs))
def eu_dis_dots_within_10perc(ref_clean_dotdata):
    #eg of ref_clean_dotdata=[[2380, 471], [2381, 472], [2382, 472], [2383, 472], [2384, 472], [2385, 472], [2386, 472], [2387, 472], [2388, 472], [2389, 473], [2390, 474], [2391, 475], [2402, 486], [2403, 487], [2404, 488], [2405, 489], [2406, 490], [2407, 491], [2408, 492], [2409, 493], [2410, 494], [2411, 495], [2418, 504], [2440, 528], [2441, 529], [2442, 530], [2443, 531], [2454, 542], [2455, 543], [2456, 544], [2473, 563], [2474, 564], [2475, 565], [2476, 566], [2477, 567], [2478, 568], [2479, 569], [2480, 570], [2481, 571], [2482, 572], [2483, 573], [2484, 574], [2485, 575], [2486, 576], [2487, 577], [2488, 578], [2489, 579], [2490, 580], [2500, 589], [2501, 590], [2502, 591], [2503, 592], [2504, 593], [2523, 613], [2524, 614], [2525, 615], [2526, 616], [2527, 617], [2528, 618], [2529, 619], [2530, 620], [2531, 621], [2532, 622], [2533, 623], [2557, 642], [2558, 643], [2559, 644], [2560, 645], [2561, 646], [2562, 647]]
    dis_10perc=[abs(float(i[0]-i[1])/float(i[0])) for i in ref_clean_dotdata if i[0]>0]
    return len([i for i in dis_10perc if i <0.1])    
def flank_length_calculate(bps):
    if int(bps[-1])-int(bps[1])<100:
        flank_length=(int(bps[-1])-int(bps[1]))
    else:
        if int(bps[-1])-int(bps[1])<500:
            flank_length=int(bps[-1])-int(bps[1])
        else:
            flank_length=500
    return flank_length
def genoCN_extract(pin):
    out=[0,0]
    rec_pos=-1
    if 'CN' in pin[8]:
        for x in pin[8].split(':'):
            rec_pos+=1
            if x=='CN': break
    geno=[i.split(':')[rec_pos] for i in pin[9:]]
    out=[0 if i=='2' else 1 for i in geno]
    return out
def genotype_extract(pin):
    out=[0,0]
    rec_pos=-1
    if 'GT' in pin[8]:
        for x in pin[8].split(':'):
            rec_pos+=1
            if x=='GT': break
    geno=[i.split(':')[rec_pos] for i in pin[9:]]
    for i in geno:
        if '/' in i:
            if geno=='./.': out.append(1)
            else:   out.append(sum([int(j) for j in i.split('/')]))
        elif '|' in i:
            if geno=='.|.': out.append(1)
            else:   out.append(sum([int(j) for j in i.split('|')]))
        elif i=='.':    out.append(1)
    return out
def INS_length_detect(pin):
    out=0
    for x in pin[7].split(';'):
        if 'SVLEN=' in x:
            out=int(x.split('=')[1])
    return out
def intersect(a, b):
    return ''.join(sorted(list(set(a) & set(b))))
def interval_dis_calcu_max(pos_check):
    #eg of pos_check=[97,98]
    if len(pos_check)>1:
        out=[pos_check[i+1]-pos_check[i] for i in range(len(pos_check)-1)]
        return max(out)
    else:
        return 'NA'
def k_means_cluster(data_list):
    if max(data_list[0])-min(data_list[0])>10 and max(data_list[1])-min(data_list[1])>10:
        array_diagnal=np.array([[data_list[0][x],data_list[1][x]] for x in range(len(data_list[0]))])
        ks = range(1,min([5,len(data_list[0])+1]))
        KMeans = [cluster.KMeans(n_clusters = i, init="k-means++").fit(array_diagnal) for i in ks]
        KMeans_predict=[cluster.KMeans(n_clusters = i, init="k-means++").fit_predict(array_diagnal) for i in ks]
        BIC=[]
        BIC_rec=[]
        for x in ks:
            if KMeans_predict[x-1].max()<x-1: continue
            else:
                BIC_i=compute_bic(KMeans[x-1],array_diagnal)
                if abs(BIC_i)<10**8:
                    BIC.append(BIC_i)
                    BIC_rec.append(x)
        #BIC = [compute_bic(kmeansi,array_diagnal) for kmeansi in KMeans]
        #ks_picked=ks[BIC.index(max(BIC))]
        ks_picked=BIC_rec[BIC.index(max(BIC))]
        if ks_picked==1:
            return [data_list]
        else:
            out=[]
            std_rec=[scipy.std(data_list[0]),scipy.std(data_list[1])]
            whitened = whiten(array_diagnal)
            centroids, distortion=kmeans(whitened,ks_picked)
            idx,_= vq(whitened,centroids)
            for x in range(ks_picked):
                group1=[[int(i) for i in array_diagnal[idx==x,0]],[int(i) for i in array_diagnal[idx==x,1]]]
                out.append(group1)
            return out
    else:
        return [data_list]
def k_means_cluster_Predict(data_list,info):
    array_diagnal=np.array([[data_list[0][x],data_list[1][x]] for x in range(len(data_list[0]))])
    ks = range(1,len(info))
    KMeans = [cluster.KMeans(n_clusters = i, init="k-means++").fit(array_diagnal) for i in ks]
    BIC = [compute_bic(kmeansi,array_diagnal) for kmeansi in KMeans]
    ks_picked=ks[BIC.index(max(BIC))]
    if ks_picked==1:
        return [data_list]
    else:
        out=[]
        std_rec=[scipy.std(data_list[0]),scipy.std(data_list[1])]
        whitened = whiten(array_diagnal)
        centroids, distortion=kmeans(whitened,ks_picked)
        idx,_= vq(whitened,centroids)
        for x in range(ks_picked):
            group1=[[int(i) for i in array_diagnal[idx==x,0]],[int(i) for i in array_diagnal[idx==x,1]]]
            out.append(group1)
        return out
def key_modify(key):
    if 'R' in key:
            key=key.replace('R','N')
    if 'r' in key:
            key=key.replace('r','n')
    if 'Y' in key:
            key=key.replace('Y','N')
    if 'y' in key:
            key=key.replace('y','n')
    if 'S' in key:
            key=key.replace('S','N')
    if 's' in key:
            key=key.replace('s','n')
    if 'W' in key:
            key=key.replace('W','N')
    if 'w' in key:
            key=key.replace('w','n')
    if 'K' in key:
            key=key.replace('K','N')
    if 'k' in key:
            key=key.replace('k','n')
    if 'M' in key:
            key=key.replace('M','N')
    if 'm' in key:
            key=key.replace('m','n')
    if 'B' in key:
            key=key.replace('B','N')
    if 'b' in key:
            key=key.replace('b','n')
    if 'D' in key:
            key=key.replace('D','N')
    if 'd' in key:
            key=key.replace('d','n')
    if 'H' in key:
            key=key.replace('H','N')
    if 'h' in key:
            key=key.replace('h','n')
    if 'V' in key:
            key=key.replace('V','N')
    if 'v' in key:
            key=key.replace('v','n')
    return key
def kmerhits(seq1, seq2, kmerlen, nth_base=1, inversions=False):
    # hash table for finding hits
    lookup = {}
    # store sequence hashes in hash table
    #print "hashing seq1..."
    seq1len = len(seq1)
    for i in xrange(seq1len - kmerlen + 1):
        key = seq1[i:i+kmerlen]
        for subkey in subkeys(key, nth_base, inversions):
            lookup.setdefault(subkey, []).append(i)
    # match every nth base by look up hashes in hash table
    #print "hashing seq2..."
    hits = []
    for i in xrange(len(seq2) - kmerlen + 1):
        key = seq2[i:i+kmerlen]
        # only need to specify inversions for one seq
        for subkey in subkeys(key, nth_base, False):
            subhits = []
            if kmerlen>40:
                print 'Window size:' +str(kmerlen)
                for k1 in lookup.keys():
                    if edit_dis_setup(k1,subkey)<int(kmerlen/10)+1:
                        subhits+=lookup[k1]
            else:
                subhits=lookup.get(subkey,[])
            if subhits != []:
                # store hits to hits list
                for hit in subhits:
                    hits.append((i, hit))
                # break out of loop to avoid doubly counting
                # exact matches
                break
    return hits
def let_to_block_info(let,let_hash,chromos):
    #eg of let='ab'; eg of let_hash={'a': ['chrY', '10818935', '10819073'], 'b': ['chrY', '10819073', '10926507'], '+': ['chrY', '10926507', '10927007'], '-': ['chrY', '10818435', 10818935]}
    out=[]
    for i in let:
        if not i=='^':
            out+=let_hash[i]
    return(block_modify(out,chromos))
def letter_subgroup(k2_hap):
    #eg of k2_hap='ac^b^'
    inverted_sv=[]
    for x in k2_hap:
        if not x=='^':  inverted_sv.append(x)
        else:   inverted_sv[-1]+='^'
    inverted_sv_2=[]
    for x in inverted_sv:
        if inverted_sv_2==[]: inverted_sv_2.append(x)
        else:
            if not '^' in inverted_sv_2[-1] and not '^' in x and ord(x)-ord(inverted_sv_2[-1][-1])==1:  inverted_sv_2[-1]+=x
            elif '^' in inverted_sv_2[-1] and '^'  in x and ord(x[0])-ord(inverted_sv_2[-1][-2])==-1:   inverted_sv_2[-1]+=x
            else:   inverted_sv_2.append(x)
    inverted_sv_3=[]
    for i in inverted_sv_2:
        if not '^' in i:    inverted_sv_3.append(i)
        else:
            inverted_sv_3.append(i.replace('^','')[::-1]+'^')
    return inverted_sv_3
def letter_split(let):
    #eg of let='c^ba'
    out=[]
    for x in let:
        if not x=='^':  out.append(x)
        else:   out[-1]+=x
    return out
def list_unify(list):
    out=[]
    for i in list:
        if not i in out:    out.append(i)
    return out
def make_event_figure_1(valor_score_list,best_read_rec,window_size,ref_seq,alt_seq,out_figure_name):
    if not best_read_rec=='':
        if not best_read_rec==[]:
            dotdata_record=[dotdata(window_size,ref_seq,ref_seq),    dotdata(window_size,alt_seq,alt_seq),    dotdata(window_size,best_read_rec[0],ref_seq[best_read_rec[1]:]),    dotdata(window_size,best_read_rec[0],alt_seq[best_read_rec[1]:])]
            grdevices = importr('grDevices')
            graphics = importr("graphics")
            if  len(out_figure_name.split('/')[-1])>150:
                out_figure_name='/'.join(out_figure_name.split('/')[:-1])+'/'+out_figure_name.split('/')[-1][:140]+'.'+out_figure_name.split('.')[-1]
            print out_figure_name
            grdevices.png(file=out_figure_name, width=512, height=512)
            graphics.par(mfrow=IntVector([2,2]))
            #penal 1
            graphics.plot([min([i[0] for i in dotdata_record[0]]),max([i[0] for i in dotdata_record[0]])],[min([i[1] for i in dotdata_record[0]]),max([i[1] for i in dotdata_record[0]])],type='n',xlab='oroginal_reference',ylab='oroginal_reference')
            graphics.points([i[0] for i in dotdata_record[0]],[i[1] for i in dotdata_record[0]],col='red',pch='+')
            #penal 2
            graphics.plot([min([i[0] for i in dotdata_record[1]]),max([i[0] for i in dotdata_record[1]])],[min([i[1] for i in dotdata_record[1]]),max([i[1] for i in dotdata_record[1]])],type='n',xlab='altered_reference',ylab='altered_reference')
            graphics.points([i[0] for i in dotdata_record[1]],[i[1] for i in dotdata_record[1]],col='red',pch='+')
            #penal 3
            graphics.plot([min([i[0] for i in dotdata_record[2]]),max([i[0] for i in dotdata_record[2]])],[min([i[1] for i in dotdata_record[2]]),max([i[1] for i in dotdata_record[2]])],type='n',ylab='PacBio read',xlab='oroginal_reference')
            graphics.points([i[0] for i in dotdata_record[2]],[i[1] for i in dotdata_record[2]],col='red',pch='+')
            #penal 4
            graphics.plot([min([i[0] for i in dotdata_record[3]]),max([i[0] for i in dotdata_record[3]])],[min([i[1] for i in dotdata_record[3]]),max([i[1] for i in dotdata_record[3]])],type='n',ylab='PacBio read',xlab='altered_reference')
            graphics.points([i[0] for i in dotdata_record[3]],[i[1] for i in dotdata_record[3]],col='red',pch='+')
            grdevices.dev_off()
def minimize_pacbio_read_list(pacbio_read_info,ideal_list_length=20):
    #eg of pacbio_read_info=chop_pacbio_read_by_pos(bam_in_new,chrom,start,end,flank_length)
    if len(pacbio_read_info)>ideal_list_length:
        out=[]
        temp_hash={}
        for x in pacbio_read_info:
            if not x[1] in temp_hash.keys():    temp_hash[x[1]]=[x]
            else:    temp_hash[x[1]]+=[x]
        for x in sorted(temp_hash.keys()):
            if len(out)<ideal_list_length:    out+=temp_hash[x]
            return out[:ideal_list_length]
    else: return pacbio_read_info
def reverse(seq):
    return seq[::-1]
def path_mkdir(path):
    if not os.path.isdir(path):
        os.system(r'''mkdir %s'''%(path))
def path_modify(path):
    if not path[-1]=='/':
        path+='/'
    return path
def polarity_detect(pin):
    out='+'
    for x in pin[7].split(';'):
        if 'MEIINFO=' in x:
            out=x.split(',')[-1]
    return out
def qual_check_repetitive_region(dotdata_qual_check):
    #qual_check_R_2(dotdata_qual_check,txt_file)  #not necessary for validator
    diagnal=0
    other=[[],[]]
    for x in dotdata_qual_check:
        if x[0]==x[1]:
            diagnal+=1 
        else:
            if x[0]>x[1]:
                other[0].append(x[0])
                other[1].append(x[1])
    if len(dotdata_qual_check)>0 and float(len(other[0]))/float(len(dotdata_qual_check))>0.1 and float(len(other[0]))/float(len(dotdata_qual_check))<0.5:
        other_cluster=X_means_cluster_reformat(other)
        range_cluster=cluster_range_decide(other_cluster)
        size_cluster=cluster_size_decide(range_cluster)
    else:
        size_cluster=[0]
    return [float(diagnal)/float(len(dotdata_qual_check)),size_cluster]
def simple_del_diploid_decide(k1,k2):
    #eg of k1='ab/ab'   ; eg of k2='a/a'
    k2_haps=k2.split('/')
    k1_hap=k1.split('/')[0]
    out=[]
    for x in k2_haps:
        if x==k1_hap: out.append('NA')
        else:
            out.append(simple_del_haploid_decide(k1_hap,x))
    return out
def simple_del_haploid_decide(k1_hap,k2_hap):
    #eg of k1_hap='ab'  ;    eg of k2_hap='b'
    if k1_hap==k2_hap: return 'FALSE'   #no alt
    if k2_hap=='': return [i for i in k1_hap]
    if '^' in k2_hap:   return 'FALSE'  #check if inv included
    dup_test=[k2_hap.count(x) for x in k2_hap]
    if max(dup_test)>1:     return 'FALSE'  #check if dup included
    if len(k2_hap)==1 and len(k1_hap)>1:    return letter_subgroup(''.join([i for i in k1_hap if not i in k2_hap]))   #del
    pos_compare=[ord(k2_hap[i+1])-ord(k2_hap[i]) for i in range(len(k2_hap)-1)]
    if min(pos_compare)<1: return 'FALSE'
    return letter_subgroup(''.join([i for i in k1_hap if not i in k2_hap]))
def simple_inv_diploid_decide(k1,k2):
    #eg of k1='ab/ab'   ; eg of k2='ab^/ab'
    k2_haps=k2.split('/')
    k1_hap=k1.split('/')[0]
    out=[]
    for x in k2_haps:
        if x==k1_hap: out.append('NA')
        else:
            out.append(simple_inv_haploid_decide(k1_hap,x))
    return out
def simple_inv_haploid_decide(k1_hap,k2_hap):
    #eg of k1_hap='ab'  ;    eg of k2_hap='b^a^'
    if not '^' in k2_hap:   return 'FALSE'      #if not block inverted
    if len(k2_hap.replace('^',''))==1 and len(k1_hap)==1:   return [i for i in k1_hap]
    dup_test=[k2_hap.count(i) for i in k2_hap if not i=='^']
    if max(dup_test)>1: return 'FALSE'
    inverted_sv_new=letter_subgroup(k2_hap)
    if ''.join([i.replace('^','') for i in inverted_sv_new])==k1_hap: return [i[:-1] for i in inverted_sv_new if '^' in i]
    else:   return 'FALSE'
def simple_tandup_diploid_decide(k1,k2):
    #eg of k1='ab/ab'   ; eg of k2='abb/ab'
    k2_haps=k2.split('/')
    k1_hap=k1.split('/')[0]
    out=[]
    for x in k2_haps:
        if x==k1_hap: out.append('NA')
        else:
            out.append(simple_tandup_haploid_decide(k1_hap,x))
    return out
def simple_tandup_haploid_decide(k1_hap,k2_hap):
    if '^' in k2_hap:   return 'FALSE'
    dup_count=[k2_hap.count(i) for i in k1_hap]
    if min(dup_count)<1 or max(dup_count)<2:    return 'FALSE'  #deletion structure inside
    out=[]
    temp1=[]
    for x in k2_hap:
        if temp1==[]:   temp1.append(x)
        elif ord(x)-ord(temp1[-1][-1])==1:  temp1[-1]+=x
        else:   temp1.append(x)
    overlap_portion=[]
    overlap_count=[]
    for x in temp1:
        if out==[]:
            out.append(x)
        else:
            overlap=intersect(out[-1],x)
            if not len(overlap) >len(out[-1]) and not len(overlap)>len(x):
                if out[-1][-len(overlap):]==x[:len(overlap)]:
                    out[-1]+=x[len(overlap):]
                    if not overlap in overlap_portion:
                        overlap_portion.append(overlap)
                        overlap_count.append(2)
                    else:
                        overlap_count[overlap_portion.index(overlap)]+=1
                else:
                    out.append(x)
            else:
                out.append(x)
    if ''.join(out)==k1_hap:
        return [overlap_portion,overlap_count]
    return 'FALSE'
def simple_disdup_diploid_decide(k1,k2):
    #eg of k1='ab/ab'   ; eg of k2='bab/ab'
    k2_haps=k2.split('/')
    k1_hap=k1.split('/')[0]
    out=[]
    for x in k2_haps:
        if x==k1_hap: out.append('NA')
        else:
            out.append(simple_disdup_haploid_decide(k1_hap,x))
    return out
def simple_disdup_haploid_decide(k1_hap,k2_hap):
    #eg of k1_hap='abcd'    ;   eg of k2_hap='babdcd'
    if not '^' in k2_hap:
        if simple_tandup_haploid_decide(k1_hap,k2_hap)=='FALSE': 
            dup_dis=letter_subgroup(k2_hap)
            overlap=[intersect(dup_dis[i],dup_dis[i+1]) for i in range(len(dup_dis)-1)]
            if len(list_unify(overlap))==len(overlap):
                dup_count=[k2_hap.count(i) for i in k1_hap]
                if not min(dup_count)<1 and not max(dup_count)<2:    #deletion structure inside
                    dup_block=[k1_hap[i] for i in range(len(dup_count)) if dup_count[i]>1]
                    dup_block_combined=dup_block_combine(dup_block,k1_hap,k2_hap)
                    dis_dup_check=[]
                    no_dup_block=[]
                    for x in k2_hap: 
                        if not x in dup_block:
                            no_dup_block.append(k2_hap.index(x))
                    for x in dup_block_combined:
                        dis_dup_check.append([])
                        for y in range(len(k2_hap)-len(x)+1):
                            if k2_hap[y:(y+len(x))]==x:
                                dis_dup_check[-1].append(y)
                    original_pos=[]
                    for x in itertools.product(*dis_dup_check):
                        x_modify_new=x_to_x_modify_new(x,dup_block_combined)
                        temp_structure=[k2_hap[i] for i in sorted(x_modify_new+no_dup_block)]
                        if ''.join(temp_structure)==k1_hap:
                            original_pos+=list(x)
                    if len(original_pos)>0:
                        insert_pos=[]
                        for i in dis_dup_check:
                            for j in i:
                                if not j in original_pos:
                                    insert_pos.append(j)
                        k2_hap_new=['-']+[i for i in k2_hap]+['+']
                        insert_block=[]
                        pos_rec=-1
                        for i in insert_pos:
                            pos_rec+=1
                            if len(dup_block_combined[pos_rec])==1:
                                insert_block.append([k2_hap_new[i],k2_hap_new[i+1],k2_hap_new[i+2]])
                            else:
                                insert_block.append([k2_hap_new[i]]+k2_hap_new[(i+1):(i+len(dup_block_combined[pos_rec])+2)])
                        #insert_block=[[k2_hap_new[i],k2_hap_new[i+1],k2_hap_new[i+2]] for i in insert_pos]
                        return [dup_block_combined,insert_block]
    return 'FALSE'
def simple_del_chop_pacbio_read_simple_short(bam_in,sv_info,flank_length):
    #eg of info = ['a/a','/','chr1', 101553562, 101553905]
    #pick read around left breakpoints
    #block_length={}
    #for x in range(len(info)-4):
    #    block_length[chr(97+x)]=int(info[x+4])-int(info[x+3])
    bam_in_new_list=bam_in_decide(bam_in,sv_info)
    if bam_in_new_list=='': return [[],[],[]]
    pacbio_read_info=[]
    for bam_in_new in bam_in_new_list:
        pacbio_read_info+=chop_pacbio_read_by_pos(bam_in_new, sv_info[0],int(sv_info[1])-flank_length,int(sv_info[1])+flank_length,flank_length)
    pacbio_read_info=minimize_pacbio_read_list(pacbio_read_info)
    return pacbio_read_info
def simple_chop_pacbio_read_simple_short(bam_in,sv_info,flank_length):
    #eg of info = ['a/a','/','chr1', 101553562, 101553905]
    #eg of sv_info=['chr1', 101553562, 101553905]
    bam_in_new_list=bam_in_decide(bam_in,sv_info)
    if bam_in_new_list=='': return [[],[],[]]
    pacbio_read_info=[]
    for bam_in_new in bam_in_new_list:
        pacbio_read_info+=chop_pacbio_read_by_pos(bam_in_new, sv_info[0],int(sv_info[1])-flank_length,int(sv_info[-1])+flank_length,flank_length)
    pacbio_read_info=minimize_pacbio_read_list(pacbio_read_info)
    return pacbio_read_info
def del_inv_Valor(bam_in,ref,sv_info,out_figure_name):
    sv_block=[sv_info[0][0],sv_info[0][1],sv_info[-1][2]]
    flank_length=flank_length_calculate(sv_block)
    valor_score_list=[]
    best_read_rec=''
    if sv_info[1][1]-sv_info[0][2]<100:
        if sv_block[2]-sv_block[1]<default_max_sv_test: #only try to read in all reads with sv <100K; else: try breakpoints ; 
            ref_seq=ref_seq_readin(ref,sv_block[0],sv_block[1]-flank_length,sv_block[2]+flank_length)
            [window_size,window_size_qc]=window_size_refine(ref_seq)
            if not window_size=='Error':
                alt_seq=ref_seq[:flank_length]
                for x in sv_info:
                    if x[-1]=='del':continue
                    elif x[-1]=='inv':  alt_seq+=reverse(complementary(ref_seq_readin(ref,x[0],x[1],x[2])))
                alt_seq+=ref_seq[-flank_length:]
                [window_size,window_size_qc]=window_size_refine(alt_seq)
                if not window_size=='Error':
                    all_reads=simple_chop_pacbio_read_simple_short(bam_in,sv_block,flank_length)
                    if len(all_reads)>10:
                        for x in all_reads:
                            valor_single_read_score=calcu_valor_single_read_score_abs_dis_m1b(ref_seq,alt_seq,x,window_size)
                            if not valor_single_read_score==[1,1] and not valor_single_read_score[0]==0:
                                valor_score_list.append(1-float(valor_single_read_score[1])/float(valor_single_read_score[0]))
                                if valor_score_list[-1]==max(valor_score_list):best_read_rec=x
                        make_event_figure_1(valor_score_list,best_read_rec,window_size,ref_seq,alt_seq,out_figure_name)
                    else:
                        if len(sv_info)==2 and [i[-1] for i in sv_info]==['del','inv']:
                            valor_score_list=long_del_inv(bam_in,ref,sv_info,out_figure_name)
        else:
            if len(sv_info)==2 and [i[-1] for i in sv_info]==['del','inv']:
                valor_score_list=long_del_inv(bam_in,ref,sv_info,out_figure_name)
    else:
        for sub_sv_info in sv_info:
            if 'del' in sub_sv_info:    valor_score_list+=simple_del_Valor(bam_in,ref,sub_sv_info[:-1],'.'.join(out_figure_name.split('.')[:-1])+'_'.join([str(i) for i in sub_sv_info])+'.'+out_figure_name.split('.')[-1])
            elif 'inv' in sub_sv_info:    valor_score_list+=simple_inv_Valor(bam_in,ref,sub_sv_info[:-1],'.'.join(out_figure_name.split('.')[:-1])+'_'.join([str(i) for i in sub_sv_info])+'.'+out_figure_name.split('.')[-1])
    return valor_score_list
def dup_inv_ValoR(bam_in,ref,sv_info,out_figure_name):
    #eg of sv_info=['chr1', 114103333, 114103408, 'chr1', 114111746]
    sv_info[1:3]=[int(i) for i in sv_info[1:3]]
    dup_block=sv_info[:3]
    ins_point=[sv_info[3],int(sv_info[4])]
    flank_length=flank_length_calculate(dup_block)
    valor_score_list=[]
    best_read_rec=''
    bp_info=sorted(sv_info[1:3]+[sv_info[4]])
    run_flag=0
    if sv_info[0]==sv_info[3] and max(bp_info)-min(bp_info)<default_max_sv_test: #try to test on the whole region
        ref_seq=ref_seq_readin(ref,sv_info[0],min(bp_info)-flank_length,max(bp_info)+flank_length)
        [window_size,window_size_qc]=window_size_refine(ref_seq)
        if not window_size=='Error':
            run_flag+=1
            ref_structure=[i for i in 'ab']
            if sv_info[4]>sv_info[2]:           alt_structure=['a','b','a^']
            elif sv_info[4]<sv_info[1]:         alt_structure=['b^','a','b']
            else:                               alt_structure=['a','a^']
            all_reads=simple_chop_pacbio_read_simple_short(bam_in,[sv_info[0]]+bp_info+[bp_info[-1]+sv_info[2]-sv_info[1]],flank_length)
            if len(all_reads)>0:
                alt_seq=ref_seq_readin(ref,sv_info[0],min(bp_info)-flank_length,min(bp_info))
                a_seq=ref_seq_readin(ref,sv_info[0],bp_info[0],bp_info[1])
                b_seq=ref_seq_readin(ref,sv_info[0],bp_info[1],bp_info[2])
                for x in alt_structure:
                    if x=='a':  alt_seq+=a_seq
                    elif x=='a^':   alt_seq+=reverse(complementary(a_seq))
                    elif x=='b':    alt_seq+=b_seq
                    elif x=='b^':   alt_seq+=reverse(complementary(b_seq))
                alt_seq+=ref_seq_readin(ref,sv_info[0],max(bp_info),max(bp_info)+flank_length)
                [window_size,window_size_qc]=window_size_refine(alt_seq)
                if not window_size=='Error':
                    for x in all_reads:
                        valor_single_read_score=calcu_valor_single_read_score_directed_dis_m1b(ref_seq,alt_seq,x,window_size)
                        if not valor_single_read_score==[1,1] and not valor_single_read_score[0]==0:
                            valor_score_list.append(1-float(valor_single_read_score[1])/float(valor_single_read_score[0]))
                            if valor_score_list[-1]==max(valor_score_list):best_read_rec=x
                    make_event_figure_1(valor_score_list,best_read_rec,window_size,ref_seq,alt_seq,out_figure_name)
    if run_flag==0:
        if max(bp_info)-min(bp_info)<default_max_sv_test: #only try to read in all reads with sv <100K; else: try breakpoints ; 
            ref_seq=ref_seq_readin(ref,ins_point[0],ins_point[1]-flank_length,ins_point[1]+flank_length)
            [window_size,window_size_qc]=window_size_refine(ref_seq)
            if not window_size=='Error':
                all_reads=simple_del_chop_pacbio_read_simple_short(bam_in,ins_point,flank_length)
                if len(all_reads)>0:
                    alt_seq=ref_seq[:flank_length]+reverse(complementary(ref_seq_readin(ref,dup_block[0],dup_block[1],dup_block[2])))+ref_seq[-flank_length:]
                    [window_size,window_size_qc]=window_size_refine(alt_seq)
                    if not window_size=='Error':
                        for x in all_reads:
                            valor_single_read_score=calcu_valor_single_read_score_abs_dis_m1b(ref_seq,alt_seq,x,window_size)
                            if not valor_single_read_score==[1,1] and not valor_single_read_score[0]==0:
                                valor_score_list.append(1-float(valor_single_read_score[1])/float(valor_single_read_score[0]))
                                if valor_score_list[-1]==max(valor_score_list):best_read_rec=x
                        make_event_figure_1(valor_score_list,best_read_rec,window_size,ref_seq,alt_seq,out_figure_name)
        else:
            ref_seq=ref_seq_readin(ref,ins_point[0],ins_point[1]-flank_length,ins_point[1]+flank_length)
            [window_size,window_size_qc]=window_size_refine(ref_seq)
            if not window_size=='Error':
                all_reads=simple_del_chop_pacbio_read_simple_short(bam_in,ins_point,flank_length)
                if len(all_reads)>0:
                    alt_seq=ref_seq[:flank_length]+reverse(complementary(ref_seq_readin(ref,dup_block[0],dup_block[2]-flank_length,dup_block[2])))
                    [window_size,window_size_qc]=window_size_refine(alt_seq)
                    if not window_size=='Error':
                        for x in all_reads:
                            valor_single_read_score=calcu_valor_single_read_score_within_10Perc_m1b(ref_seq,alt_seq,x,window_size)
                            if not valor_single_read_score==[1,1] and not valor_single_read_score[0]==0:
                                valor_score_list.append(1-float(valor_single_read_score[1])/float(valor_single_read_score[0]))
                                if valor_score_list[-1]==max(valor_score_list):best_read_rec=x
                        make_event_figure_1(valor_score_list,best_read_rec,window_size,ref_seq,alt_seq,out_figure_name)
    return valor_score_list
def CANNOT_CLASSIFY_ValoR(bam_in,ref,sv_info,out_figure_name):
    ref_sv=sv_info[0].split('_')
    alt_sv=list_unify([i for i in sv_info[1].split('_') if not i in ref_sv])
    chromos=chromos_readin(ref)
    bp_info=block_subsplot(sv_info[2:],chromos)
    flank_length=max([flank_length_calculate(i) for i in bp_info])
    valor_score_list=[]
    best_read_rec=''
    run_flag=0
    if len(bp_info)==1: #inter-chromosome event not considered here
        if bp_info[0][-1]-bp_info[0][1]<default_max_sv_test:
            ref_seq=ref_seq_readin(ref,bp_info[0][0],bp_info[0][1]-flank_length,bp_info[0][-1]+flank_length)
            [window_size,window_size_qc]=window_size_refine(ref_seq)
            if not window_size=='Error':
                all_reads=simple_chop_pacbio_read_simple_short(bam_in,bp_info[0],flank_length)
                bp_let_hash=bp_to_chr_hash(bp_info[0],chromos,flank_length)
                if len(all_reads)>10:
                    run_flag+=1
                    bp_let_seq={}
                    for i in bp_let_hash.keys():
                        bp_let_seq[i]=ref_seq_readin(ref,bp_let_hash[i][0],int(bp_let_hash[i][1]),int(bp_let_hash[i][-1]))
                    for alt_allele in alt_sv:
                        alt_seq=ref_seq[:flank_length]
                        for i in letter_split(alt_allele):
                            if not '^' in i:    alt_seq+=bp_let_seq[i]
                            else:               alt_seq+=reverse(complementary(bp_let_seq[i[0]]))
                        alt_seq+=ref_seq[-flank_length:]
                        [window_size,window_size_qc]=window_size_refine(alt_seq)
                        if not window_size=='Error':
                            for x in all_reads:
                                if not max([alt_allele.count(i) for i in alt_allele]+[0])>1:    valor_single_read_score=calcu_valor_single_read_score_abs_dis_m1b(ref_seq,alt_seq,x,window_size)
                                else:                                                       valor_single_read_score=calcu_valor_single_read_score_directed_dis_m1b(ref_seq,alt_seq,x,window_size)
                                if not valor_single_read_score==[1,1] and not valor_single_read_score[0]==0:
                                    valor_score_list.append(1-float(valor_single_read_score[1])/float(valor_single_read_score[0]))
                                    if valor_score_list[-1]==max(valor_score_list):best_read_rec=x
                            make_event_figure_1(valor_score_list,best_read_rec,window_size,ref_seq,alt_seq,'.'.join(out_figure_name.split('.')[:-1]+[ref_sv[0]+'.vs.'+alt_allele,out_figure_name.split('.')[-1]]))
        if run_flag==0:
            for alt_allele in alt_sv:
                alt_juncs=block_around_check(alt_allele,ref_sv[0])
                bp_let_hash=bp_to_chr_hash(bp_info[0],chromos,flank_length)
                for alt_jun in alt_juncs:
                    print alt_jun
                    if not '^' in alt_jun[0]:                        ref_seq_a=ref_seq_readin(ref,bp_let_hash[alt_jun[0][0]][0],bp_let_hash[alt_jun[0][0]][2]-flank_length,bp_let_hash[alt_jun[0][0]][2]+flank_length)
                    else:                      ref_seq_a=reverse(complementary(ref_seq_readin(ref,bp_let_hash[alt_jun[0][0]][0],bp_let_hash[alt_jun[0][0]][1]-flank_length,bp_let_hash[alt_jun[0][0]][1]+flank_length)))
                    if not '^' in alt_jun[1]:                        ref_seq_b=ref_seq_readin(ref,bp_let_hash[alt_jun[1][0]][0],bp_let_hash[alt_jun[1][0]][1]-flank_length,bp_let_hash[alt_jun[1][0]][1]+flank_length)
                    else:                      ref_seq_b=reverse(complementary(ref_seq_readin(ref,bp_let_hash[alt_jun[1][0]][0],bp_let_hash[alt_jun[1][0]][2]-flank_length,bp_let_hash[alt_jun[1][0]][2]+flank_length)))
                    [window_size,window_size_qc]=window_size_refine(ref_seq_a+ref_seq_b)
                    if not window_size=='Error':
                        alt_seq=ref_seq_a[-flank_length:]+ref_seq_b[:flank_length]
                        [window_size,window_size_qc]=window_size_refine(alt_seq)
                        if not window_size=='Error':
                            if not '^' in alt_jun[0]:                   all_reads_a=simple_del_chop_pacbio_read_simple_short(bam_in,[bp_let_hash[alt_jun[0][0]][0],bp_let_hash[alt_jun[0][0]][2]],flank_length)
                            else:                                       all_reads_a=simple_del_chop_pacbio_read_simple_short(bam_in,[bp_let_hash[alt_jun[0][0]][0],bp_let_hash[alt_jun[0][0]][1]],flank_length)
                            #if not '^' in alt_jun[1]:                   all_reads_b=simple_del_chop_pacbio_read_simple_short(bam_in,[bp_let_hash[alt_jun[0][0]][0],bp_let_hash[alt_jun[0][0]][1]],flank_length)
                            #else:                                       all_reads_b=simple_del_chop_pacbio_read_simple_short(bam_in,[bp_let_hash[alt_jun[0][0]][0],bp_let_hash[alt_jun[0][0]][2]],flank_length)
                            if len(all_reads_a)>0:
                                for x in all_reads_a:
                                    valor_single_read_score=calcu_valor_single_read_score_within_10Perc_m1b(ref_seq_a,alt_seq,x,window_size)
                                    if not valor_single_read_score==[1,1] and not valor_single_read_score[0]==0:
                                        valor_score_list.append(1-float(valor_single_read_score[1])/float(valor_single_read_score[0]))
                            #if len(all_reads_b)>0:
                            #    for x in all_reads_a:
                            #        valor_single_read_score=calcu_valor_single_read_score_within_10Perc_m1b(ref_seq_a,alt_seq,x,window_size)
                            #        if not valor_single_read_score==[1,1] and not valor_single_read_score[0]==0:
                            #            valor_score_list.append(1-float(valor_single_read_score[1])/float(valor_single_read_score[0]))
    return valor_score_list
def long_del_inv(bam_in,ref,sv_info,out_figure_name):
    #eg of sv_info=[['chr19', 46275941, 46314150, 'del'], ['chr19', 46314150, 46314312, 'inv']]
    valor_score_list=[]
    best_read_rec=''
    flank_length=500
    ref_seq=ref_seq_readin(ref,sv_info[0][0],sv_info[0][1]-flank_length,sv_info[1][1]+flank_length) 
    [window_size,window_size_qc]=window_size_refine(ref_seq)
    if not window_size=='Error':
        alt_seq=ref_seq[:flank_length]+reverse(complementary(ref_seq_readin(ref,sv_info[1][0],sv_info[1][2]-flank_length,sv_info[1][2])))
        [window_size,window_size_qc]=window_size_refine(alt_seq)
        if not window_size=='Error':
            all_reads=simple_del_chop_pacbio_read_simple_short(bam_in,sv_info[0],flank_length)
            if len(all_reads)>0:
                for x in all_reads:
                    valor_single_read_score=calcu_valor_single_read_score_within_10Perc_m1b(ref_seq,alt_seq,x,window_size)
                    if not valor_single_read_score==[1,1] and not valor_single_read_score[0]==0:
                        valor_score_list.append(1-float(valor_single_read_score[1])/float(valor_single_read_score[0]))
                        if valor_score_list[-1]==max(valor_score_list):best_read_rec=x
                make_event_figure_1(valor_score_list,best_read_rec,window_size,ref_seq,alt_seq,out_figure_name)
    return valor_score_list
def simple_del_Valor(bam_in,ref,sv_info,out_figure_name):
    #eg of sv_info=['chr1', 101553562, 101553905]
    flank_length=flank_length_calculate(sv_info)
    valor_score_list=[]
    best_read_rec=''
    if sv_info[2]-sv_info[1]<default_max_sv_test: #only try to read in all reads with sv <100K; else: try breakpoints ; 
        all_reads=simple_del_chop_pacbio_read_simple_short(bam_in,sv_info,flank_length)
        if len(all_reads)>0:
            ref_seq=ref_seq_readin(ref,sv_info[0],sv_info[1]-flank_length,sv_info[2]+flank_length)
            [window_size,window_size_qc]=window_size_refine(ref_seq)
            if not window_size=='Error':
                alt_seq=ref_seq[:flank_length]+ref_seq[-flank_length:]
                for x in all_reads:
                    valor_single_read_score=calcu_valor_single_read_score_abs_dis_m1b(ref_seq,alt_seq,x,window_size)
                    if not valor_single_read_score==[1,1] and not valor_single_read_score[0]==0:
                        valor_score_list.append(1-float(valor_single_read_score[1])/float(valor_single_read_score[0]))
                        if valor_score_list[-1]==max(valor_score_list):best_read_rec=x
                make_event_figure_1(valor_score_list,best_read_rec,window_size,ref_seq,alt_seq,out_figure_name)
    else:
        #change the way of read in ref seq
        all_reads=simple_del_chop_pacbio_read_simple_short(bam_in,sv_info,flank_length)
        if len(all_reads)>0:
            ref_seq=ref_seq_readin(ref,sv_info[0],sv_info[1]-flank_length,sv_info[1]+flank_length)
            [window_size,window_size_qc]=window_size_refine(ref_seq)
            if not window_size=='Error':
                alt_seq=ref_seq_readin(ref,sv_info[0],sv_info[1]-flank_length,sv_info[1])+ref_seq_readin(ref,sv_info[0],sv_info[2],sv_info[2]+flank_length)
                [window_size,window_size_qc]=window_size_refine(alt_seq)
                if not window_size=='Error':
                    for x in all_reads:
                        valor_single_read_score=calcu_valor_single_read_score_within_10Perc_m1b(ref_seq,alt_seq,x,window_size)
                        if not valor_single_read_score==[1,1] and not valor_single_read_score[0]==0:
                            valor_score_list.append(1-float(valor_single_read_score[1])/float(valor_single_read_score[0]))
                            if valor_score_list[-1]==max(valor_score_list):best_read_rec=x
                    make_event_figure_1(valor_score_list,best_read_rec,window_size,ref_seq,alt_seq,out_figure_name)
    return valor_score_list
def simple_disdup_Valor(bam_in,ref,sv_info,out_figure_name):
    sv_info[1:3]=[int(i) for i in sv_info[1:3]]
    dup_block=sv_info[:3]
    ins_point=[sv_info[3],int(sv_info[4])]
    flank_length=flank_length_calculate(dup_block)
    valor_score_list=[]
    best_read_rec=''
    bp_info=sorted([int(i) for i in sv_info[1:3]+[sv_info[4]]])
    run_flag=0
    if sv_info[0]==sv_info[3] and max(bp_info)-min(bp_info)<default_max_sv_test: #try to test on the whole region
        ref_seq=ref_seq_readin(ref,sv_info[0],min(bp_info)-flank_length,max(bp_info)+flank_length)
        [window_size,window_size_qc]=window_size_refine(ref_seq)
        if not window_size=='Error':
            all_reads=simple_chop_pacbio_read_simple_short(bam_in,[sv_info[0]]+bp_info+[int(bp_info[-1])+sv_info[2]-sv_info[1]],flank_length)
            if len(all_reads)>0:
                run_flag+=1
                ref_structure=[i for i in 'ab']
                if sv_info[4]>sv_info[2]:   alt_structure=['a','b','a']
                elif sv_info[4]<sv_info[1]: alt_structure=['b','a','b']
                alt_seq=ref_seq_readin(ref,sv_info[0],min(bp_info)-flank_length,min(bp_info))
                a_seq=ref_seq_readin(ref,sv_info[0],bp_info[0],bp_info[1])
                b_seq=ref_seq_readin(ref,sv_info[0],bp_info[1],bp_info[2])
                for x in alt_structure:
                    if x=='a':  alt_seq+=a_seq
                    elif x=='b':    alt_seq+=b_seq
                alt_seq+=ref_seq_readin(ref,sv_info[0],max(bp_info),max(bp_info)+flank_length)
                [window_size,window_size_qc]=window_size_refine(alt_seq)
                if not window_size=='Error':
                    for x in all_reads:
                        valor_single_read_score=calcu_valor_single_read_score_directed_dis_m1b(ref_seq,alt_seq,x,window_size)
                        if not valor_single_read_score==[1,1] and not valor_single_read_score[0]==0:
                            valor_score_list.append(1-float(valor_single_read_score[1])/float(valor_single_read_score[0]))
                            if valor_score_list[-1]==max(valor_score_list):best_read_rec=x
                    make_event_figure_1(valor_score_list,best_read_rec,window_size,ref_seq,alt_seq,out_figure_name)
    if run_flag==0:
        if max(bp_info)-min(bp_info)<default_max_sv_test: #only try to read in all reads with sv <100K; else: try breakpoints ; 
            all_reads=simple_del_chop_pacbio_read_simple_short(bam_in,ins_point,flank_length)
            if len(all_reads)>0:
                ref_seq=ref_seq_readin(ref,ins_point[0],ins_point[1]-flank_length,ins_point[1]+flank_length)
                [window_size,window_size_qc]=window_size_refine(ref_seq)
                if not window_size=='Error':
                    alt_seq=ref_seq[:flank_length]+ref_seq_readin(ref,dup_block[0],dup_block[1],dup_block[2])+ref_seq[-flank_length:]
                    [window_size,window_size_qc]=window_size_refine(alt_seq)
                    if not window_size=='Error':
                        for x in all_reads:
                            valor_single_read_score=calcu_valor_single_read_score_abs_dis_m1b(ref_seq,alt_seq,x,window_size)
                            if not valor_single_read_score==[1,1] and not valor_single_read_score[0]==0:
                                valor_score_list.append(1-float(valor_single_read_score[1])/float(valor_single_read_score[0]))
                                if valor_score_list[-1]==max(valor_score_list):best_read_rec=x
                        make_event_figure_1(valor_score_list,best_read_rec,window_size,ref_seq,alt_seq,out_figure_name)
        else:
            all_reads=simple_del_chop_pacbio_read_simple_short(bam_in,ins_point,flank_length)
            if len(all_reads)>0:
                ref_seq=ref_seq_readin(ref,ins_point[0],ins_point[1]-flank_length,ins_point[1]+flank_length)
                [window_size,window_size_qc]=window_size_refine(ref_seq)
                if not window_size=='Error':
                    alt_seq=ref_seq[:flank_length]+ref_seq_readin(ref,dup_block[0],dup_block[1],dup_block[1]+flank_length)
                    [window_size,window_size_qc]=window_size_refine(alt_seq)
                    if not window_size=='Error':
                        for x in all_reads:
                            valor_single_read_score=calcu_valor_single_read_score_within_10Perc_m1b(ref_seq,alt_seq,x,window_size)
                            if not valor_single_read_score==[1,1] and not valor_single_read_score[0]==0:
                                valor_score_list.append(1-float(valor_single_read_score[1])/float(valor_single_read_score[0]))
                                if valor_score_list[-1]==max(valor_score_list):best_read_rec=x
                        make_event_figure_1(valor_score_list,best_read_rec,window_size,ref_seq,alt_seq,out_figure_name)
    return valor_score_list
def simple_inv_Valor(bam_in,ref,sv_info,out_figure_name):
    #eg of sv_info=['chr1', 101553562, 101553905]
    flank_length=flank_length_calculate(sv_info)
    valor_score_list=[]
    best_read_rec=''
    if sv_info[2]-sv_info[1]<default_max_sv_test: #only try to read in all reads with sv <100K; else: try breakpoints ; 
        ref_seq=ref_seq_readin(ref,sv_info[0],sv_info[1]-flank_length,sv_info[2]+flank_length)
        [window_size,window_size_qc]=window_size_refine(ref_seq)
        if not window_size=='Error':
            alt_seq=ref_seq[:flank_length]+reverse(complementary(ref_seq[flank_length:(-flank_length)]))+ref_seq[-flank_length:]
            [window_size,window_size_qc]=window_size_refine(alt_seq)
            if not window_size=='Error':
                all_reads=simple_chop_pacbio_read_simple_short(bam_in,sv_info,flank_length)
                if len(all_reads)>10:
                    for x in all_reads:
                        valor_single_read_score=calcu_valor_single_read_score_abs_dis_m1b(ref_seq,alt_seq,x,window_size)
                        if not valor_single_read_score==[1,1] and not valor_single_read_score[0]==0:
                            valor_score_list.append(1-float(valor_single_read_score[1])/float(valor_single_read_score[0]))
                            if valor_score_list[-1]==max(valor_score_list):best_read_rec=x
                    make_event_figure_1(valor_score_list,best_read_rec,window_size,ref_seq,alt_seq,out_figure_name)
                    return valor_score_list
    ref_seq=ref_seq_readin(ref,sv_info[0],sv_info[1]-flank_length,sv_info[1]+flank_length)
    [window_size,window_size_qc]=window_size_refine(ref_seq)
    if not window_size=='Error':
        alt_seq=ref_seq[:flank_length]+ref_seq_readin(ref,sv_info[0],sv_info[2]-flank_length,sv_info[2],'TRUE')
        [window_size,window_size_qc]=window_size_refine(alt_seq)
        if not window_size=='Error':
            all_reads=simple_del_chop_pacbio_read_simple_short(bam_in,sv_info,flank_length)
            if len(all_reads)>0:
                for x in all_reads:
                    valor_single_read_score=calcu_valor_single_read_score_within_10Perc_m1b(ref_seq,alt_seq,x,window_size)
                    if not valor_single_read_score==[1,1] and not valor_single_read_score[0]==0:
                        valor_score_list.append(1-float(valor_single_read_score[1])/float(valor_single_read_score[0]))
                        if valor_score_list[-1]==max(valor_score_list):best_read_rec=x
                make_event_figure_1(valor_score_list,best_read_rec,window_size,ref_seq,alt_seq,out_figure_name)
    return valor_score_list
def simple_ins_Valor(bam_in,ref,ins_pos,ins_seq,out_figure_name,POLARITY):
    #eg of ins_pos='chr1_83144055'
    #eg of bam_in='/nfs/turbo/remillsscr/scratch2_trans/datasets/1000genomes/vol1/ftp/data_collections/hgsv_sv_discovery/PacBio/alignment/HG00512.XXX.bam'
    #eg of ref='/scratch/remills_flux/xuefzhao/reference/GRCh38.1KGP/GRCh38_full_analysis_set_plus_decoy_hla.fa'
    if POLARITY=='+':
        ins_seq_2=ins_seq
    elif POLARITY=='-':
        ins_seq_2=reverse(complementary(ins_seq))
    flank_length=default_flank_length if len(ins_seq)>default_flank_length else len(ins_seq)
    valor_score_list=[]
    best_read_rec=''
    all_reads=simple_chop_pacbio_read_simple_short(bam_in,['_'.join(ins_pos.split('_')[:-1]),ins_pos.split('_')[-1]]+[int(ins_pos.split('_')[-1])+len(ins_seq)],flank_length)
    if len(all_reads)>0:
        ref_seq=ref_seq_readin(ref,['_'.join(ins_pos.split('_')[:-1]),ins_pos.split('_')[-1]][0],int(ins_pos.split('_')[-1])-flank_length,int(ins_pos.split('_')[-1])+flank_length+len(ins_seq),reverse_flag='FALSE')
        [window_size,window_size_qc]=window_size_refine(ref_seq)
        if not window_size=='Error':
            alt_seq=ref_seq[:flank_length]+ins_seq_2+ref_seq[flank_length:(2*flank_length)]
            [window_size,window_size_qc]=window_size_refine(alt_seq)
            if not window_size=='Error':
                if len(all_reads)>10:
                    for x in all_reads:
                        if float(x[0].count('N')+x[0].count('n'))/float(len(x[0]))<0.1:
                            valor_single_read_score=calcu_valor_single_read_score_abs_dis_m1b(ref_seq,alt_seq,x,window_size)
                            if not valor_single_read_score==[1,1] and not valor_single_read_score[0]==0:
                                valor_score_list.append(1-float(valor_single_read_score[1])/float(valor_single_read_score[0]))
                                if valor_score_list[-1]==max(valor_score_list):best_read_rec=x
                else:
                    all_reads=simple_chop_pacbio_read_simple_short(bam_in,['_'.join(ins_pos.split('_')[:-1]),ins_pos.split('_')[-1]]+[int(ins_pos.split('_')[-1])],flank_length)
                    for x in all_reads:
                        if float(x[0].count('N')+x[0].count('n'))/float(len(x[0]))<0.1:
                            valor_single_read_score=calcu_valor_single_read_score_within_10Perc_m1b(ref_seq,alt_seq,x,window_size)
                            if not valor_single_read_score==[1,1] and not valor_single_read_score[0]==0:
                                valor_score_list.append(1-float(valor_single_read_score[1])/float(valor_single_read_score[0]))
                                if valor_score_list[-1]==max(valor_score_list):best_read_rec=x
                if ins_seq_2.count('N')==len(ins_seq_2):        make_event_figure_1(valor_score_list,best_read_rec,window_size,ref_seq,ref_seq[2:flank_length],out_figure_name)
                else:            make_event_figure_1(valor_score_list,best_read_rec,window_size,ref_seq,alt_seq,out_figure_name)
    return valor_score_list
def subkeys(key, nth_base, inversions):
    subkeys_info = []
    key=key_modify(key)
    keylen = len(key)
    # speed tip from:
    # http://wiki.python.org/moin/PythonSpeed/PerformanceTips#String_Concatenation
    if nth_base == 1:
        subkeys_info = [key]
    elif nth_base != 0:
        for k in range(nth_base):
            substr_list = [key[j] for j in range(keylen) if (j % nth_base == k)]
            subkeys_info.append("".join(substr_list))
    else:
        # nth_base = 0 is a special case for third base mismatches
        # for every codon, only include the first 2 bases in the hash
        subkeys_info = ["".join([key[i] for i in range(len(key)) if i % 3 != 2])]
    if inversions:
        for i in range(len(subkeys_info)):
            subkeys_info.append("".join([invert_base[c] for c in reversed(subkeys_info[i])]))
    return subkeys_info
def svtype_extract(pin):
    svtype=''
    for x in pin[7].split(';'):
        if 'SVTYPE' in x:
            svtype=x.split('=')[1]
    if svtype=='':
        svtype=pin[4]
    return svtype
def sv_len_extract(pin):
    svtype=''
    for x in pin[7].split(';'):
        if 'SVLEN' in x:
            svtype=x.split('=')[1]
    if svtype=='':
        svtype=0
    return svtype
def sv_insert_point_define(pin):
    svtype=[0,0]
    for x in pin[7].split(';'):
        if 'insert_point=' in x:
            svtype=x.split('=')[1].split(':')
    return svtype
def chr_start_end_extract(pin):
    out=[pin[0],int(pin[1])]
    for x in pin[7].split(';'):
        if 'END=' in x and x.split('=')[0]=='END':
            out.append(int(x.split('=')[1]))
    return out
def ref_seq_readin(ref,chrom,start,end,reverse_flag='FALSE'):
    #reverse=='TRUE': return rev-comp-seq   ; if not specified, default as 'FALSE'
    #else: return original seq
    fref=os.popen(r'''samtools faidx %s %s:%d-%d'''%(ref,chrom,int(start),int(end)))
    fref.readline().strip().split()
    seq=''
    while True:
            pref=fref.readline().strip().split()
            if not pref: break
            seq+=pref[0]
    fref.close()
    if not reverse_flag=='FALSE':
        return seq
    else:
        return reverse(complementary(seq))
def result_organize_ins(info_list):
    #eg of info_list=[key_event,valor_score_event]=['chr2_82961201', [-9.228366096827557, -106.46718851834126, -667.0858781654538, -38.56838396416415, -64.87185751169045, -147.77261544769615, -28.29536680099185, -25.378519434143666, -17.23542013374081, -113.00564782332029, -64.53043553409316]]
    if len(info_list[1])>0:
        pos_values=[i for i in info_list[1] if float(i)>0]
        neg_values=[i for i in info_list[1] if not float(i)>0]
        geno_value=float(len(pos_values))/float(len(pos_values)+len(neg_values))
        if not pos_values==[]:
            qual_value=np.mean(pos_values)
        else:
            qual_value=0
        return [info_list[0]]+[qual_value,geno_value,','.join([str(round(float(i),2)) for i in info_list[1]])]
    else:
        return [info_list[0]]+['NA' for i in range(3)]
def vcf_rec_hash_modify(vcf_rec_hash):
    out={}
    for k1 in vcf_rec_hash.keys():
        if not vcf_rec_hash[k1] in out.keys():  out[vcf_rec_hash[k1]]=[]
        out[vcf_rec_hash[k1]].append(k1)
    return out
def vcf_valor_modify(vcf_input,vcf_rec_hash_new):
    valor_input=vcf_input+'.valor'
    valor_rec={}
    fin=open(vcf_input)
    vcf_info_hash={}
    rec=-1
    for line in fin:
        rec+=1
        pin=line.strip().split()
        if not pin[0][0]=='#':
            if pin[6]=='PASS':
                vcf_info_hash[rec]=pin
    fin.close()
    fin=open(valor_input)
    keep_rec=[]
    for line in fin:
        pin=line.strip().split()
        if pin[0] in vcf_rec_hash_new.keys():   
            sv_rec_label=vcf_rec_hash_new[pin[0]]
            for y in sv_rec_label:
                vcf_info_hash[y]+=[round(float(i),2) if not i=='NA' else i for i in pin[1:3]]
                keep_rec.append(y)
    fin.close()
    fo=open(valor_input,'w')
    print >>fo, '\t'.join(['#CHROM','POS','ID','REF','ALT','QUAL','FILTER','INFO','FORMAT','SAMPLE','VaLoR_value','VaLoR_geno'])
    for k1 in sorted(vcf_info_hash.keys()):
        if k1 in keep_rec:
            print>>fo, '\t'.join([str(i) for i in vcf_info_hash[k1]])
    fo.close()
def window_size_refine(seq2,region_QC_Cff=0.4):
    window_size=10
    if not seq2.count('N')+seq2.count('n')>100:
        dotdata_qual_check=dotdata(window_size,seq2,seq2)
        if len(dotdata_qual_check)>0:
            region_QC=qual_check_repetitive_region(dotdata_qual_check)
            while True:
                if window_size>30: break
                if region_QC[0]>region_QC_Cff or sum(region_QC[1])/float(len(seq2))<0.3: break
                else:
                    window_size+=10
                    dotdata_qual_check=dotdata(window_size,seq2,seq2)
                    region_QC=qual_check_repetitive_region(dotdata_qual_check)
            return [window_size,region_QC]
        else:   return ['Error','Error']
    else:       return ['Error','Error']
def write_output_initiate(out_name):
    fo=open(out_name,'w')
    print >>fo, '\t'.join(['chr','start','end','SV_description','VaLoR_quality_score','VaLoR_genotype_score','other'])
    fo.close()
def write_output_main(out_name,out_list):
    fo=open(out_name,'a')
    print >>fo, '\t'.join([str(i) for i in out_list])
    fo.close()
def ref_seq_readin(ref,chrom,start,end,reverse_flag='FALSE'):
    #reverse=='TRUE': return rev-comp-seq    ; if not specified, default as 'FALSE'
    #else: return original seq
    end-=1
    fref=os.popen(r'''samtools faidx %s %s:%d-%d'''%(ref,chrom,int(start),int(end)))
    fref.readline().strip().split()
    seq=''
    while True:
            pref=fref.readline().strip().split()
            if not pref: break
            seq+=pref[0]
    fref.close()
    if reverse_flag=='FALSE':
        return seq
    else:
        return reverse(complementary(seq))
def unify_list(list):
    out=[]
    for x in list:
        if not x in out:
            out.append(x)
    return out
def x_to_x_modify_new(x,dup_block_combined):
    x_modify=[[i] for i in list(x)]
    block_rec=-1
    for y in dup_block_combined:
        block_rec+=1
        if len(y)>1:
            x_modify[block_rec]+=[x_modify[block_rec][0]+1+i for i in range(len(y)-1)]
    x_modify_new=[]
    for y in x_modify:  x_modify_new+=y
    return x_modify_new
def X_means_cluster(data_list):
    temp_result=[i for i in k_means_cluster(data_list) if not i==[[],[]]]
    if temp_result==[data_list]:
        return temp_result[0]
    else:
        out=[]
        for i in temp_result:
            out+=X_means_cluster(i)
        return out
def X_means_cluster_reformat(data_list):
    out=X_means_cluster(data_list)
    out2=[]
    for y in range(len(out)/2):
        out2.append([out[2*y],out[2*y+1]])
    return out2

def calcu_eu_dis_svelter(global_list,km):
    [lenght_cff,dots_num_cff,clu_dis_cff, point_dis_cff, simple_dis_cff, invert_base, dict_opts, out_path, out_file_Cannot_Validate, sample_name, start, delta, bam_in, ref, chromos, region_QC_Cff, min_length, min_read_compare, case_number, qc_file_name]=global_list
    [k1,k2,k3]=km
    simple_del_test=simple_del_diploid_decide(k1,k2)
    out_hash={}
    if not 'FALSE' in simple_del_test:
        #eg of [k1,k2,k3]=['a/a','/a',['chr1', 909055, 909641]]
        chr_let_hash=bp_to_chr_hash(k3,chromos)
        del_block=[let_to_block_info(i,chr_let_hash,chromos) for i in simple_del_test if not i=='NA']
        if len(unify_list(del_block))==1: #homo- event
            SV_rec=0
            for sv_block in del_block[0]:
                SV_rec+=1
                sv_info=[sv_block[0]]+[int(i) for i in sv_block[1:]]
                valor_score_event=simple_del_Valor(bam_in,ref,sv_info,out_path+sample_name+'.'+':'.join(k3)+'.DEL.png')
                out_hash[SV_rec]=result_organize_ins(valor_score_event)
        else:
            SV_rec=0
            for x in del_block:
                for sv_block in x:
                    SV_rec+=1
                    sv_info=[sv_block[0]]+[int(i) for i in sv_block[1:]]
                    valor_score_event=simple_del_Valor(bam_in,ref,sv_info,out_path+sample_name+'.'+':'.join(k3)+'.DEL.png')
                    out_hash[SV_rec]=result_organize_ins(valor_score_event)
    else:
        simple_inv_test=simple_inv_diploid_decide(k1,k2)
        if not 'FALSE' in simple_inv_test:
            #eg of km=['a/a', 'a^/a', ['chr1', '236755681', '237282111']]
            chr_let_hash=bp_to_chr_hash(k3,chromos)
            inv_block=[let_to_block_info(i,chr_let_hash,chromos) for i in simple_inv_test if not i=='NA']
            SV_rec=0
            for x in unify_list(inv_block):
                for sv_block in x:
                    SV_rec+=1
                    sv_info=[sv_block[0]]+[int(i) for i in sv_block[1:]]
                    valor_score_event=simple_inv_Valor(bam_in,ref,sv_info,out_path+sample_name+':'.join(k3)+'.INV.png')
                    out_hash[SV_rec]=result_organize_ins(['_'.join([k1,k2]+k3),valor_score_event])
        else: 
            simple_disdup_test=simple_disdup_diploid_decide(k1,k2)
            if not 'FALSE' in simple_disdup_test:
                SV_rec=0
                #eg of km=['abc/abc', 'abc/babc',['chr2', '65591889', '65592631', '65593270', '65594027']]
                chr_let_hash=bp_to_chr_hash(k3,chromos)
                dup_block=[ [let_to_block_info(j[1],chr_let_hash,chromos),[chr_let_hash[j[0]][0],chr_let_hash[j[0]][2]] ]    for j in i[1] for i in simple_disdup_test if not i=='NA']
                for x in dup_block:
                    for y in x[0]:
                        SV_rec+=1
                        dup_info=y+x[1]
                        valor_score_event=simple_disdup_Valor(bam_in,ref,dup_info,out_path+sample_name+':'.join([str(i) for i in dup_info])+'.DISDUP.png')
                        out_hash[SV_rec]=result_organize_ins(valor_score_event)
    return out_hash






