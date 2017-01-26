def chop_pacbio_read_simple_short(bam_in,info,flank_length):
	bps=info[2:]
	block_length={}
	for x in range(len(info[2:])-2):
		block_length[chr(97+x)]=int(info[x+4])-int(info[x+3])
	alA_len=np.sum([block_length[x] for x in info[1].split('/')[0] if not x=='^'])
	alB_len=np.sum([block_length[x] for x in info[1].split('/')[1] if not x=='^'])
	alRef_len=int(info[-1])-int(info[3])
	len_cff=max([alA_len,alB_len,alRef_len])
	bam_in_new_list=bam_in_decide(bam_in,bps)
	if bam_in_new_list=='': return [[],[],[]]
	out=[]
	out2=[]
	out3=[]
	tandem_test=[]
	for bam_in_new in bam_in_new_list:
		fbam=os.popen(r'''samtools view %s %s:%d-%d'''%(bam_in_new,bps[0],int(bps[1])-flank_length,int(bps[-1])+flank_length))
		for line in fbam:
			pbam=line.strip().split()
			if not pbam in tandem_test:
				tandem_test.append(pbam)
				#out3.append(int(pbam[3])-int(bps[1])+flank_length)
				if not pbam[0]=='@': 
					if int(pbam[3])<int(bps[1])-flank_length+1:
						align_info=cigar2alignstart(pbam[5],int(pbam[3]),bps,flank_length)
						align_start=align_info[0]
						miss_bp=align_info[1]+1
						if not miss_bp>flank_length/2:
							align_pos=int(pbam[3])
							target_read=pbam[9][align_start:]
							if len(target_read)>flank_length+len_cff:
								out.append(target_read[:max([alA_len,alB_len,alRef_len])+2*flank_length])
								out2.append(miss_bp)
								out3.append(pbam[0])
		fbam.close()
	return [out,out2,out3]
def chop_pacbio_read_simple_insert(bam_in,info,flank_length,insert_seq):
	bps=info[2:]
	block_length={}
	for x in range(len(info[2:])-2):
		block_length[chr(97+x)]=int(info[x+4])-int(info[x+3])
	len_cff=len(insert_seq)+flank_length
	bam_in_new_list=bam_in_decide(bam_in,bps)
	if bam_in_new_list=='': return [[],[],[]]
	out=[]
	out2=[]
	out3=[]
	tandem_test=[]
	for bam_in_new in bam_in_new_list:
		fbam=os.popen(r'''samtools view %s %s:%d-%d'''%(bam_in_new,bps[0],int(bps[1])-flank_length,int(bps[1])+flank_length))
		for line in fbam:
			pbam=line.strip().split()
			if not pbam in tandem_test:
				tandem_test.append(pbam)
				if not pbam[0]=='@': 
					if int(pbam[3])<int(bps[1])-flank_length+1:
						align_info=cigar2alignstart(pbam[5],int(pbam[3]),bps,flank_length)
						print(align_info)
						align_start=align_info[0]
						miss_bp=align_info[1]+1
						if not miss_bp>flank_length/2:
							align_pos=int(pbam[3])
							target_read=pbam[9][align_start:]
							if len(target_read)>flank_length+len_cff:
								out.append(target_read[:len_cff+flank_length])
								out2.append(miss_bp)
								out3.append(pbam[0])
		fbam.close()
	return [out,out2,out3]
def chop_pacbio_read_simple_insert(bam_in,info,flank_length,insert_seq):
	bps=info[2:]
	block_length={}
	for x in range(len(info[2:])-2):
		block_length[chr(97+x)]=int(info[x+4])-int(info[x+3])
	len_cff=len(insert_seq)+flank_length
	bam_in_new_list=bam_in_decide(bam_in,bps)
	if bam_in_new_list=='': return [[],[],[]]
	out=[]
	out2=[]
	out3=[]
	tandem_test=[]
	for bam_in_new in bam_in_new_list:
		fbam=os.popen(r'''samtools view %s %s:%d-%d'''%(bam_in_new,bps[0],int(bps[1])-flank_length,int(bps[1])+flank_length))
		for line in fbam:
			pbam=line.strip().split()
			if not pbam in tandem_test:
				tandem_test.append(pbam)
				if not pbam[0]=='@': 
					if int(pbam[3])<int(bps[1])-flank_length+1:
						align_info=cigar2alignstart(pbam[5],int(pbam[3]),bps,flank_length)
						print(align_info)
						align_start=align_info[0]
						miss_bp=align_info[1]+1
						if not miss_bp>flank_length/2:
							align_pos=int(pbam[3])
							target_read=pbam[9][align_start:]
							if len(target_read)>flank_length+len_cff:
								out.append(target_read[:len_cff+flank_length])
								out2.append(miss_bp)
								out3.append(pbam[0])
		fbam.close()
	return [out,out2,out3]
def chop_pacbio_read_insert_short(bam_in,bps,flank_length,insert_seq):
	bam_in_new_list=bam_in_decide(bam_in,bps)
	if bam_in_new_list=='': return [[],[],[]]
	out=[]
	out2=[]
	out3=[]
	tandem_test=[]
	for bam_in_new in bam_in_new_list:
		fbam=os.popen(r'''samtools view %s %s:%d-%d'''%(bam_in_new,bps[0],int(bps[1])-2*flank_length,int(bps[1])+flank_length))
		for line in fbam:
			pbam=line.strip().split()
			if not pbam in tandem_test:
				tandem_test.append(pbam)
				if not pbam[0]=='@': 
					if int(pbam[3])<int(bps[1])-flank_length+1 and len(pbam[9])>2*flank_length+len(insert_seq):
						align_info=cigar2alignstart_left(pbam[5],int(pbam[3]),bps,flank_length)
						align_start=align_info[0]
						miss_bp=align_info[1]+1
						#print [align_start,miss_bp]
						target_read=pbam[9][align_start:]
						if len(target_read)>2*flank_length+len(insert_seq) and miss_bp<flank_length:
							out.append(target_read[:len(insert_seq)+2*flank_length])
							out2.append(miss_bp)
							out3.append(pbam[0])
		fbam.close()
	return [out,out2,out3]
def chop_pacbio_read_insert_long(bam_in,bps,flank_length,insert_seq):
	bam_in_new_list=bam_in_decide(bam_in,bps)
	if bam_in_new_list=='': return [[],[],[]]
	out=[]
	out2=[]
	out3=[]
	tandem_test=[]
	for bam_in_new in bam_in_new_list:
		fbam=os.popen(r'''samtools view %s %s:%d-%d'''%(bam_in_new,bps[0],int(bps[1])-2*flank_length,int(bps[1])+flank_length))
		for line in fbam:
			pbam=line.strip().split()
			if not pbam in tandem_test:
				tandem_test.append(pbam)
				if not pbam[0]=='@': 
					if int(pbam[3])<int(bps[1])-flank_length+1 and len(pbam[9])>2*flank_length:
						align_info=cigar2alignstart_left(pbam[5],int(pbam[3]),bps,flank_length)
						align_start=align_info[0]
						miss_bp=align_info[1]+1
						#print [align_start,miss_bp]
						target_read=pbam[9][align_start:]
						if len(target_read)>2*flank_length and miss_bp<flank_length:
							out.append(target_read[:2*flank_length])
							out2.append(miss_bp)
							out3.append(pbam[0])
		fbam.close()
	return [out,out2,out3]
def chop_pacbio_read_left(bam_in,bps,flank_length):
	#get reads align starts bps[1]-flank_length, and ends bps[1]+flank
	bam_in_new_list=bam_in_decide(bam_in,bps)
	if bam_in_new_list=='': return [[],[],[]]
	out=[]
	out2=[]
	out3=[]
	tandem_test=[]
	for bam_in_new in bam_in_new_list:
		fbam=os.popen(r'''samtools view %s %s:%d-%d'''%(bam_in_new,bps[0],int(bps[1])-2*flank_length,int(bps[1])+flank_length))
		for line in fbam:
			pbam=line.strip().split()
			if not pbam in tandem_test:
				tandem_test.append(pbam)
				if not pbam[0]=='@': 
					if int(pbam[3])<int(bps[1])-flank_length+1:
						align_info=cigar2alignstart_left(pbam[5],int(pbam[3]),bps,flank_length)
						align_start=align_info[0]
						miss_bp=align_info[1]+1
						#print [align_start,miss_bp]
						target_read=pbam[9][align_start:]
						if len(target_read)>2*flank_length and miss_bp<flank_length:
							out.append(target_read[:2*flank_length])
							out2.append(miss_bp)
							out3.append(pbam[0])
		fbam.close()
	return [out,out2,out3]
def chop_pacbio_read_right(bam_in,bps,flank_length):
	bam_in_new_list=bam_in_decide(bam_in,bps)
	if bam_in_new_list=='': return [[],[],[]]
	out=[]
	out2=[]
	out3=[]
	tandem_test=[]
	for bam_in_new in bam_in_new_list:
		fbam=os.popen(r'''samtools view %s %s:%d-%d'''%(bam_in_new,bps[0],int(bps[2])-flank_length,int(bps[2])+flank_length))
		for line in fbam:
			pbam=line.strip().split()
			if not pbam in tandem_test:
				tandem_test.append(pbam)
				if not pbam[0]=='@': 
					if int(pbam[3])+int(cigar2reaadlength(pbam[5]))>int(bps[2])+flank_length-1:
						align_info=cigar2alignstart_right(flank_length,pbam[5],int(pbam[3]),bps)
						align_end=align_info[0]
						miss_bp=align_info[1]+1
						target_read=pbam[9][:align_end]
						if len(target_read)>2*flank_length:
							out.append(target_read[::-1][:2*flank_length])
							out2.append(miss_bp)
							out3.append(pbam[0])
		fbam.close()
	return [out,out2,out3]
def chop_pacbio_read_left_vcf(bps,flank_length,bam_in):
	#chop reads going from bps[1]-flank_length 
	bps_new=[bps[0],bps[2],bps[2]+flank_length]
	bam_in_new_list=bam_in_decide(bam_in,bps)
	if bam_in_new_list=='': return [[],[],[]]
	out=[]
	out2=[]
	out3=[]
	for bam_in_new in bam_in_new_list:
		fbam=os.popen(r'''samtools view %s %s:%d-%d'''%(bam_in_new,bps_new[0],int(bps_new[1])-flank_length,int(bps_new[1])+flank_length))
		for line in fbam:
			pbam=line.strip().split()
			if not pbam[0]=='@': 
				if int(pbam[3])<int(bps_new[1])-flank_length+1:
					align_info=cigar2alignstart_left(pbam[5],int(pbam[3]),bps_new,flank_length)
					align_start=align_info[0]
					miss_bp=align_info[1]+1
					if align_start<flank_length/2:
					#print [align_start,miss_bp]
						target_read=pbam[9][align_start:]
						if len(target_read)>2*flank_length:
							out.append(target_read[:2*flank_length])
							out2.append(miss_bp)
							out3.append(pbam[0])
		fbam.close()
	return [out,out2,out3]

def calcu_eu_dis_simple_long(global_list,svelter_list,plt_figure_index,PacVal_score_hash,info,SV_index):
	#info=[k1,k2]+k3
	[lenght_cff,dots_num_cff,clu_dis_cff, point_dis_cff, simple_dis_cff, invert_base, dict_opts, out_path, out_file_Cannot_Validate, sample_name, start, delta, bam_in, ref, chromosomes, region_QC_Cff, min_length, min_read_compare, case_number, qc_file_name]=global_list
	[PacVal_file_in,PacVal_file_out]=svelter_list
	window_size=10
	[k1,k2,k3,chrom]=[info[0],info[1],info[2:],info[2]]
	case_name='_'.join([str(i) for i in k3])
	txt_file=out_path+case_name+'.txt'
	flank_length=500
	bp_let_hash=bp_let_to_hash(info,flank_length)
	letters=let_to_letters(k2)
	for x_let in letters:
		for y in range(len(x_let)-1):
			junction=x_let[y]+x_let[y+1]
			if not junc_check(junction,k1)==1: #Normal junctions
				bps_temp=decide_bp_temp(chrom,junction,bp_let_hash)
				all_reads_left=chop_pacbio_read_left(bam_in,bps_temp,flank_length)
				all_reads_right=[[],[],[]]#for now, need a better idea how to deal with this
				read_hash_left=all_reads_left[0]
				miss_hash_left=all_reads_left[1]
				name_hash_left=name_hash_modify(all_reads_left[2])
				read_hash_right=all_reads_right[0]
				miss_hash_right=all_reads_right[1]
				name_hash_right=name_hash_modify(all_reads_right[2])
				seqs=Pacbio_prodce_ref_alt_junc(ref,flank_length,info,junction,bp_let_hash)
				seq2=seqs[0]
				dotdata_qual_check=dotdata(window_size,seq2,seq2)
				region_QC_a=qual_check_repetitive_region(dotdata_qual_check)
				seq3=seqs[1]
				dotdata_qual_check=dotdata(window_size,seq3,seq3)
				region_QC_b=qual_check_repetitive_region(dotdata_qual_check)
				if region_QC_a[0]>region_QC_Cff or sum(region_QC_a[1])/float(len(seq2))<0.3:
					if region_QC_b[0]>region_QC_Cff or sum(region_QC_b[1])/float(len(seq3))<0.3:
						seq4=seqs[2]
						rsquare_ref1=[]
						rsquare_ref2=[]
						rsquare_alt=[]
						seq1=''
						seq_length_limit=min([len(seq2),len(seq3),len(seq4)])
						if not read_hash_left==[]:
							miss_rec=-1
							rec_len=0
							rec_len_rev=1
							new_all_reads_left=read_hash_minimize(all_reads_left)
							read_hash_left=new_all_reads_left[0]
							miss_hash_left=new_all_reads_left[1]
							name_hash_left=name_hash_modify(new_all_reads_left[2])
							for x in read_hash_left:
								miss_rec+=1
								y=x
								y2=miss_hash_left[miss_rec]
								y3=name_hash_left[miss_rec]
								dotdata_ref1=dotdata(window_size,y,seq2[y2:])
								dotdata_ref2=dotdata(window_size,y,seq3[y2:])
								dotdata_alt=dotdata(window_size,y,seq4[y2:])	
								[rsquare_ref1,rsquare_ref2,rsquare_alt]=eu_dis_calcu_simple_long(dotdata_ref1,dotdata_ref2,dotdata_alt,rsquare_ref1,rsquare_ref2,rsquare_alt)	
								if not rsquare_ref2==[]:
									if float(max([rsquare_ref1[-1],rsquare_ref2[-1]]))==0:
										rsquare_alt[-1]+=1
										rsquare_ref1[-1]+=1
										rsquare_ref2[-1]+=1
									if float(rsquare_alt[-1]-max([rsquare_ref1[-1],rsquare_ref2[-1]]))/float(max([rsquare_ref1[-1],rsquare_ref2[-1]]))>rec_len:
										if miss_hash_left[miss_rec]<seq_length_limit:
											rec_len=float(rsquare_alt[-1]-max([rsquare_ref1[-1],rsquare_ref2[-1]]))/float(max([rsquare_ref1[-1],rsquare_ref2[-1]]))
											rec_name=name_hash_left[miss_rec]
											rec_start=miss_hash_left[miss_rec]
											seq1=y
											dotdata_for_record=[dotdata_ref1,dotdata_ref2,dotdata_alt]
									if float(min([rsquare_ref1[-1],rsquare_ref2[-1]]))==0:
										rsquare_alt[-1]+=1
										rsquare_ref1[-1]+=1
										rsquare_ref2[-1]+=1
									if float(rsquare_alt[-1]-min([rsquare_ref1[-1],rsquare_ref2[-1]]))/float(min([rsquare_ref1[-1],rsquare_ref2[-1]]))<rec_len_rev:
										if miss_hash_left[miss_rec]<seq_length_limit:
											rec_len_rev=float(rsquare_alt[-1]-min([rsquare_ref1[-1],rsquare_ref2[-1]]))/float(min([rsquare_ref1[-1],rsquare_ref2[-1]]))
											rec_name_rev=name_hash_left[miss_rec]
											rec_start_rev=miss_hash_left[miss_rec]
											seq1_rev=y
											dotdata_for_record_rev=[dotdata_ref1,dotdata_ref2,dotdata_alt]
							if rec_len>0:	
								if len(dotdata_for_record[0])>0:
									dotplot_subfigure_simple_long(plt_figure_index,window_size,seq1,seq2[rec_start:],seq3[rec_start:],seq4[rec_start:],out_path+sample_name+'.'+case_name+'.alt.'+rec_name+'.png')
							if rec_len_rev<1:
								if len(dotdata_for_record_rev[0])>0:
									dotplot_subfigure_simple_long(plt_figure_index,window_size,seq1_rev,seq2[rec_start_rev:],seq3[rec_start_rev:],seq4[rec_start_rev:],out_path+sample_name+'.'+case_name+'.ref.'+rec_name_rev+'.png')
						if not read_hash_right==[]:
							miss_rec=-1
							rec_len=0
							new_all_reads_right=read_hash_minimize(all_reads_right)
							read_hash_right=new_all_reads_right[0]
							miss_hash_right=new_all_reads_right[1]
							name_hash_right=name_hash_modify(new_all_reads_right[2])
							for x in read_hash_right:
								miss_rec+=1
								y=x
								y2=miss_hash_right[miss_rec]
								y3=name_hash_right[miss_rec]
								dotdata_ref1=dotdata(window_size,y,seq2[::-1])
								dotdata_ref2=dotdata(window_size,y,seq3[::-1])
								dotdata_alt=dotdata(window_size,y,seq4[::-1])							
								[rsquare_ref1,rsquare_ref2,rsquare_alt]=eu_dis_calcu_simple_long(dotdata_ref1,dotdata_ref2,dotdata_alt,rsquare_ref1,rsquare_ref2,rsquare_alt)	
								if not len(rsquare_ref1)==len(rsquare_ref2)==len(rsquare_alt):
									min_len=min([len(rsquare_ref1),len(rsquare_ref2),len(rsquare_alt)])
									rsquare_ref1=rsquare_ref1[:min_len]
									rsquare_ref2=rsquare_ref2[:min_len]
									rsquare_alt=rsquare_alt[:min_len]
								if not rsquare_ref2==[]:
									if max([rsquare_ref2[-1],rsquare_alt[-1]])-rsquare_ref1[-1]>rec_len:
										rec_len=max([rsquare_ref2[-1],rsquare_alt[-1]])-rsquare_ref1[-1]
										rec_name=name_hash_right[miss_rec]
										rec_start=miss_hash_right[miss_rec]
										dotdata_for_record=[dotdata_ref1,dotdata_ref2,dotdata_alt]
							if rec_len>0:
								if len(dotdata_for_record[0])>0:
									dotplot_subfigure_simple_long(plt_figure_index,window_size,seq1,seq2[rec_start:],seq3[rec_start:],seq4[rec_start:],out_path+sample_name+'.'+case_name+rec_name+'.png')
						PacVal_score_hash=write_pacval_individual_stat_simple_long(PacVal_score_hash,txt_file,rsquare_ref1,rsquare_ref2,rsquare_alt,SV_index)
						remove_files_short(txt_file)
				else:
					fo=open(out_file_Cannot_Validate,'a')
					print(' '.join([str(i) for i in info+junction]), file=fo)
					fo.close()
	return PacVal_score_hash
def calcu_fuzzy_dis_complex_short(global_list,svelter_list,plt_figure_index,PacVal_score_hash,info,SV_index):
	[lenght_cff,dots_num_cff,clu_dis_cff, point_dis_cff, simple_dis_cff, invert_base, dict_opts, out_path, out_file_Cannot_Validate, sample_name, start, delta, bam_in, ref, chromosomes, region_QC_Cff, min_length, min_read_compare, case_number, qc_file_name]=global_list
	[PacVal_file_in,PacVal_file_out]=svelter_list
	window_size=10
	[k1,k2,k3]=[info[0],info[1],info[2:]]
	case_name='_'.join([str(i) for i in k3])
	txt_file=out_path+case_name+'.txt'
	[ref_sv,alt_sv,chrom,bps]=[info[0],info[1],info[2],info[2:]]
	flank_length=flank_length_calculate(bps)
	if bps_check(bps,chromosomes)==0:
		bl_len_hash=bl_len_hash_calculate(bps,ref_sv)
		#end_point=end_point_calculate(alt_sv,bl_len_hash)
		all_reads=chop_pacbio_read_simple_short(bam_in,info,flank_length)
		if not all_reads[0]==[]:
			all_reads_new=read_hash_minimize(all_reads)
			[read_hash,miss_hash,name_hash]=[all_reads_new[0],all_reads_new[1],name_hash_modify(all_reads_new[2])]
			[rsquare_ref,rsquare_alt1,rsquare_alt2,rec_len,rec_start,rec_name]=[[],[],[],0,0,'0']
			if not read_hash==[]:
				seqs=Pacbio_prodce_ref_alt_short(ref,flank_length,info)
				seq2=seqs[0]
				[window_size,region_QC]=window_size_refine(region_QC_Cff,window_size,seq2,flank_length)
				if region_QC[0]>region_QC_Cff or sum(region_QC[1])/float(len(seq2))<0.3:
					[seq1,seq3,seq4,miss_rec,dotdata_for_record]=['',seqs[1],seqs[2],-1,[[],[],[]]]
					for x in read_hash:
						miss_rec+=1
						y=x
						y2=miss_hash[miss_rec]
						y3=name_hash[miss_rec]
						dotdata_ref=dotdata(window_size,y,seq2)
						dotdata_alt1=dotdata(window_size,y,seq3)
						dotdata_alt2=dotdata(window_size,y,seq4)							
						dotdata_ref_2=keep_diagnal_and_anti_diag_for_csv(dotdata_ref,clu_dis_cff,dots_num_cff)
						dotdata_alt1_2=keep_diagnal_and_anti_diag_for_csv(dotdata_alt1,clu_dis_cff,dots_num_cff)
						dotdata_alt2_2=keep_diagnal_and_anti_diag_for_csv(dotdata_alt2,clu_dis_cff,dots_num_cff)
						eu_dis_calcu_check_dup(dotdata_ref_2,dotdata_alt1_2,dotdata_alt2_2,rsquare_ref,rsquare_alt1,rsquare_alt2,y2,delta,info)
						if not len(rsquare_ref)==len(rsquare_alt1)==len(rsquare_alt2):
							min_len=min([len(rsquare_ref),len(rsquare_alt1),len(rsquare_alt2)])
							rsquare_ref=rsquare_ref[:min_len]
							rsquare_alt1=rsquare_alt1[:min_len]
							rsquare_alt2=rsquare_alt2[:min_len]
						if not rsquare_alt1==[]:
							if min([rsquare_alt1[-1],rsquare_alt2[-1]])-rsquare_ref[-1]<rec_len:
								rec_len=min([rsquare_alt1[-1],rsquare_alt2[-1]])-rsquare_ref[-1]
								rec_start=y2
								rec_name=y3
								seq1=y
								dotdata_for_record=[dotdata_ref,dotdata_alt1,dotdata_alt2]
							else:
								seq1=y
					if rec_len==0:
						dotdata_for_record=[dotdata_ref,dotdata_alt1,dotdata_alt2]
						rec_name=y3
						rec_start=y2
						seq1=y
					if len(dotdata_for_record[0])>0:
						dotplot_subfigure_simple_short(plt_figure_index,window_size,seq1,seq2[rec_start:],seq3[rec_start:],seq4[rec_start:],out_path+sample_name+'.'+case_name+rec_name+'.png')
						PacVal_score_hash=write_pacval_individual_stat_simple_short(PacVal_score_hash,txt_file,rsquare_ref,rsquare_alt1,rsquare_alt2,SV_index)
				else:
					fo=open(out_file_Cannot_Validate,'a')
					print(' '.join([str(i) for i in info]), file=fo)
					fo.close()
		else:
			write_null_individual_stat(txt_file)
			remove_files_short(txt_file)
	return PacVal_score_hash
def calcu_eu_dis_ins_short(global_list,svelter_list,plt_figure_index,PacVal_score_hash,info,SV_index,insert_seq):
	#info=[k1,k2]+k3
	[lenght_cff,dots_num_cff,clu_dis_cff, point_dis_cff, simple_dis_cff, invert_base, dict_opts, out_path, out_file_Cannot_Validate, sample_name, start, delta, bam_in, ref, chromosomes, region_QC_Cff, min_length, min_read_compare, case_number, qc_file_name]=global_list
	[PacVal_file_in,PacVal_file_out]=svelter_list
	window_size=10
	[k1,k2,k3]=[info[0],info[1],info[2:]]
	case_name='_'.join([str(i) for i in k3])
	txt_file=out_path+case_name+'.txt'
	[ref_sv,alt_sv,chrom,bps]=[info[0],info[1],info[2],info[2:]]
	flank_length=flank_length_calculate(bps)
	insert_seq=insert_seq_construct(insert_seq,info)
	if bps_check(bps,chromosomes)==0:
		bl_len_hash=bl_len_hash_calculate(bps,ref_sv)
		all_reads=chop_pacbio_read_insert_short(bam_in,bps,flank_length,insert_seq)
		if not all_reads[0]==[]:
			all_reads_new=read_hash_minimize(all_reads)
			[read_hash,miss_hash,name_hash]=[all_reads_new[0],all_reads_new[1],name_hash_modify(all_reads_new[2])]
			[rsquare_ref,rsquare_alt1,rsquare_alt2,rec_len,rec_start,rec_name]=[[],[],[],0,0,'0']
			if not read_hash==[]:
				seqs=Pacbio_prodce_ref_alt_ins_short(ref,flank_length,info,insert_seq)
				seq2=seqs[0]
				[window_size,region_QC]=window_size_refine(region_QC_Cff,window_size,seq2,1)
				if region_QC[0]>region_QC_Cff or sum(region_QC[1])/float(len(seq2))<0.3:
					seq3=seqs[1]
					seq4=seqs[2]
					miss_rec=-1
					dotdata_for_record=[[],[],[]]
					seq1=''
					for y in read_hash:
						miss_rec+=1
						y2=miss_hash[miss_rec]
						y3=name_hash[miss_rec]
						dotdata_ref=dotdata(window_size,y,seq2)
						dotdata_alt1=dotdata(window_size,y,seq3)
						dotdata_alt2=dotdata(window_size,y,seq4)
						dotdata_ref_2=keep_diagnal_and_anti_diag_for_csv(dotdata_ref,clu_dis_cff,dots_num_cff)
						dotdata_alt1_2=keep_diagnal_and_anti_diag_for_csv(dotdata_alt1,clu_dis_cff,dots_num_cff)
						dotdata_alt2_2=keep_diagnal_and_anti_diag_for_csv(dotdata_alt2,clu_dis_cff,dots_num_cff)
						[rsquare_ref,rsquare_alt1,rsquare_alt2]=eu_dis_calcu_check_dup(dotdata_ref_2,dotdata_alt1_2,dotdata_alt2_2,rsquare_ref,rsquare_alt1,rsquare_alt2,y2,delta,info)
						if not len(rsquare_ref)==len(rsquare_alt1)==len(rsquare_alt2):
							min_len=min([len(rsquare_ref),len(rsquare_alt1),len(rsquare_alt2)])
							rsquare_ref=rsquare_ref[:min_len]
							rsquare_alt1=rsquare_alt1[:min_len]
							rsquare_alt2=rsquare_alt2[:min_len]
						if not rsquare_alt1==[]:
							if min([rsquare_alt1[-1],rsquare_alt2[-1]])-rsquare_ref[-1]<rec_len:
								rec_len=min([rsquare_alt1[-1],rsquare_alt2[-1]])-rsquare_ref[-1]
								rec_start=y2
								rec_name=y3
								seq1=y
								dotdata_for_record=[dotdata_ref,dotdata_alt1,dotdata_alt2]
						else:
							seq1=y
					if rec_len==0:
						dotdata_for_record=[dotdata_ref,dotdata_alt1,dotdata_alt2]
						rec_name=y3
						rec_start=y2
						seq1=y
					if len(dotdata_for_record[0])>0:
						dotplot_subfigure_simple_short(plt_figure_index,window_size,seq1,seq2[rec_start:],seq3[rec_start:],seq4[rec_start:],out_path+sample_name+'.'+case_name+rec_name+'.png')
						PacVal_score_hash=write_pacval_individual_stat_simple_short(PacVal_score_hash,txt_file,rsquare_ref,rsquare_alt1,rsquare_alt2,SV_index)
				else:
					fo=open(out_file_Cannot_Validate,'a')
					print(' '.join([str(i) for i in info]), file=fo)
					fo.close()
		else:
			write_null_individual_stat(txt_file)
			remove_files_short(txt_file)
	return PacVal_score_hash
def calcu_eu_dis_ins_long(global_list,svelter_list,plt_figure_index,PacVal_score_hash,info,SV_index,insert_seq):
	#info=[k1,k2]+k3
	[lenght_cff,dots_num_cff,clu_dis_cff, point_dis_cff, simple_dis_cff, invert_base, dict_opts, out_path, out_file_Cannot_Validate, sample_name, start, delta, bam_in, ref, chromosomes, region_QC_Cff, min_length, min_read_compare, case_number, qc_file_name]=global_list
	[PacVal_file_in,PacVal_file_out]=svelter_list
	window_size=10
	[k1,k2,k3]=[info[0],info[1],info[2:]]
	case_name='_'.join([str(i) for i in k3])
	txt_file=out_path+case_name+'.txt'
	[ref_sv,alt_sv,chrom,bps]=[info[0],info[1],info[2],info[2:]]
	flank_length=flank_length_calculate(bps)
	if insert_seq in ['homo','het','']:
		insert_seq=''.join(['n' for i in range(info[4]-info[3])])
	if bps_check(bps,chromosomes)==0:
		bl_len_hash=bl_len_hash_calculate(bps,ref_sv)
		all_reads=chop_pacbio_read_insert_long(bam_in,bps,flank_length,insert_seq)
		if not all_reads[0]==[]:
			all_reads_new=read_hash_minimize(all_reads)
			[read_hash,miss_hash,name_hash]=[all_reads_new[0],all_reads_new[1],name_hash_modify(all_reads_new[2])]
			[rsquare_ref1,rsquare_ref2,rsquare_alt,rec_len,rec_start,rec_name]=[[],[],[],0,0,'0']
			if not read_hash==[]:
				seqs=Pacbio_prodce_ref_alt_ins_short(ref,flank_length,info,insert_seq)
				seq2=seqs[0]
				[window_size,region_QC]=window_size_refine(region_QC_Cff,window_size,seq2,1)
				if region_QC[0]>region_QC_Cff or sum(region_QC[1])/float(len(seq2))<0.3:
					seq3=seqs[1]
					seq4=seqs[2]
					miss_rec=-1
					dotdata_for_record=[[],[],[]]
					seq1=''
					for y in read_hash:
						miss_rec+=1
						y2=miss_hash[miss_rec]
						y3=name_hash[miss_rec]
						[dotdata_ref,dotdata_alt1,dotdata_alt2]=[dotdata(window_size,y,seq2),dotdata(window_size,y,seq3),dotdata(window_size,y,seq4)]
						[dotdata_ref_2,dotdata_alt1_2,dotdata_alt2_2]=[keep_diagnal_and_anti_diag_for_csv(dotdata_ref,clu_dis_cff,dots_num_cff),keep_diagnal_and_anti_diag_for_csv(dotdata_alt1,clu_dis_cff,dots_num_cff),keep_diagnal_and_anti_diag_for_csv(dotdata_alt2,clu_dis_cff,dots_num_cff)]
						[rsquare_ref,rsquare_alt1,rsquare_alt2]=eu_dis_calcu_check_dup(dotdata_ref_2,dotdata_alt1_2,dotdata_alt2_2,rsquare_ref,rsquare_alt1,rsquare_alt2,y2,delta,info)
						if not len(rsquare_ref)==len(rsquare_alt1)==len(rsquare_alt2):
							min_len=min([len(rsquare_ref),len(rsquare_alt1),len(rsquare_alt2)])
							rsquare_ref=rsquare_ref[:min_len]
							rsquare_alt1=rsquare_alt1[:min_len]
							rsquare_alt2=rsquare_alt2[:min_len]
						if not rsquare_alt1==[]:
							if min([rsquare_alt1[-1],rsquare_alt2[-1]])-rsquare_ref[-1]<rec_len:
								rec_len=min([rsquare_alt1[-1],rsquare_alt2[-1]])-rsquare_ref[-1]
								rec_start=y2
								rec_name=y3
								seq1=y
								dotdata_for_record=[dotdata_ref,dotdata_alt1,dotdata_alt2]
						else:
							seq1=y
					if rec_len==0:
						dotdata_for_record=[dotdata_ref,dotdata_alt1,dotdata_alt2]
						rec_name=y3
						rec_start=y2
						seq1=y
					if len(dotdata_for_record[0])>0:
						dotplot_subfigure_simple_short(plt_figure_index,window_size,seq1,seq2[rec_start:],seq3[rec_start:],seq4[rec_start:],out_path+sample_name+'.'+case_name+rec_name+'.png')
						PacVal_score_hash=write_pacval_individual_stat_simple_short(PacVal_score_hash,txt_file,rsquare_ref,rsquare_alt1,rsquare_alt2,SV_index)
				else:
					fo=open(out_file_Cannot_Validate,'a')
					print(' '.join([str(i) for i in info]), file=fo)
					fo.close()
		else:
			write_null_individual_stat(txt_file)
			remove_files_short(txt_file)
	return PacVal_score_hash
def calcu_eu_dis_csv_long_vcf(global_list,svelter_list,plt_figure_index,PacVal_score_hash,info,SV_index):
	#eg of info:['22', '50934615', '22', 50567289, 50567789, '', '22', 50934615, 50935115, '.', '.']
	[lenght_cff,dots_num_cff,clu_dis_cff, point_dis_cff, simple_dis_cff, invert_base, dict_opts, out_path, out_file_Cannot_Validate, sample_name, start, delta, bam_in, ref, chromosomes, region_QC_Cff, min_length, min_read_compare, case_number, qc_file_name]=global_list
	[PacVal_file_in,PacVal_file_out]=svelter_list
	flank_length=500
	window_size=10
	case_name=str(SV_index)
	txt_file=out_path+case_name+'.txt'
	all_reads_left=chop_pacbio_read_left_vcf(info[2:5],flank_length,bam_in)
	read_hash_left=all_reads_left[0]
	miss_hash_left=all_reads_left[1]
	name_hash_left=name_hash_modify(all_reads_left[2])
	read_hash_right=[]
	seqs=Pacbio_prodce_ref_alt_multi_chrom(out_path,ref,flank_length,info)
	seq2=seqs[0]
	seq2_QC=seqs[3]
	[window_size_a,region_QC_a]=window_size_refine(region_QC_Cff,window_size,seq2,flank_length)
	#dotdata_qual_check=dotdata(window_size,seq2_QC,seq2_QC)
	#region_QC_a=qual_check_repetitive_region(dotdata_qual_check)
	seq3=seqs[1]
	seq3_QC=seqs[4]
	[window_size_b,region_QC_b]=window_size_refine(region_QC_Cff,window_size,seq2,flank_length)
	#dotdata_qual_check=dotdata(window_size,seq3_QC,seq3_QC)
	#region_QC_b=qual_check_repetitive_region(dotdata_qual_check)
	window_size=max([window_size_a,window_size_b])
	if region_QC_a[0]>region_QC_Cff or sum(region_QC_a[1])/float(len(seq2))<0.3 or region_QC_b[0]>region_QC_Cff or sum(region_QC_b[1])/float(len(seq3))<0.3:
			seq4=seqs[2]
			rsquare_ref1=[]
			rsquare_ref2=[]
			rsquare_alt=[]
			if not read_hash_left==[]:
				new_read_hash_left=read_hash_minimize(all_reads_left)
				read_hash_left=new_read_hash_left[0]
				miss_hash_left=new_read_hash_left[1]
				name_hash_left=new_read_hash_left[2]
				miss_rec=-1
				rec_len=0
				for x in read_hash_left:
					miss_rec+=1
					y=x
					y2=miss_hash_left[miss_rec]
					y3=name_hash_left[miss_rec]
					dotdata_ref1=dotdata(window_size,y,seq2[y2:])
					dotdata_ref2=dotdata(window_size,y,seq3[y2:])
					dotdata_alt=dotdata(window_size,y,seq4[y2:])							
					[rsquare_ref1,rsquare_ref2,rsquare_alt]=eu_dis_calcu_simple_long(dotdata_ref1,dotdata_ref2,dotdata_alt,rsquare_ref1,rsquare_ref2,rsquare_alt)	
					if not rsquare_alt==[]:
						if rsquare_alt[-1]-max([rsquare_ref1[-1],rsquare_ref2[-1]])>rec_len:
							rec_len=rsquare_alt[-1]-max([rsquare_ref1[-1],rsquare_ref2[-1]])
							rec_name=name_hash_left[miss_rec]
							rec_start=miss_hash_left[miss_rec]
							dotdata_for_record=[dotdata_ref1,dotdata_ref2,dotdata_alt]
							seq1=y
				if rec_len>0:
					if len(dotdata_for_record[0])>0:
						dotplot_subfigure_simple_long(plt_figure_index,window_size,seq1,seq2[rec_start:],seq3[rec_start:],seq4[rec_start:],out_path+sample_name+'.'+case_name+rec_name+'.png')
			write_pacval_individual_stat_simple_long(txt_file,rsquare_ref1,rsquare_ref2,rsquare_alt,SV_index)
			remove_files_short(txt_file)
	else:
		fo=open(out_file_Cannot_Validate,'a')
		print(info)
		print(' '.join([str(i) for i in info]), file=fo)
		fo.close()

def Pacbio_prodce_ref_alt_short(ref,flank_length,info):
	#fin=open(txt_file)
	#pin=fin.readline().strip().split()
	#fin.close()
	ref_sv=info[0]
	ref_sv='/'.join(['1'+ref_sv.split('/')[0]+'2','1'+ref_sv.split('/')[0]+'2'])
	alt_sv=info[1]
	alt_sv='/'.join(['1'+alt_sv.split('/')[0]+'2','1'+alt_sv.split('/')[1]+'2'])
	chrom=info[2]
	bps=info[2:]
	bps=[bps[0]]+[str(int(bps[1])-flank_length)]+bps[1:]+[str(int(bps[-1])+flank_length)]
	ref_hash={}
	rec=0
	for x in ref_sv.split('/')[0]:
		rec+=1
		fref=os.popen(r'''samtools faidx %s %s:%d-%d'''%(ref,chrom,int(bps[rec]), int(bps[rec+1])))
		fref.readline().strip().split()
		seq=''
		while True:
				pref=fref.readline().strip().split()
				if not pref: break
				seq+=pref[0]
		fref.close()
		ref_hash[x]=seq
	ref_hash=ref_hash_modi(ref_hash)
	alt_sv_list=alt_sv_to_list(alt_sv)
	ref_seq=''.join([ref_hash[x] for x in ref_sv.split('/')[0]])
	alt_1_seq=''.join([ref_hash[x] for x in alt_sv_list[0]])
	alt_2_seq=''.join([ref_hash[x] for x in alt_sv_list[1]])
	return [upper_string(ref_seq),upper_string(alt_1_seq),upper_string(alt_2_seq)]
def Pacbio_prodce_ref_alt_ins_short(ref,flank_length,info,insert_seq):
	ref_sv='12/12'
	alt_sv='1'+info[1].split('/')[0]+'2/1'+info[1].split('/')[1]+'2'
	chrom=info[2]
	bps=info[2:]
	bps=[bps[0],str(int(bps[1])-flank_length),bps[1],str(int(bps[1])+flank_length)]
	ref_hash={}
	rec=0
	for x in ref_sv.split('/')[0]:
		rec+=1
		fref=os.popen(r'''samtools faidx %s %s:%d-%d'''%(ref,chrom,int(bps[rec]), int(bps[rec+1])))
		fref.readline().strip().split()
		seq=''
		while True:
				pref=fref.readline().strip().split()
				if not pref: break
				seq+=pref[0]
		fref.close()
		ref_hash[x]=seq
	ref_hash['a']=insert_seq
	ref_hash=ref_hash_modi(ref_hash)
	alt_sv_list=alt_sv_to_list(alt_sv)
	ref_seq=''.join([ref_hash[x] for x in ref_sv.split('/')[0]])
	alt_1_seq=''.join([ref_hash[x] for x in alt_sv_list[0]])
	alt_2_seq=''.join([ref_hash[x] for x in alt_sv_list[1]])
	return [upper_string(ref_seq),upper_string(alt_1_seq),upper_string(alt_2_seq)]
def Pacbio_prodce_ref_alt_ins_long(ref,flank_length,info,insert_seq):
	ref_sv='12/12'
	alt_sv='1'+info[1].split('/')[0]+'2/1'+info[1].split('/')[1]+'2'
	chrom=info[2]
	bps=info[2:]
	bps=[bps[0],str(int(bps[1])-flank_length),bps[1],str(int(bps[1])+flank_length)]
	ref_hash={}
	rec=0
	for x in ref_sv.split('/')[0]:
		rec+=1
		fref=os.popen(r'''samtools faidx %s %s:%d-%d'''%(ref,chrom,int(bps[rec]), int(bps[rec+1])))
		fref.readline().strip().split()
		seq=''
		while True:
				pref=fref.readline().strip().split()
				if not pref: break
				seq+=pref[0]
		fref.close()
		ref_hash[x]=seq
	ref_hash['a']=insert_seq
	ref_hash=ref_hash_modi(ref_hash)
	alt_sv_list=alt_sv_to_list(alt_sv)
	ref_seq=ref_hash['1']+ref_hash['2']
	alt_1_seq=ref_hash['1']+insert_seq[:flank_length]
	alt_2_seq=insert_seq[-flank_length:]+ref_hash['2']
	return [upper_string(ref_seq),upper_string(alt_1_seq),upper_string(alt_2_seq)]
def Pacbio_prodce_ref_alt_junc(ref,flank_length,info,junction,bp_let_hash):
	chrom=info[2]
	ref_hash={}
	if not '^' in junction[0] and not '^' in junction[1]:
		#['a','c']
		bps=[chrom,bp_let_hash[junction[0]][1]-flank_length,bp_let_hash[junction[0]][1],bp_let_hash[junction[1]][0],bp_let_hash[junction[1]][0]+flank_length]
		bps_left=[chrom,bp_let_hash[junction[0]][1]]
		bps_right=[chrom,bp_let_hash[junction[1]][0]]
		ref_hash[junction[0]+'_ri']=ref_seq_readin(ref,chrom,int(bps[1]),int(bps[2]),'FALSE')
		ref_hash[junction[1]+'_le']=ref_seq_readin(ref,chrom,int(bps[3]),int(bps[4]),'FALSE')
	elif '^' in junction[0] and not '^' in junction[1]:
		#['a^','c']
		bps=[chrom,bp_let_hash[junction[0]][0],bp_let_hash[junction[0]][0]+flank_length,bp_let_hash[junction[1]][0],bp_let_hash[junction[1]][0]+flank_length]
		bps_left=[chrom,bp_let_hash[junction[0]][0]]
		bps_right=[chrom,bp_let_hash[junction[1]][0]]
		ref_hash[junction[0]+'_ri']=ref_seq_readin(ref,chrom,int(bps[1]),int(bps[2]),'TRUE')
		ref_hash[junction[1]+'_le']=ref_seq_readin(ref,chrom,int(bps[3]),int(bps[4]),'FALSE')
	elif not '^' in junction[0] and '^' in junction[1]:
		#['a','c^']
		bps=[chrom,bp_let_hash[junction[0]][1]-flank_length,bp_let_hash[junction[0]][1],bp_let_hash[junction[1]][1]-flank_length,bp_let_hash[junction[1]][1]]
		bps_left=[chrom,bp_let_hash[junction[0]][1]]
		bps_right=[chrom,bp_let_hash[junction[1]][1]]
		ref_hash[junction[0]+'_ri']=ref_seq_readin(ref,chrom,int(bps[1]),int(bps[2]),'FALSE')
		ref_hash[junction[1]+'_le']=ref_seq_readin(ref,chrom,int(bps[3]),int(bps[4]),'TRUE')
	elif '^' in junction[0] and '^' in junction[1]:
		#['a^','c^']
		bps=[chrom,bp_let_hash[junction[0]][0],bp_let_hash[junction[0]][0]+flank_length,bp_let_hash[junction[1]][1]-flank_length,bp_let_hash[junction[1]][1]]			
		bps_left=[chrom,bp_let_hash[junction[0]][0]]
		bps_right=[chrom,bp_let_hash[junction[1]][1]]
		ref_hash[junction[0]+'_ri']=ref_seq_readin(ref,chrom,int(bps[1]),int(bps[2]),'TRUE')
		ref_hash[junction[1]+'_le']=ref_seq_readin(ref,chrom,int(bps[3]),int(bps[4]),'TRUE')
	ref_hash['ref_le']=ref_seq_readin(ref,chrom,bps_left[1]-flank_length,bps_left[1]+flank_length,'FALSE')
	ref_hash['ref_ri']=ref_seq_readin(ref,chrom,bps_right[1]-flank_length,bps_right[1]+flank_length,'FALSE')
	return [upper_string(ref_hash['ref_le']),upper_string(ref_hash['ref_ri']),upper_string(ref_hash[junction[0]+'_ri']+ref_hash[junction[1]+'_le'])]
	#return [ref1,ref2,alt]
def Pacbio_prodce_ref_alt_multi_chrom(out_path,ref,flank_length,k1):
	#eg of k1:['22', '50934615', 'chr22', 50567289, 50567789, '', '22', 50934615, 50935115, '.', '.']
	chrom=k1[0]
	ref_hash={}
	alt_a=k1[2:5]
	alt_b=k1[6:9]
	reverse_flag=[]
	for i in k1[9:]:
		if i=='.':
			reverse_flag.append('FALSE')
		elif i=='^':
			reverse_flag.append('TRUE')
		else:
			print('Error: wrong k1 info, reverse not properly stated.')
	if reverse_flag==['FALSE','FALSE']: 
		ref_seq_a=ref_seq_readin(ref,chrom,int(alt_a[2])-flank_length,int(alt_a[2])+flank_length,'FALSE')
		ref_seq_b=ref_seq_readin(ref,chrom,int(alt_b[1])-flank_length,int(alt_b[1])+flank_length,'FALSE')
		ref_seq_a2=ref_seq_readin(ref,chrom,int(alt_a[2])-flank_length,int(alt_a[2]),'FALSE')
		ref_seq_b2=ref_seq_readin(ref,chrom,int(alt_b[1]),int(alt_b[1])+flank_length,'FALSE')
	elif reverse_flag==['FALSE','TRUE']:
		ref_seq_a=ref_seq_readin(ref,chrom,int(alt_a[2])-flank_length,int(alt_a[2])+flank_length,'FALSE')
		ref_seq_b=ref_seq_readin(ref,chrom,int(alt_b[2])-flank_length,int(alt_b[2])+flank_length,'FALSE')
		ref_seq_a2=ref_seq_readin(ref,chrom,int(alt_a[2])-flank_length,int(alt_a[2]),'FALSE')
		ref_seq_b2=ref_seq_readin(ref,chrom,int(alt_b[2])-flank_length,int(alt_b[2]),'FALSE')
	elif reverse_flag==['TRUE','FALSE']:
		ref_seq_a=ref_seq_readin(ref,chrom,int(alt_a[1])-flank_length,int(alt_a[1])+flank_length,'FALSE')
		ref_seq_b=ref_seq_readin(ref,chrom,int(alt_b[1])-flank_length,int(alt_b[1])+flank_length,'FALSE')
		ref_seq_a2=ref_seq_readin(ref,chrom,int(alt_a[1]),int(alt_a[1])+flank_length,'FALSE')
		ref_seq_b2=ref_seq_readin(ref,chrom,int(alt_b[1]),int(alt_b[1])+flank_length,'FALSE')
	alt_seq=ref_seq_readin(ref,chrom,int(alt_a[1]),int(alt_a[2]),reverse_flag[0])+ref_seq_readin(ref,chrom,int(alt_b[1]),int(alt_b[2]),reverse_flag[1])
	return [ref_seq_a,ref_seq_b,alt_seq,ref_seq_a2,ref_seq_b2]

def eu_dis_calcu_simple_long_single(fo_ref):
	#return number of dots falling within 10# dis to diagnal
	temp_data1=[[],[]]
	for x in fo_ref:
		temp_data1[0].append(int(x[0]))
		temp_data1[1].append(int(x[1]))
	temp_data2=[[],[]]
	for x in range(len(temp_data1[0])):
		if np.abs(temp_data1[0][x]-temp_data1[1][x])<float(temp_data1[0][x])/10.0:
			temp_data2[0].append(temp_data1[0][x])
			temp_data2[1].append(temp_data1[1][x])			
	return len(temp_data2[0])
def eu_dis_calcu_simple_long(fo1_ref,fo2_ref,fo_alt,rsquare_ref1,rsquare_ref2,rsquare_alt):
	fo1_ref_score=eu_dis_calcu_simple_long_single(fo1_ref)
	fo2_ref_score=eu_dis_calcu_simple_long_single(fo2_ref)
	fo_alt_score=eu_dis_calcu_simple_long_single(fo_alt)
	if not fo1_ref_score==0 and not fo2_ref_score==0 and not fo_alt_score==0:
		rsquare_ref1.append(fo1_ref_score)
		rsquare_ref2.append(fo2_ref_score)
		rsquare_alt.append(fo_alt_score)
	return [rsquare_ref1,rsquare_ref2,rsquare_alt]
def eu_dis_calcu_add_opposite(fo_ref,fo1_alt,fo2_alt,rsquare_ref,rsquare_alt1,rsquare_alt2,y2,delta,info):
	#calculate total distance of all dots to diagnal; add symmetric dot simultaneously
	#Calculate total distance of all dots to diagnal; smaller=better.
	structure1=info[0].split('/')[0]
	structure2=info[1].split('/')[0]
	structure3=info[1].split('/')[1]
	if sum([int(dup_decide(structure1)),int(dup_decide(structure2)),int(dup_decide(structure3))])>0:
			rsquare_ref=eu_dis_calcu_1(fo_ref,rsquare_ref,y2,delta)
			rsquare_alt1=eu_dis_calcu_1(fo1_alt,rsquare_alt1,y2,delta)
			rsquare_alt2=eu_dis_calcu_1(fo2_alt,rsquare_alt2,y2,delta)
	else:
			rsquare_ref=eu_dis_calcu_2(fo_ref,rsquare_ref,y2,delta)
			rsquare_alt1=eu_dis_calcu_2(fo1_alt,rsquare_alt1,y2,delta)
			rsquare_alt2=eu_dis_calcu_2(fo2_alt,rsquare_alt2,y2,delta)
	return [rsquare_ref,rsquare_alt1,rsquare_alt2]	
def eu_dis_calcu_check_dup(fo_ref,fo1_alt,fo2_alt,rsquare_ref,rsquare_alt1,rsquare_alt2,y2,delta,info):
	#Calculate total distance of all dots to diagnal; smaller=better.
	structure1=info[0].split('/')[0]
	structure2=info[1].split('/')[0]
	structure3=info[1].split('/')[1]
	if sum([int(dup_decide(structure1)),int(dup_decide(structure2)),int(dup_decide(structure3))])>0:	#calculate directioned dis when dup exits in alt structure
		rsquare_ref=eu_dis_calcu_1(fo_ref,rsquare_ref,y2,delta)
		rsquare_alt1=eu_dis_calcu_1(fo1_alt,rsquare_alt1,y2,delta)
		rsquare_alt2=eu_dis_calcu_1(fo2_alt,rsquare_alt2,y2,delta)
	if '^' in structure2 or '^' in structure3:	#calculate directioned dis when inversion do not show
		rsquare_ref=eu_dis_calcu_1(fo_ref,rsquare_ref,y2,delta)
		rsquare_alt1=eu_dis_calcu_1(fo1_alt,rsquare_alt1,y2,delta)
		rsquare_alt2=eu_dis_calcu_1(fo2_alt,rsquare_alt2,y2,delta)		
	else:
		rsquare_ref=eu_dis_calcu_2(fo_ref,rsquare_ref,y2,delta)
		rsquare_alt1=eu_dis_calcu_2(fo1_alt,rsquare_alt1,y2,delta)
		rsquare_alt2=eu_dis_calcu_2(fo2_alt,rsquare_alt2,y2,delta)
	return [rsquare_ref,rsquare_alt1,rsquare_alt2]

def dotplot_subfigure_simple_short(plt_figure_index,kmerlen,seq1,seq2,seq3,seq4,figurename):
	nth_base = 1
	inversions = True
	hits_ref_ref = kmerhits(seq2, seq2, kmerlen, nth_base, inversions)
	hits_ref = kmerhits(seq1, seq2, kmerlen, nth_base, inversions)
	hits_alt1 = kmerhits(seq1, seq3, kmerlen, nth_base, inversions)
	hits_alt2 = kmerhits(seq1, seq4, kmerlen, nth_base, inversions)
	if float(len(hits_ref))/float(len(hits_ref_ref))<0.1:
		kmerlen_new=int(0.8*kmerlen)
		hits_ref = kmerhits(seq1, seq2, kmerlen, nth_base, inversions)
		hits_alt1 = kmerhits(seq1, seq3, kmerlen, nth_base, inversions)
		hits_alt2 = kmerhits(seq1, seq4, kmerlen, nth_base, inversions)
	elif float(len(hits_ref))/float(len(hits_ref_ref))>10:
		kmerlen_new=int(1.2*kmerlen)
		hits_ref = kmerhits(seq1, seq2, kmerlen, nth_base, inversions)
		hits_alt1 = kmerhits(seq1, seq3, kmerlen, nth_base, inversions)
		hits_alt2 = kmerhits(seq1, seq4, kmerlen, nth_base, inversions)
	fig=plt.figure(plt_figure_index)
	makeDotplot_subfigure(hits_ref_ref,'ref vs. ref',221)
	makeDotplot_subfigure(hits_ref,'read vs. ref',222)
	makeDotplot_subfigure(hits_alt1,'read vs. allele1',223)
	makeDotplot_subfigure(hits_alt2,'read vs. allele2',224)
	plt.savefig(figurename)
	#plt.show()
	plt.close(fig)
def dotplot_subfigure_simple_long(plt_figure_index,kmerlen,seq1,seq2,seq3,seq4,figurename):
	nth_base = 1
	inversions = True
	hits_ref_ref = kmerhits(seq2, seq2, kmerlen, nth_base, inversions)
	hits_ref1 = kmerhits(seq1, seq2, kmerlen, nth_base, inversions)
	hits_ref2 = kmerhits(seq1, seq3, kmerlen, nth_base, inversions)
	hits_alt = kmerhits(seq1, seq4, kmerlen, nth_base, inversions)
	if float(len(hits_ref1))/float(len(hits_ref_ref))<0.1:
		kmerlen_new=int(0.8*kmerlen)
		hits_ref1 = kmerhits(seq1, seq2, kmerlen, nth_base, inversions)
		hits_ref2 = kmerhits(seq1, seq3, kmerlen, nth_base, inversions)
		hits_alt = kmerhits(seq1, seq4, kmerlen, nth_base, inversions)
	elif float(len(hits_ref1))/float(len(hits_ref_ref))>10:
		kmerlen_new=int(1.2*kmerlen)
		hits_ref1 = kmerhits(seq1, seq2, kmerlen, nth_base, inversions)
		hits_ref2 = kmerhits(seq1, seq3, kmerlen, nth_base, inversions)
		hits_alt = kmerhits(seq1, seq4, kmerlen, nth_base, inversions)
	fig=plt.figure(plt_figure_index)
	makeDotplot_subfigure(hits_ref_ref,'ref vs. ref',221)
	makeDotplot_subfigure(hits_alt,'read vs. alt_junction',222)
	makeDotplot_subfigure(hits_ref1,'read vs. ref_left',223)
	makeDotplot_subfigure(hits_ref2,'read vs. ref_right',224)
	plt.savefig(figurename)
	#plt.show()
	plt.close(fig)
def write_list(list,fileout):
	fin=open(fileout,'w')
	for line in list:
		print(' '.join([str(i) for i in line]), file=fin)
	fin.close()
def qc_ref_region_check(qc_file_name,window_size,ref,out_path,sample_name,info):
	k1=info[0]
	k2=info[1]
	k3=info[2:]
	case_name='_'.join([str(i) for i in k3+[k1.replace('/','-'),k2.replace('/','-')]])
	txt_file=out_path+case_name+'.txt'
	#fo=open(txt_file,'w')
	#print >>fo, ' '.join([str(i) for i in [k1,k2]+k3])
	#fo.close()
	ref_sv=info[0]
	alt_sv=info[1]
	chrom=info[2]
	bps=info[2:]
	#delta=delta_calculate(bps)
	if bps[2]-bps[1]<5000:
		flank_length=flank_length_calculate(bps)
		seqs=Pacbio_prodce_ref_alt_short(ref,flank_length,info)
		seq2=seqs[0]
		dotdata_qual_check=dotdata(window_size,seq2[flank_length:-flank_length],seq2[flank_length:-flank_length])
		region_QC=qual_check_repetitive_region(dotdata_qual_check)
		QC_plot_name=out_path+sample_name+'.'+'.'.join([str(i) for i in bps])+'.png'
		p = makeDotplot(QC_plot_name, dotdata_qual_check, len(seq2), len(seq2))
		qc_ref_file_write(qc_file_name,sample_name,bps,region_QC,'ab',seq2)		
	else:
		flank_length=500
		info_a=info[:3]+[info[3]-500,info[3]+500]
		info_b=info[:3]+[info[4]-500,info[4]+500]
		#left bp flank
		seqs=Pacbio_prodce_ref_alt_short(ref,flank_length,info_a)
		seq2=seqs[0]
		dotdata_qual_check=dotdata(window_size,seq2,seq2)
		region_QC=qual_check_repetitive_region(dotdata_qual_check)
		QC_plot_name=out_path+sample_name+'.'+'.'.join([str(i) for i in bps])+'.a.png'
		p = makeDotplot(QC_plot_name, dotdata_qual_check, len(seq2), len(seq2))
		qc_ref_file_write(qc_file_name,sample_name,bps,region_QC,'a',seq2)		
		#right bp flank
		seqs=Pacbio_prodce_ref_alt_short(ref,flank_length,info_b)
		seq2=seqs[0]
		dotdata_qual_check=dotdata(window_size,seq2,seq2)
		region_QC=qual_check_repetitive_region(dotdata_qual_check)
		QC_plot_name=out_path+sample_name+'.'+'.'.join([str(i) for i in bps])+'.b.png'
		p = makeDotplot(QC_plot_name, dotdata_qual_check, len(seq2), len(seq2))
		qc_ref_file_write(qc_file_name,sample_name,bps,region_QC,'b',seq2)		



'vapor svelter --sv-input /scratch/remills_flux/xuefzhao/SV_discovery_index/download/SVelter.version8/results/HG00514.chr9.bam.svelter --pacbio-input /scratch/remills_flux/xuefzhao/SV_discovery_index/smrt_temp/alignment/HG00514.chr9.XXX.bam --reference /scratch/remills_flux/xuefzhao/reference/GRCh38.1KGP/GRCh38_full_analysis_set_plus_decoy_hla.fa --output-path /scratch/remills_flux/xuefzhao/SV_discovery_index/download/SVelter.version8//VaPoR_Vali/HG00514.chr9.bam'




data1=read.table('ref.ref')
data2=read.table('read.ref')
data3=read.table('alt.alt')
data4=read.table('read.alt')
par(mfrow=c(2,2))
plot(data1[,1],data1[,2],pch=22)
plot(data2[,1],data2[,2],pch=22)
plot(data3[,1],data3[,2],pch=22)
plot(data4[,1],data4[,2],pch=22)








