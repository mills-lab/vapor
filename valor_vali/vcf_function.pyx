import os
def case_hash_readin(filein,flank_length):
	out=[]
	fin=open(filein)
	rec=-1
	for line in fin:
		pin=line.strip().split()
		if not pin[0][0]=='#': 
			rec+=1
			SVtype=SVtype_decide(pin)
			if SVtype=='DEL' or 'DEL' in pin[4]: 
				gt_temp=Genotype_caller(pin)
				ka='a/a'
				kb=DEL_alt_make(gt_temp)
				#if ka==kb: continue
				sv_start=int(pin[1])
				sv_end=SV_end_read(pin)
				if sv_end>sv_start:
					out.append([pin[0],sv_start,sv_end,ka,kb,rec])					
			elif SVtype=='INV':
				gt_temp=Genotype_caller(pin)
				ka='a/a'
				kb=INV_alt_make(gt_temp)
				#if ka==kb: continue
				sv_start=int(pin[1])
				sv_end=SV_end_read(pin)
				if sv_end>sv_start:
					out.append([pin[0],sv_start,sv_end,ka,kb,rec])
			elif SVtype=='DUP' or 'DUP' in pin[4]:
				gt_temp=Genotype_caller(pin)
				ka='a/a'
				kb=DUP_TAN_alt_make(gt_temp)
				#if ka==kb: continue
				sv_start=int(pin[1])
				sv_end=SV_end_read(pin)
				if sv_end>sv_start:
					out.append([pin[0],sv_start,sv_end,ka,kb,rec])
			elif 'INS' in SVtype or 'INS' in pin[4] or SVtype in ['ALU','LINE1','SVA']:
				ka='/'
				gt_temp=Genotype_caller(pin)
				kb=INS_alt_make(gt_temp)
				sv_start=int(pin[1])
				SV_length=SVlength_decide(pin)
				out.append([pin[0],sv_start,SV_length+sv_start,ka,kb,rec])
			elif SVtype=='CNV':
				gt_temp=Genotype_caller(pin)
				ka='a/a'
				kb=CNV_alt_make(gt_temp,pin)
				#if ka==kb: continue
				sv_start=int(pin[1])
				sv_end=SV_end_read(pin)
				if sv_end>sv_start:
					out.append([pin[0],sv_start,sv_end,ka,kb,rec])
			elif '[' in pin[4] or ']' in pin[4]:
				if ']' in pin[4] and pin[4].split(']')[0]=='':
					#eg:']13:123456]AGTNNNNNCAT'
					info=[pin[0],pin[1],pin[4].split(']')[1].split(':')[0],
							int(pin[4].split(']')[1].split(':')[1])-flank_length,
							int(pin[4].split(']')[1].split(':')[1]),
							pin[4].split(']')[2][:-1],
							pin[0],int(pin[1]),int(pin[1])+flank_length,'.','.']
					out.append(info+[rec])
				elif '[' in pin[4] and pin[4].split('[')[2]=='':
					#eg:'CAGTNNNNNCA[2:321682['
					info=[pin[0],pin[1],pin[0],int(pin[1])-flank_length,int(pin[1]),
							pin[4].split('[')[0][1:],
							pin[4].split('[')[1].split(':')[0],
							int(pin[4].split('[')[1].split(':')[1]),
							int(pin[4].split('[')[1].split(':')[1])+flank_length,'.','.']
					out.append(info+[rec])
				elif ']' in pin[4] and pin[4].split(']')[2]=='':
					info=[pin[0],pin[1],pin[0],int(pin[1])-flank_length,int(pin[1]),
							pin[4].split(']')[0][1:],
							pin[4].split(']')[1].split(':')[0],
							int(pin[4].split(']')[1].split(':')[1])-flank_length,
							int(pin[4].split(']')[1].split(':')[1]),'.','^']
					out.append(info+[rec])
				elif '[' in pin[4] and pin[4].split('[')[0]=='':
					info=[pin[0],pin[1],pin[4].split('[')[1].split(':')[0],
							int(pin[4].split('[')[1].split(':')[1]),
							int(pin[4].split('[')[1].split(':')[1])+flank_length,
							pin[4].split('[')[2][:-1],
							pin[0],int(pin[1]),int(pin[1])+flank_length,'^','.']
					out.append(info+[rec])
				else:
					print 'unconsidered situation'
					print pin
			elif pin[4]=='<INS>':continue
			elif pin[4]=='<CNV>':continue #for now. modify later
			elif pin[4]=='<DUP>':continue #for now. modify later
			elif pin[4]=='<INS:ME>':continue #for now. modify later
			#else:
			#SV_type=SVtype_decide(pin)
			else:
					print 'unconsidered situation'
					print pin
	fin.close()
	return out
def case_hash_unify(case_hash):
	out=[]
	rec=[]
	for x in case_hash:
		if not '.' in x and not '^' in x:
			if not x in out:
				out.append(x)
		else:
			if not x[2:-1] in rec:
				rec.append(x[2:-1])
				out.append(x)
	return out 
def cigar2alignstart_left(cigar,start,bps,flank_length):
	#eg cigar2alignstart(pbam[5],int(pbam[3]),bps)
	import re
	pcigar=re.compile(r'''(\d+)([MIDNSHP=X])''')
	cigars=[]
	for m in pcigar.finditer(cigar):
		cigars.append((m.groups()[0],m.groups()[1]))
	read_rec=0
	align_rec=start
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
		if align_rec>int(bps[1])-flank_length: break
	return [read_rec,int(align_rec)-int(bps[1])+flank_length]
def CNV_alt_make(gt_temp,pin):
	alt_candidates=['<CN1>']+pin[4].split(',')
	gt_temp=Genotype_caller(pin)
	gt_info=[alt_candidates[x] for x in gt_temp]
	copynumber_info=[int(i.replace('<CN','').replace('>','')) for i in gt_info]
	kb='/'.join([''.join(['a' for i in range(copynumber_info[0])]),''.join(['a' for i in range(copynumber_info[1])])])
	return kb
def DEL_alt_make(genotype):
	#eg of genotype:[0,0],[0,1],[1,1]
	if type(genotype[0])==type(1) and type(genotype[1])==type(1):
		if sum(genotype)==0:
			kb='a/a'
		elif sum(genotype)==1:
			kb='/a'
		elif sum(genotype)==2:
			kb='/'
		else:
			kb='error'
		return kb
	else:
		print 'wrong genotype!'
		return 'error'
def DUP_TAN_alt_make(genotype):
	#eg of genotype:[0,0],[0,1],[1,1]
	if type(genotype[0])==type(1) and type(genotype[1])==type(1):
		if sum(genotype)==0:
			kb='a/a'
		elif sum(genotype)==1:
			kb='aa/a'
		elif sum(genotype)==2:
			kb='aa/aa'
		else:
			kb='error'
		return kb
	else:
		print 'wrong genotype!'
		return 'error'
def Genotype_caller(pin):
	gt_pos=pin[8].split(':').index('GT')
	gt_out=pin[9].split(':')[gt_pos]
	if '/' in gt_out:
		if gt_out.split('/')[0] in ['0','1'] and gt_out.split('/')[1] in ['0','1']:
			return sorted([int(i) for i in gt_out.split('/')])
		else:
			#if GT look like './.', no idea how to interprete for now.
			return [0,0]
	elif '|' in gt_out:
		if gt_out.split('|')[0] in ['0','1'] and gt_out.split('|')[1] in ['0','1']:
			return sorted([int(i) for i in gt_out.split('|')])
		else:
			#if GT look like './.', no idea how to interprete for now.
			return [0,0]
def genotype_decide(pin,sample_pos):
	gt=pin[sample_pos]
	if '/' in gt:
		gt_new=gt.split('/')
	elif '|' in gt:
		gt_new=gt.split('|')
	else:
		gt_new=[0,0]
	gt_new=[int(i) for i in gt_new]
	return gt_new
def INS_alt_make(genotype):
	if type(genotype[0])==type(1) and type(genotype[1])==type(1):
		if sum(genotype)==0:
			kb='/'
		elif sum(genotype)==1:
			kb='a/'
		elif sum(genotype)==2:
			kb='a/a'
		else:
			kb='error'
		return kb
	else:
		print 'wrong genotype!'
		return 'error'		
def Insertion_ref_recognize(ref_letter):
	#return 0 if ref letter, 1 else
	#eg of ref_letter: abc/abc ,  ac/ac
	#eg of normal ref: abc/abc
	#eg of abnormal ref: ac/ab
	out=1
	if ref_letter.split('/')[0]==ref_letter.split('/')[1]:
		test_a=[ord(x) for x in ref_letter.split('/')[0]]
		test_b=[test_a[x+1]-test_a[x] for x in range(len(test_a)-1)]
		test_c=[i for i in test_b if i>1]
		if test_c==[]:
			out=0
	return out
def INV_alt_make(genotype):
	#eg of genotype:[0,0],[0,1],[1,1]
	if type(genotype[0])==type(1) and type(genotype[1])==type(1):
		if sum(genotype)==0:
			kb='a/a'
		elif sum(genotype)==1:
			kb='a^/a'
		elif sum(genotype)==2:
			kb='a^/a^'
		else:
			kb='error'
		return kb
	else:
		print 'wrong genotype!'
		return 'error'
def name_hash_modify(name_hash):
	out=[]
	for x in name_hash:
		if '/' in x:
			y=x.replace('/','-')
		else:
			y=x
		out.append(y)
	return out
def PacVal_score_hash_modify(PacVal_score_hash):
	out={}
	for k1 in PacVal_score_hash.keys():
		key1=str(k1)
		if not int(key1.split('.')[0]) in out.keys():
			out[int(key1.split('.')[0])]={}
		if not k1 in out[int(key1.split('.')[0])].keys():
			out[int(key1.split('.')[0])][k1]=PacVal_score_hash[k1]
	return out
def SVtype_decide(pin):
	sv_type=''
	for x in pin[7].split(';'):
		if x.split('=')[0]=='SVTYPE':
			sv_type=x.split('=')[1]
	return sv_type
def SVend_decide(pin):
	end=int(pin[1])
	for x in pin[7].split(';'):
		if x.split('=')[0]=='END':
			end=int(x.split('=')[1])
	if end>int(pin[1]):
		return end
	else:
		print 'error: SVend < SVstart'
		return end
def SVlength_decide(pin):
	length=0
	for x in pin[7].split(';'):
		if x.split('=')[0]=='SVLEN':
			length=int(x.split('=')[1])
	if length>0:
		return length
	else: 
		print 'error: SVlength = 0'
		return length
def SV_end_read(pin):
	end=0
	for x in pin[7].split(';'):
		if x.split('=')[0]=='END':
			end=int(x.split('=')[1])
	return end
def SVtype_decide(pin):
	sv_type=''
	for x in pin[7].split(';'):
		if x.split('=')[0]=='SVTYPE':
			sv_type=x.split('=')[1]
	return sv_type
def write_PacVal_score_hash(PacVal_file_out,PacVal_file_in,PacVal_score_hash):
	PacVal_score_new=PacVal_score_hash_modify(PacVal_score_hash)
	fo=open(PacVal_file_out,'w')
	rec=-1
	fin=open(PacVal_file_in)
	for line in fin:
		pin=line.strip().split()
		if not pin[0]=='chr': 
			rec+=1
			if rec in PacVal_score_new.keys():
				if len(PacVal_score_new[rec].keys())>0:
					for k2 in PacVal_score_new[rec].keys():
						pin2=pin+[str(PacVal_score_new[rec][k2])]
						print >>fo, ' '.join([str(i) for i in pin2])
				else:
					pin2=pin+['-1']
					print >>fo, ' '.join([str(i) for i in pin2])
			else:
				pin2=pin+['-1'] #SV locates in repetitive regions, cannot validate
				print >>fo, ' '.join([str(i) for i in pin2])
	fin.close()
	fo.close()