def svelter_read_in(filein):
	out={}
	fin=open(filein)
	pin=fin.readline().strip().split()
	rec=-1
	for line in fin:
		pin=line.strip().split()
		rec+=1
		if not pin[4] in list(out.keys()):
			out[pin[4]]={}
		if not pin[5] in list(out[pin[4]].keys()):
			out[pin[4]][pin[5]]=[]
		if not pin[3].split(':') in out[pin[4]][pin[5]]:
			out[pin[4]][pin[5]].append(pin[3].split(':')+[rec])
	fin.close()
	return out
def PacVal_score_hash_modify(PacVal_score_hash):
	out={}
	for k1 in list(PacVal_score_hash.keys()):
		key1=str(k1)
		if not int(key1.split('.')[0]) in list(out.keys()):
			out[int(key1.split('.')[0])]={}
		if not k1 in list(out[int(key1.split('.')[0])].keys()):
			out[int(key1.split('.')[0])][k1]=PacVal_score_hash[k1]
	return out
def write_PacVal_score_hash(PacVal_file_in,PacVal_file_out,PacVal_score_hash_new,PacVal_score_hash):
	PacVal_score_hash_2=PacVal_score_hash_modify(PacVal_score_hash)
	PacVal_score_hash_new_2=PacVal_score_hash_modify(PacVal_score_hash_new)
	fo=open(PacVal_file_out,'w')
	print(' '.join(['chr','start','end','bp_info','ref_structure','alt_structure','svelter_score','vapor_score','vapor_genotype']), file=fo)
	rec=-1
	fin=open(PacVal_file_in)
	for line in fin:
		pin=line.strip().split()
		if not pin[0]=='chr': 
			rec+=1
			if rec in list(PacVal_score_hash_2.keys()) and rec in list(PacVal_score_hash_new_2.keys()):
				if len(list(PacVal_score_hash_2[rec].keys()))>0 and len(list(PacVal_score_hash_new_2[rec].keys()))>0:
					for k2 in list(PacVal_score_hash_2[rec].keys()):
						if k2 in list(PacVal_score_hash_new_2[rec].keys()):
							pin2=pin+[str(PacVal_score_hash_new_2[rec][k2]),str(PacVal_score_hash_2[rec][k2])]
							print(' '.join([str(i) for i in pin2]), file=fo)
				else:
					pin2=pin+['-1','-1']
					print(' '.join([str(i) for i in pin2]), file=fo)
			else:
				pin2=pin+['-1','-1'] #SV locates in repetitive regions, cannot validate
				print(' '.join([str(i) for i in pin2]), file=fo)
	fin.close()
	fo.close()

