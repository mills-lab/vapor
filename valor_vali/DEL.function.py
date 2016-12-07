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

