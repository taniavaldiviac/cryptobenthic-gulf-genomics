import sys,re
PR={'MiFishU':('GTCGGTAAAACTCGTGCCAGC','CATAGTGGGGTATCTAATCCCAGTTTG'),'teleo':('ACACCGCCCGTCACTCT','CTTCCGGTACACTTACCATG'),
'Riaz':('ACTGGGATTAGATACCCC','TAGAACAGGCTCCTCTAG'),'MarVer1':('CGTGCCAGCCACCGCG','GGGTATCTAATCCYAGTTTG'),
'Fish16S':('GACCCTATGGAGCTTTAGAC','CGCTGTTATCCCTADRGTAACT'),'LerayCOI':('GGWACWGGWTGAACWGTWTAYCCYCC','TAIACYTCIGGRTGICCRAARAAYCA')}
IU={'A':'A','C':'C','G':'G','T':'T','R':'AG','Y':'CT','S':'CG','W':'AT','K':'GT','M':'AC','B':'CGT','D':'AGT','H':'ACT','V':'ACG','N':'ACGT','I':'ACGT'}
comp=str.maketrans('ACGTRYSWKMBDHVNI','TGCAYRSWMKVHDBNI')
rc=lambda s:s.translate(comp)[::-1]
def best(seq,p):
    L=len(p);b=(99,-1)
    for i in range(len(seq)-L+1):
        mm=0
        for j in range(L):
            if seq[i+j] not in IU[p[j]]:
                mm+=1
                if mm>=b[0]:break
        if mm<b[0]:b=(mm,i)
    return b
def read(f):
    s={};k=None
    for l in open(f):
        l=l.strip()
        if l.startswith('>'):k=l[1:];s[k]=''
        else:s[k]+=l.upper()
    return s
seqs=read(sys.argv[1]);tot=sum(len(v) for v in seqs.values());seqs={k:(v+v[:1500] if len(v)>15000 else v) for k,v in seqs.items()};
print(f"{len(seqs)} secuencias, {tot} pb")
out={}
for name,(f,r) in PR.items():
    res='no'
    for sid,s in seqs.items():
        for S in (s,rc(s)):
            mf,pf=best(S,f); mr,pr=best(S[pf:pf+800],rc(r)); pr+=pf
            if mf<=4 and mr<=4 and 0<pr-pf<800:
                amp=S[pf+len(f):pr]; res=f"{len(amp)} pb (mm F{mf}/R{mr})"; out[name]=amp;break
        if res!='no':break
    print(f"  {name}: {res}")
if len(sys.argv)>2:
    with open(sys.argv[2],'w') as h:
        for k,v in out.items(): h.write(f">{k}\n{v}\n")
