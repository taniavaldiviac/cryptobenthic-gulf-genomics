import sys,subprocess,glob,itertools
# Usage: python validate_annotation.py ref_proteins.fa ref_rna.fa  (run from the folder containing asm/)
REF_PROT,REF_RNA=sys.argv[1],sys.argv[2]
B='TCAG'; aa2='FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSS**VVVVAAAADDEEGGGG'
code={a+b+c:aa2[i] for i,(a,b,c) in enumerate(itertools.product(B,B,B))}
comp=str.maketrans('ACGTN','TGCAN'); rc=lambda s:s.translate(comp)[::-1]
def read(f):
    s={};k=None
    for l in open(f):
        l=l.strip()
        if l.startswith('>'):k=l[1:].split()[0];s[k]=''
        else:s[k]+=l.upper()
    return s
def bl(cmd):
    return [l.split('\t') for l in subprocess.run(cmd,capture_output=True,text=True).stdout.strip().split('\n') if l]
print('ensamblaje\tPCG_13\tPCG_sin_stop_interno\tPCG_completos(>=90%)\trRNA\ttRNA_22\tgenes_con_problemas')
for f in sorted(glob.glob('asm/*/*path_sequence.fasta')):
    t=f.split('/')[1]
    if not any(k in t for k in ['_SRR']): continue
    seqs=read(f)
    P=bl(['tblastn','-query',REF_PROT,'-subject',f,'-db_gencode','2','-evalue','1e-10','-outfmt','6 qseqid sseqid qstart qend sstart send qlen pident','-max_hsps','1'])
    genes={};prob=[]
    for q,s,qs,qe,ss,se,ql,pid in P:
        cov=(int(qe)-int(qs)+1)/int(ql)
        if q in genes and genes[q][0]>=cov: continue
        a,b=sorted((int(ss),int(se))); nt=seqs[s][a-1:b]; nt=nt if int(ss)<int(se) else rc(nt)
        prot=''.join(code.get(nt[i:i+3],'X') for i in range(0,len(nt)-2,3))
        genes[q]=(cov,prot[:-1].count('*'))
    ok_stop=sum(1 for g in genes.values() if g[1]==0); full=sum(1 for g in genes.values() if g[0]>=0.9)
    prob=[f'{q}(stop{v[1]})' for q,v in genes.items() if v[1]>0]
    R=bl(['blastn','-query',REF_RNA,'-subject',f,'-task','blastn','-evalue','1e-3','-outfmt','6 qseqid length qlen','-max_hsps','1'])
    rr=len({r[0] for r in R if r[0].startswith('rRNA') and int(r[1])>=0.6*int(r[2])}); tr=len({r[0] for r in R if r[0].startswith('tRNA')})
    print(f"{t}\t{len(genes)}\t{ok_stop}\t{full}\t{rr}/2\t{tr}/22\t{','.join(prob)}")
