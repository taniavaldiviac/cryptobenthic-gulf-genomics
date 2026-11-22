import subprocess,glob,json,re
exec(open('eval_markers.py').read().split("seqs=read")[0].replace("import sys,re","import re"))
SP={'Apogon_retrosella_SRR30361239':'Apogon retrosella','Coryphopterus_urospilus_SRR30564698':'Coryphopterus urospilus','Gillellus_semicinctus_SRR16897090':'Gillellus semicinctus','Lythrypnus_dalli_SRR5170315':'Lythrypnus dalli','Thalassoma_lucasanum_SRR30361224':'Thalassoma lucasanum','Myripristis_leiognathus_SRR30718875':'Myripristis leiognathus','Neoniphon_suborbitalis_SRR30553915':'Neoniphon suborbitalis','Poeciliopsis_latidens_SRR13299717':'Poeciliopsis latidens','Citharichthys_xanthostigma_SRR21844004':'Citharichthys xanthostigma'}
out={}
for t,sp in SP.items():
    f=f'annot/norm/{sp.replace(" ","_")}.fa'; s=max(read(f).values(),key=len); L=len(s); circ=t not in ('Gillellus_semicinctus_SRR16897090',)
    open('/tmp/q.fa','w').write('>q\n'+s+s[:2000]+'\n')
    feats=[]
    P=subprocess.run(['tblastn','-query','annot/ref_prot.fa','-subject','/tmp/q.fa','-db_gencode','2','-evalue','1e-8','-outfmt','6 qseqid sstart send qstart qend qlen bitscore','-max_hsps','1'],capture_output=True,text=True).stdout.split('\n')
    R=subprocess.run(['blastn','-query','annot/ref_rna.fa','-subject','/tmp/q.fa','-task','blastn','-evalue','1e-3','-outfmt','6 qseqid sstart send length qlen bitscore','-max_hsps','1'],capture_output=True,text=True).stdout.split('\n')
    for l in P:
        if not l: continue
        q,a,b,qs,qe,ql,bs=l.split('\t'); a,b=int(a),int(b); strand='+' if a<b else '-'; a,b=min(a,b),max(a,b)
        if a>L: continue
        feats.append(dict(name=q,type='CDS',start=a,end=min(b,L+2000),strand=strand))
    for l in R:
        if not l: continue
        q,a,b,ln,ql,bs=l.split('\t'); a,b=int(a),int(b); strand='+' if a<b else '-'; a,b=min(a,b),max(a,b)
        if a>L: continue
        typ,nm=q.split('_',1)
        if typ=='rRNA' and int(ln)<0.6*int(ql): continue
        feats.append(dict(name=nm.replace('trn','').replace('rrn12','12S').replace('rrn16','16S'),type=typ,start=a,end=b,strand=strand))
    # dedupe tRNA by name keep first
    seen=set(); F=[]
    for x in sorted(feats,key=lambda x:x['start']):
        k=(x['type'],x['name'])
        if k in seen: continue
        seen.add(k); F.append(x)
    amps=[]
    S2=s+s[:2000]
    for n,(Fp,Rp) in PR.items():
        lim=9 if n=='MiFishU' else 4
        mf,pf=best(S2,Fp); mr,pr=best(S2[pf:pf+900],rc(Rp)); pr+=pf
        if mf<=lim and mr<=4 and 0<pr-pf<900: amps.append(dict(name=n,start=pf+1,end=pr+len(Rp),mmF=mf,mmR=mr))
    # rotate so 12S starts near 0? keep raw; record 12S start for alignment
    r12=[x for x in F if x['name']=='12S']
    off=(r12[0]['start']-1) if r12 else 0
    out[sp]=dict(length=L,offset=off,features=F,amplicons=amps,circular=circ)
json.dump(out,open('annot/coords.json','w'),indent=0)
for sp,v in out.items(): print(sp,v['length'],'feat',len(v['features']),'amps',[a['name'] for a in v['amplicons']],'off',v['offset'])
