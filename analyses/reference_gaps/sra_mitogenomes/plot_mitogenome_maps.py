import json, matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, Wedge, Patch
D=json.load(open('annot/coords.json'))
order=['Apogon retrosella','Coryphopterus urospilus','Lythrypnus dalli','Gillellus semicinctus','Thalassoma lucasanum','Myripristis leiognathus','Neoniphon suborbitalis','Poeciliopsis latidens','Citharichthys xanthostigma']
crypto={'Apogon retrosella','Coryphopterus urospilus','Lythrypnus dalli','Gillellus semicinctus'}
def gcol(f):
    n=f['name']
    if f['type']=='rRNA': return '#C44E52'
    if f['type']=='tRNA': return '#4D4D4D'
    if n.startswith('ND'): return '#4C72B0'
    if n.startswith('COX'): return '#DD8452'
    if n.startswith('ATP'): return '#55A868'
    if n=='CYTB': return '#8172B3'
    return '#999999'
MK={'MiFishU':('MiFish-U','#1B9E77'),'MarVer1':('MarVer1','#66A61E'),'Riaz':('Riaz 12S-V5','#E7298A'),'teleo':('teleo','#E6AB02'),'Fish16S':('Fish16S','#7570B3'),'LerayCOI':('Leray COI','#A6761D')}
mk_order=['MiFishU','MarVer1','Riaz','teleo','Fish16S','LerayCOI']
def norm(sp):
    v=D[sp]; sh=v['offset']-70; L=v['length']
    rot=lambda x:(x-sh-1)%L+1
    F=[]
    for f in v['features']:
        a,b=rot(f['start']),rot(min(f['end'],f['start']+ (f['end']-f['start'])))
        if b<a: b=a+(f['end']-f['start'])
        F.append({**f,'start':a,'end':b})
    A=[{**a,'start':rot(a['start']),'end':rot(a['start'])+(a['end']-a['start'])} for a in v['amplicons']]
    return L,F,A
fig,(ax,axz)=plt.subplots(1,2,figsize=(15,7.2),gridspec_kw={'width_ratios':[3.2,1.3]})
for i,sp in enumerate(order):
    y=len(order)-1-i; L,F,A=norm(sp)
    for a_ in (ax,axz):
        a_.add_patch(Rectangle((1,y-0.04),L,0.08,color='#D9D9D9',lw=0))
    for f in F:
        h=0.32 if f['type']!='tRNA' else 0.22
        y0=y+0.05 if f['strand']=='+' else y-0.05-h
        for a_ in (ax,axz):
            a_.add_patch(Rectangle((f['start'],y0),f['end']-f['start'],h,color=gcol(f),lw=0))
        if f['type']!='tRNA' and f['end']-f['start']>500:
            ax.text((f['start']+f['end'])/2,y0+h/2,f['name'],ha='center',va='center',fontsize=6,color='white')
        if f['type']=='rRNA': axz.text((f['start']+f['end'])/2,y0+h/2,f['name'],ha='center',va='center',fontsize=7,color='white')
    cr_start=max(f['end'] for f in F if f['start']<L-50)
    ax.text(L+150,y,f"{L/1000:.1f} kb",va='center',fontsize=7,color='#555555')
    for a_ in A:
        k=mk_order.index(a_['name']); yy=y+0.42+k*0.075
        axz.plot([a_['start'],a_['end']],[yy,yy],color=MK[a_['name']][1],lw=2.2,solid_capstyle='butt')
        if a_['name']=='LerayCOI': ax.plot([a_['start'],a_['end']],[y+0.42,y+0.42],color=MK[a_['name']][1],lw=2.2)
        else: ax.plot([a_['start'],a_['end']],[y+0.42+k*0.06,y+0.42+k*0.06],color=MK[a_['name']][1],lw=2.2)
    lab=sp.split(); txt=f"$\\it{{{lab[0][0]}.\\ {lab[1]}}}$"+(' *' if sp in crypto else '')
    ax.text(-300,y,txt,ha='right',va='center',fontsize=9)
    if sp=='Gillellus semicinctus': ax.text(L+1300,y,'lineal',va='center',fontsize=7,color='#C44E52')
ax.set_xlim(-4200,20500); ax.set_ylim(-0.7,len(order)-0.1); ax.axis('off')
ax.set_title('a) Mitogenomas recuperados de SRA, alineados en tRNA-Phe',loc='left',fontsize=11)
for x in range(0,20000,2000): ax.text(x,-0.62,f'{x//1000}',ha='center',fontsize=7,color='#777777')
ax.text(10000,-0.85,'Posición (kb)',ha='center',fontsize=8,color='#555555')
axz.set_xlim(0,3200); axz.set_ylim(-0.7,len(order)-0.1+0.4); axz.axis('off')
axz.set_title('b) Acercamiento a 12S-16S y amplicones',loc='left',fontsize=11)
for x in range(0,3201,500): axz.text(x,-0.62,f'{x}',ha='center',fontsize=7,color='#777777')
axz.text(1600,-0.85,'Posición (pb)',ha='center',fontsize=8,color='#555555')
leg1=[Patch(color=c,label=l) for l,c in [('ND (complejo I)','#4C72B0'),('COX','#DD8452'),('ATP','#55A868'),('CYTB','#8172B3'),('rRNA','#C44E52'),('tRNA','#4D4D4D')]]
leg2=[Patch(color=MK[k][1],label=MK[k][0]) for k in mk_order]
fig.legend(handles=leg1+leg2,loc='lower center',ncol=12,fontsize=7.5,frameon=False,bbox_to_anchor=(0.5,-0.01))
fig.text(0.01,0.06,'* criptobentónico. Genes de la hebra pesada arriba de la línea; ND6 y tRNA de la hebra ligera abajo. Anotación por transferencia desde Apogon cf. pseudomaculatus (QB112149).',fontsize=7,color='#555555')
plt.tight_layout(rect=(0,0.07,1,1))
for ext in ('png','svg','pdf'): plt.savefig(f'/mnt/user-data/outputs/figuras/mapas_mitogenomas_lineales.{ext}',dpi=300,bbox_inches='tight')
# circular
sp='Apogon retrosella'; L,F,A=norm(sp)
fig2,axc=plt.subplots(figsize=(7.5,7.5)); axc.set_aspect('equal'); axc.axis('off')
ang=lambda p:90-360*(p-1)/L
axc.add_patch(Wedge((0,0),1.0,0,360,width=0.012,color='#BBBBBB'))
for f in F:
    r,w=(1.0,0.14) if f['strand']=='+' else (0.86,0.14)
    if f['type']=='tRNA': w=0.08; r=r if f['strand']=='+' else 0.92
    axc.add_patch(Wedge((0,0),r+(0.0 if f['strand']=='+' else 0),ang(f['end']),ang(f['start']),width=w,color=gcol(f),lw=0))
    if f['type']!='tRNA':
        import math; m=math.radians(ang((f['start']+f['end'])/2)); rr=1.36 if f['strand']=='+' else 0.62
        axc.text(rr*math.cos(m),rr*math.sin(m),f['name'],ha='center',va='center',fontsize=8)
for a_ in A:
    k=mk_order.index(a_['name']); r=1.04+0.03*k
    axc.add_patch(Wedge((0,0),r+0.02,ang(a_['end']),ang(a_['start']),width=0.02,color=MK[a_['name']][1],lw=0))
axc.text(0,0.08,'Apogon retrosella',ha='center',fontsize=12,style='italic'); axc.text(0,-0.07,f'{L:,} pb · circular',ha='center',fontsize=9,color='#555555')
axc.text(0,-0.2,'SRR30361239 · 40× cobertura media',ha='center',fontsize=8,color='#777777')
axc.set_xlim(-1.6,1.6); axc.set_ylim(-1.6,1.6)
axc.legend(handles=leg1+leg2,loc='lower center',ncol=4,fontsize=7,frameon=False,bbox_to_anchor=(0.5,-0.08))
for ext in ('png','svg','pdf'): fig2.savefig(f'/mnt/user-data/outputs/figuras/mapa_circular_Apogon_retrosella.{ext}',dpi=300,bbox_inches='tight')
print('ok')
