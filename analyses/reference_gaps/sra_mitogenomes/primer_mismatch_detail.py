import sys; sys.argv=[sys.argv[0]]+sys.argv[1:]
exec(open('eval_markers.py').read().split("seqs=read")[0])
seqs=read(sys.argv[1]); s=max(seqs.values(),key=len)
for name,(f,r) in PR.items():
    for S,st in ((s,'+'),(rc(s),'-')):
        mf,pf=best(S,f); mr,pr=best(S,rc(r))
        if 0<pr-pf<900:
            fs=S[pf:pf+len(f)]; rs=rc(S[pr:pr+len(r)])
            mark=lambda p,t:''.join('.' if t[i] in IU[p[i]] else t[i] for i in range(len(p)))
            print(f"{name:9s} F mm{mf}: {mark(f,fs)}   R mm{mr}: {mark(r,rs)}  amp {pr-pf-len(f)}"); break
