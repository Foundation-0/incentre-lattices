"""Triangle centres: where incentre (= medial-axis branch point = Chebyshev centre = pole of inaccessibility) and centroid diverge. stdlib only."""
import math
def cen(A,B,C):
    a=math.dist(B,C);b=math.dist(A,C);c=math.dist(A,B);s=a+b+c
    inc=tuple((a*A[i]+b*B[i]+c*C[i])/s for i in (0,1))
    cg=tuple((A[i]+B[i]+C[i])/3 for i in (0,1))
    return inc,cg
tris={"equilateral":((0,0),(1,0),(.5,math.sqrt(3)/2)),
"right isosceles (user)":((0,0),(2,0),(1,1)),
"tall isosceles 1:3":((0,0),(1,0),(.5,3)),
"flat isosceles 1:0.2":((0,0),(1,0),(.5,.2)),
"obtuse 120deg apex":((0,0),(2,0),(1,math.tan(math.radians(30)))),
"scalene 3-4-5":((0,0),(4,0),(0,3))}
rows=[];svg=""
for k,(n,(A,B,C)) in enumerate(tris.items()):
    inc,cg=cen(A,B,C); per=math.dist(A,B)+math.dist(B,C)+math.dist(C,A)
    d=math.dist(inc,cg); rows.append((n,d/per*100,d))
    sc=90/ max(max(p[0] for p in (A,B,C)),max(p[1] for p in (A,B,C))); ox=(k%3)*120+15; oy=(k//3)*120+105
    f=lambda p:(ox+p[0]*sc,oy-p[1]*sc)
    svg+='<polygon points="'+" ".join("%.1f,%.1f"%f(p) for p in (A,B,C))+'" fill="#efe3c8" stroke="#333"/>'
    svg+='<circle cx="%.1f" cy="%.1f" r="3.5" fill="#c0392b"/><circle cx="%.1f" cy="%.1f" r="3.5" fill="#1f6f8b"/>'%(*f(inc),*f(cg))
    svg+=f'<text x="{ox}" y="{oy+10}" font-size="8" font-family="sans-serif">{n}</text>'
open("09_triangle_centres.svg","w").write('<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 380 250" width="760" height="500"><rect width="380" height="250" fill="#fbf7ee"/>'+svg+'<text x="10" y="240" font-size="9" font-family="sans-serif"><tspan fill="#c0392b">red = incentre</tspan>  <tspan fill="#1f6f8b">blue = centroid</tspan></text></svg>')
for n,p,d in rows: print(f"{n:28s} |incentre-centroid| = {p:5.2f}% of perimeter")
