"""Generate SVGs: user's drawing tiled, 4.8.8 overlay, Hankin PIC variants. Pure stdlib."""
import math, os, itertools
OUT=os.path.dirname(os.path.abspath(__file__))
R2=math.sqrt(2); T=R2-1   # tan 22.5
def svg(name,w,h,body,vb):
    open(f"{OUT}/{name}.svg","w").write(f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="{vb}" width="{w}" height="{h}"><rect x="-1e3" y="-1e3" width="3e3" height="3e3" fill="#fbf7ee"/>{body}</svg>')
def L(a,b,c="#222",w=0.03): return f'<line x1="{a[0]:.5f}" y1="{-a[1]:.5f}" x2="{b[0]:.5f}" y2="{-b[1]:.5f}" stroke="{c}" stroke-width="{w}" stroke-linecap="round"/>'
def P(p,r=0.06,c="#c0392b"): return f'<circle cx="{p[0]:.5f}" cy="{-p[1]:.5f}" r="{r}" fill="{c}"/>'
def poly(pts,fill="none",c="#222",w=0.03): return '<polygon points="'+" ".join(f"{x:.5f},{-y:.5f}" for x,y in pts)+f'" fill="{fill}" stroke="{c}" stroke-width="{w}" stroke-linejoin="round"/>'

def incentre(A,B,C):
    a=math.dist(B,C);b=math.dist(A,C);c=math.dist(A,B);s=a+b+c
    return tuple((a*A[i]+b*B[i]+c*C[i])/s for i in (0,1))
def cell(i,j):
    """User's drawing in square [2i,2i+2]x[2j,2j+2]: returns segments and dots (VERIFIED incentre by formula)."""
    x,y=2*i,2*j; C=(x+1,y+1); cs=[(x,y),(x+2,y),(x+2,y+2),(x,y+2)]
    segs=[(cs[k],cs[(k+1)%4],'sq') for k in range(4)]+[(cs[k],C,'x') for k in range(4)]
    dots=[]
    for k in range(4):
        A,B=cs[k],cs[(k+1)%4]; d=incentre(A,B,C); dots.append(d)
        segs+= [(d,A,'sp'),(d,B,'sp'),(d,C,'sp')]
    return segs,dots
N=3
def drawing(show_tiling=False):
    body="";alld=[]
    for i in range(N):
        for j in range(N):
            s,d=cell(i,j);alld+=d
            for a,b,k in s:
                col={'sq':'#222','x':'#222','sp':'#1f6f8b'}[k]
                body+=L(a,b,col,0.035 if k!='sp' else 0.03)
    if show_tiling:
        # join dots: edges of 4.8.8 = pairs of dots at distance 2(sqrt2-1)
        e=2*T
        for a,b in itertools.combinations(alld,2):
            if abs(math.dist(a,b)-e)<1e-6: body+=L(a,b,"#c0392b",0.05)
    for d in alld: body+=P(d)
    return body,alld
b,d=drawing(); svg("01_user_drawing_tiled",700,700,b,f"-0.2 {-2*N-0.2} {2*N+0.4} {2*N+0.4}")
b,d=drawing(True); svg("02_user_drawing_with_4.8.8",700,700,b,f"-0.2 {-2*N-0.2} {2*N+0.4} {2*N+0.4}")

# verification: all dots at distance sqrt2-1 from the 3 triangle edges; 4.8.8 edge counts
s,dots=cell(0,0)
d=dots[0]; print("incentre dot",d,"expected",(1,T))
e=2*T
alld=[]
for i in range(N):
    for j in range(N): alld+=cell(i,j)[1]
deg={}
for a,b in itertools.combinations(alld,2):
    if abs(math.dist(a,b)-e)<1e-9: deg[a]=deg.get(a,0)+1; deg[b]=deg.get(b,0)+1
print("interior dot degrees (should be 3):",sorted(set(deg[k] for k in deg if 1<k[0]<2*N-1 and 1<k[1]<2*N-1)))

# ---- Hankin PIC ----
def regular(n,cx,cy,r,rot):
    return [(cx+r*math.cos(rot+2*math.pi*k/n),cy+r*math.sin(rot+2*math.pi*k/n)) for k in range(n)]
def rayint(p,u,q,v):
    det=u[0]*(-v[1])-u[1]*(-v[0])
    if abs(det)<1e-12: return None
    dx,dy=q[0]-p[0],q[1]-p[1]
    t=(dx*(-v[1])-dy*(-v[0]))/det; s=(u[0]*dy-u[1]*dx)/det
    return (t,s)
def pic(poly_pts,theta):
    """Hankin: from each edge midpoint two rays at angle theta (deg) from the edge, into polygon. Each ray is cut where it meets the ray from the adjacent edge. Returns segments."""
    n=len(poly_pts); mids=[];dirs=[]
    cx=sum(p[0] for p in poly_pts)/n;cy=sum(p[1] for p in poly_pts)/n
    segs=[]; th=math.radians(theta+0.001)
    rays=[]  # per edge: (mid, dir toward next vertex, dir toward prev vertex)
    for k in range(n):
        A,B=poly_pts[k],poly_pts[(k+1)%n]; m=((A[0]+B[0])/2,(A[1]+B[1])/2)
        e=(B[0]-A[0],B[1]-A[1]); l=math.hypot(*e); e=(e[0]/l,e[1]/l)
        # inward normal
        nrm=(-e[1],e[0]); 
        if nrm[0]*(cx-m[0])+nrm[1]*(cy-m[1])<0: nrm=(-nrm[0],-nrm[1])
        fw=(e[0]*math.cos(th)+nrm[0]*math.sin(th), e[1]*math.cos(th)+nrm[1]*math.sin(th))
        bw=(-e[0]*math.cos(th)+nrm[0]*math.sin(th), -e[1]*math.cos(th)+nrm[1]*math.sin(th))
        rays.append((m,fw,bw))
    for k in range(n):
        m1,f1,_=rays[k]; m2,_,b2=rays[(k+1)%n]
        r=rayint(m1,f1,m2,b2)
        if r and r[0]>0 and r[1]>0:
            q=(m1[0]+f1[0]*r[0],m1[1]+f1[1]*r[0]); segs+=[(m1,q),(m2,q)]
    return segs
def tiling488(nx,ny):
    """regular 4.8.8: octagons at 2*(i,j) centred on user's square corners, squares at centres."""
    ro=1/math.cos(math.pi/8)  # circumradius for inradius 1
    polys=[]
    for i in range(nx+1):
        for j in range(ny+1): polys.append(regular(8,2*i,2*j,ro,math.pi/8))
    rs=(2*T)/R2  # small square circumradius: side 2T
    for i in range(nx):
        for j in range(ny): polys.append(regular(4,2*i+1,2*j+1,rs,math.pi/4*0+math.pi/2*0+math.pi/4*0+math.pi/4*0))
    return polys
# small square: dots at (1,T),(2-T,1)... axis-aligned diamond => rotation 90 deg offset: vertices at angles 270,0,... => rot=0 with circumradius=1-T
def tiling488(nx,ny):
    ro=1/math.cos(math.pi/8); polys=[]
    for i in range(nx+1):
        for j in range(ny+1): polys.append(regular(8,2*i,2*j,ro,math.pi/8))
    for i in range(nx):
        for j in range(ny): polys.append(regular(4,2*i+1,2*j+1,1-T,0))
    return polys
def tetrakis(nx,ny):
    polys=[]
    for i in range(nx):
        for j in range(ny):
            x,y=2*i,2*j;c=(x+1,y+1);cs=[(x,y),(x+2,y),(x+2,y+2),(x,y+2)]
            for k in range(4): polys.append([cs[k],cs[(k+1)%4],c])
    return polys
def pic_svg(name,polys,theta,nx,ny,show=True,col="#8e2b1f"):
    body=""
    if show:
        for p in polys: body+=poly(p,"#efe3c8" if len(p)==8 else ("#d9e6ea" if len(p)==4 else "#f3ecdc"),"#999",0.02)
    seen=set()
    for p in polys:
        for a,b in pic(p,theta):
            key=tuple(sorted([(round(a[0],4),round(a[1],4)),(round(b[0],4),round(b[1],4))]))
            if key in seen: continue
            seen.add(key); body+=L(a,b,col,0.04)
    svg(name,600,600,body,f"0 {-2*ny} {2*nx} {2*ny}")
for th in (22.5,30,45,60,67.5):
    pic_svg(f"03_hankin_4.8.8_theta{th}",tiling488(3,3),th,3,3)
for th in (22.5,45,60,67.5):
    pic_svg(f"04_hankin_tetrakis_theta{th}",tetrakis(3,3),th,3,3)
# user's own spokes on 4.8.8: vertex-to-centre radii of each polygon, plus 4.8.8 edges
body=""
for p in tiling488(3,3):
    n=len(p);cx=sum(q[0] for q in p)/n;cy=sum(q[1] for q in p)/n
    body+=poly(p,"none","#c0392b",0.05)
    for q in p: body+=L((cx,cy),q,"#1f6f8b",0.03)
svg("05_488_with_vertex_spokes",600,600,body,f"0 -6 6 6")
print("done")
# zoomed single-cell views for checking
for th in (30,45,67.5):
    pic_svg(f"06_zoom_4.8.8_theta{th}",tiling488(2,2),th,2,2)
def sheet(name,mk,thetas,nx=2,ny=2,cols=4):
    rows=(len(thetas)+cols-1)//cols; W=2*nx+0.5
    body=""
    for idx,th in enumerate(thetas):
        ox=(idx%cols)*W; oy=(idx//cols)*(2*ny+0.8)
        g=""; seen=set()
        for p in mk(nx,ny):
            for a,b in pic(p,th):
                key=tuple(sorted([(round(a[0],4),round(a[1],4)),(round(b[0],4),round(b[1],4))]))
                if key in seen: continue
                seen.add(key); g+=L(a,b,"#8e2b1f",0.05)
        body+=f'<g transform="translate({ox},{oy})"><clipPath id="c{idx}"><rect x="0" y="{-2*ny}" width="{2*nx}" height="{2*ny}"/></clipPath><g clip-path="url(#c{idx})">{g}</g><rect x="0" y="{-2*ny}" width="{2*nx}" height="{2*ny}" fill="none" stroke="#999" stroke-width="0.03"/><text x="0.1" y="0.5" font-size="0.4" font-family="sans-serif">theta={th}</text></g>'
    svg(name,1200,300*rows,body,f"-0.2 {-2*ny-0.2} {cols*W} {rows*(2*ny+0.8)+0.2}")
sheet("07_sheet_4.8.8",tiling488,[15,22.5,30,35,40,50,55,60,67.5,72,78,85])
sheet("08_sheet_tetrakis",tetrakis,[15,22.5,30,35,40,50,55,60,67.5,72,78,85])
