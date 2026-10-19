/* fold-core.js: exact finite counter-rotation mechanism of the dot/ball contact framework.
   Port of research/mechanics/mechanism.py (build_state, positions) and final_states.py (contacts, density).
   Clusters (squares in 2D, octahedra in 3D) sit on a checkerboard; even clusters are rotated by R, odd by R^T.
   Cross bars keep length L = 2 - 2t; the lattice steps are chosen on the axes (alpha_i). */
(function(root){
"use strict";
const S2=Math.SQRT2, TSTAR=2-S2;
function rot2(th){const c=Math.cos(th),s=Math.sin(th);return [[c,-s],[s,c]]}
function rodrigues(axis,th){ // rotation by th about unit axis (right-handed), as scipy Rotation.from_rotvec
  const n=Math.hypot(axis[0],axis[1],axis[2]),x=axis[0]/n,y=axis[1]/n,z=axis[2]/n,c=Math.cos(th),s=Math.sin(th),C=1-c;
  return [[c+x*x*C,x*y*C-z*s,x*z*C+y*s],[y*x*C+z*s,c+y*y*C,y*z*C-x*s],[z*x*C-y*s,z*y*C+x*s,c+z*z*C]];
}
function T(M){return M[0].map((_,j)=>M.map(r=>r[j]))}
function mv(M,v){return M.map(r=>r.reduce((a,x,i)=>a+x*v[i],0))}
/* steps: lattice spacing along each axis. returns null if the mechanism cannot close at this angle */
function steps(n,Rm,t){
  const L=2-2*t;
  // column i of t(R+R^T): s_i[k]=t*(R[k][i]+R[i][k])
  const al=[];
  for(let i=0;i<n;i++){
    let perp2=0,si=0;
    for(let k=0;k<n;k++){const v=t*(Rm[k][i]+Rm[i][k]);if(k===i)si=v;else perp2+=v*v}
    const need=L*L-perp2;if(need<-1e-12)return null;
    al.push(si+Math.sqrt(Math.max(0,need)));
  }
  return al;
}
/* build the cluster block. n=2 or 3. half: clusters at indices -half..half */
function build(n,Rm,t,half){
  const al=steps(n,Rm,t);if(!al)return null;
  const Rt=T(Rm),L=2-2*t,balls=[],clusters=[],idx={};
  const rng=[];for(let i=-half;i<=half;i++)rng.push(i);
  const cells=[];
  if(n===2)for(const i of rng)for(const j of rng)cells.push([i,j]);else for(const i of rng)for(const j of rng)for(const k of rng)cells.push([i,j,k]);
  for(const c of cells){
    const par=((c.reduce((a,b)=>a+b,0)%2)+2)%2,Rc=par?Rt:Rm,ctr=c.map((v,i)=>v*al[i]);
    clusters.push({c,par,ctr});
    for(let i=0;i<n;i++)for(const sg of [1,-1]){
      const v=new Array(n).fill(0);v[i]=sg*t;const w=mv(Rc,v);
      idx[c.join(',')+'|'+i+'|'+sg]=balls.length;
      balls.push({c,par,axis:i,sg,p:ctr.map((x,k)=>x+w[k])});
    }
  }
  // bonds: inside a cluster (all non-antipodal pairs) and cross bars (+e_i of c to -e_i of c+e_i)
  const bonds=[];
  for(const c of cells){
    const ids=[];for(let i=0;i<n;i++)for(const sg of [1,-1])ids.push([i,idx[c.join(',')+'|'+i+'|'+sg]]);
    for(let a=0;a<ids.length;a++)for(let b=a+1;b<ids.length;b++)if(ids[a][0]!==ids[b][0])bonds.push([ids[a][1],ids[b][1],'in']);
    for(let i=0;i<n;i++){const c2=c.slice();c2[i]+=1;const k2=c2.join(',')+'|'+i+'|-1';if(k2 in idx)bonds.push([idx[c.join(',')+'|'+i+'|1'],idx[k2],'x'])}
  }
  // ball radius: half the shorter bar, sqrt2 t inside a cluster or L across (as lock2 in results/fold-2d/fold2d.py)
  return {n,al,t,L,rho:Math.min(S2*t,L)/2,balls,clusters,bonds,idx,half};
}
function dist(a,b){let s=0;for(let k=0;k<a.length;k++)s+=(a[k]-b[k])*(a[k]-b[k]);return Math.sqrt(s)}
/* contacts of the central cluster's balls (ball-ball distance equals 2*rho within tol), and the density */
function analyse(st,tol){
  tol=tol||1e-9;const n=st.n,zero=new Array(n).fill(0).join(',');
  const cen=[];for(let i=0;i<n;i++)for(const sg of [1,-1])cen.push(st.idx[zero+'|'+i+'|'+sg]);
  const counts=cen.map(a=>{let c=0;for(let b=0;b<st.balls.length;b++){if(b===a)continue;if(Math.abs(dist(st.balls[a].p,st.balls[b].p)-2*st.rho)<tol)c++}return c});
  let maxErr=0;for(const [a,b] of st.bonds){const e=Math.abs(dist(st.balls[a].p,st.balls[b].p)-st.L);if(e>maxErr)maxErr=e}
  const vol=st.al.reduce((a,b)=>a*b,1);
  const ballVol=n===2?Math.PI*st.rho*st.rho:4/3*Math.PI*Math.pow(st.rho,3);
  const strain=st.al.map(a=>a/2-1);
  return {counts,maxErr,density:2*n*ballVol/vol,strain,volRatio:vol/Math.pow(2,n)};
}
const api={S2,TSTAR,rot2,rodrigues,steps,build,analyse,dist,T};
if(typeof module!=='undefined'&&module.exports)module.exports=api;else root.FOLD=api;
})(typeof window!=='undefined'?window:this);
