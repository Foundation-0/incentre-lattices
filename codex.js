/* codex.js: shared shell for the geometry site. Offline, no dependencies.
   Provides the index panel, minimisable panels, H to hide the HUD, scanlines/grain/backdrop,
   and small helpers for canvas figures (CX.canvas, CX.pointer, CX.cam3). */
(function(){
"use strict";
const PAGES=[
 ['lattice.html','LATTICE','hub'],
 ['origin.html','THE DRAWING','origin'],
 ['mirrors.html','MIRRORS','why'],
 ['matter.html','MATTER','crystals'],
 ['motion.html','MOTION','the fold'],
 ['mind.html','PERCEPTION','mind'],
 ['open.html','LEDGER','claims']
];
const EXTRA=[['cube-six.html','the cube'],['deeper.html','the full story']];
const Q=new URLSearchParams(location.search);
const here=(location.pathname.split('/').pop()||'lattice.html');
const C={g:'#37ff9a',c:'#35e6ff',a:'#ffb23a',r:'#ff4a3d',dim:'#2f8f68',bg:'#02090c',ink:'#eafff4'};
const TAU=Math.PI*2;
const CX=window.CX={C,TAU,Q};

function el(tag,attrs,html){const e=document.createElement(tag);if(attrs)for(const k in attrs)e.setAttribute(k,attrs[k]);if(html!=null)e.innerHTML=html;return e}

/* ---- decorative layers ---- */
function decor(){
  const R=Math.SQRT2-1;
  const s=el('div',{id:'sig'});
  const pt=(r,i)=>[(r*Math.cos(i*TAU/8-Math.PI/8)).toFixed(1),(r*Math.sin(i*TAU/8-Math.PI/8)).toFixed(1)];
  let h='<svg xmlns="http://www.w3.org/2000/svg" viewBox="-500 -500 1000 1000" width="100%" height="100%"><g class="rot" fill="none" stroke="#37ff9a" stroke-width="1">';
  [460,460*R,460*R*R,460*R*R*R].forEach((r,i)=>{h+='<circle r="'+r.toFixed(1)+'" opacity="'+[.35,.3,.25,.2][i]+'"/>'});
  let o='',s8='';for(let i=0;i<8;i++){o+=(i?'L':'M')+pt(460,i).join(' ');s8+=(i?'L':'M')+pt(460,(i*3)%8).join(' ')}
  h+='<path d="'+o+'Z" opacity=".3"/><path d="'+s8+'Z" stroke="#35e6ff" opacity=".18"/></g></svg>';
  s.innerHTML=h;document.body.appendChild(s);
  document.body.appendChild(el('div',{id:'scan'}));
  const g=el('div',{id:'grain'}),cc=document.createElement('canvas');cc.width=cc.height=128;
  const c=cc.getContext('2d'),im=c.createImageData(128,128);
  for(let i=0;i<im.data.length;i+=4){const v=Math.random()*255|0;im.data[i]=v*.4;im.data[i+1]=v;im.data[i+2]=v*.7;im.data[i+3]=255}
  c.putImageData(im,0,0);g.style.backgroundImage='url('+cc.toDataURL()+')';document.body.appendChild(g);
}

/* ---- index panel ---- */
function nav(){
  const p=el('nav',{'class':'panel hud',id:'p-nav'});
  let b='';
  for(const [u,n,d] of PAGES)b+='<a href="'+u+'"'+(u===here?' class="cur"':'')+'>'+(u===here?'&gt; ':'')+n+'<span>'+d+'</span></a>';
  b+='<hr><div class="foot">CLASSIC // simple, clean reading</div>';
  for(const [u,n] of EXTRA)b+='<a href="'+u+'"'+(u===here?' class="cur"':'')+'>'+n+'<span>classic</span></a>';
  b+='<div class="foot">[H] hide HUD</div>';
  p.innerHTML='<header><span class="tabtitle">INDEX // R=&radic;2&minus;1</span><button class="mb" title="minimise">&ndash;</button></header><div class="body">'+b+'</div>';
  document.body.insertBefore(p,document.body.firstChild);
  if((Q.get('min')||'').split(',').includes('index')||Q.get('min')==='all')p.classList.add('min');
}

function panels(){
  const mins=(Q.get('min')||'').split(',');
  for(const p of document.querySelectorAll('.panel')){
    const btn=p.querySelector('.mb');if(!btn||btn._w)continue;btn._w=1;
    const sync=()=>{btn.innerHTML=p.classList.contains('min')?'+':'&ndash;'};
    if(Q.get('min')==='all'||(p.id&&mins.includes(p.id.replace(/^p-/,''))))p.classList.add('min');
    sync();btn.onclick=()=>{p.classList.toggle('min');sync()};
  }
}
CX.panels=panels;

addEventListener('keydown',e=>{
  if(e.ctrlKey||e.metaKey||e.altKey)return;
  const t=e.target&&e.target.tagName;if(t==='INPUT'&&e.target.type==='text')return;
  if(e.key==='h'||e.key==='H')document.body.classList.toggle('hud-off');
});

/* ---- canvas helper ---- */
CX.canvas=function(canvas,draw,opt){
  opt=opt||{};const ctx=canvas.getContext('2d');
  const o={W:300,H:200,dpr:1,t:0,visible:true,ctx};
  function size(){
    const w=Math.max(260,canvas.parentElement.clientWidth);
    let a=typeof opt.aspect==='function'?opt.aspect(w):(opt.aspect||0.6);
    const h=Math.round(w*a);
    o.dpr=Math.min(window.devicePixelRatio||1,2);o.W=w;o.H=h;
    canvas.style.height=h+'px';canvas.width=Math.round(w*o.dpr);canvas.height=Math.round(h*o.dpr);
  }
  size();addEventListener('resize',()=>{size();o.redraw()});
  let last=performance.now();
  o.redraw=function(){ctx.setTransform(o.dpr,0,0,o.dpr,0,0);draw(ctx,o.W,o.H,o.t,0)};
  if('IntersectionObserver' in window){new IntersectionObserver(es=>{o.visible=es[0].isIntersecting}).observe(canvas)}
  function frame(now){
    requestAnimationFrame(frame);
    let dt=Math.min(.05,(now-last)/1000);last=now;if(dt<=0)dt=1/60;
    if(!o.visible&&!opt.always)return;
    o.t+=dt;ctx.setTransform(o.dpr,0,0,o.dpr,0,0);
    try{draw(ctx,o.W,o.H,o.t,dt)}catch(err){CX.fail(err)}
  }
  requestAnimationFrame(frame);
  return o;
};
CX.fail=function(e){console.error(e);let d=document.getElementById('err');
  if(!d){d=el('div',{id:'err'});d.style.cssText='position:fixed;left:8px;bottom:8px;z-index:99;color:#ff6a5a;font:11px monospace;white-space:pre-wrap;background:#000a';document.body.appendChild(d)}
  d.textContent+=(e&&e.stack||e)+'\n'};
addEventListener('error',e=>CX.fail(e.message));

CX.pointer=function(cv,h){
  const pos=e=>{const r=cv.getBoundingClientRect();return [e.clientX-r.left,e.clientY-r.top]};
  let down=false;
  cv.addEventListener('pointerdown',e=>{down=true;cv.setPointerCapture(e.pointerId);const p=pos(e);h.down&&h.down(p[0],p[1],e)});
  cv.addEventListener('pointermove',e=>{const p=pos(e);if(down)h.drag&&h.drag(p[0],p[1],e);else h.hover&&h.hover(p[0],p[1],e)});
  const up=e=>{if(!down)return;down=false;h.up&&h.up(e)};
  cv.addEventListener('pointerup',up);cv.addEventListener('pointercancel',up);
};

/* ---- tiny 3D camera: rotate (yaw about y, pitch about x), orthographic-ish with mild perspective ---- */
CX.cam3=function(){
  const k={yaw:.7,pitch:.45,zoom:1,persp:0.12};
  k.rot=function(x,y,z){
    const cy=Math.cos(k.yaw),sy=Math.sin(k.yaw),cp=Math.cos(k.pitch),sp=Math.sin(k.pitch);
    const x1=cy*x+sy*z,z1=-sy*x+cy*z;
    const y2=cp*y-sp*z1,z2=sp*y+cp*z1;
    return [x1,y2,z2];
  };
  k.proj=function(x,y,z,cx,cy,s){
    const v=k.rot(x,y,z),w=1/(1-k.persp*v[2]*0.3);
    return [cx+v[0]*s*k.zoom*w,cy-v[1]*s*k.zoom*w,v[2],w];
  };
  return k;
};
CX.drag3=function(cv,cam,onchange){
  let lx=0,ly=0;
  CX.pointer(cv,{down:(x,y)=>{lx=x;ly=y;cv.classList.add('grab')},
    drag:(x,y)=>{cam.yaw+=(x-lx)*.008;cam.pitch=Math.max(-1.5,Math.min(1.5,cam.pitch+(y-ly)*.008));lx=x;ly=y;onchange&&onchange()},
    up:()=>cv.classList.remove('grab')});
  cv.addEventListener('wheel',e=>{if(!e.shiftKey&&!e.ctrlKey)return;e.preventDefault();cam.zoom=Math.max(.4,Math.min(3,cam.zoom*Math.exp(-e.deltaY*.0015)));onchange&&onchange()},{passive:false});
};
CX.glowLine=function(ctx,x0,y0,x1,y1,col,w,blur){ctx.strokeStyle=col;ctx.lineWidth=w||1;ctx.shadowColor=col;ctx.shadowBlur=blur==null?6:blur;ctx.beginPath();ctx.moveTo(x0,y0);ctx.lineTo(x1,y1);ctx.stroke();ctx.shadowBlur=0};
CX.dot=function(ctx,x,y,r,col,blur){ctx.fillStyle=col;ctx.shadowColor=col;ctx.shadowBlur=blur==null?10:blur;ctx.beginPath();ctx.arc(x,y,r,0,TAU);ctx.fill();ctx.shadowBlur=0};
CX.text=function(ctx,s,x,y,col,size,align){ctx.font=(size||11)+'px ui-monospace,Menlo,Consolas,monospace';ctx.fillStyle=col||C.g;ctx.textAlign=align||'left';ctx.fillText(s,x,y)};
CX.shade=function(ctx,x,y,r,col,alpha){ // lit sphere via radial gradient
  const g=ctx.createRadialGradient(x-r*.35,y-r*.4,r*.1,x,y,r);
  g.addColorStop(0,'rgba(255,255,255,'+(0.55*alpha)+')');g.addColorStop(.25,col.replace('X',alpha));g.addColorStop(1,col.replace('X',alpha*.18));
  ctx.fillStyle=g;ctx.beginPath();ctx.arc(x,y,r,0,TAU);ctx.fill();
  ctx.strokeStyle=col.replace('X',alpha*.9);ctx.lineWidth=1;ctx.stroke();
};
CX.fmt=(x,d)=>(+x).toFixed(d==null?4:d);
CX.kv=function(node,rows){node.innerHTML=rows.map(r=>'<span>'+r[0]+'</span><span'+(r[2]?' class="'+r[2]+'"':'')+'>'+r[1]+'</span>').join('')};
CX.ease=t=>t<0?0:t>1?1:t*t*(3-2*t);

function ready(){
  decor();nav();
  // footer prev/next
  const i=PAGES.findIndex(p=>p[0]===here);
  const m=document.querySelector('main');
  if(m&&i>=0&&!m.querySelector('.pn')){
    const prev=PAGES[i-1],next=PAGES[i+1];
    m.appendChild(el('div',{'class':'pn hud'},(prev?'<a href="'+prev[0]+'">&larr; '+prev[1]+'</a>':'<span></span>')+(next?'<a href="'+next[0]+'">'+next[1]+' &rarr;</a>':'<a href="lattice.html">LATTICE &rarr;</a>')));
  }
  panels();
  if(Q.get('hud')==='0')document.body.classList.add('hud-off');
  if(Q.get('y'))setTimeout(()=>scrollTo({top:+Q.get('y'),behavior:'instant'}),30);
  if(location.hash){const t=document.getElementById(location.hash.slice(1));if(t)setTimeout(()=>t.scrollIntoView({behavior:'instant'}),30)}
}
CX.ready=ready;
CX.start=function(fn){
  const go=()=>{try{ready();fn&&fn()}catch(e){CX.fail(e)}};
  if(document.readyState==='loading')document.addEventListener('DOMContentLoaded',go);else go();
};
})();
