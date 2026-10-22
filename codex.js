/* codex.js: shared shell for the geometry site. Offline, no dependencies.
   Provides the index panel, minimisable panels, H to hide the HUD, scanlines/grain/backdrop,
   and small helpers for canvas figures (CX.canvas, CX.pointer, CX.cam3). */
(function(){
"use strict";
const PAGES=[
 ['index.html','ABSTRACT','overview'],
 ['lattice.html','LATTICE','the visual'],
 ['origin.html','THE DRAWING','origin'],
 ['mirrors.html','MIRRORS','why'],
 ['matter.html','MATTER','crystals'],
 ['motion.html','MOTION','the fold'],
 ['mind.html','PERCEPTION','mind'],
 ['drift.html','DRIFT','wandering mind'],
 ['seek.html','SEARCH','drift + check'],
 ['open.html','LEDGER','claims'],
 ['research.html','RESEARCH','updates']
];
const EXTRA=[['cube-six.html','The cube, hands-on'],['deeper.html','Atoms, tilings and dimensions']];
const Q=new URLSearchParams(location.search);
const here=(location.pathname.split('/').pop()||'index.html');
const C={g:'#37ff9a',c:'#35e6ff',a:'#ffb23a',r:'#ff4a3d',dim:'#68aa8e',bg:'#02090c',ink:'#eafff4'};
const TAU=Math.PI*2;
const CX=window.CX={C,TAU,Q};
CX.reduce=matchMedia('(prefers-reduced-motion: reduce)').matches; // figures start still and do not pulse

function el(tag,attrs,html){const e=document.createElement(tag);if(attrs)for(const k in attrs)e.setAttribute(k,attrs[k]);if(html!=null)e.innerHTML=html;return e}

/* ---- decorative layers ---- */
function decor(){
  const R=Math.SQRT2-1;
  const s=el('div',{id:'sig','aria-hidden':'true'});
  const pt=(r,i)=>[(r*Math.cos(i*TAU/8-Math.PI/8)).toFixed(1),(r*Math.sin(i*TAU/8-Math.PI/8)).toFixed(1)];
  let h='<svg xmlns="http://www.w3.org/2000/svg" viewBox="-500 -500 1000 1000" width="100%" height="100%"><g class="rot" fill="none" stroke="#37ff9a" stroke-width="1">';
  [460,460*R,460*R*R,460*R*R*R].forEach((r,i)=>{h+='<circle r="'+r.toFixed(1)+'" opacity="'+[.35,.3,.25,.2][i]+'"/>'});
  let o='',s8='';for(let i=0;i<8;i++){o+=(i?'L':'M')+pt(460,i).join(' ');s8+=(i?'L':'M')+pt(460,(i*3)%8).join(' ')}
  h+='<path d="'+o+'Z" opacity=".3"/><path d="'+s8+'Z" stroke="#35e6ff" opacity=".18"/></g></svg>';
  s.innerHTML=h;document.body.appendChild(s);
  document.body.appendChild(el('div',{id:'scan','aria-hidden':'true'}));
  const g=el('div',{id:'grain','aria-hidden':'true'}),cc=document.createElement('canvas');cc.width=cc.height=128;
  const c=cc.getContext('2d'),im=c.createImageData(128,128);
  for(let i=0;i<im.data.length;i+=4){const v=Math.random()*255|0;im.data[i]=v*.4;im.data[i+1]=v;im.data[i+2]=v*.7;im.data[i+3]=255}
  c.putImageData(im,0,0);g.style.backgroundImage='url('+cc.toDataURL()+')';document.body.appendChild(g);
}

/* ---- index panel ---- */
function nav(){
  const p=el('nav',{'class':'panel hud',id:'p-nav','aria-label':'Site index'});
  let b='';
  for(const [u,n,d] of PAGES)b+='<a href="'+u+'"'+(u===here?' class="cur" aria-current="page"':'')+'>'+(u===here?'&gt; ':'')+n+'<span>'+d+'</span></a>';
  b+='<hr><div class="foot">CLASSIC // plain light pages</div>';
  for(const [u,n] of EXTRA)b+='<a href="'+u+'"'+(u===here?' class="cur" aria-current="page"':'')+'>'+n+'<span>classic</span></a>';
  b+='<div class="foot">[H] hide HUD</div>';
  p.innerHTML='<button class="burger" aria-expanded="false" aria-controls="p-nav-body">&equiv; INDEX<span>R=&radic;2&minus;1</span></button><header><span class="tabtitle">INDEX // R=&radic;2&minus;1</span><button class="mb" title="minimise">&ndash;</button></header><div class="body" id="p-nav-body">'+b+'</div>';
  document.body.insertBefore(p,document.body.firstChild);
  const bg=p.querySelector('.burger');
  bg.onclick=()=>{p.classList.toggle('open');p._syncPanel&&p._syncPanel()};
  if(!matchMedia('(max-width:700px)').matches&&((Q.get('min')||'').split(',').includes('index')||Q.get('min')==='all'))p.classList.add('min');
}

function panels(){
  const mins=(Q.get('min')||'').split(',');
  for(const [i,p] of [...document.querySelectorAll('.panel')].entries()){
    const btn=p.querySelector('.mb');if(!btn||btn._w)continue;btn._w=1;
    const body=p.querySelector('.body'),name=p.querySelector('.tabtitle')?.textContent.trim()||'controls';
    if(body&&!body.id)body.id=(p.id||'panel-'+i)+'-body';
    btn.type='button';if(body)btn.setAttribute('aria-controls',body.id);
    const sync=()=>{
      const mobile=p.id==='p-nav'&&matchMedia('(max-width:700px)').matches;
      const collapsed=mobile?!p.classList.contains('open'):p.classList.contains('min');
      btn.innerHTML=collapsed?'+':'&ndash;';
      btn.setAttribute('aria-expanded',String(!collapsed));
      btn.setAttribute('aria-label',(collapsed?'Expand ':'Collapse ')+name+' panel');
      btn.title=(collapsed?'Expand ':'Collapse ')+name+' panel';
      if(body)body.hidden=collapsed;
      const burger=p.querySelector('.burger');
      if(burger){burger.setAttribute('aria-expanded',String(!collapsed));burger.setAttribute('aria-label',(collapsed?'Open':'Close')+' site index')}
    };
    if(Q.get('min')==='all'||(p.id&&mins.includes(p.id.replace(/^p-/,''))))p.classList.add('min');
    p._syncPanel=sync;sync();btn.onclick=()=>{p.classList.toggle('min');sync()};
    if(p.id==='p-nav')addEventListener('resize',sync);
  }
  for(const [i,input] of [...document.querySelectorAll('input[type="range"]')].entries()){
    if(input.labels?.length||input.hasAttribute('aria-label')||input.hasAttribute('aria-labelledby'))continue;
    let label=input.closest('.sl')?.previousElementSibling;
    if(!label?.classList.contains('lab'))label=input.closest('.grp')?.querySelector('.lab');
    if(label){if(!label.id)label.id=(input.id||'range-'+i)+'-label';input.setAttribute('aria-labelledby',label.id)}
  }
}
CX.panels=panels;

addEventListener('keydown',e=>{
  if(e.ctrlKey||e.metaKey||e.altKey)return;
  const target=e.target;if(target?.isContentEditable||target?.closest('input,select,textarea'))return;
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
  let last=performance.now(),raf=0;
  const active=()=>!document.hidden&&o.visible&&canvas.isConnected!==false;
  // a draw that keeps throwing is reported once, not on every frame
  function paint(dt){ctx.setTransform(o.dpr,0,0,o.dpr,0,0);try{draw(ctx,o.W,o.H,o.t,dt);o.err=null}catch(err){const m=String(err&&err.stack||err);if(m!==o.err){o.err=m;CX.fail(err)}}}
  o.redraw=function(){if(active())paint(0)};
  function schedule(){
    if(!active()){if(raf)cancelAnimationFrame(raf);raf=0;return}
    if(opt.animate===false){o.redraw();return}
    if(!raf){last=performance.now();raf=requestAnimationFrame(frame)}
  }
  function frame(now){
    raf=0;if(!active())return;
    let dt=Math.min(.05,(now-last)/1000);last=now;if(dt<=0)dt=1/60;
    o.t+=dt;paint(dt);raf=requestAnimationFrame(frame);
  }
  function viewport(){const r=canvas.getBoundingClientRect();o.visible=r.bottom>0&&r.right>0&&r.top<innerHeight&&r.left<innerWidth}
  size();viewport();
  addEventListener('resize',()=>{size();viewport();o.redraw();schedule()});
  document.addEventListener('visibilitychange',schedule);
  if('IntersectionObserver' in window){new IntersectionObserver(es=>{o.visible=es[0].isIntersecting;schedule()}).observe(canvas)}
  else addEventListener('scroll',()=>{viewport();schedule()},{passive:true});
  schedule();
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

/* ---- signature badge (stays visible when the HUD is hidden) ---- */
function sign(){
  if(document.getElementById('sign'))return;
  const L=(u,t)=>'<a href="'+u+'" target="_blank" rel="noopener">'+t+'</a>';
  document.body.appendChild(el('div',{id:'sign'},'Ahraz Arifuddin from '+L('https://digital.askaconsult.com','Aska Digital')+'<br>'+L('https://www.anthropic.com/claude','Claude Opus 5.5')+' from Anthropic'));
}

function ready(){
  decor();nav();sign();
  // footer prev/next
  const i=PAGES.findIndex(p=>p[0]===here);
  const m=document.querySelector('main');
  if(m){if(!m.id)m.id='main-content';m.tabIndex=-1;document.body.insertBefore(el('a',{'class':'skip-link',href:'#'+m.id},'Skip to main content'),document.body.firstChild)}
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
