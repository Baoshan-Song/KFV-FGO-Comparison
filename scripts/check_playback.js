/* Browser-independent checks of the actual replay script. Node or JavaScriptCore. */
const readText=typeof require==='function' ? p=>require('fs').readFileSync(p,'utf8') : readFile;
const log=typeof print==='function' ? print : console.log;
const template=readText('kfv_fgo/comparison.html');
function assert(ok,message){if(!ok)throw Error(message);}
function element(value='') {
  return {value,style:{},dataset:{},listeners:{},classList:{toggle(){}},textContent:'',innerHTML:'',
    addEventListener(type,fn){this.listeners[type]=fn;},setAttribute(k,v){this[k]=v;},remove(){this.removed=true;}};
}
function check(real) {
  const n=24,time=Array.from({length:n},(_,i)=>i===n-1?30:i);
  const truth=[time.map(x=>Math.cos(x/10)*10),time.map(x=>Math.sin(x/10)*10)];
  if(real)truth.push(time.map(()=>0));
  const trajectories=[truth.map(row=>row.map(v=>v+.01)),truth.map(row=>row.map(v=>v-.02))];
  const residuals=Object.fromEntries(['postfit','prefit','truth'].map(k=>[k,[[time.map(()=>.1)],[time.map(()=>.2)]]]));
  const payload={real,frames:n,names:['iEKF','Re-FGO'],kinds:['kfv','fgo'],dt:1,time,truth,
    trajectories,residuals,metrics:[{rmse:.02,up_rmse:.01},{rmse:.03,up_rmse:.02}],has_truth:true,
    errors:[0,1].map(()=>({horizontal:time.map(()=>.1),up:time.map(()=>-.2),'3d':time.map(()=>.3)}))};
  const controls=new Map(),tabs=['trajectory','error','residual','iterations','difference','satellites','runtime'].map(view=>{
    const el=element();el.dataset.view=view;controls.set(`[data-view=${view}]`,el);return el;
  });
  const get=s=>{if(!controls.has(s))controls.set(s,element(s==='.span'?'8':s==='.epoch'||s==='.frame'?'15':s==='.speed'?'1':s==='.residual'?'postfit':s==='.error-component'?'horizontal':'linear'));return controls.get(s);};
  const root={isConnected:true,querySelector:get,querySelectorAll:s=>s==='[data-view]'?tabs:[]};
  const doc={hidden:false,getElementById:()=>root,listeners:{},addEventListener(type,fn){this.listeners[type]=fn;}};
  let clock=0,nextId=0;const frames=new Map();
  const request=fn=>{frames.set(++nextId,fn);return nextId;},cancel=id=>frames.delete(id);
  function tick(ms){clock+=ms;const pending=[...frames.values()];frames.clear();pending.forEach(fn=>fn(clock));}
  const source=template.split('<script>')[1].split('</script>')[0].replaceAll('__ID__','test').replace('__PAYLOAD__',JSON.stringify(payload)).replaceAll('__INTERVAL__','80');
  new Function('document','requestAnimationFrame','cancelAnimationFrame','performance',source)(doc,request,cancel,{now:()=>clock});
  assert(get('.figure').innerHTML.includes('Complete trajectory'),'original full trajectory');
  assert(get('.figure').innerHTML.includes('stroke-dasharray="7 5"'),'method styles preserved');
  assert(!get('.figure').innerHTML.includes('NaN'),'finite SVG');
  get('.play').listeners.click();tick(80);
  assert(Number(get('.epoch').value)===16,'play advances');
  assert(get('.play').textContent.includes('Pause'),'pause button');
  get('.play').listeners.click();tick(500);
  assert(Number(get('.epoch').value)===16&&frames.size===0,'pause stops');
  get('.play').listeners.click();get('.frame').value=3;get('.frame').listeners.input();tick(16);
  assert(Number(get('.epoch').value)===3&&get('.play').textContent.includes('Play'),'scrub pauses');
  get('.speed').value='2';get('.play').listeners.click();tick(40);
  assert(Number(get('.epoch').value)===4,'speed doubles');
  get('.play').listeners.click();get('.epoch').value=n-1;get('.epoch').listeners.input();tick(16);
  get('.play').listeners.click();assert(Number(get('.epoch').value)===1,'restart at end');
  get('.play').listeners.click();get('.epoch').value=n-2;get('.epoch').listeners.input();tick(16);
  get('.play').listeners.click();tick(40);
  assert(Number(get('.epoch').value)===n-1&&frames.size===0,'auto pause at last sample');
  assert(get('.play-time').textContent==='30 s / 30 s','recorded time displayed');
  tabs.find(x=>x.dataset.view==='error').listeners.click();
  if(real){get('.error-component').value='up';get('.error-component').listeners.input();tick(16);}
  assert(get('.figure').innerHTML.includes(real?'Up error':'2D position error'),'error view');
  tabs.find(x=>x.dataset.view==='residual').listeners.click();
  assert(get('.figure').innerHTML.includes('Residual norm'),'residual view');
  assert(!get('.caption').textContent.startsWith('Cannot open'),'no renderer error');
  log('PASS '+(real?'real':'simulation')+' original styles, Play/Pause, speed, seek, endpoint, errors and residuals');
}
check(false);check(true);
