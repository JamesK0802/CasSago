/** Cancel unobservable transportation cycles while preserving EVERY
 * reference/window/length marginal exactly. All candidate paths stay in the
 * dictionary and can re-enter. This only chooses a sparse numerical basis
 * for one or two windows, the currently supported analysis contract. */
export function sparseLengthFlowBasis(paths:{reference:number;sizes:number[]}[],input:number[]){
 if(paths.some(p=>p.sizes.length>2))return input.slice();
 const weights=input.slice(),nodes=new Map<string,number>(),parent:number[]=[],adjacency:Map<number,number>[]=[];
 const node=(key:string)=>{let id=nodes.get(key);if(id===undefined){id=nodes.size;nodes.set(key,id);parent.push(id);adjacency.push(new Map());}return id;},find=(id:number):number=>parent[id]===id?id:parent[id]=find(parent[id]);
 const endpoints=paths.map(p=>[node(p.reference+':'+(p.sizes.length===1?'root':'first:'+p.sizes[0])),node(p.reference+':last:'+p.sizes.at(-1))]),order=input.flatMap((w,i)=>w>0?[i]:[]).sort((a,b)=>input[b]-input[a]),add=(i:number)=>{const[a,b]=endpoints[i];adjacency[a].set(i,b);adjacency[b].set(i,a);};
 for(const i of order){const[a,b]=endpoints[i],ra=find(a),rb=find(b);if(ra!==rb){parent[ra]=rb;add(i);continue;}
  const previous=new Map<number,{node:number;edge:number}>([[a,{node:a,edge:-1}]]),queue=[a];for(const v of queue){if(v===b)break;for(const[e,u]of adjacency[v])if(!previous.has(u)){previous.set(u,{node:v,edge:e});queue.push(u);}}
  const cycle:number[]=[];for(let v=b;v!==a;){const p=previous.get(v)!;cycle.push(p.edge);v=p.node;}cycle.reverse();let step=Infinity,blocker=-1;for(let j=0;j<cycle.length;j+=2)if(weights[cycle[j]]<step){step=weights[cycle[j]];blocker=cycle[j];}
  weights[i]+=step;for(let j=0;j<cycle.length;j++)weights[cycle[j]]+=j%2?step:-step;
  const[x,y]=endpoints[blocker];adjacency[x].delete(blocker);adjacency[y].delete(blocker);add(i);
 }
 return weights;
}
