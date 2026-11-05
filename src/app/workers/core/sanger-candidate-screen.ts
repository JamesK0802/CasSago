/** Cache exact fluorescence differences inside a physical length/target family.
 * It changes candidate screening cost, not the downstream mixture objective. */
export class CandidateSignalScreen {
 private unitAnchors=new Map<number,number[]>();private unitSparse:{at:number[];v:number[];norm:number;anchor:number}[]=[];
 private anchors=new Map<string,number>();private sparse:{at:number[];v:number[];norm:number;anchor:number;family:string}[]=[];private rhs:number[]=[];
 constructor(private columns:number[][]){}
 update(families:string[]){for(let id=this.sparse.length;id<this.columns.length;id++){const family=families[id];if(!this.anchors.has(family))this.anchors.set(family,id);const anchor=this.anchors.get(family)!,at:number[]=[],v:number[]=[];let norm=0;for(let k=0;k<this.columns[id].length;k++){const d=this.columns[id][k]-this.columns[anchor][k];if(Math.abs(d)>1e-12){at.push(k);v.push(d);norm+=d*d;}}this.sparse.push({at,v,norm,anchor,family});}}
 correlate(residual:number[]){const base=new Map<number,number>();for(const id of this.anchors.values())base.set(id,this.columns[id].reduce((s,v,k)=>s+v*residual[k],0));this.rhs=this.sparse.map(c=>base.get(c.anchor)!+c.at.reduce((s,k,j)=>s+c.v[j]*residual[k],0));return this.rhs;}

 private updateUnit(){for(let id=this.unitSparse.length;id<this.columns.length;id++){const anchor=this.sparse[id].anchor;let unitAnchor=this.unitAnchors.get(anchor);if(!unitAnchor){const col=this.columns[anchor];unitAnchor=col.map((value,k)=>{const at=k-k%4,total=col[at]+col[at+1]+col[at+2]+col[at+3];return total?value/total:0;});this.unitAnchors.set(anchor,unitAnchor);}const ua:number[]=[],uv:number[]=[];let unorm=0;for(let k=0;k<this.columns[id].length;k+=4){const col=this.columns[id],total=col[k]+col[k+1]+col[k+2]+col[k+3];for(let b=0;b<4;b++){const value=total?col[k+b]/total:0,d=value-unitAnchor[k+b];unorm+=value*value;if(Math.abs(d)>1e-12){ua.push(k+b);uv.push(d);}}}this.unitSparse.push({at:ua,v:uv,norm:unorm,anchor});}}
 /** Screening score for adding a candidate to the current mixture. Unit
  * waveforms make the derivative scale comparable across candidate lengths.
  * Full measured-dye fitting still decides acceptance. */
 additionScores(observed:number[],prediction:number[]){this.updateUnit();
  const residual=observed.map((v,k)=>v-prediction[k]),dot=(vector:number[])=>{const base=new Map<number,number>();for(const [id,col]of this.unitAnchors)base.set(id,col.reduce((sum,v,k)=>sum+v*vector[k],0));return this.unitSparse.map(c=>base.get(c.anchor)!+c.at.reduce((sum,k,j)=>sum+c.v[j]*vector[k],0));},rhs=dot(residual),cross=dot(prediction),mean=prediction.reduce((sum,v,k)=>sum+v*residual[k],0),norm=prediction.reduce((sum,v)=>sum+v*v,0);
  return this.unitSparse.map((c,id)=>{const gradient=rhs[id]-mean,curvature=Math.max(0,c.norm-2*cross[id]+norm);if(gradient<=0||curvature<=1e-14)return 0;const step=Math.min(1,gradient/curvature);return 2*step*gradient-step*step*curvature;});
 }
 sameFamily(a:number,b:number){return this.sparse[a].family===this.sparse[b].family;}
 replacementGain(id:number,old:number,weight:number){const a=this.sparse[id],b=this.sparse[old];if(a.anchor!==b.anchor)throw Error('Candidate transfer requires the same physical family.');let i=0,j=0,cross=0;while(i<a.at.length&&j<b.at.length){if(a.at[i]===b.at[j])cross+=a.v[i++]*b.v[j++];else if(a.at[i]<b.at[j])i++;else j++;}return 2*weight*(this.rhs[id]-this.rhs[old])-weight*weight*(a.norm+b.norm-2*cross);}
}
