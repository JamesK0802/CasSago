import {SangerAnalysisResult} from '../models/sanger-analysis';
import {sangerCurationView} from './sanger-curation';
interface PdfOptions {curated:boolean;excludedGroups:string[];bounds:{left:number;right:number};threshold:number;summary:{name:string;target:string;share:number;wt:number|null;inFrame:number|null;outOfFrame:number|null}[]}
/** Portable vector report; PDF output follows the active original/curated selection. */
export async function downloadSangerPdf(result:SangerAnalysisResult,options:PdfOptions){
 const [{jsPDF},{default:autoTable}]=await Promise.all([import('jspdf'),import('jspdf-autotable')]);
 const doc=new jsPDF({orientation:'landscape',unit:'mm',format:'a4'}),width=doc.internal.pageSize.getWidth();
 const pct=(v:number|null)=>v===null?'Not evaluated':(100*v).toFixed(2)+'%';
 doc.setFontSize(21);doc.setTextColor(90,92,119);doc.text('CasSago',14,18);doc.setFontSize(12);doc.text(options.curated?'Curated Sanger analysis':'Sanger analysis',14,27);doc.setTextColor(45,55,72);doc.setFontSize(9);
 doc.text(doc.splitTextToSize('Sample: '+result.sampleName+'   |   Control: '+(result.controlName||'not recorded'),width-28),14,35);
 doc.text('Window: '+options.bounds.left+' bp left / '+options.bounds.right+' bp right   |   '+new Date().toISOString(),14,43);
 if(options.curated)doc.text('Selected group fractions are normalized within each reference and target; original fit retained in the JSON report.',14,50);
 autoTable(doc,{startY:options.curated?55:48,head:[['Reference','Target','Signal share','WT','In-frame','Out-of-frame']],body:options.summary.map(row=>[row.name,row.target,pct(row.share),pct(row.wt),pct(row.inFrame),pct(row.outOfFrame)]),styles:{fontSize:8,cellPadding:2.5},headStyles:{fillColor:[124,129,158]}});
 for(const [index,a] of (result.annotations||[]).entries())for(let target=0;target<(a.cuts?.length||1);target++){
  const view=sangerCurationView(result,index,target,options.excludedGroups,options.bounds);if(!view)continue;
  const rows=options.curated?view.rows:view.view.rows,displayRows=rows.filter(row=>!row.edited||row.fraction*100>=options.threshold);
  doc.addPage();doc.setFontSize(13);doc.setTextColor(90,92,119);doc.text(a.name+' - T'+(target+1),14,18);doc.setTextColor(45,55,72);doc.setFontSize(9);
  doc.text('Reference share: '+pct(a.share)+'   |   '+a.direction+'   |   '+(a.estimable?'':'QC review'),14,25);
  autoTable(doc,{startY:31,head:[['Group','Within-reference fraction','Type','Net indel','Sequence / position state']],body:displayRows.map((row,i)=>['Group '+(i+1),pct(row.fraction),row.classification,String(row.size)+' bp',row.sequenceState+'; '+row.positionLabel]),styles:{fontSize:8,cellPadding:2.3},headStyles:{fillColor:[124,129,158]}});
  // Wrap complete visible sequence windows into fixed-length monospaced lines.
  const sequenceRows:any[][]=[['REFERENCE',view.view.reference]];
  displayRows.forEach((row,i)=>sequenceRows.push(['GROUP '+(i+1),row.tokens.map(t=>t.type==='delete'?'-'.repeat(t.val.length):t.type==='insert'?'['+t.val+']':t.val).join('')]));
  const startY=(doc as any).lastAutoTable.finalY+8;
  autoTable(doc,{startY,head:[['Group','Displayed sequence window (insertions in brackets; deletions as dashes)']],body:sequenceRows.map(([label,sequence])=>[label,sequence.match(/.{1,115}/g)?.join('\n')||'Unresolved']),styles:{font:'courier',fontSize:7,cellPadding:2,overflow:'linebreak'},columnStyles:{0:{cellWidth:29},1:{cellWidth:width-57}},headStyles:{font:'helvetica',fillColor:[124,129,158]},rowPageBreak:'avoid'});
 }
 const pages=doc.getNumberOfPages();for(let page=1;page<=pages;page++){doc.setPage(page);doc.setFont('helvetica');doc.setFontSize(8);doc.setTextColor(125);doc.text('CasSago | '+(options.curated?'Curated':'Original')+' | '+page+' / '+pages,14,doc.internal.pageSize.getHeight()-8);}
 doc.save(options.curated?'CasSago_curated.pdf':'CasSago_result.pdf');
}
