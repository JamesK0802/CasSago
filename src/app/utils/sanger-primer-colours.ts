import {CommonPartner} from '../models/sanger-design';

const OPPOSITE_COLOURS=['#d39423','#8b5fc0','#347cb6','#d46e67','#667b36','#b65b91','#527f8c','#986548'];
/** Identical oligos share a colour; distinct oligos within a set stay distinct. */
export function oppositePrimerColour(partner:CommonPartner|null|undefined,index:number):string {
 if(!partner)return OPPOSITE_COLOURS[0];
 const sequences=partner.perReference?.map(p=>p.sequence)||[partner.sequence];
 const unique=[...new Set(sequences)];
 return OPPOSITE_COLOURS[Math.max(0,unique.indexOf(sequences[index]||sequences[0]))%OPPOSITE_COLOURS.length];
}
