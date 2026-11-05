import {Injectable} from '@angular/core';
import {SangerReference,SangerDesignSettings} from '../models/sanger-design';
export interface PrimerHandoff {references:SangerReference[];workflow:SangerDesignSettings['workflow'];guide:string}
@Injectable({providedIn:'root'})
export class DesignHandoffService {private pending:PrimerHandoff|null=null;stage(input:PrimerHandoff){this.pending=structuredClone(input);}take(){const pending=this.pending;this.pending=null;return pending;}}
