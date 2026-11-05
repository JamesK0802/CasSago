import {Component,HostListener,ViewEncapsulation,ViewChildren,QueryList} from '@angular/core';
import {CommonModule} from '@angular/common';
import {SangerDesignPageComponent} from './pages/sanger-design-page/sanger-design-page.component';
import {SangerAnalysisPanelComponent} from './pages/sanger-analysis-panel/sanger-analysis-panel.component';
import {LandingPageComponent} from './pages/landing-page/landing-page.component';
import {GrnaDesignPageComponent} from './pages/grna-design-page/grna-design-page.component';
import {ContactPageComponent} from './pages/contact-page/contact-page.component';
type Page='home'|'analysis'|'viewer'|'workspace'|'guides'|'contact';
@Component({selector:'app-root',standalone:true,imports:[CommonModule,SangerDesignPageComponent,SangerAnalysisPanelComponent,LandingPageComponent,GrnaDesignPageComponent,ContactPageComponent],encapsulation:ViewEncapsulation.None,styleUrls:['./app.css'],template:`
<div class="app-shell"><nav class="top-nav" aria-label="Main navigation"><div class="nav-left"><button class="nav-brand cassago-brand" aria-label="CasSago home" (click)="showPage('home')" [disabled]="locked"><img class="nav-brand-logo" src="cassago-logo.png" alt="CasSago"></button><div class="cassago-links"><button class="nav-tab btn-tab" *ngFor="let tab of navigation" [class.active]="page===tab.id" [disabled]="locked" (click)="showPage(tab.id)">{{tab.label}}</button></div></div><div class="nav-right cassago-nav-right"><button type="button" class="contact-nav-link" [class.active]="page==='contact'" [disabled]="locked" (click)="showPage('contact')">Contact</button><button type="button" class="btn-guide-trigger" [disabled]="locked" title="Guide coming soon" aria-label="Guide (coming soon)"><span class="guide-btn-text">Guide</span></button></div></nav>
<main class="app-content"><app-landing-page *ngIf="page==='home'" (openTool)="showPage($event)" />
<div class="container crispr-body cassago-workspace" [hidden]="page!=='analysis'">
<div class="analysis-tabs-bar"><div class="tabs-scroll-container"><div class="analysis-tab-item" *ngFor="let session of sessions" [class.active]="session===selected" [class.locked]="locked"><button type="button" [disabled]="locked" (click)="selected=session">Analysis {{session}}</button><button type="button" class="tab-close-btn" *ngIf="sessions.length>1" [disabled]="locked" (click)="close(session)" aria-label="Close analysis">×</button></div></div><button class="btn-add-tab" [disabled]="locked" (click)="add()" aria-label="New analysis">+</button></div>
<app-sanger-analysis-panel *ngFor="let session of sessions" [hidden]="session!==selected" [sessionId]="session"  />
</div><div class="container crispr-body cassago-workspace" [hidden]="page!=='viewer'"><app-sanger-analysis-panel [resultsOnly]="true" [sessionId]="0" /></div>
<app-sanger-design-page [hidden]="page!=='workspace'" (templateRequested)="applyTemplate($event)" (syntheticRequested)="applySynthetic($event)" />
<section class="cassago-construction" *ngIf="page==='guides'&&!grnaUnlocked"><img src="cassago-construction.png" alt="CasSago wearing a construction hat behind a barricade"><h1>Under construction</h1><p>gRNA Design</p><form (submit)="$event.preventDefault();unlockGrna(grnaPassword)"><label for="grna-password">Development password</label><div><input id="grna-password" #grnaPassword type="password" autocomplete="off" required><button type="submit">Enter</button></div><p class="gate-error" role="alert" *ngIf="grnaError">{{grnaError}}</p></form></section>
<app-grna-design-page [hidden]="page!=='guides'||!grnaUnlocked" (primerRequested)="showPage('workspace')" />
<app-contact-page *ngIf="page==='contact'" (goHome)="showPage('home')" /></main></div>`})
export class App{
 @ViewChildren(SangerAnalysisPanelComponent) panels?:QueryList<SangerAnalysisPanelComponent>;
 @ViewChildren(SangerDesignPageComponent) designers?:QueryList<SangerDesignPageComponent>;
 @ViewChildren(GrnaDesignPageComponent) guides?:QueryList<GrnaDesignPageComponent>;
 grnaUnlocked=false;grnaError='';
 unlockGrna(input:HTMLInputElement){if(input.value==='001026'){this.grnaUnlocked=true;this.grnaError='';}else{this.grnaError='Incorrect password.';input.focus();}input.value='';}
 page:Page='home';sessions=[1];selected=1;private next=2;
 readonly navigation:{id:Page;label:string}[]=[{id:'analysis',label:'Sanger Analysis'},{id:'viewer',label:'Results Viewer'},{id:'workspace',label:'Primer Design'},{id:'guides',label:'gRNA Design'}];
 get locked(){return !!(this.panels?.some(p=>p.busy||p.loading)||this.designers?.some(p=>p.busy)||this.guides?.some(p=>p.busy));}
 constructor(){this.readRoute();}
 @HostListener('window:hashchange') readRoute(){const route=window.location.hash.replace(/^#?\/?/,'');const page:Page=route==='sanger-analysis'?'analysis':route==='results'?'viewer':route==='sanger-design'?'workspace':route==='grna-design'?'guides':route==='contact'?'contact':'home';if(this.locked){window.location.hash=this.route(this.page);return;}this.page=page;}
 private route(page:Page){return {home:'/',analysis:'/sanger-analysis',viewer:'/results',workspace:'/sanger-design',guides:'/grna-design',contact:'/contact'}[page];}
 showPage(page:Page){if(this.locked)return;this.page=page;window.location.hash=this.route(page);}
 add(){if(this.locked)return;this.selected=this.next++;this.sessions=[...this.sessions,this.selected];this.showPage('analysis');}
 close(id:number){if(this.locked||this.sessions.length===1)return;this.sessions=this.sessions.filter(s=>s!==id);if(this.selected===id)this.selected=this.sessions[0];}
 private get active(){return this.panels?.toArray()[this.sessions.indexOf(this.selected)];}
 applyTemplate(value:unknown){this.active?.applyTemplate(value);this.showPage('analysis');}
 applySynthetic(data:any){this.active?.applySynthetic(data.template,data.control,data.edit,{recipe:data.recipe,truth:data.truth,mixing:data.mixing,cycles:data.cycles});this.showPage('analysis');}
}
