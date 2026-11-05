import {Component,HostListener,ViewEncapsulation} from '@angular/core';
import {CommonModule} from '@angular/common';
import {SangerDesignPageComponent} from './pages/sanger-design-page/sanger-design-page.component';
import {GrnaDesignPageComponent} from './pages/grna-design-page/grna-design-page.component';
import {ContactPageComponent} from './pages/contact-page/contact-page.component';
@Component({selector:'app-root',standalone:true,imports:[CommonModule,SangerDesignPageComponent,GrnaDesignPageComponent,ContactPageComponent],encapsulation:ViewEncapsulation.None,styleUrls:['./app.css'],template:`
 <div class="cassago-shell"><nav class="cassago-nav" aria-label="Main navigation">
  <a class="cassago-brand" href="#/sanger-design" aria-label="CasSago home">Cas<span>Sago</span></a>
  <div class="cassago-links"><button type="button" [class.active]="page==='workspace'" (click)="showPage('workspace')">Sanger tools</button><button type="button" [class.active]="page==='guides'" (click)="showPage('guides')">gRNA Design</button><button type="button" [class.active]="page==='contact'" (click)="showPage('contact')">Contact</button></div>
 </nav><main class="cassago-content">
  <section class="cassago-gate" *ngIf="page!=='contact'&&!developmentUnlocked" aria-labelledby="development-title"><h1 id="development-title">CasSago</h1><p>This feature is under development. Enter the development password to continue.</p><form (submit)="$event.preventDefault();unlockDevelopment(password)"><label for="development-password">Password</label><input #password id="development-password" type="password" inputmode="numeric" autocomplete="off" required><button type="submit">Enter</button><p role="alert" *ngIf="developmentError">{{developmentError}}</p></form></section>
  <app-sanger-design-page *ngIf="page==='workspace'&&developmentUnlocked" />
  <app-grna-design-page *ngIf="page==='guides'&&developmentUnlocked" (primerRequested)="showPage('workspace')" />
  <app-contact-page *ngIf="page==='contact'" (goHome)="showPage('workspace')" />
 </main></div>`,styles:[`
 .cassago-shell{min-height:100dvh;background:var(--color-bg,#f7f8fa);color:#202b3d}.cassago-nav{display:flex;align-items:center;gap:36px;min-height:64px;padding:10px clamp(20px,4vw,64px);background:#fff;border-bottom:1px solid #e5e9ed}.cassago-brand{font-size:23px;font-weight:700;letter-spacing:-.7px;text-decoration:none;color:#273646}.cassago-brand span{color:#168264}.cassago-links{display:flex;align-items:center;gap:8px;flex-wrap:wrap}.cassago-links button{border:0;background:transparent;border-radius:7px;padding:10px 14px;color:#65758a;font:inherit;font-size:13px;font-weight:600;cursor:pointer}.cassago-links button.active{background:#eaf6f0;color:#16724f}.cassago-links button:hover{background:#f1f5f3}.cassago-content{padding:12px clamp(12px,2vw,32px) 32px}.cassago-gate{max-width:460px;margin:64px auto;padding:30px;border:1px solid #e2e8f0;border-radius:12px;background:#fff}.cassago-gate h1{font-size:26px;margin-bottom:14px}.cassago-gate p{line-height:1.6;color:#64748b;font-size:14px}.cassago-gate form{display:grid;gap:10px;margin-top:24px}.cassago-gate input,.cassago-gate button{font:inherit;padding:11px;border:1px solid #cbd5e1;border-radius:7px}.cassago-gate button{background:#168264;color:#fff;cursor:pointer;border:0}.cassago-gate [role=alert]{color:#b91c1c}.cassago-nav button:focus-visible,.cassago-brand:focus-visible{outline:2px solid #d99a30;outline-offset:3px}@media(max-width:600px){.cassago-nav{gap:12px;align-items:flex-start;flex-direction:column;padding:14px 18px}.cassago-links{gap:2px}.cassago-links button{padding:9px 10px}.cassago-content{padding:10px 10px 24px}}
 `]})
export class App{
 page:'workspace'|'guides'|'contact'='workspace';developmentUnlocked=false;developmentError='';
 constructor(){this.readRoute();}
 @HostListener('window:hashchange') readRoute(){const route=window.location.hash.replace(/^#/,'');this.page=route==='/grna-design'?'guides':route==='/contact'?'contact':'workspace';}
 showPage(page:'workspace'|'guides'|'contact'){this.page=page;window.location.hash=page==='guides'?'/grna-design':page==='contact'?'/contact':'/sanger-design';}
 unlockDevelopment(input:HTMLInputElement){if(input.value==='001026'){this.developmentUnlocked=true;this.developmentError='';input.value='';}else{this.developmentError='Incorrect password.';input.value='';input.focus();}}
}
