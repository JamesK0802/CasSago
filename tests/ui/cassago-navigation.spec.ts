import '@angular/compiler';
import{it,expect,vi,afterEach}from'vitest';
import{App}from'../../src/app/app';
afterEach(()=>vi.unstubAllGlobals());
it('routes the public tools and landing page and defaults to the landing page',()=>{const location={hash:''};vi.stubGlobal('window',{location});const app=new App();expect(app.page).toBe('home');expect(app.navigation).toHaveLength(4);for(const tab of app.navigation){app.showPage(tab.id);app.readRoute();expect(app.page).toBe(tab.id);}});
it('keeps sessions independent and blocks switching while a worker is busy',()=>{vi.stubGlobal('window',{location:{hash:''}});const app=new App();app.add();expect(app.sessions).toEqual([1,2]);app.panels={some:(fn:any)=>fn({busy:true}),toArray:()=>[]}as any;app.showPage('contact');app.add();app.close(2);expect(app.page).toBe('analysis');expect(app.sessions).toEqual([1,2]);});
