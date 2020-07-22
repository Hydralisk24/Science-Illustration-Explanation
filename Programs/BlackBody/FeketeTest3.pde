// Szem tulajdonsagai (csapok)
// http://people.inf.elte.hu/kabrabi/szinelmelet/Cones_SMJ2_E.png
float Rm=580e-9;    
float Gm=540e-9;    
float Bm=440e-9;
float Rsz=110e-9;    
float Gsz=090e-9;    
float Bsz=050e-9;
 
// T: (Nap homerseklete)
float T=5780;
 
float x=0.01; 
int n=4;
 
 
 
 
 
 
 
 
 
 
 
// negyzet    
float h(float x){
    return x*x;
    }
 
 
// Fekete Test Sugarzas, F_l is   W/m^2/m
// Vonalas szinkep
float FTS(float l){
    float B=0;
    float c=2;
    // c = vegtelen eset a 2. tagbol 1-t csinal: 1-exp() lesz. ha hozza adjuk a masik akkor az osszes 1, ami kell hogy legyen
    // c = 1 teljesen 0 lesz, ez felel meg a amximalis vesztesegnek
    // 1/c tehat a veszteseg, es 1 es vegtelen kozott lehet csak erteke
    
    
    // spektrum mod
    if(mod==0){
    B= 0.00464*0.00464*  3.1415* 2*3e8*3e8*6.626e-34* pow(l,-5.) /( exp(3e8*6.626e-34/l/T/1.38e-23)-1); 
    //if(szoras==1) B = B * (pow(l*1e7,-4.));  // szoras: az eg
    //if(szoras==2) B = B * (1-pow(l*1e7,-4.));  // szoras: a test
    
    if(szoras==1) B = B * (-exp( -x/(pow(l*4e7,n)) ) + exp( -x/((pow(l*4e7,n))*c) )) / (-(pow(l*4e7,n))*(1./((pow(l*4e7,n))*c)-1./(pow(l*4e7,n))));
    if(szoras==2) B = B * (exp(  -x/(pow(l*4e7,n)) ));  // szoras: a test  
    
    //if(szoras==1) B = B * (-exp( -x/(pow(l*1e7,4.)) ) + exp( -x/((pow(l*1e7,4.))*c) )) / (-(pow(l*1e7,4.))*(1./((pow(l*1e7,4.))*c)-1./(pow(l*1e7,4.))));
    //if(szoras==2) B = B * (exp(  -x/(pow(l*1e7,4.)) ));  // szoras: a test  
    
    
  }
    
    // vonalas mod
    if(mod==1){  
    if(l>T-(3e-9*x/150.) && l<T+(3e-9*x/150.) ){ B=1e10; }   
    
    /*
    // H    http://www.vernier.com/images/cache/figure.vsp-em.lp._hydrogen._phys-abm-21.001.800.500.png
    if(l>410e-9-3e-9 && l<410e-9+3e-9){ B=1e9; }
    if(l>434e-9-3e-9 && l<434e-9+3e-9){ B=1e10; }
    if(l>486e-9-3e-9 && l<486e-9+3e-9){ B=7e10; }
    if(l>656e-9-3e-9 && l<656e-9+3e-9){ B=1e11; }
    */
    
    }
    
    
    if(mod==2){
    B= 0.00464*0.00464*  3.1415* 2*3e8*3e8*6.626e-34* pow(l,-5.) /( exp(3e8*6.626e-34/l/ 5780 /1.38e-23)-1);
    //if(l>T-70e-9 && l<T+70e-9){ B=0; } 
    //if(l>T-T/8 && l<T+T/8){ B=0; } 
    //if(l>T-(T/8*x/150.) && l<T+(T/8*x/150.) ){ B=0; } 
    if(l>T-(60e-9*x/150.) && l<T+(60e-9*x/150.) ){ B=0; } 
    
    }
 
 
    
    return B;
    }
    
    
// Gauss fuggveny, m:max helye, sz:szoras
float Gauss(float x, float m, float sz){
     float y;
     
     y=1/sqrt(2*PI)/sz * exp( -h(x-m)/sz/sz/2 );
 
     
     return y;
     }
     
     
// which is max, legnagyobb erteket adja vissza 
float maxs(float a, float b, float c){
     float x;
     
     x=a;
     if(b>x) x=b;
     if(c>x) x=c;
     
     return x;
     }
     
    
 
 
// Integralas, RungeKutta modszer
float integral(float s, float tig, float m, float sz){    
 float v=0;
 float a;
 float ce;
 float ck;
 float ch;
 float cn;
 float i;
 
 float dt=0.1e-9;
  
  tig=tig/dt;
  s=s/dt;
 
  //RK ciklus
  for(i=s;i<tig;i++){
 
     a=FTS(i*dt)*Gauss(i*dt,m,sz);
     ce=dt*a*dt;
 
     a=FTS((i+0.5)*dt)*Gauss((i+0.5)*dt,m,sz);
     ck=dt*a*dt;
 
     a=FTS((i+0.5)*dt)*Gauss((i+0.5)*dt,m,sz);
     ch=dt*a*dt;
 
     a=FTS((i+1.0)*dt)*Gauss((i+1.0)*dt,m,sz);
     cn=dt*a*dt;
 
       if(v*0==0) v=v+(1.0/6.0)*ce+(1.0/3.0)*ck+(1.0/3.0)*ch+(1.0/6.0)*cn;  // RungeKutta modszer
        
       
       }  
 
 return v;    
}    
 
 
 
 
  float R; // piros   
  float G; // zold   
  float B; // kek
  
  float R2; // piros   
  float G2; // zold   
  float B2; // kek
  
  float szog=0; // szog a forgashoz
  float flux; // max fluxus
  float flux2; // max fluxus
  
  int mod=0; // mod
  int szoras=0; // szoras
 
 
 
void setup() {
  size(500, 500);
  //size(500, 500, P3D);
  
  //frameRate(10);
  
  textSize(32);
  
  // mod=, 0 FTS, 1 vonalas, 2 FTS-vonalas
  // szoras=, 0 nincs, 1 eg, 2 test
  mod=0;
  szoras=0;
  
  
  
 
}
 
void draw() {
 
  if(mod>2) mod=mod-3;
  if(szoras>2) szoras=szoras-3;
  
  
if(mod==0){
  T=30*mouseX+800;
  if (mousePressed  && mod==0 && szoras==0) {
    T=5780;
  }
  
  
}
  
if(mod==1 || mod==2){
    T=1.2e-9*mouseX+300e-9;
  if (mousePressed) {
    //T=550e-9;
  }
}


  //x=0.001*mouseY;
  if(mod==0) x=15*exp(0.03*(mouseY-175));
  if(mod==2) x=15*exp(0.0130*(mouseY+50));
  //if(mod==2) x=15*exp(0.0125*(mouseY-150));
  
  
  //T=40000;
 
  if(szoras!=0) szoras=1;
  if(mod!=0) mod=1;
  // Integralas:
  R=integral(100e-9,1000e-9,Rm,Rsz);
  G=integral(100e-9,1000e-9,Gm,Gsz);
  B=integral(100e-9,1000e-9,Bm,Bsz);
  
  if(szoras!=0){ szoras=2;
  // Integralas 2:
  R2=integral(100e-9,1000e-9,Rm,Rsz);
  G2=integral(100e-9,1000e-9,Gm,Gsz);
  B2=integral(100e-9,1000e-9,Bm,Bsz);
  
  R2=R2/0.0158752172;
  G2=G2/0.0165422747;
  B2=B2/0.0164163495;
  flux2=23*log(maxs(R2,G2,B2)/0.00464/0.00464)/log(10)+92+30;
  if(flux2>255) flux2=255; if(flux2<0) flux2=0;
  }
  
  
  if(mod!=0){ mod=2;
  // Integralas 2:
  R2=integral(100e-9,1000e-9,Rm,Rsz);
  G2=integral(100e-9,1000e-9,Gm,Gsz);
  B2=integral(100e-9,1000e-9,Bm,Bsz);
  
  R2=R2/0.0158752172;
  G2=G2/0.0165422747;
  B2=B2/0.0164163495;
  flux2=23*log(maxs(R2,G2,B2)/0.00464/0.00464)/log(10)+92+30;
  if(flux2>255) flux2=255; if(flux2<0) flux2=0;
  }
 
 
// Feher egyensuly
  R=R/0.0158752172;
  G=G/0.0165422747;
  B=B/0.0164163495;
  
  flux=23*log(maxs(R,G,B)/0.00464/0.00464)/log(10)+92+30;
  //flux=log(maxs(R,G,B)/1e-10)/log(10)*19.32;
  if(flux>255) flux=255; if(flux<0) flux=0;
  
  //if(szoras==0 || mod!=0){ R2=R; G2=G; B2=B; flux2=flux; }
  if(szoras==0 && mod==0){ R2=R; G2=G; B2=B; flux2=flux; }
  
  

  
  
  
  
  
  
  background(100);
  
  
  ellipse(250,250,400,400);
  
  //szog=szog+0.01;
  //translate(250, 250, 0);
  //rotateY(szog);
  //sphere(200);
  
  //Belso kor!
  //fill(R,G,B);
  //fill(int(255*R/maxs(R,G,B)+0.5) ,int(255*G/maxs(R,G,B)+0.5) , int(255*B/maxs(R,G,B)+0.5));
  fill(int(flux2*R2/maxs(R2,G2,B2)+0.5) ,int(flux2*G2/maxs(R2,G2,B2)+0.5) , int(flux2*B2/maxs(R2,G2,B2)+0.5));
  
  
  //Kulso kor!
  ellipse(250,250,150,150);
  fill(int(flux*R/maxs(R,G,B)+0.5) ,int(flux*G/maxs(R,G,B)+0.5) , int(flux*B/maxs(R,G,B)+0.5));
  
  
  
  
  
  if(mod==1 || mod==2){ T=T*1e9; }
  
  //println(int(T));
  //println(int(flux));
  
  //6e-9*x/150.
  
  text(int(T), 10, 30, 0); 
  if(szoras!=0 && mod==0) text(x/15000., 40, 60, 0);
  if(mod==1)              text(int(6*x/150.), 10, 60, 0);
  if(mod==2)              text(int(120*x/150.), 10, 60, 0);
  if(szoras!=0 && mod==0) text(int(n), 10, 90, 0);
  if(mod==0)              text("K", 120, 30, 0); 
  if(mod==1 || mod==2)  { text("nm", 90, 30, 0);  text("nm", 90, 60, 0); }
  if(szoras!=0 && mod==0) text("X=", 5, 60, 0);
  
  if(szoras==1 && mod==0) text("Szoras: eg", 280, 30, 0); 
  if(szoras==2 && mod==0) text("Szoras: test", 280, 30, 0); 
  if(mod!=0)              text("(Kozponti, vastagsag)", 160, 30, 0);
  if(szoras!=0 && mod==0 && n==4) text("(Rayleigh)", 8, 460, 0);
  if(szoras!=0 && mod==0 && n<4 && n!=0) text("(Mie)", 8, 460, 0);
  if(szoras!=0 && mod==0 && n==0) text("(Allando)", 8, 460, 0);
  if(szoras!=0 && mod==0 ) text("(0 - 4)", 380, 460, 0);
  
  if(mod==0) text("Fekete test sugarzas (space, s)", 8, 490, 0); 
  if(mod==2){ text("Monokromatikus (space)", 10, 490, 0); }
  //if(mod==2){ text("Feher - Mono (space)", 10, 490, 0); }
  
 
 
  
}
 
 
 
 
// valt a modok kozott SPACE-el, szoras kozott s-el
public void keyPressed()
{
     
    if (key == ' ')
    {
      mod=mod+1;
    }
    
    if (key == 's')
    {
      szoras=szoras+1;
    }
    
    
    
    
            if (key == '0')
    {
      n=0;
    }
    
            if (key == '1')
    {
      n=1;
    }
    
            if (key == '2')
    {
      n=2;
    }
    
            if (key == '3')
    {
      n=3;
    }
    
            if (key == '4')
    {
      n=4;
    }
    
}