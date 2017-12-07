%--------------------------------------------------------------------------
% Replicando: Garcia-Cicco y Kawamura (2015)
% Dealing with Dutch Disease: Fiscal Rules and Macro-Prudential Policies
% Codigo elaborado por: Carlos Rojas Quiroz
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
var lambda cR cNR hR hNR w p rL c cN cT cX cM pN pT yX z hX kX yN hN kN uX 
uN qX qN nX nN lX lcN invX invN pI inv xN xM g dgstar rev drstar lR dstar 
imp expt tb gdp rer rstar gdpm sNs sX sC RP
aX aN pX yCo pCo rW;
%--------------------------------------------------------------------------
varexo epsaX epsaN epspX epsyCo epspCo epsrW; 
%--------------------------------------------------------------------------
% Parametros estructurales
parameters theta uve xi alphaX alphaN delta gamma phiD kappa sCog sCostar 
etaR nu psi mu epsilon varepsilonX varepsilonN phiX phiN  etaC etarev;
%--------------------------------------------------------------------------
% Meta de estado estacionario
parameters stb sCo sg sN lev;
%--------------------------------------------------------------------------
% Estado estacionario
parameters tauCo_ss dgstarss rp_ss rWss pCoss pXss aNss pNss hss hNss hXss hNRss hRss
beta rstarss rLss pIss pTss qXss qNss uXss uNss kNss yNss wss kXss aXss yXss pss
zss invXss invNss invss xNss xMss gdpmss gss yCoss tbss cNss cMss cXss cTss impss
exptss phi css gdpss rerss dstar_ss dbar nXss nNss lXss lcNss iotaN iotaX sCoR
varsigma cNRss cRss lambdass rev_ss tau_ss drstarss lRss sNss sXss sCss RPss;
%--------------------------------------------------------------------------
% Parametros exogenos
parameters rhoaX rhoaN rhopX rhoyCo rhopCo rhorW;
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
% Calibracion de parametros estructurales
%--------------------------------------------------------------------------
theta       =   2.0;
uve         =   1.6;
xi          =   0.5145;
alphaX      =   0.36;
alphaN      =   0.65;
delta       =   0.015;
gamma       =   0.40;
phiD        =   0.001;
kappa       =   0.50;
sCog        =   0.40;
sCostar     =   0.60;
etaR        =   0.001;
nu          =   0.97;
epsilon     =   0.9813;
psi         =   0.3416;
mu          =   0.5921;
varepsilonX =   0.0169; 
varepsilonN =   0.1150;
phiX        =   5.4315;
phiN        =   5.4315;
etarev      =   0.00;
rhoaX       =   0.75; 
rhoaN       =   0.75;
rhopX       =   0.75;
rhoyCo      =   0.75;
rhopCo      =   0.8399;
rhorW       =   0.75;
%--------------------------------------------------------------------------
% Niveles meta de estado estacionario
%--------------------------------------------------------------------------
stb         =   0.04;
sCo         =   0.10;
sg          =   0.11;
sN          =   0.50;
lev         =   2.05;
etaC        =   0.00;
%--------------------------------------------------------------------------
% Estado estacionario
%--------------------------------------------------------------------------
tauCo_ss    =   0.35;
dgstarss    =   0.00;
rp_ss       =   1.0123;
rWss        =   1.0148;
pCoss       =   1;
pXss        =   1; 
aNss        =   1;
pNss        =   1;
hss         =   0.3;
hNss        =   0.15;
hXss        =   hss-hNss;
hNRss       =   0.3;
hRss        =   0.3;
beta        =   1/(rWss);
rstarss     =   rWss;
rLss        =   rstarss;
pIss        =   pNss^gamma;
pTss        =   pXss^xi;
qXss        =   pIss;
qNss        =   pIss;
uXss        =   qXss*((rLss)*rp_ss-1+delta);
uNss        =   qNss*((rLss)*rp_ss-1+delta);
kNss        =   (uNss/(pNss*(1-alphaN)*aNss))^(-1/alphaN)*hNss;
yNss        =   aNss*(hNss)^alphaN*kNss^(1-alphaN);
wss         =   pNss*alphaN*yNss/hNss;
kXss        =   (wss*(1-alphaX-psi)/(uXss*alphaX))*hXss;
aXss        =   (wss/(pXss*alphaX))^(1-psi)*(hXss/kXss)^(1-alphaX-psi);
yXss        =   aXss^(1/(1-psi))*hXss^(alphaX/(1-psi))*kXss^((1-alphaX-psi)/(1-psi));
zss         =   yXss;
invXss      =   delta*kXss;
invNss      =   delta*kNss;
invss       =   invXss+invNss;
xNss        =   gamma*(pIss/pNss)*invss;
xMss        =   (1-gamma)*pIss*invss;
gdpmss      =   (pXss*yXss+pNss*yNss)/(1-sCo);
gss         =   sg*gdpmss;
yCoss       =   sCo*gdpmss/pCoss;
tbss        =   stb*gdpmss;
cNss        =   yNss-xNss-gss;
cMss        =   (1-xi)*(pXss*yXss+pCoss*yCoss-xMss-tbss);
cXss        =   xi/(1-xi)*cMss/pXss;
cTss        =   (cXss/xi)^xi*(cMss/(1-xi))^(1-xi);
impss       =   cMss+xMss;
exptss      =   pXss*(yXss-cXss)+pCoss*yCoss;
phi         =   (1+(pTss/pNss)^epsilon*cTss/cNss)^(-1);
css         =   (phi^(1/epsilon)*cNss^(1-1/epsilon)+(1-phi)^(1/epsilon)*(cTss)^(1-1/epsilon))^(epsilon/(epsilon-1));
pss         =   pNss*(cNss/(phi*css))^(1/epsilon);
gdpss       =   gdpmss/pss;
rerss       =   1/pss;
dstar_ss    =   (tbss-pCoss*yCoss*sCostar*(1-tauCo_ss))/(rstarss-1);
drstarss    =   (dstar_ss-dgstarss)/(1-kappa);
dbar        =   dstar_ss/gdpmss;
nXss        =   qXss*kXss/lev;
nNss        =   qNss*kNss/lev;
lXss        =   (qXss*kXss-nXss)/pss;
lcNss       =   (qNss*kNss-nNss)/pss;
lRss        =   (lXss+lcNss)/(1-kappa);
iotaX       =   nXss - nu*((uXss+(1-delta)*qXss)*kXss-pss*lXss*(rLss));
iotaN       =   nNss - nu*((uNss+(1-delta)*qNss)*kNss-pss*lcNss*(rLss));
sCoR        =   1-sCog-sCostar;
rev_ss      =   pNss*gss;
tau_ss      =   (rev_ss-pCoss*yCoss*(tauCo_ss*(sCoR+sCostar)+sCog))/(pXss*yXss+pNss*yNss);
varsigma    =   (1-tau_ss)*wss/pss*1/(hNRss^uve);
cNRss       =   (1-tau_ss)*wss*hNRss/pss;
cRss        =   (css-kappa*cNRss)/(1-kappa);
lambdass    =   (cRss-varsigma*hRss^(1+uve)/(1+uve))^(-theta);
sNss        =   pNss*yNss/(pXss*yXss+pNss*yNss);
sXss        =   pXss*yXss/(pXss*yXss+pNss*yNss);
sCss        =   pss*css/(pXss*yXss+pNss*yNss);
RPss        =   rp_ss;
%--------------------------------------------------------------------------
model;
% Hogares ricardianos (4)
exp(lambda)                          =   (exp(cR)-varsigma*exp(hR)^(1+uve)/(1+uve))^(-theta);
exp(lambda)/exp(p)*(1-tau_ss)*exp(w) =   (exp(cR)-varsigma*exp(hR)^(1+uve)/(1+uve))^(-theta)*varsigma*exp(hR)^uve;
exp(lambda)/exp(p)                   =   beta*(exp(rstar))*exp(lambda(+1))/exp(p(+1));
exp(lambda)                          =   beta*(exp(rL))*exp(lambda(+1));
% Hogares no ricardianos (2)
(exp(cNR)-varsigma*exp(hNR)^(1+uve)/(1+uve))^(-theta)*varsigma*exp(hNR)^uve  =  (exp(cNR)-varsigma*exp(hNR)^(1+uve)/(1+uve))^(-theta)/exp(p)*(1-tau_ss)*exp(w);
exp(p)*exp(cNR)                      =   (1-tau_ss)*exp(w)*exp(hNR);
% Consumo agregado (6)
exp(c)                               =   (phi^(1/epsilon)*(exp(cN))^(1-1/epsilon)+(1-phi)^(1/epsilon)*(exp(cT))^(1-1/epsilon))^(epsilon/(epsilon-1));
exp(cT)                              =   (exp(cX)/xi)^xi*(exp(cM)/(1-xi))^(1-xi);
exp(cN)                              =   phi*(exp(p)/exp(pN))^epsilon*exp(c);
exp(cT)                              =   (1-phi)*(exp(p)/exp(pT))^epsilon*exp(c);
exp(cX)                              =   xi*(exp(pT)/exp(pX))*exp(cT);
exp(cM)                              =   (1-xi)*exp(pT)*exp(cT);
% Produccion de transables (4)
exp(yX)                              =   exp(aX)*exp(z)^psi*exp(hX)^alphaX*exp(kX(-1))^(1-alphaX-psi);
exp(z)                               =   exp(z(-1))^mu*(exp(yX(-1)))^(1-mu);
exp(w)                               =   exp(pX)*alphaX*exp(yX)/exp(hX);
exp(uX)                              =   exp(pX)*(1-alphaX-psi)*exp(yX)/exp(kX(-1));
% Produccion de no transables (3)
exp(yN)                              =   exp(aN)*exp(hN)^alphaN*exp(kN(-1))^(1-alphaN);
exp(w)                               =   exp(pN)*alphaN*exp(yN)/exp(hN);
exp(uN)                              =   exp(pN)*(1-alphaN)*exp(yN)/exp(kN(-1));
% Empresarios (6)
exp(qX)*exp(kX)                                 =   exp(nX)+exp(p)*exp(lX);
(exp(uX(+1))+(1-delta)*exp(qX(+1)))/exp(qX)     =   (exp(rL))*rp_ss*(exp(qX)*exp(kX)/exp(nX)*1/lev)^(varepsilonX);
exp(nX)                                         =   nu*((exp(uX)+(1-delta)*exp(qX))*exp(kX(-1))-exp(p)*exp(lX(-1))*(exp(rL(-1))))+iotaX;
exp(qN)*exp(kN)                                 =   exp(nN)+exp(p)*exp(lcN);
(exp(uN(+1))+(1-delta)*exp(qN(+1)))/exp(qN)     =   (exp(rL))*rp_ss*(exp(qN)*exp(kN)/exp(nN)*1/lev)^(varepsilonN);
exp(nN)                                         =   nu*((exp(uN)+(1-delta)*exp(qN))*exp(kN(-1))-exp(p)*exp(lcN(-1))*(exp(rL(-1))))+iotaN;
% Capital e Inversion (7)
exp(kX)                              =   (1-delta)*exp(kX(-1))+(1-phiX/2*(exp(invX)/exp(invX(-1))-1)^2)*exp(invX);
exp(pI)                              =   exp(qX)*(1-phiX/2*(exp(invX)/exp(invX(-1))-1)^2-phiX*(exp(invX)/exp(invX(-1))-1)*exp(invX)/exp(invX(-1))*1/exp(invX(-1))) + beta*exp(lambda(+1))/exp(lambda)*exp(qX(+1))*phiX*(exp(invX(+1))/exp(invX)-1)*(exp(invX(+1))/exp(invX))^2*(1/exp(invX))^2;
exp(kN)                              =   (1-delta)*exp(kN(-1))+(1-phiN/2*(exp(invN)/exp(invN(-1))-1)^2)*exp(invN);
exp(pI)                              =   exp(qN)*(1-phiN/2*(exp(invN)/exp(invN(-1))-1)^2-phiN*(exp(invN)/exp(invN(-1))-1)*exp(invN)/exp(invN(-1))*1/exp(invN(-1))) + beta*exp(lambda(+1))/exp(lambda)*exp(qN(+1))*phiN*(exp(invN(+1))/exp(invN)-1)*(exp(invN(+1))/exp(invN))^2*(1/exp(invN))^2;
exp(inv)                             =   (exp(xN)/gamma)^gamma*(exp(xM)/(1-gamma))^(1-gamma);
exp(xN)                              =   gamma*(exp(pI)/exp(pN))*exp(inv);
exp(xM)                              =   (1-gamma)*exp(pI)*exp(inv);
% Politica fiscal (3)
exp(rev)+dgstar                      =   exp(pN)*exp(g)+dgstar(-1)*(exp(rstar(-1)));
exp(rev)                             =   tau_ss*(exp(pX)*exp(yX)+exp(pN)*exp(yN))+exp(pCo)*exp(yCo)*(tauCo_ss*(sCoR+sCostar)+sCog);
exp(pN)*exp(g)+dgstar(-1)*(exp(rstar(-1))-1+etaR) = etaC+rev_ss+etarev*(exp(rev)-rev_ss);
% Agregacion y limpieza de mercados (14)
exp(hX)+exp(hN)                      =   (1-kappa)*exp(hR)+kappa*exp(hNR);
exp(c)                               =   (1-kappa)*exp(cR)+kappa*exp(cNR);
dstar                                =   (1-kappa)*drstar+kappa*dgstar;
(1-kappa)*exp(lR)                    =   exp(lX)+exp(lcN);
exp(inv)                             =   exp(invN)+exp(invX);
exp(yN)                              =   exp(cN)+exp(xN)+exp(g);
exp(imp)                             =   exp(cM)+exp(xM);
exp(expt)                            =   exp(pX)*(exp(yX)-exp(cX))+exp(pCo)*exp(yCo);
exp(tb)                              =   exp(expt)-exp(imp);
exp(p)*exp(gdp)                      =   exp(pX)*exp(yX)+exp(pN)*exp(yN)+exp(pCo)*exp(yCo);
exp(rer)                             =   1/exp(p);
(dstar(-1))*(exp(rstar(-1)))         =   (dstar)+exp(tb)-exp(pCo)*exp(yCo)*sCostar*(1-tauCo_ss);
exp(rstar)                           =   exp(rW)+exp(phiD*((dstar)-dstar_ss)/dstar_ss)-1;
exp(gdpm)                            =   exp(p)*exp(gdp);
% Choques exogenos (6)
aX                                   =  (1-rhoaX)*log(aXss)+rhoaX*aX(-1)+epsaX; 
aN                                   =  (1-rhoaN)*log(aNss)+rhoaN*aN(-1)+epsaN;  
pX                                   =  (1-rhopX)*log(pXss)+rhopX*pX(-1)+epspX;  
yCo                                  =  (1-rhoyCo)*log(yCoss)+rhoyCo*yCo(-1)+epsyCo;  
pCo                                  =  (1-rhopCo)*log(pCoss)+rhopCo*pCo(-1)+epspCo;  
rW                                   =  (1-rhorW)*log(rWss)+rhorW*rW(-1)+epsrW;
% Variables auxiliares (4)
exp(sX)                              =  exp(pX)*exp(yX)/(exp(pX)*exp(yX)+exp(pN)*exp(yN));
exp(sNs)                             =  exp(pN)*exp(yN)/(exp(pX)*exp(yX)+exp(pN)*exp(yN));
exp(sC)                              =  exp(p)*exp(c)/(exp(pX)*exp(yX)+exp(pN)*exp(yN));
exp(RP)                              =  (exp(lX)*rp_ss*(exp(qX)*exp(kX)/exp(nX)*1/lev)^(varepsilonX) + exp(lcN)*rp_ss*(exp(qN)*exp(kN)/exp(nN)*1/lev)^(varepsilonN))/(exp(lX)+exp(lcN));
end;

% Estado estacionario
steady_state_model;
lambda              =log(lambdass); 
cR                  =log(cRss); 
cNR                 =log(cNRss); 
hR                  =log(hRss); 
hNR                 =log(hNRss); 
w                   =log(wss);
p                   =log(pss);
rL                  =log(rLss); 
c                   =log(css);
cN                  =log(cNss); 
cT                  =log(cTss); 
cX                  =log(cXss); 
cM                  =log(cMss); 
pN                  =log(pNss); 
pT                  =log(pTss); 
yX                  =log(yXss); 
z                   =log(zss);
hX                  =log(hXss); 
kX                  =log(kXss); 
yN                  =log(yNss); 
hN                  =log(hNss); 
kN                  =log(kNss); 
uX                  =log(uXss); 
uN                  =log(uNss); 
qX                  =log(qXss); 
qN                  =log(qNss); 
nX                  =log(nXss); 
nN                  =log(nNss); 
lX                  =log(lXss); 
lcN                 =log(lcNss); 
invX                =log(invXss); 
invN                =log(invNss);
pI                  =log(pIss); 
inv                 =log(invss); 
xN                  =log(xNss); 
xM                  =log(xMss); 
g                   =log(gss); 
dgstar              =(dgstarss); 
rev                 =log(rev_ss);
drstar              =(drstarss);
lR                  =log(lRss); 
dstar               =(dstar_ss); 
imp                 =log(impss);
expt                =log(exptss);
tb                  =log(tbss); 
gdp                 =log(gdpss); 
rer                 =log(rerss);
rstar               =log(rstarss);
gdpm                =log(gdpmss);
aX                  =log(aXss); 
aN                  =log(aNss); 
pX                  =log(pXss); 
yCo                 =log(yCoss); 
pCo                 =log(pCoss); 
rW                  =log(rWss);
sNs                 =log(sNss);
sX                  =log(sXss);
sC                  =log(sCss);
RP                  =log(RPss);
end;

shocks;
var epsaX   ; stderr 0.00; 
var epsaN   ; stderr 0.00; 
var epspX   ; stderr 0.00; 
var epsyCo  ; stderr 0.00; 
var epspCo  ; stderr 0.1113*100; 
var epsrW   ; stderr 0.00;
end;
 
resid;
steady;
check;

stoch_simul(order = 1, irf=24) pCo sX sNs rer sC RP;
