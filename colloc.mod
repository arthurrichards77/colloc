# 1 "colloc.mot"
# 1 "<command-line>"
# 1 "colloc.mot"



## WARNING - don't edit the MOD file directly







param Na; ## number of agents
param Nc; ## number of collocation points
param Ne; ## number of constraint evaluation points
param Nf; ## number of finite elements
param Nz; ## number of outputs

param D{1..Nc,1..Nc}; ## differentiator matrix
param D2{1..Nc,1..Nc}; ## double differentiator matrix
param E{1..Ne,1..Nc}; ## evaluation matrix
param Ed{1..Ne,1..Nc}; ## evaluation matrix for derivatives
param Edd{1..Ne,1..Nc}; ## evaluation matrix for double derivatives

param zi{1..Nz,1..Na}; ## initial states
param zf{1..Nz,1..Na}; ## terminal states
param zdi{1..Nz,1..Na}; ## initial state derivatives
param zdf{1..Nz,1..Na}; ## terminal state derivatives

param zguess{1..Nz,1..Nc*Nf*Na};
param tauguess{1..Na};
param telguess{1..Na,1..Nf};
param yguess{1..(1+Nz)};

param Nu=3;
param ulo{1..Nu};
param uhi{1..Nu};

var zc{ii in 1..Nz, cc in 1..Nc, ff in 1..Nf, aa in 1..Na} := zguess[ii,cc+Nc*(ff-1)+Nc*Nf*(aa-1)]; ## outputs
var tau{aa in 1..Na} >= 0, := tauguess[aa];
var tel{aa in 1..Na, ff in 1..Nf} >=0, := telguess[aa,ff]; ## element start time

##var slk >=0, := 1;

param R = 3;
param H = 1;

param Np >=0; ## number of avoidance pairs to enforce
param ap{1..(Np+1),1..2}; ## pairs for which to enforce avoidance (one more than Np rows to avoid empty matrix if no pairs enforced)

param Npold >=0; ## number of avoidance pairs to load
param guessye{pp in 1..Npold, e1 in 2..Ne, f1 in 1..Nf, e2 in 2..Ne, f2 in 1..Nf, ii in 1..4};

var ye{pp in 1..Np, e1 in 2..Ne, f1 in 1..Nf, e2 in 2..Ne, f2 in 1..Nf, ii in 1..4} := if pp<=Npold then guessye[pp,e1,f1,e2,f2,ii] else yguess[ii];

## OBJECTIVE ****************************************************

minimize time: sum{ff in 1..Nf, aa in 1..Na} tau[aa] + sum{pp in 1..Np, e1 in 2..Ne, f1 in 1..Nf, e2 in 2..Ne, f2 in 1..Nf, ii in 1..4} (0.0001/(Ne*Ne*Nf*Nf*Np*4))*ye[pp,e1,f1,e2,f2,ii]*ye[pp,e1,f1,e2,f2,ii]; ##

## CONSTRAINTS **************************************************

## initial and terminal constraints
subject to initz{ii in 1..Nz, aa in 1..Na}: (sum{cc in 1..Nc} E[1,cc]*zc[ii,cc,1,aa]) = zi[ii,aa];
subject to termz{ii in 1..Nz, aa in 1..Na}: (sum{cc in 1..Nc} E[Ne,cc]*zc[ii,cc,Nf,aa]) = zf[ii,aa];
subject to initzdot{ii in 1..Nz, aa in 1..Na}: (sum{cc in 1..Nc} Ed[1,cc]*zc[ii,cc,1,aa]) = (tau[aa]/2)*zdi[ii,aa];
subject to termzdot{ii in 1..Nz, aa in 1..Na}: (sum{cc in 1..Nc} Ed[Ne,cc]*zc[ii,cc,Nf,aa])= (tau[aa]/2)*zdf[ii,aa];

## continuity across elements
subject to contz{ii in 1..Nz, ff in 2..Nf, aa in 1..Na}: (sum{cc in 1..Nc} E[1,cc]*zc[ii,cc,ff,aa]) = (sum{cc in 1..Nc} E[Ne,cc]*zc[ii,cc,ff-1,aa]);
subject to contzdot{ii in 1..Nz, ff in 2..Nf, aa in 1..Na}: (sum{cc in 1..Nc} Ed[1,cc]*zc[ii,cc,ff,aa]) = (sum{cc in 1..Nc} Ed[Ne,cc]*zc[ii,cc,ff-1,aa]);

## element timing continuity
subject to contt{aa in 1..Na, ff in 1..Nf}: tel[aa,ff] = sum{f1 in 1..(ff-1)} tau[aa];
##subject to contt{aa in 1..Na, ff in 1..Nf}: tel[aa,ff] >= sum{f1 in 1..(ff-1)} tau[aa];

## differential flatness speed limit
subject to maxspd{ee in 1..Ne, ff in 1..Nf, aa in 1..Na}: (sum{cc in 1..Nc} Ed[ee,cc]*zc[1,cc,ff,aa])*(sum{cc in 1..Nc} Ed[ee,cc]*zc[1,cc,ff,aa]) + (sum{cc in 1..Nc} Ed[ee,cc]*zc[2,cc,ff,aa])*(sum{cc in 1..Nc} Ed[ee,cc]*zc[2,cc,ff,aa]) <= (tau[aa]/2)*(tau[aa]/2)*uhi[1]*uhi[1];
subject to minspd{ee in 1..Ne, ff in 1..Nf, aa in 1..Na}: (sum{cc in 1..Nc} Ed[ee,cc]*zc[1,cc,ff,aa])*(sum{cc in 1..Nc} Ed[ee,cc]*zc[1,cc,ff,aa]) + (sum{cc in 1..Nc} Ed[ee,cc]*zc[2,cc,ff,aa])*(sum{cc in 1..Nc} Ed[ee,cc]*zc[2,cc,ff,aa]) >= (tau[aa]/2)*(tau[aa]/2)*ulo[1]*ulo[1];

## curvature limits
subject to turn1{ee in 1..Ne, ff in 1..Nf, aa in 1..Na}: (sum{cc in 1..Nc} Edd[ee,cc]*zc[1,cc,ff,aa])*(sum{cc in 1..Nc} Ed[ee,cc]*zc[2,cc,ff,aa]) - (sum{cc in 1..Nc} Edd[ee,cc]*zc[2,cc,ff,aa])*(sum{cc in 1..Nc} Ed[ee,cc]*zc[1,cc,ff,aa]) <= uhi[2]*((sum{cc in 1..Nc} Ed[ee,cc]*zc[1,cc,ff,aa])*(sum{cc in 1..Nc} Ed[ee,cc]*zc[1,cc,ff,aa]) + (sum{cc in 1..Nc} Ed[ee,cc]*zc[2,cc,ff,aa])*(sum{cc in 1..Nc} Ed[ee,cc]*zc[2,cc,ff,aa]))^(3/2);
subject to turn2{ee in 1..Ne, ff in 1..Nf, aa in 1..Na}: (sum{cc in 1..Nc} Edd[ee,cc]*zc[1,cc,ff,aa])*(sum{cc in 1..Nc} Ed[ee,cc]*zc[2,cc,ff,aa]) - (sum{cc in 1..Nc} Edd[ee,cc]*zc[2,cc,ff,aa])*(sum{cc in 1..Nc} Ed[ee,cc]*zc[1,cc,ff,aa]) >= -uhi[2]*((sum{cc in 1..Nc} Ed[ee,cc]*zc[1,cc,ff,aa])*(sum{cc in 1..Nc} Ed[ee,cc]*zc[1,cc,ff,aa]) + (sum{cc in 1..Nc} Ed[ee,cc]*zc[2,cc,ff,aa])*(sum{cc in 1..Nc} Ed[ee,cc]*zc[2,cc,ff,aa]))^(3/2);

## ********************* AVOIDANCE STARTS HERE *************************

subject to avoid11{pp in 1..Np, e1 in 2..Ne, f1 in 1..Nf, e2 in 2..Ne, f2 in 1..Nf}: sum{ii in 1..3} (ye[pp,e1,f1,e2,f2,ii]*((sum{cc in 1..Nc} E[e1-1,cc]*zc[ii,cc,f1,(ap[pp,1])]) - (sum{cc in 1..Nc} E[e2-1,cc]*zc[ii,cc,f2,(ap[pp,2])]))) + (ye[pp,e1,f1,e2,f2,4]*((tel[(ap[pp,1]),f1] + tau[(ap[pp,1])]*(e1-1 -1)/(Ne-1)) - (tel[(ap[pp,2]),f2] + tau[(ap[pp,2])]*(e2-1 -1)/(Ne-1)))) >= 1; ## - slk;
subject to avoid12{pp in 1..Np, e1 in 2..Ne, f1 in 1..Nf, e2 in 2..Ne, f2 in 1..Nf}: sum{ii in 1..3} (ye[pp,e1,f1,e2,f2,ii]*((sum{cc in 1..Nc} E[e1,cc]*zc[ii,cc,f1,(ap[pp,1])]) - (sum{cc in 1..Nc} E[e2-1,cc]*zc[ii,cc,f2,(ap[pp,2])]))) + (ye[pp,e1,f1,e2,f2,4]*((tel[(ap[pp,1]),f1] + tau[(ap[pp,1])]*(e1-1)/(Ne-1)) - (tel[(ap[pp,2]),f2] + tau[(ap[pp,2])]*(e2-1 -1)/(Ne-1)))) >= 1; ## - slk;
subject to avoid21{pp in 1..Np, e1 in 2..Ne, f1 in 1..Nf, e2 in 2..Ne, f2 in 1..Nf}: sum{ii in 1..3} (ye[pp,e1,f1,e2,f2,ii]*((sum{cc in 1..Nc} E[e1-1,cc]*zc[ii,cc,f1,(ap[pp,1])]) - (sum{cc in 1..Nc} E[e2,cc]*zc[ii,cc,f2,(ap[pp,2])]))) + (ye[pp,e1,f1,e2,f2,4]*((tel[(ap[pp,1]),f1] + tau[(ap[pp,1])]*(e1-1 -1)/(Ne-1)) - (tel[(ap[pp,2]),f2] + tau[(ap[pp,2])]*(e2-1)/(Ne-1)))) >= 1; ## - slk;
subject to avoid22{pp in 1..Np, e1 in 2..Ne, f1 in 1..Nf, e2 in 2..Ne, f2 in 1..Nf}: sum{ii in 1..3} (ye[pp,e1,f1,e2,f2,ii]*((sum{cc in 1..Nc} E[e1,cc]*zc[ii,cc,f1,(ap[pp,1])]) - (sum{cc in 1..Nc} E[e2,cc]*zc[ii,cc,f2,(ap[pp,2])]))) + (ye[pp,e1,f1,e2,f2,4]*((tel[(ap[pp,1]),f1] + tau[(ap[pp,1])]*(e1-1)/(Ne-1)) - (tel[(ap[pp,2]),f2] + tau[(ap[pp,2])]*(e2-1)/(Ne-1)))) >= 1; ## - slk;

param Nth = 13;
param pi = 3.14159265;
subject to dual1{th in 1..20, pp in 1..Np, e1 in 2..Ne, f1 in 1..Nf, e2 in 2..Ne, f2 in 1..Nf}: R*(cos(2*pi*th/Nth)*ye[pp,e1,f1,e2,f2,1] + sin(2*pi*th/Nth)*ye[pp,e1,f1,e2,f2,2]) + H*ye[pp,e1,f1,e2,f2,3] <= 1;
subject to dual2{th in 1..20, pp in 1..Np, e1 in 2..Ne, f1 in 1..Nf, e2 in 2..Ne, f2 in 1..Nf}: R*(cos(2*pi*th/Nth)*ye[pp,e1,f1,e2,f2,1] + sin(2*pi*th/Nth)*ye[pp,e1,f1,e2,f2,2]) - H*ye[pp,e1,f1,e2,f2,3] <= 1;

##subject to sense{pp in 1..Np, e1 in 2..Ne, f1 in 1..Nf, e2 in 2..Ne, f2 in 1..Nf:pp==3}: ye[pp,e1,f1,e2,f2,4] <= 0;
