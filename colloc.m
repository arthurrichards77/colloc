clear all
close all
clc

% start the clock
tic

%% %%%%%%%%%%%%%% build the model file %%%%%%%%%%%%%%%
!runcpp colloc.mot colloc.mod

%% %%%%%%%%%%%%%% set problem data %%%%%%%%%%%%%%%%

% max speed
maxSpeed = 2.0;

% max turn rate
maxTurn = 0.25;

% max climb
maxClimb = 0.4;

% initial configuration (x,y,h,theta)
initConfig = [0 0 0 0*pi/4;
    0  8 0 0;
    0 15 0 1*pi/4;
    0 20 0 0*pi/4;
    10 20 0 -pi/2];

% final configuration
termConfig = [20 15 0 -0*pi/4;
    20 8 0 0;
    25 0 0 -0*pi/4;
    20 0 0 0;
    10 0 0 -pi/2];

% number of vehicles
nVehs = size(initConfig,1);

% or choose cut down number
nVehs = 4;

% sense constraints - N x 2, each row [a b] implies a before b
senseCons = [];
%senseCons = [3 1]; % example for with nVehs=3
%senseCons = [1 3];

% number of sense constraints
nSense = size(senseCons,1);

% and downselect
initConfig = initConfig(1:nVehs,:);
termConfig = termConfig(1:nVehs,:);

% convert to z (output) and z-dot form
initZ = initConfig(:,1:3);
initZdot = 0.9*maxSpeed*[cos(initConfig(:,4)) sin(initConfig(:,4)) 0*initConfig(:,4)];
termZ = termConfig(:,1:3);
termZdot = 0.9*maxSpeed*[cos(termConfig(:,4)) sin(termConfig(:,4)) 0*initConfig(:,4)];

%% %%%%%%%%%%%%%% set-up collocation method %%%%%%%%%%%%%%%%

% number of elements
nElems = 3;

% degree of polynominal = number of colloc points
nColloc = 4;

% collocation points - t is generalized time, [-1,1]
tColloc = cos((0:(nColloc-1))*pi/(nColloc-1));
tColloc = sort(tColloc);

% number of collocation points
nColloc = length(tColloc);

% matrix of powers
powerMatrix = (ones(nColloc,1)*(0:(nColloc-1)));

% matrix of colloc points
collocMatrix = (tColloc'*ones(1,nColloc));

% vandermonde matrix
vanDerMonde = collocMatrix.^powerMatrix;

% invert the thing
invVanDerMonde = inv(vanDerMonde);

% differential v-d-m matrix
diffVander = powerMatrix.*(collocMatrix.^(powerMatrix-1));

% differentiator matrix
diffMatrix = diffVander*invVanDerMonde;

% and for double differentiation
dDiffVander = powerMatrix.*(powerMatrix-1).*(collocMatrix.^(powerMatrix-2));
dDiffMatrix = dDiffVander*invVanDerMonde;

%% %% construct matrix to evaluate values at endpoints etc

% constraint points
tConstr = linspace(-1,1,5);

% number of constraint points
nConstr = length(tConstr);

% matrix of powers
conPowerMatrix = (ones(nConstr,1)*(0:(nColloc-1)));

% matrix of constraint evaluation points
constrMatrix = (tConstr'*ones(1,nColloc));

% vandermonde matrix
conVanDerMonde = constrMatrix.^conPowerMatrix;

% constraint evalution matrix
evalMatrix = conVanDerMonde*invVanDerMonde;

% derivative evaluation matrices
diffEvalMatrix = (conPowerMatrix.*[zeros(nConstr,1) conVanDerMonde(:,1:end-1)])*invVanDerMonde;
dDiffEvalMatrix = (conPowerMatrix.*(conPowerMatrix-1).*[zeros(nConstr,2) conVanDerMonde(:,1:end-2)])*invVanDerMonde;

%% %%%%%%%%%%%%%% generate initial solution %%%%%%%%%%%%%%%%

%%%% generate initial solution
weightsInit = linspace(0,1,nColloc*nElems);

% loop through each vehicle
for aa=1:nVehs,
    
    % initial solution for states
    guessConfig(1:3,(aa-1)*nColloc*nElems + (1:(nColloc*nElems))) = initConfig(aa,1:3)'*(1-weightsInit) + termConfig(aa,1:3)'*weightsInit;
    
    % initial time guess
    guessTau(aa) = (1.1*norm(initConfig(aa,1:2)-termConfig(aa,1:2))/maxSpeed)/nElems;
    
end

%% initial conflict record - empty

% empty conflict list
conflicts = false(nVehs);
sepLoss = zeros(nVehs);

% logical indicator of conflicts to enforce
avoidEnf = false(nVehs);

% element start times
guessTel = guessTau'*((1:nElems)-1);

% store for previous number of pairs enforced
Npold = 0;

% conflict resolution strategy
% confStrat = 0; % all pairs from the beginning
% confStrat = 1; % add against any conflicts found
confStrat = 2; % one at a time, worst first
% confStrat = 3; % constrain one at a time, least infringing first  

% comparison: all pairs in from scratch
if confStrat==0,
    for aa=1:nVehs,
        for a2=1:(aa-1),
            avoidEnf(aa,a2) = true;
        end
    end
end

% initialize list of pairs
[ix,jx]=find(avoidEnf);
avoidPairs = [ix jx];

%% %%%%%%%%%%%%%% write global AMPL file %%%%%%%%%%%%%%%%

%%%%% write it to AMPL file %%%%%%
c = 0;
fid=fopen('colloc.dat','w');

c = c + AMPLscalarint(fid,'Na',nVehs);

c = c + AMPLscalarint(fid,'Nz',3);

c = c + AMPLscalarint(fid,'Nc',nColloc);
c = c + AMPLscalarint(fid,'Ne',nConstr);
c = c + AMPLscalarint(fid,'Nf',nElems);

c = c + AMPLmatrix(fid,'D',diffMatrix);
c = c + AMPLmatrix(fid,'D2',dDiffMatrix);
c = c + AMPLmatrix(fid,'E',evalMatrix);
c = c + AMPLmatrix(fid,'Ed',diffEvalMatrix);
c = c + AMPLmatrix(fid,'Edd',dDiffEvalMatrix);

c = c + AMPLmatrix(fid,'zi',initZ');
c = c + AMPLmatrix(fid,'zdi',initZdot');
c = c + AMPLmatrix(fid,'zf',termZ');
c = c + AMPLmatrix(fid,'zdf',termZdot');

c = c + AMPLvector(fid,'ulo',[0.75*maxSpeed -maxTurn -maxClimb]);
c = c + AMPLvector(fid,'uhi',[maxSpeed maxTurn maxClimb]);

c = c + AMPLscalarint(fid,'Ns',nSense);
if nSense>0,
    c = c + AMPLmatrixint(fid,'sc',senseCons);
end

% completed writing file
fclose(fid);
sprintf('%d bytes written',c)

%%

% start the conflict search loop
for iter = 1:10,
    
    % re-setting Npold means you ignore the previous y values and use
    % default initial guesses
    % Npold = 0
    
    % re-initialize list of pairs
    
    % warning - doing it here could change the order and screw up the
    % initialization of the "y" values
    
    %     [ix,jx]=find(avoidEnf);
    %     avoidPairs = [ix jx];
    
    %% %%%%%%%%%%%%%%% write local AMPL file %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % element start times from length times
    for aa=1:nVehs,
        guessTel(aa,:) = guessTau(aa)*((1:nElems)-1);
    end
    
    % determine starting offsets to avoid collisions
    if 1<0,
        tOffset = zeros(nVehs,1);
        for pp=1:length(ix),
            tOffset(ix(pp)) = max(tOffset(ix(pp)),tOffset(jx(pp))+nElems*guessTau(jx(pp))+1);
        end
        for aa=1:nVehs,
            guessTel(aa,:) = tOffset(aa) + guessTau(aa)*((1:nElems)-1);
        end
    end
    
    %%%%% write it to AMPL file %%%%%%
    c = 0;
    fid=fopen('colloc2.dat','w');
    
    % initial guess
    c = c + AMPLmatrix(fid,'zguess',guessConfig);
    c = c + AMPLvector(fid,'tauguess',guessTau);
    c = c + AMPLmatrix(fid,'telguess',guessTel);
    c = c + AMPLvector(fid,'yguess',[0 0 0 -1]);
    
    % write'em
    Np = size(avoidPairs,1);
    c = c + AMPLscalarint(fid,'Np',Np);
    c = c + AMPLscalarint(fid,'Npold',Npold);
    % only write the list of pairs if any are enforced
    if Np>0,
        c = c + AMPLmatrixint(fid,'ap',avoidPairs);
    end
    
    % update for next time
    Npold = size(avoidPairs,1);
    
    % completed writing file
    fclose(fid);
    sprintf('%d bytes written',c)
    
    %% %%% solve
    
    !myampl colloc.run
    
    % stop the clock
    toc
    
    %%
    
    % load data files
    ze = load('ze.dat');
    zc = load('zc.dat');
    tau = load('tau.dat');
    %t0 = load('t0.dat');
    
    % reshape them to separate vehicles
    ze2 = reshape(ze',3,nElems*(nConstr-1)+1,nVehs);
    zc2 = reshape(zc',3,nElems*nColloc,nVehs);
    
    % figure for joint plot
    f1 = figure;
    
    % for each
    for aa=1:nVehs,
        
        % construct a time vector
        te(:,aa) = linspace(0,nElems*tau(aa),nElems*(nConstr-1)+1);
        
        %%%%% plot 3
        figure(f1)
        % path
        plot3(ze2(1,:,aa),ze2(2,:,aa),ze2(3,:,aa),'.-')
        hold on
        plot3(zc2(1,:,aa),zc2(2,:,aa),zc2(3,:,aa),'rx')
        %axis square
        %axis equal
        grid on
        plot3(initConfig(aa,1)+[0 2*cos(initConfig(aa,4))],initConfig(aa,2)+[0 2*sin(initConfig(aa,4))],initConfig(aa,3)+[0 0],'g-','LineWidth',2)
        plot3(initConfig(aa,1),initConfig(aa,2),initConfig(aa,3),'go','LineWidth',2)
        plot3(termConfig(aa,1)+[0 2*cos(termConfig(aa,4))],termConfig(aa,2)+[0 2*sin(termConfig(aa,4))],termConfig(aa,3)+[0 0],'m-','LineWidth',2)
        plot3(termConfig(aa,1),termConfig(aa,2),termConfig(aa,3),'mo','LineWidth',2)
        
        view(0,90)
        axis equal
        
        %%%%% plot
        figure
        % path
        subplot 221
        plot(ze2(1,:,aa),ze2(2,:,aa),'.-')
        hold on
        plot(zc2(1,:,aa),zc2(2,:,aa),'r:x')
        %axis equal
        grid on
        
        % vertical profile
        subplot 222
        plot(ze2(1,:,aa),ze2(3,:,aa),'.-')
        hold on
        plot(zc2(1,:,aa),zc2(3,:,aa),'r:x')
        
        % speed
        subplot 223
        plot(te(2:end,aa),sqrt(sum((diff(ze2(:,:,aa)').^2)'))*(nConstr-1)/tau(aa))
        
        
    end
    
    % add the initial guess
    figure(f1)
    plot(guessConfig(1,:),guessConfig(2,:),'k+')
    
    %% relative motion plot
    
    % max time
    tmax = max(tau)*nElems;
    
    % try a fine common timescale
    ttest = linspace(0,tmax,200);
    
    % load avoidance distance
    load R.dat
    
    % plot the cylinder
    figure
    plot3(R*cos((0:40)*pi/20),R*sin((0:40)*pi/20),ones(1,41),'r', ...
        R*cos((0:40)*pi/20),R*sin((0:40)*pi/20),-ones(1,41),'r')
    xlabel('x_1 - x_2')
    ylabel('y_1 - y_2')
    grid on
    hold on
    
    for aa=1:nVehs,
        
        % build my position on the finer timescale
        xtest(:,:,aa) = interp1(te(:,aa),ze2(:,:,aa)',ttest,'linear','extrap')';
        
        % plot relative to me
        for a2 = 1:(aa-1),
            
            % find relative positions
            xrel = xtest(:,:,aa) - xtest(:,:,a2);
            
            % plot
            plot3(xrel(1,:),xrel(2,:),xrel(3,:),'b')
            
            % identify 2-D conflicts
            horizSep = sqrt(sum(xrel(1:2,:).^2));
            in2DConf = horizSep < R - 0.1;
            minHoriz = min(horizSep);
            
            % plot'em too
            plot3(xrel(1,in2DConf),xrel(2,in2DConf),xrel(3,in2DConf),'r-s')
            
            % and record the conflicts
            if any(in2DConf),
                sepLoss(aa,a2) = R-minHoriz;
                conflicts(aa,a2) = true;
            else
                conflicts(aa,a2) = false;
                sepLoss(aa,a2) = 0;
            end
            
        end
        
    end
    
    axis equal
    view(0,90)
        
    %% check for conflicts
    if any(any(conflicts)),
        
        % sanity - check none that we've already constrained
        if any(any(conflicts & avoidEnf)),
            % problem
            error('Problem: detected conflict on a pair we have already constrained!!!')
        end
        
        if confStrat==1,
            
            % constrain against any detected
            avoidEnf = avoidEnf | conflicts
            
            % extract list of conflicting pairs
            [newix,newjx]=find(conflicts);
            
            % append it to running list of pairs, maintaining bottom zero
            avoidPairs = [avoidPairs;
                newix newjx]
            
        elseif confStrat==2,
            
            % alternative strategy - constrain the worst
            [smin,imin] = max(sepLoss);
            [smin,jmin]=max(smin);
            imin = imin(jmin);
            avoidEnf(imin,jmin)=true;
            avoidPairs = [avoidPairs;
                imin jmin]
            
        elseif confStrat==3,
            
            % alternative strategy - constrain the least first
            [smin,imin] = min(sepLoss + 20*(sepLoss==0));
            [smin,jmin]=min(smin);
            imin = imin(jmin);
            avoidEnf(imin,jmin)=true;
            avoidPairs = [avoidPairs;
                imin jmin]
            
        end
        
        % update the initial guesses
        if 1>0,
            guessConfig = zc';
            guessTau = tau;
        end
        
    else
        % finished if no more
        break
    end
    
end

% final stopwatch check
toc

return

%% animate!!

movFlag=false

figure

[xcyl,ycyl,zcyl]=cylinder;

if movFlag,
    aviobj = avifile('viewNew4.avi')
end

for kk=1:200,
    clf
    title(sprintf('tau = [%s] J=%.2f',sprintf(' %.2f', tau),sum(tau)))
    hold on
    
    for aa=1:nVehs
        
        plot3(xtest(1,1:kk,aa),xtest(2,1:kk,aa),xtest(3,1:kk,aa),'k-', ...
            xtest(1,kk,aa),xtest(2,kk,aa),xtest(3,kk,aa),'k+')
        
        surf(xtest(1,kk,aa)+.5*R*xcyl,xtest(2,kk,aa)+.5*R*ycyl,xtest(3,kk,aa)+.5*(2*zcyl-1))
        
        plot3(initConfig(aa,1)+[0 2*cos(initConfig(aa,4))],initConfig(aa,2)+[0 2*sin(initConfig(aa,4))],initConfig(aa,3)+[0 0],'g-','LineWidth',2)
        plot3(initConfig(aa,1),initConfig(aa,2),initConfig(aa,3),'go','LineWidth',2)
        plot3(termConfig(aa,1)+[0 2*cos(termConfig(aa,4))],termConfig(aa,2)+[0 2*sin(termConfig(aa,4))],termConfig(aa,3)+[0 0],'m-','LineWidth',2)
        plot3(termConfig(aa,1),termConfig(aa,2),termConfig(aa,3),'mo','LineWidth',2)
        
    end
    
    %axis vis3d
    axis equal
    axis([-5 30 -5 25 -5 4])
    view(0,90)
    
    pause(0.02)
    
    if movFlag,
        F = getframe(gcf);
        aviobj = addframe(aviobj,F);
    end
    
end

if movFlag,
    aviobj = close(aviobj)
end