% Levi Manring, Duke University
% 2021
%
% This m-file calculates the curefitted version of the ExactSolnModel
% dynamics models. The result is a curves structure that includes 
% interpolated functions for the parameters needed for the dynamic curves model.
%
% Run this function to calculate the curves prior to running the curves
% models. This function generally can take a significant amount of time to
% run (perhaps a day), depending on the refinement of the interpolation mesh.

%% Initialize

% Generate the paths for the RWD ehicle dynamics model
% addpath(genpath('/RL Car/VehicleDynamicsModelAWD'));

% load VehicleParameters if these are not loaded previously
% load('VehicleParameters'); 

%% Initialize interpolation range

% Set the interpolation range for x and mu (note a finer mesh selection is
% chosen for the region of the curve that varies the most)
xlow = (-20:0.1:-7.1)';
xmid = (-7:0.001:7)';
xhigh = (7.1:0.1:20)';
xrange = [xlow; xmid; xhigh];

muArange = (-1:0.05:1)';
muBrange = (-1:0.05:1)';

Lx = length(xrange);
LmuA = length(muArange);
LmuB = length(muBrange);

TotalSteps = Lx*LmuA*LmuB; % compute number of iterations necessary given this 3D matrix

%% Get psiB curve and psiBx curve
guess_Dx = params.dim.l;
h = params.dim.R/20;
psiB = zeros(Lx,1);
psiBx = zeros(Lx,1);
% Take numerical derivative of psiB to get psiBx
for k = 1:Lx
    xA = xrange(k,1);
    Dvec = [xA-2*h, xA-h, xA, xA+h, xA+2*h]';
    psiBsub = zeros(size(Dvec));
    for m = 1:length(Dvec)
        xD = Dvec(m,1);
        psiBsub(m,1) = psiBfun(xD,guess_Dx,params);
    end
    m = 3;
    psiB(k,1) = psiBsub(m,1);
    
    % Compute finite difference derivative
    psiBx(k,1) = (-psiBsub(m+2,1) + 8*psiBsub(m+1,1) - 8*psiBsub(m-1,1) + psiBsub(m-2,1))/(12*h);
    
    % Give progress report
    fprintf('psiB, psiBx progress: k %3.0f out of %3.0f \n',k,Lx);
end

% Fit curves for psiB and psiBx and store in curves structure
curves.psiB = fit(xrange,psiB,'linearinterp');
curves.psiBx = fit(xrange,psiBx,'linearinterp');

%% Get dynamic parameters for the noslip case

% Initialize dynamic parameters
noslip.x_coeff.Hx = zeros(Lx,1);
noslip.x_coeff.Jx = zeros(Lx,1);
noslip.x_coeff.TAx = zeros(Lx,1);
noslip.x_coeff.TBx = zeros(Lx,1);
noslip.thetaA_coeff.HthetaA = zeros(Lx,1);
noslip.thetaA_coeff.JthetaA = zeros(Lx,1);
noslip.thetaA_coeff.TAthetaA = zeros(Lx,1);
noslip.thetaA_coeff.TBthetaA = zeros(Lx,1);
noslip.thetaB_coeff.HthetaB = zeros(Lx,1);
noslip.thetaB_coeff.JthetaB = zeros(Lx,1);
noslip.thetaB_coeff.TAthetaB = zeros(Lx,1);
noslip.thetaB_coeff.TBthetaB = zeros(Lx,1);
noslip.FNA_coeff.HNA = zeros(Lx,1);
noslip.FNA_coeff.JNA = zeros(Lx,1);
noslip.FNA_coeff.TANA = zeros(Lx,1);
noslip.FNA_coeff.TBNA = zeros(Lx,1);
noslip.FNB_coeff.HNB = zeros(Lx,1);
noslip.FNB_coeff.JNB = zeros(Lx,1);
noslip.FNB_coeff.TANB = zeros(Lx,1);
noslip.FNB_coeff.TBNB = zeros(Lx,1);

% Initialize vehicle angle
tM = zeros(Lx,1);

% Calculate noslip dynamic parameters
guess_Dx = params.dim.l;
for k = 1:Lx
    x = xrange(k,1);
    [x_coeffk, thetaA_coeffk, thetaB_coeffk, FNA_coeffk, FNB_coeffk, tMk, dx] = noslip_eom(x,guess_Dx,params);
    
    tM(k,1) = tMk;
    
    % Repackage the dynamic parameters
    noslip.x_coeff.Hx(k,1) = x_coeffk.Hx;
    noslip.x_coeff.Jx(k,1) = x_coeffk.Jx;
    noslip.x_coeff.TAx(k,1) = x_coeffk.TAx;
    noslip.x_coeff.TBx(k,1) = x_coeffk.TBx;
    noslip.thetaA_coeff.HthetaA(k,1) = thetaA_coeffk.HthetaA;
    noslip.thetaA_coeff.JthetaA(k,1) = thetaA_coeffk.JthetaA;
    noslip.thetaA_coeff.TAthetaA(k,1) = thetaA_coeffk.TAthetaA;
    noslip.thetaA_coeff.TBthetaA(k,1) = thetaA_coeffk.TBthetaA;
    noslip.thetaB_coeff.HthetaB(k,1) = thetaB_coeffk.HthetaB;
    noslip.thetaB_coeff.JthetaB(k,1) = thetaB_coeffk.JthetaB;
    noslip.thetaB_coeff.TAthetaB(k,1) = thetaB_coeffk.TAthetaB;
    noslip.thetaB_coeff.TBthetaB(k,1) = thetaB_coeffk.TBthetaB;
    noslip.FNA_coeff.HNA(k,1) = FNA_coeffk.HNA;
    noslip.FNA_coeff.JNA(k,1) = FNA_coeffk.JNA;
    noslip.FNA_coeff.TANA(k,1) = FNA_coeffk.TANA;
    noslip.FNA_coeff.TBNA(k,1) = FNA_coeffk.TBNA;
    noslip.FNB_coeff.HNB(k,1) = FNB_coeffk.HNB;
    noslip.FNB_coeff.JNB(k,1) = FNB_coeffk.JNB;
    noslip.FNB_coeff.TANB(k,1) = FNB_coeffk.TANB;
    noslip.FNB_coeff.TBNB(k,1) = FNB_coeffk.TBNB;
    
    guess_Dx = dx;
    
    % Give progress report
    fprintf('noslip parameters progress: k %3.0f out of %3.0f \n',k,Lx);
end

% Fit curves for noslip dynamic parameters and store in curves structure
curves.tM                           = fit(xrange,tM,'linearinterp');
curves.noslip.x_coeff.Hx            = fit(xrange,noslip.x_coeff.Hx,'linearinterp');
curves.noslip.x_coeff.Jx            = fit(xrange,noslip.x_coeff.Jx,'linearinterp');
curves.noslip.x_coeff.TAx           = fit(xrange,noslip.x_coeff.TAx,'linearinterp');
curves.noslip.x_coeff.TBx           = fit(xrange,noslip.x_coeff.TBx,'linearinterp');
curves.noslip.thetaA_coeff.HthetaA  = fit(xrange,noslip.thetaA_coeff.HthetaA,'linearinterp');
curves.noslip.thetaA_coeff.JthetaA  = fit(xrange,noslip.thetaA_coeff.JthetaA,'linearinterp');
curves.noslip.thetaA_coeff.TAthetaA = fit(xrange,noslip.thetaA_coeff.TAthetaA,'linearinterp');
curves.noslip.thetaA_coeff.TBthetaA = fit(xrange,noslip.thetaA_coeff.TBthetaA,'linearinterp');
curves.noslip.thetaB_coeff.HthetaB  = fit(xrange,noslip.thetaB_coeff.HthetaB,'linearinterp');
curves.noslip.thetaB_coeff.JthetaB  = fit(xrange,noslip.thetaB_coeff.JthetaB,'linearinterp');
curves.noslip.thetaB_coeff.TAthetaB = fit(xrange,noslip.thetaB_coeff.TAthetaB,'linearinterp');
curves.noslip.thetaB_coeff.TBthetaB = fit(xrange,noslip.thetaB_coeff.TBthetaB,'linearinterp');
curves.noslip.FNA_coeff.HNA         = fit(xrange,noslip.FNA_coeff.HNA,'linearinterp');
curves.noslip.FNA_coeff.JNA         = fit(xrange,noslip.FNA_coeff.JNA,'linearinterp');
curves.noslip.FNA_coeff.TANA        = fit(xrange,noslip.FNA_coeff.TANA,'linearinterp');
curves.noslip.FNA_coeff.TBNA        = fit(xrange,noslip.FNA_coeff.TBNA,'linearinterp');
curves.noslip.FNB_coeff.HNB         = fit(xrange,noslip.FNB_coeff.HNB,'linearinterp');
curves.noslip.FNB_coeff.JNB         = fit(xrange,noslip.FNB_coeff.JNB,'linearinterp');
curves.noslip.FNB_coeff.TANB        = fit(xrange,noslip.FNB_coeff.TANB,'linearinterp');
curves.noslip.FNB_coeff.TBNB        = fit(xrange,noslip.FNB_coeff.TBNB,'linearinterp');

%% Get dynamic parameters for the rearstickfrontslip case

% Initialize dynamic parameters
rearstickfrontslip.x_coeff.Hx = zeros(Lx,LmuB);
rearstickfrontslip.x_coeff.Jx = zeros(Lx,LmuB);
rearstickfrontslip.x_coeff.TAx = zeros(Lx,LmuB);
rearstickfrontslip.x_coeff.TBx = zeros(Lx,LmuB);
rearstickfrontslip.thetaA_coeff.HthetaA = zeros(Lx,LmuB);
rearstickfrontslip.thetaA_coeff.JthetaA = zeros(Lx,LmuB);
rearstickfrontslip.thetaA_coeff.TAthetaA = zeros(Lx,LmuB);
rearstickfrontslip.thetaA_coeff.TBthetaA = zeros(Lx,LmuB);
rearstickfrontslip.thetaB_coeff.HthetaB = zeros(Lx,LmuB);
rearstickfrontslip.thetaB_coeff.JthetaB = zeros(Lx,LmuB);
rearstickfrontslip.thetaB_coeff.TAthetaB = zeros(Lx,LmuB);
rearstickfrontslip.thetaB_coeff.TBthetaB = zeros(Lx,LmuB);
rearstickfrontslip.FNA_coeff.HNA = zeros(Lx,LmuB);
rearstickfrontslip.FNA_coeff.JNA = zeros(Lx,LmuB);
rearstickfrontslip.FNA_coeff.TANA = zeros(Lx,LmuB);
rearstickfrontslip.FNA_coeff.TBNA = zeros(Lx,LmuB);
rearstickfrontslip.FNB_coeff.HNB = zeros(Lx,LmuB);
rearstickfrontslip.FNB_coeff.JNB = zeros(Lx,LmuB);
rearstickfrontslip.FNB_coeff.TANB = zeros(Lx,LmuB);
rearstickfrontslip.FNB_coeff.TBNB = zeros(Lx,LmuB);

% Calculate rearstickfrontslip dynamic parameters
guess_Dx = params.dim.l;
for k = 1:Lx
    x = xrange(k,1);
    for m = 1:LmuB
        muB = muBrange(m,1);
        [x_coeffk, thetaA_coeffk, thetaB_coeffk, FNA_coeffk, FNB_coeffk, tMk, dx] = rearstickfrontslip_eom(x,guess_Dx,muB,params);
        
        % Repackage the dynamic parameters
        rearstickfrontslip.x_coeff.Hx(k,m) = x_coeffk.Hx;
        rearstickfrontslip.x_coeff.Jx(k,m) = x_coeffk.Jx;
        rearstickfrontslip.x_coeff.TAx(k,m) = x_coeffk.TAx;
        rearstickfrontslip.x_coeff.TBx(k,m) = x_coeffk.TBx;
        rearstickfrontslip.thetaA_coeff.HthetaA(k,m) = thetaA_coeffk.HthetaA;
        rearstickfrontslip.thetaA_coeff.JthetaA(k,m) = thetaA_coeffk.JthetaA;
        rearstickfrontslip.thetaA_coeff.TAthetaA(k,m) = thetaA_coeffk.TAthetaA;
        rearstickfrontslip.thetaA_coeff.TBthetaA(k,m) = thetaA_coeffk.TBthetaA;
        rearstickfrontslip.thetaB_coeff.HthetaB(k,m) = thetaB_coeffk.HthetaB;
        rearstickfrontslip.thetaB_coeff.JthetaB(k,m) = thetaB_coeffk.JthetaB;
        rearstickfrontslip.thetaB_coeff.TAthetaB(k,m) = thetaB_coeffk.TAthetaB;
        rearstickfrontslip.thetaB_coeff.TBthetaB(k,m) = thetaB_coeffk.TBthetaB;
        rearstickfrontslip.FNA_coeff.HNA(k,m) = FNA_coeffk.HNA;
        rearstickfrontslip.FNA_coeff.JNA(k,m) = FNA_coeffk.JNA;
        rearstickfrontslip.FNA_coeff.TANA(k,m) = FNA_coeffk.TANA;
        rearstickfrontslip.FNA_coeff.TBNA(k,m) = FNA_coeffk.TBNA;
        rearstickfrontslip.FNB_coeff.HNB(k,m) = FNB_coeffk.HNB;
        rearstickfrontslip.FNB_coeff.JNB(k,m) = FNB_coeffk.JNB;
        rearstickfrontslip.FNB_coeff.TANB(k,m) = FNB_coeffk.TANB;
        rearstickfrontslip.FNB_coeff.TBNB(k,m) = FNB_coeffk.TBNB;
    end
    guess_Dx = dx;
    
    % Give progress report
    fprintf('rearstickfrontslip parameters progress: k %3.0f out of %3.0f \n',k,Lx);
end

% Fit curves for rearstickfrontslip dynamic parameters and store in curves structure
[muBmesh,xmesh] = meshgrid(muBrange,xrange);
curves.rearstickfrontslip.x_coeff.Hx            = fit([xmesh(:), muBmesh(:)],rearstickfrontslip.x_coeff.Hx(:),'linearinterp'); fprintf('1 done\n');
curves.rearstickfrontslip.x_coeff.Jx            = fit([xmesh(:), muBmesh(:)],rearstickfrontslip.x_coeff.Jx(:),'linearinterp'); fprintf('1 done\n');
curves.rearstickfrontslip.x_coeff.TAx           = fit([xmesh(:), muBmesh(:)],rearstickfrontslip.x_coeff.TAx(:),'linearinterp'); fprintf('1 done\n');
curves.rearstickfrontslip.x_coeff.TBx           = fit([xmesh(:), muBmesh(:)],rearstickfrontslip.x_coeff.TBx(:),'linearinterp'); fprintf('1 done\n');
curves.rearstickfrontslip.thetaA_coeff.HthetaA  = fit([xmesh(:), muBmesh(:)],rearstickfrontslip.thetaA_coeff.HthetaA(:),'linearinterp'); fprintf('1 done\n');
curves.rearstickfrontslip.thetaA_coeff.JthetaA  = fit([xmesh(:), muBmesh(:)],rearstickfrontslip.thetaA_coeff.JthetaA(:),'linearinterp'); fprintf('1 done\n');
curves.rearstickfrontslip.thetaA_coeff.TAthetaA = fit([xmesh(:), muBmesh(:)],rearstickfrontslip.thetaA_coeff.TAthetaA(:),'linearinterp'); fprintf('1 done\n');
curves.rearstickfrontslip.thetaA_coeff.TBthetaA = fit([xmesh(:), muBmesh(:)],rearstickfrontslip.thetaA_coeff.TBthetaA(:),'linearinterp'); fprintf('1 done\n');
curves.rearstickfrontslip.thetaB_coeff.HthetaB  = fit([xmesh(:), muBmesh(:)],rearstickfrontslip.thetaB_coeff.HthetaB(:),'linearinterp'); fprintf('1 done\n');
curves.rearstickfrontslip.thetaB_coeff.JthetaB  = fit([xmesh(:), muBmesh(:)],rearstickfrontslip.thetaB_coeff.JthetaB(:),'linearinterp'); fprintf('1 done\n');
curves.rearstickfrontslip.thetaB_coeff.TAthetaB = fit([xmesh(:), muBmesh(:)],rearstickfrontslip.thetaB_coeff.TAthetaB(:),'linearinterp'); fprintf('1 done\n');
curves.rearstickfrontslip.thetaB_coeff.TBthetaB = fit([xmesh(:), muBmesh(:)],rearstickfrontslip.thetaB_coeff.TBthetaB(:),'linearinterp'); fprintf('1 done\n');
curves.rearstickfrontslip.FNA_coeff.HNA         = fit([xmesh(:), muBmesh(:)],rearstickfrontslip.FNA_coeff.HNA(:),'linearinterp'); fprintf('1 done\n');
curves.rearstickfrontslip.FNA_coeff.JNA         = fit([xmesh(:), muBmesh(:)],rearstickfrontslip.FNA_coeff.JNA(:),'linearinterp'); fprintf('1 done\n');
curves.rearstickfrontslip.FNA_coeff.TANA        = fit([xmesh(:), muBmesh(:)],rearstickfrontslip.FNA_coeff.TANA(:),'linearinterp'); fprintf('1 done\n');
curves.rearstickfrontslip.FNA_coeff.TBNA        = fit([xmesh(:), muBmesh(:)],rearstickfrontslip.FNA_coeff.TBNA(:),'linearinterp'); fprintf('1 done\n');
curves.rearstickfrontslip.FNB_coeff.HNB         = fit([xmesh(:), muBmesh(:)],rearstickfrontslip.FNB_coeff.HNB(:),'linearinterp'); fprintf('1 done\n');
curves.rearstickfrontslip.FNB_coeff.JNB         = fit([xmesh(:), muBmesh(:)],rearstickfrontslip.FNB_coeff.JNB(:),'linearinterp'); fprintf('1 done\n');
curves.rearstickfrontslip.FNB_coeff.TANB        = fit([xmesh(:), muBmesh(:)],rearstickfrontslip.FNB_coeff.TANB(:),'linearinterp'); fprintf('1 done\n');
curves.rearstickfrontslip.FNB_coeff.TBNB        = fit([xmesh(:), muBmesh(:)],rearstickfrontslip.FNB_coeff.TBNB(:),'linearinterp'); fprintf('1 done\n');

%% Get dynamic parameters for the rearslipfrontstick case

% Initialize dynamic parameters
rearslipfrontstick.x_coeff.Hx = zeros(Lx,LmuA);
rearslipfrontstick.x_coeff.Jx = zeros(Lx,LmuA);
rearslipfrontstick.x_coeff.TAx = zeros(Lx,LmuA);
rearslipfrontstick.x_coeff.TBx = zeros(Lx,LmuA);
rearslipfrontstick.thetaA_coeff.HthetaA = zeros(Lx,LmuA);
rearslipfrontstick.thetaA_coeff.JthetaA = zeros(Lx,LmuA);
rearslipfrontstick.thetaA_coeff.TAthetaA = zeros(Lx,LmuA);
rearslipfrontstick.thetaA_coeff.TBthetaA = zeros(Lx,LmuA);
rearslipfrontstick.thetaB_coeff.HthetaB = zeros(Lx,LmuA);
rearslipfrontstick.thetaB_coeff.JthetaB = zeros(Lx,LmuA);
rearslipfrontstick.thetaB_coeff.TAthetaB = zeros(Lx,LmuA);
rearslipfrontstick.thetaB_coeff.TBthetaB = zeros(Lx,LmuA);
rearslipfrontstick.FNA_coeff.HNA = zeros(Lx,LmuA);
rearslipfrontstick.FNA_coeff.JNA = zeros(Lx,LmuA);
rearslipfrontstick.FNA_coeff.TANA = zeros(Lx,LmuA);
rearslipfrontstick.FNA_coeff.TBNA = zeros(Lx,LmuA);
rearslipfrontstick.FNB_coeff.HNB = zeros(Lx,LmuA);
rearslipfrontstick.FNB_coeff.JNB = zeros(Lx,LmuA);
rearslipfrontstick.FNB_coeff.TANB = zeros(Lx,LmuA);
rearslipfrontstick.FNB_coeff.TBNB = zeros(Lx,LmuA);
 
% Calculate rearslipfrontstick dynamic parameters
guess_Dx = params.dim.l;
for k = 1:Lx
    x = xrange(k,1);
    for m = 1:LmuA
        muA = muArange(m,1);
        [x_coeffk, thetaA_coeffk, thetaB_coeffk, FNA_coeffk, FNB_coeffk, tMk, dx] = rearslipfrontstick_eom(x,guess_Dx,muA,params);
        
        % Repackage the dynamic parameters
        rearslipfrontstick.x_coeff.Hx(k,m) = x_coeffk.Hx;
        rearslipfrontstick.x_coeff.Jx(k,m) = x_coeffk.Jx;
        rearslipfrontstick.x_coeff.TAx(k,m) = x_coeffk.TAx;
        rearslipfrontstick.x_coeff.TBx(k,m) = x_coeffk.TBx;
        rearslipfrontstick.thetaA_coeff.HthetaA(k,m) = thetaA_coeffk.HthetaA;
        rearslipfrontstick.thetaA_coeff.JthetaA(k,m) = thetaA_coeffk.JthetaA;
        rearslipfrontstick.thetaA_coeff.TAthetaA(k,m) = thetaA_coeffk.TAthetaA;
        rearslipfrontstick.thetaA_coeff.TBthetaA(k,m) = thetaA_coeffk.TBthetaA;
        rearslipfrontstick.thetaB_coeff.HthetaB(k,m) = thetaB_coeffk.HthetaB;
        rearslipfrontstick.thetaB_coeff.JthetaB(k,m) = thetaB_coeffk.JthetaB;
        rearslipfrontstick.thetaB_coeff.TAthetaB(k,m) = thetaB_coeffk.TAthetaB;
        rearslipfrontstick.thetaB_coeff.TBthetaB(k,m) = thetaB_coeffk.TBthetaB;
        rearslipfrontstick.FNA_coeff.HNA(k,m) = FNA_coeffk.HNA;
        rearslipfrontstick.FNA_coeff.JNA(k,m) = FNA_coeffk.JNA;
        rearslipfrontstick.FNA_coeff.TANA(k,m) = FNA_coeffk.TANA;
        rearslipfrontstick.FNA_coeff.TBNA(k,m) = FNA_coeffk.TBNA;
        rearslipfrontstick.FNB_coeff.HNB(k,m) = FNB_coeffk.HNB;
        rearslipfrontstick.FNB_coeff.JNB(k,m) = FNB_coeffk.JNB;
        rearslipfrontstick.FNB_coeff.TANB(k,m) = FNB_coeffk.TANB;
        rearslipfrontstick.FNB_coeff.TBNB(k,m) = FNB_coeffk.TBNB;
    end
    guess_Dx = dx;
    
    % Give progress report
    fprintf('rearslipfrontstick parameters progress: k %3.0f out of %3.0f \n',k,Lx);
end

% Fit curves for rearslipfrontstick dynamic parameters and store in curves structure
[muAmesh,xmesh] = meshgrid(muArange,xrange);
curves.rearslipfrontstick.x_coeff.Hx            = fit([xmesh(:), muAmesh(:)],rearslipfrontstick.x_coeff.Hx(:),'linearinterp'); fprintf('1 done\n');
curves.rearslipfrontstick.x_coeff.Jx            = fit([xmesh(:), muAmesh(:)],rearslipfrontstick.x_coeff.Jx(:),'linearinterp'); fprintf('1 done\n');
curves.rearslipfrontstick.x_coeff.TAx           = fit([xmesh(:), muAmesh(:)],rearslipfrontstick.x_coeff.TAx(:),'linearinterp'); fprintf('1 done\n');
curves.rearslipfrontstick.x_coeff.TBx           = fit([xmesh(:), muAmesh(:)],rearslipfrontstick.x_coeff.TBx(:),'linearinterp'); fprintf('1 done\n');
curves.rearslipfrontstick.thetaA_coeff.HthetaA  = fit([xmesh(:), muAmesh(:)],rearslipfrontstick.thetaA_coeff.HthetaA(:),'linearinterp'); fprintf('1 done\n');
curves.rearslipfrontstick.thetaA_coeff.JthetaA  = fit([xmesh(:), muAmesh(:)],rearslipfrontstick.thetaA_coeff.JthetaA(:),'linearinterp'); fprintf('1 done\n');
curves.rearslipfrontstick.thetaA_coeff.TAthetaA = fit([xmesh(:), muAmesh(:)],rearslipfrontstick.thetaA_coeff.TAthetaA(:),'linearinterp'); fprintf('1 done\n');
curves.rearslipfrontstick.thetaA_coeff.TBthetaA = fit([xmesh(:), muAmesh(:)],rearslipfrontstick.thetaA_coeff.TBthetaA(:),'linearinterp'); fprintf('1 done\n');
curves.rearslipfrontstick.thetaB_coeff.HthetaB  = fit([xmesh(:), muAmesh(:)],rearslipfrontstick.thetaB_coeff.HthetaB(:),'linearinterp'); fprintf('1 done\n');
curves.rearslipfrontstick.thetaB_coeff.JthetaB  = fit([xmesh(:), muAmesh(:)],rearslipfrontstick.thetaB_coeff.JthetaB(:),'linearinterp'); fprintf('1 done\n');
curves.rearslipfrontstick.thetaB_coeff.TAthetaB = fit([xmesh(:), muAmesh(:)],rearslipfrontstick.thetaB_coeff.TAthetaB(:),'linearinterp'); fprintf('1 done\n');
curves.rearslipfrontstick.thetaB_coeff.TBthetaB = fit([xmesh(:), muAmesh(:)],rearslipfrontstick.thetaB_coeff.TBthetaB(:),'linearinterp'); fprintf('1 done\n');
curves.rearslipfrontstick.FNA_coeff.HNA         = fit([xmesh(:), muAmesh(:)],rearslipfrontstick.FNA_coeff.HNA(:),'linearinterp'); fprintf('1 done\n');
curves.rearslipfrontstick.FNA_coeff.JNA         = fit([xmesh(:), muAmesh(:)],rearslipfrontstick.FNA_coeff.JNA(:),'linearinterp'); fprintf('1 done\n');
curves.rearslipfrontstick.FNA_coeff.TANA        = fit([xmesh(:), muAmesh(:)],rearslipfrontstick.FNA_coeff.TANA(:),'linearinterp'); fprintf('1 done\n');
curves.rearslipfrontstick.FNA_coeff.TBNA        = fit([xmesh(:), muAmesh(:)],rearslipfrontstick.FNA_coeff.TBNA(:),'linearinterp'); fprintf('1 done\n');
curves.rearslipfrontstick.FNB_coeff.HNB         = fit([xmesh(:), muAmesh(:)],rearslipfrontstick.FNB_coeff.HNB(:),'linearinterp'); fprintf('1 done\n');
curves.rearslipfrontstick.FNB_coeff.JNB         = fit([xmesh(:), muAmesh(:)],rearslipfrontstick.FNB_coeff.JNB(:),'linearinterp'); fprintf('1 done\n');
curves.rearslipfrontstick.FNB_coeff.TANB        = fit([xmesh(:), muAmesh(:)],rearslipfrontstick.FNB_coeff.TANB(:),'linearinterp'); fprintf('1 done\n');
curves.rearslipfrontstick.FNB_coeff.TBNB        = fit([xmesh(:), muAmesh(:)],rearslipfrontstick.FNB_coeff.TBNB(:),'linearinterp'); fprintf('1 done\n');

%% Choose to save current curves structure in case of memory failure
save('CurvesMostModels','curves','-v7.3');

%% Get dynamic parameters for the allslip case
% NOTE: this section requires a lot of memory. Possibility of memory failure

% Initialize dynamic parameters
allslip.x_coeff.Hx = zeros(Lx,LmuA,LmuB);
allslip.x_coeff.Jx = zeros(Lx,LmuA,LmuB);
allslip.x_coeff.TAx = zeros(Lx,LmuA,LmuB);
allslip.x_coeff.TBx = zeros(Lx,LmuA,LmuB);
allslip.thetaA_coeff.HthetaA = zeros(Lx,LmuA,LmuB);
allslip.thetaA_coeff.JthetaA = zeros(Lx,LmuA,LmuB);
allslip.thetaA_coeff.TAthetaA = zeros(Lx,LmuA,LmuB);
allslip.thetaA_coeff.TBthetaA = zeros(Lx,LmuA,LmuB);
allslip.thetaB_coeff.HthetaB = zeros(Lx,LmuA,LmuB);
allslip.thetaB_coeff.JthetaB = zeros(Lx,LmuA,LmuB);
allslip.thetaB_coeff.TAthetaB = zeros(Lx,LmuA,LmuB);
allslip.thetaB_coeff.TBthetaB = zeros(Lx,LmuA,LmuB);
allslip.FNA_coeff.HNA = zeros(Lx,LmuA,LmuB);
allslip.FNA_coeff.JNA = zeros(Lx,LmuA,LmuB);
allslip.FNA_coeff.TANA = zeros(Lx,LmuA,LmuB);
allslip.FNA_coeff.TBNA = zeros(Lx,LmuA,LmuB);
allslip.FNB_coeff.HNB = zeros(Lx,LmuA,LmuB);
allslip.FNB_coeff.JNB = zeros(Lx,LmuA,LmuB);
allslip.FNB_coeff.TANB = zeros(Lx,LmuA,LmuB);
allslip.FNB_coeff.TBNB = zeros(Lx,LmuA,LmuB);

% Calculate allslip dynamic parameters
guess_Dx = params.dim.l;
for k = 1:Lx
    x = xrange(k,1);
    for m = 1:LmuA
        muA = muArange(m,1);
        for n = 1:LmuB
            muB = muBrange(n,1);
            [x_coeffk, thetaA_coeffk, thetaB_coeffk, FNA_coeffk, FNB_coeffk, tMk, dx] = allslip_eom(x,guess_Dx,muA,muB,params);
            
            % Repackage the dynamic parameters
            allslip.x_coeff.Hx(k,m,n) = x_coeffk.Hx;
            allslip.x_coeff.Jx(k,m,n) = x_coeffk.Jx;
            allslip.x_coeff.TAx(k,m,n) = x_coeffk.TAx;
            allslip.x_coeff.TBx(k,m,n) = x_coeffk.TBx;
            allslip.thetaA_coeff.HthetaA(k,m,n) = thetaA_coeffk.HthetaA;
            allslip.thetaA_coeff.JthetaA(k,m,n) = thetaA_coeffk.JthetaA;
            allslip.thetaA_coeff.TAthetaA(k,m,n) = thetaA_coeffk.TAthetaA;
            allslip.thetaA_coeff.TBthetaA(k,m,n) = thetaA_coeffk.TBthetaA;
            allslip.thetaB_coeff.HthetaB(k,m,n) = thetaB_coeffk.HthetaB;
            allslip.thetaB_coeff.JthetaB(k,m,n) = thetaB_coeffk.JthetaB;
            allslip.thetaB_coeff.TAthetaB(k,m,n) = thetaB_coeffk.TAthetaB;
            allslip.thetaB_coeff.TBthetaB(k,m,n) = thetaB_coeffk.TBthetaB;
            allslip.FNA_coeff.HNA(k,m,n) = FNA_coeffk.HNA;
            allslip.FNA_coeff.JNA(k,m,n) = FNA_coeffk.JNA;
            allslip.FNA_coeff.TANA(k,m,n) = FNA_coeffk.TANA;
            allslip.FNA_coeff.TBNA(k,m,n) = FNA_coeffk.TBNA;
            allslip.FNB_coeff.HNB(k,m,n) = FNB_coeffk.HNB;
            allslip.FNB_coeff.JNB(k,m,n) = FNB_coeffk.JNB;
            allslip.FNB_coeff.TANB(k,m,n) = FNB_coeffk.TANB;
            allslip.FNB_coeff.TBNB(k,m,n) = FNB_coeffk.TBNB;
        end
    end
    guess_Dx = dx;
    
    % Give progress report
    fprintf('allslip parameters progress: k %3.0f out of %3.0f \n',k,Lx);
end

% Fit curves for rearslipfrontstick dynamic parameters and store in curves structure
curves.allslip.x_coeff.Hx            = @(x,muA,muB) interpn(xrange,muArange,muBrange,allslip.x_coeff.Hx,x,muA,muB,'linear'); fprintf('1 done\n');
curves.allslip.x_coeff.Jx            = @(x,muA,muB) interpn(xrange,muArange,muBrange,allslip.x_coeff.Jx,x,muA,muB,'linear'); fprintf('1 done\n');
curves.allslip.x_coeff.TAx           = @(x,muA,muB) interpn(xrange,muArange,muBrange,allslip.x_coeff.TAx,x,muA,muB,'linear'); fprintf('1 done\n');
curves.allslip.x_coeff.TBx           = @(x,muA,muB) interpn(xrange,muArange,muBrange,allslip.x_coeff.TBx,x,muA,muB,'linear'); fprintf('1 done\n');
curves.allslip.thetaA_coeff.HthetaA  = @(x,muA,muB) interpn(xrange,muArange,muBrange,allslip.thetaA_coeff.HthetaA,x,muA,muB,'linear'); fprintf('1 done\n');
curves.allslip.thetaA_coeff.JthetaA  = @(x,muA,muB) interpn(xrange,muArange,muBrange,allslip.thetaA_coeff.JthetaA,x,muA,muB,'linear'); fprintf('1 done\n');
curves.allslip.thetaA_coeff.TAthetaA = @(x,muA,muB) interpn(xrange,muArange,muBrange,allslip.thetaA_coeff.TAthetaA,x,muA,muB,'linear'); fprintf('1 done\n');
curves.allslip.thetaA_coeff.TBthetaA = @(x,muA,muB) interpn(xrange,muArange,muBrange,allslip.thetaA_coeff.TBthetaA,x,muA,muB,'linear'); fprintf('1 done\n');
curves.allslip.thetaB_coeff.HthetaB  = @(x,muA,muB) interpn(xrange,muArange,muBrange,allslip.thetaB_coeff.HthetaB,x,muA,muB,'linear'); fprintf('1 done\n');
curves.allslip.thetaB_coeff.JthetaB  = @(x,muA,muB) interpn(xrange,muArange,muBrange,allslip.thetaB_coeff.JthetaB,x,muA,muB,'linear'); fprintf('1 done\n');
curves.allslip.thetaB_coeff.TAthetaB = @(x,muA,muB) interpn(xrange,muArange,muBrange,allslip.thetaB_coeff.TAthetaB,x,muA,muB,'linear'); fprintf('1 done\n');
curves.allslip.thetaB_coeff.TBthetaB = @(x,muA,muB) interpn(xrange,muArange,muBrange,allslip.thetaB_coeff.TBthetaB,x,muA,muB,'linear'); fprintf('1 done\n');
curves.allslip.FNA_coeff.HNA         = @(x,muA,muB) interpn(xrange,muArange,muBrange,allslip.FNA_coeff.HNA,x,muA,muB,'linear'); fprintf('1 done\n');
curves.allslip.FNA_coeff.JNA         = @(x,muA,muB) interpn(xrange,muArange,muBrange,allslip.FNA_coeff.JNA,x,muA,muB,'linear'); fprintf('1 done\n');
curves.allslip.FNA_coeff.TANA        = @(x,muA,muB) interpn(xrange,muArange,muBrange,allslip.FNA_coeff.TANA,x,muA,muB,'linear'); fprintf('1 done\n');
curves.allslip.FNA_coeff.TBNA        = @(x,muA,muB) interpn(xrange,muArange,muBrange,allslip.FNA_coeff.TBNA,x,muA,muB,'linear'); fprintf('1 done\n');
curves.allslip.FNB_coeff.HNB         = @(x,muA,muB) interpn(xrange,muArange,muBrange,allslip.FNB_coeff.HNB,x,muA,muB,'linear'); fprintf('1 done\n');
curves.allslip.FNB_coeff.JNB         = @(x,muA,muB) interpn(xrange,muArange,muBrange,allslip.FNB_coeff.JNB,x,muA,muB,'linear'); fprintf('1 done\n');
curves.allslip.FNB_coeff.TANB        = @(x,muA,muB) interpn(xrange,muArange,muBrange,allslip.FNB_coeff.TANB,x,muA,muB,'linear'); fprintf('1 done\n');
curves.allslip.FNB_coeff.TBNB        = @(x,muA,muB) interpn(xrange,muArange,muBrange,allslip.FNB_coeff.TBNB,x,muA,muB,'linear'); fprintf('1 done\n');

save('CurvesAllModels','curves','-v7.3');
