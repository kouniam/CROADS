% CROADS 
% this is the carbon and temperature part of C-ROADS. 
% emissions are prescribed by the user. hint: currently, emissions are
% around 10Gt(C) [Gigaton of carbon] per year. 

%

close all; %kill all previous plots
clear all; %delete all previous variables

%Tfinal = zeros(1,500);
%histo = zeros(1,4000);
% LAMBDA = zeros(1,500);

for l=1:1:100

%-----------------------------------------------------------------------%

Emissource = 3;     %emission scenario 1: specified; 2: historical; 3: an RCP scenario
RCPscenario = 85;   %26, 45, 60, 85 stand for RCP2.6, RCP4.5, RCP6.0, RCP8.5

%-----------------------------------------------------------------------%

if Emissource == 1 %specify emissions
    
    dt = 1;                 % time step (in years) 
    yr1 = 1;                % first year
    yr2 = 10000;            % last year 
    timvec = yr1:dt:yr2;    % times
    nt = length(timvec);    % nr of time steps (which are years)

    Emis = zeros(nt,1); % all anthrop. emissions (landuse change included), in Gt(C)/year
    % Emis(2) = 800; Emis(2000)=800; %Experiment "Pulse"
    % Emis(2:500) = 8;

elseif Emissource == 2 %reading in the emissions 
    
    load('Carbonemission_historical.mat'); 
    dt = tvec(2)-tvec(1);               % adjust to data we read in... (tvec is times we read in)
    yr1 = tvec(1); 
    yr2 = tvec(1)+1000;
    timvec = yr1:dt:yr2; 
    nt = length(timvec);                % timvec is the times we eventually use, starting in same year as tvec
    Emis = zeros(nt,1);                 % all anthrop. emissions (landuse change included), in Gt(C)/year
    Emis(1:length(Cemhist)) = Cemhist;  %loaded from the file

elseif Emissource == 3 % reading in the emissions 
    
    % here we first read in the historical data, then the corresponding
    % scenario
    
load('Carbonemission_historical.mat'); 
dt = tvec(2)-tvec(1);               % adjust to data we read in... (tvec is times we read in)
yr1 = tvec(1); 
yr2 = tvec(1)+2000;
timvec = yr1:dt:yr2; 
nt = length(timvec);                % timvec is the times we eventually use, starting in same year as tvec
Emis = zeros(nt,1);                 % all anthrop. emissions (landuse change included), in Gt(C)/year
Emis(1:length(Cemhist)) = Cemhist;  % loaded from the file

load('Carbonemission_RCP_scenarios.mat'); 
% now we expand Emis with the RCP scenario data
yr3 = tvec(1);                              % first year of the RCP scenario
findyr = find(timvec == yr3); 
varname = ['Cem' num2str(RCPscenario)];     % name of the emission time series used

Emis(findyr:findyr+length(tvec)-1) = eval(varname); 

end % end if Emissource 

Emis(1) = 0; % for consistency, since the model only really starts in step 2

%--------------------------- model parameters ---------------------------%

% OCEAN LAYERS
laydep = [100,300,300,1300,1800];   % depth of ocean layers (meter); the first is the mixed layer.
nl = length(laydep);                % nr of layers

% PRE-INDUSTRIAL CARBON 
Cat0 = 590;             % pre-ind. atmospheric carbon [Gt(C)]]
Con00 = 10.2373;        % pre-ind. oceanic carbon per meter layer depth
Con0 = Con00*laydep;    % vector of pre-ind. carbon in each layer

% PARAMETERS FOR OCEANIC PROCESSES
rev0 = 9.7;     % standard revelle or buffer factor 
revC = 3.92;    % impact C conc. on revelle factor 
Cm0 = Con0(1);  % Gt(C)/m ocean depth  * depth mixed-layer (reference ocean carbon of mixed layer)
CmT = 0.003;    % K^-1 (dependence of Cm on temperature)
eddydif = 4400; % m^2/year  eddy diffusion coefficient (ocean mixing between layers)

% BIOSPHERE PROCESSES
NPP0 = 85.1771;     % reference NPP (Net Primary Production, i.e. photosynthesis - plant respiration) in Gt(C)/year.  
humfrac = 0.428;    % humification fraction (how much of rotting plants goes into humus. rest goes as CO2 into atm.)
biotime = 10.6;     % residence time Carbon in biosphere (years) (Croads guide 2016, tab 3-32)
humtime = 27.8;     % residence time Carbon in humus (years)
fert = 0.42;        % carbon fertilisation effect (NPP increase with atm. CO2 conc.)
heatstress = -0.01; % reduction NPP by heat; K^-1;

% values for no NPP dependency on carbon concentration and temperature:
%fert = 0;
%heatstress = 0;

% PRE-INDUSTRIAL BIOSPHERE CARBON RESERVOIRS
% Values choosen such that equilibrium will be achieved for pre-industrial conditions. 
biom0 = NPP0*biotime;                   % pre-industrial carbon stock in biosphere [Gt(C)];              
humm0 = humfrac*biom0*humtime/biotime;  % pre-industrial carbon stock in humus Gt(C); 
  
% RADIATIVE FORCING AND HEAT
forCO2 = 5.35*1.0;          % W/m^2 (rad. forcing for log(Cat/Cat0)=1, i.e. for e-folding CO2 )
outrad = 1.23;              % W/m^2/K (outgoing rad. for 1K surface warming; 1.23 -> climate sens.=3K / doubling CO2)
seasurf = 0.708;            % fraction of earth's surface that is sea
heatcap = 4.23e6;           % J/K m^3. heat capacity of water
secyr = 365.25*24*3600;     % seconds per year 


%--------------------------- define variables ---------------------------%

Cat = zeros(nt,1);      % total atmospheric carbon , in Gt(C)
Con = zeros(nt,nl);     % oceanic carbon, in Gt(C) 
dCat = Cat;             % corresponding time derivative

Hon = Con;  % heat content per m^2 of sea surface
Ton = Con;  % temperature of ocean layers

NPP = Cat;    % Net Primary Production by plants/biosphere
biom = Cat;   % biomass
humm = Cat;   % humus mass

humrelease = Cat; biorelease = Cat; % fluxes from biomass / humus into atm

oceanuptake = Cat;  % flux from atm to mixed layer
revelle = Cat;      % revelle factor (depends on Con... and T?)

T = Cat;        % global surface warming in K
Tequil = Cat;   % surf. temp. that would be in equilibrium with the current atm. carbon. not needed, but nice to plot.
FCO2 = Cat;     % forcing from CO2
Fout = Cat;     % longwave forcing due to warming


%----------------- initialisation, external time series -----------------%
% Cat, Con, humm, biom, T   (NPP... not needed)
% Emissions are specified on the top

CumEmis = Emis; %cumulated emissions (just for comparison

for t=2:nt;  
    CumEmis(t) = CumEmis(t-1)+Emis(t); 
end

% now, the initialisation. Here we start from pre-industrial values defined above.
Cat(1) = Cat0;
Con(1,:) = Con0;
biom(1) = biom0; 
humm(1) = humm0;

Hon(1,:) =0;    % heat content / m^2 ocean surface in each layer
Ton(1,:) =0;    % temp of the layers
T(1) = 0;       % surf. temp

% Randomly select a value for the climate sensitivity based on the PDF

r = rand;

if (r<0.05)
    lambda = 2;
    outrad = 1.8542;
elseif (0.05<=r)&&(r<0.15)
    lambda = 2.5;
    outrad = 1.4833;
elseif (0.15<=r)&&(r<0.40)
    lambda = 3;
    outrad = 1.23;
elseif (0.40<=r)&&(r<0.65)
    lambda = 3.5;
    outrad = 1.0595;
elseif (0.65<=r)&&(r<0.85)
    lambda = 4;
    outrad = 0.9271;
elseif (0.85<=r)&&(r<0.95)
    lambda = 4.5;
    outrad = 0.8241;
else
    lambda = 5;
    outrad = 0.7417;
end

% LAMBDA(l) = lambda;

%------------------------------------------------------------------------%
%------------------------------- THE MODEL ------------------------------%
%------------------------------------------------------------------------%

% REMARKS ON THE CARBON PART OF C-ROADS
% first, compute biosphere and ocean mixing processes (both in an
% euler-foreward fashion, based on values from last time step). 

% then, we know the total carbon in atmosphere and mixedlayer, and
% equilibrate those using iterative procedure. 

for t=2:nt
    
%--------------------------- CARBON MODELLING ---------------------------%  

    % I first check how much carbon is taken into/out of the atm. by
    % bio-processes plus human emissions. 
    % Next, I check how much carbon is taken into/out of atm. by mixing with the
    % deep ocean.
    % This tells me how much carbon is in the combined system of atmosphere and
    % ocean mixed layer (namely Cat1 + Con1)
    % Finally, I distribute this amount between atm. and mixedlayer such that
    % the equilibrium is achieved. In other words, the total amount of carbon in
    % atm + mixed layer is still Cat1+Con1, and in addition, eq. 2 and 3 of
    % Sterman 2013 hold. 
    % I do that by means of an iteration procedure, which may be slightly overkill but at least it works decently. 

    % BIOSPHERE PROCESSES
    NPP(t) = NPP0*(1+fert*log(Cat(t-1)/Cat0))*(1+heatstress*T(t-1))*dt; %atm CO2 taken up by plants
    humrelease(t) = humm(t-1)/humtime*dt;                               %CO2 released into atm. by decaying humus
    biorelease(t) = biom(t-1)/biotime*(1-humfrac)*dt;                   %CO2 released to atm. by rotting plants (rest of plant CO2 goes to humus)

    dCat1 = Emis(t)*dt+humrelease(t)+biorelease(t)-NPP(t);  %disregards oceanuptake(t)
    Cat1 = Cat(t-1)+dCat1;                                  %preliminary guess Cat (disregarding ocean uptake) 

    biom(t) = biom(t-1)+NPP(t)-biom(t-1)/biotime*dt; 
    humm(t) = humm(t-1)+biom(t-1)/biotime*humfrac*dt-humm(t-1)/humtime*dt; 

    % OCEAN MIXING
    
    fluxC = zeros(1,nl-1);
    
    for k=1:nl-1; 
        fluxC(k) = eddydif*(Con(t-1,k)/laydep(k)-Con(t-1,k+1)/laydep(k+1))/((laydep(k+1)+laydep(k))/2)*dt;
        %the amount of carbon that goes from layer k into layer k+1
    end

    Con1 = Con(t-1,1)-fluxC(1); %preliminary carbon content of mixed layer
    
    for k=2:nl-1
        Con(t,k) = Con(t-1,k)+fluxC(k-1)-fluxC(k);
    end
    
    Con(t,nl) = Con(t-1,nl)+fluxC(nl-1); 

    % atm-mixedlayer equilibration
    % i do somethin similar to a newton iteration, but i don't feel like
    % differentiating the function fC below. 

    Ctot1 = Con1+Cat1; % total carbon in mixedlayer and atmosphere. 

    % the following function must become zero, if a suitable Ca is chosen. 
    % it is: 
    % total carbon - atm. carbon - "mixlayer carbon which would be in equil. with atm" 
    
    fC = @(Ctot,Ca) Ctot-Ca-Cm0*(1-CmT*T(t-1))*(Ca/Cat0)^(1/(rev0+revC*log(Ca/Cat0)));

    % initial guesses (x0=Cat1; xx0 = somewhere nearby; -2< slope of fC <-1)
    niter = 5; 
    x = zeros(niter,1);
    y = x; 
    xx = x; 
    yy = x; 
    a = x;

    x(1) = Cat1; 
    y(1) = fC(Ctot1,x(1)); 
    xx(1) = x(1)+1.5*y(1); 
    yy(1) = fC(Ctot1,xx(1));            % a guess on the "other side" of fC=0
    a(1) = (yy(1)-y(1))/(xx(1)-x(1));   % the estimated gradient
    suppress = 0;
    
    for k=1:niter-1 %iterations to find correct Ca
        
        x(k+1) = x(k)-y(k)/a(k); 
        y(k+1) = fC(Ctot1,x(k+1));
        
        if y(k+1)<1e-10; x(niter)=x(k+1); suppress=1; break; end
        
        xx(k+1) = x(k)-2*y(k)/a(k); 
        yy(k+1) = fC(Ctot1,xx(k+1));    
        a(k+1) = (yy(k+1)-y(k+1))/(xx(k+1)-x(k+1));
    end %iterations

    if (x(niter)-x(niter-1))/x(niter)>0.0001 && suppress==0
        ['iteration worked badly in time step ' num2str(t)]
    end

    Cat(t) = x(niter);
    Con(t,1) = Ctot1-Cat(t); 
    
%------------------------- TEMPERATURE MODELLING -------------------------%  

    % the surface temperature Hon(t,1), i.e. in mixed layer, changes due to:
    % irradiation FCO2, outgoing radiation Fout [both computed relative to
    % pre-industrial; Fout depends on current temperature), and heat flux to
    % deep ocean. 

    FCO2(t) = 1.3*forCO2*log(Cat(t)/Cat0); % radiative forcing (w.r.t pre-ind) by CO2. in W/m^2
    
    % FCO2(t) = forCO2*log(2); % forced atmospheric carbon doubling (exercise 1.2)
    
    % outrad = FCO2(t)/lambda; % forced outrad setting for climate sensitivity

    Fout(t) = outrad*(T(t-1)); % outgoing longwave rad. (w.r.t. pre-ind) due to warming

    Tequil(t) = FCO2(t)/outrad;

    %heat fluxes. these are per m^2 of ocean surface. 
    fluxH = zeros(1,nl-1);
    
    for k=1:nl-1
        fluxH(k) = eddydif*(Hon(t-1,k)/laydep(k)-Hon(t-1,k+1)/laydep(k+1))/((laydep(k+1)+laydep(k))/2)*dt;
        %the amount of carbon that goes from layer k into layer k+1
    end

    Hon(t,1) = Hon(t-1,1)-fluxH(1)+(FCO2(t)-Fout(t))*secyr/seasurf*dt; %heat content upper layer
    % the factor secyr serves to compute the energy/m^2 collected in a whole
    % year. The factor 1/seasurf means that radiation is collected over the
    % whole earth surface while it ends up only in the sea (ca 70% of total
    % earth surface). 

    for k=2:nl-1
        Hon(t,k) = Hon(t-1,k)+fluxH(k-1)-fluxH(k);
    end
    
    Hon(t,nl) = Hon(t-1,nl)+fluxH(nl-1);

    Ton(t,:) = Hon(t,:)./laydep/heatcap; %divide by heat capacity to get temperature. 
    T(t,:) = Ton(t,1); %for now, surface warming = mixed-layer warming

    %if Ton(t,1) 
    %histo(t) =  
    
end %time steps

    %a few auxiliary variables are not computed in the first time step. for
    % %plotting purposes, set these to second timestep.
    NPP(1) = NPP(2); 
    humrelease(1) = humrelease(2); 
    biorelease(1) = biorelease(2); 
    FCO2(1) = FCO2(2); 
    Fout(1)= Fout(2); 

    fonty = 19;
    colvec1 = ['r','m','b','c','g','k'];

    figure(1); set(gca,'FontSize',fonty); hold on
    title('Temperature')
    plot(timvec,Ton(:,1),'r','Linewidth',2);
    xlabel('time [year]')
    ylabel('temp. change ocean surface [K]')
    
    axis([2000 2100 0 6])


%Tfinal(l) = Tequil(10000);

end

two = 2*ones(1,2001);
plot(timvec,two,'k--','Linewidth',1.5);

%histogram(Ton(:,1))
% Tav = sum(Tfinal)/500

% 
% 
% %a few auxiliary variables are not computed in the first time step. for
% %plotting purposes, set these to second timestep.
% NPP(1) = NPP(2); 
% humrelease(1) = humrelease(2); 
% biorelease(1) = biorelease(2); 
% FCO2(1) = FCO2(2); 
% Fout(1)= Fout(2); 
% 
% fonty = 19;
% colvec1 = ['r','m','b','c','g','k'];
% 
% % this figure serves to check the percentage of cumulative emissions
% % remaining in the atmosphere
% 
% figure; set(gca,'FontSize',fonty); hold on
% title('Atmospheric carbon and cumulative emission')
% plot(timvec,Cat,'Linewidth',2.5)
% plot(timvec,CumEmis+Cat(1),'r','Linewidth',2.5)
% xlabel('time [year]')
% ylabel('amount of Carbon [Gt(C)]')
% legend('atm. Carbon','Pre-ind + cumul. emission')
% 
% %axis([1751 2100 500 3050])
% 
% %this figure serves to understand the size of biosphere processes. In
% %absence of human emissions (i.e. pre-industrial conditions) NPP should be
% %equal to biorelease + humrelease. these quantities are much bigger than
% %emissions (but nearly cancel).
% figure; set(gca,'FontSize',fonty); hold on
% title('Biological processes')
% plot(timvec,NPP,'b','Linewidth',2.5)
% plot(timvec,biorelease,'g','Linewidth',2.5)
% plot(timvec,humrelease,'r','Linewidth',2.5)
% plot(timvec,Emis,'m','Linewidth',2.5)
% plot(timvec,biorelease+humrelease,'c--','Linewidth',2.5)
% xlabel('time [year]')
% ylabel('carbon flux [Gt(C)/year]')
% legend('NPP','biosph. release','humus release','emission','sum bio+hum. release')
% 
% %axis([1751 2100 0 125])
% 
% %the amount of carbon ending up in the ocean. black line: total ocean,
% %other lines: only the upper X layers (i.e. the lowest line is the carbon
% %in the mixed layer)
% Contot=sum(Con,2); Con0tot=sum(Con0); 
%   Conuptaketot=Contot-Con0tot; Conuptakecum=Conuptaketot;
% 
% figure; set(gca,'FontSize',fonty); hold on
% title('Oceanic carbon uptake')
% Contot=sum(Con,2); Con0tot=sum(Con0); 
%   Conuptaketot=Contot-Con0tot; Conuptakecum=Conuptaketot;
% plot(timvec,Conuptaketot,'k','Linewidth',4.5)
% for l=1:nl-1
%     Conuptakecum=Conuptakecum-Con(:,nl-l+1)+Con0(nl-l+1);
%     plot(timvec,Conuptakecum,colvec1(mod(l,length(colvec1))),'Linewidth',2.5);
% end
% xlabel('time [year]')
% ylabel('carbon uptake by the ocean]')
% if nl==5
% legend('total ocean',      ['top ' num2str(sum(laydep(1:4))) 'm'],['top ' num2str(sum(laydep(1:3))) 'm'],...
%        ['top ' num2str(sum(laydep(1:2))) 'm'],['top ' num2str(sum(laydep(1:1))) 'm'])
% elseif nl==2
% legend('total ocean',      ['top ' num2str(sum(laydep(1:1))) 'm']) 
% else legend('total ocean','down to Xth layer')
% end
% 
% %axis([1751 2100 0 525])
%     
% %This plot is the amount of carbon added to the reservoirs. the total
% %should exactly equal the cumulative emissions at each time step. 
% figure; set(gca,'FontSize',fonty); hold on
% title('Overall carbon uptake')
% Ctot=Cat-Cat0;
% plot(timvec,Ctot,'r','Linewidth',2.5);
% Ctot=Ctot+Conuptaketot;
% plot(timvec,Ctot,'b','Linewidth',2.5);
% Ctot=Ctot+biom-biom0;
% plot(timvec,Ctot,'g','Linewidth',2.5);
% Ctot=Ctot+humm-humm0;
% plot(timvec,Ctot,'k','Linewidth',2.5);
% plot(timvec(2:end),CumEmis(1:end-1),'m--','Linewidth',3.5);
% xlabel('time [year]')
% ylabel('change carbon reservoir since pre-ind. [Gt(c)]')
% legend('Atm. carbon','...+ Ocean','...+ Biosphere','...+ humus','Cum. emission')
% 
% %axis([1751 2100 0 2450])
% 
% %this is the temperature (rel. to pre-industrial) of the ocean layers. 
% %the deep ocean continues to warm for a long time! the mixed layer temperature 
% %is taken as "surface" temperature. 
% figure; set(gca,'FontSize',fonty); hold on
% title('Temperature')
% plot(timvec,Tequil,'k--','Linewidth',3.5);
% for l=1:nl
% plot(timvec,Ton(:,l),colvec1(mod(l,length(colvec1))),'Linewidth',2);
% end
% xlabel('time [year]')
% ylabel('temp. change since pre-ind. [K]]')
% if nl==5
% legend('equil. temp.',      'Temp ocean layer 1','Temp ocean layer 2',...
%        'Temp ocean layer 3','Temp ocean layer 4','Temp ocean layer 5')
% elseif nl==2
% legend('equil. temp.','Temp ocean layer 1','Temp ocean layer 2')  
% else legend('equil. temp','temp. ocean layers (top warms fastest)')
% end
% 
% %axis([1751 2100 0 7])
% 
