%simulate full cohort dynamics after creating virtual cohort
close all
clearvars

format long

p = load_parameters();

p.eta_F_M = 5.4*1e-1; 
p.eta_L_M = 0.4498*1e-2;%0.01872*1e2;
p.eps_F_I = 2*1e-4;%0.9*1e-4*2;
p.p_F_I = 2.8235*1e4*1e-4;
p.eta_F_I = 0.0011164*1e1;
p.p_F_M = 3.5600*1e4*1e-4; 
p.p_L_M = 72560/1e2*0.05;
p.p_L_MPhi = 1872;
p.eps_L_MPhi = 1102.9/1e5;
p.eps_G_MPhi = 2664.5/1e5;
p.eta_L_I = 0.7;
p.p_L_I = 1188.7/1e2;
p.eta_L_MPhi = 1e-5;
p.eps_V_MPhi = 905.22/1e3;
p.IC_50_N = 4.7054*1e-2;
p.eps_L_T = 0.3*1e-3;
p.p_T_I = 0.01*0.8;
p.del_I_T = 238;
p.eta_F_MPhi = 1e-5;
p.p_F_MPhi = 1.3; 
p.eps_F_T = 1e-3*1.5;
p.d_V_spec = 5.5;
p.a_I_MPhi = 1100; 
p.eps_I_M = 0.11;
p.del_V_N = 768*3;
p.del_V_MPhi = 768*100;
p.p_G_M = 1.234*1e3*1e-1;
p.d_I = 0.144*0.1;
p.p_MPhi_I_L =  0.42*4;
p.p_MPhi_I_G =  0.42*4;

%---------------------------------------------------
p.V0 = 4.5;
p.phat = 394;
p.beta = 0.3;
p.d_I = 0.1;
%p.d_V = 8.4;
p.d_V_spec = 0;
p.del_V_MPhi = 76800/200;
p.del_V_N = 2304/2.5;
%---------------------------------------------------

p.eps_L_T = 1.5*1e-5;
p.p_T_I = 0.008*2;
p.del_I_T = 238*0.5;

p = Homeostasis_calculations(p);

%---------------------------------------------------
%Load the generated patient cohort


%load('Immunosuppressed.Cohort.mat');
%load('Reference.Cohort.mat');
load('Cancer.Cohort.mat'); %changing initial neutrophils - p.N0 = 0.00526 to 0.00426 (line 69) in the 'load_parameters.m' code

%Loop over the number of virtual patients you want to simulate

for nn = 1:size(patient,1)
    nn
   

    %Reset each of the generated patient values
    %p.beta = patient(nn,1);
    p.p_MPhi_I_L = patient(nn,2);
    p.p_L_MPhi = patient(nn,3);
    p.p_F_I = patient(nn,4);
    %p.eta_F_I = patient(nn,5);
    %p.eps_L_T = patient(nn,6);
    p.p_M_I = patient(nn,7);
    p.eta_F_MPhi = patient(nn,8);
    %p.tau_T = patient(nn,9);
    p.eps_F_I = patient(nn,10);
    p.p_F_M = patient(nn,11);
    %p = Homeostasis_calculations(p);
    %-----------------------------------------------------------------------
    estimated_params = [p.p_L_M  p.L_B_star  p.MPhi_I_star  p.p_L_MPhi  p.eta_G_MPhi  p.p_G_MPhi_I  p.T_prod_star p.T_M_prod_star p.p_G_M  p.M_prod_star p.N_prod_star p.eta_F_MPhi];

    if isempty(find(estimated_params<0))==0
    disp('Negative parameter')
    elseif isempty(find(estimated_params>1e9))==0
    disp('Extremely large parameter')
    end
    %----------------------------------------------------------------------

    %% MODEL WITH VIRUS
    tspan = [0 20];%set the duration of the simulation

    [time,sol,solstruc] = COVID_IMMUNE_MODEL(p,tspan);%call the solver to get the soluation
    
    time_deval = linspace(tspan(1),tspan(2),1000);%set up the time mesh
    sol_deval = deval(solstruc,time_deval);

    Tcells(nn,:)=sol_deval(10,:)*1e9;
 

    monocytes(nn,:)=sol_deval(8,:)*1e9;


    IFN(nn,:)=sol_deval(17,:);
   
     
    IL6(nn,:)=sol_deval(11,:);
  

    GCSF(nn,:)=sol_deval(15,:)*1000;
   

    GMCSF(nn,:)=sol_deval(13,:);
   
    
    Neutrophils(nn,:)=sol_deval(9,:)*1e9;
    

    if isreal(sol_deval)==1
        sol_virus(nn,:) = sol_deval(1,:);
        minimum_tissue(nn) = min(sol_deval(2,:)+sol_deval(4,:));
        max_IL6(nn) = max(sol_deval(11,:));
        max_bound_IL6(nn) = max(sol_deval(12,:));
        max_dead_tissue(nn) = max(sol_deval(5,:));
        max_inflam_macs(nn) = max(sol_deval(7,:));
        max_GCSF(nn) = max(sol_deval(15,:));
        max_IFN(nn)=max(sol_deval(17,:));
        
       
        tgrid_vec = [];
        for i = 1:length(sol_deval(17,:))-1
           tgrid_vec(i) = time_deval(i+1)-time_deval(i); 
        end
        IFN_exposure(nn) = sum(sol_deval(17,1:end-1).*tgrid_vec);
        peak_loc = find(sol_deval(17,:)==max(sol_deval(17,:)));
        IFN_peak(nn) = time_deval(peak_loc); %time of IFN peak
        days14 = find(time>3,1);
        max_T_cells(nn) = max(sol_deval(10,1:days14));
        max_T_cells(nn) = max(sol_deval(10,:));
        max_infected_cells(nn) = max(sol_deval(3,:));
        max_inflam_macs(nn) = max(sol_deval(7,:));
        max_neutrophils(nn) = max(sol_deval(9,:));
        max_monocytes(nn) = max(sol_deval(8,:));
        max_GCSF(nn) = max(sol_deval(15,:));
    else
        minimum_tissue(nn) = NaN;%+sol(3,:)
        max_IL6(nn) = NaN;
        max_bound_IL6(nn) = NaN;
        max_dead_tissue(nn) = NaN;
        IFN_exposure(nn) = NaN;
        IFN_peak(nn) = NaN; %time of IFN peak
        max_T_cells(nn) = NaN;    
        max_infected_cells(nn) = NaN;
        max_inflam_macs(nn) = NaN;
        max_neutrophils(nn) = NaN;
        max_monocytes(nn) = NaN;
        max_GCSF(nn) = NaN;
    end

end

tspan=[0 20];
time_grid = linspace(tspan(1),20,1000);
p.Smax = 0.16;
%[max_IL6, max_neutrophils, minimum_tissue, max_inflam_macs, IFN_exposure, IFN_peak, max_T_cells]=reevaluate_patients(time_grid,sol_virus,sol_tissue,sol_infected,sol_dead,sol_macs_res,sol_macs_inflam,sol_neutrophils, sol_monocytes,sol_Tcells, sol_IL6, sol_IFN, sol_GCSF,time_grid);
%Psi = (max_IL6)/(mean(max_IL6))+(max_neutrophils)/mean(max_neutrophils)+(p.Smax-minimum_tissue)/mean(p.Smax-minimum_tissue);
Psi = (max_IL6)/30.23+(max_neutrophils)/0.006325488688718+(p.Smax-minimum_tissue)/0.121176119035704; %adjusted inflamation marker
[sorted,where] = sort(real(Psi));

%Fitting patients dynamics to the data (examples)

figure(17)
options.handle     = figure(17);
%Cancer
options.color_area = hex2rgb('#AE445A');    % Cancer theme
options.color_line ='#662549';
%Immunosuppressed
%options.color_area = hex2rgb('#427D9D');    % Immunosuppressed theme
%options.color_line =hex2rgb('#164863');
%Reference
%options.color_area = hex2rgb('#78D6C6');    % Reference theme
%options.color_line =hex2rgb('#419197');
options.alpha      = 0.5;
options.line_width = 2;
options.x_axis = time_deval;
options.error      = 'std';
plot_areaerrorbar(Tcells,options)
hold on
%Cancer
e=errorbar(8,1.9E6,0.15E6); %(day of measurements, mean value, standard deviation)
f=errorbar(14,1E6,0.35E6);
%Only Covid
%e=errorbar(7,1.7E6,0.6E6);
%f=errorbar(14,2.5E6,0.95E6);
e.Color='#7E2F8E';
e.LineWidth=2;
e.Marker='^';
e.MarkerFaceColor='#7E2F8E';
f.Color='#7E2F8E';
f.LineWidth=2;
f.Marker='^';
f.MarkerFaceColor='#7E2F8E';
set(gca,'FontSize',18)
ylim([0 3.5e6])
xlabel('Time (days)')
ylabel('T cells (cells/ml)')
saveas(gcf,'Fig_IM1.fig');

figure(19)
options.handle     = figure(19);
%Cancer
options.color_area = hex2rgb('#AE445A');    
options.color_line = '#662549';
%Immunosuppressed
%options.color_area = hex2rgb('#427D9D');   
%options.color_line =hex2rgb('#164863');
%Reference
%options.color_area = hex2rgb('#78D6C6');  
%options.color_line =hex2rgb('#419197');
options.alpha      = 0.5;
options.line_width = 2;
options.x_axis = time_deval;
options.error      = 'std';
plot_areaerrorbar(IFN,options)
hold on
%Immunosuppressed
%e=errorbar(5,0.24,0.2);
%f=errorbar(7,0.1,0.06);
%g=errorbar(11,0.08,0.06);
%h=errorbar(16,0.16,0.13);
%Cancer
e=errorbar(4,0.32,0.17);
f=errorbar(7,0.13,0.06);
g=errorbar(11,0.09,0.07);
h=errorbar(13,0.11,0.10);
%Only Covid
%e=errorbar(1,0.28,0.2);
%f=errorbar(4,0.23,0.22);
%g=errorbar(9,0.12,0.11);
%h=errorbar(11,0.14,0.13);
%i=errorbar(13,0.08,0.07);
e.Color='#7E2F8E';
e.LineWidth=2;
e.Marker='^';
e.MarkerFaceColor='#7E2F8E';
f.Color='#7E2F8E';
f.LineWidth=2;
f.Marker='^';
f.MarkerFaceColor='#7E2F8E';
g.Color='#7E2F8E';
g.LineWidth=2;
g.Marker='^';
g.MarkerFaceColor='#7E2F8E';
h.Color='#7E2F8E';
h.LineWidth=2;
h.Marker='^';
h.MarkerFaceColor='#7E2F8E';
%i.Color='#7E2F8E';
%i.LineWidth=2;
%i.Marker='^';
%i.MarkerFaceColor='#7E2F8E';
set(gca,'FontSize',18)
xlabel('Time (days)')
ylabel('IFN (pg/ml)')
ylim([0 inf])
saveas(gcf,'Fig_IM4.fig');

%Correlation plot - example

%max_dead_tissue=linspace(min(max_dead_tissue),max(max_dead_tissue),10);


cStart = [100,201,207]./255;%define the starting colour as an RGB triplet
cEnd = [255,93,93]./255;%define the ending colour as an RGB triplet
c = interp1([1;10],[cStart;cEnd],(1:10)');
%Immunosuppressed
%bins=linspace(0.4e-3,1.3e-3,10)';
%Covid
bins=linspace(0.001,0.0026,10)';
%Cancer
%bins=linspace(0.5e-3,1.8e-3,10)';
figure
hold on 
for k = 1:length(patient)
    [ind,~]=find(bins<max_dead_tissue(where(k)),1,'last');
    if isempty(ind)
        ind=1;
    end
 plot(k,max_IFN(where(k)),'.','Color',c(ind,:),'MarkerSize',30)%plot each patient's IFN according to their order in the 'sort' vector 
end
  
xlabel('Patients ordered by inflammation \Psi','FontSize',18)
ylabel('Max IFN (pg/ml)','FontSize',18)
ylim([0 0.43])
set(gca,'FontSize',18)
[R,P] = corrcoef(max_dead_tissue(where),max_IFN(where));
text(2,0.41,['R = ' num2str(R(1,2)),', p =' num2str(P(1,2))],'FontSize',18)
colormap(c)%use the colourmap defined above as the colour
colorBounds = [1 2.6];%define the bounds based on the min and max IL-6 concentrations
    caxis(colorBounds)
    a=colorbar();  
a.Label.String = "Max damaged tissue (cells x 10^6/ml)";
a.Label.FontSize = 15;
a.Label.Rotation = 270;
a.Label.VerticalAlignment = "bottom";


