%This code takes raw Pb-210 data for a core and mean Ra-226 values and
%calculates interval size, aerial dry mass, inventory, and age models for
%the three common models (CFCS, CF, & CA as described in Sanchez-Cabeza &
%Ruiz-Fernandez 2012). It requires a CSV file with some Pb-210 values, the
%year of the core, at least one Ra-226 value, and dry bulk density for
%every centimeter with Pb-210 data. If there are missing Pb-210 values,
%this code exponentially interpolates a middle value between existing
%intervals using the exponential interpolation described in Appleby 2001.
%If there are missing DBD values, the code linearly interpolates between
%real values for every missing centimeter. This code works with complete
%data sets as well.

% To run this code, simply define the folder where you have saved this code
% and the import_pb_csv function, the folder where your CSV file is saved,
% where you want the figures produced to be saved, and which format you
% want figures to be saved as. Hit run and when prompted make sure the
% folder shown is correct, all remaining data and four figure will then be
% created and saved to that folder

%Written by Tyler Souza, September 2022, tsouza@unc.edu


clear; close;

codefolder = 'C:\Users\souza\Dropbox\CoosBay\Alpha\Modeling'; %Where the code is save
datafolder = 'C:\Users\souza\Dropbox\CoosBay\Alpha\Modeling\BC_Model'; %Where your CSV data lives
outputfolder = 'C:\Users\souza\Dropbox\CoosBay\Alpha\Modeling\BC_Model'; %Where you want your data output to go
site = cellstr("BC");
addpath(codefolder);
addpath(datafolder);
figures_fileType = 'png';       % [ENTER 'png', 'jpg', 'pdf', etc.]
%% Load the Pb-210 base data
import_pb_csv;
load('model.mat')
%% Assign non-missing values to variables
temp = find(~isnan(depth)); %find the location of values that are not NaN
depth = depth(1:temp(end)); %assign non missing depths
ratotal = ratotal(1);

year = year(1);
ura = ura(1:temp(end));
upb = upb(1:temp(end));
dbd = dbd(1:temp(end));
udbd = udbd(1:temp(end));
pbtotal = pbtotal(1:temp(end));
deltaz(1:length(depth)) = 1;
clear temp  
lambda = .03118071; %Pb-210 decay constant
ulambda = .00017;
dbd(dbd<0) = NaN;
dbd(dbd==2.5) = NaN;

%% Find the index of missing and non-missing intervals

ndata = find(~isnan(pbtotal));
nandata = find(isnan(pbtotal));

%% Find the middle of the interval
for i = 2:length(ndata)
    if ndata(i)-ndata(i-1)>1
    miss_depth(i-1) = (ndata(i)+ndata(i-1))/2; %find the middle of all the missing intervals
    
    else
    miss_depth(i-1) = NaN;
    end        
end
if ~isempty(nandata)
    if nandata(1) == 1 % if the top of the core is missing the code will using a linear estimate for the surface
        miss_depth(2:length(miss_depth)+1) = miss_depth;
        miss_depth(1) = 1;
    end
    if length(ndata)>length(miss_depth)
        miss_depth(length(miss_depth):length(ndata)) = NaN;
    end
else
    miss_depth(length(ndata)) = NaN;
end
%%
mid_depth = [rot90(ndata);miss_depth;]; %assign the middle depth in core of all intervals both interpolated and real
mid_depth = mid_depth(:);
mid_depth = sort(mid_depth);
mid_depth = rmmissing(mid_depth);
[a,b,i_interp] = intersect(miss_depth, mid_depth); %find the index of points in mid_depth where values will be interpolated
[a,b,i_real] = intersect(ndata, mid_depth); %find the index of points in mid_depth where values are real
clear a b
if miss_depth(1) ==1
mid_depth(1) = (mid_depth(2)-mid_depth(1))/2;
end

%% Find the interval thickness
for i = 2:length(mid_depth)-1
    Zinext(i) = (mid_depth(i+1)-mid_depth(i))/2; %half the distance to the next interval
    Ziprev(i) = (mid_depth(i)-mid_depth(i-1))/2; %half the distance to the previous interval
    Zi(i) = Zinext(i)+Ziprev(i); %half the distance to the previous interval and half the distance to the next interval
end
if ~isempty(nandata)
if nandata(1) ==1
Zi(1) = 2*mid_depth(1);
Ziprev(1) = .5*Zi(1);
Zinext(1) = .5*Zi(1);
Zi(2) = (mid_depth(2)-Zi(1))/2+(mid_depth(3)-mid_depth(2))/2;
Ziprev(2) = .5;
end
end
if ndata(1) == 1
    Ziprev(1) = .5;
    Zinext(1) = (mid_depth(2)-mid_depth(1))/2;
    Zi(1) = Ziprev(1) + Zinext(1);
end

Zinext(length(mid_depth)) = 0.5;
Ziprev(length(mid_depth)) =(mid_depth(end)-mid_depth(end-1))/2;
Zi(length(mid_depth)) = length(dbd)- mid_depth(end)+(mid_depth(end)-mid_depth(end-1))/2+.5;

%% find the dbd for the cumulative layers
dbd_interp = dbd;
if length(rmmissing(dbd)) < length(dbd) %interpolate dbd if necessary
    
    for i = 1:length(nandata)
        temp1 = ndata(find(ndata<nandata(i)));
        if temp1 > 0
        xstart(i) = temp1(end);
        temp2 = ndata(find(ndata>nandata(i)));
        xend(i) = temp2(1);    
        end

    end
    for i = 1:length(nandata)
      y1 = dbd(xstart(i));
      x1 = depth(xstart(i));
      y2 = dbd(xend(i));
      x2 = depth(xend(i));
      x = depth(nandata(i));
      y = y1 +(((x-x1)*(y2-y1))/(x2-x1));
      dbd_interp(nandata(i)) = y;
      clear y1 x1 y2 x2 x y
    end 
   
end

j = 0;
for i = 1:length(Zi)
        intstart(i) = j;   
        j = j+Zi(i);
        intend(i) = j;
        dbdsum(i) = sum(dbd_interp(ceil(intstart(i)+1):floor(intend(i)))) + ((ceil(intstart(i))-(intstart(i)))*dbd_interp(ceil(intstart(i))+1))+((intend(i)-floor(intend(i)))*dbd_interp(ceil(intend(i))));
end

adm = dbdsum; %.*Zi; %aerial dry mass

%% calculate mass depth
Mil(1) = 0;
for i = 2:length(mid_depth)+1
    Mil(i) = Mil(i-1)+adm(i-1);
end
for i = 2:length(mid_depth)+1
    Mis(i-1) = (Mil(i-1)+Mil(i))/2;
end

for i = 1:length(i_interp)
Mis_interp(i) = Mis(i_interp(i));
end
for i = 1:length(i_real)
Mis_real(i) = Mis(i_real(i));
end

%% Assign and interpolate Pb-210 values
for i = 1:length(i_real)
    Tot_pb(i_real(i)) = pbtotal(mid_depth(i_real(i))); %Assign total Pb-210 values to real intervals in mid_depth
    Real_pb(i) = pbtotal(mid_depth(i_real(i))); % Assign total Pb-210 values to only real intervals
    Pb_error(i_real(i)) = upb(mid_depth(i_real(i))); % Assign error to real intervals
end

%interpolate exponentially to find missing values
Interp_pb = [];
if ~isempty(nandata)
if i_interp(1) == 1
    n = 2;
else
    n=1;
end

for i = n:length(i_interp)
    Tot_pb(i_interp(i)) = (Tot_pb(i_interp(i)-1)-Tot_pb(i_interp(i)+1))/log(Tot_pb(i_interp(i)-1)/Tot_pb(i_interp(i)+1)); %Exponential interpolation from Appleby 2001 (B = (A-C)/ln(A/C))
    Interp_pb(i) = (Tot_pb(i_interp(i)-1)-Tot_pb(i_interp(i)+1))/log(Tot_pb(i_interp(i)-1)/Tot_pb(i_interp(i)+1));
    Pb_error(i_interp(i)) = (Pb_error(i_interp(i)-1)-Pb_error(i_interp(i)+1))/log(Pb_error(i_interp(i)-1)/Pb_error(i_interp(i)+1));
end
end
%if surface is missing, assign the top Pb value to be a linear estimate of the
%top 4 real points
if ~isempty(nandata)
if i_interp(1) ==1
    y = rot90(log(Real_pb(1:4)));
    x = rot90(Mis_real(1:4));
    X = [ones(length(x),1) x];
   linest = X\y;
   Tot_pb(1) = exp(linest(1));
   Interp_pb(1) = Tot_pb(1);
   C0 = Tot_pb(1);
end
end

% if the top of the core has a real value calculate C0
if i_real(1) ==1
    y = rot90(log(Real_pb(1:4)));
    x = rot90(Mis_real(1:4));
    X = [ones(length(x),1) x];
   linest = X\y;
   C0 = exp(linest(1));
end

%% User input on background
scatter(Tot_pb,mid_depth)
hold on
xline(ratotal,'linewidth',2)
set(gca,'Ydir','reverse')
xlabel('Activity (dpm/g)')
ylabel('Depth (cm)')
title('Do you want to use the Ra value for supported Pb-210?')
findRa = input("Do you want to use the Ra value for supported Pb-210? (y/n)","s");
close;
%% Calculate excess Pb-210
if findRa == 'y'
[r,c, V] = findnearest(ratotal,pbtotal);
background = V;
clear r c V
findAj = 'y';
end
if findRa == 'n'
    findAj = input("Do you want to calculate Ra from the bottom of the core? (y/n)","s");
    if findAj == 'y'
    scatter(Tot_pb,mid_depth)
    set(gca,'Ydir','reverse')
    xlabel('Activity (dpm/g)')
    ylabel('Depth (cm)')
    title('Select the points you want to use as background')
    [bottom, bottom_y] = ginput;
    background = mean(bottom);
    close
    end
    if findAj == 'n'
        background = input("What is the background value for this core");
    end
end

XS_pb = Tot_pb - background;
XS_pb_real = Real_pb-background;
if ~isempty(nandata)
XS_pb_interp = Interp_pb - background;
else
    XS_pb_interp = [];
end
if findRa == 'n' && findAj == 'y'
XS_pb(mid_depth >=min(bottom_y)) = 0;
end
XS_pb=max(XS_pb,0);

%% CFCS
y = rot90(XS_pb,3);
ln_Ci = log(y);
ln_Ci(XS_pb == 0) = [];
x1 = mid_depth(1:length(ln_Ci));
X = [ones(length(x1),1) x1];
slope = X\ln_Ci;
CFCS_SAR = -lambda/slope(2);
x2 = rot90(Mis,3);
x2 = x2(1:length(x1));
X = [ones(length(x2),1) x2];
slope = X\ln_Ci;
CFCS_MAR = -lambda/slope(2);

for i = 1:length(mid_depth)
    CFCS_age(i) = (mid_depth(i)-mid_depth(1))/CFCS_SAR;
    CFCS_year(i) = year-CFCS_age(i);
end

if findRa == 'n' && findAj == 'n'
    Aj = (XS_pb(end)*CFCS_MAR)/lambda;
end
% hold on
% scatter(Mis,ln_Ci)
% y_cfcs = slope(2)*X(1:end,2)+slope(1);
% plot(Mis(1:end), y_cfcs);
% Rsq2 = 1 - sum((ln_Ci(1:end) - y_cfcs).^2)/sum((ln_Ci - mean(ln_Ci)).^2);
%% calculate inventory
inv = XS_pb.*adm; %inventory is excess * aerial dry mass

for i = 1:length(XS_pb)
uam(i) = dbdsum(i)*sqrt(sumsqr([udbd(i)/dbd(i); .02]));
uinv(i) = inv(i)*sqrt(sumsqr([upb(i)/XS_pb(i); uam(i)/dbd(i)]));
end
invbel = zeros(length(mid_depth)+1,1);
uinvbel = zeros(length(mid_depth)+1,1);

for i = length(inv):-1:1
    if findAj == 'n' && i == length(inv)
        invbel(length(inv)+1) = Aj;
    end
   invbel(i) = invbel(i+1)+inv(i);
   uinvbel(i) = sqrt(sumsqr([uinvbel(i+1) uinv(i)]));
end


    


%% CF

for i = 1:length(XS_pb)
    CF_age(i) = log(invbel(1)/invbel(i))/lambda;
    CF_year(i) = year-CF_age(i);
    CF_MAR(i) = lambda * invbel(i)/XS_pb(i);
    CF_uMAR(i) = CF_MAR(i)*sqrt(sumsqr([ulambda/lambda; uinvbel(i)/invbel(i); upb(i)/XS_pb(i)]));
    CF_SAR(i) = CF_MAR(i)/dbdsum(i);
    CF_uSAR(i) = CF_SAR(i)*sqrt(sumsqr([CF_uMAR(i)/CF_MAR(i); udbd(i)/dbdsum(i)]));
end

CF_MAR(CF_MAR ==inf) = nan;
CF_SAR(CF_SAR ==inf) = nan;
fig_start = (CF_year(CF_year>0));
fig_start = fig_start(end);
%% CA
CA_age_layer(1) = 0;
for i = 1:length(XS_pb)
    CA_age(i) = log(C0/XS_pb(i))/lambda;
    CA_year(i) = year-CA_age(i);
end

for i = 2:length(CA_age)
       CA_age_layer(i) = (CA_age(i)+CA_age(i-1))/2;
end

for i = 2:length(CA_age_layer)
    delta_t(i) = CA_age_layer(i)-CA_age_layer(i-1);
end

for i = 1:length(CA_age)
    CA_SAR(i) = (intstart(i)-intend(i))/delta_t(i);
end

for i = 2:length(CA_age)+1
    CA_MAR(i-1) = (Mil(i)-Mil(i-1))/delta_t(i-1);
end

%% Summary Figure
fig_1 = figure(1);
subplot(131)
errorbar(XS_pb, mid_depth, Pb_error,'o', 'horizontal');
xlabel('Excess Pb-210 Activity (dpm)')
ylabel('Depth (cm)')
set(gca, 'fontsize', 12);
set(gca, 'YDir','reverse')
xlim([-0.5 max(Tot_pb)+.5])
subplot(132)

errorbar(CF_MAR, CF_year,CF_uMAR, 'horizontal', 'linewidth', 2);
ylabel('Year', 'fontsize', 12)
xlabel('Mass Accumulation Rate (g/cm^2/yr)', 'fontsize', 12)
grid on
set(gca, 'fontsize', 12);
%xlim([0 4])
ylim([fig_start year])

subplot(133)
errorbar(CF_SAR, CF_year, CF_uSAR, 'horizontal', 'linewidth', 2);
ylabel('Year', 'fontsize', 12)
xlabel('Sed Accumulation Rate (cm/yr)', 'fontsize', 12)
%xlim([0 4])
ylim([fig_start year])

grid on
set(gca, 'fontsize', 12);
set(figure(1),'color','w','position',[100 100 1200 600])
saveas(fig_1, strcat('PBex.', figures_fileType));

%% MAR figure
fig_2 = figure(2);
plot1 = errorbar(CF_year, CF_MAR, CF_uMAR, 'Vertical', 'linewidth', 2);
xlabel('Year', 'fontsize', 16)
ylabel('Mass Accumulation Rate (g/cm^2/yr)', 'fontsize', 16)
xlim([fig_start year])
%ylim([0 3])
grid on
set(gca, 'fontsize', 12);
ax = plot1.Parent;   
set(ax, 'XTick', round(fig_start,-1):10:year)
xtickangle(45)
set(figure(2),'color','w','position',[100 100 600 300]);
saveas(fig_2, strcat('MAR.', figures_fileType));

%% SAR Figure
fig_3 = figure(3);
plot1 = errorbar(CF_year, CF_SAR, CF_uSAR, 'vertical', 'linewidth', 2);
xlabel('Year', 'fontsize', 16)
ylabel('Sed Accumulation Rate (cm/yr)', 'fontsize', 16)
xlim([fig_start year])
%ylim([0 3])
grid on
set(gca, 'fontsize', 12);
ax = plot1.Parent;   
set(ax, 'XTick', round(fig_start,-1):10:year)
xtickangle(45)
set(figure(3),'color','w','position',[100 100 600 300]);
saveas(fig_3, strcat('SAR.', figures_fileType));

%% Comparison Figure
fig_4 = figure(4);
hold on
plot(CA_age,mid_depth,'linewidth',2)
plot(CF_age,mid_depth,'linewidth',2)
plot(CFCS_age,mid_depth,'linewidth',2)
set(gca, 'fontsize', 12);
set(gca, 'YDir','reverse')
legend('CA','CF','CFCS')
xlabel('Age (years)')
ylabel('Depth in Core (cm)')
title('Age vs Depth For The Three Models')
saveas(fig_4, strcat('Comparison.', figures_fileType));
%%
fig5 = figure(5);
scatter(Tot_pb,mid_depth)
hold on
xline(background,'linewidth',2)
set(gca,'Ydir','reverse')
xlabel('Activity (dpm/g)')
ylabel('Depth (cm)')

%% plot the real and interpolated values in log space
% if ~isempty(nandata)
% hold on
% scatter(log(XS_pb_real),log(Mis_real),'filled')
% scatter(log(XS_pb_interp),log(Mis_interp), 'filled')
% set(gca,'ydir','reverse')
% end

%% Plot the interval widths

% fig6 = figure(6);
% CM = cbrewer('qual','Accent',7);
% CM2 = repmat(CM,length(intstart)/length(CM),1);
% for i = 2:length(CF_year)
%     width(i) = CF_year(i)-CF_year(i-1);
% end
% 
% for i = 1:length(intstart)
% hold on
% xs = [-5 5 5 -5];
% ys = [-intstart(i) -intstart(i) -intend(i) -intend(i)];
% patch(xs,ys,CM2(i,:),'edgecolor','black')
% end
% 
% xlim([-20,20])

%% plot the real and interpolated values in real space
fig7 =figure(7);
subplot(1,2,1)
if ~isempty(nandata)
hold on
% scatter(XS_pb_real,Mis_real,'filled')
% scatter(XS_pb_interp,Mis_interp, 'filled')
scatter(XS_pb_interp, miss_depth,'filled')
scatter(XS_pb_real,real_depth,'filled')
set(gca,'ydir','reverse')
end
subplot(1,2,2)
hold on
scatter(dbdsum,mid_depth)
scatter(dbd_interp,depth)
set(gca,'ydir','reverse')
%% Save all the data into one structure and clear the workspace
Core_data.meta.site_name = site;
Core_data.meta.core_name = csvfile;
Core_data.meta.core_length = max(mid_depth);
Core_data.meta.year = year;
Core_data.raw_data.total_pb = pbtotal;
Core_data.raw_data.uPb = upb;
Core_data.raw_data.Ra226 = ratotal;
Core_data.raw_data.uRa = ura(1);
Core_data.raw_data.depth = depth;
Core_data.raw_data.dbd = dbd;
Core_data.raw_data.lambda = lambda;
Core_data.raw_data.ulambda = ulambda;
Core_data.raw_data.udbd = udbd;
Core_data.interpolation.cumulative_dbd = dbdsum;
Core_data.interpolation.background = background;
Core_data.interpolation.C0 = C0;
Core_data.interpolation.interpolated_index = i_interp;
Core_data.interpolation.interpolated_Pb = Interp_pb;
Core_data.interpolation.real_index = i_real;
Core_data.interpolation.interval_start = intstart;
Core_data.interpolation.interval_end = intend;
Core_data.interpolation.inventory = inv;
Core_data.interpolation.inventory_below = invbel;
Core_data.interpolation.ln_Ci = ln_Ci;
Core_data.interpolation.mass_depth_layer = Mil;
Core_data.interpolation.mass_depth_section = Mis;
if ~isempty(nandata)
Core_data.interpolation.mass_depth_interpolated = Mis_interp;
end
Core_data.interpolation.mass_depth_real = Mis_real;
Core_data.interpolation.missing_depths = miss_depth;
Core_data.interpolation.middle_of_interval = mid_depth;
Core_data.interpolation.original_real_index = ndata;
Core_data.interpolation.uPb = Pb_error;
Core_data.interpolation.real_Pb = Real_pb;
Core_data.interpolation.total_Pb = Tot_pb;
Core_data.interpolation.uADM = uam;
Core_data.interpolation.uInventory = uinv;
Core_data.interpolation.uInventory_below = uinvbel;
Core_data.interpolation.excess_Pb = XS_pb;
Core_data.interpolation.real_excess_Pb = XS_pb_real;
Core_data.interpolation.interpolate_excess_Pb = XS_pb_interp;
Core_data.interpolation.interval_width = Zi;
if ~isempty(nandata)&& findAj == n
Core_data.interpolation.missing_inventory = Aj;
end
Core_data.CFCS.age = CFCS_age;
Core_data.CFCS.year = CFCS_year;
Core_data.CFCS.MAR = CFCS_MAR;
Core_data.CFCS.SAR = CFCS_SAR;
Core_data.CA.age = CA_age;
Core_data.CA.age_layer = CA_age_layer;
Core_data.CA.MAR = CA_MAR;
Core_data.CA.SAR = CA_SAR;
Core_data.CA.year = CA_year;
Core_data.CF.age = CF_age;
Core_data.CF.year = CF_year;
Core_data.CF.MAR = CF_MAR;
Core_data.CF.Mean_MAR = nanmean(CF_MAR);
Core_data.CF.SAR = CF_SAR;
Core_data.CF.Mean_SAR = nanmean(CF_SAR);
Core_data.CF.present_MAR = nanmean(CF_MAR(CF_year>1990));
Core_data.CF.past_MAR = nanmean(CF_MAR(CF_year<1940));
Core_data.CF.delta_MAR = Core_data.CF.present_MAR-Core_data.CF.past_MAR;
Core_data.CF.present_SAR = nanmean(CF_SAR(1:5));
Core_data.CF.past_SAR = nanmean(CF_SAR(CF_year<1940));
Core_data.CF.delta_SAR = Core_data.CF.present_SAR-Core_data.CF.past_SAR;
Core_data.CF.uMAR = CF_uMAR;
Core_data.CF.uSAR = CF_uSAR;

 clearvars -except Core_data
 save('Core_data');
% % 

