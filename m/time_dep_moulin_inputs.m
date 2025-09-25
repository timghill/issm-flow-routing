% Example of time-dependent moulin inputs
addpath(genpath([getenv('ISSM_DIR'), '/bin']))
addpath(genpath([getenv('ISSM_DIR'), '/lib']))


% Load previously mapped moulins
moulin_map = load('moulins.mat');
moulins = moulin_map.moulins;
catchments = moulin_map.catchments;

% Compute melt: seasonal sinusoid x lapse rate
% Load ISSM model and set surface and melt fields
% loadmodel('ASE_2300_ks_1e3_GlaDS_Steady_State_ks.mat')

% Compute melt by surface elevation
z = md.geometry.surface;
ela = 1500;
melt = 2/365/86400*(ela-phi)/ela;
melt(melt<0) = 0;

% Get melt on elements
elmelt = mean(melt(md.mesh.elements), 2);
disp('elmelt')
size(elmelt)

% Areas
area = GetAreas(md.mesh.elements, md.mesh.x, md.mesh.y);

% Seasonal part
tt = 0:365;
seasonal = 0.5 - cos(2*pi*(tt/365));
seasonal_melt = elmelt*seasonal;
seasonal_melt(seasonal_melt<0) = 0;
size(seasonal_melt)


% Accumulate into moulins
nmoulin = length(moulins);
nt = length(tt);
moulin_inputs = zeros(nmoulin, nt);
for i=1:nmoulin
    elements = catchments{moulins(i)};
    moulin_inputs(i, :) = sum(area(elements).*seasonal_melt(elements, :), 1);
end

figure;
plot(tt, moulin_inputs)
xlabel('Day of year')
ylabel('Moulin discharge (m^3/s)')
xlim([0, 365])
grid on
print('moulin_discharge', '-dpng', '-r400')
