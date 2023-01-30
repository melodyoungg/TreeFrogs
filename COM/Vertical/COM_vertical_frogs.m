% nicholas flaim

% this script will take an input from a subfolder either in the form of
% 8channel .au2 files, Fxyz files, Axyz data or positional data, crops each stride
% within each trial and export all data related to energy recovery in
% different excel files


close all
% clear
clc

% inputs
write_data = false;
subject_mass = 0.05425;                       % [kg]  (will be overwritten in force trials)
fps_DLC = 125;
cutoffFreq = 15;
barlength = 0.1 ;                       % distance between 2 borders [m] (only for positional data)
orientation_grav = [1,0,0];             % [x,y,z] (put the 1 in the axis that aligns with gravity)
orientation_vel = [1,0,0];              % [x,y,z]
substrate = 'Treadwall';

% orientation cell
orientation_grav_labels = {'Vertical','Vertical - Mediolateral','Upright'};
        
% creates folders if they do not exist
if ~isfolder('Positional Data')
    mkdir('Positional Data')
end
cd('Positional Data')
if ~isfolder('COM Data')
    mkdir('COM Data')
end
cd ..;

if ~isfolder('Force8')
    mkdir('Force8')
end
cd('Force8')
if ~isfolder('COM Data')
    mkdir('COM Data')
end
cd ..;

if ~isfolder('Accel Data')
    mkdir('Accel Data')
end
cd('Accel Data')
if ~isfolder('COM Data')
    mkdir('COM Data')
end
cd ..;

if ~isfolder('Force3')
    mkdir('Force3')
end
cd('Force3')
if ~isfolder('COM Data')
    mkdir('COM Data')
end
cd ..;

% asking user to select file
disp('Select File for analysis')
[filename,pathname,~] = uigetfile('*');
fullfilename = fullfile(pathname,filename);

% determining type file type (force or positional)
% if count(fullfilename,'Force8') >= 1
%     filetype = 'force8';
% elseif count(fullfilename,'Positional Data') >= 1
%     filetype = 'positional';
% elseif count(fullfilename, 'Force3') >= 1
%     filetype = 'force3';
% elseif count(fullfilename, 'Accel Data') >= 1
%     filetype = 'accel3';
% else
%     error(strcat("Unsure what type of file this is: '",filename,"'"));
% end

filetype = 'force8';
% if file is Force8, creates .csv file
if strcmpi(filetype,'force8')
    cd('Force8')
    fullerfilename = strcat(pathname,'ExtractedData/',strrep(filename,'.au2','.csv'));
    if isfile(fullerfilename)
        cd ExtractedData
    else
        extractau2revised()     % note: running this func puts you in ExtractedData folder
    end
    pathname = pwd;
    filename = strrep(filename,'.au2','.csv');
    fullfilename = fullfile(pathname,filename);
    
    % saving force file info bc og variables will be overwritten later
    filename_force = filename;
    fullfilename_force = fullfilename;
    
    % creates force file folder if it does not exist and returns to
    % original folder
    if ~isfolder('Force Files')
        mkdir('Force Files')
    end
    cd ..; cd ..;
end

% if force8 data, converts 8channel into Fxyz
% if force3 data, just reads that data
if strcmpi(filetype,'force8')
    [force_xyz, cal_t_bounds] = Raw8ChannelData_to_ForceXYZ(fullfilename);
elseif strcmpi(filetype,'force3')
    force_xyz = readmatrix(fullfilename);
end

% if force data, exports data (only if force8), then does low pass filter,
% calcs Fnet, accel_xyz and gets positional data file
if contains(filetype,'force') == 1
    % crop useful data
    fig1 = figure(1);
    plot(force_xyz(:,1),force_xyz(:,2),force_xyz(:,1),force_xyz(:,3),...
        force_xyz(:,1),force_xyz(:,4))
    title({'Please select two points between which the values will be kept.',...
        'All data outside the points selected will be ignored'});
    [startstop,~] = ginput(2);
    startstop = sort(round2actual(force_xyz(:,1),startstop));
    close(fig1)
    force_xyz_export = force_xyz(force_xyz(:,1) >= startstop(1) &...
        force_xyz(:,1) <= startstop(2),:);
    
    % exports forceXYZ data to 'Force Files' folder if force8 data
    if strcmpi(filetype,'force8')
        cd('Force8')
        cd ExtractedData;
        cd('Force Files');
        newfilename = strcat('Force_Data_',filename);
        csvwrite(newfilename,force_xyz_export);
        cd ..; cd ..; cd ..;
        fprintf('%s created and placed in ''Force Files'' folder.\n', newfilename);
    end
    
    % getting stride information
    cd Force8;
    stride_file = 'StrideFrames_Frogs_COM.xlsx';
    stride_data = readcell(stride_file);
    cd ..;
    event_name = extractBefore(filename,'-Camera');
    stride_labels = stride_data(1,:);
    stride_data = stride_data(strcmpi(stride_data(:,2),event_name),:);
    if size(stride_data,1) == 0
        error(strcat(event_name,'_is not found in_',stride_file))
    end
    
    % number of strides and avg velocity
    num_strides = size(stride_data,1);
    try
        velocity = cell2mat(stride_data(:,13));
    catch
        error('no velocity for this trial')
    end
    if ismissing(stride_data{1,9})  % if no mid stride data is there
        error(strcat('missing stride data needed for stance & swing phase, check:',...
            stride_file,'_ to confirm'));
    end
    stride_times = cell2mat(stride_data(:,7:9))/120;
    
%     overwriting old mass
    plot(force_xyz(:,1),force_xyz(:,2))
    title({'select two points between which'},{'represent the subjects body weight'});
    [subject_mass_bounds, ~] = ginput(2);
    subject_mass_bounds = sort(subject_mass_bounds);
    subject_mass_bounds = round2actual(force_xyz(:,1),subject_mass_bounds);
    subject_mass = force_xyz(force_xyz(:,1) >= subject_mass_bounds(1) & ...
        force_xyz(:,1) <= subject_mass_bounds(2),2);
    subject_mass = mean(subject_mass)/9.807;
    
    % remove offsets and drift 
    force_old = force_xyz;
    %force_xyz = drift_correction(force_xyz,cal_t_bounds);
    %force_xyz = offset_correction(force_xyz,cal_t_bounds);
    force_xyz = force_xyz(force_xyz(:,1) >= startstop(1) &...
        force_xyz(:,1) <= startstop(2),:);
    
    % low pass fourier transform
    force_xyz(:,2) = real(lowpassFourier(force_xyz(:,1),force_xyz(:,2),cutoffFreq));
    force_xyz(:,3) = real(lowpassFourier(force_xyz(:,1),force_xyz(:,3),cutoffFreq));
    force_xyz(:,4) = real(lowpassFourier(force_xyz(:,1),force_xyz(:,4),cutoffFreq));
    
    
    % backing out subject body weight from x-forces
    force_net = force_xyz;
    force_net(:,2:4) = force_net(:,2:4) - orientation_grav.*subject_mass*9.807;
    
    % keeping force_og to crop each stride for COM data
    force_og = force_net;
    
    % getting acceleration on COM
    accel_xyz = force_net;
    accel_xyz(:,2:4) = accel_xyz(:,2:4)./subject_mass;
    
    
end

% if accel3 file, accel is read from user file, remove offset from
% appropriate columns
if strcmpi(filetype,'accel3')
    accel_xyz = readmatrix(fullfilename);
    accel_xyz = accel_xyz(~isnan(accel_xyz(:,1)),:);
    accel_xyz(:,2:4) = accel_xyz(:,2:4) - 2048.*orientation_grav;
    accel_xyz(:,2:4) = accel_xyz(:,2:4)./2048*9.807;    % now in [m/s2]
    
    % cropping data
    fig1 = figure(1);
    plot(accel_xyz(:,1), accel_xyz(:,3))
    title('Please crop useful data')
    [startstop,~] = ginput(2);
    startstop = sort(round2actual(accel_xyz(:,1),startstop));
    close(fig1)
    accel_xyz = accel_xyz(accel_xyz(:,1) >= startstop(1) &...
        accel_xyz(:,1) <= startstop(2),:);
    
    % marking strides using mediolateral data
    stride_times = selecting_strides(accel_xyz(:,[1,3]));
    num_strides = (length(stride_times)-3)/2+1;
    
    % removes offset by ensuring avg accel in all 3 dimensions = 0
    accel_xyz(:,2:4) = accel_xyz(:,2:4) - mean(accel_xyz(:,2:4));
end

% if positional data grabbing pos data, remove outliers, and smoothing
if contains(filetype,'positional') == 1
    posdata = readmatrix(fullfilename);
    posdata(:,1) = posdata(:,1)./fps_DLC;
    
    % removing outliers and smoothing position a bit
    posdata(:,8) = removing_outliers(posdata(:,8));
    posdata(:,8) = movmean(posdata(:,8),min(50,round(length(posdata(:,8))*0.03)));
    posdata(:,9) = removing_outliers(posdata(:,9));
    posdata(:,9) = movmean(posdata(:,9),min(50,round(length(posdata(:,9))*0.03)));
end

% if positional data, cropping posdata to fit within time bounds of the force/accel data
if strcmpi(filetype,'positional')
    posdata = posdata(posdata(:,1) >= accel_xyz(1,1) & ...
        posdata(:,1) <= accel_xyz(end,1),:);
    
    %plot(posdata(:,1),posdata(:,8),posdata(:,1),posdata_smooth)
    
    % velocity conversion factor
    top = mean(posdata(:,5));
    bot = mean(posdata(:,2));
    conversion = barlength / (top - bot);   % [m / DLC unit]
    
    % getting velocity from positional data using CFD and smoothing
    vel_xz_from_pos = zeros(length(posdata(:,8)),3);
    vel_xz_from_pos(:,1:2) = central_finite_difference(posdata(:,1),posdata(:,8),5);
    vel_xz_from_pos(:,[1,3]) = central_finite_difference(posdata(:,1),posdata(:,9),5);
    vel_xz_from_pos(:,2:3) = vel_xz_from_pos(:,2:3) .* conversion;
    vel_xz_from_pos(:,2:3) = movmean(vel_xz_from_pos(:,2:3),min(50,round(length(posdata(:,8))*0.03)));
    
    
    % plotting position vs velocity from pos
    figure(1)
    subplot(2,2,1)
    plot(posdata(:,1),posdata(:,8))
    subplot(2,2,3)
    plot(vel_xz_from_pos(:,1),vel_xz_from_pos(:,2))
    
    subplot(2,2,2)
    plot(posdata(:,1),posdata(:,9))
    subplot(2,2,4)
    plot(vel_xz_from_pos(:,1),vel_xz_from_pos(:,3))
end

% initializing values for step loop below
accel_og = accel_xyz;
tic
all_data_scaled = cell([1,num_strides]);
all_data_spaciotemp = cell([num_strides,1]);
enef_all = zeros(1,num_strides);

% for each stride, must find v, pos, energies, powers, perrec, and export
for k = 1:num_strides
    
    % getting current stride times
    if strcmpi(filetype,'accel3')
        stride_times_current = [stride_times(2*k-1),stride_times(2*k),stride_times(2*k+1)];
    elseif contains(filetype,'force')
        stride_times_current = [stride_times(k,1),stride_times(k,3),stride_times(k,2)];
    end
    accel_xyz = accel_og(accel_og(:,1) >= stride_times_current(1) & ...
        accel_og(:,1) <= stride_times_current(3),:);
    % getting velocity and position from acceleration data
    if ~strcmpi(filetype,'positional')
        % if strcmpi(filetype,'accel3') POSSIBLY NOT NEEDED
            % grabbing velocity offset only if first time through
            if k == 1 &&  strcmpi(filetype,'accel3')
                cd('Accel Data');
                velocity_file = readcell('Velocity.xlsx');
                event_name = extractBefore(filename,'.');
                velocity_offset = cell2mat(velocity_file(contains(velocity_file(:,1),event_name),3));
                subject_mass = cell2mat(velocity_file(contains(velocity_file(:,1),event_name),10));
                if isempty(velocity_offset)
                    error(strcat('Cannot find matching trial in Velocity.xlsx for: ',filename));
                end
                cd ..;
            elseif strcmpi(filetype, 'force8')
                velocity_offset = velocity(k);
            end
            
            velocity_offset = velocity_offset.*orientation_vel;
            
            % integrating acc to get vel, then adjusting for initial velocity
            % and offsets
            
            %detrend the acceleration:
            accel_trend = accel_xyz;
            accel_xyz(:,2:4) = detrend(accel_xyz(:,2:4));
            
            vel_xyz_from_force = numerical_integration(accel_xyz(:,1),accel_xyz(:,2:4));
            
            mean_true_xvel = mean(vel_xyz_from_force(:,2));
            mean_calc_xvel = velocity(k);
            offset = mean_calc_xvel-mean_true_xvel;
            mean_true_xvel = mean_true_xvel+offset;
            vel_xyz_from_force(:,2) = vel_xyz_from_force(:,2)+offset;

%             vel_xyz_from_force(:,2:4) = vel_xyz_from_force(:,2:4) + velocity_offset;
            velocity_mean = (~orientation_vel).*mean(vel_xyz_from_force(:,2:4));
%             vel_xyz_from_force(:,2:4) = vel_xyz_from_force(:,2:4) - velocity_mean;
%         else
%             vel_xyz_from_force = numerical_integration(accel_xyz(:,1),accel_xyz(:,2:4));
%             vel_xyz_from_force(:,2) = vel_xyz_from_force(:,2) - mean(vel_xz_from_pos(:,2));             % assume avg v from force and from pos must be same
%             vel_xyz_from_force(:,3:4) = vel_xyz_from_force(:,3:4) - mean(vel_xyz_from_force(:,3:4));    % assume avg Vy and Vz must be zero
%         end
        
        % integrating vel to get position
        pos_xyz_from_force = numerical_integration(vel_xyz_from_force(:,1),vel_xyz_from_force(:,2:4));
        
        % plots Ax, Vx, and Sx vs Time
        figure(k)
        subplot(4,1,1)
        plot(accel_xyz(:,1),accel_xyz(:,2))
        title("Accel X [m/s^2]")
        subplot(4,1,2)
        plot(vel_xyz_from_force(:,1),vel_xyz_from_force(:,2))
        title('Vel X [m/s]')
        subplot(4,1,3)
        plot(pos_xyz_from_force(:,1),pos_xyz_from_force(:,2))
        title('Pos X [m]')
        
        subplot(4,1,4)
        plot(pos_xyz_from_force(:,1),pos_xyz_from_force(:,4))
        title('Pos Z [m]')
    end
    
    % grabbing the relevant velocity and position data
    if ~strcmpi(filetype,'positional')
        Vxyz = vel_xyz_from_force;
        Pxyz = pos_xyz_from_force;
    else
        Vxyz(:,[1,2,4]) = vel_xz_from_pos;
        Vxyz(:,3) = zeros(length(Vxyz),1);  % if pos data, we don't have y data
        Pxyz = [posdata(:,1),posdata(:,8),...
            zeros(length(posdata(:,9)),1),posdata(:,9)];
        Pxyz(:,2:4) = Pxyz(:,2:4).*conversion;
    end
    
    % calculating percent recovery normally and % rec without KEx, KEy, KEz
    % if this is a horizontal trial
    perrec = zeros(1,4);
    for j = 1:4
        % calculating kinetic energies [J/kg] (all y data is zero if only have DLC data)
        EKxyz = [Vxyz(:,1), 0.5*Vxyz(:,2:4).^2];
        if j <= 3
            EKtot = [Vxyz(:,1), EKxyz(:,2) + EKxyz(:,3) + EKxyz(:,4) - EKxyz(:,j+1)];
        else
            EKtot = [Vxyz(:,1), EKxyz(:,2) + EKxyz(:,3) + EKxyz(:,4)];
        end
        % calculating potential energy [J/kg]
        EP = [Pxyz(:,1), 9.807.*(Pxyz(:,2)-Pxyz(1,2))];
        
        % total mechanical energy [J/kg]
        Etot = [Vxyz(:,1), EKtot(:,2) + EP(:,2)];
        
        % total horizontal and mechanical energies [J/kg]
        Ev = [Vxyz(:,1), EP(:,2) + EKxyz(:,2)];
        Eh = [Vxyz(:,1), EKxyz(:,3) + EKxyz(:,4)];
        
        % calculating percent recovery
        deltaET = zeros(length(Etot),1);
        for i = 2:length(Etot)
            deltaET(i-1) = Etot(i,2)-Etot(i-1,2);
        end
        deltaET(deltaET<0) = 0;
        
        deltaEK = zeros(length(EKtot),1);
        for i = 2:length(EKtot)
            deltaEK(i-1) = EKtot(i,2)-EKtot(i-1,2);
        end
        deltaEK(deltaEK<0) = 0;
        
        deltaEP = zeros(length(EP),1);
        for i = 2:length(EP)
            deltaEP(i-1) = EP(i,2)-EP(i-1,2);
        end
        deltaEP(deltaEP<0) = 0;
        
        ETs = sum(deltaET);
        EKs = sum(deltaEK);
        EPs = sum(deltaEP);
        
        % percent recovery [w/o X, w/o Y, w/o Z, with all]
        perrec(j) = ((EKs + EPs - ETs)/(EKs + EPs))*100;
    end
    
    % calculating climbing efficiency of this step
    enef_force = [accel_xyz(:,1), abs(subject_mass * (accel_xyz(:,2)+9.807)), ...
        abs(subject_mass*accel_xyz(:,3)),abs(subject_mass*accel_xyz(:,4))];  % fxyz = [m(ax+g),m*ay,m*az]
    enef_force(end,:) = [];
    enef_disp = zeros(length(enef_force),4);
    enef_disp(:,1) = enef_force(:,1);
    for col = 2:4
        for i = 1:length(enef_force)
            enef_disp(i,col) = abs(pos_xyz_from_force(i+1,col) - pos_xyz_from_force(i,col));
        end
    end
    enef_energy_in_xyz = [transpose(enef_force(:,2))*enef_disp(:,2),...
        transpose(enef_force(:,3))*enef_disp(:,3), ...
        transpose(enef_force(:,4))*enef_disp(:,4)];
    enef_energy_in = sum(enef_energy_in_xyz);
    enef_energy_out = subject_mass * (9.807*pos_xyz_from_force(end,2) + 0.5 * vel_xyz_from_force(end,2)^2 ...
        + 0.5 * vel_xyz_from_force(end,3)^2 + 0.5 * vel_xyz_from_force(end,4)^2) - ...
        subject_mass * (9.807*pos_xyz_from_force(1,2) + 0.5 * vel_xyz_from_force(1,2)^2 ...
        + 0.5 * vel_xyz_from_force(1,3)^2 + 0.5 * vel_xyz_from_force(1,4)^2);
    
    enef_all(k) = enef_energy_out / enef_energy_in * 100;
    
    % calculating all relevant power info by doing P = dE/dt
    PKxyz = central_finite_difference(EKxyz(:,1),EKxyz(:,2:end),1);
    PKtot = central_finite_difference(EKtot(:,1),EKtot(:,2),1);
    PP = central_finite_difference(EP(:,1),EP(:,2),1);
    Ptot = central_finite_difference(Etot(:,1),Etot(:,2),1);
    
    % calculating sum of positive power increments [Px,Py,Pz,PK,PP,Ptot]
    power_matrix = [EKxyz(:,2:4),EKtot(:,2),EP(:,2),Etot(:,2)];
    sum_power = zeros(size(power_matrix));
    for j = 1:size(power_matrix,2)
        for i = 2:size(power_matrix,1)
            sum_power(i,j) = power_matrix(i,j) - power_matrix(i-1,j);
        end
    end
    sum_power(sum_power < 0) = 0;
    sum_power = sum(sum_power)/(EKxyz(end,1)-EKxyz(2,1));
    
    % NaN force values for pos and acc data, cropping data for force trials 
    if strcmpi(filetype,'positional') || strcmpi(filetype,'accel3')
        force_xyz = NaN(length(Vxyz),4);
        if strcmpi(filetype,'positional')
            accel_xyz = NaN(length(Vxyz),4);
        end
    else        
        force_xyz = force_og(force_og(:,1) >= accel_xyz(1,1) & ...
            force_og(:,1) <= accel_xyz(end,1),:);
    end
    
    % putting together COM data
    COMdata = [Vxyz(:,1),force_xyz(:,2:4),accel_xyz(:,2:4),Vxyz(:,2:4),...
        Pxyz(:,2:4),EKxyz(:,2:4),EKtot(:,2),EP(:,2),Eh(:,2),Ev(:,2),Etot(:,2),...
        PKxyz(:,2:4),PKtot(:,2),PP(:,2),Ptot(:,2)];
    if contains(filetype,'force')
        COMlabels = {'Time (s)','Fx (N)', 'Fy (N)', 'Fz (N)',...
            'Ax (m s^-2)',  'Ay (m s^-2)',  'Az (m s^-2)',...
            'Vx (m s^-1)',  'Vy (m s^-1)',  'Vz (m s^-1)',...
            'Pos X (m)', 'Pos Y (m)', 'Pos Z (m)',...
            'EK X (J)', 'EK Y (J)', 'EK Z (J)', 'EK Total (J)',...
            'EP (J)', 'EHorizontal (J)', 'EVertical (J)','E Total (J)',...
            'PK X (W)','PK Y (W)','PK Z (W)','PK Total (W)','PP (W)','P Total (W)'};
    elseif strcmpi(filetype,'accel3')
        COMlabels = {'Time (s)','Fx (N)', 'Fy (N)', 'Fz (N)',...
            'Ax (m s^-2)',  'Ay (m s^-2)',  'Az (m s^-2)',...
            'Vx (m s^-1)',  'Vy (m s^-1)',  'Vz (m s^-1)',...
            'Pos X (m)', 'Pos Y (m)', 'Pos Z (m)',...
            'EK X (J kg^-1)', 'EK Y (J kg^-1)', 'EK Z (J kg^-1)', 'EK Total (J kg^-1)',...
            'EP (J kg^-1)', 'EHorizontal (J kg^-1)', 'EVertical (J kg^-1)','E Total (J kg^-1)',...
            'PK X (W kg^-1)','PK Y (W kg^-1)','PK Z (W kg^-1)','PK Total (W kg^-1)','PP (W kg^-1)','P Total (W kg^-1)'};
    else
        error('check filetype')
    end
    COMexport = [COMlabels;num2cell(COMdata)];
    
    % scaling data
    COMdata_scaled = scaling_data(COMdata,100);
    all_data_scaled{k} = COMdata_scaled;
    COMexport_scaled = [COMlabels;num2cell(COMdata_scaled)];
    
    % spaciotemp data
    if contains(filetype,'force')
        if orientation_grav(3) == 1         % if horizontal
            spaciotemp_labels = {'Filename','Stride No.','Orientation','Substrate',...
                'Individual','Mass [kg]','Avg Foreaft V [m/s]','Recovery (no X) [%]',...
                'Recovery (no Y) [%]','Recovery (no Z) [%]','Recovery Total [%]',...
                'Phase [deg]','Duty Factor','Stance Time [s]','Swing Time [s]','Stride Time [s]',...
                'Stride Freq [Hz]','Stride Length [m]','sumPx [W]','sumPy [W]','sumPz [W]',...
                'sumPK [W]','sumPP [W]','sumPT [W]'};
        else                                % if vertical
            spaciotemp_labels = {'Filename','Stride No.','Orientation','Substrate',...
                'Individual','Mass [kg]','Avg Foreaft V [m/s]','Wout [J]','Win [J]','Energy Efficiency [%]',...
                'Phase [deg]','Duty Factor','Stance Time [s]','Swing Time [s]','Stride Time [s]',...
                'Stride Freq [Hz]','Stride Length [m]','sumPx [W]','sumPy [W]','sumPz [W]',...
                'sumPK [W]','sumPP [W]','sumPT [W]','sumEK [J]','sumEP [J]','sumET [J]'};
        end
    elseif strcmpi(filetype,'accel3')
        if orientation_grav(3) == 1         % if horizontal
            spaciotemp_labels = {'Filename','Stride No.','Orientation','Substrate',...
                'Individual','Mass [kg]','Avg Foreaft V [m/s]','Recovery (no X) [%]',...
                'Recovery (no Y) [%]','Recovery (no Z) [%]','Recovery Total [%]',...
                'Phase [deg]','Duty Factor','Stance Time [s]','Swing Time [s]','Stride Time [s]',...
                'Stride Freq [Hz]','Stride Length [m]','sumPx [W/kg]','sumPy [W/kg]','sumPz [W/kg]',...
                'sumPK [W/kg]','sumPP [W/kg]','sumPT [W/kg]'};
        else                                % if vertical
            spaciotemp_labels = {'Filename','Stride No.','Orientation','Substrate',...
                'Individual','Mass [kg]','Avg Foreaft V [m/s]','Wout [J/kg]','Win [J/kg]','Energy Efficiency [%]',...
                'Phase [deg]','Duty Factor','Stance Time [s]','Swing Time [s]','Stride Time [s]',...
                'Stride Freq [Hz]','Stride Length [m]','sumPx [W/kg]','sumPy [W/kg]','sumPz [W/kg]',...
                'sumPK [W/kg]','sumPP [W/kg]','sumPT [W/kg]','sumEK [J/kg]','sumEP [J/kg]','sumET [J/kg]'};
        end
    else
        error('filetype not yet supported')
    end
    phase = {'N/A'};
    stance_time = {stride_times_current(2) - stride_times_current(1)};
    swing_time = {stride_times_current(3) - stride_times_current(2)};
    stride_time = {stride_times_current(3) - stride_times_current(1)};
    duty_factor = {stance_time{1}/stride_time{1} * 100};
    stride_freq = {1/stride_time{1}};
    velocity_foreaft = {velocity_offset(logical(orientation_vel))};
    stride_length = {stride_time{1} * velocity_foreaft{1}};
    orientation_cell = orientation_grav_labels{logical(orientation_grav)};
    if strcmpi(filetype,'accel3')
        individual = {extractBefore(filename,'_')};
        enef_energy_in = enef_energy_in / subject_mass;
        enef_energy_out = enef_energy_out / subject_mass;
        EKs = EKs / subject_mass;
        EPs = EPs / subject_mass;
        ETs = ETs / subject_mass;
    elseif contains(filetype,'force8')
        individual = stride_data(k,6);
    else
        error('check me out')
    end
    
    if orientation_grav(3) == 1         % if horizontal
        spaciotemp_data = [{filename},k,{orientation_cell},{substrate},...
            individual,subject_mass,velocity_foreaft,num2cell(perrec),...
            phase,duty_factor,stance_time,swing_time,stride_time,...
            stride_freq,stride_length,num2cell(sum_power)];
    else                                % if vertical
        spaciotemp_data = [{filename},k,{orientation_cell},{substrate},...
            {extractBefore(filename,'_')},subject_mass,velocity_foreaft,...
            {enef_energy_out},{enef_energy_in},{enef_all(k)},...
            phase,duty_factor,stance_time,swing_time,stride_time,...
            stride_freq,stride_length,num2cell(sum_power),EKs,EPs,ETs];
    end
    
    spaciotemp_export = [spaciotemp_labels;spaciotemp_data];
    all_data_spaciotemp{k,1} = spaciotemp_data;
        
        % exporting data
        % goes into proper COM Data folder
        if strcmpi(filetype,'accel3')
            cd('Accel Data')
        elseif strcmpi(filetype,'force8')
            cd('Force8')
        else
            error('need to add a case for other file types')
        end
        cd('COM Data')
        
        % if file exists prior to trial, overwrite that file completely
        if k == 1
            files_to_delete = dir;
            files_to_delete = {files_to_delete.name};
            for a = 1:length(files_to_delete)
                if contains(files_to_delete{a},event_name)
                    delete(files_to_delete{a})
                end
            end
        end
        
        % filenames
        COM_filename = strcat(cd,filesep,event_name,'_COM.xlsx');
        scaled_filename = strcat(cd,filesep,event_name,'_Scaled.xlsx');
        spaciotemp_filename = strcat(cd,filesep,event_name,'_Spaciotemp.xlsx');
        
        % write cell
        if write_data
            sheet_name = strcat('Stride_',num2str(k));
            writecell(COMexport,COM_filename,'Sheet',sheet_name);
            writecell(COMexport_scaled,scaled_filename,'Sheet',sheet_name);
            writecell(spaciotemp_export,spaciotemp_filename,'Sheet',sheet_name);
        end
        cd ..; cd ..;
        
end


% if _Scaled_Data & _Spaciotemp are not files, creates a file and fills with labels
if strcmpi(filetype,'accel3')
    cd('Accel Data')
elseif strcmpi(filetype,'force8')
    cd('Force8')
else
    error(strcat('Add new condition for this filetype:',filetype));
end

cd('COM Data')
if ~isfile('_Scaled_Data.xlsx')
    for z = 2:length(COMlabels)
        labels = [{'Event'}, num2cell(1:101)];
        %             scaled_data = [{strcat(event_name,'_Stride_',num2str(k))},num2cell(transpose(COMdata_scaled(:,z)))];
        %             all_data_scaled{k} = scaled_data;
        %             scaled_export = [labels;scaled_data];
        if write_data
            writecell(labels,'_Scaled_Data.xlsx','Sheet',COMlabels{z},...
                'AutoFitWidth',false,'UseExcel',true)
        end
    end
end
if ~isfile('_Spaciotemp_Data.xlsx') && write_data
    writecell(spaciotemp_labels,'_Spaciotemp_Data.xlsx','AutoFitWidth',false)
end

% writing to _Scaled_Data.xlsx
for j = 2:length(COMlabels)
    scaled_data = zeros(101,num_strides);
    scaled_labels = cell(1,num_strides);
    for i = 1:num_strides
        % building export matrices/cells so that I only export once per sheet
        scaled_data(:,i) = all_data_scaled{i}(:,j);
        scaled_labels{i} = strcat(event_name,'_Stride_',num2str(i));
    end
    
    scaled_export = transpose([scaled_labels;num2cell(scaled_data)]);
    
    if write_data
        writecell(scaled_export,'_Scaled_Data.xlsx','Sheet',COMlabels{j},...
            'WriteMode','append','AutoFitWidth',false,'UseExcel',true);
    end
    
    percent_done = round((j-1)/length(COMlabels)*100);
    disp(strcat('Aproximately ',num2str(percent_done),'% Complete'));
end

% writing to _Spaciotemp_Data.xlsx
all_data_spaciotemp = vertcat(all_data_spaciotemp{:,1});
if write_data
    writecell(all_data_spaciotemp,'_Spaciotemp_Data.xlsx',...
        'WriteMode','append','AutoFitWidth',false,'UseExcel',true);
end

clc
disp("All Done");
toc

% FUNCTIONS
function [coordAct] = round2actual(data, coordGinput)
% this function takes in data, which should be a column vector or matrix
% whos height is longer than its width and returns the value(s) in each
% column that is closest to coordGinput)


[mdata, ndata] = size(data);
lengGinput = length(coordGinput);
diff = zeros(mdata, ndata);
coordAct = zeros(lengGinput,1);

for j = 1:lengGinput
    for n = 1:ndata
        for m = 1:mdata
            diff(m,n) = abs(data(m,n) - coordGinput(j));
        end
    end
    minimum = min(diff);
    location = find(minimum == diff);
    if length(location) > 1
        location = location(1);
    end
    coordAct(j) = data(location,1);
end
coordAct = sort(coordAct);              % sorts values in ascending order

end
function [ForceXYZ,cal_x] = Raw8ChannelData_to_ForceXYZ(fullfilename)
    % this function will take the file passed to it and output XYZ forces


    numbersFP = readmatrix(fullfilename);       % force plate data
    
    numbersS = [0.0025060 -0.0051789 -0.0029982 0.0017707 0.0000271 -0.0532926 0.0531128 0.0000278;
    0.0002548 -0.0016677 -0.0019367 0.0043635 0.0532233 0.0001949 -0.0001942 0.0546583;
    -0.0780750 -0.0758160 -0.0752487 -0.0768511 0.0003307 0.0012481 -0.0001544 0.0000142]; % calibration data
    
    fig1 = figure(1);                                   % creates figure 1
    plot(numbersFP(:,2), numbersFP(:,3));       % plots the vectors that we grabbed from numbersFP
    grid on                                     % turns on grid
    title({'Please click two points between which the values',...
        'will be used for calibration'})
    
    %goodtime = false;       % value is only set to true once the user selects points with a low std dev.
    counter2 = 0;           % counter variable, used to change the instructions for graph input if user messes up
    
    %while goodtime == false     % while loop repeats until std dev is under specified threshold
    if counter2 == 0
        disp('Please click two coordinates on the graph, between which the y-values will be averaged to calibrate data.');
    end
    [x,~] = ginput(2);                  % user clicks two points on graph, x and y coordinates are saved
    answer = transpose(x);              % grabs x-coordinates and puts them in 1x2 matrix
    
    if answer(2) < answer(1)            % if the user clicks from left to right, we must switch the order of the answer vector
        placeholder = answer(1);        % saves the first value in a placeholder variable
        answer(1) = answer(2);          % overwrites the first value with the second
        answer(2) = placeholder;        % overwrites the second value with the placeholder variable
    end
    
    [rows, ~] = size(numbersFP);  % captures no. of rows and columns in numbersFP matrix
    
    for i = 1:rows                      % finds row containing lower desired t-value
        if numbersFP(i,2) < answer(1)   % by recording and overwriting the row number until the time value
            minX = numbersFP(i+1,1);    % check is larger than the desired value
        end
        if numbersFP(i,2) < answer(2)   % finds row containing higher t-value by recording and overwriting
            maxX = numbersFP(i+1,1);    % the row number until t-value checked is higher than desire value
        end
    end
    
    offloading = zeros(1,8);            % creating the array for offloading values
    offloadingstd = zeros(1,8);         % creating array for std dev of values / average of values
    
    for i = 1:8
        offloading(1,i) = mean(numbersFP(minX:maxX,i+2));   % assigns the average of the values to the ith position in 'offloading'
        offloadingstd(1,i) = std(numbersFP(minX:maxX,i+2))/abs(offloading(1,i));    % calcs std dev of numbers used to find the avg
    end
    
    %avgstd = mean(offloadingstd);       % takes the average of the std dev array
    
%     if avgstd < 0.5                     % if avg of std dev is acceptable low, the goodtime variable is changed to true
%         goodtime = true;                % and the code exits the while loop
%     elseif avgstd >= 0.5                % if avg of std dev is too high, the goodtime variable is kept as false
%         fprintf('\n')                   % and code asks you to choose different values
%         disp('Please try again, standard deviation of y-values is too high')
%         counter2 = 1;
%     end
        
    %end
    close(fig1)
    clc                             % clears command window of error messages
    delta = zeros(rows,8);          % creates the delta vector
    
    for i = 1:8
        delta(:,i) = numbersFP(:, i+2) - offloading(1,i);    % delta = numbersFP - offloading
    end
    
    force = delta*transpose(numbersS);  % matrix multiplication 10010x8 * 8x3 = 10010x3
    [ForceXYZ] = force*4.4;        % gives final force in units of Newton
    
    ForceXYZ = [numbersFP(:,2), ForceXYZ];
    cal_x = x;

end
function newForceXYZ = offset_correction(ForceXYZ,t_bounds)
% this function will correct x, y and z arrays for any offset by
% subtracting the average value of each column between the t_bounds.  The
% t_bounds should come from the cal_x_bounds output of the
% raw8channeldata_to_forcexyz function

if size(ForceXYZ,2) == 4
    force_crop = ForceXYZ(ForceXYZ(:,1) >= t_bounds(1) & ...
        ForceXYZ(:,1) <= t_bounds(2),2:4);
    
    ForceXYZ(:,2:4) = ForceXYZ(:,2:4) - mean(force_crop);
    newForceXYZ = ForceXYZ;
else
    error('input vector must be four columns: t,x,y,z')
end
end
function newForceXYZ = drift_correction(ForceXYZ,x_bounds)
% this function takes in a four column matrix [t,x,y,z] and removes drift
% by ensuring the slope between the 2 x_bounds is zero.
% x_bounds should be the output 'cal_x_bounds' from the function
% Raw8ChannelData_to_ForceXYZ

if max(size(x_bounds)) == 2 && min(size(x_bounds)) == 1
    
    if size(ForceXYZ,2) == 4
        
        X = round2actual(ForceXYZ(:,1),x_bounds);
        
        for j = 2:4
            Y = ForceXYZ(ForceXYZ(:,1) == X(1) | ForceXYZ(:,1) == X(2),j);
            
            rise = Y(2)-Y(1);
            run = X(2)-X(1);
            slope=rise/run;
            
            % new drift code
            subtract = zeros(length(ForceXYZ),1);
            for i = 1:length(ForceXYZ)
                subtract(i) = slope*(ForceXYZ(i,1)-ForceXYZ(1,1));
            end
            ForceXYZ(:,j) = ForceXYZ(:,j) - subtract;
        end
        newForceXYZ = ForceXYZ;
        
    else
        error('input must be 4 columns')
    end
else
    error('input x_bounds must be 2x1 or 1x2')
end
end
function newData = removing_outliers(oldData,timevar)

% READ ME:
% this function takes in either 1 or 2 1-dimensional vectors and removes
% the outliers of the first input vector, replacing them with values
% calculated using linear interpolation.  The 2nd input is optional and is
% only needed if the data in the first vector is not spaced evenly across
% time (or whatever dependent variable is being used)

% the user will click two points which form opposite vertices of a
% rectangle, anything inside will be deleted and interpolated.  When the
% user is done they must click two points outside the domain of data.

% the output of this function is only the 1-dimensional vector that had
% oiutliers removed

% note: bounds of box when cropping must be within the domain of points
% because if not we do not have a y value to interpolate with


[m,n] = size(oldData);

if ~exist('timevar','var')
    index_num = transpose(1:m);
    newData = [index_num,oldData];
else
    [m2,n2] = size(timevar);
    if ~(m2 == m && n2 == n)
        error("Dimensions of input vectors into 'removing_outliers' do not match");
    end
    newData = [timevar,oldData];
end

if n == 1
    
    % code will continue cropping bounds until user
    % double clicks outside of x domain
    counterPlot = 1;
    tbounds = [1,1];
    
    
    while all([tbounds <= newData(end,1), tbounds >= 1])
        
        % plot original data if first time through
        if counterPlot == 1
            fig1 = figure(1);
            plot(newData(:,1),newData(:,2),'b');
            title('click to remove outliers')
            hold on
        end
        
        % grab bounds of box within which we will remove outliers
        [tbounds,ybounds] = ginput(2);
        tbounds = sort(tbounds);
        ybounds = sort(ybounds);
        
        % overwriting data with Nan if it's inside the box created by user
        for i = 1:m
            if newData(i,1) >= tbounds(1) && newData(i,1) <= tbounds(2) ...
                    && newData(i,2) >= ybounds(1) && newData(i,2) <= ybounds(2)
                newData(i,2) = NaN;
            end
        end
        
        % if data is NaN, data gets overwritten by linear interpolation
        for i = 1:m
            if isnan(newData(i,2))
                x1 = newData(i-1,1);    % takes most recent non-NaN value
                y1 = newData(i-1,2);
                
                j = 0;
                while isnan(newData(i+j,2))
                    j = j + 1;
                end
                
                x2 = newData(i+j,1);    % finds next non-NaN value
                y2 = newData(i+j,2);
                xcurrent = newData(i,1);
                
                newData(i,2) = y1 + (y2-y1)/(x2-x1)*(xcurrent-x1);  % linear interpolation
                
            end
        end
        % plotting new data
        if counterPlot > 1
            delete(newDataPlot)
        end
        newDataPlot = plot(newData(:,1),newData(:,2),'r');
        % location of legend so it doesn't block data
        if abs(newData(end,2) - max(newData(:,2))) <= ...
                abs(newData(1,2) - max(newData(:,2)))
            legend('Original Data', 'With Outliers Removed','Location','northwest');
        else
            legend('Original Data', 'With Outliers Removed','Location','northeast');
        end
        counterPlot = counterPlot + 1;
    end
    close(fig1)
    newData = newData(:,2);
end
end
function [ y ] = lowpassFourier( t, x, fcutoff )

X = fft(x);
N = length(t);
f = 1 / ( t(2)-t(1) );

for i = 0:N/2
    if ( f * i/N > fcutoff)
        X(i) = 0;
        X(N - i) = 0;
    end
end

y = ifft(X);
end
function time_and_vel = central_finite_difference(time,position,width)
% this function will take in a time and position matrix of the same length
% and will perform a numerical derivative on them for each column of position.  
% For the first row it will do FFD and for the last point it will do BFD, 
% both with a step of 1.

% for all other points it will perform CFD with a width input by the user,
% however it will use a smaller width near the ends of the vector as
% necessary

% ex: first point: FFD: (pos(2)-pos(1))/(time(2)-time(1)
% ex: 3rd point: CFD: (pos(5)-pos(1))/(time(5)-time(1)
% ex: middle points: CFD: (pos(i+wid)-pos(i-wid))/(time(i+wid)-time(i-wid))
% ex: last point: BFD: (pos(m)-pos(m-1))/(time(m)-time(m-1))

% checks
if ~exist('width','var')
    width = 1;
end
if width < 1
    width = 1;
end
if ~(mod(width,1) == 0)
    error('width must be an integer values >= 1')
end

[m,n] = size(time);
[m2,n2] = size(position);

if ~n==1
    error('time input must be a column vector')
elseif ~(m == m2)
    error('both inputs must be the same length')
end

% math using central finite difference

%creating vel array
time_and_vel = zeros(m,n2+1);
time_and_vel(:,1) = time(:,1);

% does forward finite dif for first value and backward finite dif for last
for i = 2:n2+1
    time_and_vel(1,i) = (position(2,i-1)-position(1,i-1))/(time(2)-time(1));
    time_and_vel(m,i) = (position(m,i-1)-position(m-1,i-1))/(time(m)-time(m-1));
end
% does cfd for rest of values
for j = 2:n2+1
    for i = 2:m-1
        % uses smaller width when near the ends of the vector
        width_min = min([width,i-1,m-i]);
        
        time_and_vel(i,j) = (position(i+width_min,j-1) - position(i-width_min,j-1))/...
            (time(i+width_min) - time(i-width_min));
    end
end
end
function integration_result = numerical_integration(time,values)
% this function will take in a 1-dimensional 'time' column vector and a
% vector or matrix (whose length must equal the time vectors) and perform
% integration on each column with respect to the time vector

% checks
[m,n] = size(time);
[m2,n2] = size(values);

if ~n==1
    error('time input must be a column vector')
elseif ~(m == m2)
    error('both inputs must be the same length')
end

% math
integration_result = zeros(m,n2+1);
integration_result(:,1) = time;

for j = 2:n2+1
    for i = 2:m
        integration_result(i,j) = trapz(time(1:i),values(1:i,j-1));

    end
end
end
function scaled_data = scaling_data(data,num_of_points)
% this function will take in an m x n data matrix and scale the data in
% each column to the desired number of points using 1D interpolation

% checks
if ~(mod(num_of_points,1) == 0)
    error('please enter an integer value for num_of_points in scaling_data()')
end

[m,n] = size(data);
scaled_data = zeros(num_of_points+1,n);

% math
for i = 1:n
    matlength = m - 1;                          %
    scale = num_of_points/matlength;            %
    xnew = transpose(0:scale:100);              % scaled data code
    x = transpose(0:1:100);                     %
    scaled_data(:,i) = interp1(xnew,data(:,i),x);
end
end
function points = selecting_strides(data)
% this function will take in a m-by-2 matrix and spit out the points given
% by the user.  The user will select points in the stride to determine
% stride time and duty factor and will then click outside of the
% domain of points to signify they are done

% checks
[~,n] = size(data);
if n ~= 2
    error('The size of the input matrix must be m-by-2')
end

% plotting data
fig1 = figure(1);
scatter(data(:,1),data(:,2),'.','b');
title({'Select the points for the strides'},{'Click outside the domain of points to signify completion'})

% initializing
tbound = mean(data(:,1));
points = NaN(100,1);
counter = 1;

% grabs ginput from user until they click outside the domain of values
while tbound > min(data(:,1)) && tbound < max(data(:,1))
            [tbound,~] = ginput(1);
            points(counter) = tbound;
            counter = counter + 1;
end
close(fig1)

% removes NaN values
points = points(~isnan(points));

% checks length to make sure user did it correctly
if rem(length(points),2) == 1
    error('You mucked up the selection process, please ensure an odd number of points are picked before selecting outside the data');
end

% deletes the last point since that's the point outside the domain
points(size(points,1)) = [];

if ~all(points == sort(points))
    error('Please select the points in increasing order')
end
end

