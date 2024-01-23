% s_folder = 'D:\Carl Malings\CAPS MATLAB\CAPS Satellite Data Analysis\Satellite MCD19A2';
% s_folder = 'D:\Downloads\Satellite 2a';
% s_folder = 'D:\Carl Malings\CAPS MATLAB\CAPS Satellite Data Analysis\Rwanda MCD19A2 Apr-May 2018';
% s_folder = 'D:\Carl Malings\CAPS MATLAB\CAPS Satellite Data Analysis\Nigeria MCD19A2 Jun-Sep 2019';
% s_folder = 'D:\Carl Malings\CAPS MATLAB\CAPS Satellite Data Analysis\Pittsburgh MCD19A2 Jan-Jun 2019';
% s_folder = 'D:\Carl Malings\CAPS MATLAB\CAPS Satellite Data Analysis\Rwanda MCD19A2 July-November 2017';
% s_folder = 'D:\Carl Malings\CAPS MATLAB\CAPS Satellite Data Analysis\Rwanda MCD19A2 Apr2017-May2018';
% s_folder = 'D:\Carl Malings\CAPS MATLAB\CAPS Satellite Data Analysis\Pittsburgh MCD19A2 Jan-Dec 2018';
% s_folder = 'D:\Carl Malings\CAPS MATLAB\CAPS Satellite Data Analysis\Malawi MCD19A2 Aug 2017';
% s_folder = 'D:\Carl Malings\CAPS MATLAB\CAPS Satellite Data Analysis\Kinshasa MCD19A2 Mar2018-Oct2019';
% s_folder = 'D:\Carl Malings\CAPS MATLAB\CAPS Satellite Data Analysis\Malawi MCD19A2 Jun2017-Jul2018';
% s_folder = 'D:\Carl Malings\CAPS MATLAB\CAPS Satellite Data Analysis\Kampala and Addis Ababa MCD19A2 Jan-Dec 2019';
% s_folder = 'D:\Carl Malings\CAPS MATLAB\CAPS Satellite Data Analysis\Dhaka MCD19A2 Jan-Feb 2018';
% s_folder = 'D:\Carl Malings\CAPS MATLAB\CAPS Satellite Data Analysis\Dhaka MCD19A2 Mar-Apr 2018';
% s_folder = 'D:\Carl Malings\CAPS MATLAB\CAPS Satellite Data Analysis\Dhaka MCD19A2 May-Jun 2018';
% s_folder = 'D:\Carl Malings\CAPS MATLAB\CAPS Satellite Data Analysis\Dhaka MCD19A2 Jul-Aug 2018';
% s_folder = 'D:\Carl Malings\CAPS MATLAB\CAPS Satellite Data Analysis\Dhaka MCD19A2 Sep-Oct 2018';
% s_folder = 'D:\Carl Malings\CAPS MATLAB\CAPS Satellite Data Analysis\Dhaka MCD19A2 Nov-Dec 2018';
% s_folder = 'D:\Carl Malings\CAPS MATLAB\CAPS Satellite Data Analysis\Dhaka MCD19A2 Jul-Dec 2018';
% s_folder = 'D:\Carl Malings\CAPS MATLAB\CAPS Satellite Data Analysis\Dhaka MCD19A2 Jan-Jun 2018';
% s_folder = 'D:\Carl Malings\CAPS MATLAB\CAPS Satellite Data Analysis\Dhaka MCD19A2 Jan-Apr 2018';
% s_folder = 'D:\Carl Malings\CAPS MATLAB\CAPS Satellite Data Analysis\Dhaka MCD19A2 May-Aug 2018';
% s_folder = 'D:\Carl Malings\CAPS MATLAB\CAPS Satellite Data Analysis\Dhaka MCD19A2 Sep-Dec 2018';
%s_folder = 'D:\Carl Malings\CAPS MATLAB\CAPS Satellite Data Analysis\Dhaka MCD19A2 Jan-Dec 2017';
% s_folder = 'D:\Carl Malings\CAPS MATLAB\CAPS Satellite Data Analysis\Dhaka MCD19A2 Jan-Dec 2019';
%s_folder = 'D:\Carl Malings\CAPS MATLAB\CAPS Satellite Data Analysis\Dhaka MCD19A2 Jan-Dec 2013';
s_folder = 'D:\Carl Malings\CAPS MATLAB\CAPS Satellite Data Analysis\Dhaka MCD19A2_h25 Jan-Dec 2013';             %h25 files


addpath(fullfile(CAPS_Master_Definitions('Main Directory'),'CAPS Satellite Data Analysis'))

s_grid_to_load = 'grid1km';

%t_grid_info = importfile_sn10deg(fullfile(CAPS_Master_Definitions('Main Directory'),'CAPS Satellite Data Analysis','sn_bound_10deg.txt'));
t_grid_info = importfile_sn10deg_2(fullfile(CAPS_Master_Definitions('Main Directory'),'CAPS Satellite Data Analysis','sn_gring_10deg.txt'));

% Corners of v4 h11:
% 0,0:          49.995833,-108.884749
% 0,1199:       49.995833,-93.341822
% 1199,0:       40.004167,-91.378647
% 1199,1199:    40.004167,-78.334657

RESTRICT_SIZE = 0;
% % Rwanda:
% v_lat_lims = [-3 -1];
% v_lon_lims = [28 31];
% v_lat_lims = [-2.2 -1.4];
% v_lon_lims = [29.4 30.5];

% % Nigeria:
% v_lat_lims = [3 14];
% v_lon_lims = [2 15];

% % Pittsburgh:
% v_lat_lims = [40.1 40.8];
% v_lon_lims = [-80.5 -79.7];

% % Malawi:
% v_lat_lims = [-16.2 -14.0];
% v_lon_lims = [33.6 35.7];

% Kinshasa:
% v_lat_lims = [-4.4 -4.2];
% v_lon_lims = [15.2 15.4];

% % Kampala (0.3,32.591):
% v_lat_lims = [0.2 0.4];
% v_lon_lims = [32.5 32.7];

% Addis Ababa (9.058,38.761):
% v_lat_lims = [8.9 9.2];
% v_lon_lims = [38.6 38.9];

% Dhaka (23.72,90.41):
% v_lat_lims = [23.60 23.90];
% v_lon_lims = [90.30 90.60];

%Bangladesh
v_lat_lims = [20.67 27.3];
v_lon_lims = [88.00 89.68];

% s_savename = 'Rwanda_MCD19A2_Apr_May_2018.mat';
% s_savename = 'Nigeria_MCD19A2_Jun_Sep_2019.mat';
% s_savename = 'Pittsburgh_MCD19A2_Jan_Jun_2019.mat';
% s_savename = 'Rwanda_MCD19A2_Apr2017_May2018.mat';
% s_savename = 'Pittsburgh_MCD19A2_Jan_Dec_2018.mat';
% s_savename = 'Malawi_MCD19A2_Aug_2017.mat';
% s_savename = 'Kinshasa_MCD19A2_Mar2018-Oct2019.mat';
% s_savename = 'Malawi_MCD19A2_Jun2017_Jul2018.mat';
% s_savename = 'Kampala_MCD19A2_Jan_Dec_2019.mat';
% s_savename = 'Addis_Ababa_MCD19A2_Jan_Dec_2019.mat';
% s_savename = 'Dhaka_MCD19A2_Jan_Feb_2018.mat';
% s_savename = 'Dhaka_MCD19A2_Mar_Apr_2018.mat';
% s_savename = 'Dhaka_MCD19A2_May_Jun_2018.mat';
% s_savename = 'Dhaka_MCD19A2_Jul_Aug_2018.mat';
% s_savename = 'Dhaka_MCD19A2_Sep_Oct_2018.mat';
% s_savename = 'Dhaka_MCD19A2_Jul_Dec_2018.mat';
% s_savename = 'Dhaka_MCD19A2_Jan_Jun_2018.mat';
% s_savename = 'Dhaka_MCD19A2_Jan_Apr_2018.mat';
% s_savename = 'Dhaka_MCD19A2_May_Aug_2018.mat';
% s_savename = 'Dhaka_MCD19A2_Sep_Dec_2018.mat';
% s_savename = 'Dhaka_MCD19A2_Jan_Dec_2017.mat';
% s_savename = 'Dhaka_MCD19A2_Jan_Dec_2019.mat';
% s_savename = 'Dhaka_MCD19A2_Jan_Dec_2022.mat';
 s_savename = 'BD_MCD19A2_h25_Jan_Dec_2013_No_res .mat';

%% Load Data

struct_dir = dir(s_folder);

% find valid data files
v_delete = [];

for i_file = 1:length(struct_dir)
    if isempty(strfind(struct_dir(i_file).name,'MCD19A2')) || isempty(strfind(struct_dir(i_file).name,'.hdf'))
        v_delete = [v_delete; i_file];
    end
end
struct_dir(v_delete) = [];

% Determine Relevant File Information - Tiles and timestamps
c_datestamps = cell(length(struct_dir),1);
c_tiles = cell(length(struct_dir),1);
for i_file = 1:length(struct_dir)
    c_split_filename = strsplit(struct_dir(i_file).name,'.');
    c_datestamps{i_file} = c_split_filename{2};
    c_tiles{i_file} = c_split_filename{3};
end

c_datestamps = unique(c_datestamps);
c_tiles = unique(c_tiles);
n_datestamps = length(c_datestamps);
n_tiles = length(c_tiles);
c_datatypes = {'Optical_Depth_047','Optical_Depth_055'};
n_datatypes = length(c_datatypes);

c_raw_data = cell(n_datestamps,n_tiles,n_datatypes);
c_orbit_info = cell(n_datestamps,n_tiles,n_datatypes);
c_range_info = cell(n_datestamps,n_tiles,n_datatypes);
c_scale_info = cell(n_datestamps,n_tiles,n_datatypes);
c_offset_info = cell(n_datestamps,n_tiles,n_datatypes);

% Load Files
for i_file = 1:length(struct_dir)
    % Identify File
    c_split_filename = strsplit(struct_dir(i_file).name,'.');
    s_datestamp = c_split_filename{2};
    s_tile = c_split_filename{3};
    
    i_datestamp = find(strcmp(s_datestamp,c_datestamps));
    i_tile = find(strcmp(s_tile,c_tiles));
    
    % Get file info
    hdf_info = hdfinfo(fullfile(s_folder,struct_dir(i_file).name));
    
    % Extract Orbit Time Information
    id_orbit_timestamp = find(strcmp({hdf_info.Attributes.Name},'Orbit_time_stamp'));
    s_obrits = hdf_info.Attributes(id_orbit_timestamp).Value;
    c_orbits = strsplit(s_obrits,' ');
    v_orbits = nan(length(c_orbits),1);
    for i_orbit = 1:length(c_orbits)
        if ~isempty(c_orbits{i_orbit})
            n_year = str2num(c_orbits{i_orbit}(1:4));
            n_day = str2num(c_orbits{i_orbit}(5:7));
            n_hour = str2num(c_orbits{i_orbit}(8:9));
            n_minute = str2num(c_orbits{i_orbit}(10:11));
            v_orbits(i_orbit) = datenum(n_year,1,1,0,0,0) + (n_day-1) + (n_hour/24) + (n_minute/(24*60));
        end
    end
    v_orbits(isnan(v_orbits)) = [];
    for i_type = 1:n_datatypes
        c_orbit_info{i_datestamp,i_tile,i_type} = v_orbits;
    end
    
    for i_type = 1:n_datatypes
        id_grid = find(strcmp({hdf_info.Vgroup.Name},s_grid_to_load));
        id_data_fields = find(strcmp({hdf_info.Vgroup(id_grid).Vgroup.Name},'Data Fields'));
        
        id_type = find(strcmp({hdf_info.Vgroup(id_grid).Vgroup(id_data_fields).SDS.Name},c_datatypes{i_type}));
        % Extract valid range 
        id_range_info = find(strcmp({hdf_info.Vgroup(id_grid).Vgroup(id_data_fields).SDS(id_type).Attributes.Name},'valid_range'));
        c_range_info{i_datestamp,i_tile,i_type} = hdf_info.Vgroup(id_grid).Vgroup(id_data_fields).SDS(id_type).Attributes(id_range_info).Value;
        
        % Extract Scaling Info
        id_scale_info = find(strcmp({hdf_info.Vgroup(id_grid).Vgroup(id_data_fields).SDS(id_type).Attributes.Name},'scale_factor'));
        c_scale_info{i_datestamp,i_tile,i_type} = hdf_info.Vgroup(id_grid).Vgroup(id_data_fields).SDS(id_type).Attributes(id_scale_info).Value;
        
        % Extract Offset Info
        id_offset_info = find(strcmp({hdf_info.Vgroup(id_grid).Vgroup(id_data_fields).SDS(id_type).Attributes.Name},'add_offset'));
        c_offset_info{i_datestamp,i_tile,i_type} = hdf_info.Vgroup(id_grid).Vgroup(id_data_fields).SDS(id_type).Attributes(id_offset_info).Value;
        
    end
    
    % Extract Data
    for i_type = 1:n_datatypes
        c_raw_data{i_datestamp,i_tile,i_type} = hdfread(fullfile(s_folder,struct_dir(i_file).name),c_datatypes{i_type});
    
    end
    display(['File ',num2str(i_file),' of ',num2str(length(struct_dir)),' Processed.'])
end


%% Interpret Data

% Assign Grid cordinates
m_grid_bound = nan(n_tiles,8);
for i_tile = 1:n_tiles
    i_h = strfind(c_tiles{i_tile},'h');
    i_v = strfind(c_tiles{i_tile},'v');
    n_h = str2num(c_tiles{i_tile}(i_h+1:i_v-1));
    n_v = str2num(c_tiles{i_tile}(i_v+1:length(c_tiles{i_tile})));
    id_line = find((n_h == t_grid_info.ih).*(n_v == t_grid_info.iv));
    m_grid_bound(i_tile,:) = table2array(t_grid_info(id_line,3:10));
end

c_lat = cell(n_tiles);
c_lon = cell(n_tiles);

for i_tile = 1:n_tiles
    Tile_Done = 0;
    for i_date = 1:n_datestamps
        for i_type = 1:n_datatypes
            if ~isempty(c_raw_data{i_date,i_tile,i_type}) && (Tile_Done == 0)
                n_y = size(c_raw_data{i_date,i_tile,i_type},2);
                n_x = size(c_raw_data{i_date,i_tile,i_type},3);
                
                
                
                % Create Meshgrid
               
%                 v_lr = m_grid_bound(i_tile,7:8);
%                 v_ll = m_grid_bound(i_tile,1:2);
%                 v_ur = m_grid_bound(i_tile,5:6);
%                 v_ul = m_grid_bound(i_tile,3:4);
%                 
%                 v_ul = m_grid_bound(i_tile,3:4);
%                 v_ur = m_grid_bound(i_tile,1:2);
%                 v_ll = m_grid_bound(i_tile,5:6);
%                 v_lr = m_grid_bound(i_tile,7:8);
                
                % Corners of v4 h11:
                % 0,0:          49.995833,-108.884749
                % 0,1199:       49.995833,-93.341822
                % 1199,0:       40.004167,-91.378647
                % 1199,1199:    40.004167,-78.334657
                
                % 599,599:      45.004167,-91.936460
                % 599,0:        45.004167,-98.996256
                % 599,1199:     45.004167,-84.864878
                
                % FOR v4 h11 ONLY!
%                 v_ul = [-108.884749,49.995833];
%                 v_ur = [-91.378647,40.004167];
%                 v_ll = [-93.341822,49.995833];
%                 v_lr = [-78.334657,40.004167];
                
%                 v_ul = [-108.884749,49.995833];
%                 v_ll = [-91.378647,40.004167];
%                 v_ur = [-93.341822,49.995833];
%                 v_lr = [-78.334657,40.004167];
                
                % Improved reading from chart - TEST THIS OUT!!!
                v_lr = m_grid_bound(i_tile,7:8);
                v_ll = m_grid_bound(i_tile,1:2);
                v_ur = m_grid_bound(i_tile,5:6);
                v_ul = m_grid_bound(i_tile,3:4);
                delta_y_l = abs(v_ul(2) - v_ll(2))/n_y;
                delta_y_r = abs(v_ur(2) - v_lr(2))/n_y;
                delta_x_u = abs(v_ul(1) - v_ur(1))/n_x;
                delta_x_l = abs(v_ll(1) - v_lr(1))/n_x;
                v_lr = v_lr + [-delta_x_l/2 delta_y_r/2];
                v_ll = v_ll + [delta_x_l/2 delta_y_l/2];
                v_ur = v_ur + [-delta_x_u/2 -delta_y_r/2];
                v_ul = v_ul + [delta_x_u/2 -delta_y_l/2];
                
%                 v_ul = [1   -n_y];
%                 v_ur = [n_x -n_y];
%                 v_ll = [1     -1];
%                 v_lr = [n_x   -1];
                
                mesh_X = nan(1,n_y,n_x);
                mesh_Y = nan(1,n_y,n_x);
                
                for i_y = 1:n_y
                    v_l = [0 0];
                    v_r = [0 0];
                    
                    v_l(2) = v_ul(2) + ((i_y-1)/(n_y-1))*(v_ll(2) - v_ul(2));
                    v_r(2) = v_ur(2) + ((i_y-1)/(n_y-1))*(v_lr(2) - v_ur(2));
                    
                    v_l(1) = v_ul(1)*(cos((pi/180)*v_ul(2))/cos((pi/180)*v_l(2)));
                    v_r(1) = v_ur(1)*(cos((pi/180)*v_ur(2))/cos((pi/180)*v_r(2)));
                    
                    
                    if v_r(1) == v_l(1)
                        mesh_X(1,i_y,:) = reshape(v_r(1)*ones(1,n_x),1,1,n_x);
                    else
                        mesh_X(1,i_y,:) = reshape((v_l(1):((v_r(1)-v_l(1))/(n_x-1)):v_r(1)),1,1,n_x);
                    end
                    if v_r(2) == v_l(2)
                        mesh_Y(1,i_y,:) = reshape(v_r(2)*ones(1,n_x),1,1,n_x);
                    else
                        mesh_Y(1,i_y,:) = reshape((v_l(2):((v_r(2)-v_l(2))/(n_x-1)):v_r(2)),1,1,n_x);
                    end

                    
                    
                end
                
                c_lat{i_tile} = mesh_Y(1,:,:);
                c_lon{i_tile} = mesh_X(1,:,:);
                
                Tile_Done = 1;
                
            end
        end
    end
end

% Combine into one big thing???
a_data = [];
v_time = [];
m_coords = [];
c_names = {};

%c_names = {};

n_max_name = 0;
for i_tile = 1:n_tiles
    for i_date = 1:n_datestamps
        for i_type = 1:n_datatypes
            if ~isempty(c_raw_data{i_date,i_tile,i_type})
                n_orbit = size(c_raw_data{i_date,i_tile,i_type},1);
                n_x = size(c_raw_data{i_date,i_tile,i_type},2);
                n_y = size(c_raw_data{i_date,i_tile,i_type},3);
                if (n_y == 1)
                    n_orbit = 1;
                    n_x = size(c_raw_data{i_date,i_tile,i_type},1);
                    n_y = size(c_raw_data{i_date,i_tile,i_type},2);
                    SINGLEORBIT = 1;
                else
                    SINGLEORBIT = 0;
                end
                
                v_time_temp = c_orbit_info{i_date,i_tile,i_type};
                a_data_temp = nan(n_orbit,n_x*n_y,1);
                v_names_temp = (1:(n_x*n_y)) + n_max_name;
                n_max_name_temp = max(v_names_temp);
                c_names_temp = cellstr(num2str(v_names_temp'))';
                
                a_data_temp_2 = double(c_raw_data{i_date,i_tile,i_type});
                a_data_temp_2(a_data_temp_2 < c_range_info{i_date,i_tile,i_type}(1)) = NaN;
                a_data_temp_2(a_data_temp_2 > c_range_info{i_date,i_tile,i_type}(2)) = NaN;
                a_data_temp_2 = a_data_temp_2*c_scale_info{i_date,i_tile,i_type} + c_offset_info{i_date,i_tile,i_type};
                
                for i_orbit = 1:n_orbit
                    if SINGLEORBIT == 1
                        a_data_temp(i_orbit,:) = reshape(a_data_temp_2,1,n_x*n_y,1);
                    else
                        a_data_temp(i_orbit,:) = reshape(a_data_temp_2(i_orbit,:,:),1,n_x*n_y,1);
                    end
                end
                
                m_coords_temp = [reshape(c_lat{i_tile},n_x*n_y,1) reshape(c_lon{i_tile},n_x*n_y,1)];
                
                % Trimming
                if RESTRICT_SIZE == 1
                    id_indeces_to_keep = find((m_coords_temp(:,1) >= min(v_lat_lims)).*(m_coords_temp(:,1) <= max(v_lat_lims)).*(m_coords_temp(:,2) >= min(v_lon_lims)).*(m_coords_temp(:,2) <= max(v_lon_lims)));
                    a_data_temp = a_data_temp(:,id_indeces_to_keep,:);
                    m_coords_temp = m_coords_temp(id_indeces_to_keep,:);
                    c_names_temp = c_names_temp(id_indeces_to_keep);
                end
                
                
            end

            if isempty(a_data)
                a_data = a_data_temp;
                v_time = v_time_temp;
                m_coords = m_coords_temp;
                c_names = c_names_temp;
                c_types = c_datatypes(i_type);
                c_codes = {};
            else
                [v_time,m_timestamps,a_data,c_types,c_codes,c_names] = ...
                    CAPS_Combine_Arrays(v_time,a_data,c_types,c_codes,c_names,...
                                        v_time_temp,a_data_temp,c_datatypes(i_type),{},c_names_temp);
                                    
                [~,id_coords] = ismember(c_names_temp,c_names);
                
                m_coords(id_coords,:) = m_coords_temp;


            end
            display(['Processed: tile ',num2str(i_tile),'; datestamp ',num2str(i_date),'; type ',num2str(i_type)])
        end
    end
    n_max_name = max(n_max_name,n_max_name_temp);
end



figure
i_type_to_plot = 2;
hold on
scatter(m_coords(:,2),m_coords(:,1),1,nanmean(a_data(:,:,i_type_to_plot))','filled')

%% Restrict size and save
if RESTRICT_SIZE == 0
    a_data_save = a_data;
    m_coords_save = m_coords;
    id_indeces_to_keep = find((m_coords(:,1) >= min(v_lat_lims)).*(m_coords(:,1) <= max(v_lat_lims)).*(m_coords(:,2) >= min(v_lon_lims)).*(m_coords(:,2) <= max(v_lon_lims)));
    a_data = a_data(:,id_indeces_to_keep,:);
    m_coords = m_coords(id_indeces_to_keep,:);
end

%save(s_savename,'v_time','a_data','m_coords','c_types');
save(s_savename,'v_time','a_data','m_coords','c_types','-v7.3');


%%  Re-averaging

n_averaging_time = 1;
a_data_save = a_data;
v_time_save = v_time;

[v_time,m_timestamps,a_data] = CAPS_Data_Average(v_time,m_timestamps,a_data,c_types,n_averaging_time,floor(min(v_time)));

%% Plotting
% 
% for i_timestep_to_plot = 1:length(v_time)
% i_type_to_plot = 2;
% 
% 
% figure
% hold on
% %contourf(m_coords(:,2),m_coords(:,1),a_data(i_timestep_to_plot,:,i_type_to_plot)','LineStyle','none')
% scatter(m_coords(:,2),m_coords(:,1),1,a_data(i_timestep_to_plot,:,i_type_to_plot)','filled')
% 
% 
% %40.521570, -80.155726
% %40.262697, -79.788981
% 
% % -1.89, 30.04
% % -2.01, 30.17
% 
% % plot([-79.788981,-79.788981,-80.155726,-80.155726,-79.788981],[40.262697,40.521570,40.521570,40.262697,40.262697],'k-')
% plot([30.04,30.04,30.17,30.17,30.04],[-2.01,-1.89,-1.89,-2.01,-2.01],'k-')
% colorbar
% title([c_types{i_type_to_plot},' @ ',datestr(v_time(i_timestep_to_plot))])
% ylabel('Lat')
% xlabel('Lon')
% end

%% Satellite Passes
% 
% m_passes = datevec(v_time);
% histogram(m_passes(:,4),'BinEdges',-0.5:1:23.5)

%%Developed%% Zarrah. AUI, Malling. Carl.%%%
