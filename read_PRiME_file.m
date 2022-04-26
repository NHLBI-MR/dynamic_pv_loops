function data = read_PRiME_file(filename)

% import PRiME bin dump in to matlab
% filename <.bin>
%
% The binary format is saved using a default Labview mode.
% There are blocks of data, each block starts with a 32-bit integer which is the size of the block. Block size will vary.
% Then, there are 17 channels
%     - 6 pre-adaptive filter ECG,
%     - 6 post-adaptive filter ECG,
%     - X Gradient, Y Gradient, Z Gradient,
%     - IBP0, IBP1
%
% The ECG order is RA, LA, LL, V1, V2, V3. You need to combine the electrodes together to get LeadI, LeadII, etc.
% After the initial 32-bit integer describing the block size, the data is big-endian format, 16-bit signed integer.
% All data is sampled at 1kHz.
% J Kakareka, NHLBI, 2021-01-07

if nargin < 1 % Select data
    [FileName,PathName] = uigetfile('*.*','Select bin FILE');
    filename = [PathName FileName];
    clear PathName FileName;
% else
%     make_nhlbi_toolbox;
%     filename = nhlbi_toolbox.run_path_on_sys(filename);
end


disp('===================='); disp(' ');
disp('Loading raw PRiME data: ' );
disp(['<strong>' filename '</strong>']);
disp(' ');
disp('===================='); disp(' ');

% Load data 

fid = fopen(filename, 'rb');

keep_reading = true;
data = [];

while keep_reading
    
    blk_size = fread(fid, [1, 2], 'int32', 'b');

    if blk_size
        tmp = fread(fid, [blk_size(2), blk_size(1)], 'int16', 'b');
        data = [data; tmp];
    else
        keep_reading = false;
    end
end

fclose(fid);

% convert pressure to mmHg 
data = convert_pressure(data);
data = lowpass(data,499,1000); %lowpass filter for smoothing signals

t_pressure=linspace(0, 1/1000*size(data,1),size(data,1))'; %Time vector using sampling frequency 1kHz
data=[data, t_pressure]; %add one more channel with time vector
end

%% convert pressure to real mmHg units
function data = convert_pressure(data)

    % hard coded values
    % PRIME > mmHg

    lv_scale0    =  0.020166;  
    lv_offset0   = -142.3672;  
    lv_scale1    = 0.020213; 
    lv_offset1   = 51.2406;  
    
    % PRIME offsets
    ibp0_scale  =   0.77;      
    ibp0_offset =  -1795;   
    
    ibp1_scale  =  0.76;
    ibp1_offset = -1935;
    
    % PRESSURE channels
    ibp0_ind = 16;
    ibp1_ind = 17;

% correction
ibp0_data = (((data(:,ibp0_ind) - ibp0_offset) * ibp0_scale) - lv_offset0) * lv_scale0;
ibp1_data = (((data(:,ibp1_ind) - ibp1_offset) * ibp1_scale) - lv_offset1) * lv_scale1;

data(:,[ibp0_ind ibp1_ind]) = [ibp0_data ibp1_data];


end
