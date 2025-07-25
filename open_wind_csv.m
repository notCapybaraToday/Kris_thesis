function [f, Z] = open_wind_csv(file)
    % Function to read acoustical impedance data from a CSV file
    % Input:
    %   file - file name (string), e.g., 'Impedance_25degC_05Hz.csv'
    % Output:
    %   f - frequency array [Hz]
    %   Z - magnitude of normalized impedance |Z/Zc|

    % Open the file
    fid = fopen(file, 'r');
    if fid == -1
        error('Could not open file.');
    end

    % Skip the first two lines
    fgetl(fid); % Skip #Note
    fgetl(fid); % Skip header line

    % Read the remaining data
    data = textscan(fid, '%f %f %f');
    fclose(fid);

    % Extract columns
    f = data{1};                        % Frequency
    Re_Z = data{2};                    % Real part of Z/Zc
    Im_Z = data{3};                    % Imaginary part of Z/Zc
    Z = Re_Z + 1i * Im_Z; 
end
