function [Z_in_lossy, Z_in_lossless, f] = compute_input_impedance(f, L, R1, R2, t, mu)
% COMPUTE_INPUT_IMPEDANCE computes the input impedance of a conical tube section
% using the Transfer Matrix Method (TMM), for both lossless and lossy cases.
% Formulas taken from the CHAIGNE Antoine and KERGOMARD Jean. Acoustics of Musical Instruments. 1st. Modern
% Acoustics and Signal Processing. Foreword by Murray Campbell. Cham, Switzerland: Springer,
% 2021, p. 852. ISBN: 978‐3‐030‐74334‐6.

% INPUTS:
%   f   - Frequency vector [Hz]
%   L   - Length of the conical section [m]
%   R1  - Radius at the narrow end of the cone [m]
%   R2  - Radius at the wide end of the cone [m]
%   t   - Temperature in Celsius [°C]
%   mu  - Dynamic viscosity of air [Pa·s]
%
% OUTPUTS:
%   Z_in_lossy     - Complex input impedance (Z/Zc) with losses (dimensionless)
%   Z_in_lossless  - Complex input impedance (Z/Zc) without losses (dimensionless)
%   f              - Frequency vector [Hz] (same as input)

    % Geometry and physical parameters
    X1 = R1 * L / (R2 - R1);
    X2 = X1 + L;
    S = pi * R1^2;

    % Air properties based on temperature t [Celsius]
    rho = 1.292 * 273.15 / (273.15 + t);  % Density [kg/m^3]
    c = 331.45 * sqrt((273.15 + t) / 273.15);  % Speed of sound [m/s]
    nu = mu / rho; % Kinematic viscosity [not directly used]

    omega = 2 * pi * f;
    k_ll = omega / c;  % Lossless wavenumber
    Zc = rho * c / S;

    % Lossy wavenumber approximation
    k = (((1i * omega )/ c) + (((1+1i) *3e-5/((R2+R1)/2).* sqrt(f)) ))/1i;

    % Radiation impedance
    Z_R = Zc .* ( (1/4) * (k * R1).^2 + 0.6133 * 1j .* k * R1 );

    % Initialize output arrays
    Z_in_lossless = zeros(size(f));
    Z_in_lossy = zeros(size(f));

    % Transfer matrix loop over all frequencies
    for idx = 1:length(f)
        % LOSSLESS Transfer Matrix
        M = [((R2/R1)*cos(k_ll(idx) * L) - sin(k_ll(idx)*L)/k_ll(idx)/X1), ...
             (1j * rho * c * sin(k_ll(idx) * L) / pi / R1 / R2); 
             (pi * R1 * R2 / rho / c) * (1i * sin(k_ll(idx) * L) * (1 + 1 / k_ll(idx)^2 / X1 / X2) + ...
              cos(k_ll(idx) * L) * (1 / X1 - 1 / X2) / 1i / k_ll(idx)), ...
             R1 * cos(k_ll(idx) * L) / R2 + sin(k_ll(idx) * L) / k_ll(idx) / X2];
        A = M(1,1); B = M(1,2); C = M(2,1); D = M(2,2);
        Z_in_lossless(idx) = ((A * Z_R(idx) + B) / (C * Z_R(idx) + D))./Zc;

        % LOSSY Transfer Matrix
        M_lossy = [((R2/R1)*cos(k(idx) * L) - sin(k(idx)*L)/k(idx)/X1), ...
                   (1j * rho * c * sin(k(idx) * L) / pi / R1 / R2); 
                   (pi * R1 * R2 / rho / c) * (1i * sin(k(idx) * L) * (1 + 1 / k(idx)^2 / X1 / X2) + ...
                    cos(k(idx) * L) * (1 / X1 - 1 / X2) / 1i / k(idx)), ...
                   R1 * cos(k(idx) * L) / R2 + sin(k(idx) * L) / k(idx) / X2];
        A = M_lossy(1,1); B = M_lossy(1,2); C = M_lossy(2,1); D = M_lossy(2,2);
        Z_in_lossy(idx) = ((A * Z_R(idx) + B) / (C * Z_R(idx) + D))./Zc;
    end
end
