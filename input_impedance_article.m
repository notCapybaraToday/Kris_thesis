function [Z_in_lossy, Z_in_lossless, f] = input_impedance_article(f, L, R1, R2, t, mu)
% Input impedance calculation for a cone.
%  Formulas taken from: TOURNEMENNE Romain and CHABASSIER Julie. “A Comparison of a One‐Dimensional Finite
% Element Method and the Transfer Matrix Method for the Computation of Wind Music Instru‐
% ment Impedance”. In: Acta Acustica united with Acustica 105.6 (2019), pp. 1119–1129.
    % INPUTS:
    % f  - frequency vector (Hz)
    % L  - cone length (m)
    % R1 - small radius of the cone (m)
    % R2 - large radius of the cone (m)
    % t  - temperature in Celsius
    % mu - dynamic viscosity of air (Pa·s)
    %
    % OUTPUTS:
    % Z_in_lossy     - input impedance with losses
    % Z_in_lossless  - input impedance without losses
    % f              - frequency vector (Hz)

    % Derived constants
    gamma = 1.402;
    c_p = 240;
    kappa = 5.77e-3 * (1 + 0.0033 * t); % Thermal conductivity
    rho = 1.292 * 273.15 / (273.15 + t); % Air density (kg/m^3)
    c = 331.45 * sqrt((273.15 + t) / 273.15); % Speed of sound (m/s)
    nu = mu / rho;
    r_dot = 2 * R1 + R2 / 3;
    beta = (R2 - R1) / L / R1;
    
    % Cross-sectional areas
    S1 = pi * R1^2;
    S2 = pi * R2^2;
    
    omega = 2 * pi * f;

    % Complex wavenumbers for losses
    k_v = conj(sqrt(1i * omega * rho / mu));
    k_t = conj(sqrt(1i * omega * rho * c_p / kappa));

    % Lossless wave number
    k = omega / c;
    
    % Radiation impedance
    Z_R = (rho * c / S2) .* ((1/4) * (k .* R2).^2 + 0.6133 * 1i .* k * R2);

    % Zc (characteristic impedance) for lossless and lossy
    Zc = rho * c / S1;

    % Bessel function loss terms
    k_v_arg = k_v * r_dot;
    k_t_arg = k_t * r_dot;
    J_z_v = 2 * besselj(1, k_v_arg) ./ (besselj(0, k_v_arg) .* k_v_arg);
    J_z_t = 2 * besselj(1, k_t_arg) ./ (besselj(0, k_t_arg) .* k_t_arg);

    J_z_v0 = 0;
    J_z_t0 = 0;

    % Impedance including losses
    Z_crw = rho * c .* sqrt(1 ./ ((1 + (gamma - 1) * J_z_t) .* (1 - J_z_v))) / S1;
    Z_crw0 = rho * c * sqrt(1 / ((1 + (gamma - 1) * J_z_t0) * (1 - J_z_v0))) / S1;

    % Propagation constants
    GAMMA = 1i * omega .* sqrt((1 + (gamma - 1) * J_z_t) ./ (1 - J_z_v)) / c;
    GAMMA0 = 1i * omega * sqrt((1 + (gamma - 1) * J_z_t0) / (1 - J_z_v0)) / c;

    % Preallocate output vectors
    Z_in_lossy = zeros(size(f));
    Z_in_lossless = zeros(size(f));

    for idx = 1:length(f)
        % Lossy transfer matrix
        A = R2 * cosh(GAMMA(idx) * L) / R1 - beta * sinh(GAMMA(idx) * L) / GAMMA(idx);
        B = R1 * Z_crw(idx) * sinh(GAMMA(idx) * L) / R2;
        C = (1 / Z_crw(idx)) * ((R2 / R1 - (beta / GAMMA(idx))^2) * sinh(GAMMA(idx) * L) + beta^2 * L * cosh(GAMMA(idx) * L) / GAMMA(idx));
        D = R1 * (cosh(GAMMA(idx) * L) + beta * sinh(GAMMA(idx) * L) / GAMMA(idx)) / R2;
        Z_in_lossy(idx) = ((A * Z_R(idx) + B) / (C * Z_R(idx) + D))/Zc;

        % Lossless transfer matrix
        A0 = R2 * cosh(GAMMA0(idx) * L) / R1 - beta * sinh(GAMMA0(idx) * L) / GAMMA0(idx);
        B0 = R1 * Z_crw0 * sinh(GAMMA0(idx) * L) / R2;
        C0 = (1 / Z_crw0) * ((R2 / R1 - (beta / GAMMA0(idx))^2) * sinh(GAMMA0(idx) * L) + beta^2 * L * cosh(GAMMA0(idx) * L) / GAMMA0(idx));
        D0 = R1 * (cosh(GAMMA0(idx) * L) + beta * sinh(GAMMA0(idx) * L) / GAMMA0(idx)) / R2;
        Z_in_lossless(idx) = ((A0 * Z_R(idx) + B0) / (C0 * Z_R(idx) + D0))/Zc;
    end
end
