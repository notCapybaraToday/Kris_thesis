function result_table = compare_resonances(f1, z1, f2, z2, f3, zref, N)
%COMPARE_RESONANCES compares up to N resonance frequencies (phase = 0 crossings)
% between two models and a reference.
%
% INPUTS:
%   f1, z1   - frequency and impedance of model A
%   f2, z2   - frequency and impedance of model B
%   f3, zref - frequency and impedance of the reference (OpenWind)
%   N        - number of resonances to compare (max)
%
% OUTPUT:
%   result_table - table with resonance numbers, frequencies,
%                  absolute and relative errors

    % Phase extraction
    phase1 = angle(z1);
    phase2 = angle(z2);
    phase_ref = angle(zref);

    % Use enhanced zero-crossing finder with interpolation
    f_ref_all = find_crossings(f3, phase_ref);
    f_model1_all = find_crossings(f1, phase1);
    f_model2_all = find_crossings(f2, phase2);

    % Find how many are available
    min_N = min([length(f_ref_all), length(f_model1_all), length(f_model2_all), N]);

    % Warn if not enough crossings
    if min_N < N
        warning('Only %d resonances available instead of requested %d. Using %d.', ...
                 min_N, N, min_N);
    end

    % Truncate to min_N
    f_ref = f_ref_all(1:min_N);
    f_model1 = f_model1_all(1:min_N);
    f_model2 = f_model2_all(1:min_N);

    % Error calculations
    abs_err1 = abs(f_model1 - f_ref);
    rel_err1 = abs_err1 ./ f_ref;

    abs_err2 = abs(f_model2 - f_ref);
    rel_err2 = abs_err2 ./ f_ref;

    % % Create result table
    % result_table = table((1:min_N)', f_ref', f_model1', abs_err1', 100*rel_err1', ...
    %                                  f_model2', abs_err2', 100*rel_err2', ...
    %     'VariableNames', {'Resonance', 'f_openwind_Hz', 'f_article2_Hz', 'AbsErr_2_Hz', 'RelErr_2_pct', ...
    %                                      'f_book1_Hz', 'AbsErr_1_Hz', 'RelErr_1_pct'});
    % 
    % % Display result
    % disp(result_table);

    % Create result table
result_table = table((1:min_N)', f_ref', f_model1', abs_err1', 100*rel_err1', ...
                                 f_model2', abs_err2', 100*rel_err2', ...
    'VariableNames', {'Resonance', 'f_ref_Hz', 'f_a_Hz', 'AbsErr_a_Hz', 'RelErr_a_pct', ...
                                     'f_b_Hz', 'AbsErr_b_Hz', 'RelErr_b_pct'});

% Display result table in MATLAB
disp(result_table);

% --- Generate LaTeX Table ---
fprintf('\n%% ------------------- LaTeX Table -------------------\n');
fprintf('\\begin{tabular}{c|c|c|c|c|c|c|c}\n');
fprintf('Res & $f_{\\text{ref}}$ & $f_a$ & $|\\Delta f_a|$ & RelErr$_a$ & $f_b$ & $|\\Delta f_b|$ & RelErr$_b$ \\\\\n');
fprintf('    & [Hz] & [Hz] & [Hz] & [\\%%] & [Hz] & [Hz] & [\\%%] \\\\\n');
fprintf('\\hline\n');

for i = 1:min_N
    fprintf('%d & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f \\\\\n', ...
        i, f_ref(i), f_model1(i), abs_err1(i), 100*rel_err1(i), ...
           f_model2(i), abs_err2(i), 100*rel_err2(i));
end

fprintf('\\end{tabular}\n');
fprintf('%% ---------------------------------------------------\n');
end

function f_zero = find_crossings(f, phase)
% Finds zero crossings (any sign change) with linear interpolation

    f_zero = [];
    for i = 1:length(phase)-1
        if sign(phase(i)) ~= sign(phase(i+1)) && ~isnan(phase(i)) && ~isnan(phase(i+1))
            % Linear interpolation to estimate zero crossing
            f_interp = interp1(phase(i:i+1), f(i:i+1), 0);
            f_zero(end+1) = f_interp; %#ok<AGROW>
        end
    end
end
