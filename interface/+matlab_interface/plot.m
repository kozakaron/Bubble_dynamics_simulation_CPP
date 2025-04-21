function plot_simulation(data, n, show_cpar)
    % PLOT_SIMULATION Plots the results of the simulation from data.
    % Arguments:
    %   data: struct containing the simulation data.
    %   n: how long the plotted time interval should be compared to the collapse time (default: 5).
    %   show_legend: if true, the legend will be visible with every single species (default: false).
    %   show_cpar: if true, the control parameters will be printed on the plot (default: false).

    if nargin < 2, n = 5.0; end
    if nargin < 4, show_cpar = true; end

    % Extract data
    cpar = data.cpar;
    sol = data.sol;
    x = sol.x;
    t = sol.t;

    % Find collapse time
    [~, loc_min] = findpeaks(-x(:, 1)); % Find local minima of R
    if ~isempty(loc_min)
        collapse_time = t(loc_min(1)); % First local minimum
    else
        collapse_time = t(end);
    end

    % Determine time interval for plotting
    t_last = n * collapse_time;
    if t_last < 1e-7 || t(end) < t_last || n < 0 || ~sol.success
        end_index = length(t);
    else
        end_index = find(t >= t_last, 1);
    end

    % Adjust time scale
    if t(end_index) < 1e-3
        t = t(1:end_index) * 1e6; % Convert to microseconds
        t_label = 't [\mus]';
    else
        t = t(1:end_index) * 1e3; % Convert to milliseconds
        t_label = 't [ms]';
    end

    % Extract variables
    R = x(1:end_index, 1); % Radius [m]
    T = x(1:end_index, 3); % Temperature [K]
    c = x(1:end_index, 4:end-1); % Concentrations [mol/cm^3]

    % Compute molar amounts
    V = 4.0 / 3.0 * (100.0 * R).^3 * pi; % Volume [cm^3]
    n_mol = c .* V; % Molar amounts [mol]

    % Plot R and T
    figure;
    ax1 = subplot(1, 1, 1);
    yyaxis left;
    plot(t, R / cpar.R_E, 'b', 'LineWidth', 1.5);
    ylabel('R/R_E [-]', 'Color', 'b');
    yyaxis right;
    plot(t, T, 'r-.', 'LineWidth', 1.5);
    ylabel('T [K]', 'Color', 'r');
    xlabel(t_label);
    grid on;

    % Add control parameters as text box
    if show_cpar
        text_str = sprintf('Initial conditions:\n');
        text_str = [text_str, sprintf('  R_E = %.2f [\x03bcm]\n', 1e6 * cpar.R_E)];
        text_str = [text_str, sprintf('  P_{amb} = %.2f [bar]\n', 1e-5 * cpar.P_amb)];
        text_str = [text_str, sprintf('  T_{inf} = %.2f [°C]\n', cpar.T_inf - 273.15)];
        text_str = [text_str, sprintf('  P_v = %.1f [Pa]\n', cpar.P_v)];
        text_str = [text_str, sprintf('Initial content:\n  ')];
        for i = 1:length(cpar.species)
            text_str = [text_str, sprintf('%d%% %s, ', round(100 * cpar.fractions(i)), cpar.species{i})];
        end
        text_str = [text_str(1:end-2), sprintf('\nExcitation = %s:\n', replace(cpar.excitation_type, "_", " "))];
        for i = 1:length(cpar.excitation_params)
            text_str = [text_str, sprintf('  %s = %.2f [%s]\n', ...
                data.excitation.names{i}, cpar.excitation_params(i), data.excitation.units{i})];
        end
        annotation('textbox', [0.55, 0.6, 0.3, 0.3], 'String', text_str, ...
            'FitBoxToText', 'on', 'BackgroundColor', 'white');
    end

    % Plot molar concentrations
    figure;
    ax2 = subplot(1, 1, 1);
    hold on;
    colors = lines(size(n_mol, 2)); % Generate distinct colors
    for i = 1:size(n_mol, 2)
        plot(t, n_mol(:, i), 'LineWidth', 1.5, 'Color', colors(i, :));
    end
    set(gca, 'YScale', 'log');
    ylim([1e-24, max(max(n_mol)) * 5]); % Limit lower bounds of y-axis
    xlabel(t_label);
    ylabel('n_k [mol]');
    grid on;

    % Add legend
    legend(data.mechanism.species_names, 'Location', 'best');
    hold off;

    print_data(data);
end