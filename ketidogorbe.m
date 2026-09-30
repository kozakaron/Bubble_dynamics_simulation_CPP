clc; clear; close all;

% Fájlok útvonalai
file1 = 'raw_data/run_4886.csv'; % -> plogext
file2 = 'raw_data/run_4674.csv'; % -> old mechanism

% Adatok beolvasása table-be
opts = detectImportOptions(file1);
data1 = readtable(file1, opts);
data2 = readtable(file2, opts);

% Előkészítés a Top 4 anyag meghatározásához (H2, N2, NH3 nélkül)
exclude_cols = {'t', 'R', 'R_dot', 'T', 'AR', 'H', 'H2', 'N2', 'NH3', ...
                'dissipated_energy', 'p_excitation', 'p_internal'};
all_vars = data1.Properties.VariableNames;
chem_vars = setdiff(all_vars, exclude_cols, 'stable');

final_vals = zeros(length(chem_vars), 1);
for i = 1:length(chem_vars)
    val_vector = data1.(chem_vars{i});
    final_vals(i) = val_vector(end);
end
[~, sorted_idx] = sort(final_vals, 'descend');
top4_vars = chem_vars(sorted_idx(1:min(4, length(chem_vars))));

% =========================================================================
% 1. FIGURE: LOGARITMIKUS IDŐSKÁLA ('XScale', 'log')
% =========================================================================
figure('Name', 'Szimulációs Eredmények - Logaritmikus Időskála', 'Position', [50, 50, 800, 900]);

% --- 1.1. Subplot: Sugár és Hőmérséklet ---
subplot(3, 1, 1);
yyaxis left
plot(data1.t, data1.R * 1e6, 'LineWidth', 1.5, 'DisplayName', 'plogext - R');
hold on;
plot(data2.t, data2.R * 1e6, '--', 'LineWidth', 1.5, 'DisplayName', 'old mechanism - R');
ylabel('Buboréksugár, R [\mu m]');
set(gca, 'YScale', 'log');

yyaxis right
plot(data1.t, data1.T, 'LineWidth', 1.5, 'DisplayName', 'plogext - T');
plot(data2.t, data2.T, '--', 'LineWidth', 1.5, 'DisplayName', 'old mechanism - T');
ylabel('Hőmérséklet, T [K]');
set(gca, 'YScale', 'log');

title('Buboréksugár és Hőmérséklet időgörbéi (Log időskála)');
xlabel('Idő, t [s]');
set(gca, 'XScale', 'log');
legend('Location', 'best');
grid on;
hold off;

% --- 1.2. Subplot: Ammónia (NH3) ---
subplot(3, 1, 2);
plot(data1.t, data1.NH3, 'LineWidth', 1.5, 'DisplayName', 'plogext - NH_3');
hold on;
plot(data2.t, data2.NH3, '--', 'LineWidth', 1.5, 'DisplayName', 'old mechanism - NH_3');
title('Ammónia (NH_3) időgörbéi (Log időskála)');
xlabel('Idő, t [s]');
ylabel('n_i [mol]');
set(gca, 'XScale', 'log');
set(gca, 'YScale', 'log');
ylim([1e-20, inf]);
legend('Location', 'best');
grid on;
hold off;

% --- 1.3. Subplot: Top 4 anyag ---
subplot(3, 1, 3);
hold on;
colors = lines(length(top4_vars));
legend_labels = {};
for i = 1:length(top4_vars)
    var_name = top4_vars{i};
    plot(data1.t, data1.(var_name), 'Color', colors(i,:), 'LineStyle', '-', 'LineWidth', 1.5);
    legend_labels{end+1} = ['plogext - ', var_name];
    
    if ismember(var_name, data2.Properties.VariableNames)
        plot(data2.t, data2.(var_name), 'Color', colors(i,:), 'LineStyle', '--', 'LineWidth', 1.5);
        legend_labels{end+1} = ['old mechanism - ', var_name];
    end
end
title('Egyéb legfőbb 4 anyag időgörbéi (Log időskála)');
xlabel('Idő, t [s]');
ylabel('n_i [mol]');
set(gca, 'XScale', 'log');
set(gca, 'YScale', 'log');
ylim([1e-20, inf]);
legend(legend_labels, 'Location', 'southeast', 'NumColumns', 2); % Jobb alsó sarok
grid on; box on;
hold off;


% =========================================================================
% 2. FIGURE: LINEÁRIS IDŐSKÁLA ('XScale', 'linear')
% =========================================================================
figure('Name', 'Szimulációs Eredmények - Lineáris Időskála', 'Position', [880, 50, 800, 900]);

% --- 2.1. Subplot: Sugár és Hőmérséklet ---
subplot(3, 1, 1);
yyaxis left
plot(data1.t, data1.R * 1e6, 'LineWidth', 1.5, 'DisplayName', 'plogext - R');
hold on;
plot(data2.t, data2.R * 1e6, '--', 'LineWidth', 1.5, 'DisplayName', 'old mechanism - R');
ylabel('Buboréksugár, R [\mu m]');

yyaxis right
plot(data1.t, data1.T, 'LineWidth', 1.5, 'DisplayName', 'plogext - T');
plot(data2.t, data2.T, '--', 'LineWidth', 1.5, 'DisplayName', 'old mechanism - T');
ylabel('Hőmérséklet, T [K]');

title('Buboréksugár és Hőmérséklet időgörbéi (Lineáris időskála)');
xlabel('Idő, t [s]');
set(gca, 'XScale', 'linear');
legend('Location', 'best');
grid on;
hold off;

% --- 2.2. Subplot: Ammónia (NH3) ---
subplot(3, 1, 2);
plot(data1.t, data1.NH3, 'LineWidth', 1.5, 'DisplayName', 'plogext - NH_3');
hold on;
plot(data2.t, data2.NH3, '--', 'LineWidth', 1.5, 'DisplayName', 'old mechanism - NH_3');
title('Ammónia (NH_3) időgörbéi (Lineáris időskála)');
xlabel('Idő, t [s]');
ylabel('n_i [mol]');
set(gca, 'XScale', 'linear');
set(gca, 'YScale', 'log');
ylim([1e-20, inf]);
legend('Location', 'best');
grid on;
hold off;

% --- 2.3. Subplot: Top 4 anyag ---
subplot(3, 1, 3);
hold on;
legend_labels_lin = {};
for i = 1:length(top4_vars)
    var_name = top4_vars{i};
    plot(data1.t, data1.(var_name), 'Color', colors(i,:), 'LineStyle', '-', 'LineWidth', 1.5);
    legend_labels_lin{end+1} = ['plogext - ', var_name];
    
    if ismember(var_name, data2.Properties.VariableNames)
        plot(data2.t, data2.(var_name), 'Color', colors(i,:), 'LineStyle', '--', 'LineWidth', 1.5);
        legend_labels_lin{end+1} = ['old mechanism - ', var_name];
    end
end
title('Egyéb legfőbb 4 anyag időgörbéi (Lineáris időskála)');
xlabel('Idő, t [s]');
ylabel('n_i [mol]');
set(gca, 'XScale', 'linear');
set(gca, 'YScale', 'log');
ylim([1e-20, inf]);
legend(legend_labels_lin, 'Location', 'southeast', 'NumColumns', 2); % Jobb alsó sarok
grid on; box on;
hold off;


% =========================================================================
% ÖSSZEHASONLÍTÓ TÁBLÁZAT A FOLYAMAT VÉGÉN (Command Window kimenet)
% Csak NH3 és a Top 4 anyag (N2 és H2 nélkül)
% =========================================================================
table_species = [{'NH3'}, top4_vars(:)'];

Material = table_species';
OldMechanism_Yield = zeros(length(table_species), 1);
Plogext_Yield     = zeros(length(table_species), 1);
Diff_Percent       = zeros(length(table_species), 1);

for i = 1:length(table_species)
    sp = table_species{i};
    
    % Végső értékek
    val1 = data1.(sp)(end); % plogext
    val2 = data2.(sp)(end); % old mechanism (Referencia)
    
    OldMechanism_Yield(i) = val2;
    Plogext_Yield(i)     = val1;
    
    % Eltérés százalékban az old_mechanism-hez (referencia) képest
    if val2 ~= 0
        Diff_Percent(i) = ((val1 - val2) / val2) * 100;
    else
        Diff_Percent(i) = NaN;
    end
end

comparisonTable = table(Material, OldMechanism_Yield, Plogext_Yield, Diff_Percent, ...
    'VariableNames', {'Anyag', 'old_mechanism_Ref_mol', 'plogext_mol', 'Eltérés_százalék'});

disp(' ');
disp('========================================================');
disp('   VÉGSŐ HOZAM ÖSSZEHASONLÍTÁS (Referencia: old mechanism)');
disp('========================================================');
disp(comparisonTable);