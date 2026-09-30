clc;clear;close all;

% Fájlok útvonalai
file1 = 'raw_data/run_4886.csv'; % -> plogext
file2 = 'raw_data/run_4674.csv'; % -> old mechanism

% Adatok beolvasása table-be
opts = detectImportOptions(file1);
data1 = readtable(file1, opts);
data2 = readtable(file2, opts);

% -------------------------------------------------------------------------
% 1. ÁBRA: Buboréksugár (mikrométerben, bal tengely) és Hőmérséklet (jobb tengely)
% -------------------------------------------------------------------------
figure('Name', 'Szimulációs Eredmények Összehasonlítása (Teljesen Log)', 'Position', [100, 100, 900, 900]);

subplot(3, 1, 1);
% Bal tengely: Sugár átszámítása méterről mikrométerre (* 1e6)
yyaxis left
plot(data1.t, data1.R * 1e6, 'LineWidth', 1.5, 'DisplayName', 'plogext - R');
hold on;
plot(data2.t, data2.R * 1e6, '--', 'LineWidth', 1.5, 'DisplayName', 'old mechanism - R');
ylabel('Buboréksugár, R [\mu m]');
set(gca, 'YScale', 'log'); % Bal Y tengely logaritmussá tétele

% Jobb tengely: Hőmérséklet (T)
yyaxis right
plot(data1.t, data1.T, 'LineWidth', 1.5, 'DisplayName', 'plogext - T');
plot(data2.t, data2.T, '--', 'LineWidth', 1.5, 'DisplayName', 'old mechanism - T');
ylabel('Hőmérséklet, T [K]');
set(gca, 'YScale', 'log'); % Jobb Y tengely logaritmussá tétele

title('Buboréksugár és Hőmérséklet időgörbéi (Log-Log skála)');
xlabel('Idő, t [s]');
set(gca, 'XScale', 'log'); % Időtengely logaritmussá tétele
legend('Location', 'best');
grid on;
hold off;

% -------------------------------------------------------------------------
% 2. ÁBRA: Ammónia (NH3) időgörbék
% -------------------------------------------------------------------------
subplot(3, 1, 2);
plot(data1.t, data1.NH3, 'LineWidth', 1.5, 'DisplayName', 'plogext - NH_3');
hold on;
plot(data2.t, data2.NH3, '--', 'LineWidth', 1.5, 'DisplayName', 'old mechanism - NH_3');
title('Ammónia (NH_3) időgörbéi (Log-Log skála)');
xlabel('Idő, t [s]');
ylabel('n_i (mol)');
set(gca, 'XScale', 'log'); % Időtengely logaritmussá tétele
set(gca, 'YScale', 'log'); % Függőleges tengely logaritmussá tétele
legend('Location', 'best');
grid on;
hold off;

% -------------------------------------------------------------------------
% 3. ÁBRA: Top 4 anyag (kivéve H2, N2, NH3) a folyamat végén
% -------------------------------------------------------------------------
subplot(3, 1, 3);
hold on;

% Nem kívánt oszlopok listája (ezeket kihagyjuk a top 4 keresésből)
exclude_cols = {'t', 'R', 'R_dot', 'T', 'AR', 'H', 'H2', 'N2', 'NH3', ...
                'dissipated_energy', 'p_excitation', 'p_internal'};

% Minden oszlop, ami nem tartozik a fenti fix listához
all_vars = data1.Properties.VariableNames;
chem_vars = setdiff(all_vars, exclude_cols, 'stable');

% Megkeressük a plogext (data1) végén a legnagyobb értékű 4 anyagot
final_vals = zeros(length(chem_vars), 1);
for i = 1:length(chem_vars)
    val_vector = data1.(chem_vars{i});
    final_vals(i) = val_vector(end);
end

[~, sorted_idx] = sort(final_vals, 'descend');
top4_vars = chem_vars(sorted_idx(1:min(4, length(chem_vars))));

colors = lines(length(top4_vars));

legend_labels = {};
for i = 1:length(top4_vars)
    var_name = top4_vars{i};
    
    % plogext rajzolása
    plot(data1.t, data1.(var_name), 'Color', colors(i,:), 'LineStyle', '-', 'LineWidth', 1.5);
    legend_labels{end+1} = ['plogext - ', var_name];
    
    % old mechanism rajzolása (ugyanolyan szín, de szaggatott vonal)
    if ismember(var_name, data2.Properties.VariableNames)
        plot(data2.t, data2.(var_name), 'Color', colors(i,:), 'LineStyle', '--', 'LineWidth', 1.5);
        legend_labels{end+1} = ['old mechanism - ', var_name];
    end
end

title('Egyéb legfőbb 4 anyag időgörbéi (H_2, N_2, NH_3 nélkül, Log-Log skála)');
xlabel('Idő, t [s]');
ylabel('n_i (mol)');
set(gca, 'XScale', 'log'); % Időtengely logaritmussá tétele
set(gca, 'YScale', 'log'); % Függőleges tengely logaritmussá tétele
legend(legend_labels, 'Location', 'best', 'NumColumns', 2);
grid on;
hold off;

% -------------------------------------------------------------------------
% ÖSSZEHASONLÍTÓ TÁBLÁZAT A FOLYAMAT VÉGÉN (Command Window kimenet)
% -------------------------------------------------------------------------
table_species = {'NH3', 'H2', 'N2', top4_vars{1}};

Material = table_species';
Plogext_Yield     = zeros(length(table_species), 1);
OldMechanism_Yield = zeros(length(table_species), 1);
Diff_Percent       = zeros(length(table_species), 1);

for i = 1:length(table_species)
    sp = table_species{i};
    
    % Végső érték
    val1 = data1.(sp)(end); % plogext
    val2 = data2.(sp)(end); % old mechanism
    
    Plogext_Yield(i) = val1;
    OldMechanism_Yield(i) = val2;
    
    % Eltérés százalékban a plogext-hez (referencia) képest
    if val1 ~= 0
        Diff_Percent(i) = ((val2 - val1) / val1) * 100;
    else
        Diff_Percent(i) = NaN;
    end
end

comparisonTable = table(Material, Plogext_Yield, OldMechanism_Yield, Diff_Percent, ...
    'VariableNames', {'Anyag', 'plogext_Ref_mol', 'old_mechanism_mol', 'Eltérés_százalék'});

disp(' ');
disp('========================================================');
disp('   VÉGSŐ HOZAM ÖSSZEHASONLÍTÁS (Referencia: plogext)');
disp('========================================================');
disp(comparisonTable);

% % Fájlok útvonalai (szükség szerint módosíthatók)
% file1 = 'raw_data/run_4886.csv';
% file2 = 'raw_data/run_4674.csv';
% 
% % Adatok beolvasása table-be
% opts = detectImportOptions(file1);
% data1 = readtable(file1, opts);
% data2 = readtable(file2, opts);
% 
% % -------------------------------------------------------------------------
% % 1. ÁBRA: Buboréksugár (bal tengely) és Hőmérséklet (jobb tengely)
% % -------------------------------------------------------------------------
% figure('Name', 'Szimulációs Eredmények Összehasonlítása', 'Position', [100, 100, 900, 900]);
% 
% subplot(3, 1, 1);
% % Bal tengely: Sugár (R)
% yyaxis left
% plot(data1.t, data1.R, 'LineWidth', 1.5, 'DisplayName', 'run\_4886 - R');
% hold on;
% plot(data2.t, data2.R, '--', 'LineWidth', 1.5, 'DisplayName', 'run\_4674 - R');
% ylabel('Buboréksugár, R [m]');
% 
% % Jobb tengely: Hőmérséklet (T)
% yyaxis right
% plot(data1.t, data1.T, 'LineWidth', 1.5, 'DisplayName', 'run\_4886 - T');
% plot(data2.t, data2.T, '--', 'LineWidth', 1.5, 'DisplayName', 'run\_4674 - T');
% ylabel('Hőmérséklet, T [K]');
% 
% title('Buboréksugár és Hőmérséklet időgörbéi');
% xlabel('Idő, t [s]');
% legend('Location', 'best');
% grid on;
% hold off;
% 
% % -------------------------------------------------------------------------
% % 2. ÁBRA: Ammónia (NH3) időgörbék
% % -------------------------------------------------------------------------
% subplot(3, 1, 2);
% plot(data1.t, data1.NH3, 'LineWidth', 1.5, 'DisplayName', 'run\_4886 - NH_3');
% hold on;
% plot(data2.t, data2.NH3, '--', 'LineWidth', 1.5, 'DisplayName', 'run\_4674 - NH_3');
% title('Ammónia (NH_3) időgörbéi');
% xlabel('Idő, t [s]');
% ylabel('Mennyiség / Koncentráció');
% legend('Location', 'best');
% grid on;
% hold off;
% 
% % -------------------------------------------------------------------------
% % 3. ÁBRA: Top 4 anyag (kivéve H2, N2, NH3) a folyamat végén
% % -------------------------------------------------------------------------
% subplot(3, 1, 3);
% hold on;
% 
% % Nem kívánt oszlopok listája (ezeket kihagyjuk a top 4 keresésből)
% exclude_cols = {'t', 'R', 'R_dot', 'T', 'AR', 'H', 'H2', 'N2', 'NH3', ...
%                 'dissipated_energy', 'p_excitation', 'p_internal'};
% 
% % Minden oszlop, ami nem tartozik a fenti fix listához
% all_vars = data1.Properties.VariableNames;
% chem_vars = setdiff(all_vars, exclude_cols, 'stable');
% 
% % Megkeressük a run_4886 végén a legnagyobb értékű 4 anyagot
% final_vals = zeros(length(chem_vars), 1);
% for i = 1:length(chem_vars)
%     val_vector = data1.(chem_vars{i});
%     final_vals(i) = val_vector(end);
% end
% 
% [~, sorted_idx] = sort(final_vals, 'descend');
% top4_vars = chem_vars(sorted_idx(1:min(4, length(chem_vars))));
% 
% % Színek vagy vonaltípusok a megkülönböztetéshez
% line_styles = {'-', '--', '-.', ':'};
% colors = lines(length(top4_vars));
% 
% legend_labels = {};
% for i = 1:length(top4_vars)
%     var_name = top4_vars{i};
% 
%     % run_4886 rajzolása
%     p1 = plot(data1.t, data1.(var_name), 'Color', colors(i,:), 'LineStyle', '-', 'LineWidth', 1.5);
%     legend_labels{end+1} = ['4886 - ', var_name];
% 
%     % run_4674 rajzolása (ugyanolyan szín, de szaggatott vonal)
%     if ismember(var_name, data2.Properties.VariableNames)
%         p2 = plot(data2.t, data2.(var_name), 'Color', colors(i,:), 'LineStyle', '--', 'LineWidth', 1.5);
%         legend_labels{end+1} = ['4674 - ', var_name];
%     end
% end
% 
% title('Egyéb legfőbb 4 anyag időgörbéi (H_2, N_2, NH_3 nélkül)');
% xlabel('Idő, t [s]');
% ylabel('Mennyiség / Koncentráció');
% legend(legend_labels, 'Location', 'best', 'NumColumns', 2);
% grid on;
% hold off;