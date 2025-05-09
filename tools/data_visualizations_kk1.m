% Declare data
k_vals = [221, 331, 441];
basis_sets = ["cc-pvdz", "cc-pvtz", "cc-pvqz"];
basis_sets_labels = ["DZ", "TZ", "QZ"];
data_types = "CORR";
simulation_types = "deltas";

% Raw data: rows = basis sets (DZ, TZ, QZ), columns = k-points (441, 331, 221)
data_deltas_CORR = fliplr([
    -9.246836463, -13.97767559, -15.30482266;
    -12.54932537, -17.83146619, -19.3211949;
    -13.91597517, -19.42325893, -20.92500701
]);

% x-axis: inverse square of k-point size
x_vals_data = fliplr([2^-2, 3^-2, 4^-2]);  % corresponds to 441, 331, 221
x_vals_labels = fliplr(["2^{-2}", "3^{-2}", "4^{-2}"]);
x_vals = [0, x_vals_data];  % for extrapolation
x_vals_labels = ["0", x_vals_labels];

% Pseudobasis deltas
DZ_TZ = fliplr([-13.9398, -19.4541, -21.0123]);
TZ_QZ = fliplr([-14.9133, -20.5848, -22.0954]);

% Plot colors
colors = [0 0.4470 0.7410;   % DZ
          0.8500 0.3250 0.0980; % TZ
          0.4660 0.7 0.1880];  % QZ
pseudo_colors = {[0.4940 0.1840 0.5560], 'k'};

pseudo_labels = {"(DZ, TZ)", "(TZ, QZ)"};
pseudo_data = {DZ_TZ, TZ_QZ};

for simulation_type = simulation_types
    for data_type = data_types
        % Create figure
        figure("Name", simulation_type + " Structure - "+ data_type + " Energy", ...
               "Color", "w"); 
        hold on; grid on;

        % Load matrix for current data
        iter = eval(strjoin(["data" simulation_type data_type], "_"));

        %% === Plot Original Basis Set Data and Fits ===
        for i = 1:length(basis_sets)
            y_vals = iter(i, :);
            x_data = x_vals(2:end); % exclude extrapolation point

            % Fit: all points & excl. 221
            p_all = polyfit(x_data, y_vals, 1);
            p_excl = polyfit(x_data(1:end-1), y_vals(1:end-1), 1);

            % Evaluate fits
            y_all = polyval(p_all, x_vals);
            y_excl = polyval(p_excl, x_vals);

            % Plot data
            plot(x_data, y_vals, '-o', ...
                'LineWidth', 1.5, ...
                'Color', colors(i, :), ...
                'DisplayName', basis_sets_labels(i) + " data");

            % Plot fit (all)
            plot(x_vals, y_all, '--', ...
                'LineWidth', 1.2, ...
                'Color', colors(i, :), ...
                'DisplayName', basis_sets_labels(i) + " fit (all)");

            % Extrapolated point (all)
            scatter(0, y_all(1), 40, 'o', ...
                'MarkerEdgeColor', colors(i, :), ...
                'DisplayName', basis_sets_labels(i) + " extrap. (all)");

            % Fit excluding 221
            plot(x_vals(1:end-1), y_excl(1:end-1), ':', ...
                'LineWidth', 1.2, ...
                'Color', colors(i, :), ...
                'DisplayName', basis_sets_labels(i) + " fit excl. 221");

            scatter(0, y_excl(1), 60, '*', ...
                'MarkerEdgeColor', colors(i, :), ...
                'DisplayName', basis_sets_labels(i) + " extrap. excl. 221");
        end

        %% === Plot Pseudo Basis Set Differences ===
        for j = 1:2
            y_vals = pseudo_data{j};
            col = pseudo_colors{j};
            label = pseudo_labels{j};

            % Fit (all) and (excl. 221)
            p_all = polyfit(x_data, y_vals, 1);
            p_excl = polyfit(x_data(1:end-1), y_vals(1:end-1), 1);

            y_all = polyval(p_all, x_vals);
            y_excl = polyval(p_excl, x_vals(1:end-1));

            % Plot original
            plot(x_data, y_vals, '-o', ...
                'LineWidth', 1.5, ...
                'Color', col, ...
                'DisplayName', label + " data");

            % Fit line (all)
            plot(x_vals, y_all, '--', ...
                'LineWidth', 1.2, ...
                'Color', col, ...
                'DisplayName', label + " fit");

            scatter(0, y_all(1), 40, 'o', ...
                'MarkerEdgeColor', col, ...
                'DisplayName', label + " extrap.");

            % Fit line (excl. 221)
            plot(x_vals(1:end-1), y_excl, ':', ...
                'LineWidth', 1.2, ...
                'Color', col, ...
                'DisplayName', label + " fit excl. 221");

            scatter(0, y_excl(1), 60, '*', ...
                'MarkerEdgeColor', col, ...
                'DisplayName', label + " extrap. excl. 221");
        end

        %% === Labels and Export ===
        xlabel('K-points mesh $k \times k \times 1 \rightarrow k^{-2}$', 'Interpreter', 'latex');
        xticks(x_vals);
        xticklabels(x_vals_labels);
        ylabel(data_type + " Energy (kJ/mol)", 'Interpreter', 'none');
        title(data_type + " Energy vs. K-points for Different Basis Sets", 'Interpreter', 'none');
        legend('Location', 'bestoutside', 'Interpreter', 'none');
        hold off;

        % Save figure
        filename = strjoin(["data" simulation_type data_type], "_") + ".png";
        exportgraphics(gcf, filename, 'Resolution', 300);
    end
end

% % Declare data
% k_vals = [221, 331, 441];
% basis_sets = ["cc-pvdz", "cc-pvtz", "cc-pvqz"];
% basis_sets_labels = ["DZ", "TZ", "QZ"];
% data_types = "CORR";
% simulation_types = "deltas";
% 
% % Raw data: rows = basis sets (DZ, TZ, QZ), columns = k-points (441, 331, 221)
% data_deltas_CORR = fliplr([
%     -9.246836463, -13.97767559, -15.30482266;
%     -12.54932537, -17.83146619, -19.3211949;
%     -13.91597517, -19.42325893, -20.92500701
% ]);
% 
% % x-axis: inverse square of k-point size
% x_vals_data = fliplr([2^-2, 3^-2, 4^-2]);  % corresponds to 441, 331, 221
% x_vals_labels = fliplr(["2^{-2}", "3^{-2}", "4^{-2}"]);
% 
% 
% % Add x = 0 for extrapolation
% x_vals = [0, x_vals_data];
% x_vals_labels = ["0", x_vals_labels];
% 
% DZ_TZ = fliplr([-13.9398 -19.4541 -21.0123]);
% TZ_QZ = fliplr([-14.9133 -20.5848 -22.0954]);
% 
% % Plot colors for each basis set
% color = [0 0.4470 0.7410;
%          0.8500 0.3250 0.0980;
%          [0.4660 0.7 0.1880]];
% 
% for simulation_type = simulation_types
%     for data_type = data_types
%         figure("Name", simulation_type + " Structure - "+ data_type + " Energy", ...
%                "Color", "w"); % White background
%         hold on; grid on;
% 
%         iter = eval(strjoin(["data" simulation_type data_type], "_"));
% 
%         for i = 1:length(basis_sets)
%             y_vals = iter(i, :);
%             x_data = x_vals(2:end); % exclude x=0 for fitting
% 
%             % Linear fit using polyfit (1st degree)
%             p_all = polyfit(x_data, y_vals, 1);  % using all 3 points
%             p_exclude = polyfit(x_data(1:end-1), y_vals(1:end-1), 1); % exclude 221
% 
%             % Evaluate over full x_vals
%             y_trend_all = polyval(p_all, x_vals);
%             y_trend_excl = polyval(p_exclude, x_vals);
% 
%             % Plot original data
%             plot(x_data, y_vals, '-o', ...
%                 'LineWidth', 1.5, ...
%                 'DisplayName', basis_sets_labels(i) + " data", ...
%                 'Color', color(i, :));
% 
%             % Plot full fit line
%             plot(x_vals, y_trend_all, '--', ...
%                 'LineWidth', 1.2, ...
%                 'DisplayName', basis_sets_labels(i) + " fit (all)", ...
%                 'Color', color(i, :));
% 
%             % Plot extrapolated value at x = 0
%             scatter(0, y_trend_all(1), 40, 'o', ...
%                 'MarkerEdgeColor', color(i, :), ...
%                 'DisplayName', basis_sets_labels(i) + " extrap. (all)");
% 
%             % Plot excluded-221 fit line
%             plot(x_vals(1:end-1), y_trend_excl(1:end-1), ':', ...
%                 'LineWidth', 1.2, ...
%                 'DisplayName', basis_sets_labels(i) + " fit (excl. 221)", ...
%                 'Color', color(i, :));
% 
%             % Plot extrapolated value at x = 0 (excl. 221)
%             scatter(0, y_trend_excl(1), 60, '*', ...
%                 'MarkerEdgeColor', color(i, :), ...
%                 'DisplayName', basis_sets_labels(i) + " extrap. (excl. 221)");
%         end
% 
%             % colors = ['green', 'black'];
%             % pesudo_basis = ["DZ_TZ", "TZ_QZ"];
%             % for j = 1:2
%             % Linear fit using polyfit (1st degree)
%             p_DZ_TZ = polyfit(x_data, DZ_TZ, 1);  % using all 3 points
%             p_TZ_QZ = polyfit(x_data, TZ_QZ, 1);  % using all 3 points
% 
%             y_trend_DZ_TZ  = polyval(p_DZ_TZ, x_vals);
%             y_trend_TZ_QZ  = polyval(p_TZ_QZ, x_vals);
% 
%             p_DZ_TZ_ex_221 = polyfit(x_data(1:end-1),DZ_TZ(1:end-1),1);
%             p_TZ_QZ_ex_221 = polyfit(x_data(1:end-1),TZ_QZ(1:end-1),1);
% 
%             y_trend_DZ_TZ_ex_221  = polyval(p_DZ_TZ_ex_221, x_vals(1:end-1));
%             y_trend_TZ_QZ_ex_221  = polyval(p_TZ_QZ_ex_221, x_vals(1:end-1));
% 
%             % Plot original data
%             plot(x_data, DZ_TZ, '-o', ...
%                 'LineWidth', 1.5, ...
%                 'DisplayName', "(DZ,TZ) data", ...
%                 'Color', [0.4940 0.1840 0.5560]);
% 
%             % Plot excluded-221 fit line
%             plot(x_vals, y_trend_DZ_TZ, '--', ...
%                 'LineWidth', 1.2, ...
%                 'DisplayName',"(DZ, TZ) fit", ...
%                 'Color', [0.4940 0.1840 0.5560]);
% 
%             % Plot extrapolated value at x = 0 (excl. 221)
%             scatter(0, y_trend_DZ_TZ(1), 40, 'o', ...
%                 'MarkerEdgeColor', [0.4940 0.1840 0.5560], ...
%                 'DisplayName',"(DZ, TZ) extrap.");
% 
%             % Plot excluded-221 fit line
%             plot(x_vals(1:end-1), y_trend_DZ_TZ_ex_221, ':', ...
%                 'LineWidth', 1.2, ...
%                 'DisplayName',"(DZ, TZ) fit excl. 221", ...
%                 'Color', [0.4940 0.1840 0.5560]);
% 
%             % Plot extrapolated value at x = 0 (excl. 221)
%             scatter(0, y_trend_DZ_TZ_ex_221(1), 60, '*', ...
%                 'MarkerEdgeColor', [0.4940 0.1840 0.5560], ...
%                 'DisplayName',"(DZ, TZ) extrap. excl. 221");
% 
%             % Plot original data
%             plot(x_data, TZ_QZ, '-o', ...
%                 'LineWidth', 1.5, ...
%                 'DisplayName',  "(TZ,QZ)", ...
%                 'Color', 'black');
% 
%             % Plot excluded-221 fit line
%             plot(x_vals, y_trend_TZ_QZ, '--', ...
%                 'LineWidth', 1.2, ...
%                 'DisplayName', "(TZ, QZ) fit", ...
%                 'Color','black');
% 
%             % Plot extrapolated value at x = 0 (excl. 221)
%             scatter(0, y_trend_TZ_QZ(1), 40, 'o', ...
%                 'MarkerEdgeColor', 'black', ...
%                 'DisplayName',"(TZ, QZ) extrap.");
% 
%             % Plot excluded-221 fit line
%             plot(x_vals(1:end-1), y_trend_TZ_QZ_ex_221, ':', ...
%                 'LineWidth', 1.2, ...
%                 'DisplayName', "(TZ, QZ) fit", ...
%                 'Color','black');
% 
%             % Plot extrapolated value at x = 0 (excl. 221)
%             scatter(0, y_trend_TZ_QZ_ex_221(1), 60, '*', ...
%                 'MarkerEdgeColor', 'black', ...
%                 'DisplayName', "(TZ, QZ) extrap. excl. 221");
%             % end
% 
% 
% 
%         xlabel('K-points mesh (e.g. $k \times k \times 1 \rightarrow k^{-2}$)', 'Interpreter', 'latex');
%         xticks(x_vals);
%         xticklabels(x_vals_labels);
%         ylabel(data_type + " Energy (kJ/mol)");
%         title(data_type + " Energy vs K-points for Different Basis Sets");
%         legend('Location', 'bestoutside', 'Interpreter', 'none');
%         hold off;
% 
%         % Save figure
%         filename = strjoin(["data" simulation_type data_type], "_") + ".png";
%         exportgraphics(gcf, filename, 'Resolution', 300);  % High-resolution export
%     end
% end
