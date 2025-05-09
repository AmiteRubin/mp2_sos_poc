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
]');

% data_deltas_CORR = data_deltas_CORR';

% x-axis: inverse cube of k-point size
x_vals_data = fliplr([2^-3, 3^-3, 4^-3]);  % corresponds to 441, 331, 221
x_vals_labels = fliplr(["2^{-3}", "3^{-3}", "4^{-3}"]);

% Add x = 0 for extrapolation
x_vals = [0, x_vals_data];
x_vals_labels = ["0", x_vals_labels];

% Plot colors for each basis set
color = [0 0.4470 0.7410;
         0.8500 0.3250 0.0980;
         0.4940 0.1840 0.5560];

for simulation_type = simulation_types
    for data_type = data_types
        figure("Name", simulation_type + " Structure - "+ data_type + " Energy", ...
               "Color", "w"); % White background
        hold on; grid on;

        iter = eval(strjoin(["data" simulation_type data_type], "_"));

        for i = 1:length(basis_sets)
            y_vals = iter(i, :);
            x_data = x_vals(2:end); % exclude x=0 for fitting

            % Linear fit using polyfit (1st degree)
            p_all = polyfit(x_data, y_vals, 1);  % using all 3 points
            p_exclude_qz = polyfit(x_data(2:end), y_vals(2:end), 1); % exclude qz
            p_exclude_dz = polyfit(x_data(1:end-1), y_vals(1:end-1), 1); % exclude dz

            % Evaluate over full x_vals
            y_trend_all = polyval(p_all, x_vals);
            y_trend_excl_dz = polyval(p_exclude_dz, x_vals);
            y_trend_excl_qz = polyval(p_exclude_qz, x_vals);

            % Plot original data
            plot(x_data, y_vals, '-o', ...
                'LineWidth', 1.5, ...
                'DisplayName', k_vals(i) + " data", ...
                'Color', color(i, :));

            % Plot full fit line
            plot(x_vals, y_trend_all, '--', ...
                'LineWidth', 1.2, ...
                'DisplayName', k_vals(i) + " fit (all)", ...
                'Color', color(i, :));

            % Plot extrapolated value at x = 0
            scatter(0, y_trend_all(1), 60, '+', ...
                'MarkerEdgeColor', color(i, :), ...
                'DisplayName', k_vals(i) + " extrap. (all)");

            % Plot excluded-dz fit line
            plot(x_vals, y_trend_excl_dz, ':', ...
                'LineWidth', 1.2, ...
                'DisplayName', k_vals(i) + " fit (excl. dz)", ...
                'Color', color(i, :));

            % Plot extrapolated value at x = 0 (excl. dz)
            scatter(0, y_trend_excl_dz(1), 60, '*', ...
                'MarkerEdgeColor', color(i, :), ...
                'DisplayName', k_vals(i) + " extrap. (excl. dz)");

            % Plot excluded-qz fit line
            plot(x_vals, y_trend_excl_qz, ':', ...
                'LineWidth', 1.2, ...
                'DisplayName', k_vals(i) + " fit (excl. qz)", ...
                'Color', color(i, :));

            % Plot extrapolated value at x = 0 (excl. qz)
            scatter(0, y_trend_excl_qz(1), 60, '*', ...
                'MarkerEdgeColor', color(i, :), ...
                'DisplayName', k_vals(i) + " extrap. (excl. qz)");
        end

        xlabel('K-points mesh (e.g. $k \times k \times k \rightarrow k^{-3}$)', 'Interpreter', 'latex');
        xticks(x_vals);
        xticklabels(x_vals_labels);
        ylabel(data_type + " Energy (kJ/mol)");
        title(data_type + " Energy vs K-points for Different Basis Sets");
        legend('Location', 'bestoutside', 'Interpreter', 'none');
        hold off;

        % Save figure
        filename = strjoin(["data" simulation_type data_type], "_") + ".png";
        exportgraphics(gcf, filename, 'Resolution', 300);  % High-resolution export
    end
end
