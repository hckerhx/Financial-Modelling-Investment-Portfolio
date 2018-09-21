function [] = plotWeight(x,plttitle,legends)
% x is already transposed
area(x)
legend(legends, 'Location', 'eastoutside','FontSize',12);
title(plttitle, 'FontSize', 14);
ylabel('Weights','interpreter','latex','FontSize',12);
xlabel('Rebalance Period','interpreter','latex','FontSize',12);
ylim([-0.2 1.4])

% Define the plot size in inches
% set(gcf,'Units','Inches', 'Position', [0 0 6, 3]);
% pos1 = get(gcf,'Position');
% set(gcf,'PaperPositionMode','Auto','PaperUnits','Inches',...
%     'PaperSize',[pos1(3), pos1(4)]);
end