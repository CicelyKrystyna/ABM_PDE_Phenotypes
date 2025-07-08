clear all;

mean_pheno_100 = importdata('./output/constant_02_100_mean_phenotype.txt');
endval = length(mean_pheno_100)-1
x = 0:1:endval;

fig1=figure;
plot(x*1000,mean_pheno_100, 'Color', 'Black', 'LineWidth', 2, 'LineStyle', '-');
hold on;
yline(20/23, 'LineWidth',2, 'Color', 'Red', 'Linestyle', '--');
set(gca,'FontSize',20);
axis([0 10000000 0.0 1.0]);
Xl=xlabel('timestep', 'Interpreter', 'latex');
Yl=ylabel('average phenotype', 'Interpreter', 'latex');
saveas(fig1,'average_phenotype_o2_100.png');

mean_pheno_20 = importdata('./output/constant_02_20_mean_phenotype.txt');
endval = length(mean_pheno_20)-1
x = 0:1:endval;

fig2=figure;
plot(x*1000,mean_pheno_20, 'Color', 'Black', 'LineWidth', 2, 'LineStyle', '-');
hold on;
yline(20/27, 'LineWidth',2, 'Color', 'Red', 'Linestyle', '--');
set(gca,'FontSize',20);
axis([0 10000000 0.0 1.0]);
Xl=xlabel('timestep', 'Interpreter', 'latex');
Yl=ylabel('average phenotype', 'Interpreter', 'latex');
saveas(fig2,'average_phenotype_o2_20.png');

mean_pheno_6 = importdata('./output/constant_02_6_mean_phenotype.txt');
endval = length(mean_pheno_6)-1
x = 0:1:endval;

fig3=figure;
plot(x*1000,mean_pheno_6, 'Color', 'Black', 'LineWidth', 2, 'LineStyle', '-');
hold on;
yline(15/29, 'LineWidth',2, 'Color', 'Red', 'Linestyle', '--');
set(gca,'FontSize',20);
axis([0 10000000 0.0 1.0]);
Xl=xlabel('timestep', 'Interpreter', 'latex');
Yl=ylabel('average phenotype', 'Interpreter', 'latex');
saveas(fig3,'average_phenotype_o2_6.png');

mean_pheno_1 = importdata('./output/constant_02_1_mean_phenotype.txt');
endval = length(mean_pheno_1)-1
x = 0:1:endval;

fig4=figure;
plot(x*1000,mean_pheno_1, 'Color', 'Black', 'LineWidth', 2, 'LineStyle', '-');
hold on;
yline(10/61, 'LineWidth',2, 'Color', 'Red', 'Linestyle', '--');
set(gca,'FontSize',20);
axis([0 10000000 0.0 1.0]);
Xl=xlabel('timestep', 'Interpreter', 'latex');
Yl=ylabel('average phenotype', 'Interpreter', 'latex');
saveas(fig4,'average_phenotype_o2_1.png');