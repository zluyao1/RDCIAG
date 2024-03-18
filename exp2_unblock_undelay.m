function exp2_unblock_undelay(n,m,s,lambda,maxiter,num_repeats)

fprintf('Experiment A, n=%d, m=%d, s=%d, maxiter=%d,',n,m,s,maxiter)

err_05 = zeros(maxiter,num_repeats);
err_1 = zeros(maxiter,num_repeats);
err_15 = zeros(maxiter,num_repeats);
err_2 = zeros(maxiter,num_repeats);

fprintf('experiment ')
% Loop over number of repeats
rng(69462991)
for repeats = 1:num_repeats
    fprintf('# %d, ',repeats)
    % Set up instance
    A = randn(m,n);   
    xhat = sparserandn(n,s);  % true solution
    b = A*xhat;  
    
    x_05 = zeros(n,1);
    y_05 = zeros(n,m); 
    x_1 = zeros(n,1);
    y_1 = zeros(n,m); 
    x_15 = zeros(n,1);
    y_15 = zeros(n,m); 
    x_2 = zeros(n,1);
    y_2 = zeros(n,m); 
    
    S = @(z) max(abs(z)-lambda,0).*sign(z);
    project = @(x,a,b) x-(a'*x-b)*a/norm(a)^2;

    % rand
    alpha = 0.5;
    prox = @(x,a,b) x-alpha*project(x/alpha,a,b);
    for iter = 1:maxiter
        
        for i = 1:n
            x_05(i) = S(-sum(y_05(i,:)));
        end
        index = randperm(m,1);
        y_05(:,index) = prox(y_05(:,index) + alpha*x_05, A(index,:)', b(index));

        err_05(iter,repeats) = norm(x_05-xhat);
    end
    
    % alpha = 1
    alpha = 1;
    prox = @(x,a,b) x-alpha*project(x/alpha,a,b);
    for iter = 1:maxiter
        
        for i = 1:n
            x_1(i) = S(-sum(y_1(i,:)));
        end
        index = randperm(m,1);
        y_1(:,index) = prox(y_1(:,index) + alpha*x_1, A(index,:)', b(index));

        err_1(iter,repeats) = norm(x_1-xhat);
    end

    % alpha = 1.5
    alpha = 1.5;
    prox = @(x,a,b) x-alpha*project(x/alpha,a,b);
    for iter = 1:maxiter
        
        for i = 1:n
            x_15(i) = S(-sum(y_15(i,:)));
        end
        index = randperm(m,1);
        y_15(:,index) = prox(y_15(:,index) + alpha*x_15, A(index,:)', b(index));

        err_15(iter,repeats) = norm(x_15-xhat);
    end
    
    % alpha = 2
    alpha = 2;
    prox = @(x,a,b) x-alpha*project(x/alpha,a,b);
    for iter = 1:maxiter
        
        for i = 1:n
            x_2(i) = S(-sum(y_2(i,:)));
        end
        index = randperm(m,1);
        y_2(:,index) = prox(y_2(:,index) + alpha*x_2, A(index,:)', b(index));

        err_2(iter,repeats) = norm(x_2-xhat);
    end    
end

% plot
subplot(1,1,1)
% mediumblue = [0 0.4470 0.7410];
% mediumyellow = [0.8500 0.3250 0.0980];
% lightblue = [0.6 0.6 1];
% lightyellow = [1 1 0];

lightblue = [0.6 0.6 1];
mediumblue = [0.3 0.3 1];
lightyellow = [1 1 0];
mediumyellow = [1 1 0];


lightred = [1 0.9 0.9];
mediumred = [1 0.6 0.6];
lightgreen = [0.9 1 0.9];
mediumgreen = [0.6 1 0.6];

meanerr_05 = median(err_05,2);
meanerr_1 = median(err_1,2);
meanerr_15 = median(err_15,2);
meanerr_2 = median(err_2,2);

semilogy(1:maxiter,meanerr_05,'k','LineWidth',3), hold on
semilogy(1:maxiter,meanerr_1,'k','LineWidth',3), hold on
semilogy(1:maxiter,meanerr_15,'k','LineWidth',3), hold on
semilogy(1:maxiter,meanerr_2,'k','LineWidth',3), hold on

hold off
ax = get(gca);
set(gca, 'linewidth', 3, 'fontsize', 30, 'fontname', 'times');

LL = find(meanerr_05<10^(-20));
if isempty(LL) == 1
    LL = maxiter;
else
LL = LL(1);
end
nn = 1 : LL;
focus = int32(floor(linspace(0, LL-1,11))+1);
Markers = {'ks';'rp';'go';'bd';'y^';'m'};

plot(1:maxiter,log10(meanerr_05),'Color',mediumblue,'LineStyle','-','LineWidth',3),hold on
plot(nn(focus),log10(meanerr_05(focus)),Markers{1},'MarkerSize',5,'MarkerFaceColor',lightblue),hold on

plot(1:maxiter,log10(meanerr_1),'Color',mediumred,'LineStyle','-','LineWidth',3),hold on
plot(nn(focus),log10(meanerr_1(focus)),Markers{2},'MarkerSize',5,'MarkerFaceColor',lightred),hold on

plot(1:maxiter,log10(meanerr_15),'Color',mediumgreen,'LineStyle','-','LineWidth',3),hold on
plot(nn(focus),log10(meanerr_15(focus)),Markers{3},'MarkerSize',5,'MarkerFaceColor',lightgreen),hold on

plot(1:maxiter,log10(meanerr_2),'Color',mediumyellow,'LineStyle','-','LineWidth',3),hold on
plot(nn(focus),log10(meanerr_2(focus)),Markers{4},'MarkerSize',5,'MarkerFaceColor',lightyellow),hold on

legend('\alpha = 0.5','','\alpha = 1','','\alpha = 1.5','','\alpha = 2','','Location','NorthEast','Fontsize', 12);
%title('Dual objective')
xl = xlabel('iteration','fontname','Times'); 
set(xl,'Fontsize',30)
yl = ylabel('log(Error)');
set(yl,'Fontsize',30,'fontname','Times')
hold off


yticks = log10(ax.YTick);
newyticks = cell(length(yticks));
for i=1:length(yticks)
    newyticks(i)=cellstr(strcat('10^{',num2str(yticks(i))));
    newyticks(i)=cellstr(strcat(newyticks(i),'}'));
end

set(gca,'YTick',log10(10.^yticks));
set(gca,'YTicklabel',newyticks);
tit = title([sprintf('m = %d, n = %d, s=%d',m,n,s)]);
set(tit,'FontSize',30,'fontname','Times');
grid on;

end