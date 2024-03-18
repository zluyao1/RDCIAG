function exp2_unblock_delay_tau2(n,m,s,alpha,lambda,maxiter,num_repeats)

fprintf('Experiment A, n=%d, m=%d, s=%d, maxiter=%d,',n,m,s,maxiter)

err_10 = zeros(maxiter,num_repeats);
err_50 = zeros(maxiter,num_repeats);
err_100 = zeros(maxiter,num_repeats);

fprintf('experiment ')
% Loop over number of repeats
rng(69462991)
for repeats = 1:num_repeats
    fprintf('# %d, ',repeats)
    % Set up instance
    A = randn(m,n);   
    xhat = sparserandn(n,s);  % true solution
    b = A*xhat;  
    
    tau = 10;
    x_10 = zeros(n,1);
    y_10 = zeros(n,m,tau+2);  
    tau = 50;
    x_50 = zeros(n,1);
    y_50 = zeros(n,m,tau+2);
    tau = 100;
    x_100 = zeros(n,1);
    y_100 = zeros(n,m,tau+2);
    
    S = @(z) max(abs(z)-lambda,0).*sign(z);
    project = @(x,a,b) x-(a'*x-b)*a/norm(a)^2;
    prox = @(x,a,b) x-alpha*project(x/alpha,a,b);

    % rand
    tau = 10;
    for iter = 1:maxiter
        
        for i = 1:n
            value = randperm(tau+1,1);
            yever = y_10(:,:,value);            
            x_10(i) = S(-sum(yever(i,:)));
        end
        index = randperm(m,1);
        ynow = y_10(:,:,tau+1);
        ynow(:,index) = prox(ynow(:,index) + alpha*x_10, A(index,:)', b(index));

        y_10(:,:,tau+2) = ynow;
        y_10(:,:,(1:(tau+1))) = y_10(:,:,(2:(tau+2)));

        err_10(iter,repeats) = norm(x_10-xhat);
    end
    
    % tau = 50
    tau = 50;
    for iter = 1:maxiter
        
        for i = 1:n
            value = randperm(tau+1,1);
            yever = y_50(:,:,value);            
            x_50(i) = S(-sum(yever(i,:)));
        end
        index = randperm(m,1);
        ynow = y_50(:,:,tau+1);
        ynow(:,index) = prox(ynow(:,index) + alpha*x_50, A(index,:)', b(index));

        y_50(:,:,tau+2) = ynow;
        y_50(:,:,(1:(tau+1))) = y_50(:,:,(2:(tau+2)));

        err_50(iter,repeats) = norm(x_50-xhat);
    end

    % tau = 100
    tau = 100;
    for iter = 1:maxiter
        
        for i = 1:n
            value = randperm(tau+1,1);
            yever = y_100(:,:,value);            
            x_100(i) = S(-sum(yever(i,:)));
        end
        index = randperm(m,1);
        ynow = y_100(:,:,tau+1);
        ynow(:,index) = prox(ynow(:,index) + alpha*x_100, A(index,:)', b(index));

        y_100(:,:,tau+2) = ynow;
        y_100(:,:,(1:(tau+1))) = y_100(:,:,(2:(tau+2)));

        err_100(iter,repeats) = norm(x_100-xhat);
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


lightred = [1 0.9 0.9];
mediumred = [1 0.6 0.6];
lightgreen = [0.9 1 0.9];
mediumgreen = [0.6 1 0.6];

meanerr_10 = median(err_10,2);
meanerr_50 = median(err_50,2);
meanerr_100 = median(err_100,2);

semilogy(1:maxiter,meanerr_10,'k','LineWidth',3), hold on
semilogy(1:maxiter,meanerr_50,'k','LineWidth',3), hold on
semilogy(1:maxiter,meanerr_100,'k','LineWidth',3), hold on

hold off
ax = get(gca);
set(gca, 'linewidth', 3, 'fontsize', 30, 'fontname', 'times');

LL = find(meanerr_10 < 1e-20);
if isempty(LL) == 1
    LL = maxiter;
else
LL = LL(1);
end
nn = 1 : LL;
focus = int32(floor(linspace(0, LL-1,11))+1);
Markers = {'ks';'rp';'go';'bd';'y^';'m'};

plot(1:maxiter,log10(meanerr_10),'Color',mediumblue,'LineStyle','-','LineWidth',3),hold on
plot(nn(focus),log10(meanerr_10(focus)),Markers{1},'MarkerSize',5,'MarkerFaceColor',lightblue),hold on

plot(1:maxiter,log10(meanerr_50),'Color',mediumred,'LineStyle','-','LineWidth',3),hold on
plot(nn(focus),log10(meanerr_50(focus)),Markers{2},'MarkerSize',5,'MarkerFaceColor',lightred),hold on

plot(1:maxiter,log10(meanerr_100),'Color',mediumgreen,'LineStyle','-','LineWidth',3),hold on
plot(nn(focus),log10(meanerr_100(focus)),Markers{3},'MarkerSize',5,'MarkerFaceColor',lightgreen),hold on

legend('\tau = 10','','\tau = 50','','\tau = 100','','Location','NorthEast','Fontsize', 12);
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