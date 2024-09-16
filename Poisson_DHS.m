clear
pATh       = "/Users/lshjr3/Documents/DHS_MPselection/";

options    = detectImportOptions(char(pATh + "RakingOuTpUT.csv"));
for j = 1:numel(options.VariableTypes)
    if isequal(options.VariableTypes{j},'char')
        options.VariableTypes{j} = 'categorical';
    end
    if isequal(options.VariableTypes{j},'datetime')
        options.VariableOptions(1,j).InputFormat = 'dd/MM/yyyy';
    end
end
dATa       = readtable(char(pATh + "RakingOuTpUT.csv"),options);
lISt       = tabulate(dATa.survey);
lISt       = lISt(:,1);
R          = 2500;
for i = 1:numel(lISt)
    sURvEy{i}       = dATa(dATa.survey == lISt{i},:);
    woman           = tabulate(string(sURvEy{i}.cluster) + " - " + string(sURvEy{i}.caseid));
    woman           = cell2mat(woman(:,2));
    K               = cumsum([1 woman(1:end - 1)']');
    
    cluster         = tabulate(string(sURvEy{i}.cluster(K)));
    cluster         = cell2mat(cluster(:,2));
    rng(0);
    for j = 1:numel(cluster)
        S{j,1} = [(1:cluster(j))',unidrnd(cluster(j),cluster(j),R)] + sum(cluster(1:j - 1));
    end
    S               = [cell2mat(S);ones(1,R + 1)*(sum(cluster) + 1)];
    for r = 1:size(S,2)
        temp    = tabulate(S(:,r));
        w(:,r)  = repelem(temp(1:end - 1,2),woman).*sURvEy{i}.W;
        wr(:,r) = repelem(temp(1:end - 1,2),woman).*recode(sURvEy{i}.WR,NaN,0);
        clear temp
        clc;
        r/(R + 1)
    end
    wEIgHt{i}       = w;
    wEIgHtR{i}      = wr;
    exposure{i}     = recode(table2array(sURvEy{i}(:,{'exposure15','exposure20','exposure25','exposure30','exposure35','exposure40','exposure45','exposure50','exposure55'})),NaN,0);
    events{i}       = recode(table2array(sURvEy{i}(:,{'events15','events20','events25','events30','events35','events40','events45','events50','events55'})),NaN,0);     
    X{i}            = sURvEy{i}(:,{'sexS','mobile','ageG','Education','Marital','Q','Region','UR','sex_head'});
    T.pOP{i,1}      = {char(sURvEy{i}.name(1))};
    T.clusters{i,1} = [numel(cluster) NaN NaN];
    clear cluster woman S K w wr r ans
end

for i = 1:numel(lISt)
    sEL             = [ones(size(sURvEy{i}.mobile)) sURvEy{i}.mobile];
    for j = 1:size(sEL,2)
        m{i}{j}  = (events{i}'*(wEIgHt{i}.*sEL(:,j)))./(exposure{i}'*(wEIgHt{i}.*sEL(:,j)));
    end
    m{i}{3}         = (events{i}'*wEIgHtR{i})./(exposure{i}'*wEIgHtR{i});
    
    for j = 1:numel(m{i})
        q{i}{j}  = 1 - prod(1 - 5*m{i}{j}./(1 + (5 - 2.5)*m{i}{j}));
    end
    
    sex             = [(sURvEy{i}.sexS == 1) (sURvEy{i}.sexS == 2)];
    sEL             = [sex.*sEL(:,1) sex.*sEL(:,2)];
    for j = 1:size(sEL,2)
        ms{i}{j} = (events{i}'*(wEIgHt{i}.*sEL(:,j)))./(exposure{i}'*(wEIgHt{i}.*sEL(:,j)));
    end
    ms{i}{5}        = (events{i}'*(wEIgHtR{i}.*(sURvEy{i}.sexS == 1)))./(exposure{i}'*(wEIgHtR{i}.*(sURvEy{i}.sexS == 1)));
    ms{i}{6}        = (events{i}'*(wEIgHtR{i}.*(sURvEy{i}.sexS == 2)))./(exposure{i}'*(wEIgHtR{i}.*(sURvEy{i}.sexS == 2)));
    
    for j = 1:numel(ms{i})
        qs{i}{j} = 1 - prod(1 - 5*ms{i}{j}./(1 + (5 - 2.5)*ms{i}{j}));
    end
end

for i = 1:numel(lISt)    
    sEL             = find(sum(exposure{i},2) > 0);    
    w               = (exposure{i} > 0).*wEIgHt{i}(:,1);
    w               = (w./sum(w,1)).*sum(exposure{i} > 0);
    
    for j = 20:5:45
        ageG(:,j/5 - 3) = (X{i}.ageG >= j);
    end
    
    for j = 1:max(X{i}.Region) - 1
        Region(:,j) = (X{i}.Region == j + 1);
    end
    
    covariates1     = (X{i}.mobile == 1);
    covariates1     = mat2cell(kron(ones(1,size(events{i},2)),covariates1(sEL,:)),numel(sEL),ones(1,size(events{i},2))*size(covariates1,2));
    covariates2     = [X{i}.mobile == 1,X{i}.sexS == 1,X{i}.Education >= 2,X{i}.Education >= 3,X{i}.Marital == 1,X{i}.Q == 2,X{i}.UR == 1,X{i}.sex_head == 1,ageG,Region];
    covariates2     = mat2cell(kron(ones(1,size(events{i},2)),covariates2(sEL,:)),numel(sEL),ones(1,size(events{i},2))*size(covariates2,2));
        
    mo.events       = events{i}(sEL,:);
    mo.exposure     = exposure{i}(sEL,:);
    mo.w            = w(sEL,:);
    
    mo.type         = 'Poisson';
    mo.X            = covariates1;
    mo.Beta         = [NaN(size(events{i},2),1);NaN(size(mo.X{1},2),1)];
    mo.sEL          = find(isnan(mo.Beta));
    [f{i}{1},B,~,s] = CG_PoissonNB(mo.events,mo.exposure,mo.X,mo.w,1000,mo.type,mo.Beta);
    B               = B + s*norminv([0.50 0.025 0.975],0,1);
    Poisson1{i}     = mat2cell(B,[size(mo.events,2) 1],3);

    mo.X            = covariates2;
    mo.Beta         = [NaN(size(events{i},2),1);NaN(size(mo.X{1},2),1)];
    mo.sEL          = find(isnan(mo.Beta));
    [f{i}{j},B,~,s] = CG_PoissonNB(mo.events,mo.exposure,mo.X,mo.w,1000,mo.type,mo.Beta);
    B               = B + s*norminv([0.50 0.025 0.975],0,1);
    Poisson2{i}     = mat2cell(B,[size(mo.events,2) 1 1 2 1 1 1 1 size(ageG,2) size(Region,2)],3);

    T.women{i,1}    = (sURvEy{i}.k <= 1)'*wEIgHt{i};
    T.age{i,1}      = (((sURvEy{i}.k <= 1).*sURvEy{i}.age)'*wEIgHt{i})./T.women{i,1};
    T.siblings{i,1} = ((sURvEy{i}.k >= 1)'*wEIgHt{i})./T.women{i,1};
    T.SRS{i,1}      = (((sURvEy{i}.k >= 1).*(sURvEy{i}.sexS == 1))'*wEIgHt{i})./(((sURvEy{i}.k >= 1).*(sURvEy{i}.sexS == 2))'*wEIgHt{i});
    T.urban{i,1}    = (((sURvEy{i}.k <= 1).*(sURvEy{i}.UR == 1))'*wEIgHt{i})./T.women{i,1}; 
    T.mobile{i,1}   = (((sURvEy{i}.k <= 1).*(sURvEy{i}.mobile == 1))'*wEIgHt{i})./T.women{i,1};
    clear Region ageG sEL covariates1 covariates2 w j B s
    clc;
    i
end
save(char(pATh + "database.mat"),'-v7.3');

for i = 1:numel(lISt)
    sORt(i,:) = [mean(T.mobile{i}(2:end)) i];
end
sORt       = sortrows(sORt);
sORt       = sORt(:,2);

for i = 1:numel(lISt)
    bOxG{i,1}       = T.clusters{sORt(i)};
    bOxG{i,2}       = prctile(round(T.women{sORt(i)}(2:end),1),[50 2.5 97.5]);
    bOxG{i,3}       = prctile(T.age{sORt(i)}(2:end),[50 2.5 97.5]);
    bOxG{i,4}       = prctile(T.siblings{sORt(i)}(2:end),[50 2.5 97.5]);
    bOxG{i,5}       = prctile(T.SRS{sORt(i)}(2:end),[50 2.5 97.5]);
    bOxG{i,6}       = prctile(T.urban{sORt(i)}(2:end),[50 2.5 97.5])*100;
    bOxG{i,7}       = prctile(T.mobile{sORt(i)}(2:end),[50 2.5 97.5])*100;
end
    
sEt        = {{'Sample'},{'Descriptive Statistics'}};
vARs1      = {'$\textit{clusters}$','$\textit{women}$'};
vARs2      = {'$\textit{average age}$','$\textit{average siblings}$','$\textit{sex ratio}$','$\textit{urban}$ $\mathrm{(}\mathrm{\%}\mathrm{)}$','$\textit{mobile owners}$ $\mathrm{(}\mathrm{\%}\mathrm{)}$'};

foRMaT     = {'%0.0f','%0.0f','%0.2f','%0.2f','%0.2f','%0.2f','%0.2f'};
lABs       = {{1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25}};
nOTe       = {'Survey',char("$\mathrm{p50}$/$\mathit{[p2.5,p97.5],}$ Bootstrapping each survey within clusters " + string(R) + " times with replacement.")};
tABleBAyEs(sEt,{vARs1,vARs2},foRMaT,lABs,nOTe,T.pOP(sORt),cell2mat(bOxG),0.075,0.100,[]);
saveas(gcf,char(pATh + "rev Table 1.png"));


for i = 1:numel(lISt)
    bOx{i,1}        = prctile(T.mobile{sORt(i)}(2:end),[50 2.5 97.5])*100;
    bOx{i,2}        = prctile(q{sORt(i)}{1}(2:end),[50 2.5 97.5])*1000;
    bOx{i,3}        = prctile(log(q{sORt(i)}{2}(2:end)./q{sORt(i)}{1}(2:end)),[50 2.5 97.5]);
    bOx{i,4}        = Poisson1{sORt(i)}{2};
    bOx{i,5}        = Poisson2{sORt(i)}{2};
    bOx{i,6}        = prctile(log(q{sORt(i)}{3}(2:end)./q{sORt(i)}{1}(2:end)),[50 2.5 97.5]);    
end

sEt        = {{'Observed (life table estimates)'},{'Poisson (excess mortality in ln)'},{'Post-stratified Owners'}};
vARs1      = {'$\textit{mobile owners}$ $\mathrm{(}\mathrm{\%}\mathrm{)}$','$\mathrm{_4}\mathrm{_5}\mathit{q}\mathrm{_1}\mathrm{_5}$ $\mathrm{000}$ $\mathrm{(}\mathit{total}\mathrm{)}$','$\mathrm{ln[}\mathit{owners/total}\mathrm{]}$'};
vARs2      = {'$\textit{without controls}$','$\textit{with controls}$'};
vARs3      = {'$\mathrm{ln[}\mathit{owners/total}\mathrm{]}$'};

foRMaT     = {'%0.2f','%0.2f','%0.4f','%0.4f','%0.4f','%0.4f'};
nOTe       = {'Survey',char("$\mathrm{p50}$/$\mathit{[p2.5,p97.5],}$ Bootstrapping each survey within clusters " + string(R) + " times with replacement. Poisson is ML estimation and CIs are calculated from the Hessian matrix, under normality assumptions.")};
tABleBAyEs(sEt,{vARs1,vARs2,vARs3},foRMaT,lABs,nOTe,T.pOP(sORt),cell2mat(bOx),0.075,0.130,[]);
saveas(gcf,char(pATh + "rev Table 2.png"));

for i = 1:numel(lISt)
    bOxS{i,1}       = prctile(log(qs{sORt(i)}{1}(2:end)./qs{sORt(i)}{2}(2:end)),[50 2.5 97.5]);
    bOxS{i,2}       = prctile(log(qs{sORt(i)}{3}(2:end)./qs{sORt(i)}{4}(2:end)),[50 2.5 97.5]);
    bOxS{i,3}       = Poisson2{sORt(i)}{3}(1,:);
    bOxS{i,4}       = prctile(log(qs{sORt(i)}{5}(2:end)./qs{sORt(i)}{6}(2:end)),[50 2.5 97.5]);    
end

sEt        = {{'Observed $\mathrm{_4}\mathrm{_5}\mathit{q}\mathrm{_1}\mathrm{_5}$'},{'Poisson (excess mortality in ln)'},{'Post-stratified $\mathrm{_4}\mathrm{_5}\mathit{q}\mathrm{_1}\mathrm{_5}$'}};
vARs1      = {'$\mathrm{ln[}\mathit{male/female}\mathrm{]\,(total)}$','$\mathrm{ln[}\mathit{male/female}\mathrm{]\,(owners)}$'};
vARs2      = {'$\mathit{control:\,male=1}$'};
vARs3      = {'$\mathrm{ln[}\mathit{male/female}\mathrm{]\,(owners)}$'};

foRMaT     = {'%0.4f','%0.4f','%0.4f','%0.4f'};
nOTe       = {'Survey',char("$\mathrm{p50}$/$\mathit{[p2.5,p97.5],}$ Bootstrapping each survey within clusters " + string(R) + " times with replacement. Poisson is ML estimation and CIs are calculated from the Hessian matrix, under normality assumptions.")};
tABleBAyEs(sEt,{vARs1,vARs2,vARs3},foRMaT,lABs,nOTe,T.pOP(sORt),cell2mat(bOxS),0.075,0.175,[]);
saveas(gcf,char(pATh + "rev Table 3.png"));


LAB                      = {'$\textit{total}$','$\textit{mobile\,owners}$','$\textit{post-stratified}$'};
coloR                    = {[0.00 0.00 0.75],[0.95 0.00 0.95],[0.45 0.65 0.20],[0.85 0.35 0.01],[0.65 0.10 0.20],[0.00 0.55 0.65],[0.05 0.05 0.05]};
pix                      = 1/37.7952755906;
fi                       = figure('Color',[1 1 1]);
fi.Position              = [0 0 21 42]/pix;
axes1                    = axes('Parent',fi,'Position',[0.025 0.025 0.975 0.975]);
hold(axes1,'on');
TL                       = tiledlayout(5,5,'Padding','compact','TileSpacing','compact');
for i = 1:numel(T.pOP)
    nexttile(i)
    ax{i}                       = gca;
    ax{i}.FontName              = 'Times New Roman';
    ax{i}.FontSize              = 8;
    ax{i}.YScale                = 'log';
    ax{i}.XAxis.TickLabelFormat = '%.0f';
    ax{i}.YAxis.TickLabelFormat = '%.1f';
    ax{i}.XTick                 = 15:5:60;
    ax{i}.XAxis.MinorTickValues = 15:1:60;    
    xlim([15 60])
    
    if i >= 21
        xlabel('$\mathit{age}$ (in years)','Interpreter','latex','FontName','Times New Roman','FontSize',9);
    end
    if ismember(i,1:5:25)
        ylabel('$\mathrm{_5}\mathit{M}\mathrm{_x}$ (log scale)','Interpreter','latex','FontName','Times New Roman','FontSize',9);
    end
    title(char(string(i) + ". " + string(T.pOP{sORt(i)}{1})),'Interpreter','latex');
    grid on;
    box on;
    hold on;
end

xn = (15:5:55)' + 2.5;
for i = 1:numel(T.pOP)
    nexttile(i)
    for j = 1:numel(m{sORt(i)})
        plot(xn,prctile(m{sORt(i)}{j}(:,2:end)',50)','color',coloR{j},'LineWidth',1.00);
    end
    for j = 1:numel(m{sORt(i)})
        fill([xn;flip(xn)],[prctile(m{sORt(i)}{j}(:,2:end)',2.5)';flip(prctile(m{sORt(i)}{j}(:,2:end)',97.5)')],coloR{j},'FaceAlpha',.10,'EdgeAlpha',0.25,'LineStyle',':','EdgeColor',coloR{j});
    end
    if isequal(i,23)
        legend(LAB,'Interpreter','latex','FontName','Times New Roman','FontSize',9,'FontAngle','oblique','Location','southoutside','NumColumns',3,'Box','off');
    end
end
saveas(gcf,char(pATh + "rev Figure 1.png"));

fi                       = figure('Color',[1 1 1]);
fi.Position              = [0 0 21 42]/pix;
axes1                    = axes('Parent',fi,'Position',[0.025 0.025 0.975 0.975]);
hold(axes1,'on');
TL                       = tiledlayout(5,5,'Padding','compact','TileSpacing','compact');
for i = 1:numel(T.pOP)
    nexttile(i)
    ax{i}                       = gca;
    ax{i}.FontName              = 'Times New Roman';
    ax{i}.FontSize              = 8;
    ax{i}.XScale                = 'log';
    ax{i}.XAxis.TickLabelFormat = '%.0f';
    ax{i}.YAxis.TickLabelFormat = '%.0f';
    
    if i >= 21
        xlabel('$\mathrm{_4}\mathrm{_5}\mathit{q}\mathrm{_1}\mathrm{_5}$ $\mathrm{000}$ (log scale)','Interpreter','latex','FontName','Times New Roman','FontSize',9);
    end    
    if ismember(i,1:5:25)
        ylabel('$\mathit{density}$','Interpreter','latex','FontName','Times New Roman','FontSize',9);
    end
    title(char(string(i) + ". " + string(T.pOP{sORt(i)}{1})),'Interpreter','latex');
    grid on;
    box on;
    hold on;
end

for i = 1:numel(T.pOP)
    nexttile(i)
    qw = q{sORt(i)};
    for j = 1:numel(qw)
        [fi,xi] = ksdensity(log(max(qw{j}(2:end),eps)));
        xi      = exp(xi)*1000;
        plot(xi,fi,'color',coloR{j},'LineWidth',1.00);
    end
    for j = 1:numel(qw)
        [fi,xi] = ksdensity(log(max(qw{j}(2:end),eps)));
        xi      = exp(xi)*1000;
        fill(xi,fi,coloR{j},'FaceAlpha',.10,'EdgeAlpha',0.25,'LineStyle',':','EdgeColor',coloR{j});
    end
    for j = 1:numel(qw)
        xi      = prctile(qw{j}(2:end),[50 2.5 97.5]);
        [fi,xi] = ksdensity(log(max(qw{j}(2:end),eps)),log(xi));
        xi      = exp(xi)*1000;
        plot(xi([1 1]),[0 fi(1)],':','color',coloR{j},'LineWidth',0.50);
        plot(xi([2 2]),[0 fi(2)],'color',coloR{j},'LineWidth',0.50);
        plot(xi([3 3]),[0 fi(3)],'color',coloR{j},'LineWidth',0.50);
    end
    if isequal(i,23)
        legend(LAB,'Interpreter','latex','FontName','Times New Roman','FontSize',9,'FontAngle','oblique','Location','southoutside','NumColumns',3,'Box','off');
    end
end
saveas(gcf,char(pATh + "rev Figure 2.png"));


LAB                      = {'$\textit{owners/total}$','$\textit{post-stratified/total}$'};
fi                       = figure('Color',[1 1 1]);
fi.Position              = [0 0 21 42]/pix;
axes1                    = axes('Parent',fi,'Position',[0.025 0.025 0.975 0.975]);
hold(axes1,'on');
TL                       = tiledlayout(5,5,'Padding','compact','TileSpacing','compact');
for i = 1:numel(T.pOP)
    nexttile(i)
    ax{i}                       = gca;
    ax{i}.FontName              = 'Times New Roman';
    ax{i}.FontSize              = 8;
    ax{i}.XAxis.TickLabelFormat = '%.2f';
    ax{i}.YAxis.TickLabelFormat = '%.1f';
    if i >= 21
        xlabel('$\mathrm{_4}\mathrm{_5}\mathit{q}\mathrm{_1}\mathrm{_5}$ $\mathit{(ratio)}$','Interpreter','latex','FontName','Times New Roman','FontSize',9);
    end    
    if ismember(i,1:5:25)
        ylabel('$\mathit{density}$','Interpreter','latex','FontName','Times New Roman','FontSize',9);
    end
    title(char(string(i) + ". " + string(T.pOP{sORt(i)}{1})),'Interpreter','latex');
    grid on;
    box on;
    hold on;
end

for i = 1:numel(T.pOP)
    nexttile(i)
    qw = q{sORt(i)};
    mx = 0;
    for j = 2:numel(qw)
        [fi,xi] = ksdensity(log(max(qw{j}(2:end),eps)./max(qw{1}(2:end),eps)));
        xi      = exp(xi);
        mx      = max([fi,mx]);
        plot(xi,fi,'color',coloR{j},'LineWidth',1.00);
    end
    for j = 2:numel(qw)
        [fi,xi] = ksdensity(log(max(qw{j}(2:end),eps)./max(qw{1}(2:end),eps)));
        xi      = exp(xi);
        fill(xi,fi,coloR{j},'FaceAlpha',.10,'EdgeAlpha',0.25,'LineStyle',':','EdgeColor',coloR{j});
    end
    for j = 2:numel(qw)
        xi      = prctile(qw{j}(2:end)./qw{1}(2:end),[50 2.5 97.5]);
        [fi,xi] = ksdensity(log(max(qw{j}(2:end),eps)./max(qw{1}(2:end),eps)),log(xi));
        xi      = exp(xi);
        plot(xi([1 1]),[0 fi(1)],':','color',coloR{j},'LineWidth',0.50);
        plot(xi([2 2]),[0 fi(2)],'color',coloR{j},'LineWidth',0.50);
        plot(xi([3 3]),[0 fi(3)],'color',coloR{j},'LineWidth',0.50);
    end    
    plot([1 1],[0 mx*1.05],'color','k','LineWidth',0.75);
    ylim([0 mx*1.05])
    if isequal(i,23)
        legend(LAB,'Interpreter','latex','FontName','Times New Roman','FontSize',9,'FontAngle','oblique','Location','southoutside','NumColumns',3,'Box','off');
    end
end
saveas(gcf,char(pATh + "rev Figure 2B.png"));


LAB                      = {'$\textit{owners/total}$','$\textit{post-stratified/total}$'};
fi                       = figure('Color',[1 1 1]);
fi.Position              = [0 0 21 42]/pix;
axes1                    = axes('Parent',fi,'Position',[0.025 0.025 0.975 0.975]);
hold(axes1,'on');
TL                       = tiledlayout(5,5,'Padding','compact','TileSpacing','compact');
for i = 1:numel(T.pOP)
    nexttile(i)
    ax{i}                       = gca;
    ax{i}.FontName              = 'Times New Roman';
    ax{i}.FontSize              = 8;
    ax{i}.XAxis.TickLabelFormat = '%.0f';
    ax{i}.YAxis.TickLabelFormat = '%.1f';
    ax{i}.XTick                 = 15:5:60;
    ax{i}.XAxis.MinorTickValues = 15:1:60;    
    xlim([15 60])
    ylim([0.5 2.0])
    
    if i >= 21
        xlabel('$\mathit{age}$ (in years)','Interpreter','latex','FontName','Times New Roman','FontSize',9);
    end
    if ismember(i,1:5:25)
        ylabel('$\mathrm{_5}\mathit{M}\mathrm{_x}$ (ratio)','Interpreter','latex','FontName','Times New Roman','FontSize',9);
    end
    title(char(string(i) + ". " + string(T.pOP{sORt(i)}{1})),'Interpreter','latex');
    grid on;
    box on;
    hold on;
end

for i = 1:numel(T.pOP)
    nexttile(i)
    for j = 2:numel(m{sORt(i)})
        plot(xn,prctile(m{sORt(i)}{j}(:,2:end)'./m{sORt(i)}{1}(:,2:end)',50)','color',coloR{j},'LineWidth',1.00);
    end
    for j = 2:numel(m{sORt(i)})
        fill([xn;flip(xn)],[prctile(m{sORt(i)}{j}(:,2:end)'./m{sORt(i)}{1}(:,2:end)',2.5)';flip(prctile(m{sORt(i)}{j}(:,2:end)'./m{sORt(i)}{1}(:,2:end)',97.5)')],coloR{j},'FaceAlpha',.10,'EdgeAlpha',0.25,'LineStyle',':','EdgeColor',coloR{j});
    end
    plot([min(xn - 2.5) max(xn + 2.5)],[1 1],'color','k','LineWidth',0.50);
    if isequal(i,23)
        legend(LAB,'Interpreter','latex','FontName','Times New Roman','FontSize',9,'FontAngle','oblique','Location','southoutside','NumColumns',3,'Box','off');
    end
end
saveas(gcf,char(pATh + "rev Figure 3.png"));

LAB                      = {'$\textit{total (male)}$','$\textit{total (female)}$','$\textit{mobile owners (male)}$','$\textit{mobile owners (female)}$','$\textit{post-stratified (male)}$','$\textit{post-stratified (female)}$'};
coloR                    = {[0.00 0.00 0.75],[0.95 0.00 0.95],[0.45 0.65 0.20],[0.85 0.35 0.01],[0.65 0.10 0.20],[0.00 0.55 0.65],[0.05 0.05 0.05]};
pix                      = 1/37.7952755906;
fi                       = figure('Color',[1 1 1]);
fi.Position              = [0 0 21 42]/pix;
axes1                    = axes('Parent',fi,'Position',[0.025 0.025 0.975 0.975]);
hold(axes1,'on');
TL                       = tiledlayout(5,5,'Padding','compact','TileSpacing','compact');
for i = 1:numel(T.pOP)
    nexttile(i)
    ax{i}                       = gca;
    ax{i}.FontName              = 'Times New Roman';
    ax{i}.FontSize              = 8;
    ax{i}.YScale                = 'log';
    ax{i}.XAxis.TickLabelFormat = '%.0f';
    ax{i}.YAxis.TickLabelFormat = '%.1f';
    ax{i}.XTick                 = 15:5:60;
    ax{i}.XAxis.MinorTickValues = 15:1:60;    
    xlim([15 60])
    
    if i >= 21
        xlabel('$\mathit{age}$ (in years)','Interpreter','latex','FontName','Times New Roman','FontSize',9);
    end
    if ismember(i,1:5:25)
        ylabel('$\mathrm{_5}\mathit{M}\mathrm{_x}$ (log scale)','Interpreter','latex','FontName','Times New Roman','FontSize',9);
    end
    title(char(string(i) + ". " + string(T.pOP{sORt(i)}{1})),'Interpreter','latex');
    grid on;
    box on;
    hold on;
end

for i = 1:numel(T.pOP)
    nexttile(i)
    for j = 1:numel(ms{sORt(i)})
        plot(xn,prctile(ms{sORt(i)}{j}(:,2:end)',50)','color',coloR{j},'LineWidth',1.00);
    end
    for j = 1:numel(ms{sORt(i)})
        fill([xn;flip(xn)],[prctile(ms{sORt(i)}{j}(:,2:end)',2.5)';flip(prctile(ms{sORt(i)}{j}(:,2:end)',97.5)')],coloR{j},'FaceAlpha',.10,'EdgeAlpha',0.25,'LineStyle',':','EdgeColor',coloR{j});
    end
    if isequal(i,23)
        legend(LAB,'Interpreter','latex','FontName','Times New Roman','FontSize',9,'FontAngle','oblique','Location','southoutside','NumColumns',3,'Box','off');
    end
end
saveas(gcf,char(pATh + "rev Figure 4.png"));


LAB                      = {'$\textit{total}$','$\textit{mobile owners}$','$\textit{post-stratified}$'};
fi                       = figure('Color',[1 1 1]);
fi.Position              = [0 0 21 42]/pix;
axes1                    = axes('Parent',fi,'Position',[0.025 0.025 0.975 0.975]);
hold(axes1,'on');
TL                       = tiledlayout(5,5,'Padding','compact','TileSpacing','compact');
for i = 1:numel(T.pOP)
    nexttile(i)
    ax{i}                       = gca;
    ax{i}.FontName              = 'Times New Roman';
    ax{i}.FontSize              = 8;
    ax{i}.XScale                = 'log';
    ax{i}.XAxis.TickLabelFormat = '%.1f';
    ax{i}.YAxis.TickLabelFormat = '%.1f';
    %xlim([.75 3.0])
    %ylim([0 5])
    if i >= 21
        xlabel('$\mathrm{_4}\mathrm{_5}\mathit{q}\mathrm{_1}\mathrm{_5}$ $\mathit{(sex\,ratio)}$','Interpreter','latex','FontName','Times New Roman','FontSize',9);
    end    
    if ismember(i,1:5:25)
        ylabel('$\mathit{density}$','Interpreter','latex','FontName','Times New Roman','FontSize',9);
    end
    title(char(string(i) + ". " + string(T.pOP{sORt(i)}{1})),'Interpreter','latex');
    grid on;
    box on;
    hold on;
end

for i = 1:numel(T.pOP)
    nexttile(i)
    qw = qs{sORt(i)};
    for j = 1:2:numel(qw)
        [fi,xi] = ksdensity(log(max(qw{j}(2:end),eps)./max(qw{j + 1}(2:end),eps)));
        plot(exp(xi),fi,'color',coloR{j},'LineWidth',1.00);
    end
    for j = 1:2:numel(qw)
        [fi,xi] = ksdensity(log(max(qw{j}(2:end),eps)./max(qw{j + 1}(2:end),eps)));
        fill(exp(xi),fi,coloR{j},'FaceAlpha',.10,'EdgeAlpha',0.25,'LineStyle',':','EdgeColor',coloR{j});
    end
    for j = 1:2:numel(qw)
        xi      = prctile(max(qw{j}(2:end),eps)./max(qw{j + 1}(2:end),eps),[50 2.5 97.5]);
        [fi,xi] = ksdensity(log(max(qw{j}(2:end),eps)./max(qw{j + 1}(2:end),eps)),log(xi));
        xi      = exp(xi);
        plot(xi([1 1]),[0 fi(1)],':','color',coloR{j},'LineWidth',0.50);
        plot(xi([2 2]),[0 fi(2)],'color',coloR{j},'LineWidth',0.50);
        plot(xi([3 3]),[0 fi(3)],'color',coloR{j},'LineWidth',0.50);        
    end
    
    if isequal(i,23)
        legend(LAB,'Interpreter','latex','FontName','Times New Roman','FontSize',9,'FontAngle','oblique','Location','southoutside','NumColumns',3,'Box','off');
    end
end
saveas(gcf,char(pATh + "rev Figure 5.png"));


fi                       = figure('Color',[1 1 1]);
fi.Position              = [0 0 21 42]/pix;
axes1                    = axes('Parent',fi,'Position',[0.025 0.025 0.975 0.975]);
hold(axes1,'on');
TL                       = tiledlayout(5,5,'Padding','compact','TileSpacing','compact');
for i = 1:numel(T.pOP)
    nexttile(i)
    ax{i}                       = gca;
    ax{i}.FontName              = 'Times New Roman';
    ax{i}.FontSize              = 8;
    ax{i}.XAxis.TickLabelFormat = '%.0f';
    ax{i}.YAxis.TickLabelFormat = '%.1f';
    ax{i}.XTick                 = 15:5:60;
    ax{i}.XAxis.MinorTickValues = 15:1:60;    
    xlim([15 60])
    ylim([0.5 3.0])
    
    if i >= 21
        xlabel('$\mathit{age}$ (in years)','Interpreter','latex','FontName','Times New Roman','FontSize',9);
    end
    if ismember(i,1:5:25)
        ylabel('$\mathrm{_5}\mathit{M}\mathrm{_x}$ (sex ratio)','Interpreter','latex','FontName','Times New Roman','FontSize',9);
    end
    title(char(string(i) + ". " + string(T.pOP{sORt(i)}{1})),'Interpreter','latex');
    grid on;
    box on;
    hold on;
end

for i = 1:numel(T.pOP)
    nexttile(i)
    for j = 1:2:numel(ms{sORt(i)})
        plot(xn,prctile(max(ms{sORt(i)}{j}(:,2:end),eps)'./max(ms{sORt(i)}{j + 1}(:,2:end),eps)',50)','color',coloR{j},'LineWidth',1.00);
    end
    for j = 1:2:numel(ms{sORt(i)})
        fill([xn;flip(xn)],[prctile(max(ms{sORt(i)}{j}(:,2:end),eps)'./max(ms{sORt(i)}{j + 1}(:,2:end),eps)',2.5)';flip(prctile(max(ms{sORt(i)}{j}(:,2:end),eps)'./max(ms{sORt(i)}{j + 1}(:,2:end),eps)',97.5)')],coloR{j},'FaceAlpha',.05,'EdgeAlpha',0.25,'LineStyle',':','EdgeColor',coloR{j});
    end
    plot([min(xn - 2.5) max(xn + 2.5)],[1 1],'color','k','LineWidth',0.50);
    if isequal(i,23)
        legend(LAB,'Interpreter','latex','FontName','Times New Roman','FontSize',9,'FontAngle','oblique','Location','southoutside','NumColumns',3,'Box','off');
    end
end
saveas(gcf,char(pATh + "rev Figure 6.png"));