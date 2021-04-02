clc
clear all
close all
tic

% format shortg
% profile viewer

%% IIR Filter Fitness

[Hfilt Wfilt]=IIR_main();

%% Initialize Parameters

k1damp = 0.004 ;
k2damp = 0.9   ; 

k1 = 0;
k2 = 0;

% glocal = 0: Global search, 1: Local search
% localdist: maximum distance used for choosing local balls within it
stallgenlimit = 200;
TolFun = -inf;
% Na = 'testfunc';
localnum = 3;
glocal = 0;
% cont = num2str(1);
% fitnessfunc = [Na, cont];
plots = 1;
numofruns = 200;
numofdims = 9;       % case1:2;  case2:4;  case3:9;
numofballs = 30;
%         c1 = .1024;
%         c2 = .6577;
%         shift1 = 40;
%         shift2 = 47;
%         scale1 = 0.55;
%         scale2 = 0.55;

Xmininit = repmat(0, 1, numofdims);     %case1:-2; case2:-2; case3:0;
Xmaxinit = repmat(1, 1, numofdims);     %case1:2; case2:1; case3:1;    

%% Loop Main

% fitnessfunc = str2func(Fitness);
bestfit = Inf; % stores previous total best
Neval = 0; % number of evaluating fitness function% 
% expanding Xmin, Xmax to cover all balls
Xmin = repmat(Xmininit, numofballs, 1);
Xmax = repmat(Xmaxinit, numofballs, 1);

% generating initial balls
% an option to generate custom initial balls can be added later
X = Xmin + (Xmax - Xmin) .* rand(numofballs, numofdims);

% initializing balls acceleration
A = zeros(numofballs, numofdims);

% initializes some variables for displaying the results
meanfits = zeros(numofruns, 1);
bests    = zeros(numofruns, 1);
worsts   = zeros(numofruns, 1);

%%%**********************************************************************
N       = rand(size(Hfilt,1),1)    ;
heights = Fitness(X,Hfilt,Wfilt,N) ;
%%%***********************************************************************

Neval = Neval + numofballs;
[tmpbestfit, tmpbestfitidx] = min(heights);
bestfit = tmpbestfit;
bestpop = X(tmpbestfitidx, :);
stallgenctrl = 0;
t = 1;

% if plots
%     hold on;
% end
%% Main loop
while ((t <= numofruns) && (stallgenctrl <= stallgenlimit))
    % calculating the acceleration for each ball
    A(:, :) = 0;
    
    % Choosing global or local algorithm
    % glocal = 0: Global search, 1: Local search
    
    if glocal
        % local version
        for i = 1:numofballs
            dists = dist(X(i, :), X');
            [~, localind] = sort(dists);
            localind = localind(2:localnum + 1);
            
            for j = 1:localnum
                dheight = heights(localind(j)) - heights(i);
                
                % uses better balls to estimate the slope and calculate the
                % acceleration. In addition, it ensures (X(i, :) - X(j, :)) > 0
                % for all dimensions.
                if dheight < 0
                    A(i, :) = A(i, :) + sin(atan(dheight ./ (X(i, :) - X(localind(j), :))));
%                     A(i, :) = A(i, :) + sin(abs(dheight)./(sqrt(((dheight).^2)+((X(i, :) - X(j, :)).^2))));
                end
            end
        end
    else
        % global version
        for i = 1:numofballs
            for j = 1:numofballs
                dheight = heights(j) - heights(i);
                
                % uses better balls to estimate the slope and calculate the
                % acceleration. In addition, it ensures (X(i, :) - X(j, :)) > 0
                % for all dimensions.
                if dheight < 0
                    A(i, :) = A(i, :) + sin(atan(dheight ./ (X(i, :) - X(j, :))));
%                     A(i, :) = A(i, :) + sin(abs(dheight)./(sqrt(((dheight).^2)+((X(i, :) - X(j, :)).^2))));
                end
            end
        end
    end
    
    % sigmoid method for changing coefficient:
    % higher c1: faster convergence in first steps, worst local search
    % lower c1: slower convergence and better global search in first steps,
    % better local search
    % c2 results to better local search
%     k1 = c1 ./ (1 + exp((t - shift1) .* scale1));
%     k2 = c2 ./ (1 + exp(-(t - shift2) .* scale2));

%%%***********************************************************************
     k1 = k1damp*( (numofruns - t) / numofruns) ;
     k2 = k2damp*( t / numofruns)               ;
 %%%**********************************************************************   
    
    % updating balls
    besttoX = repmat(bestpop, numofballs, 1) - X;
    X = X + k1 .* rand(numofballs, numofdims) .* A + ...
        k2 .* rand(numofballs, numofdims) .* besttoX;

    % ensures that all balls lie in the problem's boundaries
    tmpmaxchk = X > Xmax;
    tmpminchk = X < Xmin;
    X = X .* ~(tmpmaxchk | tmpminchk) + Xmax .* tmpmaxchk + Xmin .* tmpminchk;
    
    % evaluates fitness of each ball
    
    %%%*******************************************************************
    heights = Fitness(X,Hfilt,Wfilt,N);
    %%%**********************************************************************
    
    Neval = Neval + numofballs;

    % finding and storing the global best ball and its fitness
    [tmpbestfit, tmpbestfitidx] = min(heights);
    
    if abs(tmpbestfit - bestfit) < TolFun
        stallgenctrl = stallgenctrl + 1;
    else
        stallgenctrl = 0;
    end
    
    if tmpbestfit < bestfit
        bestfit = tmpbestfit;
        bestpop = X(tmpbestfitidx, :);
    end
      
%     t = t + 1;
    
 % updating variables for displaying the results
    meanfits(t) = mean(heights);
    bests(t) = bestfit;
    worsts(t) = max(heights);
     
     %%
    if plots
        disp(['Iteration ' num2str(t) '  :BestCost= ' num2str(bests(t))]);
%         t
%         plot(t, bests(t), '.r','LineWidth',1);
% %       legend('bests - ipo')
%         xlabel('Iteration')
%         ylabel('Fitness')
%         plot(t, meanfits(t), '.b','LineWidth',4);
%         legend('best','mean')
%         xlabel('Iteration')
%         ylabel('Fitness')
% %         plot(t, worsts(t), '.r');
%         figure(gcf);
%         hold on
    end
    t = t + 1;
end
%% Implementation of Output Filter *******************************

[Bsoa Asoa Z_f P_f] = Matching_Coefs(bestpop)
IIR_main();
disp([ ' Best Solution = '  num2str(bestpop)])
disp([ ' Best Fitness = '  num2str(bests(t-1))])

disp([ ' Time = '  num2str(toc)])

figure(1);
plot(bests,'k','LineWidth',2);
plot(bests,'.b','LineWidth',1);
hold on
plot(meanfits,'.r','LineWidth',1);
legend('best','mean')
xlabel('Iteration')
ylabel('Fitness')
hold off
figure(2);
zplane(Z_f,P_f); %%% Displays the poles and zeros of discrete-time systems.
legend('Zero','Pole');
xlabel('Real Part');
ylabel('Imaginary Plot');
title('Pole-Zero Plot in IPO');
figure(3);
H = abs(Hfilt);
Hdb=20*log10(H);
plot(Wfilt/512,Hdb);grid
title('Magnitude response of a chebyshev I bandpass filter');
hold off

%% End of main loop
if (stallgenctrl > stallgenlimit)
    display('Halted for stallgen');
    display(t);
end
% if plots
%     hold off
% end
% display(bestfit)









