clear all %#ok<CLALL>
close all
clc

Names = {'Issam' 'Sutanu' 'Guangyi' 'Hankun' 'Ana' 'Mingyuan' 'Barron' 'Wei' 'Mohamed' 'Marcelo' 'Orlando' 'Jordan' 'Joshua' 'Da' 'Mengchen' 'Will' 'Yufeng' 'Enguang' 'John' 'GA 2017' 'Simplex'};
Participants = [20 21] ;  

times = zeros(length(Participants),3) ;
fvals = times ;
scores = zeros(length(Participants),2,4) ;
ranks = zeros(length(Participants),2,4) ;

options.maxIter = 10000 ;
options.display = 1 ;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Challenge 1: 2D Morlet Wavelet %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('Challenge 1: 2D Morlet Wavelet\n')

func = @(x) TwoD_MorletProblem(x) ;
bounds = [-10 10; -10 10] ; 

iters = 10000*ones(length(Names)) ;
iters(21) = 500 ; 
iters(19) = 50 ;

fprintf('Commencing simulations\n')
for cnt = 1:length(Participants)
  options.maxIter = iters(Participants(cnt)) ;
  if Participants(cnt) < 10
    fstring = ['max_00' num2str(Participants(cnt))] ;
  else
    fstring = ['max_0' num2str(Participants(cnt))] ;
  end
  fprintf(['' fstring ':\n'])
  %Assess computational time and accuracy
  tic
  eval(['[x_opt,x_val] = ' fstring '(func,bounds,options);'])
  times(cnt,1) = toc ;
  for cntr = 1:length(bounds)
    if x_opt(cntr) < bounds(cntr,1) 
      x_opt(cntr) = bounds(cntr,1) ;
    elseif x_opt(cntr) > bounds(cntr,2)
      x_opt(cntr) = bounds(cntr,2) ;
    end
  end
  fvals(cnt,1) = func(x_opt) ;
  if fvals(cnt,1) == Inf || isnan(fvals(cnt,1))
    fvals(cnt,1) = 10 ;
  end
  fprintf(['Time for ' fstring ' is ' num2str(times(cnt,1)) ' seconds.\n'])
  fprintf(['Maximum for ' fstring ' is ' num2str(fvals(cnt,1)) '.\n\n'])
end
