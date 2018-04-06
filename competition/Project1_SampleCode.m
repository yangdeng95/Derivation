clear all %#ok<CLALL>
% close all
clc

%%%%%%%%%%%%%%%%%
%% Challenge 1 %%
%%%%%%%%%%%%%%%%%

func1 = @(x) 3*sin(2*x)+1e6 ;
xval = 1 ;
deriv = -12*sin(2*x) ;
NTimes = 1e6 ;

year = '1' ;
% Names = {'Senyo' 'Nidish' 'Andrew' 'Nick' 'Iyabo' 'Qi' 'Mianmian' 'Scott' 'John' 'Deng' 'Kaichun' 'Sichao'} ;
% Participants = [1 2 3 4 5 6 7 8 9 10 11 12] ;

% Modify this accordingly to test out your different algorithms:
Names = {'Algorithm 1' 'Algorithm 2'} ;
Participants = [0] ;

times1 = zeros(length(Participants),1) ;
errors1 = zeros(length(Participants),1) ;
scores1 = zeros(length(Participants),2) ;
ranks1 = zeros(length(Participants),2) ;

fprintf(['Challenge 1: ' func2str(func1) '\n'])
for cnt = 1:length(Participants)
  if Participants(cnt) < 10
    fstring = ['func_' year '0' num2str(Participants(cnt))] ;
  else
    fstring = ['func_' year num2str(Participants(cnt))] ;
  end
  fprintf(['' fstring ':\n'])
  eval(['fnc = @(func,x) ' fstring '(func,x);'])
  %Assess computational time, averaged over NTimes solutions:
  tic
  for cntr = 1:NTimes
    d2f = fnc(func1,xval) ;
  end
  times1(cnt) = toc/NTimes ;
  errors1(cnt) = abs(d2f - deriv)/abs(deriv) ;
  fprintf(['Average time for ' fstring ' is ' num2str(times1(cnt)) ' seconds.\n'])
  fprintf(['Error for ' fstring ' is ' num2str(errors1(cnt)*100) '%%.\n\n'])
end

for cnt = 1:length(Participants)
  scores1(cnt,1) = -log10((errors1(cnt)/(max(errors1)) ));
end
scores1(:,1) = scores1(:,1)/max(scores1(:,1)) ;
for cnt = 1:length(Participants)
  scores1(cnt,2) = 0.5*scores1(cnt,1) + 0.5*(1-times1(cnt)/max(times1)) ;
end

ranks1(:,1) = rankdata(scores1(:,1)) ;
ranks1(:,2) = rankdata(scores1(:,2)) ;

figure()
imagesc(ranks1)
title(['Challenge 1: ' func2str(func1)])
set(gca,'XTick',[1 2],'XTickLabel',{'Accuracy' 'Weighted'})
set(gca,'YTick',linspace(1,length(Participants),length(Participants)))
set(gca,'YTickLabel',{Names{1,Participants}}) %#ok<CCAT1>
colormap('Hot')
set(gcf,'Position',[5 100 490 675])
for cntr = 1:length(Participants)
  text(1,cntr,[num2str(scores1(cntr,1)*100) '%'],'FontWeight','bold','HorizontalAlignment','center')
  text(2,cntr,[num2str(scores1(cntr,2)*100) '%'],'FontWeight','bold','HorizontalAlignment','center')
end

