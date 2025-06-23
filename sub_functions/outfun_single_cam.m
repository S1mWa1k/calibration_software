function stop=outfun_single_cam(x,optimValues,state)

% OPTIMPLOTRESNORM Plot value of the norm of residuals at each iteration.
%
%   STOP = OPTIMPLOTRESNORM(X,OPTIMVALUES,STATE) plots OPTIMVALUES.resnorm.
%
%   Example:
%   Create an options structure that will use OPTIMPLOTRESNORM as the plot
%   function
%     options = optimset('PlotFcns',@optimplotresnorm);
%
%   Pass the options into an optimization problem to view the plot
%     lsqnonlin(@(x) sin(3*x),[1 4],[],[],options);

%   Copyright 2006 The MathWorks, Inc.
%   $Revision: 1.1.6.1 $  $Date: 2006/11/11 22:46:32 $
% 
% res=resnorm./size(x,2);

persistent plotavailable plotresnorm
stop = false;

switch state
     case 'init'
         if isfield(optimValues,'resnorm')
             plotavailable = true;
         else
             plotavailable = false;
%              title('Norm of Residuals: not available','interp','none');
         end
    case 'iter'
        if plotavailable
            if optimValues.iteration == 0
                figure
                drawnow
%                   plotresnorm=waitbar(0,'Working...');
                  
%                 % The 'iter' case is  called during the zeroth iteration,
%                 % but it now has values that were empty during the 'init' case
%                 plotresnorm = area(optimValues.iteration,optimValues.resnorm/length(optimValues.residual)); hold on
                plotresnorm = plot(optimValues.iteration,optimValues.resnorm/length(optimValues.residual),'linewidth',2); hold on
                
                xlabel('Iteration','interp','none');
                set(plotresnorm,'Tag','optimplotresnorm');
                ylabel('Mean pixel error','interp','none');
                xlim([0 150])
                set(gca,'YScale','log')
                if optimValues.resnorm/length(optimValues.residual)<1.1
                    ylim([0 1.1])
                else
                    ylim([0 optimValues.resnorm/length(optimValues.residual)])
                end
                
                plot([xlim],[1 1],'r','linewidth',2)
                title(['Norm of Residuals: ',num2str(optimValues.resnorm)],'interp','none');

            else
%                   waitbar(optimValues.iteration/20);
                  
                plotresnorm = findobj(get(gca,'Children'),'Tag','optimplotresnorm');
                newX = [get(plotresnorm,'Xdata') optimValues.iteration];
                newY = [get(plotresnorm,'Ydata') optimValues.resnorm/length(optimValues.residual)];
                set(plotresnorm,'Xdata',newX, 'Ydata',newY);
                xlim([0 150])
                set(get(gca,'Title'),'String',['Norm of Residuals: ',num2str(optimValues.resnorm)]);
            end
        end
        
    case 'done'
%         stop=true;
        base_gca=get(plotresnorm,'parent');
        p=get(base_gca,'parent');
        close(p)
end
drawnow