%Copyright (C) 2016 Piotr Beben

%This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
%This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
%You should have received a copy of the GNU General Public License along with this program. If not, see http://www.gnu.org/licenses/.





classdef TestFunctions
    
    
    properties
    end
    
    methods  (Static)
        
        %---------------------------------------------------------------------------------------------
        
        function test_2D(cloud, complex, complexPts, sensitivity, iterations, stepOver, delay, cloudColor, printProgress, varargin)
            
            if isempty(varargin)
                figureType = 'default';
            else
                figureType = varargin{1};
            end;
            

            
            fh = figure(1);
            set(fh,'Units','Normalized','Position',[0.2, 0.1, 0.6, 0.75])
            
            hold off
            
            PointCloud.plot2D(cloud,cloudColor);
            
            hold on
            grid on
            
            [p0,colors] = complex.plot2D(complexPts,[]);
            set(p0,'facealpha',0.35);
            
            switch figureType
                case 'figure_fit_data'
                    m = min(cloud(:));
                    M = max(cloud(:));
                    setAxis = true;
                case 'figure_fit_complex'
                    m = min(complexPts(:));
                    M = max(complexPts(:));
                    setAxis = true;
                otherwise
                    setAxis = false;
            end
            
            if setAxis
                P = 0.1*(M-m);
                vaxis = [m-P, M+P, m-P, M+P];
                axis(vaxis)
            end
            
            set(gca,'FontSize',10);
            pause(1+delay);
            
            
            fittedPts = complexPts;
            optimalOrder = [];
            
            for i = 1:iterations
                
                for j=1:stepOver
                    [fittedPts,optimalOrder,stats]  = SMeans.next_fitting(cloud,complex,fittedPts,sensitivity,printProgress,optimalOrder);
                end
                
                
                %%{
                delete(p0)
                
                [p0,~] = complex.plot2D(fittedPts,colors);
                set(p0,'facealpha',0.35);
                %}
                
                
                
                if setAxis; axis(vaxis); end
                
                title(sprintf('iter #%d, avg. change %0.4f, max. change %0.4f,',i*stepOver,stats(1),stats(2)),'FontSize',14)
                
                
                pause(delay);
                
            end
        end
        %---------------------------------------------------------------------------------------------
        
        
        function test_3D(cloud, complex, complexPts, sensitivity, iterations, stepOver, delay, cloudColor, printProgress, varargin)
            
            if complex.get_dimension() < 3
                bdComplex = complex;
            else
                bdComplex = complex.get_boundary();
            end
            
            if isempty(varargin)
                figureType = 'default';
            else
                figureType = varargin{1};
            end;
            
            
            
            fh = figure(1);
            set(fh,'Units','Normalized','Position',[0.2, 0.1, 0.6, 0.75])
            
            
            hold off
            
            PointCloud.plot3D(cloud,cloudColor);
            
            hold on
            grid on
            
            [p0,colors] = bdComplex.plot3D(complexPts,[]);
            set(p0,'facealpha',0.2);
            
            switch figureType
                case 'figure_fit_data'
                    m = min(cloud(:));
                    M = max(cloud(:));
                    setAxis = true;
                case 'figure_fit_complex'
                    m = min(complexPts(:));
                    M = max(complexPts(:));
                    setAxis = true;
                otherwise
                    setAxis = false;
            end
            
            if setAxis
                P = 0.1*(M-m);
                vaxis = [m-P, M+P, m-P, M+P, m-P, M+P];
                axis(vaxis)
            end
            
            set(gca,'FontSize',10);
            pause(1+delay);
            
            
            fittedPts = complexPts;
            optimalOrder = [];
            
            for i = 1:iterations
                
                for j=1:stepOver
                    [fittedPts,optimalOrder,stats]  = SMeans.next_fitting(cloud,complex,fittedPts,sensitivity,printProgress,optimalOrder);
                end
                
                %%{
                delete(p0)
                
                [p0,~] = bdComplex.plot3D(fittedPts,colors);
                set(p0,'facealpha',0.2);
                %}
                
                
                if setAxis; axis(vaxis); end
                
                title(sprintf('iter #%d, avg. change %0.4f, max. change %0.4f,',i*stepOver,stats(1),stats(2)),'FontSize',14)
                
                pause(delay);
                
                
            end
        end
        
        %---------------------------------------------------------------------------------------------
        
        
        
    end
    
end

