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
                    m1 = min(cloud(1,:));
                    M1 = max(cloud(1,:));
                    m2 = min(cloud(2,:));
                    M2 = max(cloud(2,:));
                    setAxis = true;
                case 'figure_fit_complex'
                    m1 = min(complexPts(1,:));
                    M1 = max(complexPts(1,:));
                    m2 = min(complexPts(2,:));
                    M2 = max(complexPts(2,:));
                    setAxis = true;
                 case 'figure_fit_all'
                    m1 = min( min(complexPts(1,:)), min(cloud(1,:)) );
                    M1 = max( max(complexPts(1,:)), max(cloud(1,:)) );
                    m2 = min( min(complexPts(2,:)), min(cloud(2,:)) );
                    M2 = max( max(complexPts(2,:)), max(cloud(2,:)) );
                    setAxis = true;    
                otherwise
                    setAxis = false;
            end
            
            if setAxis
                P1 = 0.1*(M1-m1);
                P2 = 0.1*(M2-m2);
                if P1==0; P1 = 0.01; end;
                if P2==0; P2 = 0.01; end;
               
                vaxis = [m1-P1, M1+P1, m2-P2, M2+P2];
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
                    m1 = min(cloud(1,:));
                    M1 = max(cloud(1,:));
                    m2 = min(cloud(2,:));
                    M2 = max(cloud(2,:));
                    m3 = min(cloud(3,:));
                    M3 = max(cloud(3,:));
                    setAxis = true;
                case 'figure_fit_complex'
                    m1 = min(complexPts(1,:));
                    M1 = max(complexPts(1,:));
                    m2 = min(complexPts(2,:));
                    M2 = max(complexPts(2,:));
                    m3 = min(complexPts(3,:));
                    M3 = max(complexPts(3,:));
                    setAxis = true;
                case 'figure_fit_all'
                    m1 = min( min(complexPts(1,:)), min(cloud(1,:)) );
                    M1 = max( max(complexPts(1,:)), max(cloud(1,:)) );
                    m2 = min( min(complexPts(2,:)), min(cloud(2,:)) );
                    M2 = max( max(complexPts(2,:)), max(cloud(2,:)) );
                    m3 = min( min(complexPts(3,:)), min(cloud(3,:)) );
                    M3 = max( max(complexPts(3,:)), max(cloud(3,:)) );
                    setAxis = true;
                otherwise
                    setAxis = false;
            end
            
            if setAxis
                P1 = 0.1*(M1-m1);
                P2 = 0.1*(M2-m2);
                P3 = 0.1*(M3-m3);  
                if P1==0; P1 = 0.01; end;
                if P2==0; P2 = 0.01; end;
                if P3==0; P3 = 0.01; end;

                vaxis = [m1-P1, M1+P1, m2-P2, M2+P2, m3-P3, M3+P3];
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

