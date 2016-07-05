%Copyright (C) 2016 Piotr Beben

%This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
%This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
%You should have received a copy of the GNU General Public License along with this program. If not, see http://www.gnu.org/licenses/.



classdef PointCloud
    

    
    methods (Static)

        %--------------------------------------------------------------------------------------------- 
        function output = precision(newPrecision)
            
            persistent currentPrecision;
                 
            if nargin >= 1
                currentPrecision = newPrecision;
            elseif isempty(currentPrecision)
                currentPrecision = 'double';
            end
            
            output = currentPrecision;
        end
        
        %---------------------------------------------------------------------------------------------
        
        function points = build_grid( subdivisions, dims, center )
            
            dimension = length(subdivisions);
            
            if dimension == 0
                points = center;
            else
                d = dimension;
                section = PointCloud.build_grid( subdivisions(1:(d-1)), dims(1:(d-1)), center );
                
                if subdivisions(d)>0
                    incr = [ zeros(d-1,1); dims(d)/subdivisions(d); zeros(length(center)-d,1) ];
                    bottom = [ zeros(d-1,1); -dims(d)/2; zeros(length(center)-d,1) ];

                    
                    
                    points = NaN(length(center),prod(subdivisions+1),PointCloud.precision());
                    
                    for i = 1:(subdivisions(d)+1)
                        j = (i-1)*size(section,2);
                        points( :, (j+1):(j+size(section,2)) ) = bsxfun(@plus,section,bottom+(i-1)*incr);
                    end
                else
                    points = section;
                end
            end
            
        end
        
   
        
        %---------------------------------------------------------------------------------------------
        
        function points = build_random( dimension, numPoints, type, mode, varargin )
            
            
            switch type
                
                case 'centered'
                    
                    center = varargin{1};           %column vector
                    axisScale = varargin{2};
                    
                    switch mode
                        case 'normal'
                            points = bsxfun( @plus, diag(axisScale)*randn(dimension,numPoints), center );
                        case 'uniform'
                            points = bsxfun( @plus, diag(2*axisScale)*(rand(dimension,numPoints)-1), center );
                        otherwise
                            
                            throw( MException('Invalid:Option','Invalid Option: %s',mode) );
                    end
                    
                case 'custom'
                    
                    curve_array = varargin{1};      %cell array of parametric curves (function handles), 1 per dimension                    
                    tMin = varargin{2};             %parameter lower and upper bounds, 1 entry for each unique parameter
                    tMax = varargin{3};
                    noise = varargin{4};
                    
                    
                    switch mode
                        case 'normal'
                            rand_noise = @() noise*randn(1,numPoints,PointCloud.precision());
                        case 'uniform'
                            rand_noise = @() (2*noise)*rand(1,numPoints,PointCloud.precision()) - noise;
                        case 'none' 
                            rand_noise = @() zeros(1,numPoints,PointCloud.precision());
                        otherwise
                            throw( MException('Invalid:Option','Invalid Option: %s',mode) );
                    end
                    
       
                    
                    
                    parameters = cell( 1, length(tMin) );
                    
                    for i = 1:length(parameters)
                        parameters{i} = (tMax(i)-tMin(i))*rand(1,numPoints) + tMin(i);
                    end
                    
                                        
                    points = NaN(dimension,numPoints,PointCloud.precision());
                    
                    for i = 1:dimension
                        f = curve_array{i};
                        points(i,:) = f( parameters ) + rand_noise();
                       
                    end
                    
        
                    
                otherwise
                    throw( MException('Invalid:Option','Invalid Option: %s',type) );
                    
            end
            
        end
        
        %---------------------------------------------------------------------------------------------
        function [dims,center] = get_cloud_dimensions(cloud)
            
            dims = zeros(1,size(cloud,1)); 
            center = zeros(size(cloud,1),1);
            
            for i = 1:length(dims)
                M = max(cloud(i,:));
                m = min(cloud(i,:));
                dims(i) = M - m;
                center(i) = m + (M-m)/2;
            end
        
        end
        
        
        %---------------------------------------------------------------------------------------------
        
        
        function plot = plot2D( points, color )
            
            plot = scatter( points(1,:), points(2,:), 20, color, 'filled' );

        end
        %---------------------------------------------------------------------------------------------
        
        
        function plot = plot3D( points, color )
            
            plot = scatter3( points(1,:), points(2,:), points(3,:), 10, color, 'filled' );

        end
        
        %---------------------------------------------------------------------------------------------
        
    end
    
    
   
end

