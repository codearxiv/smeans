%Copyright (C) 2016 Piotr Beben

%This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
%This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
%You should have received a copy of the GNU General Public License along with this program. If not, see http://www.gnu.org/licenses/.




classdef SMeans
   
  
    methods (Static)
        
        %---------------------------------------------------------------------------------------------
        
        function points = baryc_to_euclid( baryCoords, simplexPts )
            
            %points =  bsxfun( @plus, simplex_pts(:,2:end)*baryCoords, simplex_pts(:,1) );
            points = simplexPts*baryCoords;
            
        end
        
        %---------------------------------------------------------------------------------------------       
        function distancesSq = distSq_pt_pts( point, points )
            
            distancesSq = sum( bsxfun(@minus,points,point).^2 , 1 );
            
        end

        %---------------------------------------------------------------------------------------------
        
        function distancesSq = distSq_pts_pts( points1, points2 )
            
            distancesSq = sum( (points2-points1).^2 , 1 );
            
        end
    
        %---------------------------------------------------------------------------------------------
        
        function estimate = size_estimate( simplexPts )
         
            centroid = sum(simplexPts,2)./size(simplexPts,2);            
            estimate = sqrt(sum( SMeans.distSq_pt_pts(centroid,simplexPts) ));
            
        end
        
       
        %---------------------------------------------------------------------------------------------
        %get indices of those points for which exists a nearest point on *interior* of simplex
        
        function [ baryCoords, indexedInterior ] = nearest_points_on_simplex_plane( points, simplexPts )
            
            
            if size(simplexPts,2) >= 2
                
                simplexVect = simplexPts(:,2:end)-repmat(simplexPts(:,1),1,size(simplexPts,2)-1);
                %simplexVect = bsxfun( @minus, simplexPts(:,2:end), simplexPts(:,1) );
                pointsVect = bsxfun( @minus, points, simplexPts(:,1) );
                
                simplexPInv = pinv(simplexVect);
                baryCoordsNorm = simplexPInv * pointsVect;                
                baryNormSum = sum(baryCoordsNorm,1);
                baryCoords = [ 1-baryNormSum; baryCoordsNorm ];
                
                indexedInterior = find( all( baryCoordsNorm >= 0, 1 ) & baryNormSum <= 1 );
                
                
                
            else
                indexedInterior = find( all( bsxfun(@eq,points,simplexPts(:,1)), 1) );
                baryCoords = ones(1,length(indexedInterior));
            end
                
         
        end
        
        %---------------------------------------------------------------------------------------------
        
        function [ baryCoords, indices ] = points_inside_simplex( points, simplexPts )
            
            
            if size(simplexPts,2) >= 2
                
                simplexVect = simplexPts(:,2:end)-repmat(simplexPts(:,1),1,size(simplexPts,2)-1);
                %simplexVect = bsxfun( @minus, simplexPts(:,2:end), simplexPts(:,1) );
                pointsVect = bsxfun( @minus, points, simplexPts(:,1) );
                
                simplexPInv = pinv(simplexVect);
                baryCoordsNorm = simplexPInv * pointsVect;
                baryNormSum = sum(baryCoordsNorm,1);
                

                indices = find( all( baryCoordsNorm >= 0, 1) & baryNormSum <= 1 & SMeans.distSq_pts_pts(pointsVect,simplexVect*baryCoordsNorm) == 0 );
                
                baryCoords = [ 1-baryNormSum(indices); baryCoordsNorm(:,indices) ];
            
            else        
                indices = find( SMeans.distSq_pt_pts(simplexPts(:,1),points) == 0 );              
                baryCoords = zeros(1,length(indices),'like',points);
            end
            
        end
        %---------------------------------------------------------------------------------------------
        
        function [ baryCoords, distances ] = distSq_pts_to_simplex( points, simplexPts )
 
          
            if size(simplexPts,2) >= 2
                
                distances = inf(1,size(points,2),'like',points);
                
                [ baryCoords, indexedInterior ] = SMeans.nearest_points_on_simplex_plane( points, simplexPts );
         
                if ~isempty(indexedInterior)
                                                        
                    distances(indexedInterior) = ...
                        SMeans.distSq_pts_pts(points(:,indexedInterior),SMeans.baryc_to_euclid(baryCoords(:,indexedInterior),simplexPts));
                 
                end
                
                
                
                
                if length(indexedInterior)<size(points,2)
                    
                    leftOverIndices = find(distances==inf);
                    
                    for i = 1:size(simplexPts,2)
                        
                        markedHalfPlane = baryCoords(i,leftOverIndices) < 0;
                        indexedHalfPlane = leftOverIndices( markedHalfPlane );
                        
                        if ~isempty(indexedHalfPlane)
                            [ bdBary, bdDist ] = ...
                                SMeans.distSq_pts_to_simplex(points(:,indexedHalfPlane),simplexPts(:,[1:(i-1),(i+1):end]));
                            
                            isSmaller = bdDist < distances(indexedHalfPlane);
                            indexedSmaller = indexedHalfPlane(isSmaller);
                            distances( indexedSmaller ) = bdDist(isSmaller);
                            baryCoords( [1:(i-1),(i+1):end], indexedSmaller ) = bdBary(:,isSmaller);
                            baryCoords( i, indexedSmaller ) = 0;
                            
                            leftOverIndices = leftOverIndices(~markedHalfPlane);
                            
                        end
                        
                    end
                    
                    
                    
                end
                

            else
                  
                distances = SMeans.distSq_pt_pts(simplexPts(:,1),points);
                baryCoords = ones(1,size(points,2),'like',points);
                         
            end
            
        end
        
   
        %---------------------------------------------------------------------------------------------
        
        function [ baryAssign, facetAssign ] = ...
                partition_nearest( cloud, complex, complexPts, printProgress, varargin)
            
                                
            
            if printProgress; progress = 0; end;
       
            
            baryAssign = cell(1,size(cloud,2));
            facetAssign = zeros(1,size(cloud,2),'uint32');
            
            leftOverIndices = 1:size(cloud,2);
            

            
            
            if ~isempty(leftOverIndices)
                
                distances = inf(1,size(leftOverIndices,2),'like',cloud);
                extCloud = cloud(:,leftOverIndices); 
                
                for i = 1:(length(complex.facets))
                    
                    
                    simplexPts = complexPts( :, complex.facets{i} );
                    
                    [ baryCoords, distCurrent ] = SMeans.distSq_pts_to_simplex(extCloud,simplexPts);
                    
                    isSmaller = distCurrent < distances;
                    distances(isSmaller) = distCurrent(isSmaller);
                    baryAssign( leftOverIndices(isSmaller) ) = num2cell( baryCoords(:,isSmaller), 1 );
                    facetAssign( leftOverIndices(isSmaller) ) = i;
                    
                    
                    %--------------
                    %Print progress
                    if printProgress
                        currentProgress = 100*i/length(complex.facets);
                        if (currentProgress - progress) >=5
                            sprintf('1st part of current fitting is %0.0f %% complete',currentProgress)
                            progress = currentProgress;
                        end
                    end
                    %--------------         
                end
                             
            end
            
            
            
        end
        
        
        %---------------------------------------------------------------------------------------------
        %closed == true uses set \bar{N}^\ell_j = {y\in\mc S | y'\in w for some w s.t. v_j\in w}
        
        function [ complexPts, varargout ] = ...
                next_fitting( cloud, complex, complexPts, sensitivity, closed, printProgress, varargin )
            
            
            if printProgress; progress = 0; end;
            
            s = sensitivity;
            
            nout = max(nargout,1) - 1;
            getStats = (nout>0);
            dist = 0;
            
            
            if closed
                if isempty(varargin)
                    vertexStars = complex.get_vertex_stars(1:size(complexPts,2));
                else
                    vertexStars = varargin{1};
                end
            end
            
            
            
            
            [ baryAssign, facetAssign ] = ...
                SMeans.partition_nearest( cloud, complex, complexPts, printProgress );
            
            markedPts = false(1,size(complexPts,2));
            fittedPts = zeros(size(complexPts,1),size(complexPts,2),'like',complexPts);
            weightCounts = zeros(1,size(complexPts,2));
            
            
            
            for i = 1:size(cloud,2)
                
                nearestFacet = complex.facets{facetAssign(i)};
                markedNonZero = baryAssign{i}>0.0;
                bdFace = nearestFacet( markedNonZero );
                bdBaryAssign = baryAssign{i}( markedNonZero );
                
                
                
                if closed && length(bdFace)<length(nearestFacet)
                    
                    star = vertexStars{bdFace(1)};
                    
                    for j = 2:length(bdFace)
                        star = intersect( star, vertexStars{bdFace(j)} );
                    end
                    
                    starVertices = complex.facets{star(1)};
                    
                    for j = 2:length(star)
                        starVertices = union( starVertices, complex.facets{star(j)} );
                    end
                    
                    [~,~,bdIndices] = intersect(bdFace,starVertices);
                    baryAssignStar = zeros(length(starVertices),1);
                    baryAssignStar( bdIndices ) = bdBaryAssign;
                    
                else
                    starVertices = bdFace;
                    baryAssignStar = bdBaryAssign;
                end
                
                
                
                unmarkedIndices = find( markedPts(starVertices)==false );
                if ~isempty(unmarkedIndices)
                    markedPts( starVertices(unmarkedIndices) ) = true;
                end
                
                
                fittedPts(:,starVertices) = fittedPts(:,starVertices) + ...
                    complexPts(:,starVertices)*diag((1-baryAssignStar)./(1+s)) + cloud(:,i)*((baryAssignStar+s)./(1+s))';
                
                weightCounts(starVertices) = weightCounts(starVertices) + 1;
                          
                %--------------
                %Get stats
                if getStats
                    dist = dist + ...
                        SMeans.distSq_pt_pts( cloud(:,i), SMeans.baryc_to_euclid(bdBaryAssign,complexPts(:,bdFace)) );
                end
                %--------------
                 
                
                %--------------
                %Print progress
                if printProgress
                    currentProgress = 100*i/size(cloud,2);
                    if (currentProgress - progress) >=5
                        sprintf('2nd part of current fitting is %0.0f %% complete',currentProgress)
                        progress = currentProgress;
                    end
                end
                %--------------
                
            end
      
            
            fittedPts(:,markedPts) = bsxfun( @(x,y) x./y , fittedPts(:,markedPts), weightCounts(markedPts) );

            %-------------- 
            %Get stats
            if getStats

                indexedPts = find(markedPts);
                
                if ~isempty(indexedPts)
                    diff = sqrt(SMeans.distSq_pts_pts(fittedPts(:,indexedPts),complexPts(:,indexedPts)));
                else
                    diff = 0;
                end
             
                varargout{1} = [mean(diff),max(diff),dist];  
                    
            end       
            %--------------
  
            
                
            complexPts(:,markedPts) = fittedPts(:,markedPts);
            
     
            
        end
        
        %---------------------------------------------------------------------------------------------

        
        
        
        
        
        
        
        
        
        
        
        %---------------------------------------------------------------------------------------------
        
        function fittedPts = fit( cloud, complex, complexPts, sensitivity, closed, printProgress, stopConditions )
            
            
            
            fittedPts = complexPts;

            if closed
                vertexStars = complex.get_vertex_stars(1:size(complexPts,2));
            else 
                vertexStars = [];
            end
            
            
            switch stopConditions{1}
                
                case 'stop_on_iterations'

                    iterations = stopConditions{2};


                    
                    for i = 1:iterations
                        [fittedPts] = ...
                            SMeans.next_fitting(cloud,complex,fittedPts,sensitivity,closed,printProgress,vertexStars);

  
                        
                        if printProgress
                            sprintf('iter #%d complete,',i)
                        end
                        
                    end
                    
                    
                case {'stop_on_avg_change','stop_on_max_change'}
                    
                    
                    switch stopConditions{1}
                        case 'stop_on_avg_change'; statsSelect=1;
                        case 'stop_on_max_change'; statsSelect=2;
                    end
                    
                    change = stopConditions{2};
                    maxSuccess = stopConditions{3};
                    
                    currentSuccess = 0;
                    i= 0;
                    
                    while currentSuccess < maxSuccess 
                        [fittedPts,stats]  = ... 
                            SMeans.next_fitting(cloud,complex,fittedPts,sensitivity,closed,printProgress,vertexStars);
                        
                        
                        if printProgress
                            sprintf('iter #%d complete, avg. change %0.4f, max. change %0.4f, sum sq. dist. %0.4f',i,stats(1),stats(2),stats(3))
                        end
                        
                        if stats(statsSelect) < change  
                            currentSuccess = currentSuccess + 1;
                        else
                            currentSuccess = 0;
                        end
                        
                        i = i + 1;
                        
                    end
                    

                    
                otherwise
                    throw( MException('Invalid:Option','Invalid Option: %s',stopConditions{1}) );
            
            
            end
            

        end

        
        %---------------------------------------------------------------------------------------------
        

        
       

        
    end
end

