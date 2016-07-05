classdef FacetComplex < handle
    %A representation of a simplicial complex (and operations on it) 
    %in terms of its top dimensional facets. To save memory. 

 
    
    
    
    properties (SetAccess = protected)
        facets 
    end
    
    
    
    
    
    methods (Static)
        
        %--------------------------------------------
        %for sorted faces
        function answer = iscontained_sorted( face1, face2 )
            answer = length(face1)<=length(face2) && all(ismembc(face1,face2));
        end
        %--------------------------------------------
        %for unsorted faces
        function answer = iscontained( face1, face2 )
            answer = length(face1)<=length(face2) && length(face1)==length(intersect(face1,face2)); 
        end
        %--------------------------------------------
        %for unsorted faces
        function answer = isequiv( face1, face2 )
            answer = length(face1)==length(face2) && length(face1)==length(intersect(face1,face2));
        end
        %--------------------------------------------
        
        function bdFaces = face_boundary( face )
            
            bdFaces = cell(1,length(face));
            for i = 1:length(face)
                bdFaces{i} = face([1:(i-1),(i+1):end]);
            end
            
        end
        %--------------------------------------------
        
        function joinFace = face_join( face1, face2, disjoint )
            
            if disjoint
                joinFace = [ face1; face2 ];
            else
                joinFace = union(face1,face2);
            end
            
        end
        %--------------------------------------------
        
        function prodFaces = face_edge_prod( face, segments, vertexShift )

            
            prodFaces = cell(1,length(face)*segments);         
            sectionFace = cell(1,segments+1);
            
            sectionFace{1} = face;    
            for j = 2:(segments+1)
                sectionFace{j} = face + (j-1)*vertexShift;
            end
            
            for j = 1:segments
                k = (j-1)*length(face);
                
                for i = 1:length(face)
                    prodFaces{k+i} = [ sectionFace{j}(1:i); sectionFace{j+1}(i:end) ];
                end
                
            end
            
            
        end
        %--------------------------------------------
        
        function [ subdivFaces, newPts ] = face_subdiv( face, newVertex, depth, complexPts )
            
            subdivFaces = {};
            newPts = {};
            
            if ~isempty(complexPts)
                newPts{1} = sum( complexPts, 2 )./length(face);
            end
            
            
            if depth>1
                
                numSubdivFaces = length(face)^(depth-1);
                curVertex = newVertex+1;
                
                for i = 1:length(face)
                   
                    newFace = face;
                    newFace(i) = newVertex;
                    newFacePts = [ complexPts(:,1:(i-1)), newPts{1}, complexPts(:,(i+1):end) ];
                    
                    [ subdivFaces( (end+1):(end+numSubdivFaces) ), subdivPts ] = ...
                        FacetComplex.face_subdiv( newFace, curVertex, depth-1, newFacePts );
                    
                    curVertex = curVertex + length(subdivPts);
                    
                    if ~isempty(subdivPts) 
                        newPts( (end+1):(end+length(subdivPts)) ) = subdivPts; 
                    end
                end
            else
                subdivFaces = cell(1,length(face));
                for i = 1:length(face)
                    subdivFaces{i} = face; 
                    subdivFaces{i}(i) = newVertex;
                    
                end
            end
            

            
            
            
        end
        %--------------------------------------------

        function newObj = intersection( obj1, obj2, varargin )
                     
            [~,indices1] = obj2.contains_faces(obj1.facets,'get_contained_indices',varargin);
            [~,indices2] = obj1.contains_faces(obj2.facets,'get_contained_indices',varargin);
            
            newObj = FacetComplex('empty');
            newObj.assign(obj1,'subcomplex',indices1);
            newObj.insert_facets( obj2.facets(indices2), true, varargin );
            
            
        end
        
        %--------------------------------------------
        %use varargin{1}=='sorted' if facet indices in obj1 and obj2 are sorted
        
        function newObj = union( obj1, obj2, varargin )

            newObj = FacetComplex('empty');
            newObj.assign(obj1,'all');
            newObj.insert_facets( obj2.facets, true, varargin );
        end
        %--------------------------------------------
        function newObj = disjoint_union( obj1, obj2, varargin )
            
            
            firstVertex = obj1.get_max_vertices();
            vertexMap = @(x) x+firstVertex;
            
            newObj = FacetComplex('empty');
            newObj.assign(obj1,'all');
            newObj.insert_facets( obj2.get_vertex_mapped(vertexMap).facets, false );
            
        end
        %--------------------------------------------
        function newObj = join( obj1, obj2, varargin )
            
   
            
            vertexShift = obj1.get_max_vertices();
            vertexMap = @(x) x+vertexShift;
            facets2 = obj2.get_vertex_mapped(vertexMap).facets;
            
            newObj = FacetComplex('empty');
            newObj.facets = cell(1,length(obj1.facets)*length(facets2));
            
            for i=1:length(obj1.facets)
   
                k=(i-1)*length(facets2);
                
                for j=1:length(facets2)
                    newObj.facets{k+j} = [obj1.facets{i};facets2{j}];
                end
            end
            
        end
        
        %--------------------------------------------
        function newObj = prod_edge( obj, segments )
            
            
            
            vertexShift = obj.get_max_vertices();
            
            newObj = FacetComplex('empty');
            f = @(x) segments*length(x);
            
            numNewFaces = sum( cellfun(f,obj.facets) );
            newObj.facets = cell(1,numNewFaces);
            
            j=0;
            for i=1:length(obj.facets)
                k = segments*length(obj.facets{i});
                newObj.facets( (j+1):(j+k) ) = FacetComplex.face_edge_prod(obj.facets{i},segments,vertexShift);
                j = j+k;
          
            end
            
        end
        
        %--------------------------------------------
        
        function newObj = difference( obj1, obj2, varargin )
            
            
            [~,containing] = obj1.contains_faces(obj2.facets,'get_containing_indices',varargin);
            [~,contained] = obj1.contains_faces(obj2.facets,'get_contained_indices',varargin);
            
            newObj = FacetComplex('empty');
            newObj.assign( obj1, 'subcomplex', setdiff( 1:length(obj1.facets), [containing,contained] ) );
            
        end
        
        %--------------------------------------------
        

        
        
        
        
        
        %--------------------------------------------
        function newObj = build_simplex( dimension )
            
            newObj = FacetComplex('custom',{ (1:(dimension+1))' },false);
            
        end
        
       
        %--------------------------------------------
        
        function newObj = build_grid( subdivisions )
            
            newObj = FacetComplex('custom',{1},false);
            
            for i = 1:length(subdivisions)
                if subdivisions(i)>0
                    newObj = FacetComplex.prod_edge(newObj,subdivisions(i)); 
                end
            end
            
       
            
        end
        %--------------------------------------------


    end
    
    
    
    
    
    
    
    
    
    methods
        
        %--------------------------------------------
        function obj = FacetComplex( type, varargin )
            
            switch type
                case 'simplex'
                    obj = FacetComplex.build_simplex( varargin{1} );
                case 'points'
                    obj.facets = num2cell(1:varargin{1});
                case 'grid'
                    obj = FacetComplex.build_grid( varargin{1} );     
                case 'assign'
                    obj.assign( varargin{1}, 'all' );
                case 'custom'
                    if length(varargin)<2; varargin{2} = false; end;
                    obj.assign_facets( varargin{1}, varargin{2} );      %(faces, checkRedundancy)
                case 'random'
                case 'ghost'
                    obj.facets = cell( 1, varargin{1} );
                case 'empty'
                    obj.facets = {};
                otherwise
                    throw( MException('Invalid:Option','Invalid Option: %s',type) );
                    
            end
        end
        %--------------------------------------------
        
        function trim_memory( obj )
            isnotempty = @(x) ~isempty(x);
            obj.facets = obj.facets( cellfun(isnotempty,obj.facets) );
        end
        
        %--------------------------------------------
        
        function clear_redundancy( obj, mergeSubFaces, trimMem )
            
            if trimMem; obj.trim_memory(); end
            
            if mergeSubFaces          %Warning: Slow.
                
                obj.sort_facet_indices();         
                oldFacets = obj.facets;
                obj.insert_facets(oldFacets,true,'sorted');
            else
               
                dimension = obj.get_dimension();
                vertexCounts = obj.get_facet_vertex_counts();
                
                newFacets = cell(1,length(obj.facets));
                numNewFacets = 0;
         
                for i = 1:(dimension+1)
                
                    reduced = unique( cell2mat( obj.facets(vertexCounts==i) )' , 'rows' )';
                    
                    if ~isempty(reduced)
                        
                        newFacets( (numNewFacets+1):(numNewFacets+size(reduced,2)) ) = num2cell(reduced,1);                  
                        numNewFacets = numNewFacets + size(reduced,2);                   
                    end
                    
                end
          
                obj.facets = newFacets(1:numNewFacets);
               
            end
        end
        
        %--------------------------------------------
        
        function sort_facet_indices( obj )
            
            obj.facets = cellfun(@sort,obj.facets,'UniformOutput',false);
        end
        %--------------------------------------------
        
        function outFacets = get_facets(obj, varargin)
       
            if isempty(varargin)
                outFacets = obj.facets;
            else
                dimension = varargin{1};    
                outFacets = obj.facets( obj.get_facet_vertex_counts()==dimension+1 );
            end
            
            
        end

        %--------------------------------------------
        
        function outFacet = get_facet(obj, index)
            
            outFacet = obj.facets{index};
            
        end
        %--------------------------------------------
        
        function num = get_num_facets( obj )
            num = length(obj.facets);
        end
        
        %--------------------------------------------
        
        function maxVertices = get_max_vertices( obj )
            maxVertices = max( cellfun(@max,obj.facets) );
        end
        
        %--------------------------------------------
        
        function numVertices = get_num_vertices( obj )          
            numVertices = length( unique( cell2mat( obj.facets ) ) );
        end
        %--------------------------------------------
        
        function dimension = get_dimension( obj )
            dimension = max( cellfun( @length, obj.facets ) ) - 1;
        end
        
        %--------------------------------------------
        
        function counts = get_facet_vertex_counts( obj )
            counts = cellfun(@length,obj.facets);
        end
        %--------------------------------------------
        
        function mappedObj = get_vertex_mapped( obj, vertexMap )
            
            if isnumeric(vertexMap)
                map = @(i) vertexMap(i);
            else
                map = vertexMap;
            end
            
            mappedObj = FacetComplex('empty');
            mappedObj.facets = cellfun( map, obj.facets, 'UniformOutput', false );
           
            
        end
        %--------------------------------------------
        
        function newObj = get_boundary( obj )
            
            newObj = FacetComplex('empty');
            
            for i = 1:length(obj.facets)
                
                if length(obj.facets{i}) >= 2                    
                    newObj.insert_facets( FacetComplex.face_boundary( obj.facets{i} ), false );
                end
            end
            
            newObj.clear_redundancy(false,false);
            
        end
        
        %--------------------------------------------
        
        function [ newObj, newComplexPts ] = get_subdivided( obj, depth, complexPts, varargin )
            
            newObj = FacetComplex('empty');
            firstVertex = obj.get_max_vertices()+1;            
            tempPts = {};
            
            
            
            
            if isempty(varargin)
                
                
                for i = 1:length(obj.facets)
                    face = obj.facets{i};
                    [ newFaces, newPts ] = FacetComplex.face_subdiv( face, firstVertex, depth, complexPts(:,face) );
                    newObj.facets( (end+1):(end+length(newFaces)) ) = newFaces;
                    firstVertex = firstVertex + length(newPts);
                    
                    if ~isempty(newPts)
                        tempPts( (end+1):(end+length(newPts)) ) = newPts;
                    end
                end
                
                
            else
                indices = varargin{1};
                markedIndices = true(1,length(obj.facets));
                markedIndices(indices) = false;                
                facets1 = obj.facets(markedIndices);
                facets2 = {};
                
                for i = 1:length(indices)
                    face = obj.facets{indices(i)};
                    [ newFaces, newPts ] = FacetComplex.face_subdiv( face, firstVertex, depth, complexPts(:,face) );
                    facets2( (end+1):(end+length(newFaces)) ) = newFaces;
                    firstVertex = firstVertex + length(newPts);
                    
                    if ~isempty(newPts)
                        tempPts( (end+1):(end+length(newPts)) ) = newPts;
                    end
                end
                newObj.facets = [facets1,facets2];
                
            end
            
            if ~isempty(complexPts)
                newComplexPts = [ complexPts, cell2mat(tempPts) ];
            else
                newComplexPts = [];
            end
        end
        %--------------------------------------------
        
        
        
                        
        
        
        
        %--------------------------------------------
        %Slow, use sort_facet_indices(), then sorted=true in loops
                               
        function [ answer, indices ] = contains_face( obj, face, sorted )
            
            if sorted
                containsFace = @(x) FacetComplex.iscontained_sorted(face,x);
            else
                containsFace = @(x) FacetComplex.iscontained(face,x);
            end
            
            indices = find( cellfun(containsFace,obj.facets ) == 1 );
            answer = ~isempty(indices);
        end
        
        
        %--------------------------------------------
        function [ answer, index ] = contains_facet( obj, face, sorted )
            
            if sorted
                equivToFace = @(x) isequal(face,x);
            else
                equivToFace = @(x) FacetComplex.isequiv(face,x);
            end
            
            index = find( cellfun(equivToFace,obj.facets ) == 1);
            answer = ~isempty(index);
        end
        
        %--------------------------------------------
        
        function [ answer, indices ] = contains_faces( obj, faces, mode, varargin )
            
            if ~iscell(faces)
                faces = num2cell(faces,1);
            end
            
            if isempty(varargin); sorted = false;
            else sorted = strcmp(varargin{1},'sorted'); end
            
            switch mode
                
                case 'get_containing_indices'
                    
                    indices = [];
                    answer = true;
                    
                    for i = 1:length(faces)
                        [contained,newIndices] = obj.contains_face( faces{i}, sorted );
                        indices = union(indices,newIndices);
                        answer = answer && contained;
                    end
                    
                case 'get_contained_indices'
                    
                    indices = Inf(1,length(faces));
                    containedCount = 0;
                    
                    for i = 1:length(faces)
                        [contained,~] = obj.contains_face( faces{i}, sorted );
                        if contained
                            containedCount = containedCount+1;
                            indices(containedCount) = i;
                        end
                        
                    end
                    
                    indices = indices(1:containedCount);
                    answer = (containedCount==length(faces));
                    
                case 'get_none'    
                
                    
                    indices = [];
                    answer = true;
                    
                    for i=1:length(faces)
                        [contained,~] = obj.contains_face( faces{i}, sorted );
                        if ~contained
                            answer = false;
                            break;
                        end
                        
                    end
                    
                otherwise
                   
                    throw( MException('Invalid:Option','Invalid Option: %s',mode) );
 
            end
        end
        
        %--------------------------------------------
        
        function insert_facets( obj, faces, checkExisting, varargin )
            
            if ~iscell(faces)
                faces = num2cell(faces,1);
            end
            
            if checkExisting   %Slow, use sort_facet_indices(), then varargin{1}='sorted' in loops
                
                if isempty(varargin); sorted = false;
                else sorted = strcmp(varargin{1},'sorted'); end
                
                
                for i = 1:length(faces)
                    
                    [exists,~] = obj.contains_face(faces{i},sorted);
                    if ~exists
                        obj.facets{end+1} = faces{i};
                    end
                end
                
            else
                
                newNumFacets = length(obj.facets) + length(faces);
                obj.facets( (length(obj.facets)+1):newNumFacets ) = faces;
                
            end
        end
        
        
        
    
        %--------------------------------------------
        function assign_facets( obj, faces, checkRedundancy, varargin )
            
            if checkRedundancy
                obj.facets = {}; 
                obj.insert_facets( faces, true, varargin );
            
            else
                obj.facets = faces;
            end
            
            
            
            
        end
        
        
        %--------------------------------------------
        
        function delete_faces( obj, mode, array, varargin )
            
            switch mode
                case 'by_indices'
                    indices = array;
                case 'by_faces'
                    [~,indices] = obj.contains_faces(array,'get_containing_indices',varargin);
                otherwise
                    throw( MException('Invalid:Option','Invalid Option: %s',mode) );
            end
            
            obj.facets(indices) = [];

        end
        
        %--------------------------------------------
        
        function delete_facet_interiors( obj, mode, array, varargin )
            
            
            switch mode
                case 'by_indices'
                    indices = array;
                case 'by_faces'
                    [~,indices] = obj.contains_faces(array,'get_containing_indices',varargin);
                otherwise
                    throw( MException('Invalid:Option','Invalid Option: %s',mode) );
            end
            
            
            faces = obj.facets(indices);        
            obj.facets(indices) = [];
          
            for i = 1:length(faces)
                obj.insert_facets( FacetComplex.face_boundary( faces{i} ), false );
            end
        end
        %--------------------------------------------
        
        
        
        
        
        
        
        
        
        %--------------------------------------------
        
        function assign( obj, obj2, mode, varargin )
            
            switch mode
                case 'all'
                    obj.facets = obj2.facets;
                case 'trim_mem'
                    obj.facets = obj2.facets;
                case 'subcomplex'
                    if nargin > 0
                        obj.facets = obj2.facets( varargin{1} );
                    end
                otherwise
                    throw( MException('InvalidOption','Invalid Option: %s',mode) );
            end
            
        end
        %--------------------------------------------
        function answer = contained_in( obj, obj2, varargin )
                    
            [answer,~] = obj.contains_faces( obj2.facets, 'get_none', varargin );
            
        end
        
        %--------------------------------------------
        
        function intersect_with( obj, obj2, varargin )
            
          
            [~,indices1] = obj2.contains_faces(obj.facets,'get_contained_indices',varargin);
            [~,indices2] = obj.contains_faces(obj2.facets,'get_contained_indices',varargin);
            
            obj.assign( obj, 'subcomplex', indices1 );
            obj.insert_facets( obj2.facets(indices2), true, varargin );
            
            
        end
        
        %--------------------------------------------
        function union_with( obj, obj2, varargin )
            
            obj.insert_facets( obj2.facets, true, varargin );
        end
        %--------------------------------------------
        function disjoint_union_with( obj, obj2 )
                       
            firstVertex = obj.get_max_vertices();         
            vertexMap = @(x) x+firstVertex;
            
            obj.insert_facets( obj2.get_vertex_mapped(vertexMap).facets, false );
        end
        %--------------------------------------------
        function join_with( obj, obj2, varargin )
            
            firstVertex = obj.get_max_vertices();         
            vertexMap = @(x) x+firstVertex;
            facets2 = obj2.get_vertex_mapped(vertexMap).facets;
            
            newFacets = cell(1,length(obj.facets)*length(facets2));
            
            for i=1:length(obj.facets)
                
                k=(i-1)*length(facets2);
                
                for j=1:length(facets2)
                    newFacets{k+j} = [obj.facets{i};facets2{j}];
                end
            end
            
            obj.facets = newFacets;
        end
        %--------------------------------------------
        
        function subtract( obj, obj2, varargin )
            
            [~,containing] = obj.contains_faces(obj2.facets,'get_containing_indices',varargin);
            [~,contained] = obj.contains_faces(obj2.facets,'get_contained_indices',varargin);
            
            
            obj.assign( obj, 'subcomplex', setdiff( 1:length(obj.facets), [containing,contained] ) );
            
        end
        
        
        %--------------------------------------------

        
        
        
        %--------------------------------------------
        
        
        function [ plot, colors ] = plot2D( obj, complexPts, colors )
            
            
            vertexCounts = obj.get_facet_vertex_counts();
            
            triangleIndices = find( vertexCounts==3 );
            lineIndices = find( vertexCounts==2 );
            pointIndices = find( vertexCounts==1 );
            
            numTriangles = length(triangleIndices);
            numLines = length(lineIndices);
            numPoints = length(pointIndices);
            
            triangles = cell2mat( obj.facets(triangleIndices) );
            lines = cell2mat( obj.facets(lineIndices) );
            points = cell2mat( obj.facets(pointIndices) );
            
            trianglesX = NaN(3,numTriangles,'single');
            trianglesY = NaN(3,numTriangles,'single');
            
            for i=1:numTriangles
                trianglesX(1,i) = complexPts( 1, triangles(1,i) );
                trianglesY(1,i) = complexPts( 2, triangles(1,i) );
                trianglesX(2,i) = complexPts( 1, triangles(2,i) );
                trianglesY(2,i) = complexPts( 2, triangles(2,i) );
                trianglesX(3,i) = complexPts( 1, triangles(3,i) );
                trianglesY(3,i) = complexPts( 2, triangles(3,i) );
            end
            
            linesX = NaN(2,numLines,'single');
            linesY = NaN(2,numLines,'single');
            
            for i=1:numLines
                linesX(1,i) = complexPts( 1, lines(1,i) );
                linesY(1,i) = complexPts( 2, lines(1,i) );
                linesX(2,i) = complexPts( 1, lines(2,i) );
                linesY(2,i) = complexPts( 2, lines(2,i) );
            end
            
            pointsX = NaN(1,numPoints,'single');
            pointsY = NaN(1,numPoints,'single');
            
            for i=1:numPoints
                pointsX(i) = complexPts( 1, points(i) );
                pointsY(i) = complexPts( 2, points(i) );
            end
   
                        
            if isempty(colors)                
                colors = rand(1,length(obj.facets)); 
            end
            
            plot = fill(trianglesX,trianglesY,repmat(colors(triangleIndices),3,1),...
                linesX,linesY,repmat(colors(lineIndices),2,1),...
                pointsX,pointsY,colors(pointIndices));

            axis equal;
            
        end
        
        %--------------------------------------------
        
        
        function [ plot, colors ] = plot3D( obj, complexPts, colors )
            

            
            vertexCounts = obj.get_facet_vertex_counts();
            
            triangleIndices = find( vertexCounts==3 );
            lineIndices = find( vertexCounts==2 );
            pointIndices = find( vertexCounts==1 );
            
            numTriangles = length(triangleIndices);
            numLines = length(lineIndices);
            numPoints = length(pointIndices);
            
            triangles = cell2mat( obj.facets(triangleIndices) );
            lines = cell2mat( obj.facets(lineIndices) );
            points = cell2mat( obj.facets(pointIndices) );
            
            trianglesX = NaN(3,numTriangles,'single');
            trianglesY = NaN(3,numTriangles,'single');
            trianglesZ = NaN(3,numTriangles,'single');
            
            for i=1:numTriangles
                trianglesX(1,i) = complexPts( 1, triangles(1,i) );
                trianglesY(1,i) = complexPts( 2, triangles(1,i) );
                trianglesZ(1,i) = complexPts( 3, triangles(1,i) );
                trianglesX(2,i) = complexPts( 1, triangles(2,i) );
                trianglesY(2,i) = complexPts( 2, triangles(2,i) );
                trianglesZ(2,i) = complexPts( 3, triangles(2,i) );
                trianglesX(3,i) = complexPts( 1, triangles(3,i) );
                trianglesY(3,i) = complexPts( 2, triangles(3,i) );
                trianglesZ(3,i) = complexPts( 3, triangles(3,i) );
            end
            
            linesX = NaN(2,numLines,'single');
            linesY = NaN(2,numLines,'single');
            linesZ = NaN(2,numLines,'single');
            
            for i=1:numLines
                linesX(1,i) = complexPts( 1, lines(1,i) );
                linesY(1,i) = complexPts( 2, lines(1,i) );
                linesZ(1,i) = complexPts( 3, lines(1,i) );
                linesX(2,i) = complexPts( 1, lines(2,i) );
                linesY(2,i) = complexPts( 2, lines(2,i) );
                linesZ(2,i) = complexPts( 3, lines(2,i) );
            end
            
            pointsX = NaN(1,numPoints,'single');
            pointsY = NaN(1,numPoints,'single');
            pointsZ = NaN(1,numPoints,'single');
            
            for i=1:numPoints
                pointsX(i) = complexPts( 1, points(i) );
                pointsY(i) = complexPts( 2, points(i) );
                pointsZ(i) = complexPts( 3, points(i) );
            end
            
            if isempty(colors)                
                colors = rand(1,length(obj.facets)); 
            end
            
            plot = fill3(trianglesX,trianglesY,trianglesZ,repmat(colors(triangleIndices),3,1),...
                linesX,linesY,linesZ,repmat(colors(lineIndices),2,1),...
                pointsX,pointsY,pointsZ,colors(pointIndices));

            axis equal;
            
        end
        
        %--------------------------------------------
        
        
    end
    
    
    
end

