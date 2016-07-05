%Copyright (C) 2016 Piotr Beben

%This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
%This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
%You should have received a copy of the GNU General Public License along with this program. If not, see http://www.gnu.org/licenses/.
  


%%



%Higher values for sensitivity lead to faster convergence, 
%though poorer fittings and less stability. 
%Values between 0 and 0.1 seem to work best.

sensitivity = 0.1;
numPts = 100;
subdivX = 5; 
dim = 5; 

f=@(x)cos(x{1});
g=@(x)sin(x{1});
l=@(x)x{1};

cloud = PointCloud.build_random(2,numPts,'custom','normal',{f,g},[-3.14,-3.14],[3.14,3.14],0.05);
cloud = [cloud, PointCloud.build_random(2,2*numPts,'custom','normal',{l,l},[-2,-2],[2,2],0.1)];


complex = FacetComplex('grid',[subdivX,subdivX]);
complexPts = PointCloud.build_grid([subdivX,subdivX],[dim,dim],[0;0]);


TestFunctions.test_2D(cloud,complex,complexPts,sensitivity,55,1,0.0,'r',false,'figure_fit_complex')




%%


sensitivity = 0.1;
numPts = 100;
subdivX = 5;
subdivY = 1;
dim = 5; 

f=@(x)cos(x{1});
g=@(x)sin(x{1});
l=@(x)x{1};

cloud = PointCloud.build_random(2,numPts,'custom','normal',{f,g},[-3.14,-3.14],[3.14,3.14],0.05);
cloud = [cloud, PointCloud.build_random(2,2*numPts,'custom','normal',{l,l},[-2,-2],[2,2],0.1)];


complex = FacetComplex('grid',[subdivX,subdivY]);
complex = complex.get_boundary();
complexPts = PointCloud.build_grid([subdivX,subdivY],[dim,dim],[0;0]);


TestFunctions.test_2D(cloud,complex,complexPts,sensitivity,55,1,0.0,'r',false,'figure_fit_complex')



%%

sensitivity = 0.1;
numPts = 800;
subdivs = [5,5,0];
dims = [3,3,0];

f=@(x)x{1};
g=@(x)x{2};
h=@(x)0.7*cos(5*x{1}).*sin(3*x{2})+0.3*(x{1}-x{2});


cloud = PointCloud.build_random(3,numPts,'custom','uniform',{f,g,h},[-1,-1],[1,1],0.02);
complex = FacetComplex('grid',subdivs);
complexPts = PointCloud.build_grid(subdivs,dims,[0;0;-1]);


TestFunctions.test_3D(cloud,complex,complexPts,sensitivity,50,1,0.0,'r',true,'figure_fit_complex')


%%

sensitivity = 0.1;
numPts = 150;
subdivX = 3;
dim = 1;

f=@(x)cos(x{1}).*sin(x{2});
g=@(x)sin(x{1}).*sin(x{2});
h=@(x)cos(x{2});
l=@(x)x{1};


cloud = PointCloud.build_random(3,numPts,'custom','normal',{f,g,h},[-3.14,-3.14],[3.14,3.14],0.02);
cloud = [cloud, PointCloud.build_random(3,100,'custom','normal',{l,l,l},[-2,-2],[2,2],0.05)];

complex = FacetComplex('grid',[subdivX,subdivX,subdivX]);
complexPts = PointCloud.build_grid([subdivX,subdivX,subdivX],[dim,dim,dim],[0;0;0]);


TestFunctions.test_3D(cloud,complex,complexPts,sensitivity,80,1,0.0,'r',true,'figure_fit_data')





%%

sensitivity = 0.1;
numPts = 800;
subdivs = [5,5,0];
dims = [3,3,0];

f=@(x)x{1};
g=@(x)x{2};
h=@(x)0.3*(x{1}-x{2});


cloud = PointCloud.build_random(3,numPts,'custom','uniform',{f,g,h},[-1,-1],[1,1],0.02);
complex = FacetComplex('grid',subdivs);
complex = complex.get_boundary();
complexPts = PointCloud.build_grid(subdivs,dims,[0;0;-1]);


TestFunctions.test_3D(cloud,complex,complexPts,sensitivity,50,1,0.0,'r',true,'figure_fit_complex')



