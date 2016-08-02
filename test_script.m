%Copyright (C) 2016 Piotr Beben

%This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
%This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
%You should have received a copy of the GNU General Public License along with this program. If not, see http://www.gnu.org/licenses/.
  


%%

%closed == true uses set 
%\bar{N}^\ell_j = {y\in\mc S | y'\in w for some w s.t. v_j\in w}

%Sensitivity is the learning rate.
%Higher values lead to faster convergence, 
%though poorer fittings and less stability. 
%Values between 0 and 0.1 seem to work best.



closed = true;
sensitivity = 0.1;
numPts = 200;
subdiv = 60; 
dim = 1; 

f=@(x)cos(x{1})./(0.2*x{1}+1);
g=@(x)sin(x{1})./(0.2*x{1}+1)+1;


cloud = PointCloud.build_random(2,numPts,'custom','normal',{f,g},0,10,0.0);

complex = FacetComplex('cube',subdiv);
complexPts = PointCloud.build_grid(subdiv,dim,[0;0]);

TestFunctions.test_2D(cloud,complex,complexPts,sensitivity,closed,200,1,0.0,'r',false,'figure_fit_all')


%%
closed = false;
sensitivity = 0.1;
numPts = 200;
subdiv = 5; 
dim = 0.5; 

f=@(x)cos(x{1})./(0.2*x{1}+1)+1;
g=@(x)sin(x{1})./(0.2*x{1}+1)+1;


cloud = PointCloud.build_random(2,numPts,'custom','normal',{f,g},0,10,0.0);



complex = FacetComplex('cube',[subdiv,subdiv]);
complex = complex.get_boundary();
%complex.delete_faces( 'by_indices', randi(complex.get_num_facets(),1,40) );

complexPts = PointCloud.build_grid([subdiv,subdiv],[dim,dim],[0;0]);


TestFunctions.test_2D(cloud,complex,complexPts,sensitivity,closed,200,1,0.0,'r',false,'figure_fit_all')




%%
closed = true;
sensitivity = 0.1;
numPts = 200;
subdiv = 4; 
dim = 0.5; 

f=@(x)cos(x{1});
g=@(x)sin(x{1});


cloud = PointCloud.build_random(2,numPts,'custom','normal',{f,g},0,10,0.0);

complex = FacetComplex('box',[subdiv,subdiv]);
complexPts = PointCloud.build_grid([subdiv,subdiv],[dim,dim],[2;0]);

TestFunctions.test_2D(cloud,complex,complexPts,sensitivity,closed,150,1,0.0,'r',false,'figure_fit_all')



%%

closed = false;
sensitivity = 0.1;
numPts = 100;
subdivX = 5;
subdivY = 1;
dim = 5; 

f=@(x)cos(x{1});
g=@(x)sin(x{1});
l=@(x)x{1};

cloud = PointCloud.build_random(2,numPts,'custom','normal',{f,g},-3.14,3.14,0.05);
cloud = [cloud, PointCloud.build_random(2,2*numPts,'custom','normal',{l,l},[-2,-2],[2,2],0.1)];


complex = FacetComplex('cube',[subdivX,subdivY]);
complex = complex.get_boundary();
complexPts = PointCloud.build_grid([subdivX,subdivY],[dim,dim],[0;0]);


TestFunctions.test_2D(cloud,complex,complexPts,sensitivity,closed,100,1,0.0,'r',false,'figure_fit_complex')




%%

closed = false;
sensitivity = 0.1;
numPts = 100;
subdivX = 50;
subdivY = 1;
dim = 3; 

f=@(x)cos(x{1});
g=@(x)sin(x{1});
l=@(x)x{1};

cloud = PointCloud.build_random(2,numPts,'custom','normal',{f,g},-3.14,3.14,0.05);
cloud = [cloud, PointCloud.build_random(2,2*numPts,'custom','normal',{l,l},[-2,-2],[2,2],0.1)];


complex = FacetComplex('cube',subdivX);
%complex = complex.get_boundary();
complexPts = PointCloud.build_grid(subdivX,dim,[0;0]);


TestFunctions.test_2D(cloud,complex,complexPts,sensitivity,closed,100,1,0.0,'r',false,'figure_fit_all')





%%

closed = false;
sensitivity = 0.1;
numPts = 150;
subdiv = 6; 
dim = 5; 

f=@(x)cos(x{1});
g=@(x)sin(x{1});
l=@(x)x{1};

cloud = PointCloud.build_random(2,numPts,'custom','normal',{f,g},-3.14,3.14,0.05);
cloud = [cloud, PointCloud.build_random(2,numPts,'custom','normal',{l,l},[-2,-2],[2,2],0.05)];


complex = FacetComplex('cube',[subdiv,subdiv]);
complexPts = PointCloud.build_grid([subdiv,subdiv],[dim,dim],[0;0]);


TestFunctions.test_2D(cloud,complex,complexPts,sensitivity,closed,200,1,0.0,'r',false,'figure_fit_all')



%%

closed = false;
sensitivity = 0.1;
numPts = 50;
subdiv = 5;
dim = 1; 


f=@(x)0.1*cos(x{1})+0.5;
g=@(x)0.1*sin(x{1})+0.5;



cloud1 = PointCloud.build_random(2,numPts,'centered','normal',[0;0],0.25);
cloud2 = PointCloud.build_random(2,2*numPts,'centered','normal',[1;0],0.05);
cloud3 = PointCloud.build_random(2,numPts,'centered','normal',[1;1],0.25);
cloud4 = PointCloud.build_random(2,3*numPts,'custom','normal',{f,g},-3.14,3.14,0.01);


cloud = [cloud1,cloud2,cloud3,cloud4];


complex1 = FacetComplex('cube',[subdiv,subdiv]);
complexPts1 = PointCloud.build_grid([subdiv,subdiv],[dim,dim],[0;0]);

complex2 = FacetComplex('cube',[subdiv,subdiv]);
complexPts2 = PointCloud.build_grid([subdiv,subdiv],[dim,dim],[2;2]);


complex = FacetComplex.disjoint_union(complex1,complex2);
complexPts = [complexPts1,complexPts2];

TestFunctions.test_2D(cloud,complex,complexPts,sensitivity,closed,150,1,0.0,'r',0.2,false,'figure_fit_all')


%%


%%
closed = true;
sensitivity = 0.1;
numPts = 1000;
subdivs = [5,5,5];
dims = [1,1,1];

f=@(x)cos(x{1}).*sin(x{2});
g=@(x)sin(x{1}).*sin(x{2});
h=@(x)cos(x{2});



cloud = PointCloud.build_random(3,numPts,'custom','normal',{f,g,h},[-3.14,-3.14],[3.14,3.14],0.0);


complex = FacetComplex('box',subdivs);
complexPts = PointCloud.build_grid(subdivs,dims,[2;0;0]);

TestFunctions.test_3D(cloud,complex,complexPts,sensitivity,closed,150,1,0.0,'r',0.7,true,'figure_fit_all')



%%
closed = false;
sensitivity = 0.1;
numPts = 800;
subdivs = [8,8,0];
dims = [1.5,1.5,0];

f=@(x)cos(x{1}).*sin(x{2});
g=@(x)sin(x{1}).*sin(x{2});
h=@(x)cos(x{2});
%l=@(x)x{1};


cloud = PointCloud.build_random(3,numPts,'custom','normal',{f,g,h},[-3.14,-3.14],[3.14,3.14],0.02);
%cloud = [cloud, PointCloud.build_random(3,0.2*numPts,'custom','normal',{l,l,l},[-2,-2],[2,2],0.05)];

complex = FacetComplex('cube',subdivs);
complexPts = PointCloud.build_grid(subdivs,dims,[0;0;-1]);



TestFunctions.test_3D(cloud,complex,complexPts,sensitivity,closed,150,1,0.0,'r',1,true,'figure_fit_data')


%%
closed = false;
sensitivity = 0.1;
numPts = 500;
subdivs = [3,3,3];
dims = [1,1,1];

f=@(x)cos(x{1}).*sin(x{2});
g=@(x)sin(x{1}).*sin(x{2});
h=@(x)cos(x{2});
l=@(x)x{1};


cloud = PointCloud.build_random(3,numPts,'custom','normal',{f,g,h},[-3.14,-3.14],[3.14,3.14],0.02);
cloud = [cloud, PointCloud.build_random(3,0.2*numPts,'custom','normal',{l,l,l},[-2,-2],[2,2],0.05)];

complex = FacetComplex('cube',subdivs);
complexPts = PointCloud.build_grid(subdivs,dims,[2;0;0]);



TestFunctions.test_3D(cloud,complex,complexPts,sensitivity,closed,100,1,0.0,'r',0.2,true,'figure_fit_all')



%%
closed = false;
sensitivity = 0.1;
numPts = 800;
subdivs = [3,3,3];
dims = [1,1,1];

f=@(x)cos(x{1}).*sin(x{2});
g=@(x)sin(x{1}).*sin(x{2});
h=@(x)cos(x{2});
l=@(x)x{1};


cloud = PointCloud.build_random(3,numPts,'custom','normal',{f,g,h},[-3.14,-3.14],[3.14,3.14],0.00);
cloud = [cloud, PointCloud.build_random(3,0.2*numPts,'custom','normal',{l,l,l},[-2,-2],[2,2],0.00)];

complex = FacetComplex('cube',subdivs);
complex = complex.get_boundary();
complex = complex.get_boundary();
complexPts = PointCloud.build_grid(subdivs,dims,[2;0;0]);


TestFunctions.test_3D(cloud,complex,complexPts,sensitivity,closed,100,1,0.0,'r',0.2,true,'figure_fit_all')


%%


closed = false;
sensitivity = 0.1;
numPts = 1000;
subdivs = [9,9,0];
dims = [3,3,0];

f=@(x)x{1};
g=@(x)x{2};
h=@(x)(1/3)*cos(5*x{1}).*sin(5*x{2})+(1/5)*(x{1}-x{2});

cloud = PointCloud.build_random(3,numPts,'custom','uniform',{f,g,h},[-1.5,-1.5],[1.5,1.5],0.02);


complex = FacetComplex('cube',subdivs);
complexPts = PointCloud.build_grid(subdivs,dims,[0;0;-1]);

TestFunctions.test_3D(cloud,complex,complexPts,sensitivity,closed,40,1,0.0,'r',0.2,true,'figure_fit_all')

%%

closed = false;
sensitivity = 0.1;
numPts = 800;
subdivs = [5,5,0];
dims = [0.5,0.5,0];

f=@(x)x{1};
g=@(x)x{2};
h=@(x)0.3*(x{1}-x{2});


cloud = PointCloud.build_random(3,numPts,'custom','uniform',{f,g,h},[-1,-1],[1,1],0.02);
complex = FacetComplex('cube',subdivs);
complex = complex.get_boundary();
complexPts = PointCloud.build_grid(subdivs,dims,[0;0;-1]);


TestFunctions.test_3D(cloud,complex,complexPts,sensitivity,closed,50,1,0.0,'r',0.2,true,'figure_fit_all')



%%


%%

sensitivity = 0.1;
numPts = 2000;


f1=@(x)cos(x{1});
f2=@(x)sin(x{1}).*cos(x{2});
f3=@(x)sin(x{1}).*sin(x{2}).*cos(x{3});
f4=@(x)sin(x{1}).*sin(x{2}).*sin(x{3});



cloud = PointCloud.build_random(4,numPts,'custom','normal',{f1,f2,f3,f4},[-3.14,-3.14,-3.14],[3.14,3.14,3.14],0.02);


subdivs = [4,4,4,0];
dims = [1,1,1,0];

complex = FacetComplex('cube',subdivs);
complexPts = PointCloud.build_grid(subdivs,dims,[2;0;0;0]);

TestFunctions.test_high_dimension(cloud,complex,complexPts,sensitivity,false,30,false)


subdivs = [3,3,3,3];
dims = [1,1,1,1];

complex = FacetComplex('box',subdivs);
complexPts = PointCloud.build_grid(subdivs,dims,[2;0;0;0]);

TestFunctions.test_high_dimension(cloud,complex,complexPts,sensitivity,true,30,false)

legend('3D mesh','4D mesh boundary')
