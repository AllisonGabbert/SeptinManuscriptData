%% Initialize ImSAnE project
% Start by installing four toolboxes from MATLAB - Deep Learning Toolbox,
% Image Processing Toolbox, Statistics and Machine Learning Toolbox, and
% Curve Fitting Toolbox. 

% Then, download the ImSAnE repo. Run the setup.m file from that to
% initialize ImSAnE. Make sure that ImSAnE, this script, and your images are
% all in the same directory. 


% We start by clearing the memory and closing all figures.
clear all; close all;

[scriptPath,~,ext] = fileparts(matlab.desktop.editor.getActiveFilename);

dataDir = cd; 
projectDir = cd;
xp = project.Experiment(projectDir, dataDir);
cd(dataDir)

fileMeta = struct();
fileMeta.dataDir         = dataDir;

% Add your file name to the following line in place of MYFILENAME. Be sure this file is in the same
% directory as this script and the ImSAnE script
fileMeta.filenameFormat  = 'MYFILENAME.tif'; %yourfilename

% When converting your file to a TIF file in ImageJ, it is easiest to reduce
% the number of channels to just the channel for the membrane marker. If you
% decide to use multiple channels, change the number in the following line
fileMeta.nChannels       = 1; %number of channels
fileMeta.timePoints      = [1];

% In ImageJ, select 'Image', then 'show info' and input the width, height, and depth into
% the following line
fileMeta.stackSize       = [960,960,162]; %[width, height, depth]

% In ImageJ, select 'Image', then 'properties' and input the pixel width and
% height and the voxel depth in the following line
fileMeta.stackResolution = [.030 .030 .160]; %[pixel width, pixel height, voxel depth]
fileMeta.swapZT          = 0;

expMeta                  = struct();
expMeta.channelsUsed     = [1 2];
expMeta.channelColor     = [1 2];
expMeta.description      = 'Channel1 Channel2';
expMeta.dynamicSurface   = 1;
expMeta.jitterCorrection = 0; % 1: Correct for sample translation
expMeta.fitTime          = fileMeta.timePoints(1); 
expMeta.detectorType     = 'surfaceDetection.ilastikDetector'; 
expMeta.fitterType       = 'surfaceFitting.meshWrapper'; 

xp.setFileMeta(fileMeta);
xp.setExpMeta(expMeta);
xp.initNew();

%% Load data for surface detection and rescale to unit aspect ratio
xp.loadTime(1);
xp.rescaleStackToUnitAspect();

%% Detect the surface

% For Single Time
% Input your file name in the myDetectOpts function below, in place of '/MYFILENAME'
myDetectOpts = struct('channel', 1, 'sigma', 2, 'ssfactor',4,...
           'rmRadialOutliers', 10,'thresh',.6,'amin', 200,'dildisc',4,...
           'fileName',[xp.fileMeta.dataDir, '/MYFILENAME'],...
           'foreGroundChannel',1,'zdim',2); 
xp.setDetectOptions(myDetectOpts);

%% Prepare subsampled file for Ilastik

xp.detector.prepareIlastik(xp.stack);

% There should now be an h5 file in your directory. Open this in Ilastik and
% follow the directions given in the Star Protocol. The command window will
% indicate a pause button; press any key in the command window to resume
% once you've finished work in Ilastik. 
pause_button = "on. Switch to Ilastik and do work there, then press a key below to continue"
pause 
%% Detect the surface

xp.detectSurface(); 
%% Inspect the point cloud in a cross section

inspectOptions= struct('dimension', 'z', 'value', 400, 'pointCloud', 'r');
xp.detector.inspectQuality(inspectOptions, xp.stack);

%% Inspect pointcloud in 3d
ssfactor = 6;
xp.detector.pointCloud.inspect(ssfactor);

%% Save pointCloud to obj format

clear OBJ
OBJ.vertices = xp.detector.pointCloud.unalignedPoints;
OBJ.objects(1).type='f';
OBJ.objects(1).data.vertices=[];

% Input your file name below in place of 'MYFILENAME'
write_wobj(OBJ, fullfile(projectDir, 'MYFILENAME.obj')); %rename to desired file name

% There should now be an OBJ file in your directory. Open this in Meshlab and
% follow the directions given in the Star Protocol. The command window will
% indicate a pause button; press any key in the command window to resume
% once you've finished work in Meshlab. 
pause_button = "on. Do work in Meshlab and then press a key to continue"
pause 
%% Read surface mesh produced by meshlab

% Input your file name below in place of 'MYFILENAME'
outputMeshFormat = 'MYFILENAME.ply';  %rename to desired file name
outputMesh = sprintf(outputMeshFormat,xp.currentTime);

% There is an axis permutation in the pointCloud,that we need to take care of 
% to have thesurface axes match the data axes. 
meshPermutation = [1,2,3]; 
                           
% Reads the output mesh, and creates a mesh structure. 
mesh = read_ply_mod(outputMesh);
mesh.v = mesh.v(:,meshPermutation); 
mesh.vn = mesh.vn(:,meshPermutation);

%% Create seeds for the centers of the charts and load into fitter

zdir = 1;
temp = find(mesh.v(:,zdir) == min(mesh.v(:,zdir))); %Finds minimum value in first column
seeds(1) = temp(1);

%% Before defining the options, we compute the data sets orientation since we prefer having a fixed frame of reference. 
points = mesh.v;
pc = surfaceDetection.PointCloud(points);
pc.determineROI(5);%5
rotation    = pc.ROI.rotation;
translation = pc.ROI.translation;

%% Set seeds points

% Defined
VorSeedsXinit = [max(mesh.v(:,1)) median(mesh.v(:,2)) median(mesh.v(:,3)) ; min(mesh.v(:,1)) median(mesh.v(:,2)) median(mesh.v(:,3))]
diskSeedsXinit = [max(mesh.v(:,1)) median(mesh.v(:,2)) median(mesh.v(:,3))]

diskSeeds = pointMatch(diskSeedsXinit, mesh.v);
VorSeeds = pointMatch(VorSeedsXinit, mesh.v);

fitOptions = struct('VorSeeds', VorSeeds, 'transitionWidth', 10,...
    'diskSeeds', diskSeeds, 'diskRadius', 40, 'makeTMaps', false);
xp.setFitOptions(fitOptions);
xp.fitSurface(mesh);


%% Visualize overlapping regions

xp.fitter.inspectMesh(1:2);
view([0 0 1]);

%% Submesh

xp.fitter.inspectMesh(1);

%% Generate the atlas

xp.generateSOI();


%% Create maps of the surface data

onionOpts = struct('nLayers', 131, 'layerDistance', 1, 'sigma', 1,...
                    'makeIP', 'MIP', 'IPonly', false);
 
xp.SOI.pullbackStack(xp.stack, [], xp.currentTime, onionOpts);

%% Calcualte g and detG

xp.SOI.NCalcInducedMetric(1);

gField = xp.SOI.getField('metric');
gtidx = xp.tIdx(1);
gdata = gField(gtidx);

for k = 1 : length(gdata.patches)
sqrtdetg{gtidx,k} = sqrt(gdata.patches{k}.apply{1,1}.*gdata.patches{k}.apply{2,2}- ...
          gdata.patches{k}.apply{1,2}.*gdata.patches{k}.apply{2,1});
end

%% Visualize in 2d

% The summed intensity projection of each region at the current time

dataField = xp.SOI.getField('data_MIP');
tidx = xp.tIdx(xp.currentTime);
data = dataField(tidx)


% Creates the color version of the map (pullback) to the plane, stored to be used 
% as texture for the 3D rendering the next block
color = {};
separation = 5
figure,
for i = 1:2   
    % the two channels
    R = mat2gray(data.patches{i}.apply{1},[0 500]); %Almost max bit depth
    G = mat2gray(data.patches{i}.apply{2},[0,500]);

    % concatenate to make color image
    color{i} = cat(3,R,G,G);

    % show the m
    subplot(ceil(numel(data.patches)/2), 2, i)
    imshow(permute(color{i}, [2 1 3]),[],'InitialMagnification',100)
end

%% Visualize in 3d

separation = 0;

figure
hold on;

% for i = 1:numel(data.patches)
for i = 1:2    
    % the 3D coordinates of a region
    X = xp.SOI.embedding(tidx).patches{i}.apply;
    
    % the mean normal direction
    d = mean(xp.fitter.fittedParam.submeshes{i}.vn);
    d = d./norm(d);
    d = separation*d;
    
    % show each region in 3D displaced along its normal by separation
    surf(X{1} + d(1),X{2} + d(2),X{3} + d(3), color{i}, 'FaceColor','texturemap');  
end

axis equal;
axis on;
shading flat
set(gcf,'color','w');
view(60, -60); % set viewing angle
colormap('gray')
hold off;


%% Curvature

FV.vertices = mesh.v;
FV.faces    = mesh.f;

getderivatives=0;
[PrincipalCurvatures,PrincipalDir1,PrincipalDir2,FaceCMatrix,VertexCMatrix,Cmagnitude]= GetCurvatures( FV ,getderivatives);

GausianCurvature=PrincipalCurvatures(1,:).*PrincipalCurvatures(2,:);
MeanCurvature = 0.5*(PrincipalCurvatures(1,:)+PrincipalCurvatures(2,:));


%% Visualize the curvature atlas in 3D

separation = 100;

figure
cmp = jet;
colormap(cmp(end:-1:1,:))

curvature_Type = MeanCurvature;

mesh_h=patch(FV,'FaceVertexCdata',curvature_Type','facecolor','interp','edgecolor','interp','EdgeAlpha',0.2);

%set some visualization properties
set(mesh_h,'ambientstrength',1);

axis on
xlabel('x')
ylabel('y')
zlabel('z')

camlight();
lighting phong
axis equal
caxis([-0.02 0.04]);




%% Save Curvature 3D model
pause_button = "on. Change the file name in lines following line 285 and then press a key to continue"
pause %Be sure to change file name

% Input your file name below in place of MYFILENAME
curvatureFileName = strcat('MYFILENAME', num2str(xp.currentTime));

savefig(curvatureFileName)

