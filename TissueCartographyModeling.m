%% Initialize ImSAnE project
%
% We start by clearing the memory and closing all figures.
%
clear all; close all;


[scriptPath,~,ext] = fileparts(matlab.desktop.editor.getActiveFilename);

dataDir = cd; 
projectDir = cd;
xp = project.Experiment(projectDir, dataDir);
cd(dataDir)

fileMeta                 = struct();
fileMeta.dataDir         = dataDir;
fileMeta.filenameFormat  = 'image.tif'; %yourfilename
fileMeta.nChannels       = 2; %number of channels
fileMeta.timePoints      = [1];
fileMeta.stackSize       = [1984,1984,136]; %size
fileMeta.stackResolution = [.0341 .0341 .180]; %imageresolution
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

%For Single Time
myDetectOpts = struct('channel', 1, 'sigma', 2, 'ssfactor',4,...
           'rmRadialOutliers', 10,'thresh',.6,'amin', 200,'dildisc',4,...
           'fileName',[xp.fileMeta.dataDir, '/3rdpnutRNAiLAGPnut1.31.21'],...
           'foreGroundChannel',1,'zdim',2); 
xp.setDetectOptions(myDetectOpts);

%% Prepare subsampled file for Ilastik

xp.detector.prepareIlastik(xp.stack);
%CHANGED ilastikDetector.m under +surfaceDetection to get xyzct etc right.
%Revert for other image types

%% Detect the surface

xp.detectSurface(); 
%% Inspect the point cloud in a cross section

inspectOptions= struct('dimension', 'z', 'value', 400, 'pointCloud', 'r');
xp.detector.inspectQuality(inspectOptions, xp.stack);

%for i = 1:20:500
%    inspectOptions= struct('dimension', 'z', 'value', i, 'pointCloud', 'r');
%    xp.detector.inspectQuality(inspectOptions, xp.stack);
%    pause(.2)
%end

%% Inspect pointcloud in 3d
ssfactor = 6;
xp.detector.pointCloud.inspect(ssfactor);

%% Save pointCloud to obj format

clear OBJ
OBJ.vertices = xp.detector.pointCloud.unalignedPoints;
OBJ.objects(1).type='f';
OBJ.objects(1).data.vertices=[];
write_wobj(OBJ, fullfile(projectDir, 'image.obj')); %rename to desired file name

%% Read surface mesh produced by meshlab
outputMeshFormat = 'E:/imsane-master/image.ply';  %rename to desired file name
outputMesh = sprintf(outputMeshFormat,xp.currentTime);

meshPermutation = [1,2,3]; % There is an axis permutation in the pointCloud, 
                           % that we need to take care of to have the
                           % surface axes match the data axes. 
                           
% read the output mesh, and crate a mesh structure. 
mesh = read_ply_mod(outputMesh);
mesh.v = mesh.v(:,meshPermutation); 
mesh.vn = mesh.vn(:,meshPermutation);

%% create seeds for the centers of the charts and load into fitter

zdir = 1;
temp = find(mesh.v(:,zdir) == min(mesh.v(:,zdir))); %Finds minimum value in first column
seeds(1) = temp(1);

%% before defining the options, we compute the data sets orientation since 
% we prefer having a fixed frame of reference. 
points = mesh.v;
pc = surfaceDetection.PointCloud(points);
pc.determineROI(5);%5
rotation    = pc.ROI.rotation;
translation = pc.ROI.translation;

%% Set seeds points

%Defined
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
% 
% We call the process of creating the 2D maps of the surface data "pulling
% back the data" (in accordance with standard mathematical terminology)
% Calling the function SurfaceOfInterest.pullbackStack generates the Field 
% objects in SOI.fields that contain these pullbacks.
% The function is called as:
%
% pullbackStack(stack, ROI, time, onionOpts)
%
% Where the arguments are:
%
% stack:    Stack object containing the data
% ROI:      RegionOfInterest object containing affine
%           transformation from embedding coordinates to stack
%           ROI can be left empty: ROI = []
% time:     the time of the pullback, for storing in the SOI
%
% onionOpts: multi-layer pullback options, structure with fields
%   - nLayers:          number of layers
%   - layerDistance:    spacing of layers
%   - sigma:            smoothing of surface for layer
%                       displacement
%   - makeIP:           intensity projection: 'MIP', 'SIP',
%                       'both'
%   - IPonly:           only IP is stored as a field
%                       WARNING: this option will delete
%                       previously stored other layers if they
%                       were made

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

% the summed intensity projection of each region at the current time

dataField = xp.SOI.getField('data_MIP');
tidx = xp.tIdx(xp.currentTime);
data = dataField(tidx)


%% 
% NOTE: For more information on the organization of the data structures, 
% see the supplemental information of the manuscript.

% the color version of the map (pullback) to the plane, stored to be used 
% as texture for the 3D rendering the next block
color = {};
separation = 5

figure,
%for i = 1:numel(data.patches)
for i = 1:2   
    % the two channels
    R = mat2gray(data.patches{i}.apply{1},[0 500]); %Almost max bit depth
    %R = mat2gray(gdata.patches{i}.apply{1},[0 4]);
    %colormap(gray)
    G = mat2gray(data.patches{i}.apply{2},[0,500]);
    %G = mat2gray(gdata.patches{i}.apply{2},[0,500]);
    %colormap(gray)
    %B = mat2gray(data.patches{i}.apply{3},[0,600]);
    % a little bit of manual adjustment of the lookup table
    %R = imadjust(R, [0 0.6]);
    %G = imadjust(G, [0 0.8]);

    % concatenate to make color image
    color{i} = cat(3,R,G,G);
    %color{i} = R
    % make the background white 
   % color{i}(color{i}==0) = 1;

    % show the m
    subplot(ceil(numel(data.patches)/2), 2, i)
    imshow(permute(color{i}, [2 1 3]),[],'InitialMagnification',100)
end

%% Visualize in 3d

separation = 0;

figure
hold on;

%for i = 1:numel(data.patches)
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

%% Record rotating Video

OptionZ.FrameRate=30;OptionZ.Duration=15;OptionZ.Periodic=true;
%CaptureFigVid([-20,10;-110,10;-190,80;-290,10;-380,10], 'WellMadeVid',OptionZ)
CaptureFigVid([-20,10;-110,10;-190,10;-290,10;-380,10], 'movie',OptionZ)

%% Loop for dynamic atlas generation

% % keep track of seed position propagation in time
% VorSeedsX{1} = mesh.v(VorSeeds,:);
% diskSeedsX{1} = mesh.v(diskSeeds,:);
% 
% % point in parametrization to hold fixed next round
% fixedPtU = {};
% fixedPtX = {};
% for k = 1:numel(xp.fitter.fittedParam.submeshes)
% %for k = 3
%     subm = xp.fitter.fittedParam.submeshes{k};
%     fixedPtUtmp = sqrt(mean(subm.u{1}.^2)/2);
%     fixedPtIdx = pointMatch(fixedPtUtmp, subm.u{1});
%     fixedPtX{k} = subm.v(fixedPtIdx,:);
%     fixedPtU{k} = subm.u{1}(fixedPtIdx,:);
% end

for t = fileMeta.timePoints([10:end]) 
    
    t
    
    tidx = xp.tIdx(t);
    
    % raw data loading
    xp.loadTime(t);
    xp.rescaleStackToUnitAspect();
    

    % Read surface mesh produced by meshlab
    outputMesh = sprintf(outputMeshFormat, t);
    mesh = read_ply_mod(outputMesh);
    mesh.v = mesh.v(:,meshPermutation); 
    mesh.vn = mesh.vn(:,meshPermutation);
%     
%     nVorSeeds = 2;
%     nDiskSeeds = 0;
%     VorSeeds = floor(rand([nVorSeeds 1])*size(mesh.v,1))+1;
%     diskSeeds = floor(rand([nDiskSeeds 1])*size(mesh.v,1))+1;

%     %%%%%% Generate seeds at anterior/posterior most positions %%%%%%%%%
%      VorSeedsXinit = [max(mesh.v(:,1)) median(mesh.v(:,2)) median(mesh.v(:,3)) ; min(mesh.v(:,1)) median(mesh.v(:,2)) median(mesh.v(:,3))]
%      diskSeedsXinit = [max(mesh.v(:,1)) median(mesh.v(:,2)) median(mesh.v(:,3))]
%  
%      diskSeeds = pointMatch(diskSeedsXinit, mesh.v);
%      VorSeeds = pointMatch(VorSeedsXinit, mesh.v);

%%%%% Generate seeds based off of predetermined list for nurse cell %%%%%%%

    VorSeedsXinit = pullbackCenters{t};
    diskSeedsXinit = [max(mesh.v(:,1)) median(mesh.v(:,2)) median(mesh.v(:,3))]
    
    diskSeeds = pointMatch(diskSeedsXinit, mesh.v);
    VorSeeds = pointMatch(VorSeedsXinit, mesh.v);

    fitOptions = struct('VorSeeds', VorSeeds, 'transitionWidth', 250,...
        'diskSeeds', diskSeeds, 'diskRadius', 40, 'makeTMaps', false);
    xp.setFitOptions(fitOptions);
    
  %  xp.setFitOptions(fitOptions);
    xp.fitSurface(mesh); 
    
    % populate SOI
    xp.fitter.populateSOI(xp.SOI, xp.currentTime);

    % Pullback the stack to the desired charts
    xp.SOI.pullbackStack(xp.stack, [], xp.currentTime, onionOpts);
    
    %Calcualte g and detG

    xp.SOI.NCalcInducedMetric(t);

    gField = xp.SOI.getField('metric');
    gtidx = xp.tIdx(xp.currentTime);
    gdata = gField(gtidx);

    for k = 1 : length(gdata.patches)
    sqrtdetg{gtidx,k} = sqrt(gdata.patches{k}.apply{1,1}.*gdata.patches{k}.apply{2,2}- ...
              gdata.patches{k}.apply{1,2}.*gdata.patches{k}.apply{2,1});
    end
    
    % Generate Figure
    dataField = xp.SOI.getField('data_MIP');
    tidx = xp.tIdx(xp.currentTime);
    data = dataField(tidx);

    % NOTE: For more information on the organization of the data structures, 
    % see the supplemental information of the manuscript.

    % the color version of the map (pullback) to the plane, stored to be used 
    % as texture for the 3D renderinging the next block
    color = {};

    figure,
    %for i = 1:numel(data.patches)
     for i = 1:2   
        % the two channels
        R = mat2gray(data.patches{i}.apply{1},[0 15000]);
        G = mat2gray(data.patches{i}.apply{2},[0 20000]);

        % a little bit of manual adjustment of the lookup table
        %R = imadjust(R, [0 0.6]);
        %G = imadjust(G, [0 0.8]);

        % concatenate to make color image
        color{i} = cat(3,G,R,G);
        %color{i} = R;
        % make the background white 
       % color{i}(color{i}==0) = 1;

        % show the m
    end

    %% Visualize in 3d

    separation = 3;

    figure
    hold on;

    %for i = 1:numel(data.patches)
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

    hold off;

    threeDFileName = strcat('3D_13pix_', num2str(xp.currentTime));
    savefig(threeDFileName)
    close all
    
    
end

%% Save the surface of interest to disc

imwriteOptions = {'tif'};
saveDir = fullfile(projectDir, '181105_13pixels');
options = struct('dir',saveDir,'imwriteOptions',{imwriteOptions},...
                    'make8bit',false);
xp.SOI.save(options)


%% Generate sqrtdetg * intensity pullback

intImages = dir('june2021image2.tif');
for i = 10
    intFileName = intImages(i).name
    temp = sqrtdetg{i,1} .* imread(intFileName);
end











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
%figure('name','Triangle Mesh Curvature Example','numbertitle','off','color','w');
cmp = jet;
colormap(cmp(end:-1:1,:))

%curvature_Type = GausianCurvature;
curvature_Type = MeanCurvature;

%caxis([min(curvature_Type) max(curvature_Type)]); % color overlay the gaussian curvature
mesh_h=patch(FV,'FaceVertexCdata',curvature_Type','facecolor','interp','edgecolor','interp','EdgeAlpha',0.2);

%set some visualization properties
set(mesh_h,'ambientstrength',1);

axis on
%view([-45,35.2]);
xlabel('x')
ylabel('y')
zlabel('z')

camlight();
lighting phong
%colorbar();
axis equal
caxis([-0.02 0.04]);
%caxis([-0.001 0.001]);



% FV is the mesh that curvature_Type is being mapped onto.
% Access FV with FV.vertices([],[]) [row],[column].  Examples
% FV.vertices(:,1) = entire first column FV.vertices(1,:) = entire first
% row

%Plot with - plot3(FV.vertices(:,1), FV.vertices(:,2), FV.vertices(:,3))

%Want to line up each xyz coordinate with its respective gausian curvature
%value.  They are both 7306 at least. 

%% Save Curvature 3D model

curvatureFileName = strcat('3Dcurvature', num2str(xp.currentTime));

savefig(curvatureFileName)

%% Make Curvature Movie

v = VideoWriter('3Dcurvature.avi')
v.FrameRate = 10;
movieRotation = -15;
movieRotation2 = 36;
open(v);
for i = 1:288
    movieRotation = movieRotation + 2.0;
    movieRotation2 = movieRotation2 - .25;
    view([movieRotation,movieRotation2])
    frame = getframe(gcf);
    writeVideo(v, frame);
end
close(v);

%% Batch process curvature figures

for t = fileMeta.timePoints([1,2,4:end]) 
    
    tidx = xp.tIdx(t);
    
    % raw data loading
    xp.loadTime(t);
    xp.rescaleStackToUnitAspect();
    

    % Read surface mesh produced by meshlab
    outputMesh = sprintf(outputMeshFormat, t);
    mesh = read_ply_mod(outputMesh);
    mesh.v = mesh.v(:,meshPermutation); 
    mesh.vn = mesh.vn(:,meshPermutation);

    %Curvature Calculation
    FV.vertices = mesh.v;
    FV.faces    = mesh.f;

    getderivatives=0;
    [PrincipalCurvatures,PrincipalDir1,PrincipalDir2,FaceCMatrix,VertexCMatrix,Cmagnitude]= GetCurvatures( FV ,getderivatives);

    GausianCurvature=PrincipalCurvatures(1,:).*PrincipalCurvatures(2,:);
    MeanCurvature = 0.5*(PrincipalCurvatures(1,:)+PrincipalCurvatures(2,:));
    figure
    
    %figure('name','Triangle Mesh Curvature Example','numbertitle','off','color','w');
    cmp = jet;
    colormap(cmp(end:-1:1,:))

    %curvature_Type = GausianCurvature;
    curvature_Type = MeanCurvature;

    %caxis([min(curvature_Type) max(curvature_Type)]); % color overlay the gaussian curvature
    mesh_h=patch(FV,'FaceVertexCdata',curvature_Type','facecolor','interp','edgecolor','interp','EdgeAlpha',0.2);

    %set some visualization properties
    set(mesh_h,'ambientstrength',1);

    axis off
    %view([-45,35.2]);
    camlight();
    lighting phong
    %colorbar();
    axis equal
    caxis([-0.02 0.04]);
    %caxis([-0.001 0.001]);
    
    %Save
    
    curvatureFileName = strcat('curvature_', num2str(xp.currentTime));
    savefig(curvatureFileName)
    close all
end







