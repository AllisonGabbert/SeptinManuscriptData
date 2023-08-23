%% Spectral Analysis on a set of PLYs from the same condition/genotype
% Here, we perform a spectral decomposition of the shapes of a set of 
% surfaces stored on disk as PLY files. 
%
% The weight of each mode is the amount of deformation amplifying a mode
% with unit norm. Note that vecnorm(V', 2, 2) == 1, so each mode has a
% Euclidean norm of 1, and the amount of spectral weight amplifies this
% deformation to match the amount of that pattern measured in the mesh
% relative to a reference sphere. The reference sphere is obtained from the
% mesh by conformalized mean curvature flow.
%
% This workflow is based on code from:
% N. P. Mitchell*, D. J. Cislo*, TubULAR: Tracking in toto deformations of
% dynamic tissues via constrained maps. [biorxiv]
% but is here adapted to use spherical harmonics

%% First change directories to where PLY meshes of cell clusters exist

% Running this line in the script grabs the path of this script
mfilePath = mfilename('fullpath');
if contains(mfilePath,'LiveEditorEvaluationHelper')
    mfilePath = matlab.desktop.editor.getActiveFilename;
end
[ scriptDir, scriptName ] = fileparts(mfilePath) ; 
cd(scriptDir)
close all
clc
clearvars

%%
% Define where the PLY surfaces are stored, as datadir
% --> CHANGE THIS LINE TO REFLECT YOUR LOCAL FILEPATH
datadir = '../../3D_surface_models/sqh_knockdown/';

% resolution of the PLYs, in um / voxel width. Make sure the PLYs are saved
% with this resolution, and in isotropic resolution (same for x, y, and z)
pix2um = 0.0302 ;  % um / pixel

% Add DECLab to the path: https://github.com/DillonCislo/DECLab
% This is a submodule of the SeptinManuscriptData repo.
addpath(genpath('../../DECLab/'))

% Also download and add gptoolbox to the path: https://github.com/alecjacobson/gptoolbox
addpath(genpath('../gptoolbox/mesh/'))

blue    = [0.0000, 0.4470, 0.7410] ; % 1
red     = [0.8500, 0.3250, 0.0980] ; % 2
yellow  = [0.9290, 0.6940, 0.1250] ; % 3
purple  = [0.4940, 0.1840, 0.5560] ; % 4 
green   = [0.4660, 0.6740, 0.1880] ; % 5 
sky     = [0.3010, 0.7450, 0.9330] ; % 6 
maroon  = [0.6350, 0.0780, 0.1840] ; % 7
gray    = [0.2000, 0.2000, 0.2000] ; % 8
brown   = [0.5400, 0.2500, 0.0900] ; % 9
teal    = [0.0000, 0.6500, 0.5200] ; % 10
pink    = [1.0000, 0.5137, 0.5137] ; % 11
dark_green =  [0.0392, 0.5059, 0.2745] ;   % 12

colors = [blue; red; yellow; purple; green; ...
          sky; maroon; gray; brown; teal; pink; dark_green] ;
        
%% Default Options
% overwrite previous results if on disk
overwrite = true ;

% The number of Laplacian eigenvectors to calculate
nModes = 1000;

fns = dir(fullfile(datadir, '*.ply')) ;

% The output directory is within the data directory 'datadir'
outdir = fullfile(datadir, 'analysis') ;
if ~exist(outdir, 'dir')
    mkdir(outdir)
end

%% Loop over each PLY file
for ii = 1:length(fns)
    fn = strrep(fns(ii).name, '.ply', '') ;
    disp(['ii = ' num2str(ii) ': ' fn])

    outfn1 = fullfile(outdir, ...
        [fn '_powerSpectrum_radialu.mat']) ;
    if ~exist(outfn1, 'file') || overwrite 

        mesh = read_ply_mod(fullfile(fns(ii).folder, [fn '.ply'])) ;
        % Convert mesh to microns
        mesh.v = mesh.v * pix2um ;
        
        %% Conformally map to the unit sphere
        clf
        Ufn = fullfile(outdir, ...
            [fn '_conformalMappingToUnitSphere.mat']) ;
        if ~exist(Ufn, 'file') || overwrite
            U = conformalized_mean_curvature_flow(mesh.v,mesh.f, 'LaplacianType', 'cotangent') ;
            Urescaled = conformalized_mean_curvature_flow(mesh.v,mesh.f, 'LaplacianType', 'cotangent', 'RescaleOutput', true) ;

            [Ur2fit,ind2fit]= farthest_points(Urescaled, 200) ;
            offset= mean(Ur2fit) ;
            ptCloud = pointCloud(Ur2fit - offset) ;
            sphereModel = pcfitsphere(ptCloud,1) ;
            sphereCenter = sphereModel.Center + offset ;
            sphereRadius = sphereModel.Radius ;
            sphereParameters = sphereModel.Parameters ;
            radii = vecnorm(mesh.v - sphereCenter, 2, 2) ;
            radii0 = vecnorm(Urescaled - sphereCenter, 2, 2) ;

            % Check if we get a more uniform radius from using offset in addition
            tmp = mean([offset; sphereCenter]) ;
            if std(radii0) > std(vecnorm(Urescaled - tmp, 2, 2))
                if std(vecnorm(Urescaled - tmp, 2, 2)) < std(vecnorm(Urescaled - tmp, 2, 2))
                    disp('Using FarthestPtSearch sampling mean only, since here it is most accurate')
                    sphereCenter = offset ;
                else
                    disp('Using 1/2 * (sphereModel.Center + FarthestPtSearch sampling mean)')
                    sphereCenter = tmp ;
                end
                sphereRadius = mean(vecnorm(Ur2fit - sphereCenter, 2, 2)) ;
            end
            radii = vecnorm(mesh.v - sphereCenter, 2, 2) ;
            radii0 = vecnorm(Urescaled - sphereCenter, 2, 2) ;

            %% plot the result
            close all
            figure('position', [0, 0, 1000, 1000])
            subplot(2, 2, 1)
            trisurf(triangulation(mesh.f, mesh.v), 'facecolor', blue, ...
                 'edgecolor', 'none', 'facealpha', 0.2); 
            hold on;
            trisurf(triangulation(mesh.f, Urescaled), 'facecolor', [0.5, 0.5, 0.5], ...
                 'edgecolor', 'none', 'facealpha', 0.2); 
            hold on; 
            view([65, 25])
            axis off

            camlight
            axis equal ;
            grid off ;
            trisurf(triangulation(mesh.f, mesh.v), 'facecolor', [0.2, 0.2, 1],...
                 'edgecolor', 'none', 'facealpha', 0.2)
            axis equal; 
            grid off ; % axis off ;
            subplot(2, 2, 2)
            trisurf(triangulation(mesh.f, Urescaled), ...
                 'edgecolor', 'none', 'facevertexCdata', radii-radii0)
            axis equal ;
            grid off ;
            title('radial position of surface relative to sphere');
            clims = caxis ;
            view([65, 25])
            caxis(max(abs(radii-radii0)) * [-1,1])
            try
                colormap(bam)
            catch
            end
            c = colorbar ;
            c.Label.String = '\delta R \equiv R-R_0' ;
            axis off

            subplot(2, 2, 3)
            trisurf(triangulation(mesh.f, mesh.v), ...
                 'edgecolor', 'k', 'facecolor', 'w')
            axis equal ;
            grid off
            caxis(clims)
            view([65, 25])
            axis off
            
            subplot(2, 2, 4)
            trisurf(triangulation(mesh.f, Urescaled), ...
                 'edgecolor', 'k', 'facecolor', 'w')
            axis equal ;
            grid off
            caxis(clims)
            view([65, 25])
            axis off

            set(gcf, 'color', 'w')
            figfn = fullfile(outdir,[ fn '_sphericalFit.png']) ;
            % saveas(gcf, figfn)
            export_fig(figfn, '-r300')
            close all

            %% save the result
            save(Ufn, 'U', 'Urescaled', 'radii', 'radii0', ...
                'sphereCenter', 'sphereRadius', 'sphereParameters')
        else
            load(Ufn, 'U', 'Urescaled', 'radii', 'radii0', ...
                'sphereCenter', 'sphereRadius', 'sphereParameters')
        end

        % At this point we have performed mean curvature flow to the surface
        % according to "Can mean curvature flow be made non-singular?" 
        % [Kazhdan et al. 2012], giving a sphere of radius R centered at 
        % r0. We found the distance of the            
        % surface vertices from r0 and subtracted R, the mapped radial
        % distance from r0. This gives us a measure of how much the
        % surface is extended beyond R or retracted from R at each
        % point on the mesh. This defines the 'corrugation' field 
        % \deltar, which can be expressed either
        % as a function of position on the cell surface or as a function of
        % position on the mapped sphere found by mean curvature flow, with 
        % values \deltar(x) mapped to \deltar(.
        % Conveniently, scalar field patterns on the sphere can be compared
        % directly across samples. Therefore, we then decompose the
        % scalar field defined on the mapped (spherical) surface into 
        % spherical harmonics. Spherical harmonics are a set of 
        % functions which can reproduce the original field when multiplied by 
        % weights and summed together. More formally, they are 
        % eigenmodes of the Laplace-Beltrami operator defined on the sphere.
        % The inner product between a given eigenmode and the measured
        % corrugation field \delta r.

        %% Construct Mesh Laplace-Beltrami Operator ===============================

        % laplaceBeltrami = construct_laplace_beltrami( face, vertex );
        laplaceBeltrami = cotmatrix( Urescaled - sphereCenter, mesh.f );

        %% Find Eigenfunctions of Laplace-Beltrami Operator =======================
        tic
        disp('computing spherical harmonics')
        [V,~] = eigs(laplaceBeltrami,nModes,0);

        toc

        %% View Results ===========================================================
        close all ;
        figure('units', 'centimeters', 'position', [0,0,10,5])

        outfn = fullfile(outdir, ...
            [fn '_powerSpectrum_radialu.mat']) ;
        if ~exist(outfn, 'file') || overwrite 
            ff =  radii - radii0 ;
            titleStr = 'Laplace-Beltrami power spectrum of \deltar' ;
            % The 'Fourier Transform' of the signal
            rawPowers = V' * ff;
            % rawPowersNormV = rawPowers ./ length(mesh.v) ;

            % Plot Results ------------------------------------------------------------
            clf
            bar(abs(rawPowers), 'FaceColor',colors(1, :),'EdgeColor','none')
            xlabel('spherical harmonic index')
            ylabel('spectral power')
            sgtitle(titleStr);
            figfn = fullfile(outdir, ...
                [fn '_powerSpectrum_radialu.png']) ;
            saveas(gcf, figfn)


            %% Sort by \ell value
            llvals = [] ;
            lls = 0:30 ;
            dmyk = 1 ;
            powers = zeros(numel(lls),1) ;
            % powersNormV = powers ;
            for Lind = 1:numel(lls)
                ll = lls(Lind) ;
                powersLs = zeros(2*ll + 1, 1) ;
                for qq = 1:(2*ll + 1)
                    llvals(dmyk)= ll ;
                    powersLs(qq) = abs(rawPowers(dmyk)) ;
                    dmyk = dmyk + 1 ;
                end
                powers(Lind) = (sum(powersLs)) ;
            end

            % Plot results of power spectrum, indexed by l, the spherical
            % harmonic index
            clf
            bar(lls, powers, 'FaceColor',colors(1, :),'EdgeColor','none')
            xlabel('degree of roughness')
            ylabel('spectral power')
            sgtitle(titleStr)
            figfn = fullfile(outdir, ...
                 [fn '_powerSpectrumSummed_radialu.pdf']) ;
            saveas(gcf, figfn)
            
            % save results
            save(outfn, 'powers', ...
                'lls', 'llvals', 'rawPowers')
        else
            load(outfn, 'powers', 'lls', 'llvals', 'rawPowers')
        end

    else
        load(outfn2, 'powers', 'lls', 'llvals', 'rawPowers')
    end
end

%% Compare all surfaces
for ii = 1:length(fns)
    fn = strrep(fns(ii).name, '.ply', '') ;
    outfn = fullfile(outdir, ...
        [fn '_powerSpectrum_radialu.mat']) ;
    load(outfn, 'powers', 'llvals')

    if ii == 1
        powersAll = powers ;
    else
        powersAll = cat(2, powersAll, powers) ;
    end
end
clf
bar(lls, powersAll, 'edgecolor', 'none') 
xlabel('degree of roughness')
ylabel('spectral power')
title('Comparison across surfaces')
saveas(gcf, fullfile(outdir, 'comparison_of_powerSpectra_radialu.pdf'))

% save stats for these powers
save(fullfile(outdir, 'radialu_spectralPowers.mat'), ...
    'powersAll',  'lls', 'dirs')

%% Compare all experiment cases 
close all
colors = [ 31, 177, 3; ...
    110,110,110; ...
    203, 41,123] ./ 255.0 ;
load(fullfile(datadir, 'analysis', 'radialu_spectralPowers.mat'), ...
    'powersAll',  'lls')

meanPower = mean(powersAll, 2) ;
stdPower = std(powersAll, [], 2) ;
stePower = stdPower ./ sqrt(size(powersAll, 2)) ;

lineProps = {'-','color', colors(1, :)} ;
means = meanPower' ;
h =shadedErrorBar(lls, means, stdPower', ...
    'lineProps', lineProps, 'patchSaturation', 0.2) ;
hold on;

ylim([0, Inf])

% Now add stes (standard errors on the mean)
lineProps = {'-','color', colors(1, :)} ;
errorbar(lls, meanPower, stePower, '.', 'color', colors(1, :))

xlabel('degree of roughness')
ylabel('spectral power')

title('Comparison across conditions')
saveas(gcf, 'comparison_of_powerSpectra_radialu.pdf')

xlim([10, 20])
axis square
    saveas(gcf, 'comparison_of_powerSpectra_radialu_zoom2.pdf')
close all
