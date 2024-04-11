%% This workflow is based on code from:
% N. P. Mitchell*, D. J. Cislo*, TubULAR: Tracking in toto deformations of
% dynamic tissues via constrained maps. [biorxiv]
% but is here adapted to use spherical harmonics

%% Navigate to the directory containing a folder for each condition.
% cd /path/to/PLY_Example_Data/

%% Prepare pathnames
close all
clc
clearvars

% First change directories to where PLY meshes of cell clusters exist
datadir = pwd ;

% Running this line in the script grabs the path of this script
mfilePath = mfilename('fullpath');
if contains(mfilePath,'LiveEditorEvaluationHelper')
    mfilePath = matlab.desktop.editor.getActiveFilename;
end
[ scriptDir, scriptName ] = fileparts(mfilePath) ; 
cd(scriptDir)

%%
% Define the subdirectories where different conditions are stored
dirs = {'knockdown', 'wildtype', 'overexpression' } ;

% Define the shorthand variables describing these conditions, in the same
% order as the previous line.
shorthand = {'KD', 'WT', 'OE'} ;

% Resolution of the PLYs, in um / voxel width. Make sure the PLYs are saved
% with this resolution, and in isotropic resolution (same for x, y, and z).
% In this script, we load a saved pixel resolution for each mesh if
% one exists on disk, but if not, then we use this global scale.
pix2um_global = 0.0302 ;  % um / pixel


% Add DECLab to the path: https://github.com/DillonCislo/DECLab
% This is a submodule of the SeptinManuscriptData repo.
addpath(genpath('../DECLab/'))

% Also download and add gptoolbox to the path.
% This is already included in the SeptinManuscriptData repo.
% https://github.com/alecjacobson/gptoolbox
addpath(genpath('../gptoolbox/mesh/'))

% Change back to data directory
cd(datadir)


        
%% Default Options
% overwrite previous results if on disk
overwrite = false ;

% The number of Laplacian eigenvectors to calculate
nModes = 1000;

% Here we decide what measure of surface topography to decompose into
% spherical harmonics. The options are 'radialu',  'HH', and 'dist'. 
% We use radialu, which is the radial distance of the mesh vertices from 
% their counterparts on the spherical surface to which they are mapped.
signalType = 'radialu' ;

%% Iterate over each condition (each condition is stored in a unique dir)
for jj = 1:length(dirs)
    fns = dir(fullfile(dirs{jj}, '*.ply')) ;

    outdir = fullfile(dirs{jj}, 'analysis') ;

    if ~exist(outdir, 'dir')
        mkdir(outdir)
    end
    
    %% Iterate over each PLY file surface
    for ii = 1:length(fns)
        fn = strrep(fns(ii).name, '.ply', '') ;
        disp(['ii = ' num2str(ii) ': ' fn])

        outfn1 = fullfile(outdir, ...
            [fn '_powerSpectrum_' signalType '.mat']) ;
        if ~exist(outfn1, 'file') || overwrite 

            mesh = read_ply_mod(fullfile(fns(ii).folder, [fn '.ply'])) ;

            % Identify an experiment-specific pixel resolution if saved to
            % disk
            if exist(fullfile(fns(ii).folder, [fn '_resolution.txt']), 'file')
                pix2um = dlmread([fn '_resolution.txt']) ;
            else
                pix2um = pix2um_global ;
            end

            % Convert mesh to microns
            mesh.v = mesh.v * pix2um ;
            
        
            % Conformally map to the unit sphere
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
                clf
                subplot(2, 2, 1)
                trisurf(triangulation(mesh.f, mesh.v), 'facecolor', [0., 0., 0], ...
                 'edgecolor', 'none', 'facealpha', 0.2); 
                hold on;
                trisurf(triangulation(mesh.f, Urescaled), 'facecolor', [0.5, 0.5, 0], ...
                 'edgecolor', 'none', 'facealpha', 0.2); 
                hold on; 

                camlight
                axis equal ;
                grid off ;
                trisurf(triangulation(mesh.f, mesh.v), 'facecolor', [0.2, 0.2, 1],...
                 'edgecolor', 'none', 'facealpha', 0.2)
                axis equal; 
                grid off ; % axis off ;
                subplot(2, 2, 2)
                trisurf(triangulation(mesh.f, Urescaled), ...
                 'edgecolor', 'none', 'facevertexCdata', radii)
                axis equal ;
                grid off ;
                title('radial position of surface');
                clims = caxis ;

                subplot(2, 2, 3)
                trisurf(triangulation(mesh.f, Urescaled), ...
                 'edgecolor', 'none', 'facevertexCdata', radii0)
                axis equal ;
                grid off
                title('radius of mapped surface');
                caxis(clims)
                
                set(gcf, 'color', 'w')
                figfn = fullfile(outdir,[ fn '_sphericalFit.png']) ;
                saveas(gcf, figfn)

                %% save the result
                save(Ufn, 'U', 'Urescaled', 'radii', 'radii0', ...
                    'sphereCenter', 'sphereRadius', 'sphereParameters')
            else
                load(Ufn, 'U', 'Urescaled', 'radii', 'radii0', ...
                    'sphereCenter', 'sphereRadius', 'sphereParameters')
            end

            Hfn = fullfile(outdir, [fn '_meanCurvature.mat']) ;
            if ~exist(Hfn, 'file') || overwrite 
                DEC = DiscreteExteriorCalculus(mesh.f, mesh.v) ;
                lapx = DEC.laplacian(mesh.v(:, 1)) ;
                lapy = DEC.laplacian(mesh.v(:, 2)) ;
                lapz = DEC.laplacian(mesh.v(:, 3)) ;
                H3d = sum(mesh.vn .* [lapx, lapy, lapz], 2) * 0.5 ;
                areas = 0.5 * doublearea(mesh.v, mesh.f) ;
                [V2F, F2V] = meshAveragingOperators(mesh.f, mesh.v) ;
                H3d_faces = (V2F * H3d) ;
                disp('saving mean curvature')
                save(Hfn, 'H3d', 'areas', 'H3d_faces')
            else
                load(Hfn, 'H3d')
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

            %% Construct Mesh Laplace-Beltrami Operator ===================
            laplaceBeltrami = cotmatrix( Urescaled - sphereCenter, mesh.f );

            %% Find Eigenfunctions of Laplace-Beltrami Operator ===========
            tic
            disp('computing spherical harmonics')
            [V,~] = eigs(laplaceBeltrami,nModes,0);

            toc

            %% View Results ===============================================
            close all ;
            figure('units', 'centimeters', 'position', [0,0,10,5])
            

            outfn = fullfile(outdir, ...
                [fn '_powerSpectrum_' signalType '.mat']) ;
            if ~exist(outfn, 'file') || overwrite 
                if strcmpi(signalType, 'HH')
                    ff = H3d ;
                    titleStr = 'Laplace-Beltrami power spectrum of H' ;
                elseif strcmpi(signalType, 'radialu')
                    ff =  radii - radii0 ;
                    titleStr = 'Laplace-Beltrami power spectrum of \deltar' ;
                elseif strcmpi(signalType, 'dist')
                    ff = vecnorm(mesh.v - Urescaled, 2, 2) ;
                    titleStr = 'Laplace-Beltrami power spectrum of \deltax' ;
                end
                % The 'Fourier Transform' of the signal
                rawPowers = V' * ff;
                % rawPowersNormV = rawPowers ./ length(mesh.v) ;

                % Plot Results ----------------------------------------
                % plot( abs(x), '.-', 'LineWidth', 0.5, 'MarkerFaceColor', 'b' );
                clf
                % subplot(1, 2, 1)
                colors = [
                    0.90    0.55    0.55
                    0.5     0.5     0.5
                    0.62    0.76    0.84 
                    ];
                
                bar(abs(rawPowers), 'FaceColor',colors(1, :),'EdgeColor','none')
                xlabel('spherical harmonic index (spectral shape index)')
                ylabel('spectral weight')

                sgtitle(titleStr);

                figfn = fullfile(outdir, ...
                    [fn '_powerSpectrum_' signalType '.png']) ;
                saveas(gcf, figfn)


                %% Sort by \ell value (which indexes the spatial frequency)
                llvals = [] ;
                lls = 0:30 ;
                dmyk = 1 ;
                powers = zeros(numel(lls),1) ;

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

                clf
                % subplot(1, 2, 1)
                bar(lls, powers, 'FaceColor',colors(1, :),'EdgeColor','none')
                xlabel('spectral shape index')
                ylabel('spectral weight')

                sgtitle(titleStr)
                figfn = fullfile(outdir, ...
                     [fn '_powerSpectrumSummed_' signalType '.pdf']) ;
                saveas(gcf, figfn)
                % save results
                save(outfn, 'powers', ...
                    'lls', 'llvals', 'rawPowers')
            else
                load(outfn, 'powers', 'lls', 'llvals', 'rawPowers')

            end
        
        else
            load(outfn1, 'powers', 'lls', 'llvals', 'rawPowers')
            
        end
    end

    %% Compare all surfaces for this batch of PLYs
    for ii = 1:length(fns)
        fn = strrep(fns(ii).name, '.ply', '') ;
        outfn = fullfile(outdir, ...
            [fn '_powerSpectrum_' signalType '.mat']) ;
        load(outfn, 'powers', 'llvals')

        if ii == 1
            powersAll = powers ;
        else
            powersAll = cat(2, powersAll, powers) ;
        end
    end
    clf
    bar(lls, powersAll, 'edgecolor', 'none') 
    xlabel('spherical harmonic index (spectral shape index)')
    ylabel('spectral weight, A_l')
    title('Comparison across surfaces')
    saveas(gcf, fullfile(outdir, ...
        ['comparison_of_powerSpectra_' signalType '.pdf']))

    % save stats for these powers
    save(fullfile(outdir, [signalType '_spectralPowers.mat']), ...
        'powersAll',  'lls', 'dirs')
    
    disp('done with this batch of PLYs')
end

pause(1)

%% Compare all experiment cases 
close all

colors = [
    0.90    0.55    0.55
    0.5     0.5     0.5
    0.62    0.76    0.84 
    ];

% Alternative colors:
% colors = [ 31, 177, 3; ...
%     110,110,110; ...
%     203, 41,123] ./ 255.0 ;

clc
meanPowers = {} ;
stdPowers = {} ;
stePowers = {} ;


% First obtain the normalization. Here, we choose to normalize by the peak
% power in the wildtype 
wtIndex = find(contains(dirs,'wildtype'));
outdir = dirs{wtIndex} ;
load(fullfile(outdir, 'analysis', [signalType '_spectralPowers.mat']), ...
    'powersAll',  'lls')
meanPower = mean(powersAll, 2) ;
normalization = max(meanPower(:)) ;


hs = {} ;
for pp = 1:length(dirs)

    outdir = dirs{pp} ;
    load(fullfile(outdir, 'analysis', [signalType '_spectralPowers.mat']), ...
        'powersAll',  'lls')


    meanPower = mean(powersAll, 2) ;
    stdPower = std(powersAll, [], 2) ;
    stePower = stdPower ./ sqrt(size(powersAll, 2)) ;

    if pp == 1
        lls0 = lls ;
    else
        assert(all(lls == lls0)) ;
    end
    meanPowers{pp} = meanPower ;
    stdPowers{pp} = stdPower ;
    stePowers{pp}= stePower ;

    lineProps = {'-','color', colors(pp, :)} ;
    means = meanPower' ;
    h =shadedErrorBar(lls, means / normalization, ...
        stdPower' / normalization, ...
        'lineProps', lineProps, 'patchSaturation', 0.2) ;
    hold on;

    ylim([0, Inf])
end

legend(strrep(dirs, '_', ' '),'AutoUpdate','off')

% Now add stes (standard error on the mean) 
for pp = 1:length(dirs)
    lineProps = {'-','color', colors(pp, :)} ;
    errorbar(lls, meanPowers{pp}/ normalization,...
        stePowers{pp}/ normalization, '.', 'color', colors(pp, :))
end

xlabel('spectral shape index')
ylabel('relative spectral weight')

title('Comparison across conditions')
set(gcf, 'Units', 'centimeters')
set(gcf, 'Position', [0,0,8,8])
saveas(gcf, ['comparison_of_powerSpectra_' signalType '.pdf'])

xlim([10, 20])
ylim([0,.13]) 
axis square
    saveas(gcf, ['comparison_of_powerSpectra_' signalType '_zoom.pdf'])
close all

% Save data
save(['statistics_' signalType '.mat'], ...
   'meanPowers', 'stdPowers', 'stePowers', 'dirs', 'lls')


%% Bin powers by their l value (which determines spatial scale of variation)
load(['statistics_' signalType '.mat'], ...
               'meanPowers', 'stdPowers', 'stePowers', 'dirs', 'lls')
ind = 2 ;

%% low mode
kd = meanPowers{1}(ind) ;
wt = meanPowers{2}(ind) ;
oe = meanPowers{3}(ind) ;
kd_unc = stePowers{1}(ind) ;
wt_unc = stePowers{2}(ind) ;
oe_unc = stePowers{3}(ind) ;

num = -kd+wt ;
denom = sqrt(kd_unc.^2 + wt_unc.^2) ;
zscore_kdwt = num / denom ;
pval_kdwt = normcdf(zscore_kdwt);

num = -kd+oe ;
denom = sqrt(kd_unc.^2 + oe_unc.^2) ;
zscore_kdoe = num / denom ;
pval_kdoe = normcdf(zscore_kdoe);

num = -wt+oe ;
denom = sqrt(oe_unc.^2 + wt_unc.^2) ;
zscore_wtoe = num / denom ;
pval_wtoe = normcdf(zscore_wtoe);

% superbar plot
hf = figure('Position', [100 100 400 400], 'units', 'centimeters');
clf;

Colors = [
    0.90    0.55    0.55
    0.5     0.5     0.5
    0.62    0.76    0.84
    ];

P = [0,  pval_kdwt, pval_kdoe; ...
    pval_kdwt, 0, pval_wtoe;...
    pval_kdoe, pval_wtoe, 0];

E = [kd_unc/wt, wt_unc/wt, oe_unc/wt] ;
superbar([1,2,3],[kd/wt,wt/wt,oe/wt], 'E', E, 'P', P,...
    'BarFaceColor', Colors) ;

xticks([1,2,3])
categories = {'KD', 'WT', 'OE'} ;
xticklabels(categories)
ylabel(sprintf('Relative spectral weight in mode l=%d', lls(ind)))
set(gcf, 'Units', 'centimeters')
set(gcf, 'Position', [0,0,8,8])
saveas(gcf, sprintf('low_order_mode_comparison_l%02d.pdf', lls(ind)))
means = [kd,wt,oe] ;
save(sprintf('low_order_mode_comparison_l%02.mat', lls(ind)), ...
    'P', 'E', 'means', 'categories')
    
%% high order modes
upperThres = 25 ;
close all
for thres = 6

    kd = sum(meanPowers{1}(lls > thres & lls < upperThres, :)) ;
    wt = sum(meanPowers{2}(lls > thres & lls < upperThres, :)) ;
    oe = sum(meanPowers{3}(lls > thres & lls < upperThres, :)) ;
    kd_unc = sqrt(sum((stePowers{1}(lls > thres & lls < upperThres, :)).^2)) ;
    wt_unc = sqrt(sum((stePowers{2}(lls > thres & lls < upperThres, :)).^2)) ;
    oe_unc = sqrt(sum((stePowers{3}(lls > thres & lls < upperThres, :)).^2)) ;

    % pvalue
    num = -abs(kd-wt) ;
    denom = sqrt(kd_unc.^2 + wt_unc.^2) ;
    zscore_kdwt = num / denom ;
    pval_kdwt = normcdf(zscore_kdwt);

    num = -abs(kd-oe) ;
    denom = sqrt(kd_unc.^2 + oe_unc.^2) ;
    zscore_kdoe = num / denom ;
    pval_kdoe = normcdf(zscore_kdoe);

    num = -abs(wt-oe) ;
    denom = sqrt(oe_unc.^2 + wt_unc.^2) ;
    zscore_wtoe = num / denom ;
    pval_wtoe = normcdf(zscore_wtoe);

    %% superbar plot
    hf = figure('Position', [100 100 400 400], 'units', 'centimeters');
    clf;

    Colors = [
        0.90    0.55    0.55
        0.5     0.5     0.5
        0.62    0.76    0.84
        ];
    
    P = [0,  pval_kdwt, pval_kdoe; ...
        pval_kdwt, 0, pval_wtoe;...
        pval_kdoe, pval_wtoe, 0];
    
    E = [kd_unc/wt, wt_unc/wt, oe_unc/wt] ;
    superbar([1,2,3],[kd/wt,wt/wt,oe/wt], 'E', E, 'P', P,...
        'BarFaceColor', Colors) ;
    
    xticks([1,2,3])
    categories = {'KD', 'WT', 'OE'} ;
    xticklabels(categories)
    ylabel(sprintf('Relative spectral weight in modes l>%d', thres))

    set(gcf, 'Units', 'centimeters')
    set(gcf, 'Position', [0,0,8,8])
    outfn = sprintf('high_order_mode_comparison_lgt%02d_lt%02d.pdf', thres, upperThres) ;
    saveas(gcf, outfn)
    
    % Save the data reflected in the plot
    outfn = sprintf('high_order_mode_comparison_lgt%02d_lt%02d.mat', thres, upperThres) ;
    save(outfn, 'P', 'E', 'means', 'categories')
    
end
