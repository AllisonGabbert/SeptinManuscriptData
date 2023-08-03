%% This workflow is based on code from:
% N. P. Mitchell*, D. J. Cislo*, TubULAR: Tracking in toto deformations of
% dynamic tissues via constrained maps. [biorxiv]
% but is here adapted to use spherical harmonics

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
dirs = {'septin_knockdown', 'wildtype_control', 'septin_overexpression' } ;

% Define the shorthand variables describing these conditions, in the same
% order as the previous line.
shorthand = {'KD', 'WT', 'OE'} ;

% Resolution of the PLYs, in um / voxel width. Make sure the PLYs are saved
% with this resolution, and in isotropic resolution (same for x, y, and z)
pix2um = 1.0 ;  % um / pixel


% Add DECLab to the path: https://github.com/DillonCislo/DECLab
% This is a submodule of the SeptinManuscriptData repo.
addpath(genpath('../../DECLab/'))

% Also download and add gptoolbox to the path.
% This is already included in the SeptinManuscriptData repo.
% https://github.com/alecjacobson/gptoolbox
addpath(genpath('../gptoolbox/mesh/'))

% Change back to data directory
cd(datadir)

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

%% To use the plotting features here, you must add shadedErrorBar to path
% https://www.mathworks.com/matlabcentral/fileexchange/26311-raacampbell-shadederrorbar 
addpath ../mnt/data/code/gut_matlab/plotting/shadedErrorBar/
        
%% Default Options
% overwrite previous results if on disk
overwrite = true ;

% The number of Laplacian eigenvectors to calculate
nModes = 1000;
signalTypes = {'radialu'} ; % 'HH', 'dist'}

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
            [fn '_powerSpectrum_' signalTypes{1} '.mat']) ;
        outfn2 = fullfile(outdir, ...
            [fn '_powerSpectrum_' signalTypes{2} '.mat']) ;
        if ~exist(outfn1, 'file') || ~exist(outfn2, 'file') || overwrite 

            mesh = read_ply_mod(fullfile(fns(ii).folder, [fn '.ply'])) ;
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
                H3d = sum(mesh.vn .* DEC.laplacian(mesh.v), 2) * 0.5 ;
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
            % r0. We found the distance of the            % surface vertices from r0 and subtracted R, the mapped radial
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
            % DEC = DiscreteExteriorCalculus( face, vertex );

            %% Find Eigenfunctions of Laplace-Beltrami Operator =======================
            tic
            disp('computing spherical harmonics')
            [V,~] = eigs(laplaceBeltrami,nModes,0);

            toc

            %% View Results ===========================================================
            close all ;
            figure('units', 'centimeters', 'position', [0,0,10,5])
            for kk = 1:length(signalTypes)
                
                outfn = fullfile(outdir, ...
                    [fn '_powerSpectrum_' signalTypes{kk} '.mat']) ;
                if ~exist(outfn, 'file') || overwrite 
                    if strcmpi(signalTypes{kk}, 'HH')
                        ff = H3d ;
                        titleStr = 'Laplace-Beltrami power spectrum of H' ;
                    elseif strcmpi(signalTypes{kk}, 'radialu')
                        ff =  radii - radii0 ;
                        titleStr = 'Laplace-Beltrami power spectrum of \deltar' ;
                    elseif strcmpi(signalTypes{kk}, 'dist')
                        ff = vecnorm(mesh.v - Urescaled, 2, 2) ;
                        titleStr = 'Laplace-Beltrami power spectrum of \deltax' ;
                    end
                    % The 'Fourier Transform' of the signal
                    rawPowers = V' * ff;
                    % rawPowersNormV = rawPowers ./ length(mesh.v) ;

                    % Plot Results ------------------------------------------------------------
                    % plot( abs(x), '.-', 'LineWidth', 0.5, 'MarkerFaceColor', 'b' );
                    clf
                    % subplot(1, 2, 1)
                    bar(abs(rawPowers), 'FaceColor',colors(1, :),'EdgeColor','none')
                    xlabel('spherical harmonic index')
                    ylabel('spectral power')
                    
                    sgtitle(titleStr);
                    
                    figfn = fullfile(outdir, ...
                        [fn '_powerSpectrum_' signalTypes{kk} '.png']) ;
                    saveas(gcf, figfn)


                    %% Sort by \ell value
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
                    xlabel('degree of roughness')
                    ylabel('spectral power')
                    
                    sgtitle(titleStr)
                    figfn = fullfile(outdir, ...
                         [fn '_powerSpectrumSummed_' signalTypes{kk} '.pdf']) ;
                    saveas(gcf, figfn)
                    % save results
                    save(outfn, 'powers', ...
                        'lls', 'llvals', 'rawPowers')
                else
                    load(outfn, 'powers', 'lls', 'llvals', 'rawPowers')

                end
            end
        
        else
            load(outfn2, 'powers', 'lls', 'llvals', 'rawPowers')
            
        end
    end

    %% Compare all surfaces for this batch of PLYs
    for qq = 1:length(signalTypes)
        signalType = signalTypes{qq} ;
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
        xlabel('degree of roughness')
        ylabel('spectral power')
        title('Comparison across surfaces')
        saveas(gcf, fullfile(outdir, ...
            ['comparison_of_powerSpectra_' signalType '.pdf']))
    
        % save stats for these powers
        save(fullfile(outdir, [signalType '_spectralPowers.mat']), ...
            'powersAll',  'lls', 'dirs')
    end
end


%% Compare all experiment cases 
close all

colors = [ 31, 177, 3; ...
    110,110,110; ...
    203, 41,123] ./ 255.0 ;

clc
meanPowers = {} ;
stdPowers = {} ;
stePowers = {} ;

for qq = 1:length(signalTypes)
    signalType = signalTypes{qq} ;

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
        h =shadedErrorBar(lls, means, stdPower', ...
            'lineProps', lineProps, 'patchSaturation', 0.2) ;
        hold on;

        ylim([0, Inf])
    end

    legend(strrep(dirs, '_', ' '),'AutoUpdate','off')

    % Now add stes (standard error on the mean) 
    for pp = 1:length(dirs)
        lineProps = {'-','color', colors(pp, :)} ;
        errorbar(lls, meanPowers{pp}, stePowers{pp}, '.', 'color', colors(pp, :))
    end

    xlabel('degree of roughness')
    ylabel('spectral power')

    title('Comparison across conditions')
    saveas(gcf, ['comparison_of_powerSpectra_' signalType '.pdf'])

    xlim([10, 20])
    % ylim([200,1200]) 
    axis square
        saveas(gcf, ['comparison_of_powerSpectra_' signalType '_zoom2.pdf'])
   close all

   %% Save data
   save(['statistics_' signalTypes{qq} '.mat'], ...
       'meanPowers', 'stdPowers', 'stePowers', 'dirs', 'lls')

end


%% Bin powers by their l value (which determines spatial scale of variation)
qq = 2 ;
load(['statistics_' signalTypes{qq} '.mat'], ...
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
    E = cat(3, kd_unc, wt_unc, oe_unc) ;

    colors = [ 31, 177, 3; ...
        110,110,110; ...
        203, 41,123] ./ 255.0 ;
    
    P = [0,  pval_kdwt, pval_kdoe; ...
        pval_kdwt, 0, pval_wtoe;...
        pval_kdoe, pval_wtoe, 0];
    
    E = [kd_unc, wt_unc, oe_unc] ;
    superbar([1,2,3],[kd,wt,oe], 'E', E, 'P', P,...
        'BarFaceColor', Colors) ;
    
    xticks([1,2,3])
    categories = {'KD', 'WT', 'OE'} ;
    xticklabels(categories)
    ylabel(sprintf('Spectral power in mode l=%d', lls(ind)))
    saveas(gcf, sprintf('low_order_mode_comparison_l%02d.pdf', lls(ind)))
    means = [kd,wt,oe] ;
    save(sprintf('low_oder_mode_comparison_l%02.mat', lls(ind)), ...
        'P', 'E', 'means', 'categories')
    
%% high order modes
upperThres = 2 ;
close all
for thres = 0
    % thres  = 4 ;

    kd = sum(meanPowers{1}(lls > thres & lls < upperThres, :)) ;
    wt = sum(meanPowers{2}(lls > thres & lls < upperThres, :)) ;
    oe = sum(meanPowers{3}(lls > thres & lls < upperThres, :)) ;
    kd_unc = sqrt(sum((stePowers{1}(lls > thres & lls < upperThres, :)).^2)) ;
    wt_unc = sqrt(sum((stePowers{2}(lls > thres & lls < upperThres, :)).^2)) ;
    oe_unc = sqrt(sum((stePowers{3}(lls > thres & lls < upperThres, :)).^2)) ;

    % pvalue
    num = kd-wt ;
    denom = sqrt(kd_unc.^2 + wt_unc.^2) ;
    zscore_kdwt = num / denom ;
    pval_kdwt = normcdf(zscore_kdwt);

    num = kd-oe ;
    denom = sqrt(kd_unc.^2 + oe_unc.^2) ;
    zscore_kdoe = num / denom ;
    pval_kdoe = normcdf(zscore_kdoe);

    num = wt-oe ;
    denom = sqrt(oe_unc.^2 + wt_unc.^2) ;
    zscore_wtoe = num / denom ;
    pval_wtoe = normcdf(zscore_wtoe);

    %% superbar plot
    hf = figure('Position', [100 100 400 400], 'units', 'centimeters');
    clf;
    E = cat(3, kd_unc, wt_unc, oe_unc) ;

    Colors = [
        0.90    0.55    0.55
        0.5     0.5     0.5
        0.62    0.76    0.84
        0.89    0.10    0.11
        0.12    0.47    0.70
        ];

    P = [0,  pval_kdwt, pval_kdoe; ...
        pval_kdwt, 0, pval_wtoe;...
        pval_kdoe, pval_wtoe, 0];
    
    E = [kd_unc, wt_unc, oe_unc] ;
    superbar([1,2,3],[kd,wt,oe], 'E', E, 'P', P,...
        'BarFaceColor', Colors) ;
    
    xticks([1,2,3])
    categories = {'KD', 'WT', 'OE'} ;
    xticklabels(categories)
    ylabel(sprintf('Spectral power in modes l>%d', thres))
    outfn = sprintf('high_order_mode_comparison_lgt%02d_lt%02d.pdf', thres, upperThres) ;
    saveas(gcf, outfn)
    
    % Save the data reflected in the plot
    outfn = sprintf('high_order_mode_comparison_lgt%02d_lt%02d.mat', thres, upperThres) ;
    save(outfn, 'P', 'E', 'means', 'categories')
    
end