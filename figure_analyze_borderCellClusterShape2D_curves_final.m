%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% First add raw MIP images to ../data/masked_septin_images/
% Also add resampled aligned profiles of border cell clusters to
% ../data/contours
%
% Then run the code below.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Check we get the same thing from Adele's data & do Hausdorff distance
close all; clear all

origpath = matlab.desktop.editor.getActiveFilename;
cd(fileparts(origpath))
addpath(genpath(fullfile(fileparts(origpath), 'external')))
cd('../')

% Change these paths to where the contour data is stored
curvdir = fullfile(pwd, 'data/contours/') ;
outdir = fullfile(pwd, 'figures/') ;
if ~exist(outdir, 'dir')
    mkdir(outdir)
end
rootdir = fullfile(pwd, 'data/masked_septin_images/') ;
dirs = {'SeptinKnockdown', 'Control', 'SeptinOverexpression'} ;
categories = {'KD', 'WT', 'OE'} ;
metrics = {'perim', 'area', 'shape'} ;


%% Get metrics from raw binary images

pix2um = 1 / 10.0974 ; % resolution in um / pix

dataFn = fullfile(outdir, ['2Dshapes_statistics.mat']) ;

if exist(dataFn, 'file')
    disp('Loading data from disk...')
    load(dataFn, ...
        'perims_avg', 'perims_std', 'perims_ste', ...
        'areas_avg', 'areas_std', 'areas_ste', 'shapes_avg', ...
        'shapes_std', 'shapes_ste', ...
        'perims', 'areas', 'shapes', 'Centroids')
else
    disp('computing...')
    shapes = {} ;
    areas = {} ;
    perims = {} ;
    Centroids = {} ;
    for ii = 1:length(dirs) 

        fns = dir(fullfile(rootdir, dirs{ii}, 'binary images', 'Threshold*.tif')) ;

        nExpts = length(fns) ;
        area = zeros(1, nExpts) ;
        perim = zeros(1, nExpts) ;
        nprotrusions = zeros(1, nExpts) ;
        curves = cell(nExpts, 1) ;
        centroids = zeros(nExpts, 2) ;

        for jj = 1:nExpts
            disp(['jj = ' num2str(jj) '/' num2str(nExpts)])
            mask = imread(fullfile(fns(jj).folder, fns(jj).name)) ;

            % check that all images are the same resolution so we can compare
            % them
            assert(all(size(mask)== 512))        

            % check that the outline is right
            props = regionprops(mask, "Area", "Perimeter", "Centroid")  ;
            [~, clusterInd] = max([props.Area]) ;

            tmp = bwboundaries(mask) ;
            curv = tmp{clusterInd} ; % I think the index should be 1, but let's just confirm by plotting it
            % plot(curv(:, 1), curv(:, 2), '.-') ; axis equal; pause(1)  % Note: we're just inspecting this contour. Should look like a cluster
            % axis equal ;
            % hold on
            curves{jj} = curv * pix2um ; 
            centroids(jj, :) = props(clusterInd).Centroid * pix2um ;

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Count protrusions
            radius = vecnorm(curv - props(clusterInd).Centroid, 2, 2) ;
            radius2 = [radius; radius] ; % repeat to avoid edge effects
            % radsm = savgol(double(radius2), double(5), ...
            %     double(round(0.3 * length(radius)))) ; % smooth by 10% perimeter
            % [pks,locs,w,p] = findpeaks(radsm) ;
            % what is radius of a circle with the same area?
            rad_thres = sqrt(props(clusterInd).Area / pi ) ; % threshold radius
            pks = peakfinder(radius2, rad_thres*0.3, rad_thres, 1, false, false) ;
            npts = size(curv, 1) ;
            pks = pks(pks > npts * 0.5 & pks < npts * 1.5) ;
            pks = mod(pks, npts) ;
            pks(pks == 0) = npts ;
            nprotrusions(jj) = length(pks) ;

            %scatter(curv(pks, 1), curv(pks, 2), 50, 'r') ; axis equal; pause(1)  % Note: we're just inspecting this contour. Should look like a cluster
            set(gcf, 'visible', 'off')
            imshow(mask) ;
            hold on
            scatter(curv(pks, 2), curv(pks, 1), 350, 'r', 'filled')
            pause(0.01)
            if ~exist(fullfile(fns(jj).folder, 'protrusionIms'), 'dir')
                mkdir(fullfile(fns(jj).folder, 'protrusionIms')) 
            end
            saveas(gcf, fullfile(fns(jj).folder, 'protrusionIms', fns(jj).name))
            clf

            % store perimeter and area
            area(jj) =  props(clusterInd).Area * pix2um^2;
            perim(jj) = props(clusterInd).Perimeter * pix2um ;
        end
        areas{ii} = area ;
        perims{ii} = perim ;
        shapes{ii} = perim ./ sqrt(area)  ;
        Centroids{ii} = centroids ; 
    end


    %% Statistics
    perims_avg = zeros(1, 3) ;
    perims_std = zeros(1, 3) ;
    perims_ste = zeros(1, 3) ;
    areas_avg = zeros(1, 3) ;
    areas_std = zeros(1, 3) ;
    areas_ste = zeros(1, 3) ;
    shapes_avg = zeros(1, 3) ;
    shapes_std = zeros(1, 3) ;
    shapes_ste = zeros(1, 3) ;
    for ii = 1:3
        perims_avg(ii) = mean(perims{ii}) ;
        perims_std(ii) = std(perims{ii}) ;
        perims_ste(ii) = std(perims{ii}) ./ sqrt(length(perims{ii})) ;
        areas_avg(ii) = mean(areas{ii}) ;
        areas_std(ii) = std(areas{ii}) ;
        areas_ste(ii) = std(areas{ii}) ./ sqrt(length(areas{ii})) ;
        shapes_avg(ii) = mean(shapes{ii}) ;
        shapes_std(ii) = std(shapes{ii}) ;
        shapes_ste(ii) = std(shapes{ii}) ./ sqrt(length(shapes{ii})) ;
    end

    save(dataFn, ...
        'perims_avg', 'perims_std', 'perims_ste', ...
        'areas_avg', 'areas_std', 'areas_ste', 'shapes_avg', ...
        'shapes_std', 'shapes_ste',...
        'perims', 'areas', 'shapes', 'Centroids')
    close all
end
disp('done')

%% Plot area, perim, and shape indices

% smooth plot
profile_step = [0.03, 0.03, 0.03] ;

faceColors = [31, 177, 3; ...
    110,110,110; ...
    203, 41,123] ./ 255.0 ;

% smooth violin plot of area, perim, shape
for qq = 1:length(metrics)
    maxy = 0 ;
    miny = Inf ;
    maxys = [0,0,0] ;
    close all
    figure('units', 'centimeters', 'position', [0, 0, 10, 10])
    
    aas = {} ;
    for ii = 1:length(categories)
        if qq == 1
            aa = perims{ii};
        elseif qq==2
            aa = areas{ii} ;
        elseif qq == 3
            aa = shapes{ii} ;
        end
        
        aas{ii} = aa ;
        
        violin(aa(:), 'x', [ii,ii+1],...
            'facecolor',faceColors(ii,:),'edgecolor','none', 'facealpha', 1, ...
            'bw',profile_step(qq)*(max(aa) - min(aa)), 'mc','k-','medc', ...
            '', 'plotlegend', false)
        
        hold on;
        maxys(ii) = max(aa(:)) ;
        miny = min(miny, min(aa(:))) ;
        maxy = max(maxy, max(aa(:))) ;
    end
    xlim([0, 4])
    yLow = miny * 0.9 ;
    yHigh = maxy * 1.1 ;
    topY = yHigh * 1.1 ;
    ylim([yLow, yHigh])
    xticks([1,2,3])
    set(gca, "xtickLabels", categories)
    if qq == 1
        title('Border cell cluster perimeter')
        ylabel(['perimeter, P [' char(181) 'm]'])
        saveas(gcf, fullfile(outdir, '2Dshapes_perim.pdf'))
    elseif qq == 2
        title('Border cell cluster areas')
        ylabel(['area, A [' char(181) 'm^2]'])
        saveas(gcf, fullfile(outdir, '2Dshapes_area.pdf'))
    elseif qq == 3
        title('Border cell cluster shape')
        ylabel(['shape index, P/\surd A'])
        saveas(gcf, fullfile(outdir, '2Dshapes_shape.pdf'))
    end
    % 
    % 
    % Add p value and save plot again
    % Add p value and save plot again
    if qq == 1    
        wt = perims_avg(2);
        kd = perims_avg(1);
        oe = perims_avg(3);
        wt_unc = perims_ste(2) ;
        kd_unc = perims_ste(1) ;
        oe_unc = perims_ste(3) ;
    elseif qq==2
        wt = areas_avg(2) ;
        kd = areas_avg(1) ;
        oe = areas_avg(3) ;
        wt_unc = areas_ste(2) ;
        kd_unc = areas_ste(1) ;
        oe_unc = areas_ste(3) ;    
    elseif qq == 3
        wt = shapes_avg(2) ;
        kd = shapes_avg(1) ;
        oe = shapes_avg(3) ;
        wt_unc = shapes_ste(2) ;
        kd_unc = shapes_ste(1) ;
        oe_unc = shapes_ste(3) ;
    else
        error('here')
    end
    % 
    % num = -abs(kd-wt) ;
    % denom = sqrt(kd_unc.^2 + wt_unc.^2) ;
    % zscore_kdwt = num / denom ;
    % pval_kdwt = normcdf(zscore_kdwt);
    % 
    % num = -abs(kd-oe) ;
    % denom = sqrt(kd_unc.^2 + oe_unc.^2) ;
    % zscore_kdoe = num / denom ;
    % pval_kdoe = normcdf(zscore_kdoe);
    % 
    % num = -abs(wt-oe) ;
    % denom = sqrt(oe_unc.^2 + wt_unc.^2) ;
    % zscore_wtoe = num / denom ;
    % pval_wtoe = normcdf(zscore_wtoe);

    %if qq == 1
    % preview it:
    % clf; ecdf(aas{1}); hold on; ecdf(aas{2})
    [~,pval_kdwt_1] = kstest2(aas{2}, aas{1}, 'Tail', 'larger') ;
    [~,pval_kdwt_2] = kstest2(aas{1}, aas{2}, 'Tail', 'larger') ;
    pval_kdwt = min(pval_kdwt_1, pval_kdwt_2) ;
    [~,pval_kdoe_1] = kstest2(aas{3}, aas{1}, 'Tail', 'larger') ;
    [~,pval_kdoe_2] = kstest2(aas{1}, aas{3}, 'Tail', 'larger') ;
    pval_kdoe = min(pval_kdoe_1, pval_kdoe_2) ;
    % clf; ecdf(aas{2}); hold on; ecdf(aas{3})
    [~,pval_wtoe_1] = kstest2(aas{2}, aas{3}, 'Tail', 'larger') ;
    [~,pval_wtoe_2] = kstest2(aas{3}, aas{2}, 'Tail', 'larger') ;
    pval_wtoe = min(pval_wtoe_1, pval_wtoe_2) ;
    
    % Draw the pairwise pvalue bars
    %       'PStarThreshold' : Values which p-values must exceed (be smaller
    %           than or equal to) to earn a star. If PStarShowGT is false,
    %           significance is indicated with a star for every value in
    %           PStarThreshold which exceeds the p-value. If PStarShowGT is
    %           true, a p-value smaller than every element in PStarThreshold is
    %           indicated with (e.g.) '>***' instead of '****', to show the
    %           maximum measured precision has been exceeded. Default is
    %           [0.05, 0.01, 0.001, 0.0001].
    
    P = [0,  pval_kdwt, pval_kdoe; ...
        pval_kdwt, 0, pval_wtoe;...
        pval_kdoe, pval_wtoe, 0];
    E = [kd_unc, wt_unc, oe_unc] ;
    
    means = [kd,wt,oe] ;
    hold on;
    if qq ~= 4
        superbar([1,2,3],maxys, 'E', E, 'P', P,...
            'BarFaceColor', 'none', 'BarEdgeColor', 'none') % , 'ErrorbarColor', 'k') ;
    else
        superbar([1,2,3],[kd,wt,oe], 'E', E, 'P', P,...
            'BarFaceColor', faceColors, 'BarEdgeColor', 'none') % , 'ErrorbarColor', 'k') ;
    end
    
    disp(['Saving figure for ' metrics{qq}])
    if qq == 1
        ylim([-Inf, 280])
        saveas(gcf, fullfile(outdir, '2Dshapes_perim_pval_kstest.pdf'))
    elseif qq == 2
        ylim([-Inf, 1050])
        saveas(gcf, fullfile(outdir, '2Dshapes_area_pval_kstest.pdf'))
    elseif qq == 3
        % draw line for circle
        plot([0, 4], 2*pi / sqrt(pi) * [1,1], 'k--')
        ylim([3.5, 12])
        saveas(gcf, fullfile(outdir, '2Dshapes_shape_pval_kstest.pdf'))
        saveas(gcf, fullfile(outdir, '2Dshapes_shape_pval_kstest.png'))
    else
        error('index not understood')
    end
    % Save the statistical analysis to disk
    data = aa ;
    outputfn_tmp = fullfile(outdir, ['Curves_2Dshapes_' metrics{qq} '_comparison.mat']) ;
    disp(['Saving data to disk: ' outputfn_tmp])
    save(outputfn_tmp, ...
        'P', 'E', 'means', 'categories', 'data')
    close all

end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PART II                                                                 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Get metrics from resampled contours
overwrite = true ;
categories = {'KD', 'WT', 'OE'} ;
metrics = {'nProtrusions', 'Hdists'} ;

pix2um = 1 / 10.0974 ; % resolution in um / pix
% wtm = dlmread(fullfile(curvdir, 'control_mean.txt')) ;
% kdm = dlmread(fullfile(curvdir, 'knockdown_mean.txt')) ;
% oem = dlmread(fullfile(curvdir, 'overexp_mean.txt')) ;


wtc = dlmread(fullfile(curvdir, 'control_shapes.txt')) ;
kdc = dlmread(fullfile(curvdir, 'knockdown_shapes.txt')) ;
oec = dlmread(fullfile(curvdir, 'overexp_shapes.txt')) ;
AllCurves = { kdc, wtc, oec} ;

outFn = fullfile(outdir, ['Curves_2Dshapes_statistics.mat']) ;
if exist(outFn, 'file') && ~overwrite
    load(outFn,  ...
        'nProtrusions', 'Hdists', 'Curves', 'HdistsIJ')
else

    % Note data is in format x1, y1, x2, y2, ... xN, yN
    Curves = {} ;
    nProtrusions = {} ;
    for ii = 1:length(categories)
        curvs = AllCurves{ii} ;
        nExpts = size(curvs, 1) ;
        area = zeros(1, nExpts) ;
        perim = zeros(1, nExpts) ;
        nprotrusions = zeros(1, nExpts) ;
        curves = cell(nExpts, 1) ;

        for jj = 1:nExpts 
            disp(['shape/counting jj = ' num2str(jj) '/' num2str(nExpts)])

            xx = curvs(jj, 1:2:end) * pix2um ;
            yy = curvs(jj, 2:2:end) * pix2um ;

            % plot(xx, yy, '.')
            % axis equal 
            % pause(0.5)

            area(jj) = polyarea(xx, yy) ;
            perim(jj) = sum(vecnorm( diff([xx(:), yy(:)] ),2,2)) ;

            curves{jj} = [xx(:), yy(:)] ; 


            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Count protrusions
            radius = vecnorm([xx(:), yy(:)], 2, 2) ;
            radius2 = [radius; radius] ; % repeat to avoid edge effects
            % radsm = savgol(double(radius2), double(5), ...
            %     double(round(0.3 * length(radius)))) ; % smooth by 10% perimeter
            % [pks,locs,w,p] = findpeaks(radsm) ;
            % what is radius of a circle with the same area?
            rad_thres = sqrt(area(jj) / pi ) ; % threshold radius
            [pks, peakMag] = peakfinder(radius2, rad_thres*0.3, rad_thres, 1, false, false) ;
            
            % dprotrusions{jj} = radius2(pks) - rad_thres ;
            % dprotrusions_peakfinder{jj} = radius2(pks) ;
            
            npts = numel(xx(:)) ;
            pks = pks(pks > npts * 0.5 & pks < npts * 1.5) ;
            pks = mod(pks, npts) ;
            pks(pks == 0) = npts ;
            nprotrusions(jj) = length(pks) ;

            set(gcf, 'visible', 'off')
            plot([xx, xx(1)], [yy, yy(1)], '-')
            axis equal
            xlabel(['x [' char(181) 'm]'])
            ylabel(['y [' char(181) 'm]'])
            hold on
            scatter(xx(pks), yy(pks), 350, 'r', 'filled')
            pause(0.01)
            if ~exist(fullfile(curvdir, 'protrusionIms'), 'dir')
                mkdir(fullfile(curvdir, 'protrusionIms'))
            end
            saveas(gcf, fullfile(curvdir, 'protrusionIms', ['protrusions_' categories{ii} sprintf('_%03d.png', jj)]))
            clf

        end
        areas{ii} = area ;
        perims{ii} = perim ;
        shapes{ii} = perim ./ sqrt(area)  ; 
        nProtrusions{ii} = nprotrusions ;
        Curves{ii} = curves ; 
    end

    % Hausdorff distance between curves

    % You have the outlines of clusters as curves in a cell array like
    % curves = {curv1, curv2, curv3, ....}, where each curve is a Nx2 vector of
    % points, for ex as obtained by taking the outline of the black-and-white
    % mask. You also have the average profile as avgcurv. Then:

    Hdists = cell(3, 1) ;
    HdistsIJ = cell(3, 1) ;
    for qq = 1:length(categories)
        curves = Curves{qq} ;
        nExperiments = numel(curves) ;
        hdist = nan(nExperiments, nExperiments) ;
        dmyk = 1;
        for ii = 1:nExperiments
            disp(['Hausdorff ii = ' num2str(ii) '/' num2str(nExperiments)])
            othercurv = curves{ii} ;
            for jj = 1:nExperiments
                if jj ~= ii
                    curv = curves{jj} ; 
                    % For each point in the curv ask how far away is the reference curve
                    dists = zeros(length(curv), 1) ;
                    for kk = 1:length(curv)
                        dists(kk) = min(vecnorm(curv(kk,:) - othercurv, 2, 2)) ;
                    end
                    % Hausdorff distance is max(min(euclidean point-wise distance))
                    hdist(ii, jj) = max(dists) ;
                end
            end
        end

        % Note that these are not all independent, so let's symmetrize and take
        % the upper triangular matrix entries. There are N(N-1)/2 of them
        HdistsIJ{qq} = hdist ;
        hsym = 0.5 * (hdist + hdist') ;
        uhsym = zeros(1, nExperiments*(nExperiments-1)*0.5) ;
        dmyk = 1 ;
        for ii = 1:nExperiments
            for jj = ii+1:nExperiments
                uhsym(dmyk) = hsym(ii, jj) ;
                dmyk = dmyk + 1 ;
            end
        end
        Hdists{qq} = uhsym ;
    end
    disp('done')

    % Now hdist is an array of Hausdorff distances for this collection of
    % experiments. Done!


    % Statistics
    nprotr_avg = zeros(1, 3) ;
    nprotr_std = zeros(1, 3) ;
    nprotr_ste = zeros(1, 3) ;
    hdist_avg = zeros(1, 3) ;
    hdist_std = zeros(1, 3) ;
    hdist_ste = zeros(1, 3) ;
    for ii = 1:3
        nprotr_avg(ii) = mean(nProtrusions{ii}) ;
        nprotr_std(ii) = std(nProtrusions{ii}) ;
        nprotr_ste(ii) = std(nProtrusions{ii}) ./ sqrt(length(nProtrusions{ii})) ;
        hdist_avg(ii) = mean(Hdists{ii}) ;
        hdist_std(ii) = std(Hdists{ii}) ;
        hdist_ste(ii) = std(Hdists{ii}) ./ sqrt(length(Hdists{ii})) ;
        
    end

    save(outFn, ...
        'nprotr_avg', 'nprotr_std', 'nprotr_ste', ...
        'hdist_avg', 'hdist_std', 'hdist_ste', ...
        'perims', 'areas', 'shapes', 'nProtrusions', 'Hdists', 'Curves', 'HdistsIJ')
    close all
end



%% plot curves
for pp = 1
    close all
    for qq = 1:length(categories)
        subplot(2, length(categories), qq) 
        curves = Curves{qq} ;

        nExpts = length(curves) ;
        for ii = 1:nExpts 
            curv = curves{ii} ;
            plot([curv(:, 1); curv(1, 1)], ...
                [curv(:, 2); curv(1, 2)], '-')
            hold on;

            if ii == 1
                allcurvs = curv ;
            else
                allcurvs = [allcurvs; curv] ;
            end
        end

        % allcurvs = [allcurvs; zeros(20*nExpts, 2)] ;
        axis equal
        axis square
        xlabel(['AP position [' char(181) 'm]'])
        ylabel(['DV position [' char(181) 'm]'])
        xlim([-30, 30])
        ylim([-30, 30])
        title(categories{qq})

        % Histcounts density plot
        subplot(2, length(categories), qq+length(categories))
        % nn = histcounts2(divv_tp(:), H2vn_tp(:), ...
        %     xedges, yedges) ;
        % imagesc(xedges, yedges, log10(nn)') ;

        gridx1 = -30:.5:30;
        gridx2 = -30:.5:30;
        [x1g,x2g] = meshgrid(gridx1, gridx2);
        x1 = x1g(:);
        x2 = x2g(:);
        xi = [x1 x2];
        [dens, xi, h3] = ksdensity(allcurvs, xi);
        xx = reshape(xi(:, 1), size(x1g)) ;
        yy = reshape(xi(:, 2), size(x1g)) ;
        densG = reshape(dens, size(x1g)) ;
        imagesc(gridx1, gridx2, densG)
        axis equal ;
        if pp == 1
            tmp = cubehelix(256,0.,1,1,0.4,[0,1],[0,1]) ;
            colormap(flipud(tmp))
        else
            tmp = cubehelix(256,0.,1,1,0.4,[0,1],[0,1]) ;
            % tmp = [1,1,1; tmp] ;
            colormap(tmp)
        end
        % tmp = cubehelix(256,1.,1,1,2.,[0,1],[0,1]) ;
        % colormap(tmp)
        xlim([-30, 30])
        ylim([-30, 30])
        set(gca, 'Ydir', 'normal')
        xlabel(['AP position [' char(181) 'm]'])
        ylabel(['DV position [' char(181) 'm]'])

        caxis([0, 0.0019]) 
        set(gcf, 'color', 'w')

    end
    saveas(gcf, fullfile(outdir, sprintf('Curves_overlaid_curves_%02d.pdf', pp)))
    clf
    subplot(2, 3, 1)
    colorbar
    if pp == 1
        colormap(flipud(tmp))
    else
        colormap(tmp)
    end
        caxis([0, 0.0019])
        set(gcf, 'color', 'w')
        axis off
    saveas(gcf, fullfile(outdir, sprintf('Curves_overlaid_curves_colorbar_%02d.pdf', pp)))
end


%% Plot Hausdorff distance indices

% smooth plot
profile_step = [ 0.03, 0.01] ;

faceColors = [31, 177, 3; ...
    110,110,110; ...
    203, 41,123] ./ 255.0 ;

% smooth violin plot of area, perim, shape
for qq = 1:length(metrics)
    maxy = 0 ;
    miny = Inf ;
    maxys = [0,0,0] ;
    close all
    figure('units', 'centimeters', 'position', [0, 0, 10, 10])
    
    aas = {} ;
    for ii = 1:length(categories)
        if qq == 1
            aa = nProtrusions{ii} ;
        elseif qq == 2
            aa = Hdists{ii} ;
        end
        
        aas{ii} = aa ;
        
        if qq ~= 1
            violin(aa(:), 'x', [ii,ii+1],...
                'facecolor',faceColors(ii,:),'edgecolor','none', 'facealpha', 1, ...
                'bw',profile_step(qq)*(max(aa) - min(aa)), 'mc','k-','medc', '', 'plotlegend', false)
        else
            ptrs = 0:10 ; % possible #protrusions
            % counts = histc(aa, ptrs) ; % count occurrence of each occuring number of protrusions 
            boxyViolinHistogram(ii, aa, max(aa)) ;
            hold on;
        end
        hold on;
        maxys(ii) = max(aa(:)) ;
        miny = min(miny, min(aa(:))) ;
        maxy = max(maxy, max(aa(:))) ;
    end
    xlim([0, 4])
    yLow = miny * 0.9 ;
    yHigh = maxy * 1.1 ;
    topY = yHigh * 1.1 ;
    ylim([yLow, yHigh])
    xticks([1,2,3])
    set(gca, "xtickLabels", categories)
    if qq == 1
        title('Border cell cluster protrusions')
        ylabel(['number of protrusions'])
        saveas(gcf, fullfile(outdir, 'Curves_2Dshapes_protrusions.pdf'))
        clf
    elseif qq == 5
        title('Hausdorff distance between profiles')
        ylabel(['distance between profiles [' char(181) 'm]'])
        saveas(gcf, fullfile(outdir, 'Curves_2Dshapes_hdists.pdf'))
    end
    % 
    % 
    % Add p value and save plot again
    if qq == 1    
        wt = nprotr_avg(2) ;
        kd = nprotr_avg(1) ;
        oe = nprotr_avg(3) ;
        wt_unc = nprotr_ste(2) ;
        kd_unc = nprotr_ste(1) ;
        oe_unc = nprotr_ste(3) ;
    elseif qq == 2
        wt = hdist_avg(2) ;
        kd = hdist_avg(1) ;
        oe = hdist_avg(3) ;
        wt_unc = hdist_ste(2) ;
        kd_unc = hdist_ste(1) ;
        oe_unc = hdist_ste(3) ;
    else
        error('Index not recognized')
    end
    % 
    % num = -abs(kd-wt) ;
    % denom = sqrt(kd_unc.^2 + wt_unc.^2) ;
    % zscore_kdwt = num / denom ;
    % pval_kdwt = normcdf(zscore_kdwt);
    % 
    % num = -abs(kd-oe) ;
    % denom = sqrt(kd_unc.^2 + oe_unc.^2) ;
    % zscore_kdoe = num / denom ;
    % pval_kdoe = normcdf(zscore_kdoe);
    % 
    % num = -abs(wt-oe) ;
    % denom = sqrt(oe_unc.^2 + wt_unc.^2) ;
    % zscore_wtoe = num / denom ;
    % pval_wtoe = normcdf(zscore_wtoe);

    %if qq == 1
    % preview it:
    % clf; ecdf(aas{1}); hold on; ecdf(aas{2})
    [~,pval_kdwt_1] = kstest2(aas{2}, aas{1}, 'Tail', 'larger') ;
    [~,pval_kdwt_2] = kstest2(aas{1}, aas{2}, 'Tail', 'larger') ;
    pval_kdwt = min(pval_kdwt_1, pval_kdwt_2) ;
    [~,pval_kdoe_1] = kstest2(aas{3}, aas{1}, 'Tail', 'larger') ;
    [~,pval_kdoe_2] = kstest2(aas{1}, aas{3}, 'Tail', 'larger') ;
    pval_kdoe = min(pval_kdoe_1, pval_kdoe_2) ;
    % clf; ecdf(aas{2}); hold on; ecdf(aas{3})
    [~,pval_wtoe_1] = kstest2(aas{2}, aas{3}, 'Tail', 'larger') ;
    [~,pval_wtoe_2] = kstest2(aas{3}, aas{2}, 'Tail', 'larger') ;
    pval_wtoe = min(pval_wtoe_1, pval_wtoe_2) ;
    
    % Draw the pairwise pvalue bars
    %       'PStarThreshold' : Values which p-values must exceed (be smaller
    %           than or equal to) to earn a star. If PStarShowGT is false,
    %           significance is indicated with a star for every value in
    %           PStarThreshold which exceeds the p-value. If PStarShowGT is
    %           true, a p-value smaller than every element in PStarThreshold is
    %           indicated with (e.g.) '>***' instead of '****', to show the
    %           maximum measured precision has been exceeded. Default is
    %           [0.05, 0.01, 0.001, 0.0001].
    
    P = [0,  pval_kdwt, pval_kdoe; ...
        pval_kdwt, 0, pval_wtoe;...
        pval_kdoe, pval_wtoe, 0];
    E = [kd_unc, wt_unc, oe_unc] ;
    
    means = [kd,wt,oe] ;
    hold on;
    if qq ~= 1
        superbar([1,2,3],maxys, 'E', E, 'P', P,...
            'BarFaceColor', 'none', 'BarEdgeColor', 'none') % , 'ErrorbarColor', 'k') ;
    else
        superbar([1,2,3],[kd,wt,oe], 'E', E, 'P', P,...
            'BarFaceColor', faceColors, 'BarEdgeColor', 'none') % , 'ErrorbarColor', 'k') ;
    end
    
    if qq == 1
        title('Border cell cluster protrusions')
        ylabel('number of protrusions')
        xticks([1,2,3])
        set(gca, "xtickLabels", categories)
        ylim([0, 3.5])
        saveas(gcf, fullfile(outdir, 'Curves_2Dshapes_protrusions_pval_kstest.pdf'))
        saveas(gcf, fullfile(outdir, 'Curves_2Dshapes_protrusions_pval_kstest.png'))
    elseif qq == 2
        ylim([-Inf, 20])
        saveas(gcf, fullfile(outdir, 'Curves_2Dshapes_hdists_pval_kstest.pdf'))               
    else
        error('code for this case here')
    end
    % Save the statistical analysis to disk
    data = aa ;
    save(fullfile(outdir, ['Curves_2Dshapes_' metrics{qq} '_comparison.mat']), ...
        'P', 'E', 'means', 'categories', 'data')
    close all

end

