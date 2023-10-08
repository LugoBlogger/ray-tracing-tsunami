% TraceContour_Arthur.m
%  Use Fermat's principle to find the trajectory that gives minimum time.
%
% Author: Henokh Lugo Hariyanto (henokh.lugo.h@mail.ugm.ac.id)
% 
% Logs
%
% [2022/04/18]
% - This is an implementation of (Arthur, 1946)
% 
% [2022/04/20]
% - Add a feature to rotate the bathymetry before selecting the area
% - Overlay significant wave height data


clear;

%% -- user-defined input
toolboxPath = 'deg_utm_conversion/';
addpath(genpath(toolboxPath));

% ---- Socotra survey bathy
bathyDir = './ieodo_digitization/';
bathyFilename = 'Ieodo_lon_lat_depth.csv';

towerCoor = [32, 7, 22.63;      % This is the Ieodo Ocean Research Station position (lat, lon)
             125, 10, 56.81];

zeroRefCoor = [32, 7, 42.74;      % This is the Socotra rock position (lat, lon)
               125, 10, 49.80];   
        
T_colormap = [...
      [214, 61, 41]/255;
      [243, 134, 79]/255;
      [238, 240, 168]/255;
      [170, 216, 232]/255;
      [69, 117, 180]/255];
x_colormap_space = [0; 62; 124; 190; 255];

khoa_colormap = interp1(x_colormap_space/255, T_colormap, linspace(1,0,255));

% --- moments data (var, skew, kurt) from Hawassi simulation
hawassiDataDir = './hawassi_simulation_data/';
%momentsDataFilename = 'var_skew_kurt_T19_x_-1600_4_959_y_-1320_4_1796.dat';   % 4 is the gridsize
momentsDataFilename = 'rotated_dom_45_var_skew_kurt_T15_x_-1600_4_959_y_-1200_4_1036.dat';


%% -- read bathymetry dataset
depthData = readBathyDataset(bathyDir, bathyFilename);


%% -- transform coordinates (lat, lon) to UTM (northing easting)
% Also interpolating to given selected area

% selected_xRange = [NaN, NaN];
% selected_yRange = [-500, 2000];

% selected_xRange = [NaN, NaN];
% selected_yRange = [NaN, NaN];

% -- the numbers are round-off
% selected_xRange = [-1600, 1000];  % for case 1 (wave influx from north)
% selected_yRange = [-1300, 1800];
% rotateAngle = 0;   % in degree

selected_xRange = [-1600, 1000];  % for rotating bathymetry
selected_yRange = [-1250, 1000];
rotateAngle = 45;   % in degree
  
Nx_selected = 512;     % the best setting to have peak of the selected line around -4.6
Ny_selected = 512;

smoothingFac = 1.5;       % 0 or -1 is no smoothin factor
                          % 1.5 is the best for this bathymetry in my
                          % opinion (less noisy contour line)

[X_selected, Y_selected, depth_selected] = transformToUTM(depthData, zeroRefCoor, ...
  selected_xRange, selected_yRange, Nx_selected, Ny_selected, smoothingFac, ...
  rotateAngle);


%% -- load moments data and get Hs
isLoadMoments = true;
momentsData = readMomentsData(hawassiDataDir, momentsDataFilename);


%% -- compute group speed
%selectedTp = 8.36;   % mean wave period of Tp_influx = 10 s
%selectedTp = 12;
selectedTp = 12.56;   % mean wave period of Tp_influx = 15 s
%selectedTp = 15;
%selectedTp = 15.9;   % mean wave period of Tp_influx = 19 s
%selectedTp = 19;


numOfLevel = 100;     % define the group speed interval here
[CgMesh, kMesh, Cg_interpolant] = computeGroupSpeed(selectedTp, ...
  X_selected, Y_selected, depth_selected);
levelLineCgMesh = linspace(0, max(CgMesh(:)), numOfLevel);


%% -- plot contour Cg
isPlot_contourCg = false;
levelLineCgMeshPlot = [0.7, 2, 4, 7, 10, 13, 13.5, 13.7, 13.9, 14];
%levelLineCgMeshPlot = levelLineCgMesh;

xlimPlot = [];
ylimPlot = [];

plot_contourCg(X_selected, Y_selected, CgMesh, levelLineCgMeshPlot, ...
  xlimPlot, ylimPlot, isPlot_contourCg)


%% -- get min path
% for non-rotating bathymetry (wave from North)
% x0_initRay = -1600:25:1000;
% y0_initRay = ones(size(x0_initRay)) * selected_yRange(2);
% init_dir = [zeros(length(x0_initRay), 1), -ones(length(x0_initRay), 1)];

% for non-rotating bathymetry (wave from North)  (dense lines)
% x0_initRay = -1600:1:1000;     % 491 secs
% y0_initRay = ones(size(x0_initRay)) * selected_yRange(2);
% init_dir = [zeros(length(x0_initRay), 1), -ones(length(x0_initRay), 1)]; 
 
% for non-rotating bathymetry (wave from South-West)
% x0_initRay = [-1600*ones(1, length(-1300:50:1800)), -1600:50:1000];
% y0_initRay = [-1300:50:1800, -1300*ones(1, length(-1600:50:1000))];
% init_dir = [ones(length(x0_initRay), 1), ones(length(x0_initRay), 1)] / sqrt(2); 

% for rotating bathymetry (wave from South-West -> wave from South)
% x0_initRay = -1600:25:1000;
% y0_initRay = -1250*ones(1, length(x0_initRay));
% init_dir = [zeros(length(x0_initRay), 1), ones(length(x0_initRay), 1)]; 

% for rotating bathymetry (wave from South-West -> wave from South)   (dense lines)
%x0_initRay = -1600:1:1000;
x0_initRay = -1600:25:1000;
y0_initRay = -1250*ones(1, length(x0_initRay));
init_dir = [zeros(length(x0_initRay), 1), ones(length(x0_initRay), 1)]; 

search_radius = 1;   % in meter(s)
delta_n = 0.1;       % in meter(s); how far the step along direction normal to the ray path
maxPathStep = 5000;                         

minTimePath = computeMinPath(X_selected, Y_selected, Cg_interpolant, ...
  x0_initRay, y0_initRay, init_dir, search_radius, delta_n, maxPathStep);


%% -- plot minimum path (ray of refraction)
isPlot_minPath = true;
%levelLineCgMeshPlot = 9:0.01:15;
%levelLineCgMeshPlot = [];
levelLineDepthMeshPlot = [];
xRange = [-500, 300];
yRange = [-700, 400];

plotUnderlay.momentsData = momentsData;
plotUnderlay.underlay_Hs = true;

plotMinPath(X_selected, Y_selected, minTimePath, depth_selected, ...
  levelLineDepthMeshPlot, search_radius, delta_n, xRange, yRange, plotUnderlay, isPlot_minPath);

%% -- plot surface created by minimum path (each ray at different z)
isPlot_surfMinPath = false;
plotSurfMinPath(X_selected, Y_selected, depth_selected, levelLineDepthMeshPlot, minTimePath, isPlot_surfMinPath)

%% -- plot histogram2 to see the distribution of points of rays
isPlot_bivarHisto = false;
plotBivarHisto(X_selected, Y_selected, minTimePath, isPlot_bivarHisto)


%% -- function declarations
function plotBivarHisto(X_selected, Y_selected, minTimePath, isPlot_bivarHisto)
  if isPlot_bivarHisto
    fprintf('   plotBivarHist: ');
    tic;
    
    % collect all points into scattered data
    scattered_xy = minTimePath(1).xy;
    scattered_xy = scattered_xy(~isnan(scattered_xy(:, 1)),:);
%     scattered_z = ones(size(scattered_xy, 1), 1);
    N_rays = length(minTimePath);
    for i = 2:N_rays
      ith_ray = minTimePath(i).xy;
      ith_ray = ith_ray(~isnan(ith_ray(:, 1)), :);
      temp = [scattered_xy; ith_ray];
      
%       temp_z = [scattered_z; ones(size(ith_ray, 1), 1) * i];
      
      scattered_xy = temp;
%       scattered_z = temp_z;
    end
    
%     F_interpolant = scatteredInterpolant(scattered_xy(:, 1), scattered_xy(:, 2), ...
%       scattered_z, 'natural', 'none');
% 
%     raysDensity = F_interpolant(X_selected, Y_selected);
    
    
    fig = figure('Color', 'white', 'Units', 'centimeters', ...
      'Position', [15, 8, 16, 14]);
      ax = axes('Parent', fig);
      
%       Nx_bins = 28;
%       Ny_bins = 32;

      Nx_bins = 512;
      Ny_bins = 512;

      
      if isempty(Nx_bins)
        histo2_vecs = histogram2(ax, scattered_xy(:, 1), scattered_xy(:, 2), ...
          'DisplayStyle', 'tile');
      else
        histo2_vecs = histogram2(ax, scattered_xy(:, 1), scattered_xy(:, 2), ...
          [Nx_bins, Ny_bins], 'DisplayStyle', 'tile');
      end
        

%       imagesc(ax, raysDensity, 'XData', X_selected(1, :), 'YData', Y_selected(:, 1));
%       set(ax, 'FontSize', 12, 'YDir', 'normal');

      %caxis(ax, [min(histo2_vecs.Data(:)), max(histo2_vecs.Data(:))]);
      cbarHandler = colorbar(ax);
      colormap(ax, 'jet');
      ylabel(cbarHandler, 'density of points along rays');
      set(cbarHandler, 'FontSize', 12);

      set(ax, 'FontSize', 12);
        
      xlabel(ax, 'easting [m]');
      ylabel(ax, 'northing [m]');
      
      daspect(ax, [1, 1, 1]);
      
      fprintf('%.2f [s]\n', toc);
  end
end

function plotSurfMinPath(X_selected, Y_selected, depthMesh, levelLineDepthMeshPlot, minTimePath, isPlot_surfMinPath)
  if isPlot_surfMinPath
    fprintf('   plotSurfMinPath: ');
    tic;
    
    matlabBlue = [0, 114, 189]/255.;
    N_rays = length(minTimePath);
    N_points = length(minTimePath(1).xy(:,1));
    zSpace = 1:N_rays;
    
    fig = figure('Color', 'white', 'Units', 'centimeters', ...
      'Position', [20, 8, 18, 14]);
      ax = axes('Parent', fig);
      
      % -- plot contour depth
      if ~isempty(levelLineDepthMeshPlot)
        [contourData, contourHandler] = contour(ax, X_selected, Y_selected, depthMesh, ...
          'LevelList', levelLineDepthMeshPlot, 'ShowText', 'off');
      else
        [contourData, contourHandler] = contour(ax, X_selected, Y_selected, depthMesh, ...
          'ShowText', 'off');
      end
      contourHandler.LineColor = [1, 1, 1]*0.5;
      contourHandler.LineStyle = '--';
      
      hold(ax, 'on');
      for i = 1:N_rays
        plot3(ax, minTimePath(i).xy(:,1), minTimePath(i).xy(:,2), ...
          ones(N_points,1)*zSpace(i), 'Color', matlabBlue);
      end
      
      view(ax, 3);
      
      set(ax, 'FontSize', 12);
      grid(ax, 'on');
      box(ax, 'off');

      xlabel(ax, 'easting [m]');
      ylabel(ax, 'northing [m]');
      zlabel(ax, 'i-th wave ray');
      
      fprintf('%.2f [s]\n', toc);
  end
end

function plotMinPath(X_selected, Y_selected, minTimePath, ...
  depthMesh, levelLineDepthMeshPlot, search_radius, delta_n,  xRange, yRange, ...
  plotUnderlay, isPlot_minPath)

  if isPlot_minPath
    fprintf('   plotMinPath: ');
    tic;
    
    matlabBlue = [0, 114, 189]/255.;
    %matlabRed = [163, 20, 46]/255;
    %matlabSky = [77, 191, 237]/255;
    %matlabBrown = [217, 83, 25]/255;

    
    fig = figure('Color', 'white', 'Units', 'centimeters', ...
      'Position', [20, 6, 18, 14]);
      ax = axes('Parent', fig);
      
      if plotUnderlay.underlay_Hs
        momentsData = plotUnderlay.momentsData;
        xSpace_Hs = momentsData.xSpace;
        ySpace_Hs = momentsData.ySpace;
        Hs = reshape(4. * sqrt(momentsData.moments(:, 1)), length(xSpace_Hs), length(ySpace_Hs));
        
        imagesc(ax, Hs', 'XData', xSpace_Hs, 'YData', ySpace_Hs, ...
          'AlphaData', 1);
  
        caxis(ax, [min(Hs(:)), max(Hs(:))]);
        cbarHandler = colorbar(ax);
        colormap(ax, 'jet');
        ylabel(cbarHandler, 'Significant wave height [m]');
        set(cbarHandler, 'FontSize', 12);
        
        hold(ax, 'on');
      end
  

      if ~isempty(levelLineDepthMeshPlot)
        [contourData, contourHandler] = contour(ax, X_selected, Y_selected, depthMesh, ...
          'LevelList', levelLineDepthMeshPlot, 'ShowText', 'off');
      else
        [contourData, contourHandler] = contour(ax, X_selected, Y_selected, depthMesh, ...
          'ShowText', 'off');
      end
      contourHandler.LineColor = [1, 1, 1]*0.5;
      contourHandler.LineStyle = '--';
%       clabel(contourData, contourHandler, 'FontSize', 12, 'Color', 'k', ...
%         'LabelSpacing', 100);
      
      
      N_init_dir = length(minTimePath);
      for i = 1:N_init_dir
        hold(ax, 'on');
        plot(ax, minTimePath(i).xy(:, 1), minTimePath(i).xy(:, 2), 'Color', matlabBlue, ...
          'LineWidth', 0.1);
      end

    if ~isempty(xRange)
      xlim(ax, xRange);
    end
    
    if ~isempty(yRange)
      ylim(ax, yRange);
    end
      
      
    set(ax, 'YDir', 'normal', 'FontSize', 12);
    xlabel(ax, 'easting [m]');
    ylabel(ax, 'northing [m]');

    title(ax, sprintf('Depth contour and Hs; h\\_s = %g [m]; h\\_n = %g [m]', ...
      search_radius, delta_n));
    
    daspect(ax, [1, 1, 1]);
    
    fprintf('%.2f [s]\n', toc);
  end
end

function minTimePath = computeMinPath(X_selected, Y_selected, ...
  Cg_interpolant, x0_initRay, y0_initRay, init_dir, search_radius, ...
  delta_n, maxPathStep)

  fprintf('   Compute minTimePath: ');
  tic;

  N_init_points = length(x0_initRay);
  minTimePath = struct('xy', repmat({[]}, 1, N_init_points));
  
  matrixRot90 = [0, -1;
                 1, 0];

  for i = 1:N_init_points
    message = sprintf(' (procees %d/%d init points)\n', i, N_init_points);
    fprintf('%s', message);
    
    path_j = ones(maxPathStep, 2) * NaN;
    x0 = x0_initRay(i);
    y0 = y0_initRay(i);
    path_j(1, :) = [x0, y0];
    
    dir_i = init_dir(i, :) / norm(init_dir(i, :));
    vec_advanced = dir_i * search_radius;
    
    x0 = vec_advanced(1) + x0;
    y0 = vec_advanced(2) + y0;
    path_j(2, :) = [x0, y0];
    
    for j = 2:maxPathStep-1
      
      normal_dirToRayPath = (matrixRot90 * dir_i')';
      xy_advanced_normal_dir = (delta_n * normal_dirToRayPath) + [x0, y0];
      Cg_advanced_normal_dir = Cg_interpolant(xy_advanced_normal_dir(1), xy_advanced_normal_dir(2));
      Cg_j = Cg_interpolant(x0, y0);
      
      rotAngle = search_radius * (-1/Cg_j) * (Cg_advanced_normal_dir - Cg_j)/delta_n;
      dir_j = ([cos(rotAngle), -sin(rotAngle);
               sin(rotAngle), cos(rotAngle)] * vec_advanced')';
      
      xy_advanced = (dir_j * search_radius) + [x0, y0];

      path_j(j+1, :) = xy_advanced;

      vec_advanced = xy_advanced - [x0, y0];
      
      % if xy_advance_min outside region, we stop the iteration
      if xy_advanced(1) < X_selected(1,1) || xy_advanced(1) > X_selected(1, end) ...
         || xy_advanced(2) > Y_selected(end, 1) || xy_advanced(2) < Y_selected(1, 1)
        break;
      end
      
      % update x0, y0
      x0 = vec_advanced(1) + x0;
      y0 = vec_advanced(2) + y0;
      
      % update direction
      dir_i = dir_j;
      
      % -- check output by plot it
      %plot_angle_search(xy_advanced_arr, x0, y0, Cg_over_perimeter, Cg_at_x0_y0)
    end
    
    
    minTimePath(i).xy = path_j;
    fprintf(repmat('\b', 1, length(message)));
  end
  
  fprintf('%.2f [s]\n', toc);
end

function plot_contourCg(Xprocess, Yprocess, CgMesh, levelLineCgMeshPlot, ...
  xlimPlot, ylimPlot, isPlot_contourCg)

  if isPlot_contourCg
    fprintf('   plot contour Cg: ');
    tic;
    
    fig = figure('Color', 'w', 'Units', 'centimeters', 'Position', [15, 3, 18, 14]);
      ax = axes('Parent', fig);


      xSpace = Xprocess(1,:);
      ySpace = Yprocess(:,1);

      imagesc(ax, CgMesh, 'XData', xSpace, 'YData', ySpace);

      hold(ax, 'on')
        [contourData, contourHandler] = contour(ax, Xprocess, Yprocess, CgMesh, ...
          'LevelList', levelLineCgMeshPlot, 'ShowText', 'on');
         contourHandler.LineColor = 'k';
         contourHandler.LineStyle = '-';
         clabel(contourData, contourHandler, 'FontSize', 12, 'Color', 'k', ...
           'LabelSpacing', 100);

      caxis(ax, [min(levelLineCgMeshPlot), max(levelLineCgMeshPlot)]);

      cbarHandler = colorbar(ax);
      ylabel(cbarHandler, 'Group velocity [m/s]');
      set(cbarHandler, 'FontSize', 12);

      if ~isempty(xlimPlot) && ~isempty(ylimPlot)    
        xlim(ax, xlimPlot);
        ylim(ax, ylimPlot);
      end

      grid(ax, 'on'); 
      set(ax, 'YDir', 'normal', 'FontSize', 12, 'GridColor', [0, 0, 0, 0.2]);
      xlabel(ax, 'easting [m]');
      ylabel(ax, 'northing [m]');

      title(ax, 'group speed');
      daspect(ax, [1, 1, 1]);
      
    fprintf('%.2f [s]\n', toc);
  end
end

function [CgMesh, kMesh, Cg_interpolant] = computeGroupSpeed(selectedTp, ...
  X_selected, Y_selected, IM_combine)


  fprintf('   Compute group speed and its interpolant: ');
  tic;
  
  freq = 2*pi/selectedTp;
  g = 9.81;

  % (Guo, 2002) - `Simple and explicit solution of wave dispersion equation`
  beta = 2.4908;   % transitional shape parameter
  kMesh = (freq.^2/g) .* (1 - exp(-(freq.*sqrt(abs(IM_combine)/g)).^beta)) .^ (-1/beta);

  positiveDepth = -IM_combine;
  positiveDepth(positiveDepth <= 0) = 0;

  % (Dingemans, 1994) `Water wave propagation over uneven bottoms` Eqs. 2.29a and 2.29c
  CpMesh = sqrt(g * tanh(kMesh .* positiveDepth) ./ kMesh);
  CgMesh = 0.5 * CpMesh .* (1 + 2*kMesh .* positiveDepth ./ sinh(2*kMesh .* positiveDepth));

  % get interolation function for Cg
  Cg_interpolant = griddedInterpolant(X_selected', Y_selected', CgMesh', 'spline');
  
  fprintf('%.2f [s]\n', toc);
end

function momentsData = readMomentsData(hawassiDataDir, momentsDataFilename)
  fprintf('   Read %s: ', momentsDataFilename);
  tic;
    momentsData.moments = load([hawassiDataDir, momentsDataFilename]);
    
  % -- get xSpace and ySpace
  [~, name, ~] = fileparts(momentsDataFilename);
  [startIdx, endIdx] = regexp(name, '-?[0-9]+_-?[0-9]+_-?[0-9]+');
  
  tempData = cellfun(@str2num,split(name(startIdx(1):endIdx(1)), '_'));
  momentsData.xSpace = tempData(1):tempData(2):tempData(3);
  
  tempData = cellfun(@str2num,split(name(startIdx(2):endIdx(2)), '_'));
  momentsData.ySpace = tempData(1):tempData(2):tempData(3);
    
  fprintf('%.2f [s]\n', toc);
end

function [X_selected, Y_selected, IM_selected] = transformToUTM(depthData, ...
  zeroRefCoor, selected_xRange, selected_yRange, Nx_selected, Ny_selected, ...
  smoothingFac, rotateAngle)

  [UTM_x, UTM_y, ~] = deg2utm(depthData.lat', depthData.lon');
  
  refPos = dms2degrees(zeroRefCoor);
  [Ref_x, Ref_y, ~] = deg2utm(refPos(1), refPos(2));
  
  UTM_x = UTM_x - Ref_x;
  UTM_y = UTM_y - Ref_y;
  
  % Rotate the data
  if abs(rotateAngle) > 1e-6
    UTM_xy = [cosd(rotateAngle), -sind(rotateAngle);
              sind(rotateAngle),  cosd(rotateAngle)] * [UTM_x' ; UTM_y'];
    UTM_x = UTM_xy(1,:)';
    UTM_y = UTM_xy(2,:)';
  end
  
  % Interpolate to selected xRange and yRange
  pointGroup = [UTM_x, UTM_y, depthData.depth'];
  depthData = unique(pointGroup, 'rows');
  F_interpolant = scatteredInterpolant(depthData(:, 1), depthData(:, 2), depthData(:, 3), ...
    'natural', 'linear');

  % -- too slow
  %F_interpolant = fit([depthData(:, 1), depthData(:, 2) ], depthData(:, 3), 'cubicinterp');

  
  if isnan(selected_xRange(1))
    selected_xRange(1) = min(UTM_x);
  end
  if isnan(selected_xRange(2))
    selected_xRange(2) = max(UTM_x);
  end
  if isnan(selected_yRange(1))
    selected_yRange(1) = min(UTM_y);
  end
  if isnan(selected_yRange(2))
    selected_yRange(2) = max(UTM_y);
  end
  
  X_space = linspace(selected_xRange(1), selected_xRange(2), Nx_selected);
  Y_space = linspace(selected_yRange(1), selected_yRange(2), Ny_selected);
  
  [X_selected, Y_selected] = meshgrid(X_space, Y_space);
  
  IM_selected = F_interpolant(X_selected, Y_selected);
  
  
  if smoothingFac > 0
    IM_selected = imgaussfilt(IM_selected, smoothingFac);
  end
end

function depthData = readBathyDataset(bathyDir, bathyFilename)
  
  fprintf('   Read %s\n', bathyFilename);
  fileId = fopen([bathyDir, bathyFilename]);
  
  rowEntry = fgetl(fileId);
  initRow = 1;
  latitudeStr  = "";
  longitudeStr = "";
  depthStr     = "";
  
  while ischar(rowEntry)
    message = sprintf('   Processing %d item(s)\n', initRow);
    fprintf('%s', message);
    
    % remove leading and trailing spaces
    rowEntry = strip(rowEntry);
    
    % split by comma into a cell
    rowEntry = strsplit(rowEntry, ',');
    
    % because splitting only scan single space, we need to 
    % mask out each element which is an empty string
    rowEntry = rowEntry(~(rowEntry == ""));
    
    % please be careful to use curly brace instead of parenthesis
    % when indexing rowEntry
    tempConcat   = [longitudeStr, [rowEntry{1}]];
    longitudeStr = tempConcat;
    
    tempConcat  = [latitudeStr, [rowEntry{2}]];
    latitudeStr = tempConcat;

    tempConcat = [depthStr, [rowEntry{3}]];
    depthStr   = tempConcat;
    
    initRow = initRow + 1;
    rowEntry = fgetl(fileId);
    
    fprintf(repmat('\b', 1, length(message)));
  end
  
  fprintf('   Total item(s): %d\n', initRow - 1);
  
  fclose(fileId);
  
  longitudeStr = longitudeStr(3:end);   % exclude the header and init str
  latitudeStr  = latitudeStr(3:end);    % exclude the header and init str
  depthStr     = depthStr(3:end);       % exclude the header and init str
  
  depthData.lon   = str2double(longitudeStr);
  depthData.lat   = str2double(latitudeStr);
  depthData.depth = -str2double(depthStr);
end