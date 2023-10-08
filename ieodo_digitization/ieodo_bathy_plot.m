% ieodo_bathy_plot.m
%   This is a plot for a bathymetry around Ieodo rock.
%   The bathymetry is a digitization from the nautical chart around Ieodo
%   rock.
%   The bathymetry data is in .csv. We avoid to use readtable() because of
%   the performance issue.

clear

%% -- user-defined input
toolboxPath = 'deg_utm_conversion/';
addpath(genpath(toolboxPath));

bathyDir = '../ieodo_digitization/';

bathyFilename = 'Ieodo_lon_lat_depth.csv';
%bathyFilename = 'Ieodo_lon_lat_depth_scatterOnly.csv';   % less accurate

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
             

%% -- read bathymetry dataset
depthData = readBathyDataset(bathyDir, bathyFilename);

%% -- transform coordinates (lat, lon) to UTM (northing easting)
% Also interpolating to given selected area

% selected_xRange = [NaN, NaN];
% selected_yRange = [-500, 2000];

% selected_xRange = [NaN, NaN];
% selected_yRange = [NaN, NaN];

% selected_xRange = [-1000, 1000];
% selected_yRange = [-2000, 1000];

% -- the numbers are round-off
selected_xRange = [-1600, 1000];      % case 1 (wave influx from north)
selected_yRange = [-1300, 1800];
rotateAngle = 0;        % in degree

% -- the numbers are round-off
% selected_xRange = [-1600, 1000];      % case 2 (wave influx from south west)
% selected_yRange = [-1200, 1000];
% rotateAngle = 45;        % in degree

Nx_selected = 512;     % the best setting to have peak of the selected line around -4.6
Ny_selected = 1024;

smoothingFac = 1.5;       % 0 or -1 is no smoothin factor
                          % 1.5 is the best for this bathymetry in my
                          % opinion (less noisy contour line)

[X_selected, Y_selected, IM_selected, depth_interpolant, towerUTM_coor] ...
  = transformToUTM(depthData, zeroRefCoor, selected_xRange, selected_yRange, ...
      Nx_selected, Ny_selected, smoothingFac, rotateAngle, towerCoor);


%% -- plot bathymetry
%levelList = -[5, 10, 20, 30, 40, 50, 60];
%levelList = -4:-2:-60;
levelList = [-5, -10, -20, -30, -40, -50, -52, -54];
isPlot_bathy = true;

plot_bathy(X_selected, Y_selected, IM_selected, levelList, khoa_colormap, ...
  smoothingFac, isPlot_bathy)


%% -- plot 3D bathymetry
isPlot_3Dbathy = false;
plot_3Dbathy(X_selected, Y_selected, IM_selected, khoa_colormap, ...
  towerUTM_coor, depth_interpolant, isPlot_3Dbathy)

%% -- plot bathymetry with straight line
isPlot_bathyWithLine = false;
% lineStartPoints = [52.35, -117.6];    % (UTM_x, UTM_y) relative to zeroRefUTM
% lineEndPoints   = [-446.1, 1123];   

% lineStartPoints = [-221, -500];    % (UTM_x, UTM_y) relative to zeroRefUTM
% lineEndPoints   = [-221, 1500];    

% lineStartPoints = [-483.24, 1500];    % (UTM_x, UTM_y) relative to zeroRefUTM
% lineEndPoints   = [322.16, -1000];  

% case 1
lineStartPoints = [-60, 1500];    % (UTM_x, UTM_y) relative to zeroRefUTM
lineEndPoints   = [-60, -1000];  

plot_bathyWithLine(X_selected, Y_selected, IM_selected, levelList, ...
  khoa_colormap, lineStartPoints, lineEndPoints, smoothingFac, ...
  towerUTM_coor, isPlot_bathyWithLine)


%% -- plot selected straight line
isPlot_selectedStraightLine = false;
Nr = 512;    % number of points along selected line
lowest_depth = -60;
plot_selectedStraightLine(lineStartPoints, lineEndPoints, ...
  depth_interpolant, Nr, lowest_depth, khoa_colormap, levelList, ...
  towerUTM_coor, isPlot_selectedStraightLine);

%% -- save data to csv
isSaved = false;
saved_bathyData(X_selected, Y_selected, IM_selected, Nx_selected, Ny_selected, isSaved)

%% -- function declarations
function plot_selectedStraightLine(lineStartPoints, lineEndPoints, ...
  depth_interpolant, Nr, lowest_depth, user_colormap, levelList, ...
  towerUTM_coor, isPlot_selectedStraightLine)

  if isPlot_selectedStraightLine
    enoughGreen = [0, 127, 0]/255.;
    
    fig = figure('Color', 'w', 'Units', 'centimeters', 'Position', [7, 7, 16, 6]);
      ax = axes('Parent', fig, 'Position', [0.11, 0.25, 0.8, 0.69]);

      x_space = linspace(lineStartPoints(1), lineEndPoints(1), Nr);
      y_space = linspace(lineStartPoints(2), lineEndPoints(2), Nr);
      
      if abs(lineStartPoints(1) - lineEndPoints(1)) < 1e-6
        r_space = y_space';
      elseif abs(lineStartPoints(2) - lineEndPoints(2)) < 1e-6
        r_space = x_space';
      else
        r_space = vecnorm([x_space', y_space'] - [x_space(1), y_space(1)], 2, 2);
      end
      
      depth_along_line = depth_interpolant(x_space, y_space);

      plot(ax, r_space, depth_along_line, 'Color', 'r', 'LineWidth', 1.5);
      
      hold(ax, 'on');
      patch(ax, [r_space; flipud(r_space)], [depth_along_line'; lowest_depth*ones(size(depth_along_line'))], ...
        [depth_along_line'; flipud(depth_along_line')], 'EdgeColor', 'None');
      
      % -- add tower location
      hold(ax, 'on');
        if abs(lineStartPoints(1) - lineEndPoints(1)) < 1e-6
          r_tower = towerUTM_coor(2);
        elseif abs(lineStartPoints(2) - lineEndPoints(2)) < 1e-6
          r_tower = towerUTM_coor(1);
        else
          r_tower = norm([towerUTM_coor(1) - x_space(1), towerUTM_coor(2) - y_space(1)]);
        end
        
        depth_at_tower = depth_interpolant(towerUTM_coor(1), towerUTM_coor(2));
        plot(ax, [r_tower, r_tower], [depth_at_tower, 0], ...
          'Color', enoughGreen, 'LineWidth', 2);
      
      colormap(ax, user_colormap);
      caxis(ax, [min(levelList), max(levelList)]);
      
      grid(ax, 'on'); 
      set(ax, 'FontSize', 12, 'GridColor', [0, 0, 0, 0.2], 'Layer', 'top');
      
      xlabel(ax, 'along line [m]');
      ylabel(ax, 'depth [m]');
      
      xlim(ax, [0, r_space(end)]);
  
  end
end

function plot_bathyWithLine(X_selected, Y_selected, IM_selected, ...
  levelList, user_colormap, lineStartPoints, lineEndPoints, smoothingFac, ...
  towerUTM_coor, isPlot_bathyWithLine)
  
  enoughGreen = [0, 127, 0]/255.;

  if isPlot_bathyWithLine
    fig = figure('Color', 'w', 'Units', 'centimeters', 'Position', [5, 5, 13, 14]);
      ax = axes('Parent', fig, 'Position', [0.14, 0.11, 0.8, 0.8]);
      
      xSpace = X_selected(1, :);
      ySpace = Y_selected(:, 1);

      imagesc(ax, IM_selected, 'XData', xSpace, 'YData', ySpace);

      hold(ax, 'on')
        [contourData, contourHandler] = contour(ax, X_selected, Y_selected, IM_selected, ...
          'LevelList', levelList, 'ShowText', 'on');
         contourHandler.LineColor = 'k';
         contourHandler.LineStyle = '-';
         clabel(contourData, contourHandler, 'FontSize', 12, 'Color', 'k', ...
           'LabelSpacing', 300); 
         
      % -- add straight line
      hold(ax, 'on');
        plot(ax, [lineStartPoints(1), lineEndPoints(1)], ...
          [lineStartPoints(2), lineEndPoints(2)], 'Color', 'r', ...
          'LineWidth', 2);
       
      % -- add tower coordinate
      hold(ax, 'on');
        plot(ax, towerUTM_coor(1), towerUTM_coor(2), 'Marker', 'o', ...
          'MarkerFaceColor', 'w', 'LineWidth', 1.5, ...
          'MarkerEdgeColor', enoughGreen);
         
      cbarHandler = colorbar(ax);
      caxis(ax, [min(levelList), max(levelList)]);
      colormap(ax, user_colormap);

      ylabel(cbarHandler, 'depth [m]', 'FontSize', 12);
      set(cbarHandler, 'FontSize', 12);

      grid(ax, 'on'); 
      set(ax, 'YDir', 'normal', 'FontSize', 12, 'GridColor', [0, 0, 0, 0.2]);
      xlabel(ax, 'easting [m]');
      ylabel(ax, 'northing [m]');
      
      if smoothingFac > 0
        title(ax, sprintf('smoothingFac = %g', smoothingFac));
      else
        title(ax, sprintf('no smoothingFac'));
      end
        
        
      daspect(ax, [1, 1, 1]);
  end
end

function plot_3Dbathy(X_selected, Y_selected, IM_selected, user_colormap, ...
  towerUTM_coor, depth_interpolant, isPlot_3Dbathy)
  
  if isPlot_3Dbathy
    enoughGreen = [0, 127, 0]/255.;
    
    fig = figure('Color', 'w', 'Units', 'centimeters', 'Position', [20, 3, 20, 12]);
      ax = axes('Parent', fig);

      surf_handler = surf(ax, X_selected, Y_selected, IM_selected);

      surf_handler.EdgeColor = 'none';

      hold(ax, 'on');
        depthAtTower = depth_interpolant(towerUTM_coor(1), towerUTM_coor(2));
        plot3(ax, [towerUTM_coor(1), towerUTM_coor(1)], ...
          [towerUTM_coor(2), towerUTM_coor(2)], ...
          [depthAtTower, 0], 'Color', enoughGreen, 'LineWidth', 2);
      
      colormap(ax, user_colormap);
      
      material(surf_handler, 'dull');
      lightangle(ax, 45, 60);
      xlabel(ax, 'easting [m]');
      ylabel(ax, 'northing [m]');
      
      set(ax, 'FontSize', 12);

      daspect(ax, [1, 1, 0.1]);
  end
end

function plot_bathy(X_selected, Y_selected, IM_selected, levelList, ...
  user_colormap, smoothingFac, isPlot_bathy)

  if isPlot_bathy
    fig = figure('Color', 'w', 'Units', 'centimeters', 'Position', [3, 3, 16, 14]);
      ax = axes('Parent', fig, 'Position', [0.14, 0.11, 0.8, 0.8]);

      xSpace = X_selected(1, :);
      ySpace = Y_selected(:, 1);

      imagesc(ax, IM_selected, 'XData', xSpace, 'YData', ySpace);

      hold(ax, 'on')
        [contourData, contourHandler] = contour(ax, X_selected, Y_selected, IM_selected, ...
          'LevelList', levelList, 'ShowText', 'on');
         contourHandler.LineColor = 'k';
         contourHandler.LineStyle = '-';
         clabel(contourData, contourHandler, 'FontSize', 12, 'Color', 'k', ...
           'LabelSpacing', 300);

      cbarHandler = colorbar(ax);
      caxis(ax, [min(levelList), max(levelList)]);
      colormap(ax, user_colormap);

      ylabel(cbarHandler, 'depth [m]', 'FontSize', 12);
      set(cbarHandler, 'FontSize', 12);

      grid(ax, 'on'); 
      set(ax, 'YDir', 'normal', 'FontSize', 12, 'GridColor', [0, 0, 0, 0.2]);
      xlabel(ax, 'easting [m]');
      ylabel(ax, 'northing [m]');
      
      if smoothingFac > 0
        title(ax, sprintf('smoothingFac = %g', smoothingFac));
      else
        title(ax, sprintf('no smoothingFac'));
      end
        
        
      daspect(ax, [1, 1, 1]);
  end
end

function saved_bathyData(X_selected, Y_selected, IM_selected, ...
  Nx_selected, Ny_selected, isSaved)

  if isSaved
    filenameOutput = sprintf('Ieodo_depth_Nx_%d_Ny_%d.txt', Nx_selected, Ny_selected);
    
    fprintf('   Saved %s\n', filenameOutput);
    % Create two columns of data
    % Create a table with the data and variable names
    T = table(X_selected(:), Y_selected(:), IM_selected(:), 'VariableNames', { 'UTM_x', 'UTM_y', 'depth'} );
    % Write data to text file
    writetable(T, filenameOutput);
  end
end

function [X_selected, Y_selected, IM_selected, F_interpolant, towerUTM_coor] ...
  = transformToUTM(depthData, zeroRefCoor, selected_xRange, selected_yRange, ...
      Nx_selected, Ny_selected, smoothingFac, rotateAngle, towerCoor)

  [UTM_x, UTM_y, ~] = deg2utm(depthData.lat', depthData.lon');
  
  refPos = dms2degrees(zeroRefCoor);
  [Ref_x, Ref_y, ~] = deg2utm(refPos(1), refPos(2));
  
  UTM_x = UTM_x - Ref_x;
  UTM_y = UTM_y - Ref_y;
  
  towerCoor = dms2degrees(towerCoor);
  [towerUTM_coor(1), towerUTM_coor(2)] = deg2utm(towerCoor(1), towerCoor(2));
  towerUTM_coor(1) = towerUTM_coor(1) - Ref_x;
  towerUTM_coor(2) = towerUTM_coor(2) - Ref_y;
  
  % Rotate the data
  if abs(rotateAngle) > 1e-6
    rotMatrix = [cosd(rotateAngle), -sind(rotateAngle);
                 sind(rotateAngle),  cosd(rotateAngle)];
    UTM_xy = rotMatrix * [UTM_x' ; UTM_y'];
    UTM_x = UTM_xy(1,:)';
    UTM_y = UTM_xy(2,:)';
    
    towerUTM_coor = rotMatrix * [towerUTM_coor(1); towerUTM_coor(2)];
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

  
  %diffX = X_space(2) - X_space(1);
  %diffY = Y_space(2) - Y_space(1);
  
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