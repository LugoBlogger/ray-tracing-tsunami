% TraceContour_netCDF.m
%  A script to do wave ray tracing with Snell's law.
%  This script requires `deg2utm` conversion and `InterX.m` for finding
%  intersection. This scripts is an improvement version of `TraceBathy_Glagah.m`
%  and with input file format netCDF (.nc)
%
% Author: 
% - Durra Hadri (Durra.Saputera@uib.no)
% - Henokh Lugo Hariyanto (henokh.lugo.h@mail.ugm.ac.id)
%
% Logs
% [2023/10/07]   
% - First release
% - In user-defined input we use DMS unit instead of the decimal to
%   increase the accuracy of the location 
%
% [2023/10/08]
% - Define current level for NaN to a high number greater than
%   max(Vg_bat(:)) to emanate ray from inland.
% - Incoming ray to inland are reflected (not refracted)


clear;

%% -- user-defined input
toolboxPath = 'deg_utm_conversion/';
addpath(genpath(toolboxPath));

filenames = {'BATNAS_cropped.nc'};
stackDir = '';  % if we have more than one .tif images

% use degrees2dms() if you are given the source in the decimal form
refPos = [-6, 6, 7.2;      % This is the source of the Tsunami
          105, 25, 22.8];   % (lat; lon);
N_init_points = 20;
        
% cut area around the reference point
% bboxAroundRefPos = [];
% bboxAroundRefPos = [105.4229722, 105.4230278, ...
%                     -6.1019722, -6.1020278];    % [xMin, xMax, yMin, yMax]  (in degrees)        
bboxAroundRefPos = [10000, 10000, ...
                    10000, 10000];    % [xShiftLeft, xShiftRight, yShiftDown, yShiftUp]  (in degrees)   
rot_angle = 0;
% define border of ray tracing (after the border, the tracing will stop)
X_border = [-8000, 8000];
Y_border = [-8000, 8000];

selectedTp = 600;
smoothingFac = 2;   % bathymetry smoothing factor (-1 is no smoothing)

levelLineBathyPlot = -2200:100:0;

showBathyPlot = 'off';
showVgPlot = 'off';


%% -- read coordinate x, y (lon, lat)
fprintf("   %s\n", 'Read coordinate from netCDF file');
[IM_combine, XX, YY, refPos_coord] = getBathyData(filenames, refPos, ...
  stackDir, smoothingFac);
x = XX(:,1);
y = YY(1,:);


%% -- convert coordinate to UTM
fprintf("   %s\n", 'Convert (lon, lat) to (UTMx, UTMy)');
[Nx ,Ny] = size(IM_combine);

[Refx, Refy, ~] = deg2utm(refPos_coord(2), refPos_coord(1));
[UTx, UTy, ~] = deg2utm(YY(:), XX(:));

UTx = reshape(UTx, Nx, Ny);
UTy = reshape(UTy, Nx, Ny);

UTx = UTx-Refx;
UTy = UTy-Refy;

xUTM = UTx(:,1);
yUTM = UTy(1,:);

if ~isempty(bboxAroundRefPos)
  idx = logical(    (xUTM >= - bboxAroundRefPos(1)) ...
                 .* (xUTM <=  bboxAroundRefPos(2)));
  idy = logical(    (yUTM >= - bboxAroundRefPos(3)) ...
                 .* (yUTM <=  bboxAroundRefPos(4)));

  UTx = UTx(idx, idy);
  UTy = UTy(idx, idy);
  IM_combine = IM_combine(idx, idy);
  
  [Nx ,Ny] = size(IM_combine);
end      

Xsave = UTx(:);
Ysave = UTy(:);


%% -- plot bathy
fprintf("   %s [%s]\n", 'Plot Bathy', showBathyPlot);
fig = figure('Color', 'w', 'Visible', showBathyPlot);
  ax = axes('Parent', fig);
  
  contour(ax, UTx, UTy, IM_combine, 'LevelList', levelLineBathyPlot);
  
  hold(ax, 'on');
    plot(ax, 0, 0, '*', 'MarkerSize', 16);
    caxis(ax, [-2500, 500]);
  hold(ax, 'off');
  title(ax, 'Plot Bathy');

% IM_save = IM_combine(:);
% IM_save(IM_save<-50) = -50;
% dlmwrite('bathy2D_BATNAS_cropped.txt',[Xsave,Ysave,IM_save])


%% -- compute input data for ray tracing (wavenumber, Vg_bat, X_process, Y_process)
fprintf("   %s\n", 'Calcuate ray tracing params');
% calculating Vg in the domain
g = 9.81;      % acceleration of gravity
freq = 2*pi/selectedTp;   

% (Guo, 2002) - `Simple and explicit solution of wave dispersion equation`
beta = 2.4908;   % transitional shape parameter
wavenumber = (freq.^2/g) .* (1 - exp(-(freq.*sqrt(abs(IM_combine)/g)).^beta)) .^ (-1/beta);

positiveDepth = -IM_combine;
positiveDepth(positiveDepth <= 0) = 0;

% (Dingemans, 1994) `Water wave propagation over uneven bottoms` Eqs. 2.29a and 2.29c
Vp_bat = sqrt(g * tanh(wavenumber .* positiveDepth) ./ wavenumber);
Vg_bat = 0.5 * Vp_bat .* (1 + 2*wavenumber .* positiveDepth ./ sinh(2*wavenumber .* positiveDepth));

% numOfLevel = 100;     % define the group velocity interval here
numOfLevel = 20;     % define the group velocity interval here
levelLineVg_bat = linspace(0, max(Vg_bat(:)), numOfLevel); 
% levelLineVg_bat = linspace(5, 50, 50); 

%Add coordinate rotation
XY_Rot = [cosd(rot_angle), -sind(rot_angle);
          sind(rot_angle) cosd(rot_angle)] ...
          * [Xsave';
             Ysave'];
Xprocess = reshape(XY_Rot(1,:), Nx, Ny);
Yprocess = reshape(XY_Rot(2,:), Nx, Ny);

% contour(Xprocess,Yprocess,reshape(IM_save,Nx,Ny),'ShowText','on','LevelList',[-100:10:-10,-5,-2.5])

dVLevel = levelLineVg_bat(2) - levelLineVg_bat(1);


%% -- plot group speed
fprintf("   %s [%s]\n", 'Plot Group Speed', showVgPlot);
fig = figure('Color', 'w', 'Visible', showVgPlot);
  ax = axes('Parent', fig);
  cmatrixVg_bat = contour(ax, Xprocess, Yprocess, Vg_bat, 'LevelList', levelLineVg_bat);
  title(ax, 'Plot Vg\_bat')

  hold(ax, 'on');
    plot(ax, 0, 0, '*', 'MarkerSize', 16);
  hold(ax, 'off');
  
clineVg_bat = getContourLine(cmatrixVg_bat);

% replot clineVg_bat to verify getContourLine function
% you have to copy manually the contour line made by this plot to the
% previous contour plot.
% fig = figure('Color', 'w');
%   ax = axes('Parent', fig);
%  
%   for i = 1:size(clineVg_bat, 2)
%     if i ~= 1
%       hold(ax, 'on');
%     end
%     xData = clineVg_bat(i).xy(1,:);
%     yData = clineVg_bat(i).xy(2,:);
%     plot(ax, xData, yData);
%   end


%% -- Start the animation of ray tracing
fig = figure('Color', 'w', 'Units', 'centimeters', 'Position', [3, 3, 16, 16]);

  ax = axes('Parent', fig);
  
  contour(ax, Xprocess, Yprocess, IM_combine, ...
    'LevelList', levelLineBathyPlot, 'ShowText','on');
  
  hold(ax, 'on')
  contour(ax, Xprocess, Yprocess, Vg_bat);
  
  caxis(ax, [min(levelLineBathyPlot), max(levelLineBathyPlot)]);
  
  cbarHandler = colorbar(ax);
  ylabel(cbarHandler, 'Group velocity');
  set(cbarHandler, 'FontSize', 12);
  
  set(ax, 'FontSize', 12);
  xlabel(ax, 'easting [m]');
  ylabel(ax, 'northing [m]');
  hold(ax, 'on');


XY0 = [zeros(1, N_init_points);
       zeros(1, N_init_points)];
contour_flag = 0;

R = [cosd(90), -sind(90); 
     sind(90), cosd(90)];
   
Realized_lvl = zeros(size(clineVg_bat));
for ii = 1:length(clineVg_bat)
  Realized_lvl(ii) = clineVg_bat(ii).Level;
end
Realized_lvl = unique(Realized_lvl);
   
F_interpolant = scatteredInterpolant(Xprocess(:), Yprocess(:), double(Vg_bat(:)));

for idx_XY = 1:size(XY0, 2)   % iterate over number of starting points
    
  % define direction of the vector at starting point
  % [Xstart, Xend, Ystart, Yend]
  dtheta = -2.0*pi/(N_init_points);
  radius = 5000;   % <--- set this
  xy_line = [XY0(1, idx_XY), radius*cos((idx_XY-1)*dtheta);
             XY0(2, idx_XY), radius*sin((idx_XY-1)*dtheta)];     

  Txy = 0;                            % Rotation Degree minus means clockwise
  Rxy = [cosd(Txy), -sind(Txy);
         sind(Txy), cosd(Txy)];
  dvector = xy_line(:,2) - xy_line(:,1);
  xy_line(:,2) = xy_line(:,1) + Rxy*dvector;

  % level at the starting point
  curr_level = F_interpolant(xy_line(1,1), xy_line(2,1)); 
  if isnan(curr_level)  % if the ray emanates from the inland, set the source to be high enough.
    curr_level = max(Realized_lvl) + 1;
  end
    
  % interpolate curr_level based on Realized_lvl
  if curr_level >= max(Realized_lvl)
    curr_level = max(Realized_lvl) + dVLevel/2;
  elseif curr_level < min(Realized_lvl)
    curr_level = min(Realized_lvl) - dVLevel/2;
  else
    id1 = find(Realized_lvl <= curr_level, 1, 'Last');
    id2 = find(Realized_lvl >= curr_level, 1, 'First');
    curr_level = mean(Realized_lvl([id1, id2]));
  end
    
  jj = 1;    % iterator for storing each intersection between path and contour line
  while true
    Pmin = [];
    idmin = [];   % this will store which index corresponds to cline_bat.Level

    for ii = 1:length(clineVg_bat)  % iterate over all contour line to find any intersection
      P = InterX(xy_line(:,jj:jj+1), clineVg_bat(ii).xy);

      if ~isempty(P) 
        pos_diff = P-xy_line(:,jj);  
        norm_diff = sqrt(pos_diff(1,:).^2 + pos_diff(2,:).^2);

        % -- if there are more than one intersection point, select the 
        %    closest to xy_line(:, jj) 
        if size(P, 2) > 1
          min_norm_diff = norm_diff(norm_diff > 1e-6);
          min_norm_diff = min(min_norm_diff(:));
          P = P(:, abs(norm_diff - min_norm_diff) < 1e-6);
          norm_diff = min(norm_diff(norm_diff > 1e-6));
        end

        if ~isempty(Pmin)
          if norm_diff <= norm(Pmin - xy_line(:, jj)) && norm_diff > 1e-3
            Pmin = P;
            idmin = ii;
          end
        else
          if norm_diff > 1e-6
            Pmin = P;
            idmin = ii;
          end
        end
      end
    end
    
    contour_flag=1;
    %fprintf('size(Pmin) = %d, %d\n', size(Pmin));
    if ~isempty(Pmin)
      fprintf('   countourLev at Pmin: %g\n', F_interpolant(Pmin(1,1), Pmin(2,1)));
    else
      fprintf('   countourLev at Pmin: Pmin is empty\n');
    end

    if ~isempty(Pmin)
      xy_line(:, jj+1) = Pmin;                % set xy_line(:, jj+1) to the nearest intersection point
      prev_dir = xy_line(:, jj+1) - xy_line(:, jj);
      prev_dir = prev_dir / norm(prev_dir);   % normalized vector from xy_line(:, jj) to xy_line(:, jj+1)

      dxy_line = clineVg_bat(idmin).xy - Pmin;
      norm_dxy_line = sqrt(dxy_line(1,:).^2 + dxy_line(2,:).^2);
      
      % -- iterate over dxy_line to get min_dir?
      for zz = 1:size(clineVg_bat(idmin).xy,2) - 1
        min_dir = dxy_line(:, zz) / norm_dxy_line(zz);
        plus_dir = dxy_line(:, zz+1)/norm_dxy_line(zz + 1);
        dir_diff = plus_dir + min_dir;

        if round(norm(dir_diff)) == 0
          break
        end
      end

      str1 = sprintf('current level = %.2f, next level = %.2f\n', ...
        curr_level, clineVg_bat(idmin).Level);
      lvl_sgn = sign(curr_level - clineVg_bat(idmin).Level);  % +1 uphill
                                                              % -1 downhill

      % -- dVLevel is always positive; uphill -> v_ratio < 1; downhill ->
      %    v_ratio > 1; v_ratio represents ratio between refracted index to
      %    incident index. Uphill motion means moving to higher refractive
      %    index. Downhill motion means moving to lower refractive index.
      v_ratio = (curr_level - dVLevel*lvl_sgn) / curr_level;  

      curr_level = curr_level - dVLevel*lvl_sgn;

      %perpendicular vector to the contour line
      pp_vector = R * min_dir;
      pp_vector = (prev_dir'*pp_vector) * pp_vector;   % this makes sure that the pp_vector
                                                       % never reverses
                                                       % its motion

      pp_vector = pp_vector/norm(pp_vector);

      sin_t = cross([pp_vector;0], [prev_dir;0]);   % incident angle
      sin_t1 = sin_t;
      sin_t = sin_t(3)*v_ratio;                     % refracted angle

      %larger than critical angle 
      %(maybe wrong, but this is my best approximate)
      if abs(sin_t)>1
        % sin_t= sign(sin_t);
        sin_t= sin_t1(3);
      end
      cos_t = sqrt(1 - sin_t^2);
      title(ax, [str1, '', sprintf('v\\_ratio = %.2f, theta1 = %.2f', ...
        v_ratio, rad2deg(asin(sin_t1(3))))])

      R_pp = [cos_t, -sin_t; 
              sin_t cos_t];
      pp_vector = R_pp * pp_vector;  % rotate pp_vector to the direction guided by
                                     % v_ratio 

      if abs(sin_t * v_ratio)>1
        RR = sign(sin_t) *[sin_t, -cos_t; 
                           cos_t, sin_t];
        pp_vector = (RR*RR) * pp_vector;
      end

      % -- store the next 
      if pp_vector(1) > 0
        xy_line(:, jj+2) = [X_border(2);
                            Pmin(2) + (X_border(2)-Pmin(1))*pp_vector(2)/pp_vector(1)];
      elseif pp_vector(1) < 0
        xy_line(:, jj+2) = [X_border(1);
                            Pmin(2) + (X_border(1)-Pmin(1))*pp_vector(2)/pp_vector(1)];
      elseif abs(pp_vector(1)) < 1e-6
        if pp_vector(2) > 0
          xy_line(:, jj+2) = [Pmin(1);
                              Y_border(2)];
        else
          xy_line(:,jj+2) = [Pmin(1);
                             Y_border(1)];
        end
      end
      
      axis(ax, 'equal');

      xlim(ax, X_border);
      ylim(ax, Y_border);

      fprintf('   jj = %i\n',jj)
      fprintf('   Pmin: \n');
      disp(Pmin);
      
      if jj == 1
        linePlotHandler = plot(ax, xy_line(1, 1:end-1), xy_line(2, 1:end-1), ...
          'r', 'LineWidth', 3);
        pointPlotHandler = plot(ax, Pmin(1), Pmin(2), 'or');
        quivPlotHandler = quiver(ax, Pmin(1), Pmin(2), ...
          pp_vector(1)/norm(pp_vector), pp_vector(2)/norm(pp_vector), ...
          'AutoScale', 'on', 'AutoScaleFactor', 500, 'MaxHeadSize', 8, ...
          'LineWidth', 1, 'Color', 'r');
      else
        linePlotHandler.XData = xy_line(1, 1:end-1);
        linePlotHandler.YData = xy_line(2, 1:end-1);
        pointPlotHandler.XData = Pmin(1);
        pointPlotHandler.YData = Pmin(2);
        quivPlotHandler.XData = Pmin(1);
        quivPlotHandler.YData = Pmin(2);
        quivPlotHandler.UData = pp_vector(1);
        quivPlotHandler.VData = pp_vector(2);
      end

      if Pmin(1) <= X_border(1) || Pmin(1) >= X_border(2) ...
          || Pmin(2) <= Y_border(1) || Pmin(2)>= Y_border(2)
        
        fprintf('   X0,Y0 = %.2f, %.2f tracing is done\n',XY0(:,idx_XY));
        break
      end
      
      % hit the new inland (stop the ray tracing)
      radius_next = 100;  % <-- set this
      test_next_nearest = xy_line(:,jj+1) + radius_next*pp_vector;
      if isnan(F_interpolant(test_next_nearest(1), test_next_nearest(2)))
        fprintf('   X0,Y0 = %.2f, %.2f tracing is done\n',XY0(:,idx_XY));
        break
      end
      
      pause(0.001)

      jj = jj + 1;
      
    else
      fprintf('   (X0, Y0) = (%.2f, %.2f) tracing is done\n', XY0(:, idx_XY));
      break
    end
    

  end

  linePlotHandler.Color = [1, 0, 1];
  linePlotHandler.LineWidth = 1;
  pointPlotHandler.Color = [1, 0, 1];
  quivPlotHandler.Color = [1, 0, 1];
  
  fprintf('   xy_line: \n');
  disp(xy_line);
end


%% -- function declarations
function [IM_combine, XX, YY, refPos_coord] = getBathyData(filenames, ...
  refPos, stackDir, smoothingFac)

  numOfFiles = size(filenames, 2);
  
  switch numOfFiles
    case 2   % we stack horizontally
      switch stackDir
        case 'horizontal'
          % not implemented
      end
      
    case 1
      x = ncread(filenames{1}, 'lon');
      y = ncread(filenames{1}, 'lat');
      
      [XX, YY] = meshgrid(x, y);

      IM_combine = ncread(filenames{1}, 'Band1');
      IM_combine = double(IM_combine); 
      %IM_combine = IM_combine';
      
      if smoothingFac > 0
        IM_combine = imgaussfilt(IM_combine, smoothingFac);
      end
        
      refPos = dms2degrees(refPos);
      refPos_coord = [refPos(2), refPos(1)];
      
      XX = XX';
      YY = YY';
  end

end

function cline = getContourLine(cmatrix)
  % An improvement version of getContourLine by finding first the number 
  % of lines.
  % See `RadarImage.m: get_zero_crossing()
  
  cmatrix_colLength = size(cmatrix, 2);
  lineCounter = 1;
  nPointCounter = 1;
  
  % scan number of lines
  while true
    nPoint = cmatrix(2, nPointCounter);
    nPointCounter = nPointCounter + nPoint + 1;
    
    if nPointCounter > cmatrix_colLength
      break
    end
    
    lineCounter = lineCounter + 1;
  end
  
  % get cline that consists of the columns: 1) line level, 2) num. of
  % points, 3) xy-coordinates
  cline = struct('Level', repmat({[]}, 1, lineCounter), ...
                 'No', repmat({[]}, 1, lineCounter), ...
                 'xy', repmat({[]}, 1, lineCounter));
  lineCounter = 1;
  nPointCounter = 1;
  while true
    nPoint = cmatrix(2, nPointCounter);
    cline(lineCounter).Level = cmatrix(1, nPointCounter);
    cline(lineCounter).No = nPoint;
    cline(lineCounter).xy = cmatrix(:, nPointCounter+1:nPointCounter+nPoint);
    
    nPointCounter = nPointCounter + nPoint + 1;
    
    if nPointCounter > cmatrix_colLength
      break
    else
      lineCounter = lineCounter + 1;
    end
     
  end

end