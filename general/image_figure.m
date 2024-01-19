function fig = image_figure()
% This is the machine-generated representation of a Handle Graphics object
% and its children.  Note that handle values may change when these objects
% are re-created. This may cause problems with any callbacks written to
% depend on the value of the handle at the time the object was saved.
%
% To reopen this object, just type the name of the M-file at the MATLAB
% prompt. The M-file and its associated MAT-file must be on your path.

load image_figure

h0 = figure('Color',[0.8 0.8 0.8], ...
	'Colormap',mat0, ...
	'MenuBar','none', ...
	'NumberTitle','off', ...
	'PointerShapeCData',mat1, ...
	'Position',[73 27 896 672], ...
    'Resize','on', ...
	'Tag','image_figure');
h1 = uimenu('Parent',h0, ...
	'Label','File', ...
	'Tag','image_fig_file_menu');
h2 = uimenu('Parent',h1, ...
	'Callback','imageCbf save_image', ...
	'Label','Save Figure', ...
	'Tag','save_image_menu');
h2 = uimenu('Parent',h1, ...
	'Label','Export', ...
	'Tag','export_image_menu');
h3 = uimenu('Parent',h2, ...
	'Callback','imageCbf export_jpg', ...
	'Label','*.jpg', ...
	'Tag','jpg_menu');
h3 = uimenu('Parent',h2, ...
	'Callback','imageCbf export_bmp', ...
	'Label','*.bmp', ...
	'Tag','bmp_menu');
if nargout > 0, fig = h0; end
