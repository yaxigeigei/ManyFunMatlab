function  [ files, n ] = findfiles_all_subpath( filespec )
%FINDFILES_ALL_SUBPATH 
% [files n]= findfiles_all_subpath( filespec )

cmd = [ 'dir "' filespec '" /s /b'  ];
[ ~, filenames ] = dos(cmd);
files = textscan(filenames, '%s', 'Delimiter', '\n');
files = files{1};
n = length(files);

end