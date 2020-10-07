function appRegion = applicationRegion( locations, faceNormals )
%APPLICATIONREGION Construct application region struct for user function
% This undocumented function may be changed or removed in a future release.

%       Copyright 2014 The MathWorks, Inc.

appRegion.x = locations(1,:);
appRegion.y = locations(2,:);
appRegion.z = locations(3,:);
appRegion.nx = faceNormals(1,:);
appRegion.ny = faceNormals(2,:);
appRegion.nz = faceNormals(3,:);
end

