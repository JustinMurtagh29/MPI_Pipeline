function [ interfaces, raw ] = interfacesFromTrData( data )
%INTERFACESFROMTRDATA Get the interface struct from the training data
%struct.
% INPUT data: struct
%           Data struct from SynEM training cubes.
% OUTPUT interfaces: struct
%           Struct containing the fields 'surface' and 'subseg' required
%           for feature calculation.
%        raw: 3d uint8
%           The training cube raw data.
% Author: Benedikt Staffler <benedikt.staffler@brain.mpg.de>

raw = data.raw;
interfaces.surface = data.interfaceSurfaceList;
interfaces.subseg = data.subsegmentsList;

end

