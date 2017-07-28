function [ data, metadata ] = calculateContacts( knossos_conf, segmentation_conf, bbox, rinclude, useInGUI )
%CALCULATECONTACTS Calculate contact surfaces in a given volume.

fprintf('%s INFO loading raw data and segmentation\n',datestr(now,'yyyy-mm-dd HH:MM:SS,FFF'));
raw = SynEM.Paper.TrainingSet.readKnossosRoi(knossos_conf.parfolder,knossos_conf.fileprefix,bbox,'uint8',[128 128 128],0);
if all(~raw)
    error('The raw file read from %s is empty.',knossos_conf.parfolder);
end
walls = SynEM.Paper.TrainingSet.readKnossosRoi(segmentation_conf.parfolder,segmentation_conf.fileprefix,bbox,'uint16',[128 128 128],10);
if all(~walls)
    error('The walls file read form %s is empty.',segmentation_conf.parfolder);
end
segments = bwlabeln(walls,26);
rinclude = sort(rinclude);

%calculate contacts
[contactSurfaceList, subsegmentsList, neighborIDs ] = SynEM.Paper.TrainingSet.contactCalculationPipeline( segments,~walls,rinclude,bbox,useInGUI );

%create data for gui
data.raw = raw;
if max(segments(:)) < 65537
    data.segments = uint16(segments);
    data.neighborIDs = uint16(neighborIDs);
else
    data.segments = uint32(segments);
    data.neighborIDs = uint32(neighborIDs);
end
data.contactSurfaceList = contactSurfaceList;
data.subsegmentsList = subsegmentsList;


%create metadata
metadata.experiment = knossos_conf.parfolder;
metadata.boundingBox = bbox;
if useInGUI
    metadata.innerBoundingBox = [100,bbox(1,2)-bbox(1,1)-99;100,bbox(2,2)-bbox(2,1)-99;40,bbox(3,2)-bbox(3,1)-39];
else
    metadata.innerBoundingBox = [40,bbox(1,2)-bbox(1,1)-39;40,bbox(2,2)-bbox(2,1)-39;14,bbox(3,2)-bbox(3,1)-13];
end

metadata.includedSubsegmentDistance = rinclude;

fprintf('%s INFO finished contact calculation\n',datestr(now,'yyyy-mm-dd HH:MM:SS,FFF'));

end

