
function onlyFwdPass3DonKnossosFolder(cnetLocation, gpuSwitch, input, result, bbox, normFunc)

load(cnetLocation,'cnet');
cnet.run.savingPath = [fileparts(cnetLocation) '/'];

fwdPass3Dfaster(cnet,input,result,bbox,normFunc);


end


