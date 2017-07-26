function seg2stack(p)


stack.folder = fullfile(p.saveFolder,'globalSeg');

f = load(fullfile(p.saveFolder,'numEl.mat'));

stack.section = struct('sectionID','section1','bbox', p.bbox-1, 'resolutions',1);
stack.layer = struct('largestValue', f.numElTotalUpper(end,end,end),'typ','segmentation','class','uint32');   % largest value might be wrong

WK.prepareStack(stack)

