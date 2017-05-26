%load('skel')
%addpath('/gaba/u/kboerg/code/skel2graph3d-matlab/')
w = size(skel,1);
l = size(skel,2);
h = size(skel,3);
tic
% initial step: condense, convert to voxels and back, detect cells
THR = 0;
[~,node,link] = Skel2Graph3D(skel,THR);
toc
% total length of network
wl = sum(cellfun('length',{node.links}));
tic
skel2 = Graph2Skel3D(node,link,w,l,h);
if wl
[~,node2,link2] = Skel2Graph3D(skel2,THR);
else
    node2=node;
end
toc
% calculate new total length of network
wl_new = sum(cellfun('length',{node2.links}));

% iterate the same steps until network length changed by less than 0.5%
while(wl_new~=wl)
    'once again'
    wl = wl_new;
    if (wl==0)
        break;
    end
     skel2 = Graph2Skel3D(node2,link2,w,l,h);
     [A2,node2,link2] = Skel2Graph3D(skel2,THR);

     wl_new = sum(cellfun('length',{node2.links}));

end;
