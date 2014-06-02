function [imfeat] = filter3d(parameter,I, idx)

%filter3d: gives back the filtered image(s) imfeats of the input image I
%the applied filter is specified by type and the siz parameter (two in the
%case of the structure tensor), the filters are created with fspecial3
%
%Currently used 'gaussiansmoothedgradmagnitude' idx(2), 'sortedeigenvalueshessian' idx(1,4), 'intensitygaussiansmoothed' idx(3),
%'laplaceofgaussian' idx(5), 'differenceofgaussians' idx(6)

% profile on;
% 
type = parameter.filter{idx}{1};
sizs = parameter.filter{idx}{2};
siz2s = parameter.filter{idx}{3};
if ~isa(I, 'single')  %uniform data type
    I=single(I);
end

imfeats=cell(length(sizs),1);
for u=1:length(sizs)
	siz=sizs(u);	
	switch type
	    case 'intensitygaussiansmoothed'
             h=fspecial3('gaussianFWHMhalfSize',siz);             
             imfeatx = convn(I,h{1},'same');
             imfeaty = convn(imfeatx,h{2},'same');
             imfeat = convn(imfeaty,h{3},'same');
             clear imfeaty imfeatx
            
	    case 'gaussiansmoothedgradmagnitude'
	        hx=fspecial3('gaussgradient1',siz,1);
	        hy=fspecial3('gaussgradient1',siz,2);
	        hz=fspecial3('gaussgradient1',siz,3);
	        
	        Ix_x=convn(I,hx{1},'same');
            Ix_y=convn(Ix_x,hx{2},'same');
            Ix_z=convn(Ix_y,hx{3},'same');
            clear Ix_x Ix_y Ix_0
            Iy_x=convn(I,hy{1},'same');
            Iy_y=convn(Iy_x,hy{2},'same');
            Iy_z=convn(Iy_y,hy{3},'same');
            clear Iy_x Iy_y Iy_0
            Iz_x=convn(I,hz{1},'same');
            Iz_y=convn(Iz_x,hz{2},'same');
            Iz_z=convn(Iz_y,hz{3},'same');
            
	        imfeat=sqrt(Ix_z.^2+Iy_z.^2+Iz_z.^2);
            clear Iz_x Iz_y Iz_0 Iz_z
          
	    case 'laplaceofgaussian'   
		clear h  
	        h=fspecial3('log',siz); 
               	imfeat_x=convn(I,h{1}, 'same');
		imfeat_y=convn(imfeat_x,h{2},'same');
		imfeat=convn(imfeat_y,h{3},'same');
             
	    case 'differenceofgaussians'
	         h=fspecial3('dog',siz);
             
             imfeatx = convn(I,h{1},'same');
             imfeaty = convn(imfeatx,h{2},'same');
             imfeat1 = convn(imfeaty,h{3},'same');
             
             imfeatx2 = convn(I,h{4},'same');
             imfeaty2 = convn(imfeatx2,h{5},'same');
             imfeat2 = convn(imfeaty2,h{6},'same');
            
             imfeat = imfeat1 - imfeat2;
	        
	    case 'sortedeigenvalueshessian'
	        
	        %start by creating 6 filter for each entry in the Hessian
	        hxx=fspecial3('gaussgradient2',siz,11);
	        hyy=fspecial3('gaussgradient2',siz,22);
	        hzz=fspecial3('gaussgradient2',siz,33);
	        hxy=fspecial3('gaussgradient2',siz,12);
	        hxz=fspecial3('gaussgradient2',siz,13);
	        hyz=fspecial3('gaussgradient2',siz,23);
	        
	        %apply the filter to the image, ignoring boundaries, i.e. intensity outside
	        %the image is set to zero
	        Ixx_x=convn(I,hxx{1},'same');
            Ixx_y=convn(Ixx_x,hxx{2},'same');
            Ixx=convn(Ixx_y,hxx{3},'same');
            
            Iyy_x=convn(I,hyy{1},'same');
            Iyy_y=convn(Iyy_x,hyy{2},'same');
            Iyy=convn(Iyy_y,hyy{3},'same');
            clear Iyy_x Iyy_y Ixx_x Ixx_y
            
            Izz_x=convn(I,hzz{1},'same');
            Izz_y=convn(Izz_x,hzz{2},'same');
            Izz=convn(Izz_y,hzz{3},'same');
            
            Ixy_x=convn(I,hxy{1},'same');
            Ixy_y=convn(Ixy_x,hxy{2},'same');
            Ixy=convn(Ixy_y,hxy{3},'same');
            clear Izz_x Izz_y Ixy_x Ixy_y
            
            Ixz_x=convn(I,hxz{1},'same');
            Ixz_y=convn(Ixz_x,hxz{2},'same');
            Ixz=convn(Ixz_y,hxz{3},'same');
            
            Iyz_x=convn(I,hyz{1},'same');
            Iyz_y=convn(Iyz_x,hyz{2},'same');
            Iyz=convn(Iyz_y,hyz{3},'same');
            clear Ixz_x Ixz_y Iyz_x Iyz_y
 	       
	        %Calculate eigenvectors and eigenvalues of the Hessian for each point
	        [a,b,c]=size(I);
	        imfeat=cell(3,1);
	        newSize = [1 1 a*b*c];
	        Ieigen= eig3([reshape(Ixx, newSize) reshape(Ixy, newSize) reshape(Ixz, newSize); ...
	            reshape(Ixy, newSize) reshape(Iyy, newSize) reshape(Iyz, newSize); ...
	            reshape(Ixz, newSize) reshape(Iyz, newSize) reshape(Izz, newSize)]);
	        Ieigen = sort(Ieigen,2);
	        imfeat{1}=reshape(Ieigen(:,1), [a b c]);
	        imfeat{2}=reshape(Ieigen(:,2), [a b c]);
      	 	imfeat{3}=reshape(Ieigen(:,3), [a b c]); 
        
        	% normalize gradient
        	maxValue = max(imfeat{3}(:));
        	imfeat{1} = imfeat{1} ./ maxValue;
        	imfeat{2} = imfeat{2} ./ maxValue;
        	imfeat{3} = imfeat{3} ./ maxValue;  

	    case 'sortedeigenvaluesstructure'
        	
        	if isempty(siz2s)
        	    error('The structure tensor needs two size parameters!');
        	else
        	    siz2=siz2s(u);
        	    %calculate gradients with fspecial3 using size siz2
        	    hx=fspecial3('gradient',siz2,1);
        	    hy=fspecial3('gradient',siz2,2);
        	    hz=fspecial3('gradient',siz2,3);
                
                Ix_x=convn(I,hx{1},'same');
                Ix_y=convn(Ix_x,hx{2},'same');
		clear Ix_x
                Ix=convn(Ix_y,hx{3},'same');
		clear Ix_y
                Iy_x=convn(I,hy{1},'same');
		clear Iy_x
                Iy_y=convn(Iy_x,hy{2},'same');
		clear Iy_x
                Iy=convn(Iy_y,hy{3},'same');
		clear Iy_y
                Iz_x=convn(I,hz{1},'same');
                Iz_y=convn(Iz_x,hz{2},'same');
		clear Iz_x
                Iz=convn(Iz_y,hz{3},'same');
                
                clear Iz_y hx hy hz
        	    
            	%calculate elements of the structure tensor
                h=fspecial3('gaussianFWHMhalfSize',siz2); %window function
  
            	Ixx1=convn(Ix.^2,h{1},'same');
                Ixx2=convn(Ixx1,h{2},'same');
                Ixx=convn(Ixx2,h{3},'same');
                clear Ixx1 Ixx2
                
                Ixy1=convn(Ix.*Iy,h{1},'same');
                Ixy2=convn(Ixy1,h{2},'same');
                Ixy=convn(Ixy2,h{3},'same');
                clear Ixy1 Ixy2
                
                Ixz1=convn(Ix.*Iz,h{1},'same');
                Ixz2=convn(Ixz1,h{2},'same');
                Ixz=convn(Ixz2,h{3},'same');
                clear Ixz1 Ixz2
                
                Iyy1=convn(Iy.^2,h{1},'same');
                Iyy2=convn(Iyy1,h{2},'same');
                Iyy=convn(Iyy2,h{3},'same');
                clear Iyy1 Iyy2
                
                Iyz1=convn(Iy.*Iz,h{1},'same');
                Iyz2=convn(Iyz1,h{2},'same');
                Iyz=convn(Iyz2,h{3},'same');
                clear Iyz1 Iyz2
                
                Izz1=convn(Iz.^2,h{1},'same');
                Izz2=convn(Izz1,h{2},'same');
                Izz=convn(Izz2,h{3},'same');               
                clear Izz1 Izz2 Iz Ix Iy
                
            	%calculate the eigenvalues of the structure tensor
            	[a,b,c]=size(I);
            	imfeat=cell(3,1);
            	newSize = [1 1 a*b*c];
            	Ieigen = eig3([reshape(Ixx, newSize) reshape(Ixy, newSize) reshape(Ixz, newSize); ...
            	    reshape(Ixy, newSize) reshape(Iyy, newSize) reshape(Iyz, newSize); ...
            	    reshape(Ixz, newSize) reshape(Iyz, newSize) reshape(Izz, newSize)]);
            	Ieigen = sort(Ieigen,1);
            	imfeat{1}=reshape(Ieigen(:,1), [a b c]);
            	imfeat{2}=reshape(Ieigen(:,2), [a b c]);
            	imfeat{3}=reshape(Ieigen(:,3), [a b c]);            
           
        	end
        
        	% normalize gradient
        	maxValue = max(imfeat{3}(:));
        	imfeat{1} = imfeat{1} ./ maxValue;
        	imfeat{2} = imfeat{2} ./ maxValue;
        	imfeat{3} = imfeat{3} ./ maxValue;
   	 otherwise 
   	     error('Unknown feature type!');
	end 

imfeats{u}=imfeat;
end

if strcmp(type,'sortedeigenvaluesstructure') || strcmp(type, 'sortedeigenvalueshessian')
	imfeats = cat(1,imfeats{:});
end

%save([parameter.feature.root type '.mat'], 'imfeats', '-v7.3');
%save([parameter.feature.root input type '.mat'], 'imfeats', '-v7.3');

%profile viewer;
%profile off;
%profsave(profile('info'), ['P:\Filter' input type  '/']);

end

