path = 'test_image';
num = 3;

for i = 1:1:num  
	file = [path,'/image',int2str(i),'.jpg']; 
	img=imread(file);

% simple linear iterativeclustering;
	[l, Am, C] = slic(img,200, 10, 1, 'median');

%show the superpixels results
	figure(i);
	I1=drawregionboundaries(l, img, [0 0 0]);
	imshow(I1);
	hold on
    imwrite(I1,[path,'/superpixels/super_',strcat(int2str(i),'.jpg')]); 

end
