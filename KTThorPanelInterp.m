function img_interp = KTThorPanelInterp(img)
% up-to-date as of 08-23-2021

img_interp = img;

vert_right_top = [129, 257, 385, 513, 641, 769, 897];
vert_left_top = [128, 256, 384, 512, 640, 768, 896];

vert_right_bottom = [109, 237, 365, 493, 621, 749, 877, 1005];
vert_left_bottom = [108, 236, 364, 492, 620, 748, 876, 1004];

% % nearest neighbor interp
% img_interp(256,:) = img(255,:);
% img_interp(257,:) = img(258,:);
% 
% img_interp(1:256,vert_right_top) = img(1:256,vert_right_top+1);
% img_interp(1:256,vert_left_top) = img(1:256,vert_left_top-1);
% 
% img_interp(257:512,vert_right_bottom) = img(257:512,vert_right_bottom+1);
% img_interp(257:512,vert_left_bottom) = img(257:512,vert_left_bottom-1);

% linear interp
img_interp(256,:) = 2/3.*img(255,:)+1/3.*img(258,:);
img_interp(257,:) = 2/3.*img(258,:)+1/3.*img(255,:);

img_interp(1:256,vert_right_top) = 2/3.*img(1:256,vert_right_top+1)+1/3.*img(1:256,vert_right_top-2);
img_interp(1:256,vert_left_top) = 2/3.*img(1:256,vert_left_top-1)+1/3.*img(1:256,vert_right_top+2);

img_interp(257:512,vert_right_bottom) = 2/3.*img(257:512,vert_right_bottom+1)+1/3.*img(257:512,vert_right_bottom-2);
img_interp(257:512,vert_left_bottom) = 2/3.*img(257:512,vert_left_bottom-1)+1/3.*img(257:512,vert_right_bottom+2);
        
end