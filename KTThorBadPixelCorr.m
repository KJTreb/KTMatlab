function img_3d_corr = KTThorBadPixelCorr(img_3d, numviews)
% img_3d = MgThorBadPixelCorr(img_3d)
%input and output projection data (512, 1024, 810)
% correct Thor detector bad pixels

cols = [466; 256; 750; 751; 787; 115; 321; 288; 239; 256; 236; 275; 267; 753; 754; 247; 248; 497; 253; 250; 262; 519; 256; 108; 534; 352; 532; 229; 154; 100; 280; 225; 236]; 
rows = [97; 172; 200; 200; 223; 228; 231; 257; 261; 261; 262; 262; 272; 278; 278; 285; 285; 297; 297; 305; 314; 331; 348; 357; 284; 362; 277; 247; 242; 239; 238; 235; 253];

img_3d_corr=img_3d;

cols=[cols;cols+1;cols-1;cols+2;cols-2];
rows=[rows;rows;rows;rows;rows];
for view = 1:numviews
    for i=1:length(cols)
        row=rows(i);
        col=cols(i);
        img_3d_corr(row, col, view) = (img_3d(row+1, col, view) + img_3d(row+1, col+1, view) + img_3d(row+1, col-1, view) + img_3d(row, col+1, view) + img_3d(row, col-1, view) + img_3d(row-1, col+1, view) + img_3d(row-1, col, view) + img_3d(row-1, col-1, view))./16 + ...
            (img_3d(row+2, col, view) + img_3d(row+2, col+1, view) + img_3d(row+2, col+2, view) + img_3d(row+2, col-1, view) + img_3d(row+2, col-2, view) + img_3d(row+1, col+2, view) + img_3d(row+1, col-2, view) + img_3d(row, col-2, view) + img_3d(row, col+2, view) + img_3d(row-1, col+2, view) + img_3d(row-1, col-2, view) + img_3d(row-2, col, view) + img_3d(row-2, col-1, view)+ img_3d(row-2, col+1, view)+ img_3d(row-2, col-2, view)+ img_3d(row-2, col+2, view))./32;
        %img_3d_corr(row,col,view) = img_3d_corr(row,col,view)/2;
    end
end
end

