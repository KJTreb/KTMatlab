function img_3dnew = KTThorBadPixelCorr(img_3d, frames)
% img_3d = MgThorBadPixelCorr(img_3d)
%input and output projection data (512, 1024, frames)
% correct Thor detector bad pixels
% up-to-date as of 08-23-2021
img_3dnew = img_3d;

cols = [466; 256; 750; 751; 787; 115; 321; 288; 239; 256; 236; 275; 267; 753; 754; 247; 248; 497; 253; 250; 262; 519; 256; 108; 534; 352; 532; 229; 154; 100; 280; 225; 236]; 
rows = [97; 172; 200; 200; 223; 228; 231; 257; 261; 261; 262; 262; 272; 278; 278; 285; 285; 297; 297; 305; 314; 331; 348; 357; 284; 362; 277; 247; 242; 239; 238; 235; 253];

% more pixels in HE bin
colsnew = [157; 147; 151; 166; 176; 204; 204; 236; 236; 250; 255; 132; 146; 187; 199; 251; 249; 246; 237; 370; 376; 567; 39; 24; 24; 877; 767; 971; 351; 352; 351; 352; 582; 876; 877; 878; 876; 877; 878; 610; 611; 315]+1; %x
rowsnew = [225; 245; 247; 249; 247; 244; 250; 244; 256; 234; 238; 217; 238; 211; 244; 217; 198; 194; 195; 152; 269; 392; 403; 449; 448; 86; 254; 256; 360; 360; 362; 362; 371; 85; 85; 85; 87; 87; 87; 296; 296; 463]+1; %y

cols=[cols; colsnew];
rows=[rows; rowsnew];

cols=[cols;cols+1;cols-1];
rows=[rows;rows;rows];
for view = 1:frames
    for i=1:length(cols)
        row=rows(i);
        col=cols(i);
%         img_3dnew(row, col, view) = (img_3d(row+1, col, view) + img_3d(row+1, col+1, view) + img_3d(row+1, col-1, view) + img_3d(row, col+1, view) + img_3d(row, col-1, view) + img_3d(row-1, col+1, view) + img_3d(row-1, col, view) + img_3d(row-1, col-1, view))./8 + ...
%             (img_3d(row+2, col, view) + img_3d(row+2, col+1, view) + img_3d(row+2, col+2, view) + img_3d(row+2, col-1, view) + img_3d(row+2, col-2, view) + img_3d(row+1, col+2, view) + img_3d(row+1, col-2, view) + img_3d(row, col-2, view) + img_3d(row, col+2, view) + img_3d(row-1, col+2, view) + img_3d(row-1, col-2, view) + img_3d(row-2, col, view) + img_3d(row-2, col-1, view)+ img_3d(row-2, col+1, view)+ img_3d(row-2, col-2, view)+ img_3d(row-2, col+2, view))./16;
        
          % nearest neighbor  
%         img_3dnew(row, col, view) = img_3d(row+7, col+7, view);
        img_3dnew(row, col, view) = 1/4.*img_3d(row+4, col, view)+1/4.*img_3d(row, col+4, view)+ 1/4.*img_3d(row-4, col, view)+1/4.*img_3d(row, col-4, view);
    end
end
end

