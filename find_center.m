function [tic_x,tic_y] = find_center(temp)
imshow(temp);
temp_gray=rgb2gray(temp);%figure;imshow(temp_gray);
[row,col] = find(temp_gray > 0);
tic_x = sum(col)/length(col)-1;
tic_y = sum(row)/length(row)-1;
end

