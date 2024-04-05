row_vec = [1 4 5 7 8];
col_vec = row_vec'; % transpose

mat23 = [1 3 4;3 4 6];
mat32 = [1 3; 4 5; 6 7]; % initialize
mul22 = mat23 * mat32; % matrix multiplication

%[m,s] = multiply3(mat23(1,1),mat23(1,2),mat23(1,3))

m = mat23(1,1)*mat23(1,2)*mat23(1,3);
s = mat23(1,1)+mat23(1,2)+mat23(1,3);