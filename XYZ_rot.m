function out = XYZ_rot(A, gira)
%% this function change the 3d matrix to make it easy to be used
% A is the matrix, gira should specify the target dimension design, like:
% gira = [1,2,3] means no change, gira = [2,3,1] means XYZ->YZX.
if gira == [1,2,3]
    out = A;
else
    s = size(A);
    out = zeros(s(gira));
    if gira == [2,1,3]
        for k = 1:s(3)
            out(:,:,k) = A(:,:,k)';
        end
    elseif gira == [3,1,2]
        for k = 1:s(2)
            B(:,:,1) = A(:,k,:);
            out(:,:,k) = B';
        end
    elseif gira == [3,2,1]
        for k = 1:s(1)
            B(:,:,1) = A(k,:,:);
            out(:,:,k) = B';
        end
    elseif gira == [1,3,2]
        for k = 1:s(2)
            out(:,:,k) = A(:,k,:);
        end
    elseif gira == [2,3,1]
        for k = 1:s(1)
            out(:,:,k) = A(k,:,:);
        end
    end
end
end