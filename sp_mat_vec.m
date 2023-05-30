function b = sp_mat_vec(A,s)

supp = abs(s)>0;
b = A(:,supp)*s(supp);