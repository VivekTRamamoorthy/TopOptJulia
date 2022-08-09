#  AN TOPOLOGY OPTIMIZATION CODE IN JULIA %%%%
using SparseArrays
using SuiteSparse
using LinearAlgebra

# % function top88(nelx,nely,volfrac,penal,rmin,ft)

# % top88(100,100,0.5,1,1,2)
# clc;
# nelx=20;nely=20;volfrac=0.5;penal=2;rmin=1;ft=2;
nelx=5;nely=5;volfrac=0.5;penal=2;rmin=1;ft=2;
# %% MATERIAL PROPERTIES
# E0 = 1;
# Emin = 1e-9;
# nu = 0.3;
E0 = 1;
Emin = 1e-9;
nu = 0.3;
# %% PREPARE FINITE ELEMENT ANALYSIS
# A11 = [12 3 -6 -3; 3 12 3 0; -6 3 12 -3; -3 0 -3 12];
# A12 = [-6 -3 0 3; -3 -6 -3 -6; 0 -3 -6 3; 3 -6 3 -6];
# B11 = [-4 3 -2 9; 3 -4 -9 4; -2 -9 -4 -3; 9 4 -3 -4];
# B12 = [ 2 -3 4 -9; -3 2 9 -2; 4 9 2 3; -9 -2 3 2];

# %% PREPARE FINITE ELEMENT ANALYSIS
A11 = [12.0 3 -6 -3; 3 12 3 0; -6 3 12 -3; -3 0 -3 12];
A12 = [-6.0 -3 0 3; -3 -6 -3 -6; 0 -3 -6 3; 3 -6 3 -6];
B11 = [-4.0 3 -2 9; 3 -4 -9 4; -2 -9 -4 -3; 9 4 -3 -4];
B12 = [ 2.0 -3 4 -9; -3 2 9 -2; 4 9 2 3; -9 -2 3 2];

# KE = 1/(1-nu^2)/24*([A11 A12;A12' A11]+nu*[B11 B12;B12' B11]);
KE = 1/(1-nu^2)/24*([A11 A12;A12' A11]+nu*[B11 B12;B12' B11]);

# nodenrs = reshape(1:(1+nelx)*(1+nely),1+nely,1+nelx);
# edofVec = reshape(2*nodenrs(1:end-1,1:end-1)+1,nelx*nely,1);

nodenrs = reshape(1.0:1:(1+nelx)*(1+nely),1+nely,1+nelx);
edofVec = reshape(2*nodenrs[1:end-1,1:end-1].+1,nelx*nely,1);


function repmat(matrix,nrows::Int,ncols::Int)
    rescol=matrix;
    for i=1:nrows-1
        rescol = [rescol; matrix];
    end
    res = rescol;
    for j=1:ncols-1
        res=[res rescol];
    end
    return res
end
            

edofMat = Int.( repmat(edofVec,1,8)+repmat([0 1 2*nely.+[2 3 0 1] -2 -1],nelx*nely,1));


# iK = reshape(kron(edofMat,ones(8,1))',64*nelx*nely,1);
# jK = reshape(kron(edofMat,ones(1,8))',64*nelx*nely,1);

iK = vec(kron(edofMat,ones(8,1))');
jK = vec(kron(edofMat,ones(1,8))');


# % DEFINE LOADS AND SUPPORTS (HALF MBB-BEAM)
# F = sparse(2,1,-1,2*(nely+1)*(nelx+1),1);
# U = zeros(2*(nely+1)*(nelx+1),1);

F = sparse([2],[1],[1.0],2*(nely+1)*(nelx+1),1);
U = zeros(2*(nely+1)*(nelx+1),1);


# fixeddofs = union([1:2:2*(nely+1)],[2*(nelx+1)*(nely+1)]);
fixeddofs = union(1:2:2*(nely+1),[2*(nelx+1)*(nely+1)]);

# alldofs = [1:2*(nely+1)*(nelx+1)];
alldofs = 1:2*(nely+1)*(nelx+1);


# freedofs = setdiff(alldofs,fixeddofs);
freedofs = setdiff(alldofs,fixeddofs);


# %% PREPARE FILTER
# iH = ones(nelx*nely*(2*(ceil(rmin)-1)+1)^2,1);
# jH = ones(size(iH));
# sH = zeros(size(iH));
iH = ones(nelx*nely*(2*(ceil(rmin)-1)+1)^2);
jH = ones(size(iH));
sH = zeros(size(iH));


# k = 0;

# for i1 = 1:nelx
#     for j1 = 1:nely
#         e1 = (i1-1)*nely+j1;
#         for i2 = max(i1-(ceil(rmin)-1),1):min(i1+(ceil(rmin)-1),nelx)
#             for j2 = max(j1-(ceil(rmin)-1),1):min(j1+(ceil(rmin)-1),nely)
#                 e2 = (i2-1)*nely+j2;
#                 k = k+1;
#                 iH(k) = e1;
#                 jH(k) = e2;
#                 sH(k) = max(0,rmin-sqrt((i1-i2)^2+(j1-j2)^2));
#             end
#         end
#     end
# end
ik=0;

for i1 = 1:nelx
    for j1 = 1:nely
        e1 = (i1-1)*nely+j1;
        for i2 = max(i1-(ceil(rmin)-1),1):min(i1+(ceil(rmin)-1),nelx)
            for j2 = max(j1-(ceil(rmin)-1),1):min(j1+(ceil(rmin)-1),nely)
                e2 = (i2-1)*nely+j2;
                 global ik = ik+1;
                 iH[ik] = e1;
                 jH[ik] = e2;
                 sH[ik] = max(0,rmin-sqrt((i1-i2)^2+(j1-j2)^2));
            end
        end
    end
end




# H = sparse(iH,jH,sH);
H = sparse(iH,jH,sH);


# Hs = sum(H,2);
Hs = sum(H,dims = 2);


# %% INITIALIZE ITERATION
# x = repmat(volfrac,nely,nelx);
x = repmat(volfrac,nely,nelx);


# xPhys = x;
xPhys = x;

# loop = 0;
loop = 0;

function colon(matrixA)
    colonA=zeros(size(matrixA,1)*size(matrixA,2));
    index=1
    for col = 1:size(matrixA,1)
        for row = 1:size(matrixA,2)
            colonA[index] = matrixA[row,col];
            index = index+1;
        end
    end
    return colonA;
end

# change = 1;
change = 1;

# %% START ITERATION
# for loop=1:2
#     % FE-ANALYSIS
#     sK = reshape(KE(:)*(Emin+xPhys(:)'.^penal*(E0-Emin)),64*nelx*nely,1);
    sK = vec(KE[:].*(Emin.+colon(xPhys)'.^penal*(E0-Emin))) #64*nelx*nely,1);
    print(size(sK))
# 
#     K = sparse(iK,jK,sK); K = (K+K')/2;

    K = sparse(iK,jK,sK); 
    K = (K+K')/2;
    Kfree=K[freedofs,freedofs]
    Ffree=F[freedofs]

    Ufree=sparsevec([],[],length(Ffree));
#     U(freedofs) = K(freedofs,freedofs)\F(freedofs);
    Kfree=Matrix(Kfree);
    # L=factorize(Symmetric(Kfree))
    # CHOLMOD.Factor(Kfree)
    Ufree = Kfree\Ffree;

    # Ffree=F(freedofs);
#     Kfree=full(K(freedofs,freedofs));
#     Ufree=U(freedofs);
    U=zeros(length(F));
    U[freedofs] = Ufree;

#     % OBJECTIVE FUNCTION AND SENSITIVITY ANALYSIS
#     ce = reshape(sum((U(edofMat)*KE).*U(edofMat),2),nely,nelx);
    ce = reshape(  sum((U[edofMat]*KE).*U[edofMat],dims = 2),  (nely,nelx));

# %     display(ce)
#     c = sum(sum((Emin+xPhys.^penal*(E0-Emin)).*ce));
    c = sum((Emin.+xPhys.^penal.*(E0-Emin)).*ce);



#     disp([' compliance c = ' num2str(c)])
    print("compliance c =  $c ")
#     dc = -penal*(E0-Emin)*xPhys.^(penal-1).*ce;
    dc = -penal*(E0-Emin).*xPhys.^(penal-1).*ce;

# %     display(xPhys)
    print(xPhys)

#     dv = ones(nely,nelx);
    dv = ones(nely,nelx);

#     % FILTERING/MODIFICATION OF SENSITIVITIES
#     if ft == 1
#         dc(:) = H*(x(:).*dc(:))./Hs./max(1e-3,x(:));
#     elseif ft == 2
#         dc(:) = H*(dc(:)./Hs);
#         dv(:) = H*(dv(:)./Hs);
#     end
    if ft == 1
        # dc(:) = H*(x(:).*dc(:))./Hs./max(1e-3,x(:));
         dc[:] = H*(x[:].*dc[:])./Hs./max(1e-3*ones(size(x[:])),x[:]);

    elseif ft == 2
        # dc(:) = H*(dc(:)./Hs);
        # dv(:) = H*(dv(:)./Hs);
        dc[:] = H*(dc[:]./Hs);
        dv[:] = H*(dv[:]./Hs);
    end



#     % OPTIMALITY CRITERIA UPDATE OF DESIGN VARIABLES AND PHYSICAL DENSITIES

#     l1 = 0; l2 = 1e9; move = 0.2;
    l1 = 0; l2 = 1e9; move = 0.2;
#     whileloopcounter=0;
    whileloopcounter=0;

#     while (l2-l1)/(l1+l2) > 1e-3
#         lmid = 0.5*(l2+l1);
#         xnew = max(0,max(x-move,min(1,min(x+move,x.*sqrt(-dc./dv/lmid)))));
# %                    display(xnew)
#         if ft == 1
#             xPhys = xnew;
#         elseif ft == 2
#             xPhys(:) = (H*xnew(:))./Hs;
#  
# %             disp(['xPhy(1,1) = ' num2str(xPhys(1,1))])
#         end
#         if sum(xPhys(:)) > volfrac*nelx*nely, l1 = lmid; else l2 = lmid; end
#         disp(['while loop counter: ' num2str(whileloopcounter) ' (l2-l1)/(l1+l2): ' num2str((l2-l1)/(l1+l2))])
# %         disp(['x(1,1) ' num2str(x(1,1))])
#         whileloopcounter=whileloopcounter+1;
#         if whileloopcounter>100,break,end
# 
#     end
    while (l2-l1)/(l1+l2) > 1e-3
        global l2,l1
        lmid = 0.5*(l2+l1);
        xnew = max(0,max(x-move,min(1,min(x+move,x.*sqrt(-dc./dv/lmid)))));
# %                    display(xnew)
        if ft == 1
            xPhys = xnew;
        elseif ft == 2
            xPhys[:] = (H*xnew[:])./Hs;
 
# %             disp(['xPhy(1,1) = ' num2str(xPhys(1,1))])
        end
        if sum(xPhys(:)) > volfrac*nelx*nely; l1 = lmid; 
        else l2 = lmid; 
        end
        print("while loop counter:  $(whileloopcounter)  (l2-l1)/(l1+l2):  $((l2-l1)/(l1+l2))")
# %         disp(['x(1,1) ' num2str(x(1,1))])
        whileloopcounter=whileloopcounter+1;
        if whileloopcounter>100
            break
        end

    end
#     change = max(abs(xnew(:)-x(:)));
#     disp([ 'change = ' num2str(change)])
#     x = xnew;
#         display(xPhys)
#     %% PRINT RESULTS
# %     fprintf(' It.:%5i Obj.:%11.4f Vol.:%7.3f ch.:%7.3f\n',loop,c, mean(xPhys(:)),change);
#     %% PLOT DENSITIES
#     colormap(gray); imagesc(1-xPhys); caxis([0 1]); axis equal; axis off; drawnow;
# end







# %%%% AN 88 LINE TOPOLOGY OPTIMIZATION CODE %%%%
# % function top88(nelx,nely,volfrac,penal,rmin,ft)
# % top88(100,100,0.5,1,1,2)
# clc;
# nelx=20;nely=20;volfrac=0.5;penal=2;rmin=1;ft=2;
# %% MATERIAL PROPERTIES
# E0 = 1;
# Emin = 1e-9;
# nu = 0.3;
# %% PREPARE FINITE ELEMENT ANALYSIS
# A11 = [12 3 -6 -3; 3 12 3 0; -6 3 12 -3; -3 0 -3 12];
# A12 = [-6 -3 0 3; -3 -6 -3 -6; 0 -3 -6 3; 3 -6 3 -6];
# B11 = [-4 3 -2 9; 3 -4 -9 4; -2 -9 -4 -3; 9 4 -3 -4];
# B12 = [ 2 -3 4 -9; -3 2 9 -2; 4 9 2 3; -9 -2 3 2];
# KE = 1/(1-nu^2)/24*([A11 A12;A12' A11]+nu*[B11 B12;B12' B11]);
# nodenrs = reshape(1:(1+nelx)*(1+nely),1+nely,1+nelx);
# edofVec = reshape(2*nodenrs(1:end-1,1:end-1)+1,nelx*nely,1);
# edofMat = repmat(edofVec,1,8)+repmat([0 1 2*nely+[2 3 0 1] -2 -1],nelx*nely,1);
# iK = reshape(kron(edofMat,ones(8,1))',64*nelx*nely,1);
# jK = reshape(kron(edofMat,ones(1,8))',64*nelx*nely,1);
# % DEFINE LOADS AND SUPPORTS (HALF MBB-BEAM)
# F = sparse(2,1,-1,2*(nely+1)*(nelx+1),1);
# U = zeros(2*(nely+1)*(nelx+1),1);
# fixeddofs = union([1:2:2*(nely+1)],[2*(nelx+1)*(nely+1)]);
# alldofs = [1:2*(nely+1)*(nelx+1)];
# freedofs = setdiff(alldofs,fixeddofs);
# %% PREPARE FILTER
# iH = ones(nelx*nely*(2*(ceil(rmin)-1)+1)^2,1);
# jH = ones(size(iH));
# sH = zeros(size(iH));
# k = 0;
# for i1 = 1:nelx
#     for j1 = 1:nely
#         e1 = (i1-1)*nely+j1;
#         for i2 = max(i1-(ceil(rmin)-1),1):min(i1+(ceil(rmin)-1),nelx)
#             for j2 = max(j1-(ceil(rmin)-1),1):min(j1+(ceil(rmin)-1),nely)
#                 e2 = (i2-1)*nely+j2;
#                 k = k+1;
#                 iH(k) = e1;
#                 jH(k) = e2;
#                 sH(k) = max(0,rmin-sqrt((i1-i2)^2+(j1-j2)^2));
#             end
#         end
#     end
# end
# H = sparse(iH,jH,sH);
# Hs = sum(H,2);
# %% INITIALIZE ITERATION
# x = repmat(volfrac,nely,nelx);
# xPhys = x;
# loop = 0;
# change = 1;
# %% START ITERATION
# for loop=1:36
#     % FE-ANALYSIS
#     sK = reshape(KE(:)*(Emin+xPhys(:)'.^penal*(E0-Emin)),64*nelx*nely,1);
# 
#     K = sparse(iK,jK,sK); K = (K+K')/2;
#     U(freedofs) = K(freedofs,freedofs)\F(freedofs);
#     
#     Ffree=F(freedofs);
#     Kfree=full(K(freedofs,freedofs));
#     Ufree=U(freedofs);
#     % OBJECTIVE FUNCTION AND SENSITIVITY ANALYSIS
#     ce = reshape(sum((U(edofMat)*KE).*U(edofMat),2),nely,nelx);
# %     display(ce)
#     c = sum(sum((Emin+xPhys.^penal*(E0-Emin)).*ce));
#     disp([' compliance c = ' num2str(c)])
#     dc = -penal*(E0-Emin)*xPhys.^(penal-1).*ce;
# %     display(xPhys)
#     dv = ones(nely,nelx);
#     % FILTERING/MODIFICATION OF SENSITIVITIES
#     if ft == 1
#         dc(:) = H*(x(:).*dc(:))./Hs./max(1e-3,x(:));
#     elseif ft == 2
#         dc(:) = H*(dc(:)./Hs);
#         dv(:) = H*(dv(:)./Hs);
#     end
#     % OPTIMALITY CRITERIA UPDATE OF DESIGN VARIABLES AND PHYSICAL DENSITIES
#     l1 = 0; l2 = 1e9; move = 0.2;
#     whileloopcounter=0;
#     while (l2-l1)/(l1+l2) > 1e-3
#         lmid = 0.5*(l2+l1);
#         xnew = max(0,max(x-move,min(1,min(x+move,x.*sqrt(-dc./dv/lmid)))));
# %                    display(xnew)
#         if ft == 1
#             xPhys = xnew;
#         elseif ft == 2
#             xPhys(:) = (H*xnew(:))./Hs;
#  
# %             disp(['xPhy(1,1) = ' num2str(xPhys(1,1))])
#         end
#         if sum(xPhys(:)) > volfrac*nelx*nely, l1 = lmid; else l2 = lmid; end
#         disp(['while loop counter: ' num2str(whileloopcounter) ' (l2-l1)/(l1+l2): ' num2str((l2-l1)/(l1+l2))])
# %         disp(['x(1,1) ' num2str(x(1,1))])
#         whileloopcounter=whileloopcounter+1;
#         if whileloopcounter>100,break,end
# 
#     end
#     change = max(abs(xnew(:)-x(:)));
#     disp([ 'change = ' num2str(change)])
#     x = xnew;
#         display(xPhys)
#     %% PRINT RESULTS
# %     fprintf(' It.:%5i Obj.:%11.4f Vol.:%7.3f ch.:%7.3f\n',loop,c, mean(xPhys(:)),change);
#     %% PLOT DENSITIES
#     colormap(gray); imagesc(1-xPhys); caxis([0 1]); axis equal; axis off; drawnow;
# end
# 
# 
# 
# 
# % end
# 