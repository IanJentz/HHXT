function [G,argout] = HXGeom_Airfoil(model,buildType,L,W,H,th_sw,offsets,head_thk,stickout)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

switch buildType
   
    case 'Block'
        L_bk = L;
        W_bk = W;
        W_ch = W_bk-2*th_sw;
%         [L_mean,R_mean] = arcMeanLength(W_ch,th_sw);
        L_mean = heattransLength(W_ch,th_sw);
        L_ch = (L_bk-4*th_sw-2*W_ch) + 2*L_mean;

        
    case 'Channel'
        L_ch = L;
        W_ch = W;
        W_bk = W_ch+2*th_sw;
%         [L_mean,R_mean] = arcMeanLength(W_ch,th_sw);
        L_mean = heattransLength(W_ch,th_sw);
        L_bk = (L_ch-2*L_mean)+4*th_sw+2*W_ch;
    
end

if H~=0 %3D geometry

else %2D
    
    % use decsg to build the geometry out intersecting rectangles and
    % circles
    x1 = -0.5*L_bk;
    x2 = -0.5*L_bk + 2*th_sw + W_ch - offsets(1);
    x3 =  0.5*L_bk - 2*th_sw - W_ch + offsets(1);
    x4 =  0.5*L_bk;
    
    y1 = -0.5*W_bk;
    y2 = -0.5*W_bk + th_sw;
    y3 =  0.5*W_bk - th_sw;
    y4 =  0.5*W_bk;
    
    r1 = W_ch+th_sw - offsets(3);
    r2 = th_sw;
    
    R1 = [3,4,x1,x2,x2,x1,y4,y4,y1,y1]';
    R2 = [3,4,x3,x4,x4,x3,y4,y4,y1,y1]';
    R3 = [3,4,x2,x3,x3,x2,y4,y4,y3,y3]';
    R4 = [3,4,x2,x3,x3,x2,y3,y3,y2,y2]';
    R5 = [3,4,x2,x3,x3,x2,y2,y2,y1,y1]';
    
    centers = [x2,x2,x3,x3;y4,y1,y4,y1] + [ offsets(2)*[1,1,-1,-1] ; offsets(3)*[-1,1,-1,1] ];
    
    x5 = x2+offsets(2)-r1;
    x6 = x3-offsets(2)+r1;
    
    R6 = [3,4,x5,x2,x2,x5,y1,y1,y1+offsets(3),y1+offsets(3)]';
    R7 = [3,4,x6,x3,x3,x6,y1,y1,y1+offsets(3),y1+offsets(3)]';
    R8 = [3,4,x5,x2,x2,x5,y4,y4,y4-offsets(3),y4-offsets(3)]';
    R9 = [3,4,x6,x3,x3,x6,y4,y4,y4-offsets(3),y4-offsets(3)]';
    
    zbuff = zeros(1,6);
    C1 = [1,x2 + offsets(2),y4 - offsets(3),r1,zbuff]';
    C2 = [1,x2 + offsets(2),y1 + offsets(3),r1,zbuff]';
    C3 = [1,x3 - offsets(2),y4 - offsets(3),r1,zbuff]';
    C4 = [1,x3 - offsets(2),y1 + offsets(3),r1,zbuff]';
    C5 = [1,x2,y4,r2,zbuff]';
    C6 = [1,x2,y1,r2,zbuff]';
    C7 = [1,x3,y4,r2,zbuff]';
    C8 = [1,x3,y1,r2,zbuff]';
    
    r3 = 0.5*(x2-th_sw-x1-head_thk);
    r4 = r3+head_thk;
    
    C9  = [1,x1+r4,y1,r3,zbuff]';
    C10 = [1,x1+r4,y1,r4,zbuff]';
    C11 = [1,x4-r4,y1,r3,zbuff]';
    C12 = [1,x4-r4,y1,r4,zbuff]';
    C13 = [1,x1+r4,y4,r3,zbuff]';
    C14 = [1,x1+r4,y4,r4,zbuff]';
    C15 = [1,x4-r4,y4,r3,zbuff]';
    C16 = [1,x4-r4,y4,r4,zbuff]';
    
    R10 = [3,4,x1,x2,x2,x1,y1,y1,y1-stickout,y1-stickout]';
    R11 = [3,4,x3,x4,x4,x3,y1,y1,y1-stickout,y1-stickout]';
    R12 = [3,4,x1,x2,x2,x1,y4,y4,y4+stickout,y4+stickout]';
    R13 = [3,4,x3,x4,x4,x3,y4,y4,y4+stickout,y4+stickout]';
    
    gd = [R1,R2,R3,R4,R5,C1,C2,C3,C4,C5,C6,C7,C8,...
          R6,R7,R8,R9,...
          C9,C10,C11,C12,C13,C14,C15,C16,...
          R10,R11,R12,R13];
    ns = char('R1','R2','R3','R4','R5','C1','C2','C3','C4','C5','C6','C7','C8',...
               'R6','R7','R8','R9',...
               'C9','C10','C11','C12','C13','C14','C15','C16',...
               'R10','R11','R12','R13')';
    sf = '(C1*R1)+(C2*R1)+(C3*R2)+(C4*R2)+(C5*R1)+(C6*R1)+(C7*R2)+(C8*R2)+(R1-C1-C2)+(R2-C3-C4)+R3+R4+R5';
    sf = ['(C1*R1-R8)+(C2*R1-R6)+(C3*R2-R9)+(C4*R2-R7)+(C5*R1)+(C6*R1)+(C7*R2)+(C8*R2)+(R1-C1-C2)+(R2-C3-C4)+R3+R4+R5', ...
          '+(R6-C6)+(R7-C8)+(R8-C5)+(R9-C7)',...
          '+(C10-C9)+(C12-C11)+(C14-C13)+(C16-C15)',...
          '+(C9*R10)+(C11*R11)+(C13*R12)+(C15*R13)']; 

    [dl,bt] = decsg(gd,sf,ns);
    
    % delete superfluous internal edges
    % note boundary edges cannot be deleted
    try % try the overlapping interior lines case first
        bl = [114,113,123,112,122,121,...
              143,131,154,132,155,156,...
              185:190,196:202,235:241,225:230,207:212,218:224,257:263,247:252];

        [dl2,bt2] = csgdel(dl,bt,bl);
        
        bl = [115:117,112:114,132:134,135:137,...
              151:153,148:150,157,161,162,163,164,168,...
              10,92,91,9,36,71,72,...
              59,60,7,8,28,66,67,...
              62,61,15,16,29,98,99,...
              94,93,18,17,37,103,104];
          
        [dl2,bt2] = csgdel(dl2,bt2,bl);
        
        bl = [69,68,44,45,72,73,48,49];
        [dl2,bt2] = csgdel(dl2,bt2,bl);    
        
        bl = [8,51,52,53,...
              59,60,61,9,...
              6,41,40,39,...
              5,31,32,33];
        [dl2,bt2] = csgdel(dl2,bt2,bl);    

    catch % if that didn't work than the lines are not over lapping and the set is different
%         bl = [8,58,41,47,61,5, ...
%               44,40,48,54, ...
%               39,43,55,49, ...
%               6,59,42,56,64,3];
% 
%         [dl2,bt2] = csgdel(dl,bt,bl);
    
    end
    
    % build the geometry
    G = geometryFromEdges(model,dl2);
    argout = centers;


end



    function [L_mean,R_mean] = arcMeanLength(W_ch,th_sw)
    
    a = W_ch; b = th_sw;
    Larc = @(r) r.*acos( (2*a*b+3*b^2+r.^2) ./ ( (2*a+4*b)*r) );
    L_mean = integral(Larc,b,a+b)/a;
    fun = @(R) L_mean - Larc(R);
    R_mean = fzero(fun,(th_sw+0.5*W_ch));
        
    end

    function L = heattransLength(W_ch,th_sw)
        
       H_t = sqrt(th_sw*W_ch+0.75*W_ch^2);
       theta = asin(H_t/(th_sw+W_ch));
       C =sqrt(H_t^2+0.25*W_ch^2);
       A_sliver = 0.5*theta*(th_sw+W_ch)^2 - 0.25*C^2/tan(0.5*theta);
       A = 0.5*W_ch*H_t+2*A_sliver;
       L = A/W_ch;
    end


end

