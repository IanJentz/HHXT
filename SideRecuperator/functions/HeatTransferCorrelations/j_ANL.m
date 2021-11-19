function j = j_ANL(Re,Pr)  
    
    alpha = 80;
    cosa = 1/cos(pi*alpha/360);
    Re = cosa*Re;
    
    Re_max = 1e5;
    Re_in = Re<=Re_max;
    
    aj = 0.6+0.5*tan(pi*alpha/360);
    ajlam = ( aj*0.1341*(1300^(-0.3319))*(1300/4.1) - 1 ) / (1300+50);
    jdj0 = 1+ajlam*(Re+50);
    jdj0trans = [1+ajlam*(1700+50),1];
    jdj0(Re>=1700) = 1;
    j = interp1([1700,2300],jdj0trans.*[(4.1/1700),(aj*0.1341*(2300^(-0.3319)))],Re,'linear','extrap');
    j(Re<1700) = 4.1./Re(Re<1700);   
    j(Re>2300) = aj*0.1341.*Re(Re>2300).^(-0.3319); 
    j = jdj0.*j;
    
%     a = 0.743878969781;
%     j = a*j;
%     
%     a = 1.58562611186;
%     b = -0.0455530340805;
%     
%     a = 1.64984026007;
%     b = -0.0541261571151;
%     j = j*a.*Re.^b;


    Re_pnts = [3000,4900,6400];

    
    j = 2.9903*Re.^(-0.6648);
    j(Re<=Re_pnts(1)) = 0.0022*Re(Re<=Re_pnts(1)).^(0.1889);
    inb =  Re>Re_pnts(1) & Re<Re_pnts(3);
    w2 = (Re(inb)-Re_pnts(1))/(Re_pnts(3)-Re_pnts(1));
    w1 = 1-w2;
    j(inb) = w1.*0.0022.*Re(inb).^(0.1889) + ...
             w2.*2.9903.*Re(inb).^(-0.6648);
    
         
end