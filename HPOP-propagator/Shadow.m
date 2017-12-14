%--------------------------------------------------------------------------
%    Computation of fraction (lambda) of solar disk seen by spacecraft
%           Beebe, King, Reasonberg, Preston   June 19971
% 
%                                Vector    Distance
%        Moon wrt Earth          pccor      rpc
%        Earth wrt Sun           ccor       rc
%        Moon wrt Sun            pscor      rps   
%        Satellite wrt Earth     sbcor      rsb  
%        Satellite wrt Sun       bcor       rb 
%        Satellite wrt Moon      sbpcor     rsbp
%
% Last modified:   2015/08/12   M. Mahooti
%
%--------------------------------------------------------------------------
function [lambda,ecltyp] = Shadow(pccor,ccor,pscor,sbcor,bcor,sbpcor)

SAT_Const

% shadow computation - geometric model
% lambda = 1 - no shadow
% lambda = 0 - no sunlight
% 0 < lambda < 1 - partial shadow

% no consideration is given to the change of direction associated
% with partial shadow.
lambda=1;
ecltyp = char(1);

% Check for both eclipses of the Sun by both the Earth and the Moon
% First the Earth
ubcor = zeros(3,1);
rb = norm(bcor);
rc = norm(ccor);

if(rb<=rc)
else
    % get the unit vector of the satellite relative to the sun
    for i=1:3
        ubcor(i)=bcor(i)/rb;
    end
    
    sepp = cross(sbcor,ubcor);
    
    % rsbx is the projection of sbcor along bcor
    rsbx = dot(sbcor,ubcor);
    
    % rs, rp are apparent (from satellite) radii of sun and earth
    % sep is apparent separation of their centers
    rs=R_Sun/rb;
    rp=R_Earth/rsbx;
    sep=sqrt(sepp(1)^2+sepp(2)^2+sepp(3)^2)/rsbx;
    
    lambda = get_lambda(rs,rp,sep);
end
% If no Earth eclipse, check the Moon
if(lambda<1)
    ecltyp = 'E';
    return
else
    
    for i=1:3
        pscor(i) = pccor(i) + ccor(i);
        sbpcor(i) = sbcor(i) - pccor(i);
    end
    
    rps = sqrt(pscor(1)^2+pscor(2)^2+pscor(3)^2);
    
    if(rb<=rps)
        return
    end
    
    %   unit vector of SV wrt Sun
    for i=1:3
        ubcor(i)=bcor(i)/rb;
    end
    
    % rsbx is the projection of sbcor along bcor
    rsbx=dot(sbpcor,ubcor);
    
    % rs, rp are apparent (from satellite) radii of sun and moon
    % sep is apparent angular separation of their centers
    rs=R_Sun/rb;
    rp=R_Moon/rsbx;
    
    sepp=cross(sbpcor,ubcor);
    sep=sqrt(sepp(1)^2+sepp(2)^2+sepp(3)^2)/rsbx;
    
    lambda = get_lambda(rs,rp,sep);
    
    if( lambda<1 )
        ecltyp = 'M';
    end    
end

