function lambda = get_lambda(rs,rp,sep)

% Calculate lambda
% 
% rs : apparent radius of sun as viewed from satellite (radians)
% rp : apparent radius of eclipsing body as viewed from satellite (radians)
% sep: apparent separation of the center of the Sun and eclipsing body (radians)
% 
% lambda : fraction of Sun's disk visible (1.0 = no eclipse; 0 = total eclipse)

if (rs+rp<=sep)
    % no eclipse
    lambda = 1;
    return
elseif( rp-rs>=sep )
    % full eclipse
    lambda = 0;
    return
else
    % partial eclipse, do the calculations
    if(sep<=rs-rp)
    else
        % set r1 = smaller disc, r2 = larger
        if(rs>rp)
            r1=rp;
            r2=rs;
        else
            r1=rs;
            r2=rp;        
        end        
        % phi = 1/2 angle subtended in disc 1 by arc of intersection
        phi = acos((r1*r1+sep*sep-r2*r2)/(2*r1*sep));
        if (phi<0)
            phi = pi + phi;
        end
        if(r2/r1>5)
            hgt=sqrt(r1^2-(sep-r2)^2);
            area2=hgt*(sep-r2);
            area3=0;
        else
            % thet = 1/2 angle subtended in disc 2 by arc of intersection
            % hgt  = 1/2 linear distance between ends of arc of intersection
            hgt=r1*sin(phi);
            thet=asin(hgt/r2);
            area2=sep*hgt;
            area3=thet*r2^2;            
        end
        % one disc much bigger - treat boundary as a straight line        
        area1=(pi-phi)*r1^2;
        % ari = area of non-overlapped portion of small disc
        ari=area1+area2-area3;
        area1=pi*rs^2;
        if(rs>rp)
            area2=pi*rp^2;
            lambda=(area1+ari-area2)/area1;
            return
        else
            % sun is small disc
            lambda=ari/area1;
            return
        end
        % eclipsing body is small disc
    end
    % eclipsing body lies within sun's disc - what fraction of sun's disk is blocked
    lambda=(rs^2-rp^2)/rs^2;    
end

