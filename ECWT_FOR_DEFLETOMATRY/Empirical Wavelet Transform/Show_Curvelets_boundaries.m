function Show_Curvelets_boundaries(f,Bw,Bt,option)

%=======================================================================
% function Show_Curvelets_boundaries(f,Bw,Bt,option)
% 
% This function plots the angular sectors corresponding to the 
% detected scales and angles in the 2D spectrum.
%
% Inputs:
%   -f: input image
%   -Bw: list of Fourier boundaries (in [0,pi])
%   -Bt: list of angles (in radian)
%   -option: 1=EWTC-I;2=EWTC-II; 3=EWTC-III (must be the same as the one 
%            used in the transform)
%
% Author: Jerome Gilles
% Institution: SDSU - Department of Mathematics and Statistics
% Version: 1.0 - 2013
% Version: 2.0 - 2015
%=======================================================================
ff=log(1+fftshift(abs(fft2(f))));
figure;imshow(ff,[]);
hold on;

color='white';
wid = 1;

po=floor(size(f)/2);
[r,c]=size(f);

if option==1  %Option 1
    % We plot the scale annuli
    for n=1:length(Bw)
        a=Bw(n)*ceil(size(f,2)/(2*pi));
        b=Bw(n)*ceil(size(f,1)/(2*pi));
        EWT_drawEllipse(ceil(size(f,2)/2)+1,ceil(size(f,1)/2)+1,a,b,color);
    end
    % We plot the angles limits
    for n=1:length(Bt)
        if abs(Bt(n))<=pi/4
            hold on
            p0 = po.*(1 + Bw(1)*[cos(Bt(n)) sin(Bt(n))]/pi);
            p1 = [c ceil((r+c*tan(Bt(n)))/2)+1];
            plot([p1(1),p0(1)],[p1(2),p0(2)],'Color',color,'LineWidth',wid);   
            p0 = po.*(1 - Bw(1)*[cos(Bt(n)) sin(Bt(n))]/pi);
            p2 = [1 ceil((r-c*tan(Bt(n)))/2)+1];
            plot([p2(1),p0(1)],[p2(2),p0(2)],'Color',color,'LineWidth',wid);           
        else
            hold on
            p0 = po.*(1 - Bw(1)*[cos(Bt(n)) sin(Bt(n))]/pi);
            p1 = [ceil((c+r*cot(Bt(n)))/2)+1 r];
            plot([p1(1),p0(1)],[p1(2),p0(2)],'Color',color,'LineWidth',wid);   
            p0 = po.*(1 + Bw(1)*[cos(Bt(n)) sin(Bt(n))]/pi);
            p2 = [ceil((c-r*cot(Bt(n)))/2)+1 1];
            plot([p2(1),p0(1)],[p2(2),p0(2)],'Color',color,'LineWidth',wid);   
        end
    end
elseif option==2
    % We plot the scale annuli
    for n=1:length(Bw)
        a=Bw(n)*floor(size(f,2)/(2*pi));
        b=Bw(n)*floor(size(f,1)/(2*pi));
        EWT_drawEllipse(floor(size(f,2)/2)+1,floor(size(f,1)/2)+1,a,b,color);
    end
    % We plot the angles limits
    for s=1:length(Bw)-1
        for n=1:length(Bt{s})
            if abs(Bt{s}(n))<=pi/4
                hold on
                p0 = po + Bw(s)*floor(size(f,2)/(2*pi))*[cos(Bt{s}(n)) sin(Bt{s}(n))];
                p1 = po + Bw(s+1)*floor(size(f,2)/(2*pi))*[cos(Bt{s}(n)) sin(Bt{s}(n))];
                plot([p1(1),p0(1)],[p1(2),p0(2)],'Color',color,'LineWidth',wid);   
                p0 = po - Bw(s)*floor(size(f,2)/(2*pi))*[cos(Bt{s}(n)) sin(Bt{s}(n))]+1;
                p2 = po - Bw(s+1)*floor(size(f,2)/(2*pi))*[cos(Bt{s}(n)) sin(Bt{s}(n))]+1;
                plot([p2(1),p0(1)],[p2(2),p0(2)],'Color',color,'LineWidth',wid);           
            else
                hold on
                p0 = po - Bw(s)*floor(size(f,2)/(2*pi))*[cos(Bt{s}(n)) sin(Bt{s}(n))]+1;
                p1 = po - Bw(s+1)*floor(size(f,2)/(2*pi))*[cos(Bt{s}(n)) sin(Bt{s}(n))]+1;
                plot([p1(1),p0(1)],[p1(2),p0(2)],'Color',color,'LineWidth',wid);   
                p0 = po + Bw(s)*floor(size(f,2)/(2*pi))*[cos(Bt{s}(n)) sin(Bt{s}(n))]+1;
                p2 = po + Bw(s+1)*floor(size(f,2)/(2*pi))*[cos(Bt{s}(n)) sin(Bt{s}(n))]+1;
                plot([p2(1),p0(1)],[p2(2),p0(2)],'Color',color,'LineWidth',wid);   
            end
        end
    end
    for n=1:length(Bt{end})
        if abs(Bt{end}(n))<=pi/4
            hold on
            p0 = po + Bw(end)*floor(size(f,2)/(2*pi))*[cos(Bt{end}(n)) sin(Bt{end}(n))];
            p1 = [c floor((r+c*tan(Bt{end}(n)))/2)+1];
            plot([p1(1),p0(1)],[p1(2),p0(2)],'Color',color,'LineWidth',wid);   
            p0 = po - Bw(end)*floor(size(f,2)/(2*pi))*[cos(Bt{end}(n)) sin(Bt{end}(n))];
            p2 = [1 floor((r-c*tan(Bt{end}(n)))/2)+1];
            plot([p2(1),p0(1)],[p2(2),p0(2)],'Color',color,'LineWidth',wid);           
        else
            hold on
            p0 = po - Bw(end)*floor(size(f,2)/(2*pi))*[cos(Bt{end}(n)) sin(Bt{end}(n))];
            p1 = [floor((c+r*cot(Bt{end}(n)))/2)+1 r];
            plot([p1(1),p0(1)],[p1(2),p0(2)],'Color',color,'LineWidth',wid);   
            p0 = po + Bw(end)*floor(size(f,2)/(2*pi))*[cos(Bt{end}(n)) sin(Bt{end}(n))];
            p2 = [floor((c-r*cot(Bt{end}(n)))/2)+1 1];
            plot([p2(1),p0(1)],[p2(2),p0(2)],'Color',color,'LineWidth',wid);   
        end
    end
elseif option==3
    a=Bw{1}*ceil(size(f,2)/(2*pi));
    b=Bw{1}*ceil(size(f,1)/(2*pi));
    EWT_drawEllipse(ceil(size(f,2)/2)+1,ceil(size(f,1)/2)+1,a,b,color);
    
    % We plot the angles limits
    for n=1:length(Bt)
        if abs(Bt(n))<=pi/4
            hold on
            p0 = po.*(1 + Bw{1}*[cos(Bt(n)) sin(Bt(n))]/pi)+1;
            p1 = [c floor((r+c*tan(Bt(n)))/2)+1];
            plot([p1(1),p0(1)],[p1(2),p0(2)],'Color',color,'LineWidth',wid);   
            p0 = po.*(1 - Bw{1}*[cos(Bt(n)) sin(Bt(n))]/pi)+1;
            p2 = [1 ceil((r-c*tan(Bt(n)))/2)+1];
            plot([p2(1),p0(1)],[p2(2),p0(2)],'Color',color,'LineWidth',wid);           
        else
            hold on
            p0 = po.*(1 - Bw{1}*[cos(Bt(n)) sin(Bt(n))]/pi)+1;
            p1 = [ceil((c+r*cot(Bt(n)))/2)+1 r];
            plot([p1(1),p0(1)],[p1(2),p0(2)],'Color',color,'LineWidth',wid);   
            p0 = po.*(1 + Bw{1}*[cos(Bt(n)) sin(Bt(n))]/pi)+1;
            p2 = [ceil((c-r*cot(Bt(n)))/2)+1 1];
            plot([p2(1),p0(1)],[p2(2),p0(2)],'Color',color,'LineWidth',wid);   
        end
        end
    
    % We plot the scale limits per angular sector
    for t=1:length(Bt)-1
        for s=1:length(Bw{t+1})
            a=Bw{t+1}(s)*ceil(size(f,2)/(2*pi));
            b=Bw{t+1}(s)*ceil(size(f,1)/(2*pi));
            EWT_drawArcEllipse(ceil(size(f,2)/2)+1,ceil(size(f,1)/2)+1,a,b,Bt(t),Bt(t+1),color);
        end
    end
    
    %last angular sector
    for s=1:length(Bw{end})
        a=Bw{end}(s)*ceil(size(f,2)/(2*pi));
        b=Bw{end}(s)*ceil(size(f,1)/(2*pi));
        EWT_drawArcEllipse(ceil(size(f,2)/2)+1,ceil(size(f,1)/2)+1,a,b,Bt(end),Bt(1)+pi,color);
    end
end
hold off