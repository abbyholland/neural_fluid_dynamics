function [vortex, info] = IDvortex(x,y,Vx,Vy)

% IDvortex location finder
% Detects a vortex location and size in Vx and Vy velocity fields
%   Calls the two functions IDvortexLoc and IDvortexSize
%   Then combines the results
%   supply:
%       x,y coordinate vectors (mm) and Vx and Vy fields (m/s)
%       coordinate vectors in mm (otherwise displayed units are wrong)
%   returns:
%       vortex position in true coordinates (unit)
%       size (unit^2)
%       circulation

% Based on the paper
% Combining PIV, POD and vortex identification algorithms for the study of
% unsteady turbulent swirling flows
% Laurent Graftieaux, Marc Michard and Nathalie Grosjean
% Ecole Centrale de Lyon

% Version #: 1.3, 2015-11-25
% sebastian@email.arizona.edu

MatLabSettings
if SETTINGS.IDvortex.TroubleShooting==true
    SETTINGS.IDvortex.mute=false;
end
%% Calculations

[xLoc, yLoc]=IDvortexLoc(x,y,Vx,Vy);
[area, circulation, points]=IDvortexSize(x,y,Vx,Vy);

disp(['IDvortex: found ' num2str(length(xLoc)) ' G1 peaks.' ])
disp(['IDvortex: StDev of ' num2str(length(area)) ' found G2 peaks = ' num2str(std(area))])
info.Loc.amount = length(xLoc);
info.Size.amount = length(area);
info.Size.StdDev = std(area);

for ii=1:length(xLoc)
    for kk=1:length(area)
        if length(points{kk}(:,1))<4
            in(kk)=false;
        else
            edges = convhull(points{kk}(:,1),points{kk}(:,2));
            in(kk) = inpolygon(xLoc(ii),yLoc(ii),points{kk}(edges,1),points{kk}(edges,2)); % true if VortexLoc is inside area
        end
        if in(kk)
            vortex(ii,:) = [xLoc(ii), yLoc(ii), area(kk), circulation(kk)];
			if SETTINGS.IDvortex.mute==false
				disp(['IDvortex: Found one vortex at (' num2str(xLoc(ii)) ';' num2str(yLoc(ii)) ') with area=' ...
                    num2str(area(kk)) ' mm2 and circulation=' num2str(circulation(kk)*1e-6) ' m2/s'])
                if kk~=1
                   disp('IDvortex: Found vortex not in the largest area found.')
				   info.G1peakInLargestG2peak(ii)=false;
				else
				   info.G1peakInLargestG2peak(ii)=true;
                end
			end
            break;
        end
    end
    if ~any(in) %all 'in' are false, the location is not in one of the areas
        disp('IDvortex: A found vortex could not be matched to an area.')
        vortex(ii,:) = [xLoc(ii), yLoc(ii), NaN, NaN];
    end
end
end

