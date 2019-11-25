function [ xLoc, yLoc] = IDvortexLoc(x,y,Vx,Vy)

% IDvortex location finder
% Detects a vortex location in Vx and Vy avg velocity fields
%   supply:
%       x,y coordinate vectors (mm) and Vx and Vy fields (m/s)
%   returns:
%       vortex position in true coordinates (e.g. mm) as opposed to
%       position in the vector.

% Based on the paper
% Combining PIV, POD and vortex identification algorithms for the study of
% unsteady turbulent swirling flows
% Laurent Graftieaux, Marc Michard and Nathalie Grosjean
% Ecole Centrale de Lyon

% Version #: 1.3, 2015-11-25
% sebastian@email.arizona.edu

%% Settings
MatLabSettings

% TroubleShooting
% Plot the Gamma1 field. Good to check if the found vortex is clearly
% visible and the only one.
% Displaying qualifying locations
% In general: qualifying locations should be few / in a small area.

% Avg area side: n is an odd number
% default is 11, sometimes as high as 21 was necessary
n=11;

if SETTINGS.IDvortex.TroubleShooting==true
    SETTINGS.IDvortex.mute=false;
end

%% Calculations
done=false;
while done==false
    side = (n-1)/2;
    Gamma1=0*Vx;
    for thisx=(side+1):length(x)-(side+1)
        for thisy=(side+1):length(y)-(side+1)
            if ~(Vx(thisy,thisx)==0 && Vy(thisy,thisx)==0) % ==0 where img cropped
                CenterLocation = [x(thisx) y(thisy)];
                m=1;
                s=zeros(1,(n*n)-1);
                for AAx=thisx-side:thisx+side
                    for AAy=thisy-side:thisy+side
                        if ~(AAx==thisx && AAy==thisy)
                            mLocation = [x(AAx) y(AAy)];
                            P = CenterLocation-mLocation;
                            U = [Vx(AAy,AAx) Vy(AAy,AAx)];
                            s(m)=(P(1)*U(2)-P(2)*U(1)) / (norm(P)*norm(U));
                            m=m+1;
                        end
                    end
                end
                Gamma1(thisy,thisx)=mean(s);
            end
        end
    end
    
    [xMesh, yMesh] = meshgrid(x,y);
    
    if SETTINGS.IDvortex.TroubleShooting==true
        figure('name',['Gamma1 field, n=' num2str(n)]); surf(xMesh,yMesh,Gamma1);
        xlabel('x (mm)'); ylabel('y (mm)')
        if isfield(SETTINGS.IDvortex,'xlim')
            xlim(SETTINGS.IDvortex.xlim)
        end
        if isfield(SETTINGS.IDvortex,'ylim')
            ylim(SETTINGS.IDvortex.ylim)
        end
    end
    
    % Find all Gamma1>0.95:
    t=SETTINGS.IDvortex.LocThreshold;
    vl = find(abs(Gamma1)>t);
    if isempty(vl)
        disp(['IDvortexLoc: No Gamma1>' num2str(t) ' found. Trying Gamma1>' num2str(t-0.05)])
        t=t-0.05;
        vl = find(abs(Gamma1)>t);
        if isempty(vl)
            disp(['IDvortexLoc: No Gamma1>' num2str(t) ' found. Trying Gamma1>' num2str(t-0.05)])
            t=t-0.05;
            vl = find(abs(Gamma1)>t);
        end
    end
    if isempty(vl)
        if SETTINGS.IDvortex.quick==false
            n=n+2;
            disp(['IDvortexLoc: Changed AA size to ' num2str(n) ' by ' num2str(n)])
            if n>25
                disp(['IDvortexLoc: No Gamma1 values > ' num2str(t) '. No vortex found.'])
                xLoc = NaN;
                yLoc = NaN;
                done = true;
            end
        else % quick mode: don't try other AA sizes
            disp(['IDvortexLoc: No Gamma1 values > ' num2str(t) '. No vortex found in quick mode.'])
            xLoc = NaN;
            yLoc = NaN;
            done = true;
        end
    else
        data=[xMesh(vl) yMesh(vl) Gamma1(vl)];
        
        if length(vl)>1
            links = linkage(pdist(data(:,1:2)));
            cutoffvalue = 1.0*max([abs(max(diff(x))); abs(max(diff(y)))]); % immediately adjecent points are grouped
            data(:,4)=cluster(links,'cutoff',cutoffvalue,'criterion','distance');
            
            if SETTINGS.IDvortex.TroubleShooting==true
                figure('name',['IDvortexLoc: Dendrogram of found points, n= ' num2str(n)]);
                dendrogram(links); hold on; plot(xlim,[cutoffvalue cutoffvalue],'k')
            end
        else
            data(:,4)=1; % only one entry in data
        end
        
        for ClusterNo=unique(data(:,4))'
            dataset=data(data(:,4)==ClusterNo,:);
            
            % weighted average of vortex' location
            xLoc(ClusterNo) = sum(dataset(:,1) .* dataset(:,3))/sum(dataset(:,3));
            yLoc(ClusterNo) = sum(dataset(:,2) .* dataset(:,3))/sum(dataset(:,3));
            dA = mean(diff(x))*mean(diff(y));
            area(ClusterNo) = dA * length(dataset);
            if SETTINGS.IDvortex.mute==false
                disp(['IDvortexLoc: Averaged ' num2str(length(dataset(:,1))) ' locations within ' num2str(area(ClusterNo)) 'mm^2.'])
            end
            
            if SETTINGS.IDvortex.TroubleShooting==true
                disp('IDvortexLoc: Locations found: ')
                disp('x, y, G1, cluster#')
                disp(dataset)
                disp(['Avg x = ' num2str(xLoc(ClusterNo))])
                disp(['Avg y = ' num2str(yLoc(ClusterNo))])
            end
        end
        % Sort by 'size', area here is not the vortex size tho!
        [foo, I]=sort(area,'descend');
        xLoc = xLoc(I);
        yLoc = yLoc(I);
        done=true;
    end
end
end

