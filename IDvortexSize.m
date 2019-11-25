function [area, circulation, points] = IDvortexSize(x,y,Vx,Vy)

% IDvortex size finder
% Detects the vortex size in Vx and Vy avg velocity fields
%   supply:
%       x,y coordinate vectors (mm) and Vx and Vy fields (m/s)
%   returns:
%       area in unit^2 (e.g. mm^2)
%       optional: circulation for every area (unit^2/s e.g. mm^2/s)
%       optional: coordinates that make up the area

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
% Plot the Gamma2 field with the dendrogram of qualifying locations
% Qualifying points plotted per area.

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
    Gamma2=0*Vx;
    for thisx=(side+1):length(x)-(side+1)
        for thisy=(side+1):length(y)-(side+1)
            if ~(Vx(thisy,thisx)==0 && Vy(thisy,thisx)==0) % ==0 where img cropped
                CenterLocation = [x(thisx) y(thisy)];
                xcomps = [...
                    reshape(Vx(thisy-side:thisy-1,thisx-side:thisx+side),[],1);...
                    reshape(Vx(thisy+1:thisy+side,thisx-side:thisx+side),[],1);...
                    reshape(Vx(thisy,thisx-side:thisx-1),[],1);...
                    reshape(Vx(thisy,thisx+1:thisx+side),[],1)];
                ycomps = [...
                    reshape(Vy(thisy-side:thisy-1,thisx-side:thisx+side),[],1);...
                    reshape(Vy(thisy+1:thisy+side,thisx-side:thisx+side),[],1);...
                    reshape(Vy(thisy,thisx-side:thisx-1),[],1);...
                    reshape(Vy(thisy,thisx+1:thisx+side),[],1)];
                UP = [mean(xcomps) mean(ycomps)];
                m=1;
                s=zeros(1,(n*n)-1);
                for AAx=thisx-side:thisx+side
                    for AAy=thisy-side:thisy+side
                        if ~(AAx==thisx && AAy==thisy)
                            mLocation = [x(AAx) y(AAy)];
                            P = CenterLocation-mLocation;
                            U = [Vx(AAy,AAx) Vy(AAy,AAx)];
                            UMP = U-UP;
                            s(m)=(P(1)*UMP(2)-P(2)*UMP(1)) / (norm(P)*norm(UMP));
                            m=m+1;
                        end
                    end
                end
                Gamma2(thisy,thisx)=mean(s);
            end
        end
    end
    
    [xMesh, yMesh] = meshgrid(x,y);
    [VortField,~] = curl(xMesh,yMesh,Vx*1000,Vy*1000); %mm, mm, mm/s, mm/s
    
    % Find all Gamma2>2/pi:
    threshold = 2/pi;
    vl = find(abs(Gamma2)>threshold);
    if SETTINGS.IDvortex.TroubleShooting==true
        %     Gamma2(Gamma2<threshold)=0; % optional
        figure('name',['Gamma2 field, n=' num2str(n)]); surf(xMesh,yMesh,Gamma2);
        xlabel('x (mm)'); ylabel('y (mm)')
        if isfield(SETTINGS.IDvortex,'xlim')
            xlim(SETTINGS.IDvortex.xlim)
        end
        if isfield(SETTINGS.IDvortex,'ylim')
            ylim(SETTINGS.IDvortex.ylim)
        end
    end
    if isempty(vl)
        if SETTINGS.IDvortex.quick==false
            n=n+2;
            disp(['IDvortexSize: Changed AA size to ' num2str(n) ' by ' num2str(n)])
            if n>25
                disp('IDvortexSize: No Gamma2 found. No vortex found !!!')
                area = NaN;
                circulation = NaN;
                points{1} = NaN;
                done = true;
            end
        else % quick mode: don't try other AA sizes
            disp('IDvortexSize: No Gamma2 found. No vortex found in quick mode.')
            area = NaN;
            circulation = NaN;
            points{1} = NaN;
            done = true;
        end
    else
        data=[xMesh(vl) yMesh(vl) Gamma2(vl) VortField(vl)];
        
        links = linkage(pdist(data(:,1:2)));
        cutoffvalue = 1.0*max([abs(max(diff(x))); abs(max(diff(y)))]); % immediately adjecent points are grouped
        data(:,5)=cluster(links,'cutoff',cutoffvalue,'criterion','distance');
        
        if SETTINGS.IDvortex.TroubleShooting==true
            figure('name',['IDvortexSize: Dendrogram of found points, n= ' num2str(n)]);
            dendrogram(links); hold on; plot(xlim,[cutoffvalue cutoffvalue],'k')
        end
        
        for ClusterNo=unique(data(:,5))'
            dataset{ClusterNo}=data(data(:,5)==ClusterNo,1:2);
            
            % Every point is in the center of a dA
            dA = mean(diff(x))*mean(diff(y));
            area(ClusterNo) = dA * length(dataset{ClusterNo});
            
            % Circulation:
            circulation(ClusterNo)=sum(sum(data(data(:,5)==ClusterNo,4).*dA)); % mm^2/s
            
            if SETTINGS.IDvortex.mute==false
                disp(['IDvortexSize: Found ' num2str(length(dataset{ClusterNo}(:,1))) ' points within ' num2str(area(ClusterNo)) 'mm^2.'])
            end
            
            if SETTINGS.IDvortex.TroubleShooting==true
                figure('name',['Outline of vortex # ' num2str(ClusterNo)]);
                plot(dataset{ClusterNo}(:,1),dataset{ClusterNo}(:,2),'k.');
                title(['Integrated area = ' num2str(area(ClusterNo))])
            end
        end
        [area, I]=sort(area,'descend'); % largest vortex first
        points = dataset(I);
        circulation = circulation(I);
        done=true;
    end
end
end

