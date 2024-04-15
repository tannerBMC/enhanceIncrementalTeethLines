% find rotation around horizontal axis (X) which optimizes
% contrast of vertical 1D profiles of mean projection image
%
% [projI,bestRot,bestVal,opt0,centerVal]=optimizeProjection(roi3DI,maxRot,deltaRot,showFigs,showExtraFigs)
%
% INPUT
% -----
% roi3DI           % 3D input image
% maxRot           % for search range -maxRot:deltaRot:maxRot in [deg]
% deltaRot         % [deg]
% showFigs         % 1: show figures
% showExtraFigs    % 1: show more figures
%
% OUTPUT
% -------
% projI           % best projection image
% bestRot         % best rotation angle [deg]
% bestVal         % best optimization value
% opt0            % optimization value for 0 deg projection
% centerVal       % optimization value for center slice
%
% This function was used in 
%
% Tanner et al., Extended-field synchrotron microtomography for 
% non-destructive analysis of incremental lines in archeological human 
% teeth cementum, Proceedings of SPIE 11840 (2021) 1184019,
% DOI:10.1117/12.2595180
% 
% Christine Tanner, University of Basel, Biomaterials Science Center
%
% v1, 2024/04/15, cleaned up version

function [projI,bestRot,bestVal,opt0,centerVal]=optimizeProjection(roi3DI,maxRot,deltaRot,showFigs,showExtraFigs)

pVal = 2;   % image normalization 2%

lw=2;       % linewidth
fs=12;      % fontsize
ms=9;       % marker size
    
% optimization function
% maximize standard deviation after projection

% define region and lines
sz=size(roi3DI);
dX=60;   % distance from border
dY=30;   % distance between lines
selX=dX:sz(1)-dX;
if isempty(selX)
    % reduce distance from border
    selX=dX/2:sz(1)-dX/2;
end
selY=dY:dY:sz(2);

% optimization value for projection at 0 degree
proj0I=mean(roi3DI,3);
projI=proj0I;   
opt0=mean(std(projI(selX,selY)));
optStr = 'std';

% optimization value for center slice
centerVal = mean(std(roi3DI(selX,selY,ceil(sz(3)/2))));

p2=prctile(projI(:),2);
p98=prctile(projI(:),98);
meanVal=mean(roi3DI(:));

if showExtraFigs
    figure(4)
    show3D(roi3DI)

    figure(76)
    imagesc(projI,[p2 p98])
    axis equal, axis tight
    title('proj 0 deg')
    colormap('gray')
end

%% optimize w.r.t. std along lines

rotAV=-maxRot:deltaRot:maxRot;
noRot=length(rotAV);

noLines=length(selY);
optM=zeros(noLines,noRot);
meanOptV=zeros(noRot,1);

if showExtraFigs
    figure(208),clf
    m=floor(noLines/2);  % show central 1D profile 
end

% contrasting color map
colorM=redblueu(noRot,[rotAV(1),rotAV(end)],'black');
%
% create horizontal colorbar
%
if showFigs
    tmpI(:,1,1)=colorM(:,1);
    tmpI(:,1,2)=colorM(:,2);
    tmpI(:,1,3)=colorM(:,3);
    figure(33)
    imagesc(tmpI)
    colormap(colorM)
    cb=colorbar('horizontal');
    cb.Ticks=0:0.1:1;
    cbTickLabelsV=get(cb,'TickLabels');
    for k=1:length(cbTickLabelsV)
        % in deg
        cb.TickLabels{k}=[num2str(str2num(cbTickLabelsV{k})*2*maxRot-maxRot,'%.0f') '^o'];
    end
    set(cb,'FontSize',18)
end
for k=1:noRot
    rotA=rotAV(k);
    
    % rotate image
    tmpI=imrotate3(roi3DI,rotA,[1 0 0],'crop','cubic','FillValues',meanVal);
    % mean projection
    projI=mean(tmpI,3);
    
    % mean standard deviation along lines
    optM(:,k)=std(projI(selX,selY));
    meanOptV(k)=mean(optM(:,k));
    % show progress
    fprintf('.')
    
    if showExtraFigs
        figure(207)
        imagesc(projI,[p2 p98])
        axis equal, axis tight
        title(['proj ' num2str(rotA) 'deg'])
        colormap('gray')
    
        figure(209)
        show3D(tmpI)
    
        drawnow
        pause(1)
    
        figure(208)
        subplot(3,1,1)
        hold on
        plot(projI(selX,selY(m)))
        hold off
        axis tight
        %ylim([p2 p98])
        title([num2str(meanOptV(k)) ' ' num2str(rotA) 'deg']);
        drawnow
    end
end

if showExtraFigs
    % show result without rotation
    projI=mean(roi3DI,3);
    figure(208)
    subplot(3,1,2)
    hold on
    plot(projI(selX,selY(m)),'r-')
    hold off
    axis tight
    %ylim([p2 p98])
    figure(208)
    subplot(3,1,3)
    hold on
    plot(projI(selX,selY(2*m)),'r-')
    hold off
    axis tight
end
if showFigs
    % add mean
    optM(noLines+1,:) = mean(optM);
    % Fig. 5d in manuscript 
    figure(220)
    subplot(2,1,1)
    imagesc(optM)
    xlabel('rotations')
    ylabel('1D profiles')
    hold on
    [~,maxIdx]=max(optM');
    plot(maxIdx,1:noLines+1,'r*')
    hold off
    title([optStr ' of 1D profiles'])
end

% best rotation for all 1D lines
[bestVal,bestIdx]=max(meanOptV);
bestRot=rotAV(bestIdx);
zeroC=colorM(ceil(noRot/2),:);
bestC=colorM(bestIdx,:);
    
% best rotation
tmpI=imrotate3(roi3DI,bestRot,[1 0 0],'crop','cubic','FillValues',meanVal);
projI=mean(tmpI,3);

if showFigs
    % Fig. 5d in manuscript 
    figure(220)
    subplot(2,1,2)
    fh=plot(rotAV,meanOptV,'g-');
    set(fh,'LineWidth',lw)
    hold on
    fh=plot(0,opt0,'*');
    set(fh,'MarkerSize',ms,'LineWidth',lw,'Color',zeroC)
    fh=plot(bestRot,bestVal,'o');
    set(fh,'MarkerSize',ms,'LineWidth',lw,'Color',bestC)
    hold off
    axis tight
    grid on
    fh=legend('function',['zero: ' num2str(opt0,'%.1f')],['best: ' num2str(bestVal,'%.1f')]);
    set(fh,'FontSize',fs*2)
    xlabel([])
    xticks([])
    fh=ylabel('measure'); %['mean(' optStr '(1Dprofiles))'])
    set(fh,'FontSize',fs*2)
    set(gca,'FontSize',fs*2)
end

if showExtraFigs
    figure(210)
    show3D(tmpI)
    title(['best ' num2str(bestRot) 'deg'])

    figure(211)
    imagesc(projI,[p2 p98])
    hold on
    plot(selY,selX(1)*ones(noLines),'g*')
    plot(selY,selX(end)*ones(noLines),'g*')
    hold off
    axis equal, axis tight, axis off
    title(['best ' num2str(bestRot) 'deg'])
    colormap('gray')
end

%% show 1D profiles

noX=length(selX);
proj0SD=std(proj0I(selX,selY));
projSD=std(projI(selX,selY));
tmp=proj0I(selX,selY);
profiles1p=prctile(tmp(:),1);
profiles99p=prctile(tmp(:),99);
B=15;

for midx=1:length(selY)-1
    if showFigs
        if midx<3
            % first 2 lines for Fig. 5 in manuscript
            figure(212)
            subplot(1,5,midx)
            
            fh=plot(proj0I(selX,selY(midx)),1:noX);
            set(fh,'LineWidth',lw,'Color',zeroC)
            hold on
            fh=plot(projI(selX,selY(midx)),1:noX);
            set(fh,'LineWidth',lw,'Color',bestC)
            hold off
            axis tight
            axis ij
            xlim([profiles1p,profiles99p])
            fh=legend(num2str([proj0SD(midx); projSD(midx)],'%.1f'));
            set(fh,'FontSize',fs*1.5)
            xticks([])
            yticks([])
            if midx==1
                fh=xlabel('Intensity','fontSize',fs*2);
                fh.Position = fh.Position + [50 0 0];
            end
        end
    end
    if showExtraFigs
        if midx>length(selY)-6
            % first 5 lines
            figure(213)
            subplot(1,5,midx-length(selY)+6)
            
            fh=plot(proj0I(selX,selY(midx)),1:noX);
            set(fh,'LineWidth',lw,'Color',zeroC)
            hold on
            fh=plot(projI(selX,selY(midx)),1:noX);
            set(fh,'LineWidth',lw,'Color',bestC)
            hold off
            axis tight
            axis ij
            xlim([profiles1p,profiles99p])
            fh=legend(num2str([proj0SD(midx); projSD(midx)],'%.1f'));
            set(fh,'FontSize',fs)
            fh=title(['x=' num2str(selY(midx))]);
            set(fh,'FontSize',fs)
            xticks([])
            yticks([])
            if midx==1
                ylabel('y','fontSize',fs)
            end
            xlabel('Intensity','fontSize',fs)
        end
    end
end
if showExtraFigs
    %
    % show 0, best rotation image with overlayed lines
    %
    figure(210)
    subplot(2,1,1)
    tmpRgb(:,:,1)=normalizeImagePercentile(proj0I,pVal,1);
    tmpRgb(:,:,2)=tmpRgb(:,:,1);
    tmpRgb(:,:,3)=tmpRgb(:,:,1);
    tmpI=addborder(tmpRgb,7,zeroC,'inner');
    imagesc(tmpI)
    colormap('gray')
    axis equal, axis tight, axis off
    
    subplot(2,1,2)
    tmpRgb(:,:,1)=normalizeImagePercentile(projI,pVal,1);
    tmpRgb(:,:,2)=tmpRgb(:,:,1);
    tmpRgb(:,:,3)=tmpRgb(:,:,1);
    tmpI=addborder(tmpRgb,7,bestC,'inner');
    imagesc(tmpI)
    colormap('gray')
    axis equal, axis tight, axis off
end

if showFigs
    % Fig. 5 a,b in manuscript
    %
    % show 0, best rotation image with overlayed lines
    %
    figure(211)
    subplot(2,1,1)
    tmpRgb(:,:,1)=normalizeImagePercentile(proj0I,pVal,1);
    tmpRgb(:,:,2)=tmpRgb(:,:,1);
    tmpRgb(:,:,3)=tmpRgb(:,:,1);
    tmpI=addborder(tmpRgb,7,zeroC,'inner');
    imagesc(tmpI)
    hold on
    for midx=1:length(selY)
        fh=line([selY(midx) selY(midx)],[selX(1),selX(end)]);
        set(fh,'Color',zeroC,'LineWidth',lw)
    end
    hold off
    colormap('gray')
    axis equal, axis tight, axis off
    
    subplot(2,1,2)
    tmpRgb(:,:,1)=normalizeImagePercentile(projI,pVal,1);
    tmpRgb(:,:,2)=tmpRgb(:,:,1);
    tmpRgb(:,:,3)=tmpRgb(:,:,1);
    tmpI=addborder(tmpRgb,7,bestC,'inner');
    imagesc(tmpI)
    hold on
    for midx=1:length(selY)
        fh=line([selY(midx) selY(midx)],[selX(1),selX(end)]);
        set(fh,'Color',bestC,'LineWidth',lw)
    end
    hold off
    colormap('gray')
    axis equal, axis tight, axis off
end

