% extractCementumPatchesEnhanceIL
% 
% extracts tooth cementum patches and enhances incremental line appearance
% automatically
%
% 1) extract patches from the outer band containing the incremental teeth 
% lines from the image of the whole tooth 
%
% - projection of 3D stack
% - median filter 15x
% - downsample 15x
% - imclose
% - mask
% - find outer ring
% - get centerline
% - upsample centerline position to original positions
%
% 2) straighten patches
% - extract stripe around spline fit to centerline points
%
% 3) automatic enhancement of incremental lines
% - optimize rotation to improve SNR of projection orthogonal to incremental lines
%
% This program was used in 
%
% Tanner et al., Extended-field synchrotron microtomography for 
% non-destructive analysis of incremental lines in archeological human 
% teeth cementum, Proceedings of SPIE 11840 (2021) 1184019,
% DOI:10.1117/12.2595180
% 
% Christine Tanner, University of Basel, Biomaterials Science Center
%
% v1, 2024/04/15  cleaned up version

close all
clear all

%% load data path and parameter definitions

dataParameterDefinitionTeeth

%% process all center slices in zV

for z = zV
    minZ=max(1,z-noSladd);
    maxZ=z+noSladd;
    
    slV=minZ:maxZ;
    noSl=length(slV);
    % central slice
    cidx=ceil(noSl/2);  %
    
    % read first tif slice
    fname = sprintf('%sreco_%04d.tif',imageDir,z);
    patchOutFname = sprintf('%sreco_%s%04d_patch',outDir,projStr,z);
    
    if ~exist(outDir,'dir')
        unix(['mkdir ' outDir])
    end
    figDir = sprintf('%sreco_%s%04d_',outDir,projStr,z);
    
    % where final result is saved
    roiStr = ['_pLen' num2str(pLen) '_y' num2str(years) ];
    saveFile=[patchOutFname roiStr '_optRot.mat'];
    
    if exist(saveFile,'file') && doSingle==0
        load(saveFile)
        disp([saveFile ' loaded'])
    else
        
        %% read 1 slice for size
        %
        disp(fname)
        I=imread(fname);
        sz=size(I);
        %
        % read all slices
        %
        fullI=zeros(sz(1),sz(2),noSl);
        disp(['reading ' num2str(noSl) ' slices: ' num2str(slV(1)) '-' num2str(slV(end))])
        for k=1:noSl
            fname = sprintf('%sreco_%04d.tif',imageDir,slV(k));
            try
                I=imread(fname);
                fullI(:,:,k)=I;
            catch
                % slice not available
                disp(['reading ' fname ' failed'])
                return
            end
        end
        
        % simple projection image of full image
        I=mean(fullI,3);
        
        if showFigs
            % projection normalized to 2-98 perctile intensity
            normI=normalizeImagePercentile(I,2,1);
            % normalized first slice
            normSl=normalizeImagePercentile(fullI(:,:,1),2,1);
            
            figure(1)
            imagesc(normI), axis equal, axis tight, axis off
            colormap('gray')
            if saveFigs
                imwrite(normI,[figDir 'meanProj_normIp2.png'])
            end
            figure(32)
            imagesc(normSl), axis equal, axis tight, axis off
            colormap('gray')
            if saveFigs
                imwrite(normSl,[figDir 'normIp2.png'])
            end
            
            if doScaleBars && saveFigs
                scaleBar = ones(300,sz(2));
                % 1000um = 1mm scale bar
                xmax = round(1000/pixelDimensions(1));
                
                scaleBar(101:200,end-xmax-100:end-100)=0;
                figure(77)
                imagesc(scaleBar),axis equal, axis tight,axis off
                colormap('gray')
                imwrite(scaleBar,[figDir 'scalebar1000um.png'])
            end
        end
        
        %% median filter image and bin it by factor F
        
        fI=medfilt2(I,[F F]);
        
        if showExtraFigs
            figure(2)
            imagesc(fI)
            axis equal, axis tight, axis off
            colormap('gray')
            title(['bin' num2str(F) '-medianF'])
        end
        if saveExtraFigs
            imwrite(normalizeImagePercentile(fI,2,1),[figDir 'medfilterIp2.png'])
        end
        
            
        % bin median filtered image
        lowSz=floor(sz/F);
        lowI=imresize(fI(1:lowSz(1)*F,1:lowSz(2)*F),lowSz);
        
        if showExtraFigs
            figure(20)
            imagesc(lowI)
            axis equal, axis tight, axis off
            colormap('gray')
            title(['medianF ' num2str(F) '-binned'])
        end
        if saveExtraFigs
            imwrite(normalizeImagePercentile(lowI,2,1),[figDir 'medfilterLowIp2.png'])
        end
        
        %% segment tooth
        %
        % thresholding median filtered downsampled image
        % and set background to zero
        
        % thresholding filtered, binned image
        G=2;  % foreground/background
        lowThV=multithresh(lowI(:),G);
        
        if showExtraFigs
            figure(3)
            subplot(2,1,1)
            histogram(single(lowI(:)),100)
            hold on
            plot(lowThV,zeros(G,1),'r*')
            hold off
            axis tight
            title(['medianF ' num2str(F) '-binned'])
            if saveExtraFigs
                print(3,'-dpng',[figDir 'thresholdHisto.png'])
            end
        end
        
        % select threshold index for background
        gIdx=G-1;
        if (sampleIdx==4 && hs==4)
            bgIdx=2;
        else
            bgIdx=1;
        end
        
        % threshold to remove BG from downsampled image
        lowTh=lowThV(bgIdx);
        fgLowI=lowI;
        fgLowI(fgLowI<lowTh)=0;
            
        if showFigs
            % later used
            figure(6)
            imagesc(fgLowI)
            axis equal, axis tight, axis off
            title(['(medianF ' num2str(F) '-binned < ' num2str(lowTh,'%.0f') ')->zero'])
        end
        if saveFigs
            % Fig. 3a in manuscript
            imwrite(fgLowI>0,[figDir 'fgLowIbinary.png'])
        end
        
        %% process foreground image
        
        % do morphological image closing
        R=10;    % default radius of sphere
        SE = strel('disk',R,0);  % same in 2D
        stripe255I=imclose(fgLowI>0,SE);
        if saveExtraFigs
            imwrite(stripe255I,[figDir 'fgLowIbinaryClosed.png'])
        end
           
        % determine largest connected body
        CC = bwconncomp(stripe255I>0);
        numPixels = cellfun(@numel,CC.PixelIdxList);
        [val,idx] = max(numPixels);
        maskI = zeros(size(stripe255I));
        maskI(CC.PixelIdxList{idx})=1;
        
        if showExtraFigs
            figure(61)
            subplot(2,3,1)
            imagesc(fgLowI)
            axis equal, axis tight
            title('fgLowI')
            subplot(2,3,2)
            imagesc(maskI)
            axis equal, axis tight
            title('1) thresholded,1B')
        end
        if saveExtraFigs
            imwrite(maskI,[figDir 'fgLowIbinaryClosed1B.png'])
        end
        
        if sampleIdx==2 && hs==4 && z>1900 
            % to avoid connecting with outer ring
            doErodeFill=1;
            erodeR=1;
        end
        if doErodeFill
            disp('doing imerode/largestComponent/imdilate/imfill');
            newR=erodeR;
            if newR==0
                stripe255I=fgLowI>0;
            else
                newSE = strel('sphere',newR);
                % erode
                stripe255I=imerode(fgLowI>0,newSE);
            end
            % largest connected body
            CC = bwconncomp(stripe255I>0);
            numPixels = cellfun(@numel,CC.PixelIdxList);
            [val,idx] = max(numPixels);
            maskIa = zeros(size(stripe255I));
            maskIa(CC.PixelIdxList{idx})=1;
            if newR==0
                maskIb=maskIa>0;
            else
                % dilate by same 
                maskIb=imdilate(maskIa>0,newSE);
            end
            % fill remaining holes
            maskI=imfill(maskIb>0,'holes');
            
            if showFigs
                figure(61)
                subplot(2,3,3)
                imagesc(stripe255I)
                axis equal, axis tight
                title(['2) imerode ' num2str(newR)])
                subplot(2,3,4)
                imagesc(maskIa)
                axis equal, axis tight
                title('3) 1Body')
                subplot(2,3,5)
                imagesc(maskIb)
                axis equal, axis tight
                title('4) imdilate')
                subplot(2,3,6)
                imagesc(fgLowI)
                hold on
                imcontour(maskI,'r');
                hold off
                axis equal, axis tight
                title('5) imfill')
            end
        end
        
        %% determine rim of tooth cementum
        
        % distance map to background
        [maskDbg,idxMaskDfg]=bwdist(1-maskI);
        
        % outer and inner ring
        stripe255I=maskDbg>distR;
        ringsI=maskI-stripe255I;
        if saveExtraFigs
             imwrite(ringsI,[figDir 'bothRingMasks.png'])
        end
        
        % find outer ring, largest component
        CC = bwconncomp(ringsI);
        numPixels = cellfun(@numel,CC.PixelIdxList);
        [~,largestIdx]=max(numPixels);
        outerRingMask=zeros(lowSz);
        outerRingMask(CC.PixelIdxList{largestIdx})=1;
        
        if showExtraFigs
            figure(7)
            imagesc(maskDbg)
            axis equal, axis tight, axis off
            colorbar
            title('distance to background')
            if saveExtraFigs
                print(7,'-dpng',[figDir 'maskDbg.png'])
            end
        end
       
        if showFigs
            % show ring mask on image
            figure(6)
            hold on
            imcontour(outerRingMask>0,1,'y--')
            hold off
            axis off
            colormap('gray')
            caxis([lowTh prctile(fgLowI(:),98)])
            if saveFigs
                print(6,'-dpng',[figDir 'outerRingMask.png'])
            end
        end
        
        if doErodeFill
            % new: make smoother for 'problem' cases
            newR=erodeR;
            newSE = strel('sphere',newR*3);
            % close
            tmpI=imclose(outerRingMask,newSE);
            
            if showExtraFigs
                figure(69)
                imagesc(tmpI)
                hold on
                imcontour(outerRingMask,'r')
                hold off
                axis equal, axis tight
            end
            outerRingMask=tmpI;
        end
        
        % upsample outer ring
        upOuterRingMask=imresize(outerRingMask,[lowSz(1)*F lowSz(2)*F]);
        
        
        %% extract strightened patches from outer ring of high-resolution image
        %
        % first get centerline from down-sampled image
        % then upsample to high-resolution image
        %
        clear Vertices
        %
        % extract bottom stripe of low-res image
        %
        % centerline of down-sampled ring image
        if doErodeFill
            % new: remove branches for 'problem' cases
            CL = bwskel(outerRingMask==1,'MinBranchLength',1000);
        else
            CL = bwskel(outerRingMask==1);
        end
        [ys,xs] = find(CL);  % center line points
        % center point
        mys=mean(ys);
        mxs=mean(xs);
        
        if showExtraFigs
            figure(18)
            imagesc(outerRingMask), axis equal, axis tight
            hold on
            plot(xs,ys,'.')
            plot(mxs,mys,'ro')
            hold off
            xlabel('xs')
            ylabel('ys')
        end
        
        % ensure centerline orientation
        %
        % bwskel ordered along xs, hence jumping between top and bottom
        % need to process sector of image
        % hence get angle to center
        % theta: Angular coordinate, returned as an array.
        % theta is the counterclockwise angle in the x-y plane measured in radians from the positive x-axis.
        % The value of the angle is in the range [-pi pi].
        [theta,rho] = cart2pol(xs-mxs,ys-mys);
        [thetaSort,polIdx]=sort(theta);
        rhoSort=rho(polIdx);
        % sorted in cartesian coordinates
        [tmpX,tmpY]=pol2cart(thetaSort,rhoSort);
        xsort=tmpX+mxs;
        ysort=tmpY+mys;
        % cumulative length along centerline
        noCL=length(xsort);
        dxs = xsort(2:end)-xsort(1:end-1);
        dys = ysort(2:end)-ysort(1:end-1);
        cumLenCL=cumsum(sqrt(dxs.^2+dys.^2));
        
        if showExtraFigs
            figure(180)
            subplot(2,2,1)
            plot(thetaSort,rhoSort,'x-')
            xlabel('theta')
            ylabel('rho')
            axis tight
            subplot(2,2,2)
            imagesc(outerRingMask), axis equal, axis tight
            hold on
            plot(xsort,ysort,'r-')
            hold off
            xlabel('xsort')
            ylabel('ysort')
            axis tight
        end
        
        % extract patches of a certain size, i.e. length along CL
        noP=floor(cumLenCL(end)/pLen);
        
        if showExtraFigs
            p=1;
            %
            % get patch 1 centerline to visualize
            %
            px=xsort(1+pLen*(p-1):1+pLen*p-1);
            py=ysort(1+pLen*(p-1):1+pLen*p-1);
            
            subplot(2,2,2+p)
            imagesc(outerRingMask), axis equal, axis tight
            hold on
            plot(px,py,'r-')
            hold off
            xlabel('xsort')
            ylabel('ysort')
            title(['lowPatch ' num2str(p)])
            axis tight
        end
        
        %%  centerline for high-resolution image
        %
        % CL with many branches, upsample coordinates instead
        upXs=xsort*F;
        upYs=ysort*F;
        
        % extract straightend stripe with width w
        w=lowW*F;   % width in pixels
        stripeL = pLen*F;
        
        % store performance of patches
        bestValV = zeros(noP,1);  % optimal value
        initialValV = zeros(noP,1);  % initial Value
        cxV=zeros(noP,1);  % center position of patch
        cyV=zeros(noP,1);
        centerValV=zeros(noP,1);
        bestRotV=zeros(noP,1);
    
        %% define patches to process
        
        if doSingle
            % patch of selected slice
            pV=selP:selP;
        else
            pV=1:noP;
        end
        
        if showFigs
            % Fig. 2b in manuscript
            figure(98)
            imagesc(normalizeImagePercentile(I,2,1))
            hold on
            imcontour(upOuterRingMask>0,1,'y-')
            % full centerline
            plot(upXs,upYs,'c-')
            hold off
            colormap('gray')
            axis off, axis equal, axis tight
            fh=legend({'rim','centerline'},'Location','northwest');
            set(fh,'FontSize',16)
            if saveFigs
                print(98,'-dpng',[figDir 'withRimCL.png'])
            end
        end
        %% loop to process all patches
        
        for p=pV
            %
            % get patch p, based on length along CL
            %
            % find points belonging to this lowRes CL region
            % do a bit of overlap and then cut down to same size
            if p<noP
                % use two neighbours, later cut down
                fidx1=find(cumLenCL>(pLen*(p-1)));
                fidx2=find(cumLenCL<(pLen*(p+1)));
                fidx=intersect(fidx1,fidx2);
                % for location of region, 2nd as first part moving clockwise
                fidx3=find(cumLenCL>(pLen*(p)));
                midx=intersect(fidx2,fidx3);
            else
                fidx1=find(cumLenCL>(pLen*(p-1)));
                % include start of CL
                fidx2=find(cumLenCL>0 & cumLenCL<pLen);
                fidx=[fidx1; fidx2];
                % for location of region, last sector
                midx=fidx2;
            end
            disp([num2str(1+pLen*(p-1)) ' to ' num2str(pLen*(p)) ' length ' num2str(length(fidx))])
            px=upXs(fidx);
            py=upYs(fidx);
            % save center location
            cxV(p)=mean(upXs(midx));
            cyV(p)=mean(upYs(midx));
            
            % range of points
            rx=max(px)-min(px);
            ry=max(py)-min(py);
            
            whalb=round(w/2);
            BB=[min(px)-whalb min(py)-whalb rx+2*whalb ry+2*whalb];
            cropI=imcrop(I,BB);
            cropM=imcrop(255*upOuterRingMask,BB);
            
            if showFigs
                % for Fig. 4 middle in manuscript
                figure(34)
                imagesc(normalizeImagePercentile(cropI,2,1))
                hold on
                plot(px-min(px)+whalb,py-min(py)+whalb,'r-')
                imcontour(cropM>127,1,'y-')
                hold off
                axis equal, axis tight, axis off
                colormap('gray')
            end
            
            % determine orientation of patch
            if p==noP
                % base it on first part for p==noP
                meanTheta=mean(thetaSort(fidx1));
            else
                meanTheta=mean(thetaSort(fidx));
            end
            if rx<ry/2  
                % range in x much smaller than in y
                % vertical left or right side
                [~,ix] = unique(py);
                xsub=px(ix);
                ysub=py(ix);
                % fit 2nd order poly
                polyV = polyfit(ysub,xsub,2);
                xpoly=polyval(polyV,ysub);
                if showFigs
                    figure(34)
                    hold on
                    plot(xpoly-min(px)+whalb,ysub-min(py)+whalb,'b-')
                    hold off
                    if saveFigs
                        print(34,'-dpng',[patchOutFname '_P' num2str(p,'%02d') '.png'])
                    end
                end
                
                % extract stripe
                
                % straighten does another fit to get  pixel-wise values
                % try using doing it here
                yfine=ysub(1):ysub(end);
                xfine=polyval(polyV,yfine);
                tmpI = straighten(I',[yfine',xfine'],w,showExtraFigs);
                tmpM = straighten(255*upOuterRingMask',[yfine',xfine'],w,showExtraFigs);
                % get straightened stripe for 3D image
                tmpFullI = straighten3D(permute(fullI,[2,1,3]),[yfine',xfine'],w,showExtraFigs);
                
                
                % patches with air at the bottom, lines horizontal
                if abs(meanTheta)>pi/2
                    % left side
                    stripeI= rot90(tmpI);
                    stripeRingMaskI= rot90(tmpM);
                    % get straightened stripe for 3D image
                    fullStripeI = rot90(tmpFullI);
                else
                    % right side
                    stripeI= tmpI';
                    stripeRingMaskI= tmpM';
                    % get straightened stripe for 3D image
                    fullStripeI = permute(tmpFullI,[2,1,3]);
                end
            else
                % horizontal patches, top or bottom
                % unique values in x
                [~,ix] = unique(px);
                xsub=px(ix);
                ysub=py(ix);
                % fit 2nd order poly
                polyV = polyfit(xsub,ysub,2);
                ypoly=polyval(polyV,xsub);
                if showFigs
                    figure(34)
                    hold on
                    plot(xsub-min(px)+whalb,ypoly-min(py)+whalb,'b-')
                    hold off
                    if saveFigs
                        print(34,'-dpng',[patchOutFname '_P' num2str(p,'%02d') '.png'])
                    end
                end
                % extract stripe
                
                % straighten does another fit to get  pixel-wise values
                % try using doing it here
                xfine=xsub(1):xsub(end);
                yfine=polyval(polyV,xfine);
                tmpI = straighten(I,[xfine',yfine'],w,showExtraFigs);
                tmpM = straighten(255*upOuterRingMask,[xfine',yfine'],w,showExtraFigs);
                % get straightened stripe for 3D image
                tmpFullI = straighten3D(fullI,[xfine',yfine'],w,showExtraFigs);
                 
                % patches with lines at the bottom
                if abs(meanTheta)>(pi*3/4)
                    % left side
                    if meanTheta<0
                        stripeI= rot90(tmpI);
                        stripeRingMaskI= rot90(tmpM);
                        % get straightened stripe for 3D image
                        fullStripeI = rot90(tmpFullI);
                    else
                        stripeI= rot90(tmpI,3);
                        stripeRingMaskI= rot90(tmpM,3);
                        % get straightened stripe for 3D image
                        fullStripeI = rot90(tmpFullI,3);
                    end
                elseif abs(meanTheta)<pi*1/2
                    % right side
                    if meanTheta>0
                        stripeI= tmpI';
                        stripeRingMaskI= tmpM';
                        % get straightened stripe for 3D image
                        fullStripeI = permute(tmpFullI,[2,1,3]);
                    else
                        stripeI= rot90(tmpI);
                        stripeRingMaskI= rot90(tmpM);
                        % get straightened stripe for 3D image
                        fullStripeI = rot90(tmpFullI);
                    end
                elseif meanTheta<0
                    % top patch
                    stripeI= rot90(tmpI);
                    stripeRingMaskI= rot90(tmpM);   
                    % get straightened stripe for 3D image
                    fullStripeI = rot90(tmpFullI);
                else
                    % bottom patch
                    stripeI= rot90(tmpI,3);
                    stripeRingMaskI= rot90(tmpM,3);
                    % get straightened stripe for 3D image
                    fullStripeI = rot90(tmpFullI,3);
                end
            end
            % write full strip to illustrate
            if showFigs
                % for Fig. 4 right bottom in manuscript
                figure(35)
                imagesc(normalizeImagePercentile(stripeI,2,1))
                hold on
                imcontour(stripeRingMaskI>127,1,'y-')
                hold off
                axis equal, axis tight, axis off
                colormap('gray')
                if saveFigs
                    print(35,'-dpng',[patchOutFname '_P' num2str(p,'%02d') 'fullStripe.png'])
                end
            end
            
            if stripeL<size(stripeI,2)
                % reduce stripeI to same size
                stripeI = stripeI(:,1:stripeL);
                stripeRingMaskI= stripeRingMaskI(:,1:stripeL);
                % get straightened stripe for 3D image
                fullStripeI = fullStripeI(:,1:stripeL,:);
            else
                disp('stripe is too short!')
            end
            centerI = fullStripeI(:,:,cidx);
            if showExtraFigs
                figure(88)
                imagesc(centerI)
                hold on
                imcontour(stripeRingMaskI>128,'y-')
                hold off
                axis equal, axis tight
                title(['p=' num2str(p) ', theta=' num2str(meanTheta,'%.2f')])
                colormap('gray')
            end
            if showFigs
                % Fig. 4 left in manuscript
                figure(98)
                imagesc(normalizeImagePercentile(I,2,1))
                hold on
                imcontour(upOuterRingMask>0,1,'y-')
                % full centerline
                plot(upXs,upYs,'c-')
                % selected part
                plot(xsub,ysub,'r-')
                %plot(cxV(p),cyV(p),'bx')
                hold off
                colormap('gray')
                axis off, axis equal, axis tight
                fh=legend({'rim','centerline','selected'},'Location','northwest');
                set(fh,'FontSize',16)
                if saveFigs
                    print(98,'-dpng',[patchOutFname '_P' num2str(p,'%02d') 'selectedCL.png'])
                end
                if doScaleBars
                    if p==1
                        scaleBar100 = ones(150,size(stripeI,2));
                        % 100um scale bar
                        xmax100= round(100/pixelDimensions(1));
                        scaleBar100(51:100,end-xmax100-10:end-10)=0;
                        imwrite(scaleBar100,[patchOutFname '_scaleBar100.png'],'png')
                    end
                end
            end
            %
            % write patch images
            %
            if saveData
                % Fig. 4 right top in manuscript
                imwrite(normalizeImagePercentile(stripeI,1,1),[patchOutFname '_' num2str(p,'%02d') '.png'],'png')
                imwrite(normalizeImagePercentile(fullStripeI(:,:,cidx),1,1),[patchOutFname '_' num2str(p,'%02d') '_centre.png'],'png')
            end
            
            %%
            % determine cementum subregion
            %
            % convert to 255
            roiI=normalizeImagePercentile(stripeI,1,255);
            stdV=std(roiI,[],2);  % standard deviation horizontal rows
            szRoi=size(roiI);
            
            roiMaskI=stripeRingMaskI;
            meanMaskV=mean(roiMaskI,2);  % mean horizontal rows
            roiMaskV=meanMaskV>128;
            [~,tmpIdx]=sort(abs(meanMaskV-128));  %
            maskIdx=sort(tmpIdx(1:5));
            locsMask = [maskIdx(1) maskIdx(end)];
            
            widthPx=ceil(years*2/pixelDimensions(1));
            topPx=locsMask(end)+20-widthPx;
            bottomPx=locsMask(end)+20;
            if topPx>0
                if bottomPx>szRoi(1)
                    dX=szRoi(1)-widthPx:szRoi(1);
                else
                    dX=topPx:bottomPx;
                end
            else
                % avoid negative index
                dX=1:widthPx+1;
            end
            if showExtraFigs
                figure(79)
                imagesc(roiI(dX,:))
                axis equal, axis tight
                colormap('gray')
                title(['(final) x: ' num2str(dX(1)) ' ' num2str(dX(end))])
                drawnow
            end
            
            %% rotate for best mean projection image
            
            
            % normalize intensities to 2-98th percentile
            roi3DI = normalizeImagePercentile(fullStripeI(dX,:,:),2,255);
            sz3Droi=size(roi3DI);
            cxRoi=floor(sz3Droi/2);
            
            roiStr = ['_P' num2str(p,'%02d') '_X' num2str(dX(1)) '_' num2str(dX(end))];
            disp(['processing ROI ' roiStr])
            
            figure(100)
            show3D(roi3DI)
            title('0 deg')
            meanVal=mean(roi3DI(:));
            
            if showExtraFigs
                % show slice planes
                tmpI=flip(permute(roi3DI,[3 2 1]),3);
                cxtmp=floor(size(tmpI)/2);
                figure(203)
                fh=slice(tmpI,[floor(cxtmp(2)*0.5) cxtmp(2) floor(cxtmp(2)*1.5)],[floor(noSl/2)],[],'cubic');
                for i1=1:length(fh)
                    fh(i1).EdgeColor = 'none';
                end
                xlabel('x')
                ylabel('y')
                zlabel('z')
                axis equal, axis tight,axis off
                view(3)
                colormap('gray')
                if saveExtraFigs
                    print(203,'-dpng',[patchOutFname roiStr '_roi3DPlanes.png'])
                end
                
                % zoom in
                roi3DIzoom = roi3DI(:,max(1,cxRoi(2)-cxRoi(1)):min(sz3Droi(2),cxRoi(2)+cxRoi(1)),:);
                figure(110)
                show3D(roi3DIzoom)
                title('0 deg')
                if saveExtraFigs
                    print(110,'-dpng',[patchOutFname roiStr '_roi3Dzoom.png'])
                end
                
                % show slice planes
                tmpI=flip(permute(roi3DIzoom,[3 2 1]),3);
                figure(202)
                fh=slice(tmpI,[cxRoi(1)],[floor(noSl/2)],[],'cubic');
                fh(1).EdgeColor = 'none';
                fh(1).MarkerFaceColor = 'r';
                fh(1).MarkerEdgeColor = 'r';
                fh(2).EdgeColor = 'none';
                xlabel('x')
                ylabel('y')
                zlabel('z')
                axis equal, axis tight,axis off
                view(3)
                colormap('gray')
                if saveExtraFigs
                    print(202,'-dpng',[patchOutFname roiStr '_roi3DzoomPlanes.png'])
                end
            end
            
            % get original 3D patch
            tightBB=[min(px) min(py) rx+whalb ry];
            
            crop3DI=imcrop3(fullI,[tightBB(1:2) 1 tightBB(3:4) noSl-1]);
            cxCropRoi=floor(size(crop3DI)/2);
            
            %%
            % show slice planes
            % for Fig 2 of manuscript
            if showFigs && sampleIdx==2 && p==17
                K=10;  % size of borders
                
                for loop=1:2
                    
                    if loop==1
                        tmpI=flip(permute(crop3DI(:,50:end-100,:),[3 1 2]),3);
                    else
                        % smaller region
                        tmpI=flip(permute(crop3DI(700:end,50:end-100,:),[3 1 2]),3);
                    end
                    tmpI=normalizeImagePercentile(tmpI,2,1);
                    cxtmp=floor(size(tmpI)/2);
                    figure(204)
                    if loop==1
                        % set slices, 3 across slab, 1 orthogonal
                        fh=slice(tmpI,[floor(cxtmp(2)*0.5) cxtmp(2) floor(cxtmp(2)*1.5)],[floor(3*noSl/4)],[],'cubic');
                    else
                        %fh=slice(tmpI,[floor(cxtmp(2))],[noSl-15],[],'cubic');
                        fh=slice(tmpI,[cxtmp(2)],[noSl-15],[],'cubic');
                    end
                    % add scale bar at bottom
                    xmax = round(100/pixelDimensions(1));
                    % scale bar at the end
                    if loop==1
                        fh(end).CData(end-20-xmax:end-20,20:20+K)=0;
                    else
                        % scale bar at the beginning
                        fh(end).CData(20:20+xmax,20:20+K)=0;
                    end
                    for i1=1:length(fh)
                        fh(i1).EdgeColor = 'none';
                        % make RGB images
                        rgbI(:,:,1) = fh(i1).CData;
                        rgbI(:,:,2) = fh(i1).CData;
                        rgbI(:,:,3) = fh(i1).CData;
                        fh(i1).CData = rgbI;
                        szC=size(fh(i1).CData);
                        if i1==length(fh)
                            % make green border
                            fh(i1).CData(1:K,:,:)=zeros(K,szC(2),3);
                            fh(i1).CData(1:K,:,2)=ones(K,szC(2));
                            fh(i1).CData(:,1:K,:)=zeros(szC(1),K,3);
                            fh(i1).CData(:,1:K,2)=ones(szC(1),K);
                            fh(i1).CData(end-K+1:end,:,:)=zeros(K,szC(2),3);
                            fh(i1).CData(end-K+1:end,:,2)=ones(K,szC(2));
                            fh(i1).CData(:,end-K+1:end,:)=zeros(szC(1),K,3);
                            fh(i1).CData(:,end-K+1:end,2)=ones(szC(1),K);
                        else
                            % planes in slab direction
                            % make a thin black line at the connection
                            k1=3;
                            if loop==1
                                offs=noSl-floor(3*noSl/4);
                            else
                                offs=15;
                            end
                            fh(i1).CData(end-k1-offs+1:end-offs,:,:)=zeros(k1,szC(2),3);
                            % make colored frame
                            fh(i1).CData(1:K,:,:)=ones(K,szC(2),3);
                            fh(i1).CData(1:K,:,2)=zeros(K,szC(2));
                            fh(i1).CData(:,1:K,:)=ones(szC(1),K,3);
                            fh(i1).CData(:,1:K,2)=zeros(szC(1),K);
                            fh(i1).CData(end-K+1:end,:,:)=ones(K,szC(2),3);
                            fh(i1).CData(end-K+1:end,:,2)=zeros(K,szC(2));
                            fh(i1).CData(:,end-K+1:end,:)=ones(szC(1),K,3);
                            fh(i1).CData(:,end-K+1:end,2)=zeros(szC(1),K);
                        end
                        clear rgbI
                    end
                    xlabel('x')
                    ylabel('y')
                    zlabel('z')
                    axis equal, axis tight, axis off
                    if loop==1
                        view(-55,13)
                    else
                        view(-57.5,20)
                    end
                    colormap('jet')
                    
                    if saveFigs
                        if loop==1
                            print(204,'-dpng',[patchOutFname roiStr '_crop3DPlanes.png'])
                        else
                            print(204,'-dpng',[patchOutFname roiStr '_crop3DPlanesSmaller.png'])
                        end
                    end
                end
            end
            
            %% optimize projection
            %
            fprintf('optimizeProjection');
            [projI,bestRotV(p),bestValV(p),initialValV(p),centerValV(p)]=optimizeProjection(roi3DI,maxRot,deltaRot,showFigs,showExtraFigs);
            disp('done')
            
            if abs(bestRotV(p))>(maxRot-deltaRot/2)
                % at border of range, search larger range once
                newMaxRot=maxRot*1.5;
                newDeltaRot=deltaRot*1.5;
                fprintf('optimizeProjection range +-40');
                [projI,bestRotV(p),bestValV(p),initialValV(p),centerValV(p)]=optimizeProjection(roi3DI,newMaxRot,newDeltaRot,showFigs,showExtraFigs);
                disp('done')
                extraStr='_30_15deg';
            else
                extraStr='';
            end
            if saveFigs
                % Fig. 5 colorbar
                print(33,'-dpng',[patchOutFname roiStr '_optRotMovieColorbar.png'])
                % Fig. 5a,b left in manuscript
                print(211,'-dpng',[patchOutFname roiStr '_1DprofilesOnImage.png'])
                % Fig. 5c in manuscript
                print(212,'-dpng',[patchOutFname roiStr '_1DprofilesFirstTwo.png'])
            end
            if saveExtraFigs
                % without lines overlaid
                print(210,'-dpng',[patchOutFname roiStr '_1DprofilesOnImageWithoutLines.png'])
                % last two profiles
                print(213,'-dpng',[patchOutFname roiStr '_1DprofilesLastFive.png'])
            end
            
            % show also central slice
            roiCenterI=centerI(dX,:);   % center slice
            proj0I=mean(roi3DI,3);      % no optimization
            
            if showFigs
                % Fig. 6 in manuscript
                figure(101)
                pVal=2;
                montage({normalizeImagePercentile(roiCenterI,pVal,1),normalizeImagePercentile(proj0I,pVal,1), normalizeImagePercentile(projI,pVal,1)},'Size',[3 1],'BorderSize',5, 'BackgroundColor','k')
                pause(1)
                drawnow
                title(['center ' num2str(centerValV(p),'%.1f') ', proj 0deg: ' num2str(initialValV(p),'%.1f') ', best ' num2str(bestRotV(p)) 'deg: ' num2str(bestValV(p),'%.1f')])
                if saveFigs
                    print(101,'-dpng',[patchOutFname roiStr '_optProj.png'])
                    print(220,'-dpng',[ patchOutFname roiStr '_optProjFunction.png'])
                end
            end
            
            rotI=imrotate3(roi3DI,bestRotV(p),[1 0 0],'crop','cubic','FillValues',meanVal);
            if showFigs
                figure(105)
                show3D(rotI)
                title(['best ' num2str(bestRotV(p)) 'deg'])
                % zoom in
                szRotI=size(rotI);
                cxRoi=floor(szRotI/2);
                rotIzoom = rotI(:,max(1,cxRoi(2)-cxRoi(1)):min(szRotI(2),cxRoi(2)+cxRoi(1)),:);
                figure(106)
                show3D(rotIzoom)
                title(['best ' num2str(bestRotV(p)) 'deg'])
                if saveFigs
                    print(106,'-dpng',[patchOutFname roiStr '_optRot3Dzoom.png'])
                end
            end
            
            % show slice planes
            if showExtraFigs
                tmpI=flip(permute(rotIzoom,[3 2 1]),3);
                figure(201)
                fh=slice(tmpI,[cxRoi(1)],[floor(noSl/2)],[],'cubic');
                fh(1).EdgeColor = 'none';
                fh(2).EdgeColor = 'none';
                xlabel('x')
                ylabel('y')
                zlabel('z')
                axis equal, axis tight,axis off
                view(3)
                colormap('gray')
                if saveExtraFigs
                    print(201,'-dpng',[patchOutFname roiStr '_optRot3DzoomPlanes.png'])
                end
                
                % rendering of best rotated slab
                normRoi3DI=normalizeImagePercentile(roi3DI,2,1);
                colorM=colormap('gray');
                alphaV = linspace(0,1.0,256)';
                figure(110)
                h=volshow(normRoi3DI,'BackgroundColor','k','Renderer','MaximumIntensityProjection','Alphamap',alphaV,'ColorMap',colorM,'ScaleFactors',2.3*ones(1,3));
                % capture image
                tmpF = getframe(gcf);
                figure(111)
                imshow(tmpF.cdata)
                if saveExtraFigs
                    print(111,'-dpng',[figDir 'slab_rendering_MIP_bestRot.png'])
                end
            end
        end
        %
        if doSingle
            disp('single patch processed')
        else
            % save optimization result
            save(saveFile,'bestValV','initialValV','w','pLen','years','cxV','cyV','centerValV','bestRotV')
        end
    end
end
