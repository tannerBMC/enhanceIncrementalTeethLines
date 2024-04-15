function show3D(I1,cx,I2,I3,I4)
    %
    % show image and potentially 3 segmentations
    % e.g. show3D(X,cx,mask1,mask2,mask3)
    % 
    if nargin==1
        % use center location for orthogonal slices
        cx=round(size(I1)/2);
    end
    
    ABCD = getOrthogonalSlices(I1,cx);
    imagesc(ABCD)
    
    if nargin>2
        ABCDY = getOrthogonalSlices(I2,cx);
        hold on
        imcontour(ABCDY,1,'y')
        hold off
    end
    if nargin>3
        ABCDY = getOrthogonalSlices(I3,cx);
        hold on
        imcontour(ABCDY,1,'c')
        hold off
    end
    if nargin>4
        ABCDY = getOrthogonalSlices(I4,cx);
        hold on
        imcontour(ABCDY,1,'b')
        hold off
    end
    colormap('gray')
    axis off
    axis equal
    axis tight
end

function ABCD = getOrthogonalSlices(X,cx)
    
    if ~isempty(X)
        sz = size(X);
        if nargin==1
            cx = round(sz/2);
        end
        [~,sortDim]=sort(sz);
        
        if sortDim(1)==2
            A = squeeze(X(:,cx(2),:));
            B = X(:,:,cx(3));
            C = squeeze(X(cx(1),:,:));
            D = zeros(sz(2),sz(2));
            
        elseif sortDim(1)==1
            A = squeeze(X(cx(1),:,:))';
            B = squeeze(X(:,cx(2),:))';
            C = X(:,:,cx(3));
            D = zeros(sz(1),sz(1));
        else
            A = X(:,:,cx(3))';
            B = squeeze(X(cx(1),:,:));
            C = squeeze(X(:,cx(2),:))';
            D = zeros(sz(3),sz(3));
        end
        ABCD = [A B; C D];
    else
        ABCD = [];
    end
end
