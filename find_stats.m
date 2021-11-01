function [x,y]=find_stats(pts)
    x=unique(pts(:,1));
    y=zeros(size(x,1),3);
    for ii=1:size(x,1)
        y(ii,1)=mean(pts(pts(:,1)==x(ii),2));
        y(ii,2)=max(pts(pts(:,1)==x(ii),2));
        y(ii,3)=min(pts(pts(:,1)==x(ii),2));
    end 
end