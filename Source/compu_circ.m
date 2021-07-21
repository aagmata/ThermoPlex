function f = compu_circ(nsides,coord0)

coord_vec = coord0(1,:)-coord0(2,:);
degrot = @(nsides) [cosd((180*(nsides-2)./nsides)) -sind((180*(nsides-2)./nsides));...
        sind((180*(nsides-2)./nsides)) cosd((180*(nsides-2)./nsides))];
coord_set = zeros(nsides,2);
coord_set(1:2,:) = coord0;
    
for i = 1:nsides-2
    coord_set(i+2,:) = coord_set(i+1,:) + (degrot(nsides)*coord_vec')';
    coord_vec = coord_set(i+1,:)-coord_set(i+2,:);
end

f = coord_set;

end