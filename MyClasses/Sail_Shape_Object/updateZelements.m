function elements = updateZelements(Z, elements, n)
% Local function to update the Z elements of the mesh array for one sail
% quadrant mesh.
    ind = 0;
    for lvX = 1:n
        for lvY = 1:length(elements{lvX}(:))
            ind = ind+1;
            elements{lvX}{lvY}(3) = Z(ind);
        end
    end
end