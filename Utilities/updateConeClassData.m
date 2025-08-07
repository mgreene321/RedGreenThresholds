function new_Cone_Mat_all_marked = updateConeClassData(Cone_Mat_all_marked, newCones)

% Handle cones previously labeled as "missing" that now have types
new_Cone_Mat_all_marked = Cone_Mat_all_marked;
[~, ia, ib] = intersect(Cone_Mat_all_marked(:, 1:2), newCones(:,1:2), 'rows');
new_Cone_Mat_all_marked(ia,3) = newCones(ib,3);
% delete the labeled missing cones from newCones
newCones(ib,:) = [];
% concatenate newCones and new_Cone_mat_all_marked
new_Cone_Mat_all_marked = [new_Cone_Mat_all_marked; newCones];

end