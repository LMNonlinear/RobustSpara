function A = create_empty_objarray(sz, a)

% A(2)=a(1);
% A=repmat(A(1),sz);
A = repmat(get_empty_obj(a), sz);
end
