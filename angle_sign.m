function sign_data=angle_sign(rel_pos,v_n,r0_vec)

% Given the input relative position and velocity vector define the sign of
% the angle being measured 
i_cap =v_n;

j_cap= cross(i_cap,r0_vec);
j_cap= j_cap./vecnorm((j_cap)')'; % find  unit vector normal to the plane made by rel pos, v_n

k_cap= cross(i_cap,j_cap);
k_cap= k_cap./vecnorm((k_cap)')'; % find  unit vector normal to the plane made by rel pos, v_n

sign_data=sign(dot(rel_pos,k_cap,2));
% normal = [0,0,1];
% normal = repmat(normal, size(rel_pos,1), 1);
% sign_data = sign(dot(normal, cross(rel_pos,v_n), 2));